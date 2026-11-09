suppressPackageStartupMessages({
  library(slingshot)
  library(ggplot2)
  library(cowplot)
})

defaults <- list(
  x = NULL, n = NULL, d = NULL, clusters = NULL, starts = NULL, outdir = NULL,
  `approx-points` = 150, `plot-xy` = NULL, `plot-labels` = NULL, `plot-dir` = NULL,
  `plot-max-cells` = 30000, seed = 0, prefix = "slingshot"
)
required <- c("x", "n", "d", "clusters", "starts", "outdir")

parse_args <- function(argv, defaults, required) {
  out <- defaults
  i <- 1
  while (i <= length(argv)) {
    a <- argv[i]
    if (!startsWith(a, "--")) stop("Unexpected argument: ", a)
    a <- sub("^--", "", a)
    if (grepl("=", a, fixed = TRUE)) {
      key <- sub("=.*", "", a)
      val <- sub("^[^=]*=", "", a)
    } else {
      key <- a
      i <- i + 1
      if (i > length(argv)) stop("Missing value for --", key)
      val <- argv[i]
    }
    if (!key %in% names(defaults)) stop("Unknown option: --", key)
    out[[key]] <- val
    i <- i + 1
  }
  missing <- setdiff(required, names(Filter(Negate(is.null), out)))
  if (length(missing)) stop("Missing required: ", paste0("--", missing, collapse = ", "))
  for (key in c("n", "d", "approx-points", "plot-max-cells", "seed")) {
    out[[key]] <- as.numeric(out[[key]])
  }
  out
}

opt <- parse_args(commandArgs(trailingOnly = TRUE), defaults, required)
dir.create(opt$outdir, recursive = TRUE, showWarnings = FALSE)
out_path <- function(name) file.path(opt$outdir, name)

X <- matrix(readBin(opt$x, "double", opt$n * opt$d), nrow = opt$n, byrow = TRUE)
cl <- readLines(opt$clusters)
starts <- strsplit(opt$starts, ",", fixed = TRUE)[[1]]

xy <- NULL
if (!is.null(opt$`plot-xy`)) {
  xy <- matrix(readBin(opt$`plot-xy`, "double", opt$n * 2), nrow = opt$n, byrow = TRUE)
  colnames(xy) <- c("x", "y")
  lab <- if (is.null(opt$`plot-labels`)) cl else readLines(opt$`plot-labels`)
  plot_dir <- if (is.null(opt$`plot-dir`)) opt$outdir else opt$`plot-dir`
  dir.create(plot_dir, recursive = TRUE, showWarnings = FALSE)
  set.seed(opt$seed)
  show <- if (opt$n > opt$`plot-max-cells`) sort(sample(opt$n, opt$`plot-max-cells`)) else seq_len(opt$n)
}

plot_fit <- function(sds, pt, s) {
  emb <- embedCurves(sds, xy)
  curves <- do.call(rbind, lapply(seq_along(slingCurves(emb)), function(i) {
    cv <- slingCurves(emb)[[i]]
    d <- as.data.frame(cv$s[cv$ord, , drop = FALSE])
    colnames(d) <- c("x", "y")
    d$lineage <- paste0("L", i)
    d
  }))
  df <- data.frame(xy[show, , drop = FALSE], label = lab[show], cluster = cl[show],
                   pt = rank(pt, na.last = "keep")[show] / sum(!is.na(pt)))
  cent <- aggregate(cbind(x, y) ~ cluster, data.frame(xy, cluster = cl), median)
  arr <- arrow(length = unit(0.2, "cm"), type = "closed")
  curve_layer <- geom_path(data = curves, aes(x, y, group = lineage), color = "black",
                           linewidth = 0.8, arrow = arr, inherit.aes = FALSE)
  p_lab <- ggplot(df, aes(x, y, color = label)) +
    geom_point(size = 0.2, alpha = 0.4) + curve_layer +
    geom_text(data = cent, aes(x, y, label = cluster), color = "black", size = 2.5, inherit.aes = FALSE) +
    guides(color = guide_legend(override.aes = list(size = 2, alpha = 1))) +
    theme_classic() + labs(title = paste("labels, clusters and lineages; start", s), color = NULL)
  p_pt <- ggplot(df[order(df$pt), ], aes(x, y, color = pt)) +
    geom_point(size = 0.2, alpha = 0.6) + curve_layer +
    scale_color_viridis_c(option = "plasma", na.value = "grey88") +
    theme_classic() + labs(title = "pseudotime (rank)", color = "rank")
  ggsave(file.path(plot_dir, paste0(opt$prefix, "_start_", s, ".png")),
         plot_grid(p_lab, p_pt, ncol = 2), width = 14, height = 6, dpi = 150)
}

for (s in starts) {
  sds <- slingshot(X, clusterLabels = cl, start.clus = s, approx_points = opt$`approx-points`)
  pt <- slingAvgPseudotime(sds)
  writeBin(as.double(pt), out_path(paste0("pt_", s, ".bin")))
  lin <- slingLineages(sds)
  writeLines(paste(names(lin), vapply(lin, paste, character(1), collapse = " > ")),
             out_path(paste0("lineages_", s, ".txt")))
  if (!is.null(xy)) plot_fit(sds, pt, s)
}
