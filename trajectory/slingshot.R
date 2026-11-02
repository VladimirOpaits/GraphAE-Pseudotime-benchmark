suppressPackageStartupMessages({
  library(SingleCellExperiment)
  library(slingshot)
  library(data.table)
  library(ggplot2)
  library(cowplot)
})

defaults <- list(
  coords = NULL, labels = NULL, start = NULL, outdir = NULL,
  `label-col` = NULL, end = NULL, dims = 0, subsample = 0, seed = 67,
  `approx-points` = 150, umap = NULL, time = NULL, `time-col` = NULL,
  `min-cells` = 30, prefix = "slingshot"
)
required <- c("coords", "labels", "start", "outdir")

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
  for (key in c("dims", "subsample", "seed", "approx-points", "min-cells")) {
    out[[key]] <- as.numeric(out[[key]])
  }
  out
}

opt <- parse_args(commandArgs(trailingOnly = TRUE), defaults, required)
dir.create(opt$outdir, recursive = TRUE, showWarnings = FALSE)
out_path <- function(suffix) file.path(opt$outdir, paste0(opt$prefix, "_", suffix))

read_table <- function(path) {
  d <- fread(path, data.table = FALSE, check.names = FALSE)
  rownames(d) <- as.character(d[[1]])
  d[, -1, drop = FALSE]
}

pick_col <- function(d, col, path) {
  if (is.null(col)) col <- colnames(d)[1]
  if (!col %in% colnames(d)) stop("Column '", col, "' not found in ", path)
  setNames(d[[col]], rownames(d))
}

spearman <- function(x, y) {
  ok <- complete.cases(x, y)
  if (sum(ok) < 3 || sd(x[ok]) == 0 || sd(y[ok]) == 0) return(NA_real_)
  cor(x[ok], y[ok], method = "spearman")
}

coords_df <- read_table(opt$coords)
coords <- as.matrix(coords_df[, vapply(coords_df, is.numeric, logical(1)), drop = FALSE])
if (opt$dims > 0) {
  if (opt$dims > ncol(coords)) stop("--dims exceeds available columns: ", ncol(coords))
  coords <- coords[, seq_len(opt$dims), drop = FALSE]
}

labels_df <- read_table(opt$labels)
labels <- pick_col(labels_df, opt$`label-col`, opt$labels)

time_all <- NULL
if (!is.null(opt$time)) {
  time_df <- read_table(opt$time)
  time_all <- pick_col(time_df, opt$`time-col`, opt$time)
  if (!is.numeric(time_all)) stop("Time column must be numeric")
}

umap_all <- NULL
if (!is.null(opt$umap)) {
  umap_df <- read_table(opt$umap)
  umap_all <- as.matrix(umap_df[, c("UMAP1", "UMAP2")])
}

cells <- intersect(rownames(coords), names(labels))
cells <- cells[!is.na(labels[cells])]
if (!is.null(umap_all)) cells <- intersect(cells, rownames(umap_all))
if (!length(cells)) stop("No shared cell ids between --coords and --labels")

set.seed(opt$seed)
if (opt$subsample > 0 && opt$subsample < length(cells)) {
  cells <- sample(cells, opt$subsample)
}

labels_sub <- as.character(labels[cells])
if (!opt$start %in% labels_sub) stop("--start '", opt$start, "' not found in labels")
end_clus <- if (is.null(opt$end)) NULL else strsplit(opt$end, ",", fixed = TRUE)[[1]]
if (length(setdiff(end_clus, labels_sub))) {
  stop("--end not found in labels: ", paste(setdiff(end_clus, labels_sub), collapse = ", "))
}

sce <- SingleCellExperiment(
  assays = list(counts = matrix(0, nrow = 1, ncol = length(cells),
                                dimnames = list("x", cells))),
  reducedDims = list(LATENT = coords[cells, , drop = FALSE])
)
colData(sce)$label <- labels_sub

sce <- slingshot(sce, clusterLabels = "label", reducedDim = "LATENT",
                 start.clus = opt$start, end.clus = end_clus,
                 approx_points = opt$`approx-points`)

pt <- as.data.frame(slingPseudotime(sce))
lineage_names <- colnames(pt)
lineages <- slingLineages(sce)
lineage_path <- vapply(lineages, paste, character(1), collapse = " > ")
terminal <- vapply(lineages, function(x) x[length(x)], character(1))

pt_out <- cbind(cell_id = cells, label = labels_sub, pt)
fwrite(pt_out, out_path("pseudotime.csv"))
writeLines(paste(names(lineage_path), lineage_path), out_path("lineages.txt"))

curves <- rbindlist(lapply(seq_along(slingCurves(sce)), function(i) {
  cv <- slingCurves(sce)[[i]]
  d <- as.data.frame(cv$s[cv$ord, ])
  d$lineage <- lineage_names[i]
  d$order <- seq_len(nrow(d))
  d
}))
fwrite(curves, out_path("curves_latent.csv"))

if (!is.null(time_all)) {
  time_sub <- unname(time_all[cells])
  rows <- list()
  for (lin in lineage_names) {
    on <- !is.na(pt[[lin]]) & !is.na(time_sub)
    rows[[length(rows) + 1]] <- data.table(
      lineage = lin, terminal = terminal[[lin]], path = lineage_path[[lin]],
      scope = "all", n = sum(on), spearman = spearman(pt[[lin]][on], time_sub[on])
    )
    for (cl in sort(unique(labels_sub[on]))) {
      m <- on & labels_sub == cl
      if (sum(m) < opt$`min-cells`) next
      rows[[length(rows) + 1]] <- data.table(
        lineage = lin, terminal = terminal[[lin]], path = lineage_path[[lin]],
        scope = cl, n = sum(m), spearman = spearman(pt[[lin]][m], time_sub[m])
      )
    }
  }
  metrics <- rbindlist(rows)
  within <- metrics[scope != "all" & !is.na(spearman),
                    .(scope = "within_weighted", n = sum(n),
                      spearman = weighted.mean(spearman, n)),
                    by = .(lineage, terminal, path)]
  metrics <- rbind(metrics, within)
  setorder(metrics, lineage, scope)
  fwrite(metrics, out_path("time_metrics.csv"))
  print(as.data.frame(metrics[scope %in% c("all", "within_weighted")]), row.names = FALSE)

  time_long <- rbindlist(lapply(lineage_names, function(lin) {
    on <- !is.na(pt[[lin]]) & !is.na(time_sub)
    data.table(lineage = lin, pseudotime = pt[[lin]][on], time = time_sub[on])
  }))
  time_long[, time_bin := if (uniqueN(time) <= 10) factor(time) else cut(frank(time, ties.method = "min"), 10)]
  p_time <- ggplot(time_long, aes(time_bin, pseudotime)) +
    geom_boxplot(outlier.size = 0.2) +
    facet_wrap(~lineage, scales = "free_y") +
    labs(x = "time") +
    theme_classic() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
  ggsave(out_path("time_vs_pseudotime.png"), p_time,
         width = 5 * length(lineage_names), height = 4, dpi = 200)
}

if (!is.null(umap_all)) {
  reducedDim(sce, "UMAP") <- umap_all[cells, , drop = FALSE]
  emb <- embedCurves(sce, "UMAP")
  curves_umap <- rbindlist(lapply(seq_along(slingCurves(emb)), function(i) {
    cv <- slingCurves(emb)[[i]]
    d <- as.data.frame(cv$s[cv$ord, ])
    colnames(d) <- c("UMAP1", "UMAP2")
    d$lineage <- lineage_names[i]
    d
  }))
  fwrite(curves_umap, out_path("curves_umap.csv"))

  df <- data.frame(umap_all[cells, , drop = FALSE], label = labels_sub, pt, check.names = FALSE)
  cent <- aggregate(cbind(UMAP1, UMAP2) ~ label, df, median)
  arr <- arrow(length = unit(0.25, "cm"), type = "closed")
  p_cl <- ggplot(df, aes(UMAP1, UMAP2, color = label)) +
    geom_point(size = 0.4, alpha = 0.5) +
    geom_path(data = curves_umap, aes(UMAP1, UMAP2, group = lineage),
              color = "black", linewidth = 1, arrow = arr, inherit.aes = FALSE) +
    geom_text(data = cent, aes(UMAP1, UMAP2, label = label),
              color = "black", size = 3, fontface = "bold", inherit.aes = FALSE) +
    theme_classic() + theme(legend.position = "none") +
    labs(title = "Clusters + lineages")
  plot_pt <- function(lin) {
    d <- df[order(!is.na(df[[lin]])), ]
    ggplot(d, aes(UMAP1, UMAP2, color = .data[[lin]])) +
      geom_point(size = 0.4, alpha = 0.6) +
      geom_path(data = curves_umap[curves_umap$lineage == lin, ], aes(UMAP1, UMAP2),
                color = "black", linewidth = 1, arrow = arr, inherit.aes = FALSE) +
      scale_color_viridis_c(option = "plasma", na.value = "grey88") +
      theme_classic() + labs(title = paste(lin, "pseudotime"), color = "PT")
  }
  panels <- c(list(p_cl), lapply(lineage_names, plot_pt))
  p_all <- plot_grid(plotlist = panels, ncol = length(panels))
  ggsave(out_path("umap.png"), p_all, width = 6 * length(panels), height = 7, dpi = 200)
}
