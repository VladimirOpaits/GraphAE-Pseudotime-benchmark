data.py - loads data; clone matrix and metadata
metrics.py - sister order - main metric for Larry dataset (acc, error, mean gap)
graph.py - to not recalculate graph every single time, and we gonna test different k- values and reuse them a lot
dpt.py - cli scrpit that tests diffusion pseudo time and dependency on start root (--root-seed)
trajectory.py - runs dpt.py 