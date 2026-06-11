# GPU Nodes #
In general, you should target H200 (only 1 gpu at a time is allowed), H100, or RTX2080Ti models (RTX models are the most available for quickest job start turnaround). These have the following characteristics (to access, add the portion of the resource column to your `-l` option in the submission file):

| GPU type  | Compute capability | # of CUDA cores | GPU memory     | scheduler options         |
| --------- | -----------------: | --------------: | -------------- | ------------------------- |
| H200      |                9.0 |            8448 | 141 GB         | `-l gpu,H200,cuda=1`      |
| H100      |                9.0 |            8448 | 96 GB          | `-l gpu,H100,cuda=1`      |
| RTX2080Ti |                7.5 |            4352 | 10 GB          | `-l gpu,RTX2080Ti,cuda=1` |

If you need more information on specific GPU nodes, run `GPU_NODES_AT_A_GLANCE` or `GPU_NODES_AT_A_GLANCE_WITH_USED_GPUs` on hoffman.

# Highp Nodes #
Use these nodes when running a **CPU-only job**. You can requestion them by adding `highp` to a resource request in a submission file (e.g. `-l h_rt=23:00:00,h_data=10G,highp`)

User aliceg has access to the following nodes in highp:

n1902 n6361 n7361 n7674

These nodes have the following characteristics:

CPU-type                # nodes         #cores/node     # tot. cores     memory/core (GB)        tot memory (GB)
intel-E5-2650v4          2              24               48              2.617                   62.800
intel-gold-6140          1              36               36              5.239                   188.600
intel-gold-6240          1              36               36              5.239                   188.600

TOTALS                  4               -               120              -                       440

# Gotchas #
* If a highp cpu task is taking a long time (more than 5-15min) to start, trying running the same job but removing `highp` from the resource list
* *Never* use `highp` when requesting any GPUs