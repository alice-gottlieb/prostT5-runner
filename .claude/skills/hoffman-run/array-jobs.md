# Array Jobs

An **array job** is a single `qsub` submission that fans out into N independent
*tasks*. The scheduler can place tasks on different nodes and run them fully in
parallel — the right tool for embarrassingly-parallel work (per-chunk,
per-genome, parametric sweeps). Each task is its own job with its own log, and a
single failed task can be rerun on its own without redoing the rest.

## Syntax

Add a `-t` directive to the job script:

```bash
#$ -t lower-upper:interval
```

- `#$ -t 1-49` — 49 tasks, indices 1..49 (interval defaults to 1).
- `#$ -t 200-480:2` — tasks at 200, 202, 204, … 480 (step size 2).

## Task environment variables

Inside the job script, each task reads these to know which slice it owns:

- `$SGE_TASK_ID` — this task's index (varies across the range).
- `$SGE_TASK_FIRST` — the lower bound.
- `$SGE_TASK_LAST` — the upper bound.
- `$SGE_TASK_STEPSIZE` — the interval.

## Pattern: replace the loop with `$SGE_TASK_ID`

Take a serial script that loops over inputs and **delete the `for` loop** — let
the scheduler do the looping. Use `$SGE_TASK_ID` as the file/chunk selector:

```bash
#$ -t 1-20
CHUNK_FILE="chunks/chunk_$(printf '%02d' "$SGE_TASK_ID").txt"
uv run python my_tool.py "$CHUNK_FILE" --output "out/task_${SGE_TASK_ID}"
```

The worked example `examples/submit_3di_multi_run_array.sh` goes one step
further: each task launches `PARALLEL_RUNS` runs at once on the same GPU and
sets the `-t` **step size equal to `PARALLEL_RUNS`** so the chunk ranges handled
by adjacent tasks don't overlap.

## Per-task logs

Give each task its own log by including `$TASK_ID` in the `-o` path:

```bash
#$ -o /u/scratch/a/aliceg/logs/{run}/joblog.$JOB_ID.$TASK_ID.out
#$ -j y
```

## Monitoring

```bash
ssh hoffman2 'qstat -u aliceg'          # tasks show as jobid.taskid
ssh hoffman2 'qdel {jobid}'             # kill the whole array
ssh hoffman2 'qdel {jobid} -t {range}'  # kill a subrange of tasks
```

## Resource note

The `#$ -l` request applies **per task**, not to the array as a whole. A
20-task array with `h_data=16G,gpu_mem=16G` asks for that much for *each* task,
so sizing mistakes multiply across the array. Pick GPU/CPU flags from
`resource-spec.md`; never combine `highp` with a GPU request.

See the hoffman2 docs (linked in `hoffman-help.md`) for more on job arrays.
