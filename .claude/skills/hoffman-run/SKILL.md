---
name: hoffman-run
description: Run, test, and interact with the Hoffman cluster. Use when asked to run, drive, remote test, verify, explore, profile, or describe the action of code on Hoffman. Also use when asked to run, drive, test, or verify any code on a GPU. 
---
# Run hoffman-run

Hoffman2 is a UCLA shared cluster running Univa/Altair Grid Engine (`qsub`,
`qrsh`, `qstat`). You drive it over SSH from this machine. The one rule that
shapes everything: **never run real compute on the login node** — submit a batch
job with `qsub` and poll it, or grab a node interactively with `qrsh`.

All commands below are run from the local machine; they SSH into the cluster for
you. Paths in `{curly braces}` are placeholders to fill in.

## Connect

SSH is key-based and non-interactive — run one-shot commands and read the output:

```bash
ssh hoffman2 'hostname && whoami'   # login node, user `aliceg`
```

Key locations on the cluster:
- **Repo:** `~/prostT5-runner` (push from here with git, pull there — see Transfer).
- **Scratch:** `$SCRATCH` = `/u/scratch/a/aliceg` — large, fast, **not permanent**
  (purged periodically). Job outputs, logs, and big data live here.
- `h2compute` (an SSH alias to `localhost:7022`) is a pre-existing tunnel to a
  compute node, for interactive use when one is already running.

## Run a job (qsub — primary path)

This is the default way to run anything non-trivial. The loop:

1. **Get the code there** (see Transfer) — usually `git pull` on the cluster.
2. **Write the job script.** Copy a file from `templates/`, fill the
   `{PLACEHOLDERS}`, and set the `#$ -l` resource line from `resource-spec.md`.
3. **Submit** and capture the job id:
   ```bash
   ssh hoffman2 'cd ~/prostT5-runner && qsub {my_job.sh}'
   # -> Your job 1234567 ("my_job.sh") has been submitted
   ```
4. **Poll — do not block.** Check state and read the log named in the `#$ -o`
   line of the script:
   ```bash
   ssh hoffman2 'qstat -u aliceg'          # or: myjobs -u aliceg
   ssh hoffman2 'tail -n 50 {OUTPUT_LOG}'  # the path from `#$ -o`
   ```
   States: `qw` = queued/waiting, `r` = running, `Eqw` = error (stuck — inspect,
   then `qdel` and resubmit). The log file only appears **after** the job starts.
5. **Cancel / inspect:**
   ```bash
   ssh hoffman2 'qdel {jobid}'
   ssh hoffman2 'qacct -j {jobid}'   # accounting/exit status for a finished job
   ```

For work that splits across many nodes (per-chunk, per-genome, parametric
sweeps), use an **array job** — see `array-jobs.md`. Worked example:
`examples/submit_3di_multi_run_array.sh`.

## Run interactively (qrsh — secondary)

For quick iteration on a GPU node. This **holds a session**, so use it only for
short interactive work; use `qsub` for anything long.

```bash
ssh hoffman2   # then, on the login node:
qrsh -l gpu,RTX2080Ti,cuda=1,h_rt=2:00:00,h_data=8G
# once on the node:
module load cuda/12.3
cd ~/prostT5-runner && uv run python {script.py}
```

## Check GPUs / job health

On a compute node (inside qrsh or from within a job script):

```bash
nvidia-smi
# machine-readable, e.g. pick the GPU with most free memory:
nvidia-smi --query-gpu=index,memory.free --format=csv,noheader,nounits
```

Cluster-wide GPU availability (run on the login node): `GPU_NODES_AT_A_GLANCE`
or `GPU_NODES_AT_A_GLANCE_WITH_USED_GPUs`. See `resource-spec.md`.

## Transfer files

Prefer **git** for code; `scp`/`rsync` for data and artifacts that aren't in git.
Results and scratch data are **not** committed.

The cluster checkout's working tree is **usually dirty** (many locally-modified
and untracked files). So always pull with `--ff-only`: it fast-forwards cleanly
when the incoming commits don't touch a locally-modified file, and otherwise
**refuses without merging** rather than leaving a conflicted tree. If `--ff-only`
is rejected, stop and ask the user — do not merge, stash, or force.

```bash
# code: commit/push locally, then pull on the cluster
ssh hoffman2 'cd ~/prostT5-runner && git pull --ff-only origin main'

# pull a results dir from scratch back to repo-local tmp/ (gitignored)
scp -r hoffman2:/u/scratch/a/aliceg/{run_dir} tmp/
# or rsync for big/resumable transfers
rsync -av hoffman2:/u/scratch/a/aliceg/{run_dir}/ tmp/{run_dir}/
```

## Resources

Full GPU table, highp (CPU-only) nodes, and gotchas are in `resource-spec.md`.
Two rules you need at submit time:
- Target **`RTX2080Ti`** for the fastest queue turnaround (`H100`/`H200` for big
  GPU memory). Request a GPU with `-l gpu,RTX2080Ti,cuda=1,...`.
- **Never combine `highp` with a GPU request.** `highp` is for CPU-only jobs.

## Templates & Examples

- `templates/qsub-basic-template.sh` — minimal single-core job.
- `templates/qsub-multithread-template.sh` — multi-core (`-pe shared N`); total
  memory = `h_data` × cores.
- `examples/submit_3di_multi_run_array.sh` — real multi-run GPU **array** job
  (illustrative; adapt paths/resources).

## What to NEVER do

- **Never run compute on the login node** — always `qsub` or `qrsh`.
- **Never combine `highp` with a GPU** request.
- **Never run `rm` on the server.** Do not delete cluster files; if cleanup is
  needed, ask the user.
- **Never call `rclone`.** Backups (e.g. `scratch-backup.sh`) are the user's to
  run, not yours.
- **Never commit** `results**`/scratch data, or hardcode another user's paths.
- **Never use bare `python`** on the cluster — always `uv run python` (bare
  python lacks the locked deps).

## Gotchas

- **`h_data` is per-core.** Total job memory = `h_data` × number of cores
  (`-pe shared N`). E.g. `h_data=2G` with `-pe shared 4` = 8G total.
- **Logs appear only after the job starts running** (state `r`), not while `qw`.
- **`$SCRATCH` is purged periodically** — don't treat it as permanent storage.
- **Slow `highp` start?** If a `highp` CPU job sits queued >5–15 min, resubmit
  without `highp`.

## Help / escalation

If you need more than this skill covers, see `hoffman-help.md` (man pages on the
cluster, the hoffman2 docs URL, and when to ask the user).
