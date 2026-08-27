# Running an RNA-seq Cohort on the HPC — Student Guide

A step-by-step guide for running **one RNA-seq cohort at a time** through this pipeline
on DISCOVERY. You give it a folder of paired FASTQ files; it produces a **gene counts
matrix** plus QC and (preliminary) ancestry outputs.

---

## Golden rules (read these first)

1. **One cohort at a time.** The pipeline locks its working directory while it runs.
   Start a cohort, wait for it to finish, *then* start the next one. Launching two runs
   from the same folder at once will corrupt both.
2. **Runs live on scratch; grab your results from `deliverables/`.** The full run
   (including big BAMs/temp files) is written to your scratch space, which the cluster
   **auto-deletes** after a while. When a run succeeds, the small final outputs are
   automatically copied to `deliverables/<cohort>/` in the pipeline folder — that's the
   copy you keep and share.
3. **Human, paired-end, GRCh38.** The pipeline aligns to GRCh38 / GENCODE v48. Your
   FASTQs must be paired-end human RNA-seq.
4. **Never submit a second cohort while one is still running** (see rule 1). Check first.

---

## One-time setup (do once)

1. **Get your own copy of the pipeline** so your runs don't clash with anyone else's:
   ```bash
   cd /projects/yateslab_genomics/sam        # or wherever the lab keeps your work
   git clone https://github.com/Kishaz/ngs-pipeline.git
   cd ngs-pipeline
   ```
2. **Load the environment** (do this every time you log in, too):
   ```bash
   module load anaconda3 && source activate snakemake
   ```
3. **Confirm you can read the shared references** (you do NOT build anything — the config
   already points at the lab's GRCh38 reference):
   ```bash
   ls /projects/yateslab_genomics/sam/ngs_pipeline/resources/GRCh38_v48/hisat2_index/ | head
   ```
   If that lists files, you're set.

---

## FASTQ naming the pipeline recognizes

The pipeline auto-detects samples and pairs R1/R2 from these patterns:

| R1 file looks like        | pairs with                |
|---------------------------|---------------------------|
| `SAMPLE_R1_001.fastq.gz`  | `SAMPLE_R2_001.fastq.gz`  |
| `SAMPLE_R1.fastq.gz`      | `SAMPLE_R2.fastq.gz`      |
| `SAMPLE_1.fastq.gz`       | `SAMPLE_2.fastq.gz`       |
| `SAMPLE.R1.fastq.gz`      | `SAMPLE.R2.fastq.gz`      |

The part before the R1 tag becomes the sample name. Put **one cohort's** FASTQs in **one
folder**.

---

## What the run does (raw FASTQ → results)

Each cohort goes through the full pipeline automatically:

```
FASTQs → fastp (trim/QC) → HISAT2 (align to GRCh38) → sort → mark duplicates
       → featureCounts (gene counts) → merge into counts matrix
       → ancestry estimation → QC & ancestry reports
```

**Timing:** alignment is the slow part — roughly ~30–60 min per sample, and this
single-job setup aligns ~2 samples at a time. So expect **hours to ~a day** for a
typical cohort (say 20–40 samples). The 3-day walltime in the script gives plenty of
margin. **If a cohort is 100+ samples, tell the PI** — that's better run by fanning the
jobs across the cluster (SLURM profile), which finishes far faster.

---

## Running one cohort (repeat this whole section per cohort)

### Step 1 — Make a submit script for the cohort

From your `ngs-pipeline` folder, create `run_cohort.sh` and **edit the two values at the
top** (`COHORT` and `FASTQ_DIR`):

```bash
cat > run_cohort.sh <<'EOF'
#!/bin/bash
#SBATCH --account=yateslab_genomics
#SBATCH --partition=cpu
#SBATCH --cpus-per-task=16
#SBATCH --mem=100G
#SBATCH --time=3-00:00:00
#SBATCH --job-name=rnaseq
#SBATCH --output=%x_%j.log

# ===================== EDIT THESE TWO =====================
COHORT="CRC"                                           # short cohort name (no spaces)
FASTQ_DIR="/projects/yates_lab_hpc/.../CRC/fastqs"     # folder of THIS cohort's FASTQs
# =========================================================

module load anaconda3
source activate snakemake
cd "$SLURM_SUBMIT_DIR"                                 # the pipeline folder you submit from

CFG="input_dir=$FASTQ_DIR results_dir=$COHORT results_base=/scratch/$USER samples=config/samples_${COHORT}.tsv"

# clear any stale lock, then run the whole pipeline for this cohort
snakemake --profile profiles/local --unlock --config $CFG
snakemake --profile profiles/local all --config $CFG \
    --rerun-incomplete --rerun-triggers mtime --cores 16 --resources mem_mb=90000
EOF
```

You only edit `COHORT` and `FASTQ_DIR`. Everything else (references, output location,
publishing) is handled automatically:
- results go to **`/scratch/$USER/<COHORT>/`**
- deliverables are auto-copied to **`deliverables/<COHORT>/`** in your pipeline folder.

### Step 2 — Make sure nothing else of yours is running

```bash
squeue -u $USER
ps -u $USER -o pid,cmd | grep "[s]nakemake"
```
If a previous cohort's job or a `snakemake` process is still there, **wait for it to
finish** before continuing.

### Step 3 — Preview before you commit (dry-run)

Load the same `CFG` values and do a dry-run — it discovers your samples and shows where
output will go, without running anything:
```bash
snakemake -n --config \
  input_dir="/projects/yates_lab_hpc/.../CRC/fastqs" \
  results_dir="CRC" results_base="/scratch/$USER" samples="config/samples_CRC.tsv" \
  2>&1 | grep -iE "Total samples|results_dir="
```
Check: it reports the right number of samples and `results_dir=/scratch/$USER/CRC`.

### Step 4 — Submit

```bash
sbatch run_cohort.sh
```

### Step 5 — Monitor

```bash
squeue -u $USER                    # your job; "PD" = pending (cluster busy), "R" = running
tail -f rnaseq_*.log               # live progress ("X of N steps done")
```
It's **done** when the log prints:
```
ALL STAGES COMPLETED SUCCESSFULLY
...
Deliverables:     published to .../ngs-pipeline/deliverables/CRC
```

### Step 6 — Collect your results

```bash
ls -R deliverables/CRC/
```
You'll find:
- `counts/gene_counts_matrix.tsv` — genes × samples (for DESeq2 / edgeR)
- `qc/`, `alignment_qc/` — QC report PDFs
- `ancestry/` — ancestry proportions + barplots

Copy `deliverables/<cohort>/` somewhere permanent — the full run on scratch will be
deleted later.

### Step 7 — Next cohort

Only now, go back to **Step 1** and repeat for the next cohort.

---

## If something goes wrong (quick fixes)

| Symptom | Fix |
|---|---|
| `Directory cannot be locked` | A previous run didn't exit cleanly. First confirm nothing is running (Step 2), then the wrapper's `--unlock` line clears it on the next submit. |
| Job stuck in `PD` (pending) | Cluster is busy — just wait; it starts when a slot frees. |
| `No space left on device` | Check `df -h /scratch/$USER`. Scratch should have room; if `/projects` is full that's a separate lab issue. |
| No samples discovered | Your FASTQ names don't match the table above, or the folder path is wrong. |
| A run failed partway | Fix the cause, then just re-submit the same `sbatch run_cohort.sh` — it resumes where it stopped. |

---

## Important notes

- **RNA-seq ancestry is preliminary/low-confidence.** RNA-seq covers very few of the
  ancestry-informative markers, so the ancestry proportions and any excluded samples
  should be treated as rough — do not report them as definitive. Ask the PI before using
  them.
- **Very large cohorts (100+ samples)** will take a long time as a single batch job. If a
  cohort is that big, tell the PI — it can be run faster by fanning the jobs across the
  cluster (SLURM profile), which the PI will set up.
- **Ask before changing `config/config.yaml`** — the per-cohort settings all go through
  the `run_cohort.sh` wrapper, so you shouldn't need to touch the config.
