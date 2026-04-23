# HPC Guide & Benchmarks (HPC-01)

This document serves as a technical reference for environment management, job scheduling, and resource allocation on HPC-01.

## Environment Management (Micromamba/Mamba/Conda)
### Update Environment
To synchronize the environment with changes in `environment.yml`:
```bash
# Using Micromamba (Recommended for speed)
micromamba env update -f environment.yml

# Using Mamba
mamba env update -f environment.yml
```

### Cleaning Hidden Cache
If Home space is low (`df -h ~`), run this to clear hidden installer files:
```bash
micromamba clean --all -y
# AND/OR
rm -rf ~/.conda/pkgs
```

## RNA-Seq Pipeline: HPC Operations & Performance Guide

## 1. Slurm Batch Job Optimization (Production)

To process the full Human dataset efficiently, use the provided `_hpc/hpc_run.sbatch` template.

```bash
#!/bin/bash
#SBATCH --job-name=RNAseq_Prod_Run
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=32
#SBATCH --mem=128G
#SBATCH --time=48:00:00             # 48 hours is extremely safe
```

## Resource Allocation (STAR Core Requirements)
STAR is one of the fastest aligners but is significantly **RAM-intensive**.

*   **RAM Intensity**: 
    *   **Production (Human Genome)**: Requires **~32GB RAM**. The index is loaded entirely into memory (Suffix Array).
    *   **Test Suite (Chr21)**: Requires **~1.5GB RAM**.
    *   **Risk**: If your machine has <32GB RAM, STAR will likely crash with a "segmentation fault" or "killed" message during index loading.
*   **CPU Intensity**: 
    *   STAR is highly parallelizable. Increasing `--runThreadN` significantly reduces mapping time but does not increase RAM usage significantly.
    *   Standard recommendation: **8 to 16 threads** for a fast production run.

## Useful HPC Diagnostic Commands
### 1. Check Personal Stats
*   **Check Disk Quota**: `quota -s` (Ensures you aren't hitting the /home limit)
*   **Check Job Status**: `squeue -u $USER`
*   **Cancel a Job**: `scancel <job_id>`

### 2. Check Node Resources
*   **Check Total Node Space**: `sinfo -o "%n %c %m"` (Hostname, CPUs, Memory in MB)
*   **Detailed Node View**: `scontrol show node <nodename>`
*   **Check Storage Mounts**: `df -h /mnt/22T` (Check available total storage)

## Remote Editing (Code Server / VS Code)
### 1. Start the Server (HPC Side)
Run this inside a `screen` session so it keeps running when you disconnect:
```bash
# Start/Enter a screen session
screen -S code-server

# Start the server (using a private port like 8888 or 9010)
code-server --bind-addr 127.0.0.1:8888

# Detach from screen: Press Ctrl+A, then D
```

### 2. Connect from Local Machine
On your **local laptop terminal**, run the SSH tunnel:
```bash
# Syntax: ssh -L [LocalPort]:localhost:[RemotePort] -p [SSH_Port] [User]@[Host]
ssh -L 8888:localhost:8888 -p 6820 phongdinh@trongchinh.zapto.org
```

### 3. Access in Browser
Go to: `http://localhost:8888`

> [!TIP]
> **Password**: If prompted, find your password on the HPC using:

```bash
grep password ~/.config/code-server/config.yaml
```
## Data Transfer (HPC to Local)
To download your results (e.g., QC reports or Counts) to your local machine:

### 1. Using rsync (Recommended)
`rsync` is the standard for bioinformatics because it is fast, compresses data, and can resume if interrupted. Run this on your **local machine**:

```bash
# Syntax: rsync -avzP -e 'ssh -p [Port]' [User]@[Host]:[RemotePath] [LocalPath]
rsync -avzP -e 'ssh -p 6820' phongdinh@trongchinh.zapto.org:/mnt/22T/phongdinh/RNAseq_pipeline/_data/qc ./
```

### 2. Common Flags
*   `-a`: Archive mode (preserves permissions/timestamps).
*   `-v`: Verbose (list files).
*   `-z`: Compress (saves bandwidth).
*   `-P`: Progress bars + partial transfer (resumable).

> [!TIP]
> **Downloading everything**: To download the entire project root (excluding raw data), you can use `--exclude`:
> ```bash
> rsync -avzP -e 'ssh -p 6820' --exclude '_data/fastq' phongdinh@trongchinh.zapto.org:/mnt/22T/phongdinh/RNAseq_pipeline/ ./
> ```

## Cleaning Data (Hard Reset)
To force a re-run of specific steps, delete the corresponding output directories:
*   **STAR Index**: `rm -rf _data/index/*` (Production) or `_data/index_test/*` (Test)
*   **BAM Alignment**: `rm -rf _data/bam/*` or `_data/bam_test/*`
*   **FeatureCounts**: `rm -rf _data/featurecounts/*` or `_data/featurecounts_test/*`

> [!IMPORTANT]
> **Do NOT delete** `_data/fastq_trimmed/` unless you want to re-run Trimmomatic (slow). FASTQ quality is independent of mapping indices.
