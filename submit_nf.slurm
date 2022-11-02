#! /usr/bin/env bash
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=4
#SBATCH --time=24:00:00
#SBATCH --job-name=NF_refactor
#SBATCH --output=R-%x.%J.out
#SBATCH --error=R-%x.%J.err
#SBATCH --mail-user=username@email.com
#SBATCH --mail-type=begin
#SBATCH --mail-type=end
# --account=project_name    # <= put HPC account name here, required on Atlas

# === For example, Load Modules here
# ==  Atlas HPC (will need a local install of nextflow)
# module load singularity

# == Nova HPC
# module load gcc/7.3.0-xegsmw4 nextflow
# module load singularity

# === Set working directory and in/out variables
cd ${SLURM_SUBMIT_DIR}

# === For example, using Conda Environments
# module load miniconda
# source activate /project/isu_gif_vrsc/jenchang/dot_home/.conda/envs/polishCLR_env

# === Main Program
# Pgos
nextflow run main.nf \
  --primary_assembly "/project/ag100pest/Pgos/RawData/3-unzip/all_p_ctg.fasta" \
  --illumina_reads "/project/ag100pest/Illumina_polishing/JAMU*{R1,R2}.fastq.bz2" \
  --pacbio_reads "/project/ag100pest/Pgos/RawData/m54334U_190823_194159.subreads.bam" \
  --k "21" \
  -resume
