!/bin/bash -l
#SBATCH --cluster=genius
#SBATCH --job-name fastqc
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --time=12:00:00
#SBATCH -o fastqc.%j.out
#SBATCH -A lp_svbelleghem

module load FastQC/0.11.8-Java-1.8.0_162

fastqc POG_HiFi_reads.fastq.gz
