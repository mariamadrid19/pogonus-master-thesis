#!/bin/bash -l
#SBATCH --cluster=genius
#SBATCH --job-name isocon
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=10
#SBATCH --time=48:00:00
#SBATCH -o isocon.%j.out
#SBATCH -A lp_svbelleghem

conda activate IsoCon

IsoCon pipeline -fl_reads /scratch/leuven/357/vsc35707/pogonus/reads/Pogonus_PACBIO_RNA/POG_IsoSeq_HiFi_demux.fastq -outfolder /scratch/leuven/357/vsc35707/pogonus/reads/Pogonus_PACBIO_RNA/isocon_assembly
~                                             
