#!/bin/bash
# Sample batchscript to run a parallel python job on HPC using 10 CPU cores
 
#SBATCH --partition=bch-compute                # queue to be used
#SBATCH --time=72:00:00             # Running time (in hours-minutes-seconds)
#SBATCH --job-name=cellranger-arc            # Job name
#SBATCH --mail-type=BEGIN,END,FAIL       # send and email when the job begins, ends or fails
#SBATCH --output=cellranger-count_%j.txt         # Name of the output file
#SBATCH --nodes=1               # Number of compute nodes
#SBATCH --ntasks=16               # Number of cpu cores on one node
#SBATCH --mem=64G

source /programs/biogrids.shrc

MRO_DISK_SPACE_CHECK=disable cellranger-arc count --id=Mult_ABNCR6D_corrected \
             --reference=/home/ch228212/AK/Multiome-strand-02212022/MedGenome_12122022/refdata-cellranger-arc-GRCh38-2020-A-2.0.0 \
             --libraries=libraries.Mult_CR6D_corrected.csv \
             --localcores=16 \
             --localmem=64