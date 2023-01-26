#!/bin/bash

#SBATCH -J snakemake_bill
#SBATCH -A recombinationlandscape
#SBATCH -o /shared/projects/recombinationlandscape/ARNAUD_STAGE/bill/LOGS/P1/snakemake_"%j".log
#SBATCH --partition=fast
#SBATCH --nodes=1
#SBATCH --mem=1G
#SBATCH --mail-user=arnaud.soulier@etu.umontpellier.fr
#SBATCH --mail-type=ALL

# Script variable
path="/shared/projects/recombinationlandscape/ARNAUD_STAGE/bill"

# Loading requiered module(s)
module load snakemake/6.14.0
module load minimap2/2.17
module load seqkit/2.1.0
module load samtools/1.9
module load sniffles/1.0.11
module load python/3.7
module load deeptools/3.5.0

# Mapping 
snakemake -s SCRIPTS/Snakefile --cores 1 --keep-going --cluster-config SCRIPTS/cluster_config.yaml --cluster "sbatch -J {cluster.job-name} -A recombinationlandscape -p {cluster.queue} {cluster.nodes} {cluster.cpus} {cluster.log} {cluster.mem}" --jobs 5
