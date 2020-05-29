#!/bin/bash

#PBS -N RNA_Seq_PE
#PBS -l walltime=56:00:00
#PBS -q default
#PBS -o SNAKEMAKE.out
#PBS -e SNAKEMAKE.err
#PBS -l mem=10GB
#PBS -l nodes=1:ppn=1

#cd $PBS_O_WORKDIR

WORKDIR=$PBS_O_WORKDIR
cd  $WORKDIR

# Remove necessary to make sure the run info is correct
rm -f logs/*
mkdir -p logs

PIPELINE=$(grep -Po '"pipeline": *\K"[^"]*' conf_RNA_Seq.json | tr -d '"')

cp -r $PIPELINE/* ./
python scripts/move_fastq.py -json conf_RNA_Seq.json

snakemake --jobs 200 --latency-wait 100 --cluster-config $PIPELINE/scripts/cluster.json --cluster "qsub -l {cluster.walltime} -l {cluster.cores} -l {cluster.memory} -m n -q default -e $WORKDIR/logs/ -o $WORKDIR/logs/" --jobname 's.{rulename}.{jobid}' --stats $WORKDIR/logs/snakemake.stats >& $WORKDIR/logs/snakemake.log

rm -rf scripts