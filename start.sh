#!/bin/bash

echo $0 $@
sample=$1
reads1=$2
reads2=$3
bwa_cpu=$4
gatk_cpu=$5
output=$6
work_dir=$7

mkdir -p $work_dir
cd $work_dir
set -x
nextflow run /code/main-v2.nf -profile rcecc --sample=$sample --reads1=$reads1 --reads2=$reads2 --bwa_cpu=$bwa_cpu --gatk_cpu=$gatk_cpu --output $output

