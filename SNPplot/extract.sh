#!/bin/bash
#SBATCH --partition=highmem_p
#SBATCH --job-name=chr378
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=16
#SBATCH --time=144:00:00
#SBATCH --mem=180000
#SBATCH --output=%x.%j.out
#SBATCH --error=%x.%j.err

cd /work/kylab/mike/TestxVeg/GEM/chr3-78Mb

ml PLINK/2.00-alpha2.3-x86_64-20210920-dev

genoindir=("/scratch/mf91122/LipidsxVeg/genotypeQC/GEM2pgen/pfile")
outdir=("/scratch/mf91122/TestxVeg/GEM/prelim-test3SHBG/logSHBG/Veg2/exportA")

mkdir -p $outdir

plink2 \
--pfile $genoindir/chr3 \
--extract extractchr3-78.txt \
--export A \
--out $outdir/chr3
