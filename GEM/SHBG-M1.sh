#!/bin/bash
#SBATCH --partition=batch
#SBATCH --job-name=GEM-SHBGxVeg
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=8
#SBATCH --time=144:00:00
#SBATCH --mem=30000
#SBATCH --output=%x.%j.out
#SBATCH --error=%x.%j.err
#SBATCH --array=1-22

i=$SLURM_ARRAY_TASK_ID

cd /work/kylab/mike/TestxVeg/GEM

#ml PLINK/2.00-alpha2.3-x86_64-20210920-dev
ml GEM/1.4.1-foss-2019b


genoindir=("/scratch/mf91122/LipidsxVeg/genotypeQC/GEM2")

phenodir=("/scratch/mf91122/TestxVeg/pheno")

outdir=("/scratch/mf91122/TestxVeg/GEM/SHBG")

#phenotypes=("LDL" "HDL" "Tot_Chol" "TAGs")
#phenotypes="logSHBG"
phenotypes="SHBG"

#exposures=("Consistent_Self_Reported_Vegetarian_across_all_24hr" "Self_Reported_Vegetarian_plus_strict_initial_and24")
#exposures=("Consistent_Self_Reported_Vegetarian_across_all_24hr")
exposures=("Veg2" "Veg1")

for j in ${phenotypes[@]} 
        do

for e in ${exposures[@]} 
        do

mkdir -p $outdir/$j/$e

echo running "$j" and "$e"

GEM \
--bgen $genoindir/chr"$i".bgen \
--sample $genoindir/chr"$i".sample \
--pheno-file $phenodir/SHBGxVeg_pheno.M1.csv \
--sampleid-name IID \
--pheno-name $j \
--covar-names Age Sex Geno_batch \
PC1 PC2 PC3 PC4 PC5 PC6 PC7 PC8 PC9 PC10 \
center1 center2 center3 center4 center5 \
center6 center7 center8 center9 center10 \
center11 center12 center13 center14 center15 \
center16 center17 center18 center19 center20 \
--robust 0 \
--exposure-names "$e" \
--threads 8 \
--out $outdir/$j/$e/chr"$i"

done
done
