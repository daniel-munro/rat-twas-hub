#!/bin/sh

# Run TWAS analyses
# $1 = trait name
# $2 = integer GWAS sample size
# Run assumes all FUSION code, LDREF, and weights are in a FUSION directory
# --- for parallel jobs, your job scheduler variable goes here
CHR=${SLURM_ARRAY_TASK_ID}
# ---

TRAIT=$1
N=$2
# trait=`basename $GWAS | sed 's/\.[[:alnum:]]\+\.sumstats//'`

POS="all.models.par"
GWAS="./FUSION/GWAS/$TRAIT/$TRAIT.$CHR.sumstats"
LDREF="./FUSION/LDREF/Brain."
OUT="tmp/$TRAIT/$TRAIT.$CHR"
TOTAL=`wc -l $POS | awk '{ print $1 - 1 }'`

mkdir -p "tmp/$TRAIT"
Rscript ./FUSION/FUSION.assoc_test.R \
    --sumstats $GWAS \
    --out $OUT.dat \
    --weights $POS \
    --weights_dir ./FUSION/WEIGHTS \
    --ref_ld_chr $LDREF \
    --chr $CHR --coloc_P 0.00005 --GWASN $N --PANELN panels.par
cat $OUT.dat | awk -vt=$TOTAL 'NR == 1 || $20 < 0.05/t' > $OUT.top

Rscript ./FUSION/FUSION.post_process.R \
    --sumstats $GWAS \
    --input $OUT.top \
    --out $OUT.post \
    --minp_input 1 \
    --ref_ld_chr $LDREF \
    --chr $CHR --plot --locus_win 200e3 --report
