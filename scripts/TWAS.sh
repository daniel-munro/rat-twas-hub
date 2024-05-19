#!/bin/sh

# Run TWAS analyses
# $1 = trait name
# $2 = integer GWAS sample size
# Run assumes FUSION scripts are in scripts/ and WEIGHTS, GWAS, and LDREF are in data/.

# Inputs:
# data/all.models.par
# data/GWAS/$TRAIT/$TRAIT.$CHR.sumstats
# data/LDREF/Brain.$CHR.{bed,bim,fam}
# data/panels.par

# Outputs:
# data/twas_out/$TRAIT/$TRAIT.$CHR.dat
# data/twas_out/$TRAIT/$TRAIT.$CHR.top
# data/twas_out/$TRAIT/$TRAIT.$CHR.post.*

TRAIT=$1
N=$2
# --- Get chromosome either from slurm array task ID, if run that way, otherwise from first argument
if [ -z "$SLURM_ARRAY_TASK_ID" ]; then
    CHR=$3
else
    CHR=$SLURM_ARRAY_TASK_ID
fi
# ---
# trait=`basename $GWAS | sed 's/\.[[:alnum:]]\+\.sumstats//'`

POS="data/all.models.par"
GWAS="data/GWAS/$TRAIT/$TRAIT.$CHR.sumstats"
LDREF="data/LDREF/Brain."
OUT="data/twas_out/$TRAIT/$TRAIT.$CHR"
TOTAL=`wc -l $POS | awk '{ print $1 - 1 }'`

mkdir -p "data/twas_out/$TRAIT"
Rscript scripts/FUSION.assoc_test.R \
    --sumstats "$GWAS" \
    --out "$OUT.dat" \
    --weights "$POS" \
    --weights_dir data/WEIGHTS \
    --ref_ld_chr "$LDREF" \
    --chr $CHR --coloc_P 0.00005 --GWASN $N --PANELN data/panels.par
# cat $OUT.dat | awk -vt=$TOTAL 'NR == 1 || $20 < 0.05/t' > $OUT.top
cat $OUT.dat | awk -vt=$TOTAL 'NR == 1 || $20 < 0.05' > $OUT.top # For testing

if [ $(wc -l < "$OUT.top") -eq 1 ]; then
    printf "FILE\tCHR\tP0\tP1\tHIT.GENES\tJOINT.GENES\tBEST.TWAS.P\tBEST.SNP.P\tCOND.SNP.P\tVAR.EXP\n" > "$OUT.post.report"
else
    Rscript scripts/FUSION.post_process.R \
        --sumstats "$GWAS" \
        --input "$OUT.top" \
        --out "$OUT.post" \
        --minp_input 1 \
        --ref_ld_chr "$LDREF" \
        --chr "$CHR" --locus_win 200e3 --report
        # --chr "$CHR" --plot --locus_win 200e3 --report
fi
