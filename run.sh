#!/bin/bash

set -e

wgt_dir=data/WEIGHTS
tissue=Brain
gtf=data/Rattus_norvegicus.Rnor_6.0.99.gtf

## Add gene names since Ensembl IDs were used in TWAS
echo -e "ID\tNAME" > data/gene_names.tsv
perl -nle '/gene_id "([^"]+)".+gene_name "([^"]+)"/ && print "$1\t$2"' $gtf | sort | uniq >> data/gene_names.tsv

# Create table of all models
echo -e "WGT\tPANEL\tID\tCHR\tP0\tP1\tNSNPS\tHSQ\tHSQ.SE\tHSQ.PV\tTOP1.R2\tBLUP.R2\tENET.R2\tBSLMM.R2\tLASSO.R2\tTOP1.PV\tBLUP.PV\tENET.PV\tBSLMM.PV\tLASSO.PV" > data/all_models.par
for modality in alt_polyA alt_TSS expression isoforms splicing stability; do
    tail -n+2 "$wgt_dir/$tissue/$modality.pos" \
        | join -1 2 -2 1 -t'	' - <(tail -n+2 "$wgt_dir/$tissue/$modality.profile") \
        | awk -v OFS='	' -v tis="$tissue" -v mod="$modality" '{ print tis"/"$2,"RatGTEx."tis"."mod,$1,$3,$4,$5,$6,$7,$8,$9,$10,$11,$12,$13,$14,$15,$16,$17,$18,$19 }' \
        >> data/all_models.par
done

## Run TWAS
mkdir -p data/twas_out
cat data/traits.par | tail -n+2 | while read line; do
    trait=`echo $line | awk '{ print $2 }'`
    n=`echo $line | awk '{ print $3 }'`
    echo "Running TWAS on $trait"
    # Use either slurm jobs or GNU parallel
    # sbatch -n 1 -t 04:00:00 --mem-per-cpu=8G --array=1-20 --job-name="twas" scripts/TWAS.sh $trait $n
    parallel -j 4 --joblog data/twas_out/$trait.log sh scripts/TWAS.sh $trait $n {} ::: {1..20}
done

## Process TWAS results
bash scripts/merge_TWAS_results.sh
mkdir jekyll
cp -r jekyll_base/* jekyll/
mkdir -p jekyll/traits jekyll/genes jekyll/data

N_TRAITS=`wc -l data/traits.par | awk '{ print $1 - 1 }'`
## sbatch -n 1 -t 08:00:00 --mem-per-cpu=8G --job-name="report" --array=1-$N_TRAITS REPORT.sh
## R --slave --args ${SLURM_ARRAY_TASK_ID} < build_pages_for_trait.R
## for i in `seq 1 $N_TRAITS`; do
##     Rscript build_pages_for_trait.R $i
## done
parallel -j1 --joblog data/twas_out/report.log Rscript scripts/build_pages_for_trait.R {} ::: `seq 1 $N_TRAITS`

## sbatch --wrap="Rscript scripts/build_gene_pages.R"
Rscript scripts/build_gene_pages.R
Rscript scripts/build_index_pages.R

## Compress TWAS data for site download links
cat data/traits.par | tail -n+2 | cut -f2 | while read trait; do
    cd data/twas_out && tar -cjf ../../jekyll/data/${trait}.tar.bz2 ${trait}/ ${trait}.dat ${trait}.dat.post.report && cd ../..
done
