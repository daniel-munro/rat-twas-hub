#!/bin/bash

set -e

wgt_dir=data/WEIGHTS
tissue=Brain
gtf=data/Rattus_norvegicus.Rnor_6.0.99.gtf

# Create table of all models
echo -e "WGT\tPANEL\tID\tCHR\tP0\tP1\tNSNPS\tHSQ\tHSQ.SE\tHSQ.PV\tTOP1.R2\tBLUP.R2\tENET.R2\tBSLMM.R2\tLASSO.R2\tTOP1.PV\tBLUP.PV\tENET.PV\tBSLMM.PV\tLASSO.PV" > data/all.models.par
tail -n+2 "$wgt_dir/$tissue/expression.pos" \
    | join -1 2 -2 1 -t'	' - <(tail -n+2 "$wgt_dir/$tissue/expression.profile") \
    | awk -v OFS='	' -v tis="$tissue" '{ print tis"/"$2,"RatGTEx."tis".expression",$1,$3,$4,$5,$6,$7,$8,$9,$10,$11,$12,$13,$14,$15,$16,$17,$18,$19,$20 }' \
    >> data/all.models.par
tail -n+2 "$wgt_dir/$tissue/stability.pos" \
    | join -1 2 -2 1 -t'	' - <(tail -n+2 "$wgt_dir/$tissue/stability.profile") \
    | awk -v OFS='	' -v tis="$tissue" '{ print tis"/"$2,"RatGTEx."tis".stability",$1,$3,$4,$5,$6,$7,$8,$9,$10,$11,$12,$13,$14,$15,$16,$17,$18,$19,$20 }' \
    >> data/all.models.par

## Run TWAS
mkdir -p data/tmp
cat data/traits.par | tail -n+2 | while read line; do
    trait=`echo $line | awk '{ print $2 }'`
    n=`echo $line | awk '{ print $3 }'`
    echo "Running TWAS on $trait"
    # Use either slurm jobs or GNU parallel
    # sbatch -n 1 -t 04:00:00 --mem-per-cpu=8G --array=1-20 --job-name="twas" scripts/TWAS.sh $trait $n
    parallel -j 4 --joblog data/tmp/$trait.log sh scripts/TWAS.sh $trait $n {} ::: {1..20}
done

## Add gene names since Ensembl IDs were used in TWAS
echo -e "ID\tNAME" > data/genes.par
perl -nle '/gene_id "([^"]+)".+gene_name "([^"]+)"/ && print "$1\t$2"' $gtf | sort | uniq >> data/genes.par

## Process TWAS results
bash scripts/MERGE.sh
mkdir jekyll
cp -r jekyll_base/* jekyll/
mkdir -p jekyll/traits jekyll/genes jekyll/data

N_TRAITS=`wc -l data/traits.par | awk '{ print $1 - 1 }'`
## sbatch -n 1 -t 08:00:00 --mem-per-cpu=8G --job-name="report" --array=1-$N_TRAITS REPORT.sh
parallel -j1 --joblog data/tmp/report.log Rscript scripts/REPORT_SINGLE.R {} ::: `seq 1 $N_TRAITS`

sbatch --wrap="Rscript scripts/REPORT_GENES.R"
Rscript scripts/REPORT_INDEX.R

## Compress TWAS data for site download links
cat data/trait_list.par | cut -f1 | while read line; do
    cd data/tmp && tar -cjf ../../jekyll/data/${line}.tar.bz2 ${line}/ ${line}.dat ${line}.dat.post.report && cd ..
done

mv -i data/traits data/genes data/data data/genes.md data/models.md data/traits.md data/genes.json jekyll/
