#!/bin/bash

## Create table of all models
wgt_dir=/gpfs/home/dmunro/bulk/pantry/Brain_rn6_Pheast/intermediate/twas
echo -e "WGT\tPANEL\tID\tCHR\tP0\tP1\tNSNPS\tHSQ\tHSQ.SE\tHSQ.PV\tTOP1.R2\tBLUP.R2\tENET.R2\tBSLMM.R2\tLASSO.R2\tTOP1.PV\tBLUP.PV\tENET.PV\tBSLMM.PV\tLASSO.PV" > all.models.par
tail -n+2 "$wgt_dir/expression.pos" \
    | join -1 2 -2 1 -t'	' - <(tail -n+2 "$wgt_dir/expression.profile") \
    | awk -v OFS='	' '{ print "RatGTEx.Brain."$2,"RatGTEx.Brain.expression",$1,$3,$4,$5,$6,$7,$8,$9,$10,$11,$12,$13,$14,$15,$16,$17,$18,$19,$20 }' \
    >> all.models.par
tail -n+2 "$wgt_dir/stability.pos" \
    | join -1 2 -2 1 -t'	' - <(tail -n+2 "$wgt_dir/stability.profile") \
    | awk -v OFS='	' '{ print "RatGTEx.Brain."$2,"RatGTEx.Brain.stability",$1,$3,$4,$5,$6,$7,$8,$9,$10,$11,$12,$13,$14,$15,$16,$17,$18,$19,$20 }' \
    >> all.models.par

## Run TWAS
mkdir -p genes traits data tmp
tail -n+2 traits.par | cut -f2,3 > trait_list.par
bash run_TWAS.sh

## Add gene names since Ensembl IDs were used in TWAS
echo -e "ID\tNAME" > genes.par
perl -nle '/gene_id "([^"]+)".+gene_name "([^"]+)"/ && print "$1\t$2"' /gpfs/home/dmunro/ratgtex/ref_rn6/Rattus_norvegicus.Rnor_6.0.99.gtf | sort | uniq >> genes.par

## Process TWAS results
bash MERGE.sh
bash run_REPORT.sh
sbatch --wrap="Rscript REPORT_GENES.R"
Rscript REPORT_INDEX.R

## Compress TWAS data for site download links
cat trait_list.par | cut -f1 | while read line; do
    cd tmp && tar -cjf ../data/${line}.tar.bz2 ${line}/ ${line}.dat ${line}.dat.post.report && cd ..
done

## Move into web directory
mv -i traits genes data genes.md models.md traits.md genes.json jekyll/
