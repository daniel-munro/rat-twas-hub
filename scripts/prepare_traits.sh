set -e

# Get sample sizes from GCTA logs
rm -f data/gwas_original/sample_sizes.tsv
cat data/gwas_original/traits.tsv | tail -n+2 | while read line; do
    project=`echo $line | awk '{ print $1 }'`
    trait=`echo $line | awk '{ print $2 }'`
    cat data/gwas_original/sumstats/$project/regressedlr_${trait}_chrgwas1.log \
        | grep -P '^\d+ observations' \
        | cut -f1 -d' ' \
        | sed "s/^/$project\t$trait\t/" \
        >> data/gwas_original/sample_sizes.tsv
done

# Add sample sizes to other trait info to make traits.par
Rscript scripts/prepare_trait_table.R
