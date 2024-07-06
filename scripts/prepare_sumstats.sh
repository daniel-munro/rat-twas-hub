set -e

# Convert to sumstats format
cat data/GWAS_original/traits_pruned.r2_50.tsv | head -n10 | while read line; do
    echo $line
    project=`echo $line | awk '{ print $1 }'`
    trait=`echo $line | awk '{ print $2 }'`
    mkdir -p data/GWAS/$trait
    for chrom in {1..20}; do
        awk 'BEGIN {OFS="\t"} NR==1 {print $2, $4, $5, "Z"} NR>1 && $8 != "inf" {print "chr"$2, $4, $5, $7/$8}' data/GWAS_original/sumstats/$project/regressedlr_${trait}_chrgwas${chrom}.mlma > data/GWAS/$trait/$trait.$chrom.sumstats
    done
    # Get sample sizes from GCTA logs
     grep -P '^\d+ observations' data/GWAS_original/sumstats/$project/regressedlr_${trait}_chrgwas1.log | cut -f1 -d' ' | sed "s/^/$project\t$trait\t/" >> data/GWAS_original/sample_sizes.tsv
done
