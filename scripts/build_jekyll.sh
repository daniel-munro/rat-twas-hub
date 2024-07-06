#!/bin/bash

set -e

# Need -p because Snakemake will create the output file's directory already
mkdir -p jekyll
cp -r jekyll_base/* jekyll/
mkdir -p jekyll/traits jekyll/genes jekyll/data

echo "Building trait pages..."
N_TRAITS=`wc -l data/traits.par | awk '{ print $1 - 1 }'`
parallel -j1 --joblog data/twas_out/report.log Rscript scripts/build_pages_for_trait.R {} ::: `seq 1 $N_TRAITS`

echo "Building gene pages..."
Rscript scripts/build_gene_pages.R

echo "Building index pages..."
Rscript scripts/build_index_pages.R

echo "Compressing data for download..."
## Compress TWAS data for site download links
cat data/traits.par | tail -n+2 | cut -f2 | while read trait; do
    cd data/twas_out && tar -cjf ../../jekyll/data/${trait}.tar.bz2 ${trait}/ ${trait}.dat ${trait}.dat.post.report && cd ../..
done

echo "Done building jekyll files."
