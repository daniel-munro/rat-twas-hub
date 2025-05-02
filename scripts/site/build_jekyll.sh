#!/bin/bash

set -e

THREADS=$1

mkdir -p jekyll/traits jekyll/genes jekyll/data jekyll/_data

echo "Building trait pages..."
N_TRAITS=`wc -l data/traits.par | awk '{ print $1 - 1 }'`
parallel -j$THREADS --joblog data/twas_out/report.log Rscript scripts/site/build_pages_for_trait.R {} ::: `seq 1 $N_TRAITS`

echo "Building gene pages..."
Rscript scripts/site/build_gene_pages.R
Rscript scripts/site/build_cross_species.R

echo "Building index pages..."
Rscript scripts/site/build_index_pages.R
cp data/projects.tsv jekyll/_data/
cp data/panels.par jekyll/_data/panels.tsv

echo "Done building jekyll files."
