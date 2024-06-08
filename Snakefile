import pandas as pd

wgt_dir = "data/WEIGHTS"
tissues = ["Brain"]
gtf = "data/Rattus_norvegicus.Rnor_6.0.99.gtf"

modalities = ["alt_polyA", "alt_TSS", "expression", "isoforms", "splicing", "stability"]

traits_df = pd.read_csv("data/traits.par", sep="\t", index_col=1)
traits = traits_df.index.tolist()

localrules:
    gene_names,
    all_models_file,
    twas_top,

rule all:
    input:
        expand("data/twas_out/{trait}/{trait}.{chrom}.post.report", trait=traits, chrom=range(1, 21)),

rule gene_names:
    """Add gene names since Ensembl IDs were used in TWAS"""
    input:
        gtf = gtf
    output:
        tsv = "data/gene_names.tsv"
    shell:
        """
        echo -e "ID\tNAME" > {output.tsv}
        perl -nle '/gene_id "([^"]+)".+gene_name "([^"]+)"/ && print "$1\t$2"' {input.gtf} | sort | uniq >> {output.tsv}
        """

rule all_models_file:
    """Create table of all models"""
    input:
        pos = expand(wgt_dir + "/{tissue}/{modality}.pos", tissue=tissues, modality=modalities),
        profile = expand(wgt_dir + "/{tissue}/{modality}.profile", tissue=tissues, modality=modalities)
    output:
        "data/all_models.par"
    params:
        wgt_dir = wgt_dir,
        tissues = " ".join(tissues),
        modalities = " ".join(modalities)
    shell:
        """
        echo -e "WGT\tPANEL\tID\tCHR\tP0\tP1\tNSNPS\tHSQ\tHSQ.SE\tHSQ.PV\tTOP1.R2\tBLUP.R2\tENET.R2\tBSLMM.R2\tLASSO.R2\tTOP1.PV\tBLUP.PV\tENET.PV\tBSLMM.PV\tLASSO.PV" > {output}
        for tissue in {params.tissues}; do
            for modality in {params.modalities}; do
                tail -n+2 "{params.wgt_dir}/$tissue/$modality.pos" \
                    | join -1 2 -2 1 -t'	' - <(tail -n+2 "{params.wgt_dir}/$tissue/$modality.profile") \
                    | awk -v OFS='	' -v tis="$tissue" -v mod="$modality" '{{ print tis"/"$2,"RatGTEx."tis"."mod,$1,$3,$4,$5,$6,$7,$8,$9,$10,$11,$12,$13,$14,$15,$16,$17,$18,$19 }}' \
                    >> {output}
            done
        done
        """

rule twas:
    """Run TWAS"""
    input:
        models = "data/all_models.par",
        sumstats = "data/GWAS/{trait}/{trait}.{chrom}.sumstats",
        panels = "data/panels.par",
        ldref = multiext("data/LDREF/Brain.{chrom}", ".bed", ".bim", ".fam"),
    output:
        dat = "data/twas_out/{trait}/{trait}.{chrom}.dat"
    params:
        wgt_dir = wgt_dir,
        ldref_prefix = "data/LDREF/Brain.",
        trait_dir = "data/twas_out/{trait}",
        gwas_n = lambda w: traits_df.loc[w.trait, "N"],
    group: "twas"
    resources:
        mem_mb = 8000,
        walltime = 4,
    shell:
        """
        mkdir -p {params.trait_dir}
        Rscript scripts/FUSION.assoc_test.R \
            --sumstats {input.sumstats} \
            --out {output.dat} \
            --weights {input.models} \
            --weights_dir {params.wgt_dir} \
            --ref_ld_chr {params.ldref_prefix} \
            --chr {wildcards.chrom} \
            --coloc_P 0.0005 \
            --GWASN {params.gwas_n} \
            --PANELN {input.panels}
        """
        
rule twas_top:
    input:
        models = "data/all_models.par",
        dat = "data/twas_out/{trait}/{trait}.{chrom}.dat",
    output:
        top = "data/twas_out/{trait}/{trait}.{chrom}.top"
    group: "twas"
    shell:
        """
        TOTAL=`wc -l {input.models} | awk '{{ print $1 - 1 }}'`
        cat {input.dat} | awk -vt=$TOTAL 'NR == 1 || $20 < 0.05/t' > {output.top}
        """

rule twas_post:
    input:
        top = "data/twas_out/{trait}/{trait}.{chrom}.top",
        sumstats = "data/GWAS/{trait}/{trait}.{chrom}.sumstats",
        ldref = multiext("data/LDREF/Brain.{chrom}", ".bed", ".bim", ".fam"),
    output:
        report = "data/twas_out/{trait}/{trait}.{chrom}.post.report",
    params:
        ldref_prefix = "data/LDREF/Brain.",
        out_prefix = "data/twas_out/{trait}/{trait}.{chrom}.post",
    group: "twas"
    resources:
        mem_mb = 8000,
    shell:
        """
        if [ $(wc -l < {input.top}) -eq 1 ]; then
            printf "FILE\tCHR\tP0\tP1\tHIT.GENES\tJOINT.GENES\tBEST.TWAS.P\tBEST.SNP.P\tCOND.SNP.P\tVAR.EXP\n" > {output.report}
        else
            Rscript scripts/FUSION.post_process.R \
                --sumstats {input.sumstats} \
                --input {input.top} \
                --out {params.out_prefix} \
                --minp_input 1 \
                --ref_ld_chr {params.ldref_prefix} \
                --chr {wildcards.chrom} \
                --locus_win 200e3 \
                --report
        fi
        """

# ## Process TWAS results
# bash scripts/merge_TWAS_results.sh
# mkdir jekyll
# cp -r jekyll_base/* jekyll/
# mkdir -p jekyll/traits jekyll/genes jekyll/data

# N_TRAITS=`wc -l data/traits.par | awk '{ print $1 - 1 }'`
# ## sbatch -n 1 -t 08:00:00 --mem-per-cpu=8G --job-name="report" --array=1-$N_TRAITS REPORT.sh
# ## R --slave --args ${SLURM_ARRAY_TASK_ID} < build_pages_for_trait.R
# ## for i in `seq 1 $N_TRAITS`; do
# ##     Rscript build_pages_for_trait.R $i
# ## done
# parallel -j1 --joblog data/twas_out/report.log Rscript scripts/build_pages_for_trait.R {} ::: `seq 1 $N_TRAITS`

# ## sbatch --wrap="Rscript scripts/build_gene_pages.R"
# Rscript scripts/build_gene_pages.R
# Rscript scripts/build_index_pages.R

# ## Compress TWAS data for site download links
# cat data/traits.par | tail -n+2 | cut -f2 | while read trait; do
#     cd data/twas_out && tar -cjf ../../jekyll/data/${trait}.tar.bz2 ${trait}/ ${trait}.dat ${trait}.dat.post.report && cd ../..
# done
