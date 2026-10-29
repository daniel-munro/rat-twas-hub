from functools import lru_cache, partial

import pandas as pd

tissues = ["Adipose", "BLA", "Brain", "Eye", "IC", "IL", "LHb", "Liver", "NAcc", "OFC", "PL", "pVTA", "RMTg"]
modalities = ["alt_polyA", "alt_TSS", "expression", "isoforms", "splicing", "stability"]
gtf = "data/GCF_036323735.1_GRCr8_genomic.chr.gtf"
wgt_dir = "data/WEIGHTS"
ldref_prefix = "data/LDREF/Brain_v4."

traits_df = pd.read_csv("data/traits.tsv", sep="\t", index_col=1)
traits = traits_df.index.tolist()


@lru_cache(maxsize=1)
def read_site_loci(manifest_path):
    # Called only after the checkpoint finishes; its manifest is fixed for this run.
    return pd.read_csv(manifest_path, sep="\t")


def trait_loci(trait):
    manifest_path = checkpoints.site_locus_manifest.get().output.tsv
    loci = read_site_loci(str(manifest_path))
    return loci.loc[loci.trait == trait]


def locus_inputs(wildcards):
    loci = trait_loci(wildcards.trait)
    prefix = loci.iloc[int(wildcards.locus) - 1]["FILE"]
    return [f"{prefix}.genes", f"{prefix}.cond.csv"]


def trait_locus_inputs(wildcards):
    return [f"{prefix}{suffix}" for prefix in trait_loci(wildcards.trait)["FILE"]
            for suffix in [".genes", ".cond.csv"]]


def site_locus_outputs(trait, wildcards):
    outputs = []
    for locus in range(1, len(trait_loci(trait)) + 1):
        outputs.extend([
            f"jekyll/traits/{trait}/{locus}.md",
            f"jekyll/traits/{trait}/{locus}.cond.csv",
            f"jekyll/_data/locus_models/{trait}/{locus}.tsv",
        ])
    return outputs


localrules:
    all_models_file,
    sumstats,
    twas_top,
    combine_chroms_assoc,
    combine_chroms_top,
    combine_chroms_post,
    traits_info,
    genes_n_assoc_info,
    genes_n_models_info,
    jekyll_copy_projects,
    site_locus_manifest,

wildcard_constraints:
    trait="[^/.]+"

rule all:
    input:
        expand("data/twas_out/{trait}.all.tsv", trait=traits),
        expand("data/twas_out/{trait}.top.tsv", trait=traits),
        expand("data/twas_out/{trait}.post.tsv", trait=traits),
        "data/traits.summary.tsv",
        "data/genes_n_models.tsv",
        "data/genes_n_assoc.tsv",
        "jekyll/genes.json",
        "jekyll/_data/panels.tsv",
        "jekyll/_data/projects.tsv",
        "jekyll/_data/stats.yml",
        "jekyll/_data/traits.tsv",
        expand("jekyll/data/{trait}.tar.bz2", trait=traits),
        expand("jekyll/traits/{trait}.md", trait=traits),
        expand("jekyll/assets/images/porcupine/{trait}.png", trait=traits),
        "jekyll/cross-species.json",
        "jekyll/_data/gene_pages.done",
        [partial(site_locus_outputs, trait) for trait in traits],

rule twas_results:
    """TWAS tables only; optionally select one project with --config project=ID."""
    input:
        expand(
            "data/twas_out/{trait}.{table}.tsv",
            trait=(traits_df.index[traits_df["PROJECT"] == config["project"]].tolist()
                   if "project" in config else traits),
            table=["all", "top", "post"],
        )

rule all_models_file:
    """Create table of all models"""
    input:
        pos = expand(wgt_dir + "/{tissue}/{modality}.pos", tissue=tissues, modality=modalities),
        profile = expand(wgt_dir + "/{tissue}/{modality}.profile", tissue=tissues, modality=modalities)
    output:
        "data/all_models.tsv"
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
                    | awk -v OFS='	' -v tis="$tissue" -v mod="$modality" '{{ print tis"/"$2,"RatGTEx."tis"."mod,$1,$3,$4,$5,$7,$8,$9,$10,$11,$12,$13,$14,$15,$16,$17,$18,$19,$20 }}' \
                    >> {output}
            done
        done
        """

rule sumstats:
    """Convert sumstats files"""
    input:
        gwas = lambda w: "data/gwas_original/sumstats/{project}/regressedlr_{trait}_chrgwas{chrom}.mlma.gz".format(
            project=traits_df.loc[w.trait, "PROJECT"],
            trait=w.trait,
            chrom=w.chrom
        )
    output:
        sumstats = "data/gwas/{trait}/{trait}.{chrom}.sumstats"
    shell:
        """
        mkdir -p data/gwas/{wildcards.trait}
        zcat {input.gwas} | awk 'BEGIN {{OFS="\t"}} \
            NR==1 {{print $2, $4, $5, "Z"}} \
            NR>1 && $8 != "inf" {{print "chr"$2, $4, $5, $7/$8}}' \
            > {output.sumstats}
        """

rule twas:
    """Run TWAS"""
    input:
        models = "data/all_models.tsv",
        sumstats = "data/gwas/{trait}/{trait}.{chrom}.sumstats",
        panels = "data/panels.tsv",
        ldref = multiext(f"{ldref_prefix}chr{{chrom}}", ".bed", ".bim", ".fam"),
    output:
        assoc = "data/twas_out/{trait}/{trait}.chr{chrom}.all.tsv"
    params:
        wgt_dir = wgt_dir,
        ldref_prefix = ldref_prefix,
        trait_dir = "data/twas_out/{trait}",
        gwas_n = lambda w: traits_df.loc[w.trait, "N"],
    group: "twas"
    resources:
        mem_mb = 16000,
        runtime = "16h",
    shell:
        """
        mkdir -p {params.trait_dir}
        Rscript scripts/twas/FUSION.assoc_test.R \
            --sumstats {input.sumstats} \
            --out {output.assoc} \
            --weights {input.models} \
            --weights_dir {params.wgt_dir} \
            --ref_ld_chr {params.ldref_prefix} \
            --chr chr{wildcards.chrom} \
            --coloc_P 0.0005 \
            --GWASN {params.gwas_n} \
            --PANELN {input.panels}
        """

rule twas_top:
    """Get significant TWAS hits on one chromosome for one trait"""
    input:
        models = "data/all_models.tsv",
        assoc = "data/twas_out/{trait}/{trait}.chr{chrom}.all.tsv",
    output:
        top = "data/twas_out/{trait}/{trait}.chr{chrom}.top.tsv"
    group: "twas"
    shell:
        """
        TOTAL=`wc -l {input.models} | awk '{{ print $1 - 1 }}'`
        cat {input.assoc} | awk -vt=$TOTAL 'NR == 1 || $20 < 0.05/t' > {output.top}
        """

rule twas_post:
    """Run post-processing on TWAS hits on one chromosome for one trait"""
    input:
        top = "data/twas_out/{trait}/{trait}.chr{chrom}.top.tsv",
        sumstats = "data/gwas/{trait}/{trait}.{chrom}.sumstats",
        ldref = multiext(f"{ldref_prefix}chr{{chrom}}", ".bed", ".bim", ".fam"),
    output:
        post = "data/twas_out/{trait}/{trait}.chr{chrom}.post.report",
    params:
        ldref_prefix = ldref_prefix,
        out_prefix = "data/twas_out/{trait}/{trait}.chr{chrom}.post",
        max_r2 = 0.85, # default is 0.9, but that was giving errors in solve() in conditional analysis
    group: "twas"
    resources:
        mem_mb = 16000,
    shell:
        """
        if [ $(wc -l < {input.top}) -eq 1 ]; then
            printf "FILE\tCHR\tP0\tP1\tHIT.GENES\tJOINT.GENES\tBEST.TWAS.P\tBEST.SNP.P\tCOND.SNP.P\tVAR.EXP\n" > {output.post}
        else
            Rscript scripts/twas/FUSION.post_process.R \
                --sumstats {input.sumstats} \
                --input {input.top} \
                --out {params.out_prefix} \
                --minp_input 1 \
                --ref_ld_chr {params.ldref_prefix} \
                --chr chr{wildcards.chrom} \
                --locus_win 200e3 \
                --max_r2 {params.max_r2} \
                --report
        fi
        """

rule combine_chroms_assoc:
    """Concatenate per-chromosome TWAS results for one trait"""
    input:
        assocs = expand("data/twas_out/{{trait}}/{{trait}}.chr{chrom}.all.tsv", chrom=range(1, 21))
    output:
        assoc = "data/twas_out/{trait}.all.tsv"
    shell:
        """
        cat {input.assocs} | awk 'NR == 1 || $1 != "PANEL"' > {output.assoc}
        """

rule combine_chroms_top:
    """Concatenate per-chromosome significant hits for one trait"""
    input:
        tops = expand("data/twas_out/{{trait}}/{{trait}}.chr{chrom}.top.tsv", chrom=range(1, 21))
    output:
        top = "data/twas_out/{trait}.top.tsv"
    shell:
        """
        cat {input.tops} | awk 'NR == 1 || $1 != "PANEL"' > {output.top}
        """

rule combine_chroms_post:
    """Concatenate per-chromosome TWAS post-processing reports for one trait"""
    input:
        posts = expand("data/twas_out/{{trait}}/{{trait}}.chr{chrom}.post.report", chrom=range(1, 21))
    output:
        post = "data/twas_out/{trait}.post.tsv"
    shell:
        """
        cat {input.posts} | awk 'NR == 1 || $1 != "FILE"' > {output.post}
        """

rule traits_info:
    """Generate traits info file"""
    input:
        traits = "data/traits.tsv",
        assocs = expand("data/twas_out/{trait}.all.tsv", trait=traits),
        posts = expand("data/twas_out/{trait}.post.tsv", trait=traits),
    output:
        tsv = "data/traits.summary.tsv"
    shell:
        """
        tail -n+2 {input.traits} | awk '{{ print $2 }}' | while read id; do
            avgchisq=`tail -n+2 data/twas_out/$id.all.tsv | awk '{{ print $19^2 }}' | awk 'BEGIN {{ n=0; s=0; }} {{ n++; s += $1; }} END {{ if ( n > 0 ) print s / n; else print "NA"; }}'`
            cat data/twas_out/$id.post.tsv | awk -v chi=$avgchisq -v id=$id 'BEGIN {{ loc=0; tothit=0; tot=0; }} $1 != "FILE" {{ tothit += $5; tot+=$6; loc++; }} END {{ print id,loc,tot,tothit,chi }}'
        done | awk 'BEGIN {{ print "ID NUM.LOCI NUM.JOINT.GENES NUM.GENES AVG.CHISQ" }} {{ print $0 }}' | tr ' ' '\t' > {output.tsv}
        """

rule genes_n_models_info:
    """Count total models per gene"""
    input:
        models = "data/all_models.tsv"
    output:
        tsv = "data/genes_n_models.tsv"
    shell:
        """
        tail -n+2 {input.models} | cut -f3 | sed -E 's/__.*$//' | sort | uniq -c | awk 'BEGIN {{OFS="\\t"}} {{ print $2,$1 }}' > {output.tsv}
        """

rule genes_n_assoc_info:
    """Count significant TWAS associations per gene"""
    input:
        traits = "data/traits.tsv",
        tops = expand("data/twas_out/{trait}.top.tsv", trait=traits),
    output:
        tsv = "data/genes_n_assoc.tsv"
    shell:
        """
        tail -n+2 {input.traits} | awk '{{ print $2 }}' | while read id; do
            tail -n+2 data/twas_out/$id.top.tsv | cut -f3 | sed -E 's/__.*$//' | uniq | sort | uniq
        done | sort | uniq -c | awk 'BEGIN {{OFS="\\t"}} {{ print $2,$1 }}' > {output.tsv}
        """

checkpoint site_locus_manifest:
    """Discover site loci once all traits' post-processing reports are ready."""
    input:
        traits = "data/traits.tsv",
        posts = expand("data/twas_out/{trait}.post.tsv", trait=traits),
    output:
        tsv = "data/site_loci.tsv",
    script:
        "scripts/site/locus_manifest.py"


rule jekyll_trait_page:
    """Build the Jekyll page for a single trait"""
    input:
        traits = "data/traits.tsv",
        panels = "data/panels.tsv",
        trait_info = "data/traits.summary.tsv",
        assoc = "data/twas_out/{trait}.all.tsv",
        post = "data/twas_out/{trait}.post.tsv",
        models = "data/all_models.tsv",
        locus_files = trait_locus_inputs,
        renderer = "scripts/site/build_pages_for_trait.R",
    output:
        page = "jekyll/traits/{trait}.md",
        loci = "jekyll/_data/trait_loci/{trait}.tsv",
        pleio = "jekyll/_data/trait_pleio/{trait}.tsv",
        panels = "jekyll/_data/trait_panels/{trait}.tsv",
    group: "site_trait_{trait}"
    shell:
        "Rscript {input.renderer:q} {wildcards.trait:q}"

rule jekyll_locus_page:
    """Render one locus; grouped with its trait page into one submission."""
    input:
        page = "jekyll/traits/{trait}.md",
        traits = "data/traits.tsv",
        panels = "data/panels.tsv",
        post = "data/twas_out/{trait}.post.tsv",
        locus_files = locus_inputs,
        renderer = "scripts/site/build_locus_page.R",
    output:
        page = "jekyll/traits/{trait}/{locus}.md",
        cond = "jekyll/traits/{trait}/{locus}.cond.csv",
        models = "jekyll/_data/locus_models/{trait}/{locus}.tsv",
    wildcard_constraints:
        locus = "[1-9][0-9]*",
    group: "site_trait_{trait}"
    resources:
        mem_mb = 512,
    shell:
        "Rscript {input.renderer:q} {wildcards.trait:q} {wildcards.locus:q}"


rule porcupine_plots:
    """Generate one trait plot; batch independent jobs with --group-components."""
    input:
        chroms = "data/sequence_report.tsv",
        all_models = "data/all_models.tsv",
        assoc = "data/twas_out/{trait}.all.tsv",
        post = "data/twas_out/{trait}.post.tsv",
    output:
        png = "jekyll/assets/images/porcupine/{trait}.png",
    group: "site_plots"
    resources:
        mem_mb = 2000,
    shell:
        "Rscript scripts/site/porcupine_plots.R {wildcards.trait:q}"

rule jekyll_gene_pages:
    """Build Jekyll pages for all genes
    
    Outputs for all genes:
    - jekyll/genes/{gene}.md
    - jekyll/_data/gene_models/{gene}.tsv
    """
    input:
        panels = "data/panels.tsv",
        all_models = "data/all_models.tsv",
        traits = "data/traits.tsv",
        traits_info = "data/traits.summary.tsv",
        gtf = gtf,
    output:
        marker = "jekyll/_data/gene_pages.done",
    resources:
        mem_mb = 32000,
    shell:
        """
        Rscript scripts/site/build_gene_pages.R {input.gtf}
        date > {output.marker}
        """

rule jekyll_summary_data:
    """Build summary data files for Jekyll (_data)."""
    input:
        traits = "data/traits.tsv",
        traits_info = "data/traits.summary.tsv",
        all_models = "data/all_models.tsv",
        genes_models = "data/genes_n_models.tsv",
        genes_assoc = "data/genes_n_assoc.tsv",
        panels = "data/panels.tsv",
        gtf = gtf,
    output:
        traits_tsv = "jekyll/_data/traits.tsv",
        genes_json = "jekyll/genes.json",
        stats_yml = "jekyll/_data/stats.yml",
        panels_tsv = "jekyll/_data/panels.tsv",
    shell:
        "Rscript scripts/site/build_summary_data.R {input.gtf}"

rule cross_species_orthologs:
    """Fetch and save the Ensembl 108 mapping independently of association counts."""
    input:
        human_models = ancient("data/cross_species/genes.models.nfo"),
        rat_models = ancient("data/all_models.tsv"),
    output:
        mapping = "data/cross_species/human_rat_orthologs.ensembl108.tsv",
    shell:
        "Rscript scripts/site/fetch_orthologs.R"

rule jekyll_cross_species:
    """Build cross-species JSON for the site"""
    input:
        mapping = "data/cross_species/human_rat_orthologs.ensembl108.tsv",
        human_models = "data/cross_species/genes.models.nfo",
        rat_models = "data/genes_n_models.tsv",
        human_assoc = "data/cross_species/genes.nfo",
        rat_assoc = "data/genes_n_assoc.tsv",
    output:
        json = "jekyll/cross-species.json",
    shell:
        "Rscript scripts/site/build_cross_species.R"

rule jekyll_copy_projects:
    """Copy projects.tsv into Jekyll _data directory"""
    input:
        tsv = "data/projects.tsv"
    output:
        tsv = "jekyll/_data/projects.tsv"
    shell:
        "cp {input.tsv} {output.tsv}"

rule compress_twas_data:
    """Compress one trait's results; batch jobs with --group-components."""
    input:
        assoc = "data/twas_out/{trait}.all.tsv",
        post = "data/twas_out/{trait}.post.tsv",
    output:
        tar = "jekyll/data/{trait}.tar.bz2",
    group: "site_archives"
    resources:
        mem_mb = 512,
    shell:
        """
        tar -cjf {output.tar:q} -C data/twas_out \
            {wildcards.trait:q}/ {wildcards.trait:q}.all.tsv {wildcards.trait:q}.post.tsv
        """
