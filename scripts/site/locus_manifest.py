"""List expected loci from the configured traits' post-processing reports."""

import pandas as pd

traits = pd.read_csv(snakemake.input.traits, sep="\t")
loci = []
for trait in traits.ID:
    frame = pd.read_csv(f"data/twas_out/{trait}.post.tsv", sep="\t")
    frame = frame.assign(chr_num=frame.CHR.str.removeprefix("chr").astype(int))
    frame = frame.sort_values(["chr_num", "P0"])
    frame = frame.assign(trait=trait, locus=range(1, len(frame) + 1))
    loci.append(frame[["trait", "locus", "FILE"]])

pd.concat(loci, ignore_index=True).to_csv(snakemake.output.tsv, sep="\t", index=False)
