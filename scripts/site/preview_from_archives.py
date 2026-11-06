"""Rebuild the configured Jekyll preview from published TWAS archives."""

import csv
import json
import re
import shutil
import subprocess
import tarfile
import tempfile
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
SITE = ROOT / "jekyll"
SUMMARY_FIELDS = ["NUM.LOCI", "NUM.JOINT.GENES", "NUM.GENES", "AVG.CHISQ"]


def read_tsv(path):
    with path.open(newline="") as stream:
        return list(csv.DictReader(stream, delimiter="\t"))


def write_tsv(path, fields, rows):
    with path.open("w", newline="") as stream:
        writer = csv.DictWriter(stream, fields, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


def restore_trait(trait, workdir):
    output = workdir / "data/twas_out"
    archive = SITE / f"data/{trait}.tar.bz2"
    tables = {f"{trait}.all.tsv", f"{trait}.post.tsv"}
    locus_files = re.compile(rf"{re.escape(trait)}/{re.escape(trait)}\.chr\d+\.post\.loc_\d+\.(genes|cond\.csv)")
    with tarfile.open(archive, "r|bz2") as members:
        for member in members:
            if member.name not in tables and not locus_files.fullmatch(member.name):
                continue
            if not member.isfile():
                raise ValueError(f"Expected a regular file in {archive}: {member.name}")
            destination = output / member.name
            destination.parent.mkdir(parents=True, exist_ok=True)
            with members.extractfile(member) as source, destination.open("xb") as target:
                shutil.copyfileobj(source, target)

    loci = read_tsv(output / f"{trait}.post.tsv")
    for locus in loci:
        prefix = locus["FILE"].removeprefix("data/twas_out/")
        for suffix in (".genes", ".cond.csv"):
            name = prefix + suffix
            if not locus_files.fullmatch(name) or not (output / name).is_file():
                raise ValueError(f"Missing locus file in {archive}: {locus['FILE']}{suffix}")

    # Match the workflow summary, including rows with missing Z-scores in the mean.
    with (output / f"{trait}.all.tsv").open(newline="") as stream:
        reader = csv.reader(stream, delimiter="\t")
        z_col = next(reader).index("TWAS.Z")
        n_rows = 0
        sum_chisq = 0.0
        for row in reader:
            n_rows += 1
            if row[z_col] not in ("", "NA", "NaN"):
                sum_chisq += float(row[z_col]) ** 2

    return {
        "ID": trait,
        "NUM.LOCI": len(loci),
        "NUM.JOINT.GENES": sum(int(row["JOINT.GENES"]) for row in loci),
        "NUM.GENES": sum(int(row["HIT.GENES"]) for row in loci),
        "AVG.CHISQ": f"{sum_chisq / n_rows:.6g}" if n_rows else "NA",
    }


def stage_site(workdir, genes):
    site = workdir / "jekyll"

    def ignore(directory, names):
        if Path(directory) == SITE:
            return {"_site", "_site-preview", ".jekyll-cache", "data", "_data", "traits", "genes"}
        if Path(directory) == SITE / "assets/images":
            return {"porcupine"}
        return set()

    shutil.copytree(SITE, site, ignore=ignore)
    for directory in (
        "traits", "genes", "assets/images/porcupine", "_data/gene_models",
        "_data/trait_loci", "_data/trait_panels", "_data/trait_pleio", "_data/locus_models",
    ):
        (site / directory).mkdir(parents=True)
    cached = ["_data/panels.tsv", "_data/projects.tsv", "_data/stats.yml"]
    for gene in genes:
        cached.extend((f"genes/{gene}.md", f"_data/gene_models/{gene}.tsv"))
    for file in cached:
        shutil.copyfile(SITE / file, site / file)
    return site


def main():
    config = json.loads(subprocess.check_output(
        ["ruby", "-ryaml", "-rjson", "-e",
         'puts JSON.generate(YAML.safe_load_file(ARGV.fetch(0)).fetch("preview"))',
         str(SITE / "_config.preview.yml")],
        text=True,
    ))
    for key in ("traits", "genes"):
        ids = config[key]
        if not isinstance(ids, list) or not ids or not all(
            isinstance(item, str) and re.fullmatch(r"[\w-][\w.-]*", item) for item in ids
        ):
            raise ValueError(f"Preview {key} must be a nonempty list of IDs")
        config[key] = list(dict.fromkeys(ids))

    traits = config["traits"]
    metadata = {row["ID"]: row for row in read_tsv(SITE / "assets/downloads/traits.tsv")}
    for trait in traits:
        if trait not in metadata or not (SITE / f"data/{trait}.tar.bz2").is_file():
            raise ValueError(f"Preview trait requires download metadata and an archive: {trait}")

    build_root = ROOT / ".snakemake"
    build_root.mkdir(exist_ok=True)
    workdir = Path(tempfile.mkdtemp(prefix="archive-preview-", dir=build_root))
    print(f"Preview workspace: {workdir}", flush=True)
    site = stage_site(workdir, config["genes"])
    data = workdir / "data"
    (data / "twas_out").mkdir(parents=True)
    for name in ("all_models.tsv", "sequence_report.tsv"):
        (data / name).symlink_to(ROOT / "data" / name)
    shutil.copyfile(site / "_data/panels.tsv", data / "panels.tsv")

    summaries = [restore_trait(trait, workdir) for trait in traits]
    write_tsv(data / "traits.summary.tsv", ["ID", *SUMMARY_FIELDS], summaries)
    fields = ["ID", "N", "PROJECT", "TAGS", "NAME", "DESCRIPTION"]
    write_tsv(data / "traits.tsv", ["OUTPUT", *fields], [
        {"OUTPUT": f"data/twas_out/{trait}.all.tsv", **metadata[trait]} for trait in traits
    ])
    write_tsv(site / "_data/traits.tsv", [*fields, "NUM_LOCI", "NUM_JOINT_GENES", "NUM_GENES"], [
        {**metadata[row["ID"]], **{key.replace(".", "_"): row[key] for key in SUMMARY_FIELDS[:3]}}
        for row in summaries
    ])

    def render(script, *args):
        subprocess.run(["Rscript", str(ROOT / "scripts/site" / script), *args], cwd=workdir, check=True)

    render("build_pleiotropy.R")
    for row in summaries:
        trait = row["ID"]
        render("build_pages_for_trait.R", trait)
        (site / "traits" / trait).mkdir()
        (site / "_data/locus_models" / trait).mkdir()
        for locus in range(1, row["NUM.LOCI"] + 1):
            render("build_locus_page.R", trait, str(locus))
        render("porcupine_plots.R", trait)

    subprocess.run(
        ["bundle", "exec", "jekyll", "build", "--config", "_config.yml,_config.preview.yml"],
        cwd=site, check=True,
    )
    print(f"Built {site / '_site-preview/index.html'}")


if __name__ == "__main__":
    main()
