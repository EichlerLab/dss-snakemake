import pandas as pd
import sys
import itertools
import os


# --------- Config --------- #
configfile: "config/config.yaml"


# --------- Globals --------- #
dss_cntr = config["dss_cntr"]
gtarget = config.get("gtarget", "autosome")
dss_mode = config.get("mode", "group")
tech = config.get("tech", "ont")

# --------- Load files --------- #
df = pd.read_table(config["manifest"], dtype=str, header=0).set_index(
    "sample_name", drop=False
)  # sample_name\tbed_file\tgroup_name


# --------- Constraints --------- #
wildcard_constraints:
    mode="group_comparison|model",
    gtarget="autosome|sex|wgs",
    sampleA="|".join(df["sample_name"]),
    sampleB="|".join(df["sample_name"]),


# --------- Input functions --------- #
def get_final_output(wildcards):
    fp = "results/{tech}/dss/group_comparison/{group_name}/{gtarget}/{sampleA}_vs_{sampleB}_DMR-prcnt_summary.tsv.gz"
    if config.get("annotations", []):
        if config["annotations"].items():
            fp = "results/{tech}/dss/group_comparison/{group_name}/{gtarget}/annotations/{sampleA}_vs_{sampleB}_DMR-annotated.tsv.gz"

    targets = []
    for g in df.group.unique():
        pairs, _ = get_group_pairs(group_name=g)
        for a, b in pairs:
            targets.append(
                fp.format(
                    tech=tech, group_name=g, gtarget=gtarget, sampleA=a, sampleB=b
                )
            )
    return targets


def get_target_chromosomes(wildcards):
    autosomes = ["chr{}".format(x) for x in list(range(1, 23))]
    sex_chr = ["X", "Y"]
    if wildcards.gtarget == "wgs":
        return "|".join(autosomes + sex_chr)
    elif wildcards.gtarget == "sex":
        return "|".join(sex_chr)
    else:
        return "|".join(autosomes)


def get_meth_bed(wildcards):
    return df.at[wildcards.sample, "bed"]


def get_dss_params(wildcards):
    param_dict = {
        "groupA": wildcards.sampleA,
        "groupB": wildcards.sampleB,
        "sample_names": [wildcards.sampleA, wildcards.sampleB],
    }
    return param_dict


def get_group_pairs(group_name) -> (list[tuple], pd.DataFrame):
    # Get the target sample for the group.
    target_sample = config["group_target"].get(group_name, "")

    grouped_df = df.query(rf"group == '{group_name}'")
    pairs = list(itertools.combinations(grouped_df.sample_name.values, 2))
    if target_sample:
        pairs = [x for x in pairs if target_sample in x]

    return pairs, grouped_df


def get_anno_dict(wildcards):
    return config["annotations"]


def get_dss_summaries(wildcards):
    # Target pattern
    fp = "results/{tech}/dss/group_comparison/{group_name}/{gtarget}/temp/{sample}-{sampleA}_vs_{sampleB}_DMR.tsv"
    targets = []

    pairs, grouped_df = get_group_pairs(group_name=wildcards.group_name)

    for entry in grouped_df.itertuples():
        for a, b in pairs:
            targets.append(
                fp.format(
                    tech=wildcards.tech,
                    group_name=entry.group,
                    gtarget=wildcards.gtarget,
                    sample=entry.sample_name,
                    sampleA=a,
                    sampleB=b,
                )
            )

    return targets
