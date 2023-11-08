# dss-snakemake
A Bioinformatics pipeline for differential methylation analysis using DSS. *Pipeline is WIP!*

## Getting started
1. Clone the repo
2. Install Snakemake version >= 7.6.0
3. Fill in the config and manifests.
4. Start the analysis!
    * Begin with a dry-run
    ```
    ./runlocal 10 -np
    ```
    * If dry-run looks good, proceed with:
    ```
    ./runlocal 10
    ```

## Pre-reqs
```shell
config/
└── config.yaml
```
## Config explanation
```shell
mode: group # group or regression model- only group is supported at the moment.
dss_cntr: "/path/to/singularity.sif"
gtarget: autosome # autosome, sex, or wgs
tech: ont # hifi or ont

manifest: config/manifest.tab # tab delimited columns: sample_name, bed, group

dss_exclude_regions: "/path/to/lc.bed.gz" # the pipeline only parses the three columns, bed3
dss_min_mod: 0 # optional- if you want to exclude any sites with zero methylation

# optional annotations but I do parse the keyword genes special as in I get the distance from the gene as well as intersection
annotations:
    cpg_island: ""
    genes: ""

# optional but please make sure all keys are strings
# each key (e.g. one) corresponds to the name in the manifest under group and purpose explained in next section 
group_target:
    "one": "sample_child"
```
The purpose group target is to target the sample within a group and exclude any pairs that does not contain target sample. For example:
```shell
$ cat config/manifest.tab
sample_name	bed	group
sample_father	sample_father_cpg-pileup.bed.gz	one
sample_mother	sample_mother_cpg-pileup.bed.gz	one
sample_child	sample_child_cpg-pileup.bed.gz	one
```
gives us, the following pairs: 
```shell
results/ont/dss/group_comparison/one/autosome/
├── sample_father_vs_sample_child_DMR-prcnt_summary.tsv.gz
└── sample_mother_vs_sample_child_DMR-prcnt_summary.tsv.gz
```
and **NOT**
```shell
├── sample_father_vs_sample_mother_DMR-prcnt_summary.tsv.gz
```
## FAQ
1. What do the results mean?
   1. `sample_father_vs_sample_child_DMR-prcnt_summary.tsv.gz` tells us the differential methylated regions when we compare father and child. In other words, relative to the father- these are the regions that are differential. Please refer to the [DSS](https://bioconductor.org/packages/devel/bioc/vignettes/DSS/inst/doc/DSS.html) manual for output explanation.
2. Must I have the keyword genes?
   1. If you're providing a bed file of genes, please use the key "genes" to get both intersection and distance to that intersection.