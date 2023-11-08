will_exclude_regions = False

if isinstance(config["dss_exclude_regions"], str):
    if os.path.exists(config["dss_exclude_regions"]):
        will_exclude_regions = True

        rule exclude_regions:
            input:
                bed=get_meth_bed,
                regions_to_exclude=config["dss_exclude_regions"],
            output:
                excluded_regions=temp(
                    "results/mCG/{tech}/filtered_input/{gtarget}_{sample}-filtered.bed.gz"
                ),
            threads: 1
            resources:
                mem=lambda wildcards, attempt: attempt * 4,
                hrs=72,
            envmodules:
                "modules",
                "modules-init",
                "modules-gs/prod",
                "modules-eichler/prod",
                "bedtools/2.29.2",
            shell:
                """
                bedtools subtract -a {input.bed} -b {input.regions_to_exclude} | gzip -c > {output.excluded_regions}
                """


rule dss_prepare_in_txt:
    input:
        bed=get_meth_bed
        if not will_exclude_regions
        else rules.exclude_regions.output.excluded_regions,
    output:
        bed=temp("results/mCG/{tech}/{gtarget}_{sample}.txt"),
    params:
        min_mod=config.get("dss_min_mod", 0),
        chromosomes=get_target_chromosomes,
    threads: 1
    resources:
        mem=lambda wildcards, attempt: attempt * 4,
        hrs=72,
    shell:
        """
        cmd=grep
        if grep -qE ".gz$" {input.bed}
        then
            cmd=zgrep
        fi

        if [ {wildcards.tech} == "hifi" ]; then
            # Columns grabbed are based on this documentation: https://github.com/PacificBiosciences/pb-CpG-tools#bed-file-format

            # chrom, start, n_valid, n_mod
            $cmd -wE "{params.chromosomes}" {input.bed} | awk '$7 >= {params.min_mod} {{print $1,$2,$6,$7}}' FS='\\t' OFS='\\t' > {output.bed}
        elif [ {wildcards.tech} == "ont" ]; then
            # Columns grabbed are based on this documentation: https://github.com/nanoporetech/modkit/#bedmethyl-column-descriptions

            # chrom, start, n_valid, n_mod
            $cmd -wE "{params.chromosomes}" {input.bed} | awk '$12 >= {params.min_mod} {{print $1,$2,$10,$12}}' FS='\\t' OFS='\\t' > {output.bed}
        else
            echo "Invalid tech wildcard: {wildcards.tech}" 1>&2; exit 1 
        fi
        """
