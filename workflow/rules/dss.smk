if dss_mode == "group":

    rule dss:
        input:
            file_names=lambda wildcards: expand(
                "results/mCG/{tech}/{gtarget}_{sample}.txt",
                sample=[wildcards.sampleA, wildcards.sampleB],
                gtarget=gtarget,
                tech=tech,
            ),
        output:
            dmr="results/{tech}/dss/group_comparison/{group_name}/{gtarget}/{sampleA}_vs_{sampleB}_DMR.tsv",
            dml="results/{tech}/dss/group_comparison/{group_name}/{gtarget}/{sampleA}_vs_{sampleB}_DML.tsv",
        params:
            get_dss_params,
            output_prefix="results/{tech}/dss/group_comparison/{group_name}/{gtarget}/{sampleA}_vs_{sampleB}",
        threads: 8
        resources:
            mem=lambda wildcards, attempt: attempt * 8,
            hrs=72,
        log:
            "results/{tech}/dss/group_comparison/{group_name}/{gtarget}/log/{sampleA}_vs_{sampleB}.log",
        container:
            dss_cntr
        script:
            "../scripts/DSS.R"

else:
    sys.exit(f"Unsupported mode ATM: {dss_mode}. Only group comparisons are supported.")


rule dss_summary_table:
    input:
        dss_out="results/{tech}/dss/group_comparison/{group_name}/{gtarget}/{sampleA}_vs_{sampleB}_DMR.tsv",
        sample_bed="results/mCG/{tech}/{gtarget}_{sample}.txt",
    output:
        sample_summary_table=temp(
            "results/{tech}/dss/group_comparison/{group_name}/{gtarget}/temp/{sample}-{sampleA}_vs_{sampleB}_DMR.tsv"
        ),
    threads: 1
    resources:
        mem=lambda wildcards, attempt: attempt * 16,
        hrs=72,
    run:
        df = pd.read_table(
            input.sample_bed,
            header=None,
            names=["chrom", "start", "n_valid", "n_mod"],
        )
        dss_df = pd.read_table(input.dss_out, header=0)

        # # Calculate the percentage methylated
        def get_percnt_methylated(row):
            n_valid_sum = df.loc[
                df.start.between(row.start, row.end), "n_valid"
            ].sum()
            n_mod_sum = df.loc[df.start.between(row.start, row.end), "n_mod"].sum()

            percnt_methylated = n_mod_sum / n_valid_sum

            return float("{:.2f}".format(percnt_methylated))

        dss_df[wildcards.sample] = dss_df.apply(
            lambda row: get_percnt_methylated(row=row), axis=1
        )
        dss_df.fillna(0.0, inplace=True)

        dss_df.to_csv(
            output.sample_summary_table, sep="\t", header=True, index=False
        )


rule merge_dss_summary_table:
    input:
        all_group_samples=get_dss_summaries,
    output:
        summary_table="results/{tech}/dss/group_comparison/{group_name}/{gtarget}/{sampleA}_vs_{sampleB}_DMR-prcnt_summary.tsv.gz",
    threads: 1
    resources:
        mem=lambda wildcards, attempt: attempt * 16,
        hrs=72,
    run:
        # Get the DSS columns for group comparison from one file
        target_columns = [
            "chr",
            "start",
            "end",
            "length",
            "nCG",
            "meanMethy1",
            "meanMethy2",
            "diff.Methy",
            "areaStat",
        ]
        dss_df = pd.read_table(
            input.all_group_samples[0], header=0, usecols=target_columns
        )

        # Get the sample columns
        samples_df = pd.concat(
            [
                pd.read_table(x, header=0, usecols=lambda c: c not in target_columns)
                for x in input.all_group_samples
            ],
            axis=1,
        )

        pd.concat([dss_df, samples_df], axis=1).to_csv(
            output.summary_table, sep="\t", header=True, index=False
        )


rule add_annotation:
    input:
        unpack(get_anno_dict),
        merged_summary="results/{tech}/dss/group_comparison/{group_name}/{gtarget}/{sampleA}_vs_{sampleB}_DMR-prcnt_summary.tsv.gz",
    output:
        added_annotation="results/{tech}/dss/group_comparison/{group_name}/{gtarget}/annotations/{sampleA}_vs_{sampleB}_DMR-annotated.tsv.gz",
    threads: 1
    resources:
        mem=lambda wildcards, attempt: attempt * 16,
        hrs=72,
    run:
        from pybedtools import BedTool

        df = pd.read_table(input.merged_summary, header=0)
        df["id"] = df.apply(
            lambda row: f"{row.chr}-{row.start}-{row.end}-{row.length}", axis=1
        )
        mod_df = df.loc[:, ["chr", "start", "end", "id"]].sort_values(
            by=["chr", "start", "end"]
        )

        starting_bed = BedTool.from_dataframe(mod_df)

        keys = [k for k in input.keys() if k != "merged_summary"]

        for key in keys:
            anno_path = input.get(key)
            anno_df = pd.read_table(
                anno_path, header=None, names=["chr", "pos", "end", key]
            ).sort_values(by=["chr", "pos", "end"])
            anno_bed = BedTool.from_dataframe(anno_df)

            if "gene" in key.lower():
                closest_genes_df = (
                    starting_bed.closest(anno_bed, D="ref")
                    .groupby(g=[4, 9], c=8, o="collapse")
                    .to_dataframe(
                        disable_auto_names=True,
                        names=["id", "distance(bp)", key],
                        header=None,
                    )
                )

                for idx, row in closest_genes_df.iterrows():
                    if row[key] == ".":
                        closest_genes_df.loc[idx, key] = "N/A"
                        closest_genes_df.loc[idx, "distance(bp)"] = "N/A"
                    else:
                        closest_genes_df.loc[idx, key] = ",".join(
                            pd.Series(row[key].split(",")).drop_duplicates().to_list()
                        )

                closest_genes_df = closest_genes_df[["id", key, "distance(bp)"]]
                closest_genes_df.rename(
                    columns={"distance(bp)": f"distance(bp)_{key}"}, inplace=True
                )

                df = df.merge(closest_genes_df, on="id", how="left")

                del closest_genes_df

            else:
                hits = starting_bed.intersect(anno_bed, wa=True, wb=True)

                try:
                    len(hits[0])

                    annotated_df = hits.groupby(g=[4], c=8, o="collapse").to_dataframe(
                        disable_auto_names=True, names=["id", key], header=None
                    )
                    annotated_df[key] = annotated_df.apply(
                        lambda x: ",".join(
                            pd.Series(x[key].split(",")).drop_duplicates().to_list()
                        ),
                        axis=1,
                    )
                    annotated_df = annotated_df[["id", key]]

                    df = df.merge(annotated_df, on="id", how="left")
                    del annotated_df
                except IndexError:
                    df[key] = "N/A"

                del hits

        df.drop(columns=["id"], inplace=True)
        df.fillna("N/A", inplace=True)

        df.to_csv(output.added_annotation, sep="\t", header=True, index=False)
