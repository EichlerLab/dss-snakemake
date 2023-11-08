#!/usr/bin/env Rscript
# Author: Mei Wu, https://github.com/projectoriented

####### -------------- Libraries -------------- #######
library("DSS")

####### -------------- Functions -------------- #######
process_files <- function(file_names, sample_names) {
    # Create an empty list to store the data frames
    data_list <- list()

    # Iterate over the file names
    for (i in seq_along(file_names)) {
        # Read the data from each file
        data <- read.table(file_names[i], header = FALSE)
        colnames(data) <-c('chr','pos', 'N', 'X')

        # Store the data frame in the list
        data_list[[i]] <- data
    }

    BSobj <- makeBSseqData(data_list, sample_names)

    return(BSobj)
}

by_group_comparison <- function(bs_object, group_1_names, group_2_names, n_threads, output_prefix) {
    dmlTest <- DMLtest(bs_object, group1=group_1_names, group2=group_2_names, ncores=n_threads, smoothing=TRUE, smoothing.span=500)
    dmls <- callDML(dmlTest, delta=0.1, p.threshold=1e-5)
    dmrs <- callDMR(dmlTest, delta=0.1, p.threshold=1e-5)
    write.table(dmrs, file=paste(output_prefix, "DMR.tsv", sep="_"),sep="\t", quote=FALSE, row.names=FALSE)
    write.table(dmls, file=paste(output_prefix, "DML.tsv", sep="_"),sep="\t", quote=FALSE, row.names=FALSE)
}


####### -------------- Analysis -------------- #######
log <- file(snakemake@log[[1]], open = "wt")
sink(log, type = "message")

file.names <- snakemake@input[["file_names"]]

print(file.names)

output.prefix <- snakemake@params[["output_prefix"]]
n.threads <- snakemake@threads

param_dict <- snakemake@params[[1]]

print(param_dict)

BSobj <- process_files(file.names, param_dict$sample_names)

by_group_comparison(bs_object = BSobj, group_1_names = param_dict$groupA, group_2_names = param_dict$groupB, n_threads = n.threads, output_prefix=output.prefix)
