#!/usr/bin/env Rscript
# salmon_tximport_tmm_normalize.R
# Usage:
# Rscript salmon_tximport_tmm_normalize.R --basedir ./ --tx2gene tx2gene.tsv --out GRN_Sc_TMM_normalized_CPM.txt

suppressPackageStartupMessages({
  # Helper: install if missing
  install_if_missing <- function(pkg, bioc=FALSE) {
    if (!requireNamespace(pkg, quietly=TRUE)) {
      message("Installing missing package: ", pkg)
      if (bioc) {
        if (!requireNamespace("BiocManager", quietly=TRUE)) {
          install.packages("BiocManager", repos="https://cloud.r-project.org")
        }
        BiocManager::install(pkg, ask=FALSE, update=FALSE)
      } else {
        install.packages(pkg, repos="https://cloud.r-project.org")
      }
    }
  }

  # Required packages
  install_if_missing("optparse", bioc=FALSE)
  install_if_missing("tximport", bioc=TRUE)
  install_if_missing("rtracklayer", bioc=TRUE)   # only needed if building tx2gene from GTF
  install_if_missing("edgeR", bioc=TRUE)

  library(optparse)
  library(tximport)
  library(rtracklayer)
  library(edgeR)
})

# ---- CLI options ----
option_list <- list(
  make_option(c("-b","--basedir"), type="character", default=".",
              help="Base directory to search for quant.sf (default = ./)"),
  make_option(c("-t","--tx2gene"), type="character", default=NULL,
              help="Path to tx2gene two-column TSV (tx_id <tab> gene_id). If omitted, provide --gtf"),
  make_option(c("-g","--gtf"), type="character", default=NULL,
              help="Path to GTF/GFF annotation to build tx2gene if --tx2gene not provided"),
  make_option(c("-o","--out"), type="character", default="normalized_tmm_cpm.txt",
              help="Output filename (samples x genes) [default: normalized_tmm_cpm.txt]"),
  make_option(c("--pattern"), type="character", default="quant.sf",
              help="Filename to look for inside sample folders (default: quant.sf)")
)

opt <- parse_args(OptionParser(option_list = option_list))

basedir <- opt$basedir
tx2gene_path <- opt$tx2gene
gtf_path <- opt$gtf
outfile <- opt$out
pattern <- opt$pattern

# ---- find quant.sf files ----
sf_files <- list.files(path = basedir, pattern = pattern, recursive = TRUE, full.names = TRUE)
if (length(sf_files) == 0) stop("No quant.sf files found under basedir. Check --basedir and --pattern.")

# derive sample names from parent folder, clean common Salmon prefixes/suffixes
sample_dirs <- basename(dirname(sf_files))
sample_names <- gsub("^Salmon_", "", sample_dirs)
sample_names <- gsub("_output$", "", sample_names)
if (any(duplicated(sample_names))) {
  warning("Duplicate sample names detected after cleaning; using parent directory names instead.")
  sample_names <- basename(dirname(sf_files))
}
names(sf_files) <- sample_names

message(length(sf_files), " quant.sf files found. Samples:")
print(names(sf_files))

# ---- build tx2gene ----
if (!is.null(tx2gene_path)) {
  if (!file.exists(tx2gene_path)) stop("Provided --tx2gene file not found: ", tx2gene_path)
  message("Reading tx2gene from: ", tx2gene_path)
  tx2gene <- read.table(tx2gene_path, header = FALSE, sep = "\t", stringsAsFactors = FALSE,
                        col.names = c("tx_id","gene_id"))
} else if (!is.null(gtf_path)) {
  if (!file.exists(gtf_path)) stop("Provided --gtf file not found: ", gtf_path)
  message("Building tx2gene from GTF: ", gtf_path)
  gtf <- import(gtf_path)
  mcols_gtf <- mcols(gtf)
  if (!all(c("transcript_id","gene_id") %in% colnames(mcols_gtf))) {
    stop("GTF missing transcript_id or gene_id attributes.")
  }
  tx2gene_df <- as.data.frame(mcols_gtf[, c("transcript_id","gene_id")], stringsAsFactors = FALSE)
  colnames(tx2gene_df) <- c("tx_id","gene_id")
  tx2gene <- unique(tx2gene_df[complete.cases(tx2gene_df), ])
  message("Built tx2gene with ", nrow(tx2gene), " mappings.")
} else {
  stop("Either --tx2gene or --gtf must be provided.")
}

# ---- run tximport ----
message("Running tximport (countsFromAbundance = 'lengthScaledTPM') ...")
txi <- tximport(files = sf_files, type = "salmon", tx2gene = tx2gene,
                countsFromAbundance = "lengthScaledTPM", ignoreTxVersion = TRUE)

message("Imported. Genes:", nrow(txi$counts), " Samples:", ncol(txi$counts))

# ---- edgeR TMM normalization ----
# Create DGEList from txi$counts (gene x sample)
# edgeR expects genes x samples matrix
gene_counts <- txi$counts
if (!is.matrix(gene_counts)) stop("txi$counts is not a matrix")

y <- DGEList(counts = gene_counts)
# calcNormFactors computes TMM factors and stores them in y$samples$norm.factors
y <- calcNormFactors(y, method = "TMM")

message("Library sizes (raw):")
print(y$samples$lib.size)
message("TMM normalization factors:")
print(y$samples$norm.factors)

# Compute normalized CPM (counts per million) using normalized library sizes
# cpm(..., normalized.lib.sizes=TRUE) uses lib.size * norm.factors internally
cpm_norm <- cpm(y, normalized.lib.sizes = TRUE, log = FALSE, prior.count = 0.25) 
# cpm_norm is genes x samples

# ---- prepare output: transpose to samples x genes, write file ----
norm_by_sample <- t(cpm_norm)
colnames(norm_by_sample) <- rownames(cpm_norm)   # gene names
rownames(norm_by_sample) <- colnames(cpm_norm)  # sample names

out_df <- data.frame(samples = rownames(norm_by_sample), norm_by_sample, check.names = FALSE, stringsAsFactors = FALSE)

write.table(out_df, file = outfile, sep = "\t", quote = FALSE, row.names = FALSE)
message("Wrote TMM-normalized CPM matrix to: ", outfile)
message("NOTE: Output units are CPM normalized by TMM (edgeR calcNormFactors).")
