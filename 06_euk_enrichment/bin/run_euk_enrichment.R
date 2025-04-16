#!/usr/bin/env Rscript

# Load libraries
suppressPackageStartupMessages({
  library(optparse)
  library(data.table)
  library(taxonomizr)
})

# Define arguments
option_list <- list(
  make_option(c("--aln"), type="character", help="Foldseek alignment .m8 or .m8.gz"),
  make_option(c("--map"), type="character", help="Uniprot50 UID-to-UniRef map"),
  make_option(c("--tax"), type="character", help="Uniprot50 taxonomy table"),
  make_option(c("--db"), type="character", help="Path to accessionTaxa.sql"),
  make_option(c("--eval_cut"), type="numeric", default=1e-3, help="E-value cutoff [default: %default]"),
  make_option(c("--out"), type="character", default="euk_enrichment_results.tsv", help="Output summary file")
)

opt <- parse_args(OptionParser(option_list=option_list))

# Read input with or without compression
read_auto <- function(path) {
  if (grepl("\\.gz$", path)) {
    fread(cmd = paste("pigz -cqd", path))
  } else {
    fread(path)
  }
}

# Load data
cat("Reading data...\n")
aln <- read_auto(opt$aln)[evalue <= opt$eval_cut]
aln[, UID := sub("-.*", "", sub("AF-", "", target))]
map <- fread(opt$map, col.names = c("UniRef50", "UID"))
tax <- fread(opt$tax)

cat("Assigning taxonomy...\n")
tax[, domain := getTaxonomy(V3, opt$db, desiredTaxa = c("superkingdom"))]
setnames(tax, "V1", "UniRef50")
map_tax <- merge(map, tax[, .(UniRef50, domain)], by = "UniRef50")

cat("Merging taxonomy...\n")
aln <- merge(aln, map_tax, by = "UID", all.x = TRUE)

# Global domain counts
domain_dist <- tax[, .N, by = domain]
all_euk <- domain_dist[domain == "Eukaryota"]$N
all_other <- sum(domain_dist[domain != "Eukaryota"]$N)

# Fisher test function
euk.fisher.test <- function(euk, other, all_euk, all_other){
  out <- NULL
  for (i in 1:length(euk)){
    test.mt <- matrix(c(other[i], all_other, euk[i], all_euk), nrow = 2, dimnames = list(c("this", "all"), c("other", "eukaryotic")))
    test.res <- fisher.test(test.mt, alternative = "less")
    out <- c(out,test.res$p.value)
  }
  return(out)
}

cat("Summarizing hits...\n")
aln[is.na(domain), domain := "other"]
aln[, p95bits := quantile(bits, 0.95), by = query]

# Create wide tables
wide_all <- dcast(aln[evalue < opt$eval_cut, .N, by = .(query, domain)], query ~ domain, fill = 0)
wide_quant <- dcast(aln[bits >= p95bits, .N, by = .(query, domain)], query ~ domain, fill = 0)

# Run enrichment tests
cat("Running Fisher tests...\n")
for (dt in list(wide_all, wide_quant)) {
  setnafill(dt, fill = 0)
  dt[, pvalue := euk.fisher.test(Eukaryota, Archaea + Bacteria + other, all_euk, all_other)]
  dt[, adj.pvalue := p.adjust(pvalue, method = "fdr")]
}

# Save
cat("Saving results...\n")
fwrite(wide_all, paste0("results/", opt$out, "_all.tsv"), sep = "\t")
fwrite(wide_quant, paste0("results/", opt$out, "_quant.tsv"), sep = "\t")

cat("Done.\n")
