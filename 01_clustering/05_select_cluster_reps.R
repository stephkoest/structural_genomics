# Load necessary libraries
library(data.table)
library(ggplot2)

# Load ESP data
ESP_dt <- fread("data/Meng_Li_2021_ESP_UTF8.txt", sep = "\t")[`asCOG ID` != ""]
setnames(ESP_dt, "asCOG ID", "cog")

# Load Asgard DB asCOG data
Asgard_DB_asCOG_safe_dt <- fread("results/Asgard_DB_asCOG_safe_v2.tsv", sep = "\t")
Asgard_DB_alpaca_hits_dt <- fread("data/cogs_alpaca_res.m8", header = TRUE, sep = "\t",
                                   col.names = c("domain", "query", "fident", "alnlen", "mismatch", "gapopen", 
                                                 "tstart", "tend", "tlen", "qstart", "qend", "qlen", "evalue", "bits"))

# Annotate hits
Asgard_DB_asCOG_safe_domains_dt <- Asgard_DB_asCOG_safe_dt[, .(domain = paste(c(target, paste(c(tstart, tend), collapse = "-")), collapse = "__")), by = .(cog = query, target, tstart, tend)]
Asgard_DB_alpaca_hits_anno_dt <- merge(Asgard_DB_asCOG_safe_domains_dt[, .(cog, domain)], Asgard_DB_alpaca_hits_dt, by = "domain", all.x = TRUE)[cog == query]
Asgard_DB_alpaca_hits_anno_dt[, maxbits := max(bits), by = cog]

# Select best hit
Asgard_DB_alpaca_hits_anno_best_dt <- unique(Asgard_DB_alpaca_hits_anno_dt[bits == maxbits])[, .SD[1], by = cog][, .(group = gsub("__.*", "", domain)), by = .(domain, cog)]

# Merge with annotation
asCOG.anno.dt <- fread("data/asCOGs.2020-10.def.tab", header = F, sep = "\t", col.names = c("query","arCOG", "FunCat","Name", "Annotation"))
Asgard_DB_alpaca_hits_anno_best_full_dt <- merge(Asgard_DB_alpaca_hits_anno_best_dt, asCOG_anno_dt[, .(cog = query, arCOG, FunCat, Name, Annotation)], by = "cog", all = TRUE)

# Mark ESPs
Asgard_DB_alpaca_hits_anno_best_full_dt[, ESP := "NA"]
Asgard_DB_alpaca_hits_anno_best_full_dt[cog %in% ESP_dt$cog, ESP := "ESP"]

# Write output
fwrite(Asgard_DB_alpaca_hits_anno_best_full_dt, file = "results/Asgard_DB_asCOG_representatives.tsv", sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE)

# Scoring de novo members
Asgard_DB_de_novo_safe_dt <- fread("results/Asgard_DB_left_20_blast_c50_cluster.tsv", sep = "\t", header = FALSE, col.names = c("cog", "target"))
Asgard_DB_de_novo_hits_dt <- fread("data/de_novo_alpaca_res.m8", header = TRUE, sep = "\t",
                                    col.names = c("target", "query", "fident", "alnlen", "mismatch", "gapopen", 
                                                  "tstart", "tend", "tlen", "qstart", "qend", "qlen", "evalue", "bits"))

# Annotate de novo hits
Asgard_DB_de_novo_hits_anno_dt <- merge(Asgard_DB_de_novo_safe_dt[, .(cog, target)], Asgard_DB_de_novo_hits_dt, by = "target", all.x = TRUE)[cog == query]
Asgard_DB_de_novo_hits_anno_dt[, maxbits := max(bits), by = cog]

# Select best de novo hit
Asgard_DB_de_novo_hits_anno_best_dt <- unique(Asgard_DB_de_novo_hits_anno_dt[bits == maxbits])[, .SD[1], by = cog][, .(group = gsub("__.*", "", target)), by = .(target, cog)]

# Write output
fwrite(Asgard_DB_de_novo_hits_anno_best_dt, file = "results/Asgard_DB_de_novo_representatives.tsv", sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE)
