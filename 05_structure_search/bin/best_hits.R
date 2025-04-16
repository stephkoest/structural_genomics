### best_hits.R
# Author: Stephan
# Description: Extracting domain specific best hits with input/output flexibility

# Load required libraries
suppressWarnings(suppressMessages({
  library(data.table)
  library(stringi)
  library(optparse)
}))

# Define argument parsing options
option_list <- list(
  make_option(c("-i", "--input"), type = "character", default = NULL, 
              help = "Path to input file (.m8 format)", metavar = "character"),
  make_option(c("-o", "--output"), type = "character", default = NULL, 
              help = "Path to output file (extended m8 format)", metavar = "character"),
  make_option(c("-t", "--threshold"), type = "double", default = 0.0001, 
              help = "E-value threshold (default: 0.0001)", metavar = "double")
)

# Parse arguments
opt <- parse_args(OptionParser(option_list = option_list))
if (is.null(opt$input) || is.null(opt$output)) {
  stop("Input and output file paths must be provided. Use -h for help.")
}

# Function to filter top hits from input data
getTophits <- function(indt){
  top.dt <- NULL
  neighbor.dt <- NULL
  indt <- indt[,.(query,target,fident,alnlen,mismatch,gapopen,qstart,qend,qlen,tstart,tend,tlen,evalue,bits)]
  indt[,tcov:= ((tend - tstart) + 1 ) / tlen]
  indt[,qcov:= ((qend - qstart) + 1 ) / qlen]
  #### add qcov filter
  indt[,topbit:= max(bits), by = .(query)]
  indt[,topstart:= 0][,topend:= 0]
  indt[,besthitcand:=0]
  indt[,besthit:=0]
  indt[,maxtcov:=0]
  indt[,maxqcov:=0]
  indt[bits == topbit,besthitcand:=1]
  indt[besthitcand ==1, maxtcov := max(tcov), by = .(query)]
  indt[besthitcand ==1, maxqcov := max(qcov), by = .(query)]
  indt[besthitcand ==1, maxalnlen:= max(alnlen), by = .(query)]
  #indt[besthitcand ==1 & tcov == maxtcov & qcov == maxqcov, besthit :=1]
  #####HERE YOU HAVE TO MAKE AN EXCEPTION TO DEAL WITH CASES WHERE THE HITS HAVE SAME
  indt[besthitcand ==1 & alnlen == maxalnlen, besthit :=1]
  indt[besthit ==1 ,topstart := tstart]
  indt[besthit ==1,topend := tend]
  indt[,topstart := max(topstart),by = .(query)]
  indt[,topend := max(topend),by = .(query)]
  indt[,topoverlap:=1]
  #identify not overlapping hits
  indt[tend < topstart & topstart > 0, topoverlap := 0]
  indt[tstart > topend, topoverlap := 0]
  top.dt <- indt[besthit ==1]
  neighbor.dt <- indt[topoverlap == 0][,.(query,target,fident,alnlen,mismatch,gapopen,qstart,qend,qlen,tstart,tend,tlen,evalue,bits,tcov, qcov)]
  i <- 0
  while (dim(neighbor.dt)[1] >0) {
    i=i+1
    print(paste(c("Iteration:",i), sep = " "))
    dtlen <- dim(neighbor.dt)[1]
    neighbor.dt[,topbit:= max(bits), by = .(query)]
    neighbor.dt[,besthitcand:=0]
    neighbor.dt[,besthit:=0]
    neighbor.dt[,maxtcov:=0]
    neighbor.dt[,maxqcov:=0]
    #assign new best hit
    neighbor.dt[ bits == topbit,besthitcand:=1]
    neighbor.dt[besthitcand == 1, maxtcov := max(tcov), by = .(query)]
    neighbor.dt[besthitcand == 1, maxqcov := max(qcov), by = .(query)]
    neighbor.dt[besthitcand == 1, maxalnlen:= max(alnlen), by = .(query)]
    #neighbor.dt[besthitcand ==1 & tcov == maxtcov & qcov == maxqcov, besthit :=1]
    neighbor.dt[besthitcand ==1 & alnlen == maxalnlen, besthit :=1]
    neighbor.dt[,topstart:= 0][,topend:= 0]
    neighbor.dt[besthit ==1,topstart := tstart]
    neighbor.dt[besthit ==1,topend := tend]
    neighbor.dt[,topstart := max(topstart),by = .(query)]
    neighbor.dt[,topend := max(topend),by = .(query)]
    neighbor.dt[,topoverlap:=1]
    neighbor.dt[tend < topstart & topstart > 0, topoverlap := 0]
    neighbor.dt[tstart > topend, topoverlap := 0]
    top.dt <- rbind(top.dt,neighbor.dt[besthit ==1])
    neighbor.dt <- neighbor.dt[topoverlap == 0][,.(query,target,fident,alnlen,mismatch,gapopen,qstart,qend,qlen,tstart,tend,tlen,evalue,bits, tcov, qcov)]
    olddtlen <- dtlen
    dtlen <- dim(neighbor.dt)[1]
    if (olddtlen == dtlen){
      break
    }
  }
  return(top.dt)
}

# Generalized pipeline function
process_file <- function(input_path, output_path, threshold) {
  cat("Processing file:", input_path, "\n")
  
  # Read input data
  if (grep(".gz$", input_path)){
    dt <- fread(cmd = paste("pigz -cqd ", input_path), header = TRUE, sep = "\t",
              col.names = c("query", "target", "theader", "fident", "alnlen", "mismatch", "gapopen",
               "qstart", "qend", "qlen", "tstart", "tend", "tlen", "evalue", "bits",
               "taxid", "taxname", "taxlineage"))
  } else{
    dt <- fread(input_path, header = TRUE, sep = "\t",
              col.names = c("query", "target", "theader", "fident", "alnlen", "mismatch", "gapopen",
               "qstart", "qend", "qlen", "tstart", "tend", "tlen", "evalue", "bits",
               "taxid", "taxname", "taxlineage"))
  }
  
  # Preprocess data
  dt[, maskc := 0]
  
  # Get top hits
  filtered_hits <- getTophits(dt[evalue < threshold])
  
  # Write to output
  fwrite(filtered_hits, file = output_path, sep = "\t", quote = FALSE, row.names = FALSE)
  cat("Results saved to:", output_path, "\n")
}

# Main execution
process_file(opt$input, opt$output, opt$threshold)
