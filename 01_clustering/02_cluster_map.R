getTophits <- function(indt){
  require(data.table)
  top.dt <- NULL
  neighbor.dt <- NULL
  if (!"DB" %in% colnames(indt)){
    indt[,DB := "NA"]
  }
  if (!"maskc" %in% colnames(indt)){
    indt[,maskc := 0]
  }
  #reorder columns
  indt <- indt[,.(query,target,fident,alnlen,mismatch,gapopen,qstart,qend,qlen,tstart,tend,tlen,evalue,bits,DB,maskc)]
  
  #calculate coverages
  indt[,tcov:= ((tend - tstart) + 1 ) / tlen]
  indt[,qcov:= ((qend - qstart) + 1 ) / qlen]
  
  #identify top hits
  indt[,topevalue:= min(evalue), by = .(target)]
  indt[,topbit:= max(bits), by = .(target)]
  indt[,topstart:= 0][,topend:= 0]
  indt[,besthitcand:=0]
  indt[,besthit:=0]
  indt[,maxtcov:=0]
  indt[,maxqcov:=0]
  #flag best hit candidates
  indt[evalue == topevalue & bits == topbit,besthitcand:=1]
  indt[besthitcand ==1, maxtcov := max(tcov), by = .(target)]
  indt[besthitcand ==1, maxqcov := max(qcov), by = .(target)]
  indt[besthitcand ==1, maxalnlen:= max(alnlen), by = .(target)]
  #####HERE YOU HAVE TO MAKE AN EXCEPTION TO DEAL WITH CASES WHERE THE HITS HAVE SAME
  indt[besthitcand ==1 & alnlen == maxalnlen, besthit :=1]
  indt[besthit ==1 ,topstart := tstart]
  indt[besthit ==1,topend := tend]
  indt[,topstart := max(topstart),by = .(target)]
  indt[,topend := max(topend),by = .(target)]
  indt[,topoverlap:=1]
  #identify not overlapping hits
  indt[tend < topstart & topstart > 0, topoverlap := 0]
  indt[tstart > topend, topoverlap := 0]
  top.dt <- indt[besthit ==1]
  #collect non overlapping hits
  neighbor.dt <- indt[topoverlap == 0][,.(query,target,fident,alnlen,mismatch,gapopen,qstart,qend,qlen,tstart,tend,tlen,evalue,bits,DB,maskc, tcov, qcov)]
  i <- 0
  #select additional hits until no more hits are found
  while (dim(neighbor.dt)[1] >0) {
    i=i+1
    print(paste(c("Iteration:",i), sep = " "))
    dtlen <- dim(neighbor.dt)[1]
    neighbor.dt[,topevalue:= min(evalue), by = .(target)]
    neighbor.dt[,topbit:= max(bits), by = .(target)]
    neighbor.dt[,besthitcand:=0]
    neighbor.dt[,besthit:=0]
    neighbor.dt[,maxtcov:=0]
    neighbor.dt[,maxqcov:=0]
    #assign new best hit
    neighbor.dt[evalue == topevalue & bits == topbit,besthitcand:=1]
    neighbor.dt[besthitcand == 1, maxtcov := max(tcov), by = .(target)]
    neighbor.dt[besthitcand == 1, maxqcov := max(qcov), by = .(target)]
    neighbor.dt[besthitcand == 1, maxalnlen:= max(alnlen), by = .(target)]
    neighbor.dt[besthitcand ==1 & alnlen == maxalnlen, besthit :=1]
    neighbor.dt[,topstart:= 0][,topend:= 0]
    neighbor.dt[besthit ==1,topstart := tstart]
    neighbor.dt[besthit ==1,topend := tend]
    neighbor.dt[,topstart := max(topstart),by = .(target)]
    neighbor.dt[,topend := max(topend),by = .(target)]
    neighbor.dt[,topoverlap:=1]
    neighbor.dt[tend < topstart & topstart > 0, topoverlap := 0]
    neighbor.dt[tstart > topend, topoverlap := 0]
    top.dt <- rbind(top.dt,neighbor.dt[besthit ==1])
    neighbor.dt <- neighbor.dt[topoverlap == 0][,.(query,target,fident,alnlen,mismatch,gapopen,qstart,qend,qlen,tstart,tend,tlen,evalue,bits,DB,maskc, tcov, qcov)]
    #check if no more hits are found
    olddtlen <- dtlen
    dtlen <- dim(neighbor.dt)[1]
    if (olddtlen == dtlen){
      break
    }
  }
  return(top.dt)
}
#default query coverage cutoff
query_coverage_cutoff <- 0.8
#read in mmseqs output file
asCOG_seq_db_asCOG_domains_seq_db.dt <- fread(cmd = "pigz -cqd results/Asgard_DB_asCOG_seq_domains_db.m8.gz", header = T, sep = "\t",
                                              col.names = c("target","query","fident", "alnlen", "mismatch","gapopen", "tstart", "tend","tlen", "qstart", "qend", "qlen", "evalue","bits"))
#read in listasgard proteins
asgard_prots.dt <- fread("data/Asgard_DB.txt", colnames = c("target")))
#add for legacy reasons
asCOG_seq_db_asCOG_domains_seq_db.dt[,DB:= "asCOG_seq"][,maskc:= 0]
#find top hits
asCOG_seq_db_asCOG_domains_seq_db.top.dt <-  getTophits(asCOG_seq_db_asCOG_domains_seq_db.dt)
asCOG_seq_db_asCOG_domains_seq_db.top.dt[,queryseq := query]
asCOG_seq_db_asCOG_domains_seq_db.top.dt[,query := gsub("__.*","", queryseq)]
asCOG_seq_db_asCOG_domains_seq_db.top.dt$queryseq <- NULL

#filter out hits with query coverage below threshold
Asgard_DB_asCOG_safe.dt <- unique(asCOG_seq_db_asCOG_domains_seq_db.top.dt[qcov >= query_coverage_cutoff])
Asgard_DB_left.dt <- asgard_prots.dt[!target %in% unique(Asgard_DB_asCOG_safe.dt$target)]

#identify proteins with multiple asCOG domains
Asgard_DB_asCOG_safe.multi.v <- unique(Asgard_DB_asCOG_safe.dt[,.N,by = .(target)][N>1]$target)
Asgard_DB_asCOG_safe.multi.dt <- Asgard_DB_asCOG_safe.dt[target %in% Asgard_DB_asCOG_safe.multi.v]
setorderv(Asgard_DB_asCOG_safe.multi.dt, c("target", "tend"))

#identify space between asCOG domains
space.min <- 60
multi.dt <- data.table("target" = character(), "start" = numeric(),"stop" = numeric())
for (prot in unique(Asgard_DB_asCOG_safe.multi.dt$target)){
  tmp.dt <- Asgard_DB_asCOG_safe.multi.dt[target == prot,.(target, tstart, tend,tlen)]
  tlength <- tmp.dt$tlen[1]
  temp.v <- c(as.vector(t(as.matrix(tmp.dt[,.(tstart,tend)]))))
  tstart <- temp.v[1]
  tend <- temp.v[length(temp.v)]
  temp.v <- temp.v[c(-1,-length(temp.v))]
  if (tlength - tend >= space.min) multi.dt <- rbind(multi.dt, data.table("target" = prot, "start" = tend + 1, "stop" = tlength))
  if (tstart - 1 >= space.min) multi.dt <- rbind(multi.dt, data.table("target" = prot, "start" = 1, "stop" = tstart - 1))
  for (i in seq(1, length(temp.v), by = 2)){
    interspace <- temp.v[i + 1] - temp.v[i] - 1
    if (interspace >= space.min){
      start <- temp.v[i] + 1
      stop <- temp.v[i +1] - 1
      multi.dt <- rbind(multi.dt,data.table("target" = prot, "start" = start, "stop" = stop))
    }
  }
}

#identify proteins with only one asCOG domain
Asgard_DB_asCOG_safe.single.v <- unique(Asgard_DB_asCOG_safe.dt[,.N,by = .(target)][N==1]$target)
Asgard_DB_asCOG_safe.single.dt <- Asgard_DB_asCOG_safe.dt[target %in% Asgard_DB_asCOG_safe.single.v]
Asgard_DB_asCOG_safe.single.dt[,tstartleft:= tstart - 1 ]
Asgard_DB_asCOG_safe.single.dt[,tendleft := tlen - tend]

#check for left over space in proteins with assigned asCOG domains
AsCOG_leftover.dt <- Asgard_DB_asCOG_safe.single.dt[tstartleft >= space.min,.("target" = target, "start" = 1, "stop" = tstart -1)]
AsCOG_leftover.dt <- rbind(AsCOG_leftover.dt, Asgard_DB_asCOG_safe.single.dt[tendleft >= space.min, .("target" = target, "start" = tend + 1, "stop"  =  tlen)])
AsCOG_leftover.dt <- rbind(AsCOG_leftover.dt, multi.dt)

#create output files
write.table(unique(Asgard_DB_asCOG_safe.dt), file = "results/Asgard_DB_asCOG_safe_v2.tsv",sep = "\t", quote = F, row.names = F, col.names = T)
write.table(unique(AsCOG_leftover.dt), file = "results/Asgard_DB_asCOG_left_partial_proteins_v2.tsv",sep = "\t", quote = F, row.names = F, col.names = T)
write.table(unique(Asgard_DB_left.dt[,.(seq = target)]), file = "results/Asgard_DB_asCOG_unmapped_v2.txt",sep = "\t", quote = F, row.names = F, col.names = F)
