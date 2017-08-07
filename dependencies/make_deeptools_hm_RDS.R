#!/usr/bin/env Rscript

######################################
ARGV <- commandArgs(trailingOnly = TRUE)

############ USAGE ##############
if (length(ARGV)!=2){
  print(length(ARGV))
  cat("        usage : Rscript script.R <in.matrix.gz> <out.RDS>") 
  cat("")
  cat("")
  stop("2 parameters are requiered")
}
################################
suppressMessages(library(data.table))

######## read ARGV #########
hm  <- ARGV[1]
out <- ARGV[2]
############################



read.deeptools.hm.DT <- function(file.pth, nrows = NULL){
  require(data.table)
  
  # load header to know parameter  
  header <- data.table(t(fread(paste0("zcat ", file.pth, "|head -n 1"), nrows = 1, sep = ",")))
  header <- header[, tstrsplit(V1, ":")]
  
  #remove wrong character
  header[, V1 := sub("@\\{", "", V1)]
  header[, V2 := sub("\\}", "", V2)]
  
  # convert to integer
  header[,V2 := as.integer(V2)]
  
  
  downstream <-header[V1 %like% "downstream", V2]
  unscaled_3_prime <- header[V1 %like% "unscaled 3 prime",V2]
  body <-  header[V1 %like% "body",V2]
  unscaled_5_prime <- header[V1 %like% "unscaled 5 prime",V2]
  upstream <- header[V1 %like% "upstream",V2]
  bin <- header[V1 %like% "bin size",V2]
  ref_point <- header[V1 %like% "ref point",V2]
  
  ## load deeptools <file>.matrix.gz
  if(is.null(nrows)){
    hm <- fread(paste0("zcat ", file.pth), skip = 1) 
  }
  else{
    hm <- fread(paste0("zcat ", file.pth, "| head -n ", nrows+1), skip = 1)
  }
    
  
  colnames(hm)[1:6] <- c("chr", "start", "end", "rn", "length", "strand")
  
  
  cname.upstream <- upstream*-1
  cname.body <- c("TSS", "TES")
  cname.downstream <- downstream
  
  if(body != 0){
    seq.upstream <- c(paste0(seq(cname.upstream+bin, -bin, bin)/1000,"kb"), "TSS")
    seq.downstream <- c("TES",paste0(seq(bin, cname.downstream, bin)/1000,"kb"))
    percent.bin <- 100/(body/bin)
    seq.body = paste0(seq(0+percent.bin, 100-percent.bin, percent.bin),"%")
    
    cname <- c(seq.upstream, seq.body, seq.downstream)  
  }
  else{
    cname <- paste0(seq(cname.upstream+bin,cname.downstream, by = bin)/1000, "kb")
  }
  
  colnames(hm)[7:ncol(hm)] <- cname
#   
#   summary(hm)
#   hm <- as.matrix(hm)
#   hm[] <- 0
#   hm <- as.data.frame(hm)
  hm
}


hm <- read.deeptools.hm.DT(hm)

saveRDS(hm, out)
