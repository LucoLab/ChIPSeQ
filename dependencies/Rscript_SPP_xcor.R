#!/usr/bin/env Rscript

#TODO: add param for cross strand corr

#################### DEseq2 ##################
ARGV <- commandArgs(trailingOnly = TRUE)

############ USAGE ##############
if (length(ARGV)!=5){
  print(length(ARGV))
  cat("        usage : Rscript script.R <IP.bam> <output.pdf> <nb.core> <max.shift> <bin>") 
  cat("")
  cat("")
  stop("5 parameters are requiered")
}
################################


######### read ARGV ##############
chip.bam  <- ARGV[1]
out.pdf   <- ARGV[2]
nb.core   <- as.numeric(ARGV[3])
max.shift <- as.numeric(ARGV[4])
bin       <- as.numeric(ARGV[5])
############## ##############

suppressMessages(library(spp));
suppressMessages(library(snow));
#~ cluster <- makeCluster(nb.core);
cluster <- makeSOCKcluster(rep("localhost",nb.core))


cat(paste0("read IP bam file:",chip.bam, "\n"))
chip.data <- read.bam.tags(chip.bam);

cat(paste0("Get binding info from cross-correlation profile", "\n"))
# get binding info from cross-correlation profile
# srange gives the possible range for the size of the protected region;
# srange should be higher than tag length; making the upper boundary too high will increase calculation time
#
# bin - bin tags within the specified number of basepairs to speed up calculation;
# increasing bin size decreases the accuracy of the determined parameters

#binding.characteristics <- get.binding.characteristics(chip.data,srange=c(-100,600),bin=5,cluster=cluster);
binding.characteristics <- get.binding.characteristics(chip.data,srange=c(0,max.shift),bin=bin,cluster=cluster);

# print out binding peak separation distance
cat(paste("binding peak separation distance =",binding.characteristics$peak$x, "\n"))


# plot cross-correlation profile
cat(paste0("Plot cross-correlation profile", "\n"))

pdf(file= out.pdf,width=5,height=5)
par(mar = c(3.5,3.5,1.0,0.5), mgp = c(2,0.65,0), cex = 0.8);
plot(binding.characteristics$cross.correlation,type='l',xlab="strand shift",ylab="cross-correlation");
abline(v=binding.characteristics$peak$x,lty=2,col=2)
legend("topright", legend = binding.characteristics$peak$x) 
dev.off();
