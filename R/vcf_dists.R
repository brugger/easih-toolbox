#!/broad/tools/apps/R-2.6.0/bin/Rscript

args <- commandArgs(TRUE)

matrix = read.table(args[1])

for ( i in ( colnames(matrix))) {
   if ( i == "Pos" || 
        i == "dbsnp" ||   
        i == "DB" ||   
        i == "AN" ||   
        i == "X2b_Chi" ||   
        i == "Alignability" || 
        i == "GenericAnnotation" ||
        i == "HW" ||   
        i == "LowMQ" ||   
        i == "MQRankSum" ||   
        i == "Qual_Adjusted_2blod") {
     next
   }

   cat(i, "\n")

   a <- hist(matrix[[i]][matrix$dbsnp=="known"])
  b <- hist(matrix[[i]][matrix$dbsnp=="unknown"])

   outfile = paste(i, ".ClusterReport.pdf", sep="")
   pdf(outfile, height=7, width=8)

   maxY = max(a$intensities, b$intensities)

#   plot(a$mids, a$intensities, ylim=c(0,maxY),type="b",col="orange",lwd=2,xlab=i,ylab="SNPs")
#   points(b$mids, b$intensities, type="b",col="blue",lwd=2)
   plot(a$mids, a$density, ylim=c(0,maxY),type="b",col="orange",lwd=2,xlab=i,ylab="SNPs")
   points(b$mids, b$density, type="b",col="blue",lwd=2)
   legend('topright', c('knowns','novels'),lwd=2,col=c("orange","blue"))
   dev.off()
}
