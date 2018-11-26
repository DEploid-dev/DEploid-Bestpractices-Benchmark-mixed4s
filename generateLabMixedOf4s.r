rm(list=ls())
library(DEploid)
samples = c("threeD7", "dd2gt.from.regression", "HB3gt.from.regression", "sevenG8gt.from.regression")

h <- function(w){
         if ( any( grepl( "gzfile connection", w) ) )
         invokeRestart( "muffleWarning" )
    }

gzf <- gzfile("wg_vcf/GB4.wg.vcf.gz", open = "rb")
numberOfHeaderLines <- 0
line <- withCallingHandlers( readLines(gzf, n = 1), warning = h)
while ( length(line) > 0 ){
    if (grepl("##", line )){
        numberOfHeaderLines <- numberOfHeaderLines + 1
    } else {
        break
    }
    line <- withCallingHandlers( readLines(gzf, n = 1), warning = h)
}
close(gzf)

full.vcf <- read.table("wg_vcf/GB4.wg.vcf.gz", skip = numberOfHeaderLines,
    header = T, comment.char = "", stringsAsFactors = FALSE,
    check.names = FALSE)


w = list()
w[[1]] = c(.1,.2,.3,.4)
w[[2]] = c(.25, .25, .25, .25)
w[[3]] = c(.2, .2, .2, .4)
w[[4]] = c(.3, .3, .3, .1)

panel = read.table("labStrains.eg.panel.txt", sep = "\t", header=T)

full.vcf.chrompos = paste(full.vcf[,1], full.vcf[,2], sep ="-")
vcf.chrompos = paste(panel[,1], panel[,2], sep ="-")
kept.idx = which(full.vcf.chrompos %in% vcf.chrompos)
haps = panel[,-c(1,2)]

n.loci = nrow(panel)
tmpCoverage = extractCoverageFromVcf("HB3.eg.vcf.gz")
coverage_profile = tmpCoverage$refCount + tmpCoverage$altCount

for ( wi in 1:4 ){
    newvcf = full.vcf[kept.idx,]
    newvcf$FORMAT = "GT:AD"
    newvcf[["GB4"]] = NULL

    p = w[[wi]]
    WSAF = 0;
    for ( i in 1:4){
        WSAF = WSAF + p[i] * haps[,i]
    }
    err = 0.01
    includeErrorWSAF1 = WSAF*(1-err)+(1-WSAF)*err
    altCount1 = rbinom(n.loci, coverage_profile, includeErrorWSAF1)
    refCount1 = coverage_profile - altCount1

    newvcf$cb_3D7_HB3_DD2_7G8 = paste("./.:", refCount1, ",", altCount1, sep = "")
    write.table(newvcf, file = paste("cb_3D7_HB3_DD2_7G8_w", wi, ".vcf", sep =""), sep ="\t", quote = F, row.names=F)

}
