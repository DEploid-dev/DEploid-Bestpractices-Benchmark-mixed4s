# DEploid-Lasso-Benchmark-mixed4s

Use
3D7(PG0051-C)
HB3(PG0052-C)
DD2(PG0008-CW)
7G8(PG0083-C)
GB4(PG0084-C)
make artificial mixtures

sample=3D7
sample=HB3
sample=DD2
sample=7G8
sample=GB4

vcftools --gzvcf ${sample}.wg.vcf.gz --positions tmp.pos --recode --out ${sample}
mv ${sample}.recode.vcf ${sample}.eg.vcf
bgzip ${sample}.eg.vcf

cb_3D7_HB3_DD2_7G8.ref









samples = c("GB4", "HB3", "DD2", "7G8")

ref = 0
alt = 0
for (sam in samples){
    tmpCoverage = extractCoverageFromVcf(paste(sam, ".eg.vcf.gz", sep =""))
    ref = ref + tmpCoverage$refCount
    alt = alt + tmpCoverage$altCount
}










newvcf = full.vcf[kept.idx,]
newvcf$FORMAT = "AD"
newvcf[["GB4"]] = NULL
newvcf$cb_GB4_HB3_DD2_7G8 = paste(ref, alt, sep = ",")
write.table(newvcf, file = "cb_GB4_HB3_DD2_7G8.vcf", sep ="\t", quote = F, row.names=F)



dEploid -vcf cb_3D7_HB3_DD2_7G8.vcf -panel ~/DEploid-Lasso-Benchmark/combined_pf3k_field_panel_withcross.txt -plaf ~/DEploid-Lasso-Benchmark/labStrains.eg.PLAF.txt -best
dEploid -vcf cb_3D7_HB3_DD2_7G8.vcf -panel ~/DEploid-Lasso-Benchmark/labStrains.eg.panel.txt -plaf ~/DEploid-Lasso-Benchmark/labStrains.eg.PLAF.txt -best
dEploid -vcf cb_GB4_HB3_DD2_7G8.vcf -panel ~/DEploid-Lasso-Benchmark/combined_pf3k_field_panel_withcross_withparent.txt -plaf ~/DEploid-Lasso-Benchmark/labStrains.eg.PLAF.txt -best



c("3D7", "HB3", "DD2", "7G8")
w1: c(.1,.2,.3,.4)
w2: c(.25, .25, .25, .25)
w3: c(.2, .2, .2, .4)
w4: c(.3, .3, .3, .1)

