#!/bin/bash
dataDir="/well/mcvean/joezhu/DEploid-Bestpractices-Benchmark-mixed4s/"
root="/well/mcvean/joezhu/DEploid-Lasso-Benchmark/"
function run_dEploid {
    mkdir ${root}"dEploidOut/"${sample}
    cd ${root}"dEploidOut/"${sample}
    #rm -r ErrFiles OutFiles
    mkdir ErrFiles OutFiles

for vqslod in $(seq 8 11); do
echo "
#!/bin/bash
#$ -cwd
#$ -V
#$ -P mcvean.prjc -q short.qc
#$ -e ErrFiles
#$ -o OutFiles
#$ -N ${sample}Ex
#$ -t 1-30


plaf=${dataDir}labStrains.eg.PLAF.txt
excludeAt=${dataDir}exclude.txt
vcf=${dataDir}in_silico_mixed_4s/${sample}.vcf.gz

panel=${root}crossPanel_noParent.txt

prefix=${sample}_crossPanel_vqslod${vqslod}_seed\${SGE_TASK_ID}k$@
common=\"-vcf \${vcf} -plaf \${plaf} -exclude \${excludeAt} -o \${prefix} -best\"
dEploidCommon=\"\${common} -seed \${SGE_TASK_ID} -nSample 500 -rate 8 -burn 0.67 -sigma 1.6 -vqslod ${vqslod}\"
rCommon=\"\${common} -dEprefix \${prefix}\"

time dEploid \${dEploidCommon} -panel \${panel} -k $@
## &> ${root}/dEploidOut/\${sample}/\${prefix}.time
interpretDEploid.r \${rCommon}

" > ${sample}k$@_vqslod${vqslod}.sh
    qsub ${sample}k$@_vqslod${vqslod}.sh
done
}

while read sample ; do
    run_dEploid 3
    run_dEploid 4
done < simSampleNames


