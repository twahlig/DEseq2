library(Rsubread)
library(DESeq2)

setwd("/Users/twahlig/Desktop/RNAseq/")

buildindex(basename="E40V_index",reference="E40V.fa")
#this builds the index file

align(index="E40V_index",readfile1="TSB-E40V-1-005.fastq",output_file="E40V_1.bam")
align(index="E40V_index",readfile1="TSB-E40V-2-004.fastq",output_file="E40V_2.bam")
align(index="E40V_index",readfile1="TSB-E40V-3-001.fastq",output_file="E40V_3.bam")
align(index="E40V_index",readfile1="TSB-E40-1-003.fastq",output_file="E40_1.bam")
align(index="E40V_index",readfile1="TSB-E40-2-002.fastq",output_file="E40_2.bam")
align(index="E40V_index",readfile1="TSB-E40-3-006.fastq",output_file="E40_3.bam")
#aligns reads to index
#this takes a while (3-5 hours)

countfile <- featureCounts(files=c("E40V_1.bam","E40V_2.bam","E40V_3.bam","E40_1.bam","E40_2.bam","E40_3.bam"),annot.ext="/Users/twahlig/Desktop/RNAseq/E40V_annotation.gtf",GTF.featureType="CDS",isGTFAnnotationFile=TRUE,useMetaFeatures=FALSE,reportReads=TRUE)
#depending on format of annotation file, you may have to modify your gff or gtf file so the header is properly recognized 9th column needs to be gene_id 10th column needs to start with either cds[number] or gene[number]


test <- as.data.frame(countfile$counts)
#this makes the dataframe called 'test' using the counts in countfile. countfile contains 
#other objects, so you have to designate $counts 


attribute <- read.delim("attributes.txt", sep = "\t", header = TRUE, row.names = 1)
#this creates the object 'atribute' using the text file attributes (which you have to make), which is tab separated (\t), has a header (header = true) 

results <-DESeqDataSetFromMatrix (countData = test, colData = attribute, design = ~Group)
#this creates the deseq dataset, 'test' was the name given in the first step, colData is the name of the columns in the attribute file

results <- estimateSizeFactors (results)

results <- estimateDispersions (results)
#estimates variability

results <- nbinomWaldTest(results)
#calculates variance using nBinomWaldTest

final_results <- counts (results, normalized = TRUE)
#this creates final results object that contains normalized read counts that were put in 'results'

testresults <- DESeq(results)

results <- results(testresults)

write.table(results, "results.txt")
