library( "DESeq2" )
library( "GenomicFeatures" )
library( "Rsamtools" )
library( "GenomicAlignments" )
library( "EnhancedVolcano" )

#SRR3734796, SRR3734797, SRR3734798 (1 mo, control) versus SRR3734820, SRR3734824, SRR3734825 (4 mo, test)
fls <- list.files( "/scratch16/lflorea1/projects/Coursera/M1/Sorted/", pattern="bam$", full=TRUE )  ## order?

hse <- makeTxDbFromGFF( "/scratch16/lflorea1/projects/Coursera/M1/Psiclass/WithGeneNames/psiclass_vote.withStrand.gtf", format="gtf" )
#hse <- makeTxDbFromGFF( "/scratch16/lflorea1/projects/Coursera/M1/gencode.vM30.annotation.gtf", format="gtf" )
exonsByGene <- exonsBy( hse, by="gene" ) ##
fls_subset <- fls[c(1:6)] ##

bamLst <- BamFileList( fls_subset, yieldSize=100000 )

se <- summarizeOverlaps(features=exonsByGene, reads=bamLst, mode="Union", singleEnd=FALSE, ignore.strand=TRUE, fragments=TRUE)

write.csv( assay(se), file="raw_counts.1mo_4mo.psiclass.csv")

sampleInfo <- read.csv( "pheno_1mo_4mo.csv" ) ##

seIdx <- match(colnames(se), sampleInfo$sample)
head( cbind( colData(se), sampleInfo[ seIdx, ] ) )
colData(se) <- cbind( colData(se), sampleInfo[ seIdx, ] )

dds <- DESeqDataSet( se, design = ~ condition)
dds$condition <- relevel( dds$condition, "1mo" )
dds <- DESeq(dds)

res <- results(dds)
res_pval0.05 <-res[ which(res$pvalue <= 0.05 ), ]
res_qval0.05 <-res[ which(res$padj <= 0.05 ), ]

write.csv( round(counts(dds, normalized=TRUE),digits=2), file="normalized_counts.1mo_4mo.psiclass.csv")

write.csv( as.data.frame(res), file="results_1mo_4mo.psiclass.csv" )
write.csv( as.data.frame(res_pval0.05), file="results_1mo_4mo.pval05.psiclass.csv" )
write.csv( as.data.frame(res_qval0.05), file="results_1mo_4mo.qval05.psiclass.csv" )

###############  Visualizations  ##################
#PCA plot
vsd <- vst(dds, blind=TRUE)
#rld <- rlog(dds, blind=FALSE)
head(assay(vsd), 3)
png('PCA_1mo_4mo_byCondition.png', width = 720, height = 720)
plotPCA(vsd, intgroup="condition")
dev.off()

vsd_mat <- assay(vsd)
pca <- prcomp(t(vsd_mat))
write.csv(pca$x, file ="1mo_4mo_pca_table.psiclass.csv")
summary(pca)

#volcano plot
gl=c("")
png("volcano_1mo_4mo.psiclass.png", width = 720, height = 720)
#ylim = c(0, 30)
myV <- EnhancedVolcano(res, lab=NA, selectLab = gl, labSize = 3.0, title = "1mo_4mo.psiclass", x = 'log2FoldChange', y = 'padj', pCutoff = 0.05, FCcutoff = 1.5, titleLabSize=24, cutoffLineWidth=0.2, ylab = bquote(~-Log[10]~italic(Padj)), legendLabSize=-1, legendIconSize=-1, xlim=c(-10,10) )
print(myV);
dev.off()
