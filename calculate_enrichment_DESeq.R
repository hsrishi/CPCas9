#######################################################
#Final version with combined technical reps

library("DESeq")

##1. Input data and preparations
#count table
CpCas9CountTable = read.csv('df_ExportForDESeq.csv', header = TRUE, row.names = 1) #file name hardcoded

CpCas9CountTable_AD = CpCas9CountTable[c('DFSBLO002A','DFSBLO002C','DFSBLO002D')]


#metadata
CpCas9Design_AD = data.frame(
	row.names = colnames(CpCas9CountTable_AD),
	condition = c("untreated", "treated", "treated"),
	libType = c("single-end","single-end","single-end") ) #metadata hardcoded

singleSamples_AD = CpCas9Design_AD$libType == "single-end"
countTable_AD = CpCas9CountTable_AD[ , singleSamples_AD]
condition_AD = CpCas9Design_AD$condition[ singleSamples_AD]

cdsAD = newCountDataSet( countTable_AD, condition_AD)


#normalization
cdsAD = estimateSizeFactors( cdsAD)
sizeFactors(cdsAD)

head(counts(cdsAD, normalized = TRUE))

##2. Variance estimation
cdsAD = estimateDispersions( cdsAD )
str( fitInfo(cdsAD) )
plotDispEsts( cdsAD )

head( fData(cdsAD) )

##3. Inference: Calling differential expression
resAD = nbinomTest( cdsAD, "untreated", "treated" )
head(resAD)
plotMA(resAD)
hist(resAD$padj, breaks=100, col="skyblue", border="slateblue", main="")
resSigAD = resAD[ resAD$padj < 0.1, ]
head( resSigAD[ order(resSigAD$pval), ] ) #list most significantly differentially expressed variants:
head( resSigAD[ order( resSigAD$foldChange, -resSigAD$baseMean ), ] ) #list most strongly down-regulated of the significant variants
head( resSigAD[ order( -resSigAD$foldChange, -resSigAD$baseMean ), ] ) #list most strongly up-regulated of the significant variants

write.table(resAD, file = "/Users/rishih91/ArkinLabLocal/CpCas9_SIRE_170720_100SRR_HS3B/resACD_allindex.csv", sep = ",", row.names = FALSE)
