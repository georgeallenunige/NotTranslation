library("DESeq2")

# the counts table contains counts for each replicate in soluble and total for WT and not1d (expression filter applied)
countsTable="pathToCountsTable"
cts2=read.table(countsTable)

# Then coldata1 is the input to construct the design matrix:

replicate  soltot  genotype   type         
 "Wt1"     "sol"  "Wt"     "single-read"
 "Wt2"     "sol"  "Wt"     "single-read"
 "Wt3"     "sol"  "Wt"     "single-read"
 "Not1p1"  "sol"  "Not1p"  "single-read"
 "Not1p2"  "sol"  "Not1p"  "single-read"
 "Not1p3"  "sol"  "Not1p"  "single-read"
 "Wt1"     "tot"  "Wt"     "single-read"
 "Wt2"     "tot"  "Wt"     "single-read"
 "Wt3"     "tot"  "Wt"     "single-read"
 "Not1p1"  "tot"  "Not1p"  "single-read"
 "Not1p2"  "tot"  "Not1p"  "single-read"
 "Not1p3"  "tot"  "Not1p"  "single-read"


# Then we perform a ratio of ratios test examining change in sol/tot between strains:

library("DESeq2")
dds <- DESeqDataSetFromMatrix(countData = cts2,
                              colData = coldata1,
                              design = ~ soltot + genotype + soltot:genotype)

ddsF <- DESeq(ddsFO, test="Wald")
resultsNames(ddsF)
res1Dn=results(ddsF,altHypothesis="less",contrast=list("soltottot.genotypeWt"))
res1Up=results(ddsF,altHypothesis="greater",contrast=list("soltottot.genotypeWt"))



