## Assignment 1 analysis script 

library(GEOquery)

## getting the decompressed GEO series 
gse <- getGEO(filename='data/GSE50697_family.soft') # data already downloaded to local file in the interest of time
names(GSMList(gse)) # getting the name of the samples in the dataset 

for (gsm in GSMList(gse)) { # this is the condition we are interested in: whether we have miR-203 or not
  print(Meta(gsm)[['characteristics_ch1']])
}

miR_treatment <- function(gsm) { 
  Meta(gsm)[['characteristics_ch1']][2]
}
sapply(GSMList(gse),miR_treatment) # this steps returns data like a loop from function 

pd <- data.frame(miR_treatment=as.factor(sapply(GSMList(gse),miR_treatment))) # this step saves what we are interested into a dataframe

pd$miR_treatment <- as.factor(pd$miR_treatment) # simplification of the info in column
levels(pd$miR_treatment) <- c("control","miR203") # this changes the entry into a readable format

# this matches the phenodata as we read the raw data 
celfiles <- paste0('data/',rownames(pd),'.CEL')
affydata <- read.affybatch(celfiles,phenoData = new("AnnotatedDataFrame",pd))
phenoData(affydata)

eset <- rma(affydata) # performing RMA preprocessing with our raw data
plotDensity(exprs(eset),xlab='log intensity',main="feature level densities after RMA",lwd=2) # plotting the data

pData(eset) # viewing the phenodata

library(limma) # loading the limma library for batch analysis 
model <- model.matrix( ~ 0 + eset$miR_treatment) # linear model dependent on 0 [deduce the intercept]; model matrix to specify the design
colnames(model) <- levels(eset$miR_treatment) 

# look at where the miR treatments differ
contrasts <- makeContrasts(miR203 - control, levels=model) # contrast matrix for stat testing (difference between control/miR203 data)

# fitting our data to the models we have built  
fit <- lmFit(eset, model) 
fitted.contrast <- contrasts.fit(fit,contrasts) 
fitted.ebayes <- eBayes(fitted.contrast) 

# creating a character matrix of the top 10 probe sets
ps <- rownames(topTable(fitted.ebayes))  

## Mapping probe names to gene names and other identifiers
# download the hug133plus2 annotation set if needed
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("hgu133plus2.db", version = "3.8")

library(hgu133plus2.db) # import library for annotation

AnnotationDbi::select(hgu133plus2.db,ps,c("SYMBOL","ENTREZID","GENENAME"),keytype="PROBEID") 
# the probe ID is a probe set ID which does not tell us any biological information (e.g. which gene does the probeID corresponds to)
# select from a database, using the first 10 probes, and gives it a vector of the column names, and pass in the key; specifying)

# Getting a table of upregulated genes 
ps2 <- topTable(fitted.ebayes,number = Inf,p.value = 0.05,lfc=1)
ps2_up <- rownames(ps2[ps2$logFC > 0,]) # positive fold changes 
upregulated <- AnnotationDbi::select(hgu133plus2.db,ps2_up,c("SYMBOL","ENTREZID","GENENAME"),keytype="PROBEID")

# Getting a table of downregulated genes 
ps2_down <- rownames(ps2[ps2$logFC < 0,]) # negative fold changes 
downregulated <- AnnotationDbi::select(hgu133plus2.db,ps2_down,c("SYMBOL","ENTREZID","GENENAME"),keytype="PROBEID")

# export as CSV file and analyse with David (KEGG pathway)
write.csv(upregulated, "Upregulated_Genes_LFC1.csv") 
write.csv(downregulated,'Downregulated_Genes_LFC1.csv')

# Generating the volcano plot with genes with more than 2-times LFC circled in red
volcanoplot(fitted.ebayes, main=sprintf("%d features pass our cutoffs",nrow(interesting_genes)))
points(interesting_genes[['logFC']],-log10(interesting_genes[['P.Value']]),col='red')

# getting the up/down-regulated genes with 2 times LFC for further narrowing down, those circled in red in the volcano plot
interesting_genes <- topTable(fitted.ebayes,number = Inf,p.value = 0.05,lfc=2)
interesting_genes_up <- rownames(interesting_genes[interesting_genes$logFC > 0,])
interesting_genes_down <- rownames(interesting_genes[interesting_genes$logFC < 0,])
interesting_genes_up_mapped <- AnnotationDbi::select(hgu133plus2.db,interesting_genes_up,c("SYMBOL","ENTREZID","GENENAME"),keytype="PROBEID")
interesting_genes_down_mapped <-  AnnotationDbi::select(hgu133plus2.db,interesting_genes_down,c("SYMBOL","ENTREZID","GENENAME"),keytype="PROBEID")

# exporting the genes into CSV files for better visualization
write.csv(interesting_genes_up_mapped,'Interesting_2_fold_genes_mapped_up.csv')
write.csv(interesting_genes_down_mapped,'Interesting_2_fold_genes_mapped_down.csv')
