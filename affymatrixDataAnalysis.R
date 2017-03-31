###########################Using affy to process the raw Affymatrix CEL file
# If you don't have the following Bioconductor packages installed, uncomment the following line:
# source("http://bioconductor.org/biocLite.R")
# biocLite("affy");
# biocLite("arrayQualityMetrics")
# biocLite("GEOquery")
# biocLite("gageData")
# biocLite("KEGGREST")


library(affy)
library(arrayQualityMetrics)
library(limma)
library(ggplot2)
library(Biobase)
library(GEOquery)
library(pathview)
library(gage)
library(gageData)
library(KEGGREST)

# read in the meta-data file
arrayMetaData <- read.table("metadata.tab", header=T, check.names=F, sep="\t")

# load the raw array data
# read CEL file and make a normalization using rma
affyRawDataset <- justRMA(filenames = arrayMetaData$filename,
               celfile.path = "./array_data/",
               sampleNames = arrayMetaData$sample.lable,
               compress=getOption("BioC")$affy$compress.cel
               )


summary(exprs(affyRawDataset))

featureNames(affyRawDataset)

boxplot(exprs(affyRawDataset), at=1:8, xlim=c(0.5, 8.5), ylim=c(1,18), log="", add = F, las =3)

# quality control for outlier detection
arrayQualityMetrics(expressionset = affyRawDataset, outdir = "./QC/Report_for_GSE67664", force = TRUE)



expressData <- exprs(affyRawDataset)
write.table(cbind(Probe_ID = rownames(affyRawDataset), expressData), sep="\t", row.names = FALSE,
            quote = FALSE, file="expression_rma_probe_matrix.txt")


# extract the differentially expressed genes
group=factor(c(rep("activated",4),rep("quiescent",4)))
design = model.matrix(~0+group)
colnames(design)=c("aHSC","qHSC")
rownames(design)=sampleNames(affyRawDataset)

fit=lmFit(affyRawDataset,design)
cont.matrix = makeContrasts(contrasts="aHSC-qHSC",levels=design)
fit2=contrasts.fit(fit,cont.matrix)
fit2=eBayes(fit2)
diff_dat=topTable(fit2,coef=1,n=Inf)
write.table(diff_dat,"diff_dat.txt",quote=F)


# Download GPL file, put it in the current directory, and load it:
# gpl <- getGEO('GPL13667', destdir=".")
# if GPL13667.soft already downloaded
gpl <- getGEO(filename = "GPL13667.soft")
dim(Table(gpl))
colnames(Table(gpl)) ## [1] 49386 probes, 43 columns
head(Table(gpl)[,c(1,20,21)]) ## you need to check this , which column do you need. In this case, three columns for probe id, ensembl id and entrez id
probe2geneID <- Table(gpl)[,c(1,20,21)]
probe2entrenz <- as.matrix(Table(gpl)[,c(1,21)])
write.csv(Table(gpl)[,c(1,20,21)],"GPL13667.csv")

gene_probe <- Table(gpl)[,c(1,20,21)]
diff_probe=rownames(diff_dat[abs(diff_dat[,1]) > 2,])
diff_gene=gene_probe[match(diff_probe,gene_probe[,1]),3]
diff_gene=na.omit(diff_gene)
diff_gene=unique(diff_gene)
length(diff_gene)

# Over-representation analysis
require(DOSE)
require(clusterProfiler)
library(org.Hs.eg.db)

gene=as.character(diff_gene)
# keep the first gene if multiple genes related to a probe
ego <- enrichGO(gene=gene,
                OrgDb         = org.Hs.eg.db,
                ont="CC",qvalueCutoff=0.01,
                readable=TRUE)
ekk <- enrichKEGG(gene         = gene,
                  organism     = 'hsa',
                  qvalueCutoff = 0.05)
      
write.csv(ekk,"KEGG-enrich.csv",row.names =F)
write.csv(ego,"GO-enrich.csv",row.names =F)


# Gene-set enrichment analysis using GAGE and pathview, and the outdated KEGG pathway database 
entrez.data <- mol.sum(mol.data = expressData, id.map = probe2entrenz, sum.method = "mean")

data(kegg.sets.hs)
data(sigmet.idx.hs)
samp.idx=1:4    # aHSC as treatment
ref.idx=5:8     # qHSC as control

kegg.1d.p <- gage(entrez.data, gset=kegg.sets.hs[sigmet.idx.hs], ref = ref.idx, samp = samp.idx, compare="unpaired")

path.up.full.info <-kegg.1d.p$greater
write.table(cbind(Pathway=rownames(path.up.full.info), path.up.full.info), sep="\t",row.names=FALSE,file="up.path.KEGG.list")

path.down.full.info <-kegg.1d.p$less
write.table(cbind(Pathway=rownames(path.down.full.info), path.down.full.info), sep="\t",row.names=FALSE,file="down.path.KEGG.list")


# plot the significantly perturbed pathways with the changes in gene expression profile
intensity.log.d <- entrez.data[,samp.idx] - rowMeans(entrez.data[,ref.idx])
up.pathway.idx <- path.up.full.info[, "q.val"] < 1E-4 & !is.na(kegg.1d.p$greater[,"q.val"])
path.up.full.ids <- rownames(path.up.full.info)[up.pathway.idx]   # 11 up-regulated pathways

# Display the official gene symbols by setting same.layer = F
for (gs in path.up.full.ids){
  up.pid = substr(gs, 1, 8)
  outname = gsub(" |:|/", "_", substr(gs,9,100))
  pathview(gene.data = intensity.log.d, 
           pathway.id = up.pid,
           limit=list(gene=2,cpd=1), 
           species = "hsa",
           kegg.native = T, 
           same.layer = F,
           kegg.dir="./KEGGmap", 
           gene.idtype="entrez",
           out.suffix = outname)
}

# The original KEGG pathway data file (.xml) and image file (.png) will be saved in this specified directory.
# The newly created KEGG pathway image file can be found at the current directory

# switch to the folder of down-regulation KEGG pathway created just now and then set it as current working directory
down.pathway.idx <- path.down.full.info[, "q.val"] < 1E-4 & !is.na(kegg.1d.p$less[,"q.val"])
path.down.full.ids <- rownames(path.down.full.info)[down.pathway.idx]    # 13 down-regulated pathways


for (gs in path.down.full.ids){
  down.pid = substr(gs, 1, 8)
  outname = gsub(" |:|/", "_", substr(gs,9,100))
  pathview(gene.data = intensity.log.d, 
           pathway.id = down.pid,
           limit=list(gene=2,cpd=1), 
           species = "hsa",
           kegg.native=T, 
           same.layer=F,
           kegg.dir="./KEGGmap", 
           gene.idtype="entrez",
           out.suffix = outname)
}



# Gene-set enrichment analysis using GAGE and pathview, and the up-to-date KEGG pathway database 
humanKEGG <- readList("humanKEGGpathway.gmt")
samp.idx=1:4    # aHSC as treatment
ref.idx=5:8     # qHSC as control

lastestkegg.1d.p <- gage(entrez.data, gset=humanKEGG, ref = ref.idx, samp = samp.idx, compare="unpaired")

path.up.full.info <-lastestkegg.1d.p$greater
write.table(cbind(Pathway=rownames(path.up.full.info), path.up.full.info), sep="\t",row.names=FALSE,file="up.latest.KEGG.list")

path.down.full.info <-lastestkegg.1d.p$less
write.table(cbind(Pathway=rownames(path.down.full.info), path.down.full.info), sep="\t",row.names=FALSE,file="down.latest.KEGG.list")


# plot the significantly perturbed pathways with the changes in gene expression profile
intensity.log.d <- entrez.data[,samp.idx] - rowMeans(entrez.data[,ref.idx])
up.pathway.idx <- path.up.full.info[, "q.val"] < 1E-4 & !is.na(lastestkegg.1d.p$greater[,"q.val"])
path.up.full.ids <- rownames(path.up.full.info)[up.pathway.idx]   # 14 up-regulated pathways

# Display the official gene symbols by setting same.layer = F
for (gs in path.up.full.ids[1:3]){
  up.pid = substr(gs, 1, 8)
  outname = gsub(" |:|/", "_", substr(gs,9,100))
  pathview(gene.data = intensity.log.d, 
           pathway.id = up.pid,
           limit=list(gene=2,cpd=1), 
           species = "hsa",
           kegg.native = T, 
           same.layer = F,
           kegg.dir="./KEGGmap", 
           gene.idtype="entrez",
           out.suffix = outname)
}

down.pathway.idx <- path.down.full.info[, "q.val"] < 1E-4 & !is.na(lastestkegg.1d.p$less[,"q.val"])
path.down.full.ids <- rownames(path.down.full.info)[down.pathway.idx]    # 27 down-regulated pathways


for (gs in path.down.full.ids[1:3]){
  down.pid = substr(gs, 1, 8)
  outname = gsub(" |:|/", "_", substr(gs,9,100))
  pathview(gene.data = intensity.log.d, 
           pathway.id = down.pid,
           limit=list(gene=2,cpd=1), 
           species = "hsa",
           kegg.native=T, 
           same.layer=F,
           kegg.dir="./KEGGmap", 
           gene.idtype="entrez",
           out.suffix = outname)
}



# unsupervised analysis for estimating sample distance
# PCA plot
expressDatat <- t(expressData)
pca <- prcomp(expressDatat, center = TRUE)
percentVar <- pca$sdev^2/sum(pca$sdev^2)
pca1 <- pca$x[,1]
pca2 <- pca$x[,2]

pcadf <- data.frame(ID=names(pca1), PCA1 = pca1, PCA2 = pca2, Cluster=group)

ggplot(pcadf, aes(x=PCA1, y=PCA2, label=ID, color=Cluster)) + 
       geom_point() + xlab(paste0("PC1: ", round(percentVar[1] * 100), "% variance")) + ylab(paste0("PC2: ", round(percentVar[2] * 100), "% variance")) + geom_text(size = 4, colour = "black", vjust = -1)




save.image("dataProcessing.RData")



