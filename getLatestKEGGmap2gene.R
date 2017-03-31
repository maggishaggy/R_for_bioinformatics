# If you don't have the following Bioconductor packages installed, uncomment the following line:
## try http:// if https:// URLs are not supported
# source("https://bioconductor.org/biocLite.R")
# biocLite("KEGGREST")

library(KEGGREST)
listDatabases()

# generate up-to-date gene-set/pathway for the given organism

# human
humanKEGGmapping = keggLink("hsa", "pathway")
write.table(c(list(names(humanKEGGmapping), humanKEGGmapping)), "humanPathwayGene.txt", 
            sep="\t", quote = FALSE, col.names=FALSE, row.names = FALSE) 

pathwayidhuman = unique(substr(names(humanKEGGmapping), 6, 13))

for (gs in pathwayidhuman){
  keggMapName = keggGet(gs)[[1]]$NAME
  write.table(c(list(gs, keggMapName)), "humanPathwayName.txt", quote = FALSE, col.names = FALSE, row.names = FALSE, append = TRUE)
}
