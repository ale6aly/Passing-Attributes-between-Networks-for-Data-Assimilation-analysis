# AUTHOR: Fahad Almsned
# email: almsned.fahad@gmail.com

library(pandaR)
library(dplyr)
library(tibble)
library(annotables)
library(edgeR)

#-------
#data(pandaToyData)
#pandaResult <- panda(pandaToyData$motif, pandaToyData$expression)
#motif = pandaToyData$motif
#expresstion=pandaToyData$expression
#-----

motifs = read.delim("motifs.tsv",header = FALSE,sep = "")

grch38 %>% dplyr::select(ensgene, symbol) %>% head


#MS network
counts = read.delim("MSmymatrix.txt",sep = "")
counts = rownames_to_column(counts) 
colnames(counts)[colnames(counts) == 'rowname'] <- 'gene'
names(counts)

counts = counts %>%
        dplyr::inner_join(grch38, by = c("gene" = "ensgene"))%>%
        dplyr::select(symbol, SRR3146479.bam, SRR3146478.bam, SRR3146477.bam, SRR3146482.bam, SRR3146481.bam, SRR3146480.bam, SRR3146485.bam, SRR3146484.bam, SRR3146483.bam, SRR3146491.bam, SRR3146490.bam, SRR3146489.bam)

counts = aggregate(. ~ symbol, counts, sum)
counts = column_to_rownames(counts, 'symbol')

dge <- DGEList(counts=counts)
isexpr <- rowSums(cpm(dge) > 10) >= 2
dge <- dge[isexpr,]
dge <- calcNormFactors(dge)
logCPM <- cpm(dge, log=TRUE, prior.count=3)

MSpandaResult <- panda(motifs, logCPM, remove.missing.genes = TRUE, remove.missing.motif = TRUE, progress = TRUE)
#MSpandaResult <- panda(motifs, logCPM, remove.missing.genes = TRUE, remove.missing.motif = TRUE, hamming = 1)
                       

#Normal Network
counts = read.delim("Nmymatrix.txt",sep = "")
counts = rownames_to_column(counts) 
colnames(counts)[colnames(counts) == 'rowname'] <- 'gene'
names(counts)

counts = counts %>%
        dplyr::inner_join(grch38, by = c("gene" = "ensgene"))%>%
        dplyr::select(symbol,SRR3146470.bam, SRR3146469.bam, SRR3146468.bam, SRR3146473.bam, SRR3146472.bam, SRR3146471.bam, SRR3146476.bam, SRR3146475.bam, SRR3146474.bam)

counts = aggregate(. ~ symbol, counts, sum)
counts = column_to_rownames(counts, 'symbol')

dge <- DGEList(counts=counts)
isexpr <- rowSums(cpm(dge) > 10) >= 2
dge <- dge[isexpr,]
dge <- calcNormFactors(dge)
logCPM <- cpm(dge, log=TRUE, prior.count=3)

NpandaResult <- panda(motifs, logCPM,remove.missing.genes = TRUE, remove.missing.motif = TRUE, progress = TRUE)
#NpandaResult <- panda(motifs, logCPM,remove.missing.genes = TRUE, remove.missing.motif = TRUE, hamming = 1)
#-----

plotZ(MSpandaResult,NpandaResult)
#----
calcDegree(MSpandaResult,type = "tf",trim = TRUE)
calcDegree(MSpandaResult,type = "gene",trim = TRUE)
calcDegree(NpandaResult,type = "tf",trim = TRUE)
calcDegree(NpandaResult,type = "gene",trim = TRUE)
#-----
calcDegreeDifference(MSpandaResult,NpandaResult, type= c("tf"),trim = TRUE)
#-----

mat = cbind(rep(1:5,each=10),rep(seq(11,20),5), sample(100,50)/100)
x = plotCommunityDetection(mat)
set.seed(1)
subset = sample(50,10)
mat[subset,3]=subset
plotCommunityDetection(mat,scaleEdge = 0.5)
#----

topMSpandaRes = topedges(MSpandaResult,3000)
subnet.pandaRes = subnetwork(topMSpandaRes, c("SP1", "NFKB1", "RELA", "TP53", "E2F1"))
plotGraph(subnet.pandaRes)

library(reshape2)
v = melt(subnet.pandaRes)
v  = v[v$value>0,]

write.csv(v, file = "v.csv")

x = plotCommunityDetection(v)

