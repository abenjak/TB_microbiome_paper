
#### R script used to process FIND (1st and 2nd batch) and CHUV samples
#### Published in paper "Multicenter analysis of sputum microbiota in tuberculosis patients"
#### Marion Leleu & Jacques Rougemont 
#### May 2020


### LIBRARIES --------
library("pheatmap")
library("RColorBrewer")
library("colorRamps") 
library("ggrepel")
library("limma")



############################ FUNCTIONS ###############################-----

myCor = function(x) {as.dist(1-cor(t(x), use="pairwise.complete.obs"))}
myLog = function(x) {y=x; y[x<=0]=NA; y[x>0]=log(x[x>0]); y}

centerScale <- function(x){
  mns = apply(x, 1, mean, na.rm=T) 
  x1 = sweep(x, 1, mns, "-")
  sds = sqrt(apply(x1^2, 1, sum, na.rm=T)) 
  x2 = as.matrix(sweep(x1, 1, sds, "/")) 
  x2[is.na(x2)] = 0
  x2
}

prcomp_centeredScaled <- function(x){
  #** scale & center the data
  x2=centerScale(x) 
  
  #** then apply prcomp 
  prcomp_x2=prcomp(x2, center=FALSE, scale=FALSE)
  prcomp_x2$percentage <- round(prcomp_x2$sdev^2/sum(prcomp_x2$sdev^2)* 100, 2)
  prcomp_x2$center=TRUE
  prcomp_x2$scale=TRUE
  
  class(prcomp_x2)="prcomp"
  return(prcomp_x2)
}


Shannon <- function(x)
{
  x[x==0]=NA #0.000001
  return(x*log(x))
}

get_Shannon <- function(indata, inlog=TRUE)
{
  if(inlog){
    res_Shannon=unlist(apply(exp(indata)/100,2,function(t){-1*sum(Shannon(t), na.rm=TRUE)}))
  }else{
    res_Shannon=unlist(apply(indata/100,2,function(t){-1*sum(Shannon(t), na.rm=TRUE)}))
  }
  return(res_Shannon)
}

getNormData_FIND <- function(indata, aslog=TRUE){
  Iselect = which(apply(indata, 1, function(x) sum(x>0))>10)
  indata.norm = indata[Iselect, order(colnames(indata))]
  if(aslog){
    meds = apply(myLog(indata.norm), 2, median, na.rm=T)
    mads = apply(myLog(indata.norm), 2, mad, na.rm=T)
    indata.norm = sweep(myLog(indata.norm), 2, meds,"-")
    indata.norm = sweep(indata.norm, 2, mads,"/")
  }else{
    indata.norm=apply(indata.norm,2,function(v){v[v<=0]=NA; return(v)})
    meds = apply(indata.norm, 2, median, na.rm=T)
    mads = apply(indata.norm, 2, mad, na.rm=T)
    indata.norm = sweep(indata.norm, 2, meds,"-")
    indata.norm = sweep(indata.norm, 2, mads,"/")
  }
  return(indata.norm)
}

########################### DATA -----
chuv = readRDS("CHUV.rds")
find = readRDS("FIND_1stAnd2ndBatch.rds")

## norm FIND data ----

dataFindNorm.l=list()
dataFindNorm.l[["FIND_1stBatch"]]=lapply(find$FIND_1stBatch_data,function(dataFind){
  return(getNormData_FIND(dataFind[,grep("FN", names(dataFind))], aslog=TRUE))
})
dataFindNorm.l[["FIND_2ndBatch"]]=lapply(find$FIND_2ndBatch_data,function(dataFind){
  return(getNormData_FIND(dataFind[,grep("FN", names(dataFind))], aslog=TRUE))
})


## norm chuv data ----
dataChuvNorm.l=list()
for(curLevel in names(chuv$data)){
  dataChuv = as.data.frame(chuv$data[[curLevel]])
  rl = sapply(strsplit(rownames(dataChuv), ";", fixed=T), "[[", 1)
  if (curLevel == "Genus") rl[2] = "Prevotella_2"
  rownames(dataChuv) = rl
  names(dataChuv) = sapply(strsplit(names(dataChuv), "_"), "[[", 2)
  annotChuv = apply(chuv$design, 1, function(x) paste(x[c(3,5,4,6)], collapse='|'))[1:98]
  names(annotChuv) = sapply(strsplit(names(annotChuv), "_"), "[[", 2)
  annotLabels = names(chuv$design)[c(3,5,4)]
  
  Is2 = which(apply(dataChuv, 2, function(x) sum(x>0))>4)
  dataChuv = dataChuv[,Is2]
  dataChuv = dataChuv[,order(names(dataChuv))]
  Iselect = which(apply(dataChuv, 1, function(x) sum(x>0))>5)
  dataChuv.norm = dataChuv[Iselect,]
  
  meds = apply(myLog(dataChuv.norm), 2, median, na.rm=T)
  mads = apply(myLog(dataChuv.norm), 2, mad, na.rm=T)
  dataChuv.norm = sweep(myLog(dataChuv.norm), 2, meds,"-")
  dataChuv.norm = sweep(dataChuv.norm, 2, mads,"/")
  dataChuvNorm.l[[curLevel]]=dataChuv.norm
}


############################ Colors + additional vectors -----
nlev = 11
myColors = rev(colorRampPalette(brewer.pal(nlev,"RdYlBu"))(nlev))
col_greens=brewer.pal(9, 'Greens')[2:9]
colors_Regions=c("lightskyblue","blue","brown","chartreuse3","chocolate2","deeppink")
colors_TB=c("seagreen1","blueviolet")
colorsPhylum <- c("blue","darkgreen",rgb(201/255, 79/255, 77/255),rgb(195/255,214/255,155/255),rgb(223/255,102/255,3/255),"light blue",rgb(255/255, 153/255, 255/255),"pink",rainbow(7))
names(colorsPhylum)=c("Actinobacteria","Bacteroidetes","Firmicutes","Fusobacteria","Proteobacteria","Spirochaetes","SR1","Tenericutes","TM7","Gemmatimonadetes","GN01","Planctomycetes","Cyanobacteria","Others","Others2")
patientColors <- primary.colors(20,step=3, no.white=TRUE)[c(1:6,8:9,12,18:19)]
colorsTP=brewer.pal(9, 'Purples')[3:9] 

myBreaks = seq(-1.5, 1.5, length.out=nlev+1)


# correspondance Phylum with levels 
colPhylum.l=list()
for(curLevel in colnames(find$FIND_2ndBatch_data$Genus[,1:5])){
  levelToPhylum=aggregate(find$FIND_2ndBatch_data[[curLevel]][,"Phylum"], by=list("curlevel"=as.factor(find$FIND_2ndBatch_data[[curLevel]][,curLevel])), unique); 
  colnames(levelToPhylum)=c(curLevel,"Phylum")
  levelToPhylum=data.frame(levelToPhylum, colorsPhylum[levelToPhylum$Phylum])
  rownames(levelToPhylum)=levelToPhylum[,curLevel]; levelToPhylum=levelToPhylum[-1]
  colnames(levelToPhylum)=c("Phylum","colorPhylum")
  colPhylum.l[[curLevel]]=levelToPhylum
}



############################  FIND1 (Figure 1)  ----
annots=find$FIND_1stBatch_annots
curLevel="Phylum"
dataFind.norm=dataFindNorm.l[["FIND_1stBatch"]][[curLevel]]

##Fig1.A: heatmap ----
row.means = apply(dataFind.norm, 1, mean, na.rm=T)
Q.norm = as.matrix(sweep(dataFind.norm, 1, row.means, "-"))
colnames(Q.norm) = paste(annots[colnames(dataFind.norm),"Region"],annots[colnames(dataFind.norm),"TB"],sep=" | ")
rownames(Q.norm) = rownames(dataFind.norm)

Rd = hclust(as.dist(myCor(Q.norm)))
Cd = hclust(as.dist(myCor(t(Q.norm))))

nclust=2; 
cut<-cutree(Cd,nclust); 

toPlot=Q.norm
colnames(toPlot)=paste0(1:ncol(toPlot),":",colnames(Q.norm))
curAnnot=data.frame("Region"=as.factor(sapply(colnames(Q.norm), function(s){unlist(strsplit(s," \\| "))[1]})),
                    "TB"=as.factor(sapply(colnames(Q.norm), function(s){unlist(strsplit(s," \\| "))[2]})),
                    "cluster"=as.factor(cut)
)
rownames(curAnnot)=colnames(toPlot)

p=pheatmap(toPlot, 
           clustering_distance_cols=myCor(t(Q.norm)), clustering_distance_rows=myCor(Q.norm),
           breaks = myBreaks,color=myColors,
           annotation_col = curAnnot,
           annotation_colors = list("Region"=setNames(colors_Regions[1:length(unique(curAnnot$Region))],unique(curAnnot$Region)),
                                    "TB"=setNames(colors_TB,c("N","Y")),
                                    "cluster"=setNames(col_greens[1:nclust], 1:nclust)
           ),
           cellwidth = 25, cellheight =30,
           cutree_cols=nclust,labels_row=rownames(Q.norm), labels_col = colnames(Q.norm),
           cluster_rows=Rd, cluster_cols=Cd,
           main=paste0(curLevel," FIND 1st batch")
)
print(p)

##Fig1.B: PCA -----
dataFind.pca = prcomp_centeredScaled(dataFind.norm)
data1=as.data.frame(dataFind.pca$rotation)
data1$features=annots[colnames(dataFind.norm),"Region"]
data1$TB=annots[colnames(dataFind.norm),"TB"]; data1$TB[data1$TB=="Y"]="TB"; data1$TB[data1$TB=="N"]="non-TB"

percentage <- dataFind.pca$percentage
percentage <- paste( colnames(data1), "(", paste( as.character(percentage), "%", " explained variance)", sep="") )

p <- ggplot(data1,aes(x=PC1,y=PC2,label=features, color=features , shape=TB)) + ggtitle(curLevel) 
p <- p + geom_point(size=4)+ theme_bw() + xlab(percentage[1]) + ylab(percentage[2])
p <- p + scale_shape_manual(values=c(18, 17))
p <- p + scale_color_manual(values=setNames(colors_Regions,levels(data1$features)))
p <- p + theme(panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
               legend.position = "top", legend.title = element_blank(),
               axis.line = element_line(colour = "black"),
               legend.text = element_text(color = "black", size = rel(0.5)))
p <- p + guides(color = guide_legend(override.aes = list(shape = 15)))
print(p)


# inset
data0 <- data.frame("percentExplained"=dataFind.pca$percentage,
                    "PC"=paste0("PC",1:length(dataFind.pca$percentage)))

p0 <- ggplot(data=data0, aes(x=PC, y=percentExplained)) 
p0 <- p0 + geom_bar(stat="identity", width=0.8, fill="grey75") + theme_bw() + xlab("") + ylab("% variance explained")
p0 <- p0 + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
print(p0)


############################  FIND2: DE TB/non-TB ----
annots=find$FIND_2ndBatch_annots
DE_FIND_2ndBatch.l=list()

for(curLevel in names(dataFindNorm.l[["FIND_2ndBatch"]])){
  dataFind.norm=dataFindNorm.l[["FIND_2ndBatch"]][[curLevel]]
  
  conds=factor(annots[colnames(dataFind.norm),"TB"], levels=c("Y","N")); names(conds)=colnames(dataFind.norm)
  design=model.matrix(~ 0+conds)         
  colnames(design) = sub("conds", "", colnames(design))
  rownames(design) = rownames(annots)
  fit <- lmFit(dataFind.norm,design)
  fit <- contrasts.fit(fit, makeContrasts(Y - N, levels = design))
  fit2 <- eBayes(fit)
  tt= topTable(fit2, coef=1, adjust="BH",number=nrow(dataFind.norm), sort.by="B") 
  DE_FIND_2ndBatch.l[[curLevel]]=data.frame("avgExp_TB"=rowMeans(dataFind.norm[rownames(tt),names(conds)[which(conds=="Y")]], na.rm=TRUE),
                                            "avgExp_nonTB"=rowMeans(dataFind.norm[rownames(tt),names(conds)[which(conds=="N")]], na.rm=TRUE),
                                            tt)
  outfile=paste0("DEtaxa_",curLevel,"_final.txt")
  toExport=data.frame("bacteria"=rownames(tt),
                      "Phylum"=colPhylum.l[[curLevel]][rownames(tt),"Phylum"],
                      "avgExp_TB"=rowMeans(dataFind.norm[rownames(tt),names(conds)[which(conds=="Y")]], na.rm=TRUE),
                      "avgExp_nonTB"=rowMeans(dataFind.norm[rownames(tt),names(conds)[which(conds=="N")]], na.rm=TRUE),
                      tt)
  write.table(toExport, file=outfile, sep="\t",row.names=FALSE, col.names=TRUE, quote=FALSE)
  
}



############################  FIND2 (Figure 2)  ----
annots=find$FIND_2ndBatch_annots
curLevel="Family"

dataFind.norm=dataFindNorm.l[["FIND_2ndBatch"]][[curLevel]]


# 2.A: heatmap ----
row.means = apply(dataFind.norm, 1, mean, na.rm=T)
Q.norm = as.matrix(sweep(dataFind.norm, 1, row.means, "-"))
colnames(Q.norm) = paste(annots[colnames(dataFind.norm),"Region"],annots[colnames(dataFind.norm),"TB"],sep=" | ")
rownames(Q.norm) = rownames(dataFind.norm)

colnames_Sequencing=annots[colnames(dataFind.norm),"Sequencing"]

Rd = hclust(as.dist(myCor(Q.norm)))
Cd = hclust(as.dist(myCor(t(Q.norm))))

nclust=2; 
cut<-cutree(Cd,nclust); 

toPlot=Q.norm 
colnames(toPlot)=paste0(1:ncol(toPlot),":",colnames(Q.norm))
curAnnot=data.frame("Region"=as.factor(sapply(colnames(Q.norm), function(s){unlist(strsplit(s," \\| "))[1]})),
                    "TB"=as.factor(sapply(colnames(Q.norm), function(s){unlist(strsplit(s," \\| "))[2]})),
                    "Sequencing"=as.factor(colnames_Sequencing),
                    "cluster"=as.factor(cut)
                    
)
rownames(curAnnot)=colnames(toPlot)
annotsColors=annotation_colors = list("Region"=setNames(colors_Regions[1:length(unique(curAnnot$Region))],unique(curAnnot$Region)),
                                      "TB"=setNames(colors_TB,c("N","Y")),
                                      "cluster"=setNames(col_greens[1:nclust], 1:nclust),
                                      "Sequencing"=setNames(brewer.pal(n=8,name = "Set1")[-3][1:length(unique(colnames_Sequencing))],levels(curAnnot$Sequencing))
)


isSign=rep("none",nrow(toPlot)); names(isSign)=rownames(toPlot)
isSign[rownames(toPlot)[which(DE_FIND_2ndBatch.l[[curLevel]][rownames(toPlot),"logFC"]>1 & DE_FIND_2ndBatch.l[[curLevel]][rownames(toPlot),"adj.P.Val"]<0.05)]]="up"
isSign[rownames(toPlot)[which(DE_FIND_2ndBatch.l[[curLevel]][rownames(toPlot),"logFC"]< -1 & DE_FIND_2ndBatch.l[[curLevel]][rownames(toPlot),"adj.P.Val"]<0.05)]]="down"


curRowAnnots=data.frame(
  #"upDown"=ifelse(DE_FIND_2ndBatch.l[[curLevel]][rownames(toPlot),"logFC"]>0,"up","down"),
  "isSign"=isSign,
  "log2FC"=DE_FIND_2ndBatch.l[[curLevel]][rownames(toPlot),"logFC"],
  "log10_padj"=-1*log10(DE_FIND_2ndBatch.l[[curLevel]][rownames(toPlot),"adj.P.Val"]),
  "Phylum"=as.factor(colPhylum.l[[curLevel]][rownames(toPlot),"Phylum"])
)


rownames(curRowAnnots)=rownames(toPlot)

p=pheatmap(toPlot, 
           clustering_distance_cols=myCor(t(Q.norm)), clustering_distance_rows=myCor(Q.norm),
           breaks = myBreaks,color=myColors,
           annotation_col = curAnnot,
           annotation_colors = list("Region"=setNames(colors_Regions[1:length(unique(curAnnot$Region))],unique(curAnnot$Region)),
                                    "TB"=setNames(colors_TB,c("N","Y")),
                                    "cluster"=setNames(col_greens[1:nclust], 1:nclust),
                                    "Sequencing"=setNames(brewer.pal(n=8,name = "Set1")[-3][1:length(unique(colnames_Sequencing))],levels(curAnnot$Sequencing)),
                                    "Phylum"=colorsPhylum,
                                    "isSign"=c(up="red",down="green",none="grey89"),
                                    "log2FC"=rev(brewer.pal(9,"RdYlBu")[c(1,5,9)])
           ),
           annotation_row = curRowAnnots,
           cutree_cols=2,labels_row=rownames(Q.norm), labels_col = colnames(Q.norm),
           cluster_rows=Rd, cluster_cols=Cd,#fontsize=5,
           main=paste0(curLevel," FIND 2nd batch\n")
)
print(p)

# 2.B: barplots TB/nonTB ----

cutTB=sapply(names(cut),function(s){unlist(strsplit(s," \\| "))[2]})
pieTB=matrix(table(interaction(cutTB,cut)), nrow=2)
colnames(pieTB)=paste0("cluster_",levels(factor(cut)),"\n(",table(cut),")")
rownames(pieTB)=c("non-TB","TB")
pieTB.pc=100*pieTB/colSums(pieTB)

barplot(pieTB.pc, ylab="%",yaxt='n',ylim=c(0,120),main=paste0(curLevel,"\nRepresentation of cluster per TB"), col=colors_TB, border="light grey")
axis(2,at=seq(0,100,by=20))
legend("top",legend=rownames(pieTB), fill=colors_TB, cex=0.8, ncol=2)

# 2.B: barplots Regions ----
cutRegion=sapply(names(cut),function(s){unlist(strsplit(s," \\| "))[1]})
pieRegions=matrix(table(interaction(cutRegion,cut)), nrow=2)
colnames(pieRegions)=paste0("cluster_",1:nclust,"\n(",table(cut),")")
rownames(pieRegions)=names(table(cutRegion))
pieRegions.pc=100*(pieRegions/colSums(pieRegions))

barplot(pieRegions.pc, ylab="%",yaxt='n',ylim=c(0,120),main=paste0(curLevel,"\nRepresentation of Region per cluster"), col=colors_Regions, border="light grey")
axis(2,at=seq(0,100,by=20))
legend("top",legend=levels(as.factor(names(table(cutRegion)))), fill=colors_Regions[1:length(levels(as.factor(names(table(cutRegion)))))], cex=0.8, ncol=2)


# 2.C: boxplot ----
dataDE=DE_FIND_2ndBatch.l[[curLevel]]  
dataDE=dataDE[order(dataDE$logFC),]
I <- which(abs(dataDE$logFC)>1 & dataDE$adj.P.Val<0.05)
dataToPlot=as.matrix(dataFind.norm[rownames(dataDE)[I],])

dataToPlot.forBoxplot = data.frame("abundance"=as.numeric(dataToPlot),
                                   "bacteria"=rep(rownames(dataToPlot), time=ncol(dataToPlot)),
                                   "sample"=rep(colnames(dataToPlot), each=nrow(dataToPlot)),
                                   "TB"= annots[rep(colnames(dataToPlot), each=nrow(dataToPlot)),"TB"]
)

isSign=data.frame("bacteria"=rownames(dataDE)[I],
                  "abundance"=apply(dataToPlot, 1, function(v){max(v,na.rm=TRUE)}),
                  "upDown"=ifelse(dataDE[I,"logFC"]>0,"up","down"))

pp <- ggplot(dataToPlot.forBoxplot, aes(x=bacteria, y=abundance, fill=TB)) + xlab(curLevel) + ylab("norm abundance")
pp <- pp + stat_boxplot(geom ='errorbar') + ggtitle(paste0("Representative phyla\n",curLevel))
pp <- pp + geom_boxplot(color="grey50",varwidth=FALSE,outlier.alpha = 0.1) 
pp <- pp + scale_fill_manual(values=setNames(colors_TB,c("Y","N")))
pp <- pp + theme_bw() + theme(axis.text.x = element_text(angle = 90, hjust = 1),panel.grid.major = element_blank(), panel.grid.minor = element_blank())

print(pp)


############################  FIND2 (Suppl. Figure 2)  ----
# S2.A: PCA + boxplots ----
dataFind.pca = prcomp_centeredScaled(dataFind.norm)

## inset: barplot %variance explained
data0 <- data.frame("percentExplained"=dataFind.pca$percentage)
data0 <- data.frame(data0,"PC"=factor(paste0("PC",1:length(dataFind.pca$percentage)), levels=paste0("PC",1:length(dataFind.pca$percentage))))
p0 <- ggplot(data=data0[1:10,], aes(x=PC, y=percentExplained)) +
  geom_bar(stat="identity", width=0.5) + theme(axis.text.x = element_text(angle = 90, hjust = 1))

print(p0)

## PCA (colored by sequencing, shape by country, filled by TB)
data1=as.data.frame(dataFind.pca$rotation)
data1$features=annots[colnames(dataFind.norm),"Region"]
data1$TB=annots[colnames(dataFind.norm),"TB"]; data1$TB[data1$TB=="Y"]="TB"; data1$TB[data1$TB=="N"]="non-TB"
data1$TBColors=colors_TB[as.numeric(as.factor(data1$TB))]
data1$Sequencing=annots[colnames(dataFind.norm),"Sequencing"]

p1 <- ggplot(data1,aes(x=PC1,y=PC2,shape=interaction(features,TB), color=Sequencing )) + ggtitle(curLevel) 
p1 <- p1 + geom_point(size=2)+ theme_bw() + xlab(percentage[1]) + ylab(percentage[2])
p1 <- p1 + scale_shape_manual(values=c(1, 2, 16, 17))
p1 <- p1 + scale_color_manual(values=brewer.pal(n=8,name = "Set1")[-3])
p1 <- p1 + theme(legend.position = "top", legend.title = element_blank(),
                   legend.text = element_text(color = "black", size = 8))

print(p1)



## boxplots
data11 = data.frame("value"=as.numeric(as.matrix(data1[,1:5])),
                    "PC"=rep(colnames(data1)[1:5], each=nrow(data1)),
                    "TB"=as.character(as.matrix(data1[,"TB"])),
                    "features"=as.character(as.matrix(data1[,"features"])),
                    "Sequencing"=as.character(as.matrix(data1[,"Sequencing"])))

p2 <- ggplot(data11, aes(x=PC, y=value, fill=features),label=rownames(data11)) + 
  stat_boxplot(geom ='errorbar') +
  geom_boxplot(color="grey50",varwidth=TRUE,outlier.alpha = 0.1) + scale_fill_manual(values=setNames(colors_Regions, levels(data1$features)))
print(p2)

p3 <- ggplot(data11, aes(x=PC, y=value, fill=TB),label=rownames(data11)) + 
  stat_boxplot(geom ='errorbar') +
  geom_boxplot(color="grey50",varwidth=FALSE,outlier.alpha = 0.1) + scale_fill_manual(values=colors_TB)
print(p3)

p4 <- ggplot(data11, aes(x=PC, y=value, fill=Sequencing),label=rownames(data11)) + 
  stat_boxplot(geom ='errorbar') +
  geom_boxplot(color="grey50",varwidth=FALSE,outlier.alpha = 0.1) + scale_fill_manual(values=brewer.pal(n=8,name = "Set1")[-3])
print(p4)

#S2. B: PC1 vs. PC2 bacteria, colored by Phylum + boxplots PC1 and PC2

data2 <- as.data.frame(dataFind.pca$x)
data2$feature <- row.names(data2)
data2$Phylum <- colPhylum.l[[curLevel]][rownames(data2),"Phylum"]
data2$PhylumColors <- colPhylum.l[[curLevel]][rownames(data2),"colorPhylum"]


q1 <- ggplot(data2,aes(x=PC1,y=PC2,label=feature,color=Phylum )) + ggtitle(curLevel) 
q1 <- q1 + scale_color_manual(values=colorsPhylum)
q1 <- q1 + geom_point()+ theme_bw() + geom_text_repel(size=3,nudge_y=0.05) + xlab(percentage[1]) + ylab(percentage[2])
q1 <- q1 + expand_limits(x=c(-1,1), y=c(-1, 1)) + theme(legend.position = "none")

print(q1)

## boxplot, PC scores per Bacteria (at Phylum level) 
curPC="PC1"
data21 <- data.frame("Family"=rownames(data2),
                     "Phylum"=colPhylum.l[[curLevel]][rownames(data2),"Phylum"],
                     "value"=data2[,curPC])
q2 <- ggplot(data21, aes(x=Phylum, y=value, fill=Phylum),label=rownames(data2)) + xlab(curPC) +
  stat_boxplot(geom ='errorbar') + theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  geom_boxplot(color="grey50",varwidth=TRUE,outlier.alpha = 0.1) + scale_fill_manual(values=colorsPhylum)
print(q2)

curPC="PC2"
data22 <- data.frame("Family"=rownames(data2),
                     "Phylum"=colPhylum.l[[curLevel]][rownames(data2),"Phylum"],
                     "value"=data2[,curPC])
q3 <- ggplot(data22, aes(x=Phylum, y=value, fill=Phylum),label=rownames(data22)) + xlab(curPC) +
  stat_boxplot(geom ='errorbar') + theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  geom_boxplot(color="grey50",varwidth=TRUE,outlier.alpha = 0.1) + scale_fill_manual(values=colorsPhylum)

print(q3)



############### Suppl. Fig1.: Shannon FIND1 and FIND2 -----

pdf("boxplots_Shannon_1stAnd2ndFIND_norm_allInOne.pdf", width=16, height=14)
par(mfrow=c(2,1))
annots=find$FIND_1stBatch_annots
H_all=do.call(rbind,lapply(dataFindNorm.l[["FIND_1stBatch"]], function(dd){get_Shannon(dd, inlog=TRUE)}))
toPlot=data.frame("Shannon"=as.numeric(H_all),
               "level"=factor(rep(rownames(H_all), times=ncol(H_all)), levels=c("Genus","Family","Order","Class","Phylum")), 
               "TB"=rep(annots[colnames(H_all),"TB"],each=nrow(H_all)), 
               "sample"=rep(colnames(H_all), each=nrow(H_all)))
res_ttests.l[["FIND_1stBacth"]]=apply(H_all,1, function(v){t.test(v~annots[names(v),"TB"], alternative="two.sided")$p.value})
res_paired.l[["FIND_1stBacth"]]=pairwise.t.test(toPlot$Shannon, g=interaction(toPlot$TB,toPlot$level))$p.value
b=boxplot(toPlot$Shannon ~ interaction(toPlot$TB,toPlot$level), las=3, ylim=c(0,2.1), pch=16, cex=0.8, 
          col=colors_TB, varwidth=TRUE, xaxt='n',
          ylab="Shannon index",main=paste0("FIND 1st batch"))
axis(1,at=seq(1.5,(2*nrow(H_all)),by=2), labels=levels(toPlot$level))
mypos=c(0.8,0.8,1.2,1.2,0.2)
mypos_txt=c(0.7,0.7,1.3,1.3,0.1)
for(i in 1:nrow(H_all)){
  j=seq(1,10,by=2)[i]
  lines( x=c(j,j+1), y=rep(mypos[i],2), lwd=2, lty=3)
  text(j+0.5, mypos_txt[i], cex=0.8, label=sprintf("%.2f",res_ttests.l[["FIND_1stBacth"]][unlist(strsplit(b$names[j],"\\."))[2]]))
}
legend("topright", fill=colors_TB, legend=c("non-TB","TB"), bty='n', cex=0.8)
annots=find$FIND_2ndBatch_annots
H_all=do.call(rbind,lapply(dataFindNorm.l[["FIND_2ndBatch"]], function(dd){get_Shannon(dd, inlog=TRUE)}))
toPlot=data.frame("Shannon"=as.numeric(H_all),
               "level"=factor(rep(rownames(H_all), times=ncol(H_all)), levels=c("Genus","Family","Order","Class","Phylum")), 
               "TB"=rep(annots[colnames(H_all),"TB"],each=nrow(H_all)), 
               "sample"=rep(colnames(H_all), each=nrow(H_all)))
res_ttests.l[["FIND_2ndBatch"]]=apply(H_all,1, function(v){t.test(v~annots[names(v),"TB"], alternative="two.sided")$p.value})
res_paired.l[["FIND_2ndBatch"]]=pairwise.t.test(toPlot$Shannon, g=interaction(toPlot$TB,toPlot$level))$p.value
b=boxplot(toPlot$Shannon ~ interaction(toPlot$TB,toPlot$level), las=3, ylim=c(0,2.1), pch=16, cex=0.8, 
          col=colors_TB, varwidth=TRUE, xaxt='n',
          ylab="Shannon index",main=paste0("FIND 2nd batch"))
axis(1,at=seq(1.5,(2*nrow(H_all)),by=2), labels=levels(toPlot$level))
mypos=c(1.2,1,0.6,0.6,0.2)
mypos_txt=c(1.1,0.9,0.5,0.5,0.1)
for(i in 1:nrow(H_all)){
  j=seq(1,10,by=2)[i]
  lines( x=c(j,j+1), y=rep(mypos[i],2), lwd=2, lty=3)
  text(j+0.5, mypos_txt[i], cex=0.8, label=sprintf("%.2f",res_ttests.l[["FIND_2ndBatch"]][unlist(strsplit(b$names[j],"\\."))[2]]))
}
legend("topright", fill=colors_TB, legend=c("non-TB","TB"), bty='n', cex=0.8)
dev.off()




############################  CHUV (Figure 3)  ----
dataToPlot=dataChuvNorm.l[[curLevel]]

TP=data.frame(1:6,
              "labels_AF"=c("A","B","C","D","E","F"),
              "labels_wm"=c("0w","2w","4w","8w","5m","6m"))
TP_wm=TP$labels_wm; names(TP_wm)=TP$labels_AF

annots=chuv$design; rownames(annots)=annots$sampleName
annots=data.frame(annots, 
                  "TP"=factor(TP_wm[annots[rownames(annots),"timePoint"]],levels=c("0w","2w","4w","8w","5m","6m")))


# Fig3. A: boxplots per time-point ----
LabelSamples=rep(colnames(dataToPlot),each = nrow(dataToPlot), len = ncol(dataToPlot)*nrow(dataToPlot))
LabelNames=rep(rownames(dataToPlot),ncol(dataToPlot))

dataToPlot.forBoxplot <- data.frame(values=as.numeric(as.matrix(dataToPlot)),
                                    LabelSamples=LabelSamples,
                                    LabelNames=LabelNames,
                                    LabelPatient=chuv$design[LabelSamples,"PatientID"],
                                    LabelTP=chuv$design[paste0("pc_",LabelSamples),"timePoint"],
                                    labelPhylum=sapply(LabelNames,function(x){tail(unlist(strsplit(x,";")),1)}),
                                    labelLevels=sapply(LabelNames,function(x){unlist(strsplit(x,";"))[1]}),
                                    stringsAsFactors=TRUE
)
curYmin=floor(min(dataToPlot.forBoxplot[,"values"], na.rm=TRUE))
curYmax=ceiling(max(dataToPlot.forBoxplot[,"values"], na.rm=TRUE))
b=boxplot(dataToPlot.forBoxplot[,"values"]~dataToPlot.forBoxplot[,"LabelTP"]*dataToPlot.forBoxplot[,"LabelNames"],las=3,col=colorsPhylum[dataToPlot.forBoxplot[,"labelPhylum"]],xaxt="n",pch='x',cex=0.5, outline=FALSE,ylim=c(curYmin,curYmax),plot=FALSE)
xlabels=unlist(lapply(b$names,function(x){substr(unlist(strsplit(x,";"))[1],3,nchar(x))}))
vectTP=as.numeric(as.factor(unlist(lapply(b$names,function(x){unlist(strsplit(x,"\\."))[1]}))))

boxplot(dataToPlot.forBoxplot[,"values"]~dataToPlot.forBoxplot[,"LabelTP"]*dataToPlot.forBoxplot[,"LabelNames"],las=3,col=colorsTP[vectTP],xaxt="n",pch=16,cex=0.5, outline=TRUE,ylim=c(curYmin,curYmax),plot=TRUE, ylab="Normalized Relative abundance")
abline(h=0, col="grey33", lwd=0.7, lty=2)
s=(seq(1,length(levels(dataToPlot.forBoxplot[,"LabelNames"]))*length(levels(dataToPlot.forBoxplot[,"LabelTP"])),by=length(levels(dataToPlot.forBoxplot[,"LabelTP"]))))+3
axis(1,labels=xlabels[s],at=s-0.5,las=3,cex.axis=0.8)
s=(seq(1,length(levels(dataToPlot.forBoxplot[,"LabelNames"]))*length(levels(dataToPlot.forBoxplot[,"LabelTP"])),by=length(levels(dataToPlot.forBoxplot[,"LabelTP"]))))+5.5
abline(v=s,lty=2,col="light grey")
abline(v=0,lty=2,col="light grey")
legend("top",legend=names(TP), fill=colorsTP, ncol=6, cex=0.9)
title(curLevel)


# Fig3.B: heatmap ----
dataToPlot=dataChuvNorm.l[[curLevel]]
row.means = apply(dataToPlot, 1, mean, na.rm=T)
Q.norm = as.matrix(sweep(dataToPlot, 1, row.means, "-"))
colnames(Q.norm) = colnames(dataToPlot)

Rd = hclust(as.dist(myCor(Q.norm)))
Cd = hclust(as.dist(myCor(t(Q.norm))))

Qnorm_letters=TP_wm[annots[colnames(Q.norm),"timePoint"]]

curAnnots=data.frame("TP"=factor(Qnorm_letters,levels=c("0w","2w","4w","8w","5m","6m")),
                     "Patient"=factor(annots[colnames(Q.norm),"PatientID"], levels=unique(annots[colnames(Q.norm),"PatientID"]))
); 
rownames(curAnnots)=colnames(Q.norm)

p2=pheatmap(Q.norm, 
            breaks = myBreaks,color=myColors,
            annotation_col = curAnnots,
            annotation_colors = list("TP"=setNames(col_greens[1:length(levels(factor(TP_wm[annots[colnames(Q.norm),"timePoint"]],levels=c("0w","2w","4w","8w","5m","6m"))))],levels(factor(TP_wm[annots[colnames(Q.norm),"timePoint"]],levels=c("0w","2w","4w","8w","5m","6m")))),
                                     "Patient"=setNames( unique(patientColors[annots[colnames(Q.norm),"PatientID"]]), unique(annots[colnames(Q.norm),"PatientID"]))),
            cellwidth = 12, cellheight =12,
            gaps_col = which(!duplicated(annots[colnames(Q.norm),"PatientID"]))-1,
            labels_row=rownames(Q.norm), labels_col = rep("",ncol(Q.norm)),
            cluster_rows=Rd, cluster_cols=FALSE,
            main=paste0(curLevel," CHUV\n(scaled by rows)")
)
print(p2)


### Figure S3: heatmaps mean and median abundances per time-point -----
toPlot_letters=TP_wm[annots[colnames(dataToPlot),"timePoint"]]
curAnnots=data.frame("TP"=factor(toPlot_letters,levels=c("0w","2w","4w","8w","5m","6m")),
                     "Patient"=factor(annots[colnames(dataToPlot),"PatientID"], levels=unique(annots[colnames(dataToPlot),"PatientID"]))
); 
rownames(curAnnots)=colnames(dataToPlot)


dataToPlot_meanPerTP=t(apply(dataToPlot,1,function(v){tapply(v,curAnnots$TP,function(x){mean(x, na.rm=TRUE)})}))
curMax=ceiling(max(dataToPlot_meanPerTP, na.rm=TRUE)); curMin=floor(min(dataToPlot_meanPerTP, na.rm=TRUE))
curLim=max(abs(curMin),curMax);curBreaks=seq(-1*curLim, curLim, length=41) 
p4=pheatmap(dataToPlot_meanPerTP,  cluster_rows=TRUE,cluster_cols=FALSE,
            cellwidth = 14, cellheight =14,
            breaks=curBreaks, color=colorRampPalette(rev(brewer.pal(n = 7, name ="RdYlBu")))(40),
            main=paste0(curLevel," CHUV\nMean per time-point"), silent=FALSE)
print(p4)

dataToPlot_medianPerTP=t(apply(dataToPlot,1,function(v){tapply(v,curAnnots$TP,function(x){median(x, na.rm=TRUE)})}))
curMax=ceiling(max(dataToPlot_medianPerTP, na.rm=TRUE)); curMin=floor(min(dataToPlot_medianPerTP, na.rm=TRUE))
curLim=max(abs(curMin),curMax);curBreaks=seq(-1*curLim, curLim, length=41) 
p5=pheatmap(dataToPlot_medianPerTP,  cluster_rows=TRUE,cluster_cols=FALSE,
            cellwidth = 14, cellheight =14,
            breaks=curBreaks, color=colorRampPalette(rev(brewer.pal(n = 7, name ="RdYlBu")))(40),
            main=paste0(curLevel," CHUV\nMedian per time-point"), silent=FALSE)
print(p5)

