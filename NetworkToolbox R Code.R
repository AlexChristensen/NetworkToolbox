###NetworkToolbox: Methods and Measures for Brain,
###Cognitive, and Psychometric Network Analysis
###Alexander P. Christensen, University of North Carolina at Greensboro
###Journal of Statistical Software

#install packages
install.packages("NetworkToolbox")

#load packages
library(NetworkToolbox)

setwd("C:/Users/apchrist/Box Sync/CREATIVITY & ARTS LAB/NetworkToolbox--JSS/data")
psy<-neoOpen
cog<-animals
brain<-restOpen
behav<-behavOpen

#missing data example
psy <- as.vector(as.matrix(psy))
rand <- sample(1:length(psy), 10)
psy[rand]<-NA
psy<-matrix(psy,ncol=48)

TMFG(psy,na.data="listwise")

#construct and plot psy network
network<-TMFG(psy)$A
facets<-c(rep("actions",8),rep("aesthetics",8),
          rep("fantasy",8),rep("feelings",8),
          rep("ideas",8),rep("values",8))
qgraph::qgraph(network, layout = "spring", vsize = 4,
       label.prop = 1, groups=facets,
       palette="ggplot2",title="NEO-PI-3")

#replicate bootgen 
networkrel<-bootgen(psy)

qgraph::qgraph(networkrel$bootmat, layout = "spring", vsize = 4,
       label.prop = 1, groups=facets,
       palette="ggplot2",title="NEO-PI-3 bootgen")

repnetworkrel<-bootgen(psy, seeds = networkrel$Seeds)

#plot bootgen
bootgenPlot(networkrel, bootmat = TRUE)

#Brain Networks
#convert to igraph
brainNetList <- convert2igraph(brain, neural = TRUE)

#brain data
brain<-restOpen  ###see <https://drive.google.com/file/d/1ugwi7nRrlHQYuGPzEB4wYzsizFrIMvKR/view?usp=sharing>
behav<-behavOpen

#cpmIV
openCPM <- cpmIV(brain, behav, method = "mean", model = "linear", shen = TRUE)
openCPM$results

#grab masks for figures
write.table(openCPM$negMask,"OpenNegMask.txt",row.names=FALSE,col.names=FALSE)
write.table(openCPM$posMask,"OpenPosMask.txt",row.names=FALSE,col.names=FALSE)


