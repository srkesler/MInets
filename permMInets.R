library(fabricatr)
library(minet)
library(igraph)
library(plotrix)
library(circlize)

#Shelli Kesler 4/22/19
#Compare group level mutual information networks 

# enter parameters here----------------------------------------
Nperm = 5000
dat1 = imp #matrix of cases x vars
dat2 = nimp
#--------------------------------------------------------------
nvars = ncol(dat1)
set.seed(Nperm)

#mutual information networks
net1 = minet(dat1, method = "aracne")
net2 = minet(dat2, method = "aracne")

#graphs
G1 = graph_from_adjacency_matrix(net1,mode="undirected",weighted=T)
G2 = graph_from_adjacency_matrix(net2,mode="undirected",weighted=T)

#absolute group differences in graph metrics
gDiffs = rep(NA,4)
gDiffs[1] = abs(mean(strength(G1, vids = V(G1))) - mean(strength(G2, vids = V(G2))))
eV1 = eigen_centrality(G1,scale=F); eV2 = eigen_centrality(G2,scale=F)
gDiffs[2] = abs(mean(eV1$vector - eV2$vector))
gDiffs[3] = abs(mean(betweenness(G1, v = V(G1))) - mean(betweenness(G2, v = V(G2))))
gDiffs[4] = abs(modularity(cluster_walktrap(G1)) - modularity(cluster_walktrap(G2)))

#differences in bootstrapped networks
x1 <- array(rep(NA, nvars*nvars*Nperm), dim=c(nvars, nvars, Nperm))
x2 <- array(rep(NA, nvars*nvars*Nperm), dim=c(nvars, nvars, Nperm))

bDiffs = matrix(data=NA, nrow=4,ncol=Nperm)

for (i in 1:Nperm) {
x1[,,i] <- minet(resample_data(dat1),method="aracne")
x2[,,i] <- minet(resample_data(dat2),method="aracne")

xG1 = graph_from_adjacency_matrix(x1[,,i],mode="undirected",weighted=T)
xG2 = graph_from_adjacency_matrix(x2[,,i],mode="undirected",weighted=T)

bDiffs[1,i] = abs(mean(strength(xG1, vids = V(xG1))) - mean(strength(xG2, vids = V(xG2))))
bV1 = eigen_centrality(xG1,scale=F); bV2 = eigen_centrality(xG2,scale=F)
bDiffs[2,i] = abs(mean(bV1$vector - bV2$vector))
bDiffs[3,i] = abs(mean(betweenness(xG1, v = V(xG1))) - mean(betweenness(xG2, v = V(xG2))))
bDiffs[4,i] = abs(modularity(cluster_walktrap(xG1)) - modularity(cluster_walktrap(xG2)))
}

#p values and 95% CIs
permP = rep(NA,4)
permP[1] <- mean(bDiffs[1,] > gDiffs[1])
permP[2] <- mean(bDiffs[2,] > gDiffs[2])
permP[3] <- mean(bDiffs[3,] > gDiffs[3])
permP[4] <- mean(bDiffs[4,] > gDiffs[4])

permCI = matrix(data=NA,nrow=4,ncol=2)
lolim = .025*Nperm
hilim = Nperm-lolim
for (i in 1:4){
  y = sort(bDiffs[i,],decreasing = F)
  permCI[i,1] = y[lolim]
  permCI[i,2] = y[hilim]
}

#results
labs = c('strength','eigen centrality','node betweenness','modularity')
res = matrix(data = NA, nrow=4,ncol=4)
rownames(res) <- labs
colnames(res) <- c('MeanDiff','95%CI LL','95%CI UL', 'pval')
res[,1] = t(gDiffs)
res[,2] = permCI[,1]
res[,3] = permCI[,2]
res[,4] = t(permP)
res
save(res,net1,net2,x1,x2,gDiffs,bDiffs,file="permMIdat.RData")
write.csv(res,"permMIresults.csv")

#plot
x = 1:4
yticks = round(seq(min(permCI[,1]),max(permCI[,2]),.2),digits=3)
plotCI(x, gDiffs, ui=permCI[,2], li=permCI[,1], gap=0,lwd=5, pch=25,
       ylab='Permutation Distribution Confidence Interval',xlab='',axes=FALSE,cex.lab=1.5,cex=2) 
axis(1, at=1:length(x),font.axis=2,cex.axis=1.5,las = 0,labels=labs) 
axis(2, at=yticks,font.axis=2,cex.axis=1.5)
box(lty = 'solid', col = 'black')
