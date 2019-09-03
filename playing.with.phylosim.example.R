library(phytools)
source("gen.phylosim.R")
source("manip.phylosim.R")
set.seed(7996)
##generate and visualize simulation##
phylosim<-gen.phylosim()
plot(phylosim)
##transform to phylo object and visualize (can even prune off extinct lineages!)##
tree<-trans.to.phylo(phylosim)
tree.extant<-trans.to.phylo(prune.to.lins(phylosim,get.extant.lins(phylosim)))
plot(tree$tree)
plot(tree.extant$tree)
##pass phylo object and associated trait values through phenogram function, note that it's identical to the simulation minus##
##anagenetic changes along branches##
phenogram(tree$tree,tree$trait,ftype='off')
phenogram(tree.extant$tree,tree.extant$trait,ftype='off')
##looking at lineage through time plots, noting that lineage accumulation seems to match up with the rate expected given the net##
##diversification rate used to simulate. The intercept doesn't seem to match well, but this is expected with the ascertainment##
##bias inherent with simulating bd trees (I purposefully set the seed to a value where the tree doesn't end up just going##
##extinct)##
##first I will do this simply with the phylo objects##
expec<-c(2,rep(1+(spec.rate-ext.rate)*res,length(test$time.pts)-1))
expec<-cumprod(expec)
ltt(tree$tree);lines(log(expec)~phylosim$time.pts,lty=2)
ltt(tree.extant$tree);lines(log(expec)~phylosim$time.pts,lty=2)
##haven't yet made simple function for getting ltt's from phylosim object directly, but here's a (relatively) easy way to do it##
spec.rate<-phylosim$params["speciation rate"];ext.rate<-phylosim$params["extinction rate"];res<-phylosim$params["resolution"]
bd.times<-matrix(phylosim$nodes[get.bd(phylosim),"time"],ncol=2)
time<-as.vector(bd.times)
indicator<-rep(c("b","d"),each=nrow(bd.times))
bord<-ifelse(indicator=="b",1,-1)
bord<-c(0,1)*rep(bord[order(time)],each=2)
no.lineage<-cumsum(bord)
time<-rep(sort(time),each=2)
no.lineage<-no.lineage[time!=max(phylosim$time.pts)]
time<-time[time!=max(phylosim$time.pts)]
no.lineage<-c(no.lineage,no.lineage[length(no.lineage)])
time<-c(time,max(phylosim$time.pts))
plot(log(no.lineage)~time,type='l')
expec<-c(1,rep(1+(spec.rate-ext.rate)*res,length(test$time.pts)-1))
expec<-cumprod(expec)
lines(log(expec)~test$time.pts,lty=2)
##note that the lineage through time plot for the phylosim object is identical to that of the phylo object, but with a##
##non-zero period where the phylogeny consisted of only 1 lineage##