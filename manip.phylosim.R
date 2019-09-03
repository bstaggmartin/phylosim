#Utility function: just gets the minimum and maximum observed trait values in a phylosim object given two bounding time points.
#This is a simple thing, but since you have to check both the trait and node matrix, as well as the trait values at each
#bounding time point, it's annoying to write out...
get.min.max.trait.vals<-function(phylosim,xlim){
  stopifnot(inherits(phylosim,"phylosim"))
  c(min(phylosim$trait[,which(phylosim$time.pts>=xlim[1]&phylosim$time.pts<=xlim[2])],
        get.trait.vals(max(0,xlim[1]),phylosim),
        get.trait.vals(min(max(phylosim$time.pts),xlim[2]),phylosim),na.rm=T),
    max(phylosim$trait[,which(phylosim$time.pts>=xlim[1]&phylosim$time.pts<=xlim[2])],
        get.trait.vals(max(0,xlim[1]),phylosim),
        get.trait.vals(min(max(phylosim$time.pts),xlim[2]),phylosim),na.rm=T))
}

#A plotting method for the simulations, in the style of a phenogram
plot.phylosim<-function(phylosim,xlim=c(0,1.05*max(phylosim$time.pts)),ylim=get.min.max.trait.vals(phylosim,xlim),
                        ylab="trait",xlab="time",lwd=3,col=rgb(0.25,0.25,0.25,0.5),
                        label=T,prefix="sp.",label.col="black",adj=c(-0.15,0.25),
                        dot=T,pch=19,dot.col="black",node.dot=F,node.pch=19,node.dot.col="gray"){
  stopifnot(inherits(phylosim,"phylosim"))
  plot(0,col="white",xlim=xlim,ylim=ylim,ylab=ylab,xlab=xlab)
  for(i in 1:nrow(phylosim$trait)){
    start.cap.i<-phylosim$nodes[which(phylosim$nodes[,"child.lin"]==i),]
    end.cap.i<-phylosim$nodes[which(phylosim$nodes[,"parent.lin"]==i&is.na(phylosim$nodes[,"child.lin"])),]
    trait.i<-c(start.cap.i["trait"],phylosim$trait[i,!is.na(phylosim$trait[i,])],end.cap.i["trait"])
    time.i<-c(start.cap.i["time"],phylosim$time.pts[!is.na(phylosim$trait[i,])],end.cap.i["time"])
    lines(trait.i~time.i,lwd=lwd,col=col)
  }
  if(node.dot){
    int.nodes<-which(!is.na(phylosim$nodes[,"child.lin"]))
    points(phylosim$nodes[int.nodes,"trait"]~phylosim$nodes[int.nodes,"time"],pch=node.pch,col=node.dot.col)
  }
  tips<-which(is.na(phylosim$nodes[,"child.lin"]))
  if(dot){
    points(phylosim$nodes[tips,"trait"]~phylosim$nodes[tips,"time"],pch=pch,col=dot.col)
  }
  if(label){
    text(paste("sp.",phylosim$nodes[tips,"tip.lab"]),x=phylosim$nodes[tips,"time"],y=phylosim$nodes[tips,"trait"],
         col=label.col,adj=adj)
  }
}

#Get the birth and death nodes of a given lineage (keep in mind lineages are defined differently in phylosim objects vs. phylo
#objects...)
get.bd<-function(phylosim,lin=1:nrow(phylosim$trait)){
  stopifnot(inherits(phylosim,"phylosim"))
  stopifnot(lin%in%1:nrow(phylosim$trait))
  lin.birth<-which(phylosim$nodes[,"child.lin"]%in%lin)
  lin.death<-which(phylosim$nodes[,"parent.lin"]%in%lin&is.na(phylosim$nodes[,"child.lin"]))
  return(cbind(lin.birth,lin.death))
}

#Try vectorizing these functions at some point; just using mapply for now
#Get the most proximal event to a certain time going backwards in time, be it a lineage splitting event, or a time point
#sampled for a trait value
get.bw.prox.event<-function(lin,phylosim,time){
  stopifnot(inherits(phylosim,"phylosim"))
  stopifnot(lin%in%1:nrow(phylosim$trait))
  lin.death<-which(phylosim$nodes[,"parent.lin"]==lin&is.na(phylosim$nodes[,"child.lin"]))
  if(phylosim$nodes[lin.death,"time"]<time){
    bw.prox.event<-rep(NA,5)
    names(bw.prox.event)<-c("parent.lin","child.lin","time","trait","tip.lab")
    return(bw.prox.event)
  }
  bw.trait.time<-floor(time/phylosim$params["resolution"])*phylosim$params["resolution"]
  potential.indices<-c(which(phylosim$nodes[,"child.lin"]==lin&phylosim$nodes[,"time"]<=time),
                       which(phylosim$nodes[,"parent.lin"]==lin&phylosim$nodes[,"time"]<=time))
  if(length(potential.indices)==0|all(phylosim$nodes[potential.indices,"time"]<bw.trait.time)){
    if(is.na(phylosim$trait[lin,as.character(bw.trait.time)])){
      bw.prox.event<-rep(NA,5)
      names(bw.prox.event)<-c("parent.lin","child.lin","time","trait","tip.lab")
      return(bw.prox.event)
    }
    bw.prox.event<-c(lin,lin,bw.trait.time,phylosim$trait[lin,as.character(bw.trait.time)],NA)
    names(bw.prox.event)<-c("parent.lin","child.lin","time","trait","tip.lab")
    return(bw.prox.event)
  }
  bw.prox.event<-phylosim$nodes[max(potential.indices),]
  return(bw.prox.event)
}
get.bw.prox.events<-function(lin,phylosim,time){
  stopifnot(inherits(phylosim,"phylosim"))
  stopifnot(lin%in%1:nrow(phylosim$trait))
  t(mapply(FUN=get.bw.prox.event,lin=lin,time=time,MoreArgs=list(phylosim=phylosim)))
}
#Same thing, but for going forward in time.
get.fw.prox.event<-function(lin,phylosim,time){
  stopifnot(inherits(phylosim,"phylosim"))
  stopifnot(lin%in%1:nrow(phylosim$trait))
  lin.birth<-which(phylosim$nodes[,"child.lin"]==lin)
  if(phylosim$nodes[lin.birth,"time"]>time){
    fw.prox.event<-rep(NA,5)
    names(fw.prox.event)<-c("parent.lin","child.lin","time","trait","tip.lab")
    return(fw.prox.event)
  }
  fw.trait.time<-ceiling(time/phylosim$params["resolution"])*phylosim$params["resolution"]
  potential.indices<-c(which(phylosim$nodes[,"child.lin"]==lin&phylosim$nodes[,"time"]>=time),
                       which(phylosim$nodes[,"parent.lin"]==lin&phylosim$nodes[,"time"]>=time))
  if(length(potential.indices)==0|all(phylosim$nodes[potential.indices,"time"]>fw.trait.time)){
    if(is.na(phylosim$trait[lin,as.character(fw.trait.time)])){
      fw.prox.event<-rep(NA,5)
      names(fw.prox.event)<-c("parent.lin","child.lin","time","trait","tip.lab")
      return(fw.prox.event)
    }
    fw.prox.event<-c(lin,lin,fw.trait.time,phylosim$trait[lin,as.character(fw.trait.time)],NA)
    names(fw.prox.event)<-c("parent.lin","child.lin","time","trait","tip.lab")
    return(fw.prox.event)
  }
  fw.prox.event<-phylosim$nodes[min(potential.indices),]
  return(fw.prox.event)
}
get.fw.prox.events<-function(lin,phylosim,time){
  stopifnot(is.list(phylosim)&names(phylosim)==c("time.pts","trait","nodes","params"))
  stopifnot(lin%in%1:nrow(phylosim$trait))
  t(mapply(FUN=get.fw.prox.event,lin=lin,time=time,MoreArgs=list(phylosim=phylosim)))
}

#Get value of trait at a given time for a given set of lineages
get.trait.vals<-function(time,phylosim,lin=1:nrow(phylosim$trait)){
  stopifnot(inherits(phylosim,"phylosim"))
  stopifnot(length(time)==1|length(time)==nrow(phylosim$trait))
  stopifnot(time>=min(phylosim$time.pts)&time<=max(phylosim$time.pts))
  stopifnot(lin%in%1:nrow(phylosim$trait))
  bw.prox.events<-get.bw.prox.events(lin,phylosim,time)
  bw.prox.time<-bw.prox.events[,"time"]
  bw.prox.trait<-bw.prox.events[,"trait"]
  fw.prox.events<-get.fw.prox.events(lin,phylosim,time)
  fw.prox.time<-fw.prox.events[,"time"]
  fw.prox.trait<-fw.prox.events[,"trait"]
  trait<-(fw.prox.trait-bw.prox.trait)/(fw.prox.time-bw.prox.time)*(time-bw.prox.time)+bw.prox.trait
  if(!all(!is.nan(trait))){
    nan.trait<-which(is.nan(trait))
    nan.time.pts<-which(phylosim$time.pts%in%bw.prox.time[nan.trait])
    trait[nan.trait]<-phylosim$trait[cbind(lin[nan.trait],nan.time.pts)]
  }
  names(trait)<-lin
  trait
}

#Crop the time at which a phylosim object ends
crop.time<-function(phylosim,end.time){
  existing.lins<-unique(phylosim$nodes[which(phylosim$nodes[,"time"]<=end.time&
                                               !is.na(phylosim$nodes[,"child.lin"])),"child.lin"])
  lins.labels<-phylosim$nodes[which(phylosim$nodes[,"parent.lin"]%in%existing.lins&
                                      is.na(phylosim$nodes[,"child.lin"])),"tip.lab"]
  nodes.new<-cbind(existing.lins,NA,end.time,get.trait.vals(end.time,phylosim,existing.lins),lins.labels)
  nodes.new<-rbind(phylosim$nodes,nodes.new)
  trait.new<-phylosim$trait[existing.lins,as.character(seq(0,end.time,phylosim$params["resolution"]))]
  nodes.new<-nodes.new[which(nodes.new[,"time"]<=end.time),]
  rownames(nodes.new)<-NULL
  return(list(time.pts=seq(0,end.time,phylosim$params["resolution"]),trait=trait.new,nodes=nodes.new,
              params=phylosim$params))
}

#Get all lineages extant at a given time (defaults to the maximum time)
get.extant.lins<-function(phylosim,time=max(phylosim$time.pts)){
  stopifnot(inherits(phylosim,"phylosim"))
  born.lins<-unique(phylosim$nodes[which(phylosim$nodes[,"time"]<time&
                                           !is.na(phylosim$nodes[,"child.lin"])),"child.lin"])
  died.lins<-unique(phylosim$nodes[which(phylosim$nodes[,"time"]<time&
                                           is.na(phylosim$nodes[,"child.lin"])),"parent.lin"])
  return(born.lins[!(born.lins%in%died.lins)])
}

#Get the most recent common ancestor of two lineages (returns the point at which one lineage branched from the other should
#there be a parent-child relationship between the two lineages)
get.MRCA<-function(nodes,lin1,lin2){
  roots1<-NULL
  lin1.death<-which(nodes[,"parent.lin"]==lin1&is.na(nodes[,"child.lin"]))
  roots1[1]<-max(which(nodes[,"child.lin"]==lin1),
                 which(nodes[,"parent.lin"]==lin1&nodes[,"time"]<nodes[lin1.death,"time"]))
  while(!(nodes[roots1[length(roots1)],"parent.lin"]%in%c(NA,lin2))&
        (nodes[roots1[length(roots1)],"child.lin"]!=lin2)){
    roots1[length(roots1)+1]<-max(which(nodes[,"child.lin"]==nodes[roots1[length(roots1)],"parent.lin"]),
                                  which(nodes[,"parent.lin"]==nodes[roots1[length(roots1)],"parent.lin"]&
                                          nodes[,"time"]<nodes[roots1[length(roots1)],"time"]))
  }
  roots2<-NULL
  lin2.death<-which(nodes[,"parent.lin"]==lin2&is.na(nodes[,"child.lin"]))
  roots2[1]<-max(which(nodes[,"child.lin"]==lin2),
                 which(nodes[,"parent.lin"]==lin2&nodes[,"time"]<nodes[lin2.death,"time"]))
  while(!(nodes[roots2[length(roots2)],"parent.lin"]%in%c(NA,lin1))&
        (nodes[roots2[length(roots2)],"child.lin"]!=lin1)){
    roots2[length(roots2)+1]<-max(which(nodes[,"child.lin"]==nodes[roots2[length(roots2)],"parent.lin"]),
                                  which(nodes[,"parent.lin"]==nodes[roots2[length(roots2)],"parent.lin"]&
                                          nodes[,"time"]<nodes[roots2[length(roots2)],"time"]))
  }
  max(roots1[roots1%in%roots2])
}

#get the phylogentic distance between two lineages given a time from which to calculate the distance
get.phylo.dist<-function(nodes,time,lin1,lin2){
  MRCA<-get.MRCA(nodes,lin1,lin2)
  lin1.death<-which(nodes[,"parent.lin"]==lin1&is.na(nodes[,"child.lin"]))
  lin2.death<-which(nodes[,"parent.lin"]==lin2&is.na(nodes[,"child.lin"]))
  time1<-min(time,nodes[lin1.death,"time"])-nodes[MRCA,"time"]
  time2<-min(time,nodes[lin2.death,"time"])-nodes[MRCA,"time"]
  time1+time2
}

#get the total branch length of a phylosim tree
get.tot.bl<-function(phylosim){
  sum(phylosim$nodes[get.bd(phylosim)[,2],"time"]-phylosim$nodes[get.bd(phylosim)[,1],"time"])
}

#these might be able to be sped up...in particular, see if you can use the while loops in 'get.MRCA' to make
#'prune.to.lins' more efficient
#Prunes a phylosim object to a set of lineages, merging kept lineages with cut-out lineages as necessary (the inverse of
#drop.tip, but with phylosim objects)
prune.to.lins<-function(phylosim,lins.to.keep){
  stopifnot(inherits(phylosim,"phylosim"))
  trait.new<-phylosim$trait
  nodes.new<-phylosim$nodes
  res<-phylosim$time.pts[2]-phylosim$time.pts[1]
  for(i in lins.to.keep){
    #initialize a vector of row numbers in node matrix corresponding to "roots" of that extant lineage
    roots<-NULL
    #first entry is most recent, proximal root
    roots<-which(nodes.new[,"child.lin"]==i)
    #get roots further and further back in phylogeny, only stopping when the other child of a given root (i.e., the
    #parental lineage, which always continues on after the node) is either non-existent (i.e., the ultimate root of the
    #entire tree) or another extant lineage
    while(!(nodes.new[roots[length(roots)],"parent.lin"]%in%c(NA,lins.to.keep))){
      roots[length(roots)+1]<-which(nodes.new[,"child.lin"]==nodes.new[roots[length(roots)],"parent.lin"])
    }
    #for all identified roots of extant lineage...
    for(j in roots){
      #identify what iteration of loop you're on
      iter<-which(roots==j)
      #get the parental lineage of the root (j)
      root.lin<-nodes.new[j,"parent.lin"]
      #if this isn't the last iteration...
      if(iter<length(roots)){
        #get the time point of this root (j)
        t2<-floor(nodes.new[j,"time"]/res)*res
        #get the time point of the next root back
        t1<-ceiling(nodes.new[roots[iter+1],"time"]/res)*res
        if(t2>=t1){
          #replace the trait values of the current extant lineage (i) with those of root's (j) parental lineage between
          #this root's (j) time point and the time point of the next root
          trait.new[i,as.character(seq(t1,t2,by=res))]<-trait.new[root.lin,as.character(seq(t1,t2,by=res))]
        }
      }
      #if the parental lineage of the root (j) isn't an extant lineage or non-existent...
      if(!(root.lin%in%c(NA,lins.to.keep))){
        #identify where the lineage occurs in node matrix as a child lineage, before or at time point of root (j)
        #(you don't want to replace the parental lineage after this time point, because you are recoding as a child
        #of the extant lineage, just in case this extinct lineage happens to give birth to any extant lineages down
        #the rood)
        child.lins.to.repl<-which(nodes.new[,"child.lin"]==nodes.new[j,"parent.lin"]&
                                    nodes.new[,"time"]<=nodes.new[j,"time"])
        #identify where the lineage occurs in node matrix as a parental lineage, before or at time point of root (j)
        #("")
        parent.lins.to.repl<-which(nodes.new[,"parent.lin"]==nodes.new[j,"parent.lin"]&
                                     nodes.new[,"time"]<=nodes.new[j,"time"])
        #replace identified child and parental occurrences with current extant lineage (i) instead
        nodes.new[child.lins.to.repl,"child.lin"]<-i
        nodes.new[parent.lins.to.repl,"parent.lin"]<-i
        #make the child lineage of the current root (j) the originally parental lineage of the root (j) instead
        nodes.new[j,"child.lin"]<-root.lin
      }
    }
  }
  #make new node matrix include only nodes that have both extant (or non-existent/NA) parental AND child lineages
  nodes.new<-nodes.new[which((nodes.new[,"child.lin"]%in%c(NA,lins.to.keep))&
                               (nodes.new[,"parent.lin"]%in%c(NA,lins.to.keep))),]
  #make vector of unique lineage numbers in new node matrix, in order of when they branched off into one another
  old.lin.nos<-unique(as.vector(t(nodes.new[,c("parent.lin","child.lin")])))
  trait.new<-trait.new[old.lin.nos[2:length(old.lin.nos)],]
  #make vector of new lineage numbers (making sure NA still correspond to NA)
  new.lin.nos<-c(NA,1:(length(old.lin.nos)-1))
  #name the elements of the old lineage number vector by their new lineage numbers
  names(new.lin.nos)<-old.lin.nos
  #use vector as "look-up table" to rename lineages in new node matrix (makes sure that lineages are numbered by order
  #in which they branch off into one another, making later transformations of node matrix easier)
  nodes.new[,"parent.lin"]<-new.lin.nos[as.character(nodes.new[,"parent.lin"])]
  nodes.new[,"child.lin"]<-new.lin.nos[as.character(nodes.new[,"child.lin"])]
  stopifnot(sum(is.na(nodes.new[,"child.lin"]))==sum(!is.na(nodes.new[,"child.lin"])))
  out<-list(time.pts=phylosim$time.pts,trait=trait.new,nodes=nodes.new,params=phylosim$params)
  class(out)<-"phylosim"
  out
}

#transform a phylosim object into a phylo object
trans.to.phylo<-function(phylosim){
  stopifnot(inherits(phylosim,"phylosim"))
  stopifnot(nrow(phylosim$trait)>1)
  nodes.new<-as.data.frame(phylosim$nodes[order(phylosim$nodes[,"child.lin"],phylosim$nodes[,"time"],na.last=F),])
  int.nodes<-nodes.new[!is.na(nodes.new$child.lin),]
  nodes.new[!is.na(nodes.new$child.lin),]<-int.nodes[order(int.nodes$time),]
  nodes.new<-nodes.new[!is.na(nodes.new$parent.lin),]
  nodes.new$node.num<-1:nrow(nodes.new)
  trait.new<-nodes.new$trait
  names(trait.new)<-c(paste("sp.",nodes.new$tip.lab[!is.na(nodes.new$tip.lab)]),
                      (nrow(phylosim$trait)+1):(2*nrow(phylosim$trait)-1))
  nodes.new<-nodes.new[order(nodes.new$time),]
  edges<-matrix(0,ncol=3,nrow=nrow(nodes.new)-1)
  counter<-1
  for(i in 1:nrow(nodes.new)){
    if(!is.na(nodes.new[i,"child.lin"])){
      child.lin1<-nodes.new[i,"parent.lin"]
      child.lin2<-nodes.new[i,"child.lin"]
      child.node1.rows<-which(nodes.new$parent.lin==child.lin1)
      child.node2.rows<-which(nodes.new$parent.lin==child.lin2)
      child.node1.row<-min(child.node1.rows[which(child.node1.rows>i)])
      child.node2.row<-min(child.node2.rows[which(child.node2.rows>i)])
      edges[counter,]<-c(nodes.new[i,"node.num"],nodes.new[child.node1.row,"node.num"],
                         nodes.new[child.node1.row,"time"]-nodes.new[i,"time"])
      counter<-counter+1
      edges[counter,]<-c(nodes.new[i,"node.num"],nodes.new[child.node2.row,"node.num"],
                         nodes.new[child.node2.row,"time"]-nodes.new[i,"time"])
      counter<-counter+1
    }
  }
  tree<-list(edge=edges[,1:2],tip.label=paste("sp.",nodes.new$tip.lab[!is.na(nodes.new$tip.lab)]),
             Nnode=nrow(phylosim$trait)-1,edge.length=edges[,3],
             node.label=(nrow(phylosim$trait)+1):(2*nrow(phylosim$trait)-1),root.edge=phylosim$nodes[2,"time"])
  class(tree)<-"phylo"
  return(list(tree=tree,trait=trait.new))
}