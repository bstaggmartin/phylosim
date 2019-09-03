#function for generating birth-death phylogeny with constantly-sampled continuous trait evolving via BM/OU process
gen.phylosim<-function(end.time=500,res=1,spec.rate=0.01,ext.rate=0.005,
                       rate=1,init.trait=0,trend=0,alpha=0,theta=0){
  #helpful function for returning integers between x1 and x2 (where x1 and x2 are equal length components of a vector
  #c(x1,x2), discarding any cases x1[i]==x2[i] or x[i]>x2[i])
  mult.seq<-function(x){
    stopifnot(length(x)%%2==0)
    output<-NULL
    for(i in 1:(length(x)/2)){
      if(x[length(x)/2+i]>=x[i]){
        output<-c(output,seq(x[i],x[length(x)/2+i]))
      }
    }
    return(output)
  }
  #generate vector of time points
  time.vec<-seq(0,end.time,by=res)
  #initialize matrix of trait values over time, having NA for non-existent/not yet sampled time point-lineage combos
  trait<-matrix(NA,nrow=1,ncol=length(time.vec))
  #name columns of the matrix according to time points
  colnames(trait)<-time.vec
  #initialize starting trait value
  trait[1,"0"]<-init.trait
  #initialize node matrix, with first node being the root of the stem
  nodes<-cbind(NA,1,0,trait[1,"0"]);colnames(nodes)<-c("parent.lin","child.lin","time","trait")
  #initialize loop for generating tree; for each time point...
  for(t in time.vec[1:(length(time.vec))-1]){
    #simulate trait values at the next time point; trait values at current time point plus...
    trait[,as.character(t+res)]<-trait[,as.character(t)]+
      #OU component plus...
      alpha*(theta-trait[,as.character(t)])+
      #BM component
      rnorm(nrow(trait),trend,rate*res)
    #define vector for the minimum start time of events for each lineage (-the current time point, t)
    event.starts<-rep(0,nrow(trait))
    #define vector of trait values for each lineage at the minimum start time for events
    event.starts.trait.vals<-trait[,as.character(t)]
    #simulate number of events occurring for each lineage by drawing from poisson distribution w/ lambda equal to
    #the sum of speciation and extinction rates (scaled by the resolution of the time steps). Store these as a vector
    #where each element corresponds to a lineage
    lin.events<-rpois(nrow(trait),(spec.rate+ext.rate)*res)
    #set the number of events for non-existent lineages (e.g., have gone extinct) to 0
    lin.events[which(is.na(trait[,as.character(t)]))]<-0
    #initialize loop for simulating additional lineage events and grafting on new lineages/ending old ones; while the 
    #total number of lineage events occurring across the phylogeny is greater than 0...
    while(sum(lin.events)>0){
      #create a vector of parent lineages for each event by replicating parent lineage numbers according to the
      #number of lineage events
      event.parent.lins<-rep(1:nrow(trait),lin.events)
      #do the same for event starting times
      event.starts<-rep(event.starts,lin.events)
      #again for trait values at event starting times...
      event.starts.trait.vals<-rep(event.starts.trait.vals,lin.events)
      #pick random times between start times and next time point (-the current time point, t)
      event.times<-runif(sum(lin.events),event.starts,res)
      #order the event times first by their parent lineage, then by the time of occurrence
      event.times<-event.times[order(event.parent.lins,event.times)]
      #pick whether each event is a speciation event (1) or an extinction event (0) by drawing from a binomial
      #distribution with 1 trial and probability equal to the probability of speciation relative to extinction
      event.types<-rbinom(sum(lin.events),1,spec.rate/(spec.rate+ext.rate))
      #initialize(/reset) a vector for recording which parents are going extinct between the current time point and
      #the next time point
      parents.going.ext<-NULL
      #if an extinction is occuring in between the current time point and the next...
      if(0%in%event.types){
        #find the earliest times at which an extinction occurs for each parental lineage
        earliest.exts<-which(event.times%in%tapply(event.times[event.types==0],event.parent.lins[event.types==0],min))
        parents.going.ext<-event.parent.lins[earliest.exts]
        #find the last event that occurs for each parental lineage
        last.events<-which(event.times%in%tapply(event.times[event.parent.lins%in%parents.going.ext],
                                                 event.parent.lins[event.parent.lins%in%parents.going.ext],max))
        #use the above information to find any events that need to be removed because they occur to a parental lineage
        #after it has gone extinct
        rem.ind<-mult.seq(c(earliest.exts+1,last.events))
        #if there events that need to be removed...
        if(length(rem.ind)>0){
          #remove them from vector of parental lineages for each event...
          event.parent.lins<-event.parent.lins[-rem.ind]
          #do the same for vector of event starting times...
          event.starts<-event.starts[-rem.ind]
          #again for vector of trait values at event starting times...
          event.starts.trait.vals<-event.starts.trait.vals[-rem.ind]
          #yet again for vector of specific event times...
          event.times<-event.times[-rem.ind]
          #and lastly for vector of event types
          event.types<-event.types[-rem.ind]
        }
      }
      #initialize vector of child.lineage identifiers, making every child 'NA' (hence every event is asummed to be
      #an extinction for now)
      event.child.lins<-rep(NA,length(event.parent.lins))
      #if the event type is speciation, alter the child lineage to be some number between the number of previously
      #exisiting lineages and the number of newly-generated lineages
      event.child.lins[event.types==1]<-nrow(trait)+1:sum(event.types)
      #create a vector of trait values for the trait value at the time of each event by linearly interpolating between
      #trait values of parental lineage at current time point and next time point
      event.trait.vals<-(trait[event.parent.lins,as.character(t+res)]-event.starts.trait.vals)*
        ((event.times-event.starts)/(res-event.starts))+event.starts.trait.vals
      #now create new node matrix to be grafted on to existing node matrix by binding together the parent lineages,
      #child lineages, absolute times of events, and trait values at the times of those events
      new.nodes<-cbind(event.parent.lins,event.child.lins,t+event.times,event.trait.vals)
      #now create new trait matrix of newly-generated lineages to be grafted on to existing trait matrix, leaving
      #every trait value an NA for now
      new.lins<-matrix(NA,ncol=ncol(trait),nrow=sum(event.types))
      #name to columns of the new trait matrix according to the time point vector
      colnames(new.lins)<-time.vec
      #simulate trait values of newly-generated lineages at the next time point; trait values at last event plus...
      new.lins[,as.character(t+res)]<-event.trait.vals[!is.na(event.child.lins)]+
        #OU component (scaled to intervening time)
        alpha*(theta-event.trait.vals[!is.na(event.child.lins)])*(res-event.times[!is.na(event.child.lins)])+
        #BM component (scaled to intervening time)
        rnorm(nrow(new.lins),
              trend*(res-event.times[!is.na(event.child.lins)]),
              rate*(res-event.times[!is.na(event.child.lins)]))
      #bind new node and trait matrices to old ones
      nodes<-rbind(nodes,new.nodes)
      trait<-rbind(trait,new.lins)
      #if there are parental lineages going extinct...
      if(!is.null(parents.going.ext)){
        #make their trait values at the next time point NA (since we are done interpolating trait values, these values
        #are no longer required for simulation and technically shouldn't exist, so...)
        trait[parents.going.ext,as.character(t+res)]<-NA
      }
      ##SIMULATING EVENTS ON NEWLY-GENERATED LINEAGES##
      #update event starting times by looking up the times of all speciation events in the node matrix
      event.starts<-nodes[which(nodes[,"child.lin"]%in%1:nrow(trait)),"time"]-t
      #update trait values at event starting times by doing the same with trait values
      event.starts.trait.vals<-nodes[which(nodes[,"child.lin"]%in%1:nrow(trait)),"trait"]
      #re-simulate the number of events occurring to each lineage before next time point, scaling rate by time
      #between event start and next time point
      lin.events<-rpois(nrow(trait),(spec.rate+ext.rate)*(res-event.starts))
      #make sure the number of events for any lineage that isn't newly-generated is set to 0, taking advantage of the
      #fact that the number of newly-generated lineages is equal to the number of rows in the grafted-on new trait
      #matrix
      lin.events[1:(nrow(trait)-nrow(new.lins))]<-0
    }
    #If the current time point is the last one and not all lineages are already extinct
    if(t==time.vec[(length(time.vec))-1]&!all(is.na(trait[,as.character(t+res)]))){
      #Find extant lineages
      extant.lins<-which(!is.na(trait[,as.character(t+res)]))
      #Generate a new node matrix to be grafted on to existing node matrix by generating "extinction" nodes for
      #extant lineages
      new.nodes<-cbind(extant.lins,rep(NA),rep(t+res),trait[extant.lins,as.character(t+res)])
      #Bind new node matrix to existing node matrix
      nodes<-rbind(nodes,new.nodes)
    }
  }
  #Intialize new tip label vector to be bound to existing node matrix (all NA's for now) 
  tip.lab<-rep(NA,nrow(nodes))
  #Find which nodes correspond to tips (i.e., child lineage is NA)
  tip.index<-which(is.na(nodes[,"child.lin"]))
  #Label tip according to its "parent lineage"
  tip.lab[tip.index]<-nodes[is.na(nodes[,"child.lin"]),"parent.lin"]
  #Bind tip label vector to existing node matrix
  nodes<-cbind(nodes,tip.lab)
  #Erase any row names for the node matrix
  rownames(nodes)<-NULL
  #Save a vector of parameters with which the phylogeny was generated
  params<-c(res,spec.rate,ext.rate,rate,init.trait,trend,alpha,theta)
  names(params)<-c("resolution","speciation rate","extinction rate","trait evolution rate","initial trait value",
                   "trait evolution trend","pull towards trait optimum","trait optimum")
  #Output phylosim 'object' as list of time vector, trait matrix, node matrix, and parameter vector
  out<-list(time.pts=time.vec,trait=trait,nodes=nodes,params=params)
  class(out)<-"phylosim"
  out
}