render.noise<-function(noise,png.name="slice",gif.name="render",png.out.dir,
                       pre.clear.pngs=T,post.clear.pngs=T,width=600,height=420,
                       y.compression=1,x.offset=0.5,y.offset=0.05,line.width=8,gray.max=0.7,alpha=0.6,fps=20){
  def.dir<-getwd()
  setwd(png.out.dir)
  if(pre.clear.pngs){
    file.remove(list.files(pattern=".png"))
  }
  png(paste(png.name,"%02d.png",sep=""),width=width,height=height)
  col.vec=gray(seq(gray.max,0,length.out=dim(noise)[2]),alpha=alpha)
  for(j in 1:dim(noise)[3]){
    plot((y.compression*noise[,1,j])~seq(1,dim(noise)[1]),type="l",xaxt='n',yaxt='n',xlab="",ylab="",lwd=line.width,col=col.vec[1],bty='n',
         ylim=c(min(noise)*y.compression-(dim(noise)[2]+1)*y.offset,max(noise)*y.compression),xlim=c(1,dim(noise)[1]+(dim(noise)[1]+1)*x.offset))
    for(i in 2:dim(noise)[2]){
      lines(y=(y.compression*noise[,i,j])-(i-1)*y.offset,x=seq(1,dim(noise)[1])+(i-1)*x.offset,lwd=line.width,col=col.vec[i])
    }
  }
  dev.off()
  slices<-paste("slice",c(paste("0",1:9,sep=""),10:dim(noise)[3]),".png",sep="")
  images<-list(NULL)
  for(i in 1:dim(noise)[3]){
    images[[i]]<-image_read(slices[i])
  }
  frames<-image_join(c(images))
  animation<-image_animate(frames,fps=fps)
  image_write(animation,paste(gif.name,".gif",sep=""))
  if(post.clear.pngs==T){
    file.remove(list.files(pattern=".png"))
  }
  setwd(def.dir)
}
