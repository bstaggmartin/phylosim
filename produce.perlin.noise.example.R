library("magick")
source("perlin.noise.gen.R")
source("render.noise.R")
set.seed(123)
axis<-seq(1,10,length.out=101)[-101]
noise<-gen.vec.field(no.dim=3)
interp.noise<-sapply(axis,function(ii) matrix(get.noise.value(vec.field=noise,points=expand.grid(axis,axis,ii)),ncol=100))
interp.noise<-array(t(interp.noise),dim=c(100,100,100))
render.noise(noise=interp.noise,png.out.dir="noise_example",gif.name="render",y.compression=2)
