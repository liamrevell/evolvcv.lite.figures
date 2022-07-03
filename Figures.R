library(phytools)
library(xtable)
library(RColorBrewer)

## Figure 1
pdf(file="revision-Figure1.pdf",width=9,height=7)
data(sunfish.tree)
data(sunfish.data)
layout(matrix(c(1,2,3,3),2,2,byrow=TRUE),heights=c(1,0.7))
cols<-setNames(c("blue","red"),c("non","pisc"))
plot(sunfish.tree,cols,ftype="i",
	mar=c(0.1,2.1,3.1,1.1),fsize=0.8)
mtext("a)",line=0,adj=0)
legend("topleft",c("non-piscivory","piscivory"),
	pch=15,col=cols,bty="n",pt.cex=1.5)
par(mar=c(5.1,4.1,3.1,1.1),cex.axis=0.8)
phylomorphospace(sunfish.tree,sunfish.data[,2:3],colors=cols,
	ftype="off",bty="n",xlab="relative gape width",
	ylab="relative buccal length",node.by.map=TRUE,
	las=1,cex.axis=0.8)
mtext("b)",line=0,adj=0)
fit<-evol.vcv(sunfish.tree,sunfish.data[,2:3])
par(mar=c(0.1,2.1,3.1,1.1))
plot(NA,xlim=c(0,1),ylim=c(0.2,1),axes=FALSE)
mtext("c)",line=0,adj=0)
par(cex=0.9)
text(x=0,y=1,"Fitted one-matrix model:",pos=4)
text(x=0,y=0.9,"      R = ",pos=4)
h<-0.05
lines(rep(0.155,2)-h,c(0.75,0.95))
lines(c(0.155,0.16)-h,rep(0.95,2))
lines(c(0.155,0.16)-h,rep(0.75,2))
R=fit$R.single
text(0.16-h,0.9,sprintf("%1.3f",R[1,1]),pos=4)
text(0.16-h,0.8,sprintf("%1.3f",R[1,2]),pos=4)
text(0.24-h,0.9,sprintf("%1.3f",R[2,1]),pos=4)
text(0.24-h,0.8,sprintf("%1.3f",R[2,2]),pos=4)
lines(rep(0.32,2)-h,c(0.75,0.95))
lines(c(0.32,0.315)-h,rep(0.95,2))
lines(c(0.32,0.315)-h,rep(0.75,2))
text(x=0.32-h,y=0.9,paste(",          ","log(L) = ",
	round(fit$logL1,2)),pos=4)
v<-0.4
text(x=0,y=1-v,"Fitted two-matrix model:",pos=4)
text(x=0,y=0.9-v,expression(paste("      ",
	R[non-piscivorous],"  =")),pos=4)
h<-0.05
lines(rep(0.155,2)+h,c(0.75,0.95)-v)
lines(c(0.155,0.16)+h,rep(0.95,2)-v)
lines(c(0.155,0.16)+h,rep(0.75,2)-v)
R=fit$R.multiple[,,1]
text(0.16+h,0.9-v,sprintf("%1.3f",R[1,1]),pos=4)
text(0.16+h,0.8-v,sprintf("%1.3f",R[1,2]),pos=4)
text(0.24+h,0.9-v,sprintf("%1.3f",R[2,1]),pos=4)
text(0.24+h,0.8-v,sprintf("%1.3f",R[2,2]),pos=4)
lines(rep(0.32,2)+h,c(0.75,0.95)-v)
lines(c(0.32,0.315)+h,rep(0.95,2)-v)
lines(c(0.32,0.315)+h,rep(0.75,2)-v)
text(x=0.32+h,y=0.9-v,expression(paste(",    ","      ",
	R[piscivorous],"  =")),pos=4)
h<-0.32+2*h
lines(rep(0.155+h,2),c(0.75,0.95)-v)
lines(c(0.155+h,0.16+h),rep(0.95,2)-v)
lines(c(0.155+h,0.16+h),rep(0.75,2)-v)
R=fit$R.multiple[,,2]
text(0.16+h,0.9-v,sprintf("%1.3f",R[1,1]),pos=4)
text(0.16+h,0.8-v,sprintf("%1.3f",R[1,2]),pos=4)
text(0.24+h,0.9-v,sprintf("%1.3f",R[2,1]),pos=4)
text(0.24+h,0.8-v,sprintf("%1.3f",R[2,2]),pos=4)
lines(rep(0.32+h,2),c(0.75,0.95)-v)
lines(c(0.32+h,0.315+h),rep(0.95,2)-v)
lines(c(0.32+h,0.315+h),rep(0.75,2)-v)
text(x=0.32+h,y=0.9-v,paste(",          ","log(L) = ",
	round(fit$logL.multiple,2)),pos=4)
dev.off()
## end Figure 1

## save workspace
save.image(file="Figures.Rdata")
## end save workspace

## Table 1
fit.all<-evolvcv.lite(sunfish.tree,sunfish.data[,2:3],
	models="all models")
Table1<-matrix(NA,8,9,dimnames=list(
	names(fit.all),
	c("Model","sig2[1,1]","sig2[1,2]","sig2[2,1]","sig2[2,2]",
	"r[1]","r[2]","log(L)","AIC")))
Table1<-as.data.frame(Table1)
for(i in 1:length(fit.all)){
	Table1[i,1]<-fit.all[[i]]$description
	if(i==1){
		Table1[i,2]<-fit.all[[i]]$R[1,1]
		Table1[i,4]<-fit.all[[i]]$R[2,2]
		Table1[i,6]<-fit.all[[i]]$R[1,2]/sqrt(fit.all[[i]]$R[1,1]*
			fit.all[[i]]$R[2,2])
	} else {
		Table1[i,2]<-fit.all[[i]]$R[[1]][1,1]
		if(fit.all[[i]]$R[[2]][1,1]!=Table1[i,2])
			Table1[i,3]<-fit.all[[i]]$R[[2]][1,1]
		Table1[i,4]<-fit.all[[i]]$R[[1]][2,2]
		if(fit.all[[i]]$R[[2]][2,2]!=Table1[i,4])
			Table1[i,5]<-fit.all[[i]]$R[[2]][2,2]
		Table1[i,6]<-fit.all[[i]]$R[[1]][1,2]/sqrt(
			fit.all[[i]]$R[[1]][1,1]*
			fit.all[[i]]$R[[1]][2,2])
		r2<-fit.all[[i]]$R[[2]][1,2]/sqrt(
			fit.all[[i]]$R[[2]][1,1]*
			fit.all[[i]]$R[[2]][2,2])
		if(r2!=Table1[i,6]) Table1[i,7]<-r2
	}
	Table1[i,8]<-fit.all[[i]]$logLik
	Table1[i,9]<-fit.all[[i]]$AIC
}
xtable(Table1)
## end Table 1

## save workspace
save.image(file="Figures.Rdata")
## end save workspace

## Figure 2
trop.data<-read.csv(file="tropidurid-data.csv",row.names=1)
trop.tree<-read.simmap(file="tropidurid-tree.tre",format="nexus",
	version=1.5)
pdf(file="revision-Figure2-v2.pdf",width=10,height=6)
cols<-setNames(c("white","black"),c("n_rock","rock"))
tmp<-trop.tree
tmp$tip.label<-paste(" ",tmp$tip.label)
plot(tmp,outline=TRUE,colors=cols,fsize=0.6,
	direction="upwards",ftype="i",offset=0.1,
	lwd=3)
legend(x="bottomright",c("non-rock-dwelling","rock-dwelling"),
	pch=22,pt.bg=c("white","black"),pt.cex=1.5,bty="n")
dev.off()
## end Figure 2

## save workspace
save.image(file="Figures.Rdata")
## end save workspace

## Table 2
trop.fits<-evolvcv.lite(trop.tree,trop.data,models="all models")
trop.aic<-sapply(trop.fits,function(x) x$AIC)	
trop.w<-aic.w(trop.aic)

ii<-order(trop.aic)[1:8]
bfits<-trop.fits[ii]
class(bfits)<-class(trop.fits)

TROP.RESULTS<-data.frame(model=names(bfits),
	sig11=sapply(bfits,function(x) if(x$description=="common rates, common correlation") x$R[1,1] else x$R[[1]][1,1]),
	sig12=sapply(bfits,function(x) if(x$description=="common rates, common correlation") x$R[1,1] else x$R[[2]][1,1]),
	sig21=sapply(bfits,function(x) if(x$description=="common rates, common correlation") x$R[2,2] else x$R[[1]][2,2]),
	sig22=sapply(bfits,function(x) if(x$description=="common rates, common correlation") x$R[2,2] else x$R[[2]][2,2]),
	r1=sapply(bfits,function(x) if(x$description=="common rates, common correlation") cov2cor(x$R)[1,2] else cov2cor(x$R[[1]])[1,2]),
	r2=sapply(bfits,function(x) if(x$description=="common rates, common correlation") cov2cor(x$R)[1,2] else cov2cor(x$R[[2]])[1,2]),
	logL=sapply(bfits,function(x) x$logLik),
	AIC=trop.aic[ii],Akaike.w=as.vector(trop.w[ii]))
rownames(TROP.RESULTS)<-NULL
xtable(TROP.RESULTS)
## end Table 2

## save workspace
save.image(file="Figures.Rdata")
## end save workspace

## SIMULATION

nsim<-40
ntaxa<-100
nmin<-20

## simulate trees
set.seed(99)
Q<-0.5*matrix(c(-2,1,1,1,-2,1,1,1,-2),3,3,dimnames=list(1:3,1:3))
trees<-list()
class(trees)<-c("multiSimmap","multiPhylo")
for(i in 1:nsim){
	tt<-pbtree(n=ntaxa,scale=1)
	mt<-sim.history(tt,Q)
	xx<-as.factor(getStates(mt,"tips"))
	while(any(summary(xx)<nmin)){
		mt<-sim.history(tt,Q)
		xx<-as.factor(getStates(mt,"tips"))
	}
	trees[[i]]<-mt
}

## Model 1 (no difference)
set.seed(11)
MODEL1<-matrix(NA,nsim,8,dimnames=list(1:nsim,paste("AIC(",
	c("1","2","2b","2c","3","3b","3c","4"),")",sep="")))
GENERATING1<-matrix(NA,nsim,3,dimnames=list(1:nsim,
	c("r","v1","v2")))
fits.model1<-list()
for(i in 1:nsim){
	r<-runif(n=1,min=-1,max=1)
	v1<-exp(rnorm(n=1))
	v2<-exp(rnorm(n=1))
	GENERATING1[i,]<-c(r,v1,v2)
	V<-matrix(c(v1,r*sqrt(v1*v2),r*sqrt(v1*v2),v2),2,2)
	X<-sim.corrs(trees[[i]],V)
	fits.model1[[i]]<-evolvcv.lite(trees[[i]],X,models="all models",
		tol=1e-4,try.iter=20)
	MODEL1[i,]<-sapply(fits.model1[[i]],function(x) x$AIC)
}
wMODEL1<-t(apply(MODEL1,1,aic.w))
dMODEL1<-t(apply(MODEL1,1,function(x) x-min(x)))

## Model 2 (different rates, common correlation)
set.seed(21)
MODEL2<-matrix(NA,nsim,8,dimnames=list(1:nsim,paste("AIC(",
	c("1","2","2b","2c","3","3b","3c","4"),")",sep="")))
GENERATING2<-matrix(NA,nsim,7,dimnames=list(1:nsim,
	c("r",
	paste("v1[",1:3,"]",sep=""),
	paste("v2[",1:3,"]",sep=""))))
fits.model2<-list()
for(i in 1:nsim){
	r<-runif(n=1,min=-1,max=1)
	v1<-exp(rnorm(n=3))
	v2<-exp(rnorm(n=3))
	GENERATING2[i,]<-c(r,v1,v2)
	V<-setNames(list(
		matrix(c(v1[1],r*sqrt(v1[1]*v2[1]),r*sqrt(v1[1]*v2[1]),v2[1]),2,2),
		matrix(c(v1[2],r*sqrt(v1[2]*v2[2]),r*sqrt(v1[2]*v2[2]),v2[2]),2,2),
		matrix(c(v1[3],r*sqrt(v1[3]*v2[3]),r*sqrt(v1[3]*v2[3]),v2[3]),2,2)),
		1:3)
	X<-sim.corrs(trees[[i]],V)
	fits.model2[[i]]<-evolvcv.lite(trees[[i]],X,models="all models",
		tol=1e-4,try.iter=30) ## increase try iterations
	MODEL2[i,]<-sapply(fits.model2[[i]],function(x) x$AIC)
}
wMODEL2<-t(apply(MODEL2,1,aic.w))
dMODEL2<-t(apply(MODEL2,1,function(x) x-min(x)))

cat("\n\n\n-------------------------------------------\n")
cat("Done Model 2.\n")
cat("-------------------------------------------\n\n\n")

## save workspace
save.image(file="Figures.Rdata")
## end save workspace

## Model 2b (different rates trait 1, common correlation)
set.seed(22)
MODEL2b<-matrix(NA,nsim,8,dimnames=list(1:nsim,paste("AIC(",
	c("1","2","2b","2c","3","3b","3c","4"),")",sep="")))
GENERATING2b<-matrix(NA,nsim,5,dimnames=list(1:nsim,
	c("r",
	paste("v1[",1:3,"]",sep=""),
	"v2")))
fits.model2b<-list()
for(i in 1:nsim){
	r<-runif(n=1,min=-1,max=1)
	v1<-exp(rnorm(n=3))
	v2<-rep(exp(rnorm(n=1)),3)
	GENERATING2b[i,]<-c(r,v1,v2[1])
	V<-setNames(list(
		matrix(c(v1[1],r*sqrt(v1[1]*v2[1]),r*sqrt(v1[1]*v2[1]),v2[1]),2,2),
		matrix(c(v1[2],r*sqrt(v1[2]*v2[2]),r*sqrt(v1[2]*v2[2]),v2[2]),2,2),
		matrix(c(v1[3],r*sqrt(v1[3]*v2[3]),r*sqrt(v1[3]*v2[3]),v2[3]),2,2)),
		1:3)
	X<-sim.corrs(trees[[i]],V)
	fits.model2b[[i]]<-evolvcv.lite(trees[[i]],X,models="all models",
		tol=1e-4,try.iter=20)
	MODEL2b[i,]<-sapply(fits.model2b[[i]],function(x) x$AIC)
}
wMODEL2b<-t(apply(MODEL2b,1,aic.w))
dMODEL2b<-t(apply(MODEL2b,1,function(x) x-min(x)))

cat("\n\n\n-------------------------------------------\n")
cat("Done Model 2b.\n")
cat("-------------------------------------------\n\n\n")

## save workspace
save.image(file="Figures.Rdata")
## end save workspace

## Model 2c (different rates trait 2, common correlation)
set.seed(23)
MODEL2c<-matrix(NA,nsim,8,dimnames=list(1:nsim,paste("AIC(",
	c("1","2","2b","2c","3","3b","3c","4"),")",sep="")))
GENERATING2c<-matrix(NA,nsim,5,dimnames=list(1:nsim,
	c("r","v1",
	paste("v2[",1:3,"]",sep=""))))
fits.model2c<-list()
for(i in 1:nsim){
	r<-runif(n=1,min=-1,max=1)
	v1<-rep(exp(rnorm(n=1)),3)
	v2<-exp(rnorm(n=3))
	GENERATING2c[i,]<-c(r,v1[1],v2)
	V<-setNames(list(
		matrix(c(v1[1],r*sqrt(v1[1]*v2[1]),r*sqrt(v1[1]*v2[1]),v2[1]),2,2),
		matrix(c(v1[2],r*sqrt(v1[2]*v2[2]),r*sqrt(v1[2]*v2[2]),v2[2]),2,2),
		matrix(c(v1[3],r*sqrt(v1[3]*v2[3]),r*sqrt(v1[3]*v2[3]),v2[3]),2,2)),
		1:3)
	X<-sim.corrs(trees[[i]],V)
	fits.model2c[[i]]<-evolvcv.lite(trees[[i]],X,models="all models",
		tol=1e-4,try.iter=40) # increase try iterations
	MODEL2c[i,]<-sapply(fits.model2c[[i]],function(x) x$AIC)
}
wMODEL2c<-t(apply(MODEL2c,1,aic.w))
dMODEL2c<-t(apply(MODEL2c,1,function(x) x-min(x)))

cat("\n\n\n-------------------------------------------\n")
cat("Done Model 2c.\n")
cat("-------------------------------------------\n\n\n")

## save workspace
save.image(file="Figures.Rdata")
## end save workspace

## Model 3 (same rates, different correlations)
set.seed(31)
MODEL3<-matrix(NA,nsim,8,dimnames=list(1:nsim,paste("AIC(",
	c("1","2","2b","2c","3","3b","3c","4"),")",sep="")))
GENERATING3<-matrix(NA,nsim,5,dimnames=list(1:nsim,
	c(paste("r[",1:3,"]",sep=""),"v1","v2")))
fits.model3<-list()
for(i in 1:nsim){
	r<-runif(n=3,min=-1,max=1)
	v1<-rep(exp(rnorm(n=1)),3)
	v2<-rep(exp(rnorm(n=1)),3)
	GENERATING3[i,]<-c(r,v1[1],v2[2])
	V<-setNames(list(
		matrix(c(v1[1],r[1]*sqrt(v1[1]*v2[1]),r[1]*sqrt(v1[1]*v2[1]),v2[1]),2,2),
		matrix(c(v1[2],r[2]*sqrt(v1[2]*v2[2]),r[2]*sqrt(v1[2]*v2[2]),v2[2]),2,2),
		matrix(c(v1[3],r[3]*sqrt(v1[3]*v2[3]),r[3]*sqrt(v1[3]*v2[3]),v2[3]),2,2)),
		1:3)
	X<-sim.corrs(trees[[i]],V)
	fits.model3[[i]]<-evolvcv.lite(trees[[i]],X,models="all models",
		tol=1e-4,try.iter=30) ## changed to 30
	MODEL3[i,]<-sapply(fits.model3[[i]],function(x) x$AIC)
}
wMODEL3<-t(apply(MODEL3,1,aic.w))
dMODEL3<-t(apply(MODEL3,1,function(x) x-min(x)))

cat("\n\n\n-------------------------------------------\n")
cat("Done Model 3.\n")
cat("-------------------------------------------\n\n\n")

## save workspace
save.image(file="Figures.Rdata")
## end save workspace

## Model 3b (different rates trait 1, different correlations)
set.seed(32)
MODEL3b<-matrix(NA,nsim,8,dimnames=list(1:nsim,paste("AIC(",
	c("1","2","2b","2c","3","3b","3c","4"),")",sep="")))
GENERATING3b<-matrix(NA,nsim,7,dimnames=list(1:nsim,
	c(paste("r[",1:3,"]",sep=""),
	paste("v1[",1:3,"]",sep=""),"v2")))
fits.model3b<-list()
for(i in 1:nsim){
	r<-runif(n=3,min=-1,max=1)
	v1<-exp(rnorm(n=3))
	v2<-rep(exp(rnorm(n=1)),3)
	V<-setNames(list(
		matrix(c(v1[1],r[1]*sqrt(v1[1]*v2[1]),r[1]*sqrt(v1[1]*v2[1]),v2[1]),2,2),
		matrix(c(v1[2],r[2]*sqrt(v1[2]*v2[2]),r[2]*sqrt(v1[2]*v2[2]),v2[2]),2,2),
		matrix(c(v1[3],r[3]*sqrt(v1[3]*v2[3]),r[3]*sqrt(v1[3]*v2[3]),v2[3]),2,2)),
		1:3)
	GENERATING3b[i,]<-c(r,v1,v2[1])
	X<-sim.corrs(trees[[i]],V)
	fits.model3b[[i]]<-evolvcv.lite(trees[[i]],X,models="all models",
		tol=1e-4,try.iter=20)
	MODEL3b[i,]<-sapply(fits.model3b[[i]],function(x) x$AIC)
}
wMODEL3b<-t(apply(MODEL3b,1,aic.w))
dMODEL3b<-t(apply(MODEL3b,1,function(x) x-min(x)))

cat("\n\n\n-------------------------------------------\n")
cat("Done Model 3b.\n")
cat("-------------------------------------------\n\n\n")

## save workspace
save.image(file="Figures.Rdata")
## end save workspace

## Model 3c (different rates trait 2, different correlations)
set.seed(33)
MODEL3c<-matrix(NA,nsim,8,dimnames=list(1:nsim,paste("AIC(",
	c("1","2","2b","2c","3","3b","3c","4"),")",sep="")))
GENERATING3c<-matrix(NA,nsim,7,dimnames=list(1:nsim,
	c(paste("r[",1:3,"]",sep=""),"v1",
	paste("v2[",1:3,"]",sep=""))))
fits.model3c<-list()
for(i in 1:nsim){
	r<-runif(n=3,min=-1,max=1)
	v1<-rep(exp(rnorm(n=1)),3)
	v2<-exp(rnorm(n=3))
	V<-setNames(list(
		matrix(c(v1[1],r[1]*sqrt(v1[1]*v2[1]),r[1]*sqrt(v1[1]*v2[1]),v2[1]),2,2),
		matrix(c(v1[2],r[2]*sqrt(v1[2]*v2[2]),r[2]*sqrt(v1[2]*v2[2]),v2[2]),2,2),
		matrix(c(v1[3],r[3]*sqrt(v1[3]*v2[3]),r[3]*sqrt(v1[3]*v2[3]),v2[3]),2,2)),
		1:3)
	GENERATING3c[i,]<-c(r,v1[1],v2)
	X<-sim.corrs(trees[[i]],V)
	fits.model3c[[i]]<-evolvcv.lite(trees[[i]],X,models="all models",
		tol=1e-4,try.iter=20) 
	MODEL3c[i,]<-sapply(fits.model3c[[i]],function(x) x$AIC)
}
wMODEL3c<-t(apply(MODEL3c,1,aic.w))
dMODEL3c<-t(apply(MODEL3c,1,function(x) x-min(x)))

cat("\n\n\n-------------------------------------------\n")
cat("Done Model 3c.\n")
cat("-------------------------------------------\n\n\n")

## save workspace
save.image(file="Figures.Rdata")
## end save workspace

## Model 4 (no common structure)
set.seed(41)
MODEL4<-matrix(NA,nsim,8,dimnames=list(1:nsim,paste("AIC(",
	c("1","2","2b","2c","3","3b","3c","4"),")",sep="")))
GENERATING4<-matrix(NA,nsim,9,dimnames=list(1:nsim,
	c(paste("r[",1:3,"]",sep=""),
	paste("v1[",1:3,"]",sep=""),
	paste("v2[",1:3,"]",sep=""))))
fits.model4<-list()
for(i in 1:nsim){
	r<-runif(n=3,min=-1,max=1)
	v1<-exp(rnorm(n=3))
	v2<-exp(rnorm(n=3))
	V<-setNames(list(
		matrix(c(v1[1],r[1]*sqrt(v1[1]*v2[1]),r[1]*sqrt(v1[1]*v2[1]),v2[1]),2,2),
		matrix(c(v1[2],r[2]*sqrt(v1[2]*v2[2]),r[2]*sqrt(v1[2]*v2[2]),v2[2]),2,2),
		matrix(c(v1[3],r[3]*sqrt(v1[3]*v2[3]),r[3]*sqrt(v1[3]*v2[3]),v2[3]),2,2)),
		1:3)
	GENERATING4[i,]<-c(r,v1,v2)
	X<-sim.corrs(trees[[i]],V)
	fits.model4[[i]]<-evolvcv.lite(trees[[i]],X,models="all models",
		tol=1e-4,try.iter=20)
	MODEL4[i,]<-sapply(fits.model4[[i]],function(x) x$AIC)
}
wMODEL4<-t(apply(MODEL4,1,aic.w))
dMODEL4<-t(apply(MODEL4,1,function(x) x-min(x)))

cat("\n\n\n-------------------------------------------\n")
cat("Done Model 4.\n")
cat("-------------------------------------------\n\n\n")

## save workspace
save.image(file="Figures.Rdata")
## end save workspace

W<-rbind(
	colMeans(wMODEL1),
	colMeans(wMODEL2),
	colMeans(wMODEL2b),
	colMeans(wMODEL2c),
	colMeans(wMODEL3),
	colMeans(wMODEL3b),
	colMeans(wMODEL3c),
	colMeans(wMODEL4))
rownames(W)<-colnames(W)<-paste("model",c("1","2","2b","2c","3","3b","3c","4"))

pdf(file="revision-Figure4-v2.pdf",width=10,height=5)
par(mar=c(1.1,1.1,1.1,1.1))
plot(NA,xlim=c(-0.15,1.15),ylim=c(0,1.1),axes=FALSE)
for(i in 1:nrow(W)){
	for(j in 1:ncol(W)){
		xx<-c(0,1/8)+(i-1)/8
		xx<-c(xx,xx[2:1])
		yy<-c(7/8,1)-(j-1)/8
		yy<-c(yy[1],yy,yy[2])
		polygon(x=xx,y=yy,
			col=rgb(colorRamp(c("white",palette()[2]))(W[j,i]/max(W)),
			maxColorValue=255))
		if(j==1) text(mean(xx),1,rownames(W)[i],pos=3)
		if(i==1) text(0,mean(yy),rownames(W)[j],pos=2)
	}
}
text(-0.15,0.5,"Generating model in simulation",cex=1.2,srt=90)
text(0.55,1.1,"Mean Akaike weight by model",cex=1.2)

cols<-rgb(colorRamp(c("white",palette()[2]))(seq(0,1,length.out=100)),
	maxColorValue=255)
LWD<-diff(par()$usr[1:2])/dev.size("px")[1]
nticks<-7
Y<-cbind(seq(0.05,0.95,length.out=nticks),
	seq(0.05,0.95,length.out=nticks))
X<-cbind(rep(1.1+LWD*10/2,nticks),
	rep(1.1+LWD*10/2+0.01,nticks))
for(i in 1:nrow(Y)) lines(X[i,],Y[i,])
add.color.bar(0.9,cols,title="",lims=NULL,
	digits=2,direction="upwards",subtitle="",lwd = 15,x=1.1,y=0.05,
	prompt=FALSE)
text(x=1.05,y=0.5,"mean Akaike weight",srt=90)
text(x=X[,2],y=Y[,2],signif(seq(0,max(W),length.out=nticks),
	2),pos=4,cex=0.9)
dev.off()

## Figure 3
## Model 3b (different rates trait 1, different correlations)
set.seed(32)
ii<-37
r<-GENERATING3b[ii,1:3]
v1<-GENERATING3b[ii,4:6]
v2<-rep(GENERATING3b[ii,7],3)
V<-setNames(list(
	matrix(c(v1[1],r[1]*sqrt(v1[1]*v2[1]),r[1]*sqrt(v1[1]*v2[1]),v2[1]),2,2),
	matrix(c(v1[2],r[2]*sqrt(v1[2]*v2[2]),r[2]*sqrt(v1[2]*v2[2]),v2[2]),2,2),
	matrix(c(v1[3],r[3]*sqrt(v1[3]*v2[3]),r[3]*sqrt(v1[3]*v2[3]),v2[3]),2,2)),
	1:3)
V
X<-sim.corrs(trees[[ii<-sample(1:nsim,1)]],V)
cols<-setNames(brewer.pal(3,"Dark2"),1:3)
pdf(file="revision-Figure3.pdf",width=12,height=7)
par(fg="#404040",mfrow=c(1,2))
plot(trees[[ii]],lwd=2,colors=cols,outline=TRUE,ftype="off",
	mar=c(2.1,4.1,3.1,1.1))
legend("bottomleft",paste("regime",1:3),pch=22,col="#404040",
	pt.bg=cols,bty="n",pt.cex=1.2,cex=0.8)
mtext("a)",adj=0,cex=1.2)
par(mar=c(5.1,4.1,3.1,1.1))
phylomorphospace(trees[[ii]],X,ftype="off",colors=cols,node.size=c(0,1),
	node.by.map=TRUE,bty="n",xlab=expression(x[1]),ylab=expression(x[2]),
	las=1)
mtext("b)",adj=0,cex=1.2)
dev.off()

bestModel<-matrix(NA,8,4,dimnames=
	list(paste("model",c("1","2","2b","2c","3","3b","3c","4")),
	c("Best","2nd best","3rd best",">=4th")))
for(i in 1:3)
	bestModel[1,i]<-mean(apply(MODEL1,1,function(x) which(order(x)==1)==i))
bestModel[1,4]<-1-sum(bestModel[1,1:3])

for(i in 1:3)
	bestModel[2,i]<-mean(apply(MODEL2,1,function(x) which(order(x)==2)==i))
bestModel[2,4]<-1-sum(bestModel[2,1:3])

for(i in 1:3)
	bestModel[3,i]<-mean(apply(MODEL2b,1,function(x) which(order(x)==3)==i))
bestModel[3,4]<-1-sum(bestModel[3,1:3])

for(i in 1:3)
	bestModel[4,i]<-mean(apply(MODEL2c,1,function(x) which(order(x)==4)==i))
bestModel[4,4]<-1-sum(bestModel[4,1:3])

for(i in 1:3)
	bestModel[5,i]<-mean(apply(MODEL3,1,function(x) which(order(x)==5)==i))
bestModel[5,4]<-1-sum(bestModel[5,1:3])

for(i in 1:3)
	bestModel[6,i]<-mean(apply(MODEL3b,1,function(x) which(order(x)==6)==i))
bestModel[6,4]<-1-sum(bestModel[6,1:3])

for(i in 1:3)
	bestModel[7,i]<-mean(apply(MODEL3c,1,function(x) which(order(x)==7)==i))
bestModel[7,4]<-1-sum(bestModel[7,1:3])

for(i in 1:3)
	bestModel[8,i]<-mean(apply(MODEL4,1,function(x) which(order(x)==8)==i))
bestModel[8,4]<-1-sum(bestModel[8,1:3])

xtable(bestModel)

## save workspace
save.image(file="Figures.Rdata")
## end save workspace

## BEGIN PARAMETER ESTIMATION

## MODEL 1
foo<-function(x) c(cov2cor(x[[1]]$R)[1,2],x[[1]]$R[1,1],x[[1]]$R[2,2])
ESTIMATED1<-t(sapply(fits.model1,foo))
dimnames(ESTIMATED1)<-dimnames(GENERATING1)

## MODEL 2
foo<-function(x) c(cov2cor(x[[2]]$R[["1"]])[1,2],
	x[[2]]$R[["1"]][1,1],x[[2]]$R[["2"]][1,1],x[[2]]$R[["3"]][1,1],
	x[[2]]$R[["1"]][2,2],x[[2]]$R[["2"]][2,2],x[[2]]$R[["3"]][2,2])
ESTIMATED2<-t(sapply(fits.model2,foo))
dimnames(ESTIMATED2)<-dimnames(GENERATING2)
plot(GENERATING2[,3],ESTIMATED2[,3],log="xy")

## MODEL 2b
foo<-function(x) c(cov2cor(x[[3]]$R[["1"]])[1,2],
	x[[3]]$R[["1"]][1,1],x[[3]]$R[["2"]][1,1],x[[3]]$R[["3"]][1,1],
	x[[3]]$R[["1"]][2,2])
ESTIMATED2b<-t(sapply(fits.model2b,foo))
dimnames(ESTIMATED2b)<-dimnames(GENERATING2b)
plot(GENERATING2b[,5],ESTIMATED2b[,5],log="xy")

## MODEL 2c
foo<-function(x) c(cov2cor(x[[4]]$R[["1"]])[1,2],
	x[[4]]$R[["1"]][1,1],
	x[[4]]$R[["1"]][2,2],x[[4]]$R[["2"]][2,2],x[[4]]$R[["3"]][2,2])
ESTIMATED2c<-t(sapply(fits.model2c,foo))
dimnames(ESTIMATED2c)<-dimnames(GENERATING2c)
plot(GENERATING2c[,2],ESTIMATED2c[,2],log="xy")

## MODEL 3
foo<-function(x) c(cov2cor(x[[5]]$R[["1"]])[1,2],cov2cor(x[[5]]$R[["2"]])[1,2],cov2cor(x[[5]]$R[["3"]])[1,2],
	x[[5]]$R[["1"]][1,1],
	x[[5]]$R[["1"]][2,2])
ESTIMATED3<-t(sapply(fits.model3,foo))
dimnames(ESTIMATED3)<-dimnames(GENERATING3)
plot(GENERATING3[,1],ESTIMATED3[,1])

## MODEL 3b
foo<-function(x) c(cov2cor(x[[6]]$R[["1"]])[1,2],cov2cor(x[[6]]$R[["2"]])[1,2],cov2cor(x[[6]]$R[["3"]])[1,2],
	x[[6]]$R[["1"]][1,1],x[[6]]$R[["2"]][1,1],x[[6]]$R[["3"]][1,1],
	x[[6]]$R[["1"]][2,2])
ESTIMATED3b<-t(sapply(fits.model3b,foo))
dimnames(ESTIMATED3b)<-dimnames(GENERATING3b)
plot(GENERATING3b[,1],ESTIMATED3b[,1])

## MODEL 3c
foo<-function(x) c(cov2cor(x[[7]]$R[["1"]])[1,2],cov2cor(x[[7]]$R[["2"]])[1,2],cov2cor(x[[7]]$R[["3"]])[1,2],
	x[[7]]$R[["1"]][1,1],
	x[[7]]$R[["1"]][2,2],x[[7]]$R[["2"]][2,2],x[[7]]$R[["3"]][2,2])
ESTIMATED3c<-t(sapply(fits.model3c,foo))
dimnames(ESTIMATED3c)<-dimnames(GENERATING3c)
plot(GENERATING3c[,1],ESTIMATED3c[,1])

## MODEL 4
foo<-function(x) c(cov2cor(x[[8]]$R[["1"]])[1,2],cov2cor(x[[8]]$R[["2"]])[1,2],cov2cor(x[[8]]$R[["3"]])[1,2],
	x[[8]]$R[["1"]][1,1],x[[8]]$R[["2"]][1,1],x[[8]]$R[["3"]][1,1],
	x[[8]]$R[["1"]][2,2],x[[8]]$R[["2"]][2,2],x[[8]]$R[["3"]][2,2])
ESTIMATED4<-t(sapply(fits.model4,foo))
dimnames(ESTIMATED4)<-dimnames(GENERATING4)
plot(GENERATING4[,1],ESTIMATED4[,1])

lnV<-function(X){
	ii<-grep("v",colnames(X))
	Y<-X
	Y[,ii]<-log(Y[,ii])
	Y
}

Table4<-matrix(NA,8,9,dimnames=list(names(fits.model1[[1]]),colnames(GENERATING4)))
Table4["model1",c(1,4,7)]<-diag(cor(lnV(GENERATING1),lnV(ESTIMATED1)))
Table4["model2",c(1,4:9)]<-diag(cor(lnV(GENERATING2),lnV(ESTIMATED2)))
Table4["model2b",c(1,4:6,7)]<-diag(cor(lnV(GENERATING2b),lnV(ESTIMATED2b)))
Table4["model2c",c(1,4,7:9)]<-diag(cor(lnV(GENERATING2c),lnV(ESTIMATED2c)))
Table4["model3",c(1:3,4,7)]<-diag(cor(lnV(GENERATING3),lnV(ESTIMATED3)))
Table4["model3b",c(1:3,4:6,7)]<-diag(cor(lnV(GENERATING3b),lnV(ESTIMATED3b)))
Table4["model3c",c(1:3,4,7:9)]<-diag(cor(lnV(GENERATING3c),lnV(ESTIMATED3c)))
Table4["model4",]<-diag(cor(lnV(GENERATING4),lnV(ESTIMATED4)))
 
xtable(Table4)

Table5<-matrix(NA,8,9,dimnames=list(names(fits.model1[[1]]),colnames(GENERATING4)))
Table5["model1",c(1,4,7)]<-colMeans(lnV(GENERATING1)-lnV(ESTIMATED1))
Table5["model2",c(1,4:9)]<-colMeans(lnV(GENERATING2)-lnV(ESTIMATED2))
Table5["model2b",c(1,4:6,7)]<-colMeans(lnV(GENERATING2b)-lnV(ESTIMATED2b))
Table5["model2c",c(1,4,7:9)]<-colMeans(lnV(GENERATING2c)-lnV(ESTIMATED2c))
Table5["model3",c(1:3,4,7)]<-colMeans(lnV(GENERATING3)-lnV(ESTIMATED3))
Table5["model3b",c(1:3,4:6,7)]<-colMeans(lnV(GENERATING3b)-lnV(ESTIMATED3b))
Table5["model3c",c(1:3,4,7:9)]<-colMeans(lnV(GENERATING3c)-lnV(ESTIMATED3c))
Table5["model4",]<-colMeans(lnV(GENERATING4)-lnV(ESTIMATED4))

xtable(Table5)

## save workspace
save.image(file="Figures.Rdata")
## end save workspace


## Supplementary Analysis

set.seed(53)

## load tropidurid data
trop.data<-read.csv(file="tropidurid-data.csv",row.names=1)
trop.tree<-read.simmap(file="tropidurid-tree.tre",format="nexus",
	version=1.5)
habitat<-getStates(trop.tree,"tips")

## model fits (setting root to non-rock)
fitER<-fitMk(trop.tree,habitat,model="ER",pi=c(1,0))
fitARD<-fitMk(trop.tree,habitat,model="ARD",pi=c(1,0))
model_01<-matrix(c(0,0,1,0),2,2)
fit01<-fitMk(trop.tree,habitat,model=model_01,pi=c(1,0))

TableS1<-cbind(
	sapply(list(fitER,fitARD,fit01),logLik),
	AIC(fitER,fitARD,fit01))
dimnames(TableS1)<-list(c("ER","ARD","directional"),
	c("logLik","df","AIC"))
xtable(TableS1)
## ARD best-supported model

## stochastic mapping
trop.trees<-make.simmap(trop.tree,habitat,model="ARD",pi=c(1,0),
	nsim=100)

## compute density map
dMap<-densityMap(trop.trees,plot=FALSE,res=500)
dMap<-setMap(dMap,c("white","black"))
## plot(dMap,outline=TRUE,direction="upwards")
tips<-dMap$tree$tip.label
dMap$tree$tip.label<-paste(" ",tips)

cols<-setNames(c("white","black"),c("n_rock","rock"))

## Figure A1
pdf(file="FigureA1.pdf",width=10,height=8)
par(mfrow=c(1,2))
plot(trop.trees[[25]],outline=TRUE,colors=cols,
	fsize=0.5,direction="rightwards",ftype="i",
	lwd=3,ylim=c(-8,76),mar=c(1.1,1.1,2.1,1.1))
hh<-legend(x=0,y=-2,c("non-rock-dwelling","rock-dwelling"),
	pch=22,pt.bg=c("white","black"),pt.cex=1.5,bty="n",cex=0.8)
mtext("a)",line=-1,adj=0)
plot(dMap,outline=TRUE,fsize=0.5,
	direction="rightwards",ftype="i",
	lwd=3,legend=FALSE,ylim=c(-8,76),
	mar=c(1.1,1.1,2.1,1.1))
add.color.bar(1,dMap$cols,title="PP(rock-dwelling)",
	prompt=FALSE,x=0,y=-6,fsize=0.8,subtitle="")
mtext("b)",line=-1,adj=0)
dev.off()
## end Figure A1

foo<-function(tree,data) evolvcv.lite(tree,data,models="all models")
TROP.FITS<-lapply(trop.trees,foo,data=trop.data)

foo<-function(x) sapply(x,function(x) x$AIC)
TROP.AIC<-sapply(TROP.FITS,foo)

## which ones to repeat
ij<-which(is.na(TROP.AIC),arr.ind=TRUE)

## repeat
while(length(ij)>0){
	for(i in 1:nrow(ij)){
		## MODEL<-strsplit(rownames(ij)[i],"model")[[1]][2]
		TROP.FITS[[ij[i,2]]][[ij[i,1]]]<-evolvcv.lite(trop.trees[[ij[i,2]]],
			trop.data,models="all models")[[ij[i,1]]]
	}
	foo<-function(x) sapply(x,function(x) x$AIC)
	TROP.AIC<-sapply(TROP.FITS,foo)
	ij<-which(is.na(TROP.AIC),arr.ind=TRUE)
}

TROP.W<-apply(TROP.AIC,2,aic.w)

rowMeans(TROP.W)

TROP.RANKED<-apply(TROP.AIC,2,function(x,nn) nn[order(x)],nn=rownames(TROP.AIC))

nn<-rownames(TROP.AIC)
BEST<-matrix(NA,length(nn),length(nn),dimnames=list(1:length(nn),nn))
for(i in 1:length(nn)){
	xx<-factor(TROP.RANKED[i,],levels=nn)
	BEST[i,]<-summary(xx)
}

BEST<-rbind(BEST/100,rowMeans(TROP.W))

xtable(BEST)

## save workspace
save.image(file="Figures.Rdata")
## end save workspace


## model average Model 3 across stochastic maps
aic.model3<-TROP.AIC["model3",]
aicw.model3<-aic.w(aic.model3)
model3.avefit<-TROP.FITS[[1]]$model3
model3.avefit$R$n_rock[]<-0
model3.avefit$R$rock[]<-0
model3.avefit$logLik<-0
model3.avefit$k<-0
model3.avefit$AIC<-0
for(i in 1:length(aicw.model3)){
	model3.avefit$R$n_rock<-model3.avefit$R$n_rock+aicw.model3[i]*
		TROP.FITS[[i]]$model3$R$n_rock
	model3.avefit$R$rock<-model3.avefit$R$rock+aicw.model3[i]*
		TROP.FITS[[i]]$model3$R$rock
	model3.avefit$logLik<-model3.avefit$logLik+aicw.model3[i]*
		TROP.FITS[[i]]$model3$logLik
	model3.avefit$k<-model3.avefit$k+aicw.model3[i]*TROP.FITS[[i]]$model3$k
	model3.avefit$AIC<-model3.avefit$AIC+aicw.model3[i]*
		TROP.FITS[[i]]$model3$AIC
}
Model3<-list(model3=model3.avefit)
class(Model3)<-"evolvcv.lite"

## model average Model 3c across stochastic maps
aic.model3c<-TROP.AIC["model3c",]
aicw.model3c<-aic.w(aic.model3c)
model3c.avefit<-TROP.FITS[[1]]$model3c
model3c.avefit$R$n_rock[]<-0
model3c.avefit$R$rock[]<-0
model3c.avefit$logLik<-0
model3c.avefit$k<-0
model3c.avefit$AIC<-0
for(i in 1:length(aicw.model3c)){
	model3c.avefit$R$n_rock<-model3c.avefit$R$n_rock+aicw.model3c[i]*
		TROP.FITS[[i]]$model3c$R$n_rock
	model3c.avefit$R$rock<-model3c.avefit$R$rock+aicw.model3c[i]*
		TROP.FITS[[i]]$model3c$R$rock
	model3c.avefit$logLik<-model3c.avefit$logLik+aicw.model3c[i]*
		TROP.FITS[[i]]$model3c$logLik
	model3c.avefit$k<-model3c.avefit$k+aicw.model3c[i]*TROP.FITS[[i]]$model3c$k
	model3c.avefit$AIC<-model3c.avefit$AIC+aicw.model3c[i]*
		TROP.FITS[[i]]$model3c$AIC
}
Model3c<-list(model3c=model3c.avefit)
class(Model3c)<-"evolvcv.lite"

## model average Model 3b across stochastic maps
aic.model3b<-TROP.AIC["model3b",]
aicw.model3b<-aic.w(aic.model3b)
model3b.avefit<-TROP.FITS[[1]]$model3b
model3b.avefit$R$n_rock[]<-0
model3b.avefit$R$rock[]<-0
model3b.avefit$logLik<-0
model3b.avefit$k<-0
model3b.avefit$AIC<-0
for(i in 1:length(aicw.model3b)){
	model3b.avefit$R$n_rock<-model3b.avefit$R$n_rock+aicw.model3b[i]*
		TROP.FITS[[i]]$model3b$R$n_rock
	model3b.avefit$R$rock<-model3b.avefit$R$rock+aicw.model3b[i]*
		TROP.FITS[[i]]$model3b$R$rock
	model3b.avefit$logLik<-model3b.avefit$logLik+aicw.model3b[i]*
		TROP.FITS[[i]]$model3b$logLik
	model3b.avefit$k<-model3b.avefit$k+aicw.model3b[i]*TROP.FITS[[i]]$model3b$k
	model3b.avefit$AIC<-model3b.avefit$AIC+aicw.model3b[i]*
		TROP.FITS[[i]]$model3b$AIC
}
Model3b<-list(model3b=model3b.avefit)
class(Model3b)<-"evolvcv.lite"

## model average Model 4 across stochastic maps
aic.model4<-TROP.AIC["model4",]
aicw.model4<-aic.w(aic.model4)
model4.avefit<-TROP.FITS[[1]]$model4
model4.avefit$R$n_rock[]<-0
model4.avefit$R$rock[]<-0
model4.avefit$logLik<-0
model4.avefit$k<-0
model4.avefit$AIC<-0
for(i in 1:length(aicw.model4)){
	model4.avefit$R$n_rock<-model4.avefit$R$n_rock+aicw.model4[i]*
		TROP.FITS[[i]]$model4$R$n_rock
	model4.avefit$R$rock<-model4.avefit$R$rock+aicw.model4[i]*
		TROP.FITS[[i]]$model4$R$rock
	model4.avefit$logLik<-model4.avefit$logLik+aicw.model4[i]*
		TROP.FITS[[i]]$model4$logLik
	model4.avefit$k<-model4.avefit$k+aicw.model4[i]*TROP.FITS[[i]]$model4$k
	model4.avefit$AIC<-model4.avefit$AIC+aicw.model4[i]*
		TROP.FITS[[i]]$model4$AIC
}
Model4<-list(model4=model4.avefit)
class(Model4)<-"evolvcv.lite"

Models<-c(Model3,Model3c,Model3b,Model4)
class(Models)<-"evolvcv.lite"
Models

Description<-paste(sapply(Models,function(x) x$description),names(Models))

TROP.AVERAGED<-data.frame(description=Description,
	sig11=sapply(Models,function(x) x$R[[1]][1,1]),
	sig12=sapply(Models,function(x) x$R[[2]][1,1]),
	sig21=sapply(Models,function(x) x$R[[1]][2,2]),
	sig22=sapply(Models,function(x) x$R[[2]][2,2]),
	r1=sapply(Models,function(x) cov2cor(x$R[[1]])[1,2]),
	r2=sapply(Models,function(x) cov2cor(x$R[[2]])[1,2]))
rownames(TROP.AVERAGED)<-NULL
xtable(TROP.AVERAGED)

## save workspace
save.image(file="Figures.Rdata")
## end save workspace

