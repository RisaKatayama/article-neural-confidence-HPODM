library(coin)
library(pracma)
library(dplyr)
library(extrafont)
font_import()
fonts()
library(loo)
library(runjags)
library(heplots)
library(hesim)
library(grDevices)

.mod.env = new.env()
sys.source("modules.R", envir=.mod.env)
attach(.mod.env)

## data loading ####
dat.folder = ""
mdlres.folder = ""
mdl.folder = "model/"
load(paste0(c(dat.folder,"bhvdata_bhv.RData"),collapse="/"))

load(paste0(c(mdlres.folder,"cwHDM_dimdl.RData"),collapse="/"))
mdlres = res
mdlres$dk.summaries = add.summary(mdlres$dk.model)
mdlres$dik.summaries = add.summary(mdlres$di.model)

exsbj = c(7,22) # deckest last quat < 0.5

ns.org = max(Bhvdata$sidx)
insbj = c(1:ns.org)
insbj = insbj[!is.element(insbj,exsbj)]
Bhvdata = Bhvdata[!is.element(Bhvdata$sidx,exsbj),]
Bhvdata$sidx.org = Bhvdata$sidx
ns = length(insbj)
for (i in 1:ns) {
   Bhvdata$sidx[Bhvdata$sidx.org==insbj[i]] = i
}

load.data = function(data){
  ns = max(data$sidx)
  nses = max(data$ses)
  nt = 16
  N = dim(data)[1]
  
  dkest = (data$deckest==6)*1  # deck 6 = 1
  dktrue = (data$decktrue==6)*1
  ev4.rel = c()
  ev6.rel = c()
  m.add = c()
  for (i in 1:ns) {
    for (s in 1:nses) {
      tmpadd = data$add[data$sidx==i&data$ses==s]
      
      z = dbinom(tmpadd,10,0.4)/(dbinom(tmpadd,10,0.4)+dbinom(tmpadd,10,0.6))
      x = c(0,(tmpadd<5)*z)
      tmpev4.rel = mapply(function(t,x){sum(x[1:t])}, 1:nt, MoreArgs=list(x=x))
      
      z = dbinom(tmpadd,10,0.6)/(dbinom(tmpadd,10,0.4)+dbinom(tmpadd,10,0.6))
      x = c(0,(tmpadd>5)*z)
      tmpev6.rel = mapply(function(t,x){sum(x[1:t])}, 1:nt, MoreArgs=list(x=x))
      
      ev4.rel = c(ev4.rel,tmpev4.rel)
      ev6.rel = c(ev6.rel,tmpev6.rel)
    }
  }

  vr = c()
  for (i in 1:ns) {
    for (s in 1:nses) {
      tmpdkest = dkest[data$sidx==i&data$ses==s]
      tmpdkest0 = (tmpdkest==0)*1
      tmpdkest1 = (tmpdkest==1)*1
      vr = c(vr, (sapply(1:16,function(j){nansum(tmpdkest0[1:j])})*tmpdkest0+sapply(1:16,function(j){nansum(tmpdkest1[1:j])})*tmpdkest1)/c(1:16))
    }
  }

  vh.cdk = 21 - (data$open + (4*(1-dkest) + 6*dkest) + (data$p.us + (1-data$p.us)*4))
  vh.udk = 21 - (data$open + (4*dkest + (1-dkest)*6) + (data$p.us + (1-data$p.us)*4))
  vh.add = 21 - (data$open + m.add + (data$p.us + (1-data$p.us)*4))

  udc = (data$p.us==0.5)*1  # face-down deck H = 1, deck L = 0
  
  include = !is.na(data$deckest*data$cf.deck*data$decision*data$cf.deci)
  res = data.frame(open=data$open[include], add=data$add[include], dkest=dkest[include], dktrue=dktrue[include], kcf=data$cf.deck[include],
                   ev4.rel=ev4.rel[include], ev6.rel=ev6.rel[include], vr=vr[include],
                   icf=data$cf.deci[include], deci=data$decision[include], vh.cdk=vh.cdk[include], vh.udk=vh.udk[include], vh.add=vh.add[include],
                   udc=udc[include], sidx=data$sidx[include], sidx.org=data$sidx.org[include])
  
  res$evd.rea = res$ev6.rel - res$ev4.rel
  
  return(res)
}

bdata = load.data(Bhvdata)

tmp.bdata = bdata[bdata$sidx==1,]
simns = 80
simN = dim(tmp.bdata)[1]*simns

inits1 = dump.format(list(.RNG.name="base::Super-Duper", .RNG.seed=111))
inits2 = dump.format(list(.RNG.name="base::Wichmann-Hill", .RNG.seed=222))
inits3 = dump.format(list(.RNG.name="base::Mersenne-Twister", .RNG.seed=333))


sim.dat = simulate.bhv(mdlres, simns, tmp.bdata)

dat1 = dump.format(list(dkest=sim.dat$dat$dkest, de=sim.dat$dat$evd.rea, vr=sim.dat$dat$vr, kcf=sim.dat$dat$kcf,
                        sidx=sim.dat$dat$sidx, N=simN, ns=simns))
vals = c("b.dk.mu","b.kcf.mu","a0.kcf.mu","a1.kcf.mu","b.dk.sig","b.kcf.sig","a0.kcf.sig",
         "b0.dk.p","b1.dk.p","b1.kcf.p","b2.kcf.p","a0.kcf.p")

sim.rea.dkmdl = run.jags(model=paste0(mdl.folder,"basic_dkmdl.txt"), monitor=vals, data=dat1, n.chains=3, inits=c(inits1,inits2,inits3),
                         plots=FALSE, method="parallel", burnin=10000, sample=1000, adapt=1000, thin=20)

summary.sim.rea.dkmdl = add.summary(sim.rea.dkmdl)

tmpbdata = data.frame(dkest=sim.dat$dat$dkest, de=sim.dat$dat$evd.rea, vr=sim.dat$dat$vr, kcf=sim.dat$dat$kcf, sidx=sim.dat$dat$sidx)

b.kcf = array(NaN, dim=c(simN,2)) 
theta.kcf = array(NaN, dim=c(simN,3))
for (i in 1:simns) {
  b.kcf[tmpbdata$sidx==i,1] = summary.sim.rea.dkmdl$summaries[paste0("b1.kcf.p[",i,"]"),"Median"]
  b.kcf[tmpbdata$sidx==i,2] = summary.sim.rea.dkmdl$summaries[paste0("b2.kcf.p[",i,"]"),"Median"]

  a0 = summary.sim.rea.dkmdl$summaries[paste0("a0.kcf.p[",i,"]"),"Median"]
  a1 = summary.sim.rea.dkmdl$summaries["a1.kcf.mu","Median"]
  theta.kcf[tmpbdata$sidx==i,1] = a0-a1
  theta.kcf[tmpbdata$sidx==i,2] = a0
  theta.kcf[tmpbdata$sidx==i,3] = a0+a1
}

X = cbind(abs(tmpbdata$de), tmpbdata$vr)
z.kcf = apply(X*b.kcf, 1, sum)
p.kcf.mat = array(NaN,dim=c(simN,4))
p.kcf.mat[,1] = 1 - ilogit(z.kcf+theta.kcf[,3])
p.kcf.mat[,2] = ilogit(z.kcf+theta.kcf[,3]) - ilogit(z.kcf+theta.kcf[,2])
p.kcf.mat[,3] = ilogit(z.kcf+theta.kcf[,2]) - ilogit(z.kcf+theta.kcf[,1])
p.kcf.mat[,4] = ilogit(z.kcf+theta.kcf[,1])
w.pred.kcf = apply(sweep(p.kcf.mat,2,1:4,"*"),1,sum)

ud = 3.4*(1-sim.dat$dat$udc) + 2.5*sim.dat$dat$udc
cdk = 4*(1-sim.dat$dat$dkest) + 6*sim.dat$dat$dkest
udk = 6*(1-sim.dat$dat$dkest) + 4*sim.dat$dat$dkest
dat1 = dump.format(list(open=sim.dat$dat$open, ud=ud, cdk=cdk, udk=udk, udc=udc,
                        kcf=(w.pred.kcf-1)/3, deci=sim.dat$dat$deci, icf.c=sim.dat$dat$icf, sidx=sim.dat$dat$sidx, N=simN, ns=simns))
vals = c("b.di.mu","b.icf.mu","a0.icf.mu","a1.icf.mu","b.di.sig","b.icf.sig","a0.icf.sig",
         "g.kcf.p","e.kcf.p","b0.di.p","b1.di.p","b2.di.p","b1.icf.p","b2.icf.p","a0.icf.p")

sim.cwhdm.dimdl = run.jags(model=model=paste0(mdl.folder,"cwHDM_dimdl.txt"), monitor=vals, data=dat1, n.chains=3, inits=c(inits1,inits2,inits3),
                           plots=FALSE, method="parallel", burnin=10000, sample=1000, adapt=1000, thin=20)

summary.sim.cwhdm.dimdl = add.summary(sim.cwhdm.dimdl)

## Fig S6a
pnames.dk = rownames(summary.sim.rea.dkmdl$summaries)
post.dk = rbind(sim.rea.dkmdl$mcmc[[1]],sim.rea.dkmdl$mcmc[[2]],sim.rea.dkmdl$mcmc[[3]])
hdi.post.dk = hdi(sim.rea.dkmdl)

scott.bw = function(x) {3.5*sd(x)/(length(x)^(1/3))}

dat.figs6.a1 = data.frame(val=post.dk[,pnames.dk=="b.dk.mu[1]"])
figs6.a1 = ggplot(data=dat.figs6.a1, aes(x=val)) +
          geom_histogram(fill="gray", alpha=0.7, binwidth=scott.bw) +
          geom_vline(xintercept=mdlres$dk.summaries$summaries["b.dk.mu[1]","Mean"]) +
          geom_vline(xintercept=hdi.post.dk[1,"b.dk.mu[1]"],lty="dotted") +
          geom_vline(xintercept=hdi.post.dk[2,"b.dk.mu[1]"],lty="dotted") +
          scale_x_continuous(limits=c(-0.4,0.2), name=NULL) +
          scale_y_continuous(name=NULL) +
          labs(title=expression(beta[dk]*"("*bias*")")) +
          theme_custom +
          theme(plot.title=element_text(hjust=0.5))
figs6.a1

dat.figs6.a2 = data.frame(val=post.dk[,pnames.dk=="b.dk.mu[2]"])
figs6.a2 = ggplot(data=dat.figs6.a2, aes(x=val)) +
          geom_histogram(fill="gray", alpha=0.7, binwidth=scott.bw) +
          geom_vline(xintercept=mdlres$dk.summaries$summaries["b.dk.mu[2]","Mean"]) +
          geom_vline(xintercept=hdi.post.dk[1,"b.dk.mu[2]"],lty="dotted") +
          geom_vline(xintercept=hdi.post.dk[2,"b.dk.mu[2]"],lty="dotted") +
          scale_x_continuous(expand = c(0,0.1), name=NULL) +
          scale_y_continuous(name=NULL) +
          labs(title=expression(beta[dk]*"("*Delta*E*")")) +
          theme_custom +
          theme(plot.title=element_text(hjust=0.5))
figs6.a2

dat.figs6.a3 = data.frame(val=post.dk[,pnames.dk=="b.kcf.mu[1]"])
figs6.a3 = ggplot(data=dat.figs6.a3, aes(x=val)) +
          geom_histogram(fill="gray", alpha=0.7, binwidth=scott.bw) +
          geom_vline(xintercept=mdlres$dk.summaries$summaries["b.kcf.mu[1]","Mean"]) +
          geom_vline(xintercept=hdi.post.dk[1,"b.kcf.mu[1]"],lty="dotted") +
          geom_vline(xintercept=hdi.post.dk[2,"b.kcf.mu[1]"],lty="dotted") +
          scale_x_continuous(limits=c(0.7,1.2), name=NULL) +
          scale_y_continuous(name=NULL) +
          labs(title=expression(beta[c[k]]*"(|"*Delta*E*"|)")) +
          theme_custom +
          theme(plot.title=element_text(hjust=0.5))
figs6.a3

dat.figs6.a4 = data.frame(val=post.dk[,pnames.dk=="b.kcf.mu[2]"])
figs6.a4 = ggplot(data=dat.figs6.a4, aes(x=val)) +
          geom_histogram(fill="gray", alpha=0.7, binwidth=scott.bw) +
          geom_vline(xintercept=mdlres$dk.summaries$summaries["b.kcf.mu[2]","Mean"]) +
          geom_vline(xintercept=hdi.post.dk[1,"b.kcf.mu[2]"],lty="dotted") +
          geom_vline(xintercept=hdi.post.dk[2,"b.kcf.mu[2]"],lty="dotted") +
          scale_x_continuous(limits=c(0,1.5), name=NULL) +
          scale_y_continuous(name=NULL) +
          labs(title=expression(beta[c[k]]*"(|"*s*"|)")) +
          theme_custom +
          theme(plot.title=element_text(hjust=0.5))
figs6.a4

dat.figs6.a5 = data.frame(val=post.dk[,pnames.dk=="a0.kcf.mu"])
figs6.a5 = ggplot(data=dat.figs6.a5, aes(x=val)) +
          geom_histogram(fill="gray", alpha=0.7, binwidth=scott.bw) +
          geom_vline(xintercept=mdlres$dk.summaries$summaries["a0.kcf.mu","Mean"]) +
          geom_vline(xintercept=hdi.post.dk[1,"a0.kcf.mu"],lty="dotted") +
          geom_vline(xintercept=hdi.post.dk[2,"a0.kcf.mu"],lty="dotted") +
          scale_x_continuous(limits=c(-3,-0.5), name=NULL) +
          scale_y_continuous(name=NULL) +
          labs(title=expression(theta[c[k]*","*th])) +
          theme_custom +
          theme(plot.title=element_text(hjust=0.5))
figs6.a5

dat.figs6.a6 = data.frame(val=post.dk[,pnames.dk=="a1.kcf.mu"])
figs6.a6 = ggplot(data=dat.figs6.a6, aes(x=val)) +
          geom_histogram(fill="gray", alpha=0.7, binwidth=scott.bw) +
          geom_vline(xintercept=mdlres$dk.summaries$summaries["a1.kcf.mu","Mean"]) +
          geom_vline(xintercept=hdi.post.dk[1,"a1.kcf.mu"],lty="dotted") +
          geom_vline(xintercept=hdi.post.dk[2,"a1.kcf.mu"],lty="dotted") +
          scale_x_continuous(limits=c(1.8,2), name=NULL) +
          scale_y_continuous(name=NULL) +
          labs(title=expression(theta[c[k]*","*sp])) +
          theme_custom +
          theme(plot.title=element_text(hjust=0.5))
figs6.a6

pnames.di = rownames(summary.summary.sim.cwhdm.dimdl$summaries)
post.di = rbind(summary.sim.cwhdm.dimdl$mcmc[[1]],summary.sim.cwhdm.dimdl$mcmc[[2]],summary.sim.cwhdm.dimdl$mcmc[[3]])
hdi.post.di = hdi(summary.sim.cwhdm.dimdl)

dat.figs6.a7 = data.frame(val=post.di[,pnames.di=="b.di.mu[1]"])
figs6.a7 = ggplot(data=dat.figs6.a7, aes(x=val)) +
          geom_histogram(fill="gray", alpha=0.7, binwidth=scott.bw) +
          geom_vline(xintercept=mdlres$di.summaries$summaries["b.di.mu[1]","Mean"]) +
          geom_vline(xintercept=hdi.post.di[1,"b.di.mu[1]"],lty="dotted") +
          geom_vline(xintercept=hdi.post.di[2,"b.di.mu[1]"],lty="dotted") +
          scale_x_continuous(limits=c(-1,-0.2), name=NULL) +
          scale_y_continuous(name=NULL) +
          labs(title=expression(beta[di]*"("*bias*")")) +
          theme_custom +
          theme(plot.title=element_text(hjust=0.5))
figs6.a7

dat.figs6.a8 = data.frame(val=post.di[,pnames.di=="b.di.mu[2]"])
figs6.a8 = ggplot(data=dat.figs6.a8, aes(x=val)) +
          geom_histogram(fill="gray", alpha=0.7, binwidth=scott.bw) +
          geom_vline(xintercept=mdlres$di.summaries$summaries["b.di.mu[2]","Mean"]) +
          geom_vline(xintercept=hdi.post.di[1,"b.di.mu[2]"],lty="dotted") +
          geom_vline(xintercept=hdi.post.di[2,"b.di.mu[2]"],lty="dotted") +
          scale_x_continuous(limits=c(1.1,1.6), name=NULL) +
          scale_y_continuous(name=NULL) +
          labs(title=expression(beta[di]*"("*V[hit]*")")) +
          theme_custom +
          theme(plot.title=element_text(hjust=0.5))
figs6.a8

dat.figs6.a9 = data.frame(val=post.di[,pnames.di=="b.di.mu[3]"])
figs6.a9 = ggplot(data=dat.figs6.a9, aes(x=val)) +
          geom_histogram(fill="gray", alpha=0.7, binwidth=scott.bw) +
          geom_vline(xintercept=mdlres$di.summaries$summaries["b.di.mu[3]","Mean"]) +
          geom_vline(xintercept=hdi.post.di[1,"b.di.mu[3]"],lty="dotted") +
          geom_vline(xintercept=hdi.post.di[2,"b.di.mu[3]"],lty="dotted") +
          scale_x_continuous(limits=c(-0.05,0.1), name=NULL) +
          scale_y_continuous(name=NULL) +
          labs(title=expression(beta[di]*"("*V[hit]*"|"*V[hit]*"|"*")")) +
          theme_custom +
          theme(plot.title=element_text(hjust=0.5))
figs6.a9

dat.figs6.a10 = data.frame(val=post.di[,pnames.di=="b.icf.mu[1]"])
figs6.a10 = ggplot(data=dat.figs6.a10, aes(x=val)) +
          geom_histogram(fill="gray", alpha=0.7, binwidth=scott.bw) +
          geom_vline(xintercept=mdlres$di.summaries$summaries["b.icf.mu[1]","Mean"]) +
          geom_vline(xintercept=hdi.post.di[1,"b.icf.mu[1]"],lty="dotted") +
          geom_vline(xintercept=hdi.post.di[2,"b.icf.mu[1]"],lty="dotted") +
          scale_x_continuous(limits=c(0.3,0.7), name=NULL) +
          scale_y_continuous(name=NULL) +
          labs(title=expression(beta[c[i]]*"(|"*SV*"|)")) +
          theme_custom +
          theme(plot.title=element_text(hjust=0.5))
figs6.a10

dat.figs6.a11 = data.frame(val=post.di[,pnames.di=="b.icf.mu[2]"])
figs6.a11 = ggplot(data=dat.figs6.a11, aes(x=val)) +
          geom_histogram(fill="gray", alpha=0.7, binwidth=scott.bw) +
          geom_vline(xintercept=mdlres$di.summaries$summaries["b.icf.mu[2]","Mean"]) +
          geom_vline(xintercept=hdi.post.di[1,"b.icf.mu[2]"],lty="dotted") +
          geom_vline(xintercept=hdi.post.di[2,"b.icf.mu[2]"],lty="dotted") +
          scale_x_continuous(limits=c(0.5,2), name=NULL) +
          scale_y_continuous(name=NULL) +
          labs(title=expression(beta[c[i]]*"("*c[k]*")")) +
          theme_custom +
          theme(plot.title=element_text(hjust=0.5))
figs6.a11

dat.figs6.a12 = data.frame(val=post.di[,pnames.di=="a0.icf.mu"])
figs6.a12 = ggplot(data=dat.figs6.a12, aes(x=val)) +
          geom_histogram(fill="gray", alpha=0.7, binwidth=scott.bw) +
          geom_vline(xintercept=mdlres$di.summaries$summaries["a0.icf.mu","Mean"]) +
          geom_vline(xintercept=hdi.post.di[1,"a0.icf.mu"],lty="dotted") +
          geom_vline(xintercept=hdi.post.di[2,"a0.icf.mu"],lty="dotted") +
          scale_x_continuous(limits=c(-2.5,-1), name=NULL) +
          scale_y_continuous(name=NULL) +
          labs(title=expression(theta[c[i]*","*th])) +
          theme_custom +
          theme(plot.title=element_text(hjust=0.5))
figs6.a12

dat.figs6.a13 = data.frame(val=post.di[,pnames.di=="a1.icf.mu"])
figs6.a13 = ggplot(data=dat.figs6.a13, aes(x=val)) +
          geom_histogram(fill="gray", alpha=0.7, binwidth=scott.bw) +
          geom_vline(xintercept=mdlres$di.summaries$summaries["a1.icf.mu","Mean"]) +
          geom_vline(xintercept=hdi.post.di[1,"a1.icf.mu"],lty="dotted") +
          geom_vline(xintercept=hdi.post.di[2,"a1.icf.mu"],lty="dotted") +
          scale_x_continuous(limits=c(2,2.3), name=NULL) +
          scale_y_continuous(name=NULL) +
          labs(title=expression(theta[c[i]*","*sp])) +
          theme_custom +
          theme(plot.title=element_text(hjust=0.5))
figs6.a13

## Fig S6b
dat.figs6.b1 = data.frame(x=sim.dat$prm$b.dk[,1],
                          y=sapply(1:simns,function(i){summary.sim.rea.dkmdl$summaries[pnames.dk==paste0("b0.dk.p[",i,"]"),"Median"]}))
figs6.b1 = ggplot(data=dat.figs6.b1, aes(x=x, y=y)) +
          geom_line(data=data.frame(x=c(-1.5,1.5), y=c(-1.5,1.5)), linewidth=0.25) +
          geom_point(shape=16, size=1.5, col="gray") +
          scale_x_continuous(limits=c(-1.5,1.5), name=NULL) +
          scale_y_continuous(name=NULL) +
          labs(title=expression(beta[dk]*"("*bias*")")) +
          theme_custom +
          theme(plot.title=element_text(hjust=0.5))
figs6.b1

dat.figs6.b2 = data.frame(x=sim.dat$prm$b.dk[,2],
                          y=sapply(1:simns,function(i){summary.sim.rea.dkmdl$summaries[pnames.dk==paste0("b1.dk.p[",i,"]"),"Median"]}))
figs6.b2 = ggplot(data=dat.figs6.b2, aes(x=x, y=y)) +
          geom_line(data=data.frame(x=c(0,3), y=c(0,3)), linewidth=0.25) +
          geom_point(shape=16, size=1.5, col="gray") +
          scale_x_continuous(limits=c(0,3), name=NULL, labels=scales::number_format(accuracy=0.1)) +
          scale_y_continuous(name=NULL, labels=scales::number_format(accuracy=0.1)) +
          labs(title=expression(beta[dk]*"("*Delta*E*")")) +
          theme_custom +
          theme(plot.title=element_text(hjust=0.5))
figs6.b2

dat.figs6.b3 = data.frame(x=sim.dat$prm$b.kcf[,1],
                          y=sapply(1:simns,function(i){summary.sim.rea.dkmdl$summaries[pnames.dk==paste0("b1.kcf.p[",i,"]"),"Median"]}))
figs6.b3 = ggplot(data=dat.figs6.b3, aes(x=x, y=y)) +
          geom_line(data=data.frame(x=c(-1,3), y=c(-1,3)), linewidth=0.25) +
          geom_point(shape=16, size=1.5, col="gray") +
          scale_x_continuous(limits=c(-1,3), name=NULL, labels=scales::number_format(accuracy=0.1)) +
          scale_y_continuous(name=NULL, labels=scales::number_format(accuracy=0.1)) +
          labs(title=expression(beta[c[k]]*"(|"*Delta*E*"|)")) +
          theme_custom +
          theme(plot.title=element_text(hjust=0.5))
figs6.b3

dat.figs6.b4 = data.frame(x=sim.dat$prm$b.kcf[,2],
                          y=sapply(1:simns,function(i){summary.sim.rea.dkmdl$summaries[pnames.dk==paste0("b2.kcf.p[",i,"]"),"Median"]}))
figs6.b4 = ggplot(data=dat.figs6.b4, aes(x=x, y=y)) +
          geom_line(data=data.frame(x=c(-4,4), y=c(-4,4)), linewidth=0.25) +
          geom_point(shape=16, size=1.5, col="gray") +
          scale_x_continuous(limits=c(-4,4), name=NULL, labels=scales::number_format(accuracy=0.1)) +
          scale_y_continuous(name=NULL, labels=scales::number_format(accuracy=0.1)) +
          labs(title=expression(beta[c[k]]*"("*s*")")) +
          theme_custom +
          theme(plot.title=element_text(hjust=0.5))
figs6.b4

dat.figs6.b5 = data.frame(x=sim.dat$prm$a0.kcf,
                          y=sapply(1:simns,function(i){summary.sim.rea.dkmdl$summaries[pnames.dk==paste0("a0.kcf.p[",i,"]"),"Median"]}))
figs6.b5 = ggplot(data=dat.figs6.b5, aes(x=x, y=y)) +
          geom_line(data=data.frame(x=c(-7.5,5), y=c(-7.5,5)), linewidth=0.25) +
          geom_point(shape=16, size=1.5, col="gray") +
          scale_x_continuous(limits=c(-7.5,5), name=NULL, labels=scales::number_format(accuracy=0.1)) +
          scale_y_continuous(name=NULL, labels=scales::number_format(accuracy=0.1)) +
          labs(title=expression(theta[c[k]*","*th])) +
          theme_custom +
          theme(plot.title=element_text(hjust=0.5))
figs6.b5

dat.figs6.b6 = data.frame(x=sim.dat$prm$b.di[,1],
                          y=sapply(1:simns,function(i){summary.sim.cwhdm.dimdl$summaries[pnames.di==paste0("b0.di.p[",i,"]"),"Median"]}))
figs6.b6 = ggplot(data=dat.figs6.b6, aes(x=x, y=y)) +
          geom_line(data=data.frame(x=c(-2,1), y=c(-2,1)), linewidth=0.25) +
          geom_point(shape=16, size=1.5, col="gray") +
          scale_x_continuous(limits=c(-2,1), name=NULL, labels=scales::number_format(accuracy=0.1)) +
          scale_y_continuous(name=NULL, labels=scales::number_format(accuracy=0.1)) +
          labs(title=expression(beta[di]*"("*bias*")")) +
          theme_custom +
          theme(plot.title=element_text(hjust=0.5))
figs6.b6

dat.figs6.b7 = data.frame(x=sim.dat$prm$b.di[,2],
                          y=sapply(1:simns,function(i){summary.sim.cwhdm.dimdl$summaries[pnames.di==paste0("b1.di.p[",i,"]"),"Median"]}))
figs6.b7 = ggplot(data=dat.figs6.b7, aes(x=x, y=y)) +
          geom_line(data=data.frame(x=c(0.5,2.5), y=c(0.5,2.5)), linewidth=0.25) +
          geom_point(shape=16, size=1.5, col="gray") +
          scale_x_continuous(limits=c(0.5,2.5), name=NULL) +
          scale_y_continuous(name=NULL) +
          labs(title=expression(beta[di]*"("*V[hit]*")")) +
          theme_custom +
          theme(plot.title=element_text(hjust=0.5))
figs6.b7

dat.figs6.b8 = data.frame(x=sim.dat$prm$b.di[,3],
                          y=sapply(1:simns,function(i){summary.sim.cwhdm.dimdl$summaries[pnames.di==paste0("b2.di.p[",i,"]"),"Median"]}))
figs6.b8 = ggplot(data=dat.figs6.b8, aes(x=x, y=y)) +
          geom_line(data=data.frame(x=c(-0.2,0.4), y=c(-0.2,0.4)), linewidth=0.25) +
          geom_point(shape=16, size=1.5, col="gray") +
          scale_x_continuous(limits=c(-0.2,0.4), name=NULL) +
          scale_y_continuous(name=NULL) +
          labs(title=expression(beta[di]*"("*V[hit]*"|"*V[hit]*"|)")) +
          theme_custom +
          theme(plot.title=element_text(hjust=0.5))
figs6.b8

dat.figs6.b9 = data.frame(x=sim.dat$prm$g.kcf,
                          y=sapply(1:simns,function(i){summary.sim.cwhdm.dimdl$summaries[pnames.di==paste0("g.kcf.p[",i,"]"),"Median"]}))
figs6.b9 = ggplot(data=dat.figs6.b9, aes(x=x, y=y)) +
          geom_line(data=data.frame(x=c(0,0.8), y=c(0,0.8)), linewidth=0.25) +
          geom_point(shape=16, size=1.5, col="gray") +
          scale_x_continuous(limits=c(0,0.8), name=NULL) +
          scale_y_continuous(name=NULL) +
          labs(title=expression(gamma)) +
          theme_custom +
          theme(plot.title=element_text(hjust=0.5))
figs6.b9

dat.figs6.b10 = data.frame(x=sim.dat$prm$e.kcf,
                          y=sapply(1:simns,function(i){summary.sim.cwhdm.dimdl$summaries[pnames.di==paste0("e.kcf.p[",i,"]"),"Median"]}))
figs6.b10 = ggplot(data=dat.figs6.b10, aes(x=x, y=y)) +
          geom_line(data=data.frame(x=c(0.5,1), y=c(0.5,1)), linewidth=0.25) +
          geom_point(shape=16, size=1.5, col="gray") +
          scale_x_continuous(limits=c(0.5,1), name=NULL) +
          scale_y_continuous(name=NULL) +
          labs(title=expression(epsilon)) +
          theme_custom +
          theme(plot.title=element_text(hjust=0.5))
figs6.b10

dat.figs6.b11 = data.frame(x=sim.dat$prm$b.icf[,1],
                          y=sapply(1:simns,function(i){summary.sim.cwhdm.dimdl$summaries[pnames.di==paste0("b1.icf.p[",i,"]"),"Median"]}))
figs6.b11 = ggplot(data=dat.figs6.b11, aes(x=x, y=y)) +
          geom_line(data=data.frame(x=c(-0.5,2), y=c(-0.5,2)), linewidth=0.25) +
          geom_point(shape=16, size=1.5, col="gray") +
          scale_x_continuous(limits=c(-0.5,2), name=NULL, labels=scales::number_format(accuracy=0.1)) +
          scale_y_continuous(name=NULL, labels=scales::number_format(accuracy=0.1)) +
          labs(title=expression(beta[c[i]]*"(|"*SV*"|)")) +
          theme_custom +
          theme(plot.title=element_text(hjust=0.5))
figs6.b11

dat.figs6.b12 = data.frame(x=sim.dat$prm$b.icf[,2],
                          y=sapply(1:simns,function(i){summary.sim.cwhdm.dimdl$summaries[pnames.di==paste0("b2.icf.p[",i,"]"),"Median"]}))
figs6.b12 = ggplot(data=dat.figs6.b12, aes(x=x, y=y)) +
          geom_line(data=data.frame(x=c(-2,6.5), y=c(-2,6.5)), linewidth=0.25) +
          geom_point(shape=16, size=1.5, col="gray") +
          scale_x_continuous(limits=c(-2,6.5), name=NULL, labels=scales::number_format(accuracy=0.1)) +
          scale_y_continuous(name=NULL, labels=scales::number_format(accuracy=0.1)) +
          labs(title=expression(beta[c[i]]*"("*c[k]*")")) +
          theme_custom +
          theme(plot.title=element_text(hjust=0.5))
figs6.b12

dat.figs6.b13 = data.frame(x=sim.dat$prm$a0.icf,
                          y=sapply(1:simns,function(i){summary.sim.cwhdm.dimdl$summaries[pnames.di==paste0("a0.icf.p[",i,"]"),"Median"]}))
figs6.b13 = ggplot(data=dat.figs6.b13, aes(x=x, y=y)) +
          geom_line(data=data.frame(x=c(-7.5,5), y=c(-7.5,5)), linewidth=0.25) +
          geom_point(shape=16, size=1.5, col="gray") +
          scale_x_continuous(limits=c(-7.5,5), name=NULL, labels=scales::number_format(accuracy=0.1)) +
          scale_y_continuous(name=NULL, labels=scales::number_format(accuracy=0.1)) +
          labs(title=expression(theta[c[i]*","*th])) +
          theme_custom +
          theme(plot.title=element_text(hjust=0.5))
figs6.b13