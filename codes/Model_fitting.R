# Copyright (C) <2024>  <Risa Katayama>
# 
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
# 
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
# 
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <http://www.gnu.org/licenses/>.

library(loo)
library(runjags)

.mod.env = new.env()
sys.source("modules.R", envir=.mod.env)
attach(.mod.env)

## data loading ####
dat.folder = ""
mdl.folder = ""
load(paste0(c(dat.folder,"bhvdata_bhv.RData"),collapse="/"))
load(paste0(c(dat.folder,"bhvdata_mri.RData"),collapse="/"))

exsbj = c(7,22) # subjects with low deck inference accuracy
ns.org = max(Bhvdata$sidx)
insbj = c(1:ns.org)
insbj = insbj[!is.element(insbj,exsbj)]
Bhvdata = Bhvdata[!is.element(Bhvdata$sidx,exsbj),]
Mridata = Mridata[!is.element(Mridata$sidx,exsbj),]
Bhvdata$sidx.org = Bhvdata$sidx
Mridata$sidx.org = Mridata$sidx
ns = length(insbj)
for (i in 1:ns) {
   Bhvdata$sidx[Bhvdata$sidx.org==insbj[i]] = i
   Mridata$sidx[Mridata$sidx.org==insbj[i]] = i
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
      
      x = c(5,tmpadd)
      tmpm.add = mapply(function(t,x){mean(x[1:t])}, 1:nt, MoreArgs=list(x=x))
      m.add = c(m.add,tmpm.add)
    }
  }

  p4.bio = c()
  p6.bio = c()
  for (i in 1:ns) {
    for (s in 1:nses) {
      tmpp4 = rep(0.5,nt)
      tmpadd = data$add[data$sidx==i&data$ses==s]
      for (t in 1:(nt-1)) {
        pri = c(tmpp4[t],1-tmpp4[t])
        lik = dbinom(tmpadd[t],10,c(0.4,0.6))
        pst = pri*lik
        pst = pst/sum(pst)
        tmpp4[t+1] = pst[1]
      }
      p4.bio = c(p4.bio,tmpp4)
      p6.bio = c(p6.bio,1-tmpp4)
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
                   ev4.rel=ev4.rel[include], ev6.rel=ev6.rel[include], m.add=m.add[include], p4.bio=p4.bio[include], p6.bio=p6.bio[include], vr=vr[include],
                   icf=data$cf.deci[include], deci=data$decision[include], open=data$open[include], 
                   vh.cdk=vh.cdk[include], vh.udk=vh.udk[include], vh.add=vh.add[include],
                   sidx=data$sidx[include], sidx.org=data$sidx.org[include])
  
  res$evd.rea = res$ev6.rel - res$ev4.rel
  res$evd.mean = res$m.add - 5
  res$pch.bio = (res$p4.bio-res$p6.bio)*(1-res$dkest) + (res$p6.bio-res$p4.bio)*res$dkest
  res$ent.bio = - (res$p4.bio*log(res$p4.bio) + res$p6.bio*log(res$p6.bio))
  
  return(res)
}

bdata = load.data(Bhvdata)
mdata = load.data(Mridata)
N = dim(bdata)[1]
N.m = dim(mdata)[1]

inits1 = dump.format(list(.RNG.name="base::Super-Duper", .RNG.seed=111))
inits2 = dump.format(list(.RNG.name="base::Wichmann-Hill", .RNG.seed=222))
inits3 = dump.format(list(.RNG.name="base::Mersenne-Twister", .RNG.seed=333))


## cw-HDM (proposed) model ############
## deck inference phase
dat.dk = dump.format(list(dkest=bdata$dkest, de=bdata$evd.rea, vr=bdata$vr, kcf=bdata$kcf, sidx=bdata$sidx, N=N, ns=ns))
vals   = c("b.dk.mu","b.kcf.mu","a0.kcf.mu","a1.kcf.mu","b.dk.sig","b.kcf.sig","a0.kcf.sig",
           "b0.dk.p","b1.dk.p","b1.kcf.p","b2.kcf.p","a0.kcf.p")

rea.dkmdl = run.jags(model=paste0(mdl.folder,"basic_dkmdl.txt"), monitor=vals, data=dat.dk, n.chains=3, inits=c(inits1,inits2,inits3),
                     plots=FALSE, method="parallel", burnin=10000, sample=1000, adapt=1000, thin=20)

summary.rea.dkmdl = add.summary(rea.dkmdl)

tmpbdata = data.frame(dkest=bdata$dkest, de=bdata$evd.rea, vr=bdata$vr, kcf=bdata$kcf, sidx=bdata$sidx)
tmpmdata = data.frame(dkest=mdata$dkest, de=mdata$evd.rea, vr=mdata$vr, kcf=mdata$kcf, sidx=mdata$sidx)

res.rea.dkmdl = model.validate(tmpbdata, tmpmdata, rea.dkmdl$mcmc, rownames(summary.rea.dkmdl$summaries), "rea.dkmdl")

## decision phase
ud = 3.4*(1-bdata$udc) + 2.5*bdata$udc
cdk = 4*(1-bdata$dkest) + 6*bdata$dkest
udk = 6*(1-bdata$dkest) + 4*bdata$dkest
dat.di = dump.format(list(open=bdata$open, ud=ud, cdk=cdk, udk=udk, deci=bdata$deci, icf=bdata$icf, kcf=(res.rea.dkmdl$fitres$w.pred.kcf-1)/3,
                          sidx=bdata$sidx, N=N, ns=ns))
vals = c("b.di.mu","b.icf.mu","a0.icf.mu","a1.icf.mu","b.di.sig","b.icf.sig","a0.icf.sig",
         "g.kcf.p","e.kcf.p","b0.di.p","b1.di.p","b2.di.p","b1.icf.p","b2.icf.p","a0.icf.p")

cwhdm.dimdl = run.jags(model=paste0(mdl.folder,"cwHDM_dimdl.txt"), monitor=vals, data=dat.di, n.chains=3, inits=c(inits1,inits2,inits3),
                    plots=FALSE, method="parallel", burnin=10000, sample=1000, adapt=1000, thin=20)

summary.cwhdm.dimdl = add.summary(cwhdm.dimdl)

ud = 3.4*(1-bdata$udc) + 2.5*bdata$udc
cdk = 4*(1-bdata$dkest) + 6*bdata$dkest
udk = 4*bdata$dkest + 6*(1-bdata$dkest)
tmpbdata = data.frame(open=bdata$open, ud=ud, cdk=cdk, udk=udk, kcf=(res.rea.dkmdl$fitres$w.pred.kcf-1)/3, deci=bdata$deci, icf=bdata$icf, sidx=bdata$sidx)
ud = 3.4*(1-mdata$udc) + 2.5*mdata$udc
cdk = 4*(1-mdata$dkest) + 6*mdata$dkest
udk = 4*mdata$dkest + 6*(1-mdata$dkest)
tmpmdata = data.frame(open=mdata$open, ud=ud, cdk=cdk, udk=udk, kcf=(res.rea.dkmdl$valres$w.pred.kcf-1)/3, deci=mdata$deci, icf=mdata$icf, sidx=mdata$sidx)

res.cwhdm.dimdl = model.validate(tmpbdata, tmpmdata, cwhdm.dimdl$mcmc, rownames(summary.cwhdm.dimdl$summaries), "cwhdm.dimdl")


## fp-HDM  model ############
## deck inference phase
# same as cw-HDM model

## decision phase
ud = 3.4*(1-bdata$udc) + 2.5*bdata$udc
modadd = apply(sweep(res.rea.dkmdl$fitres$pmat.dk, 2, c(4,6), "*"), 1, sum)
dat.di = dump.format(list(open=bdata$open, ud=ud, add=modadd, deci=bdata$deci, icf=bdata$icf, kcf=(res.rea.dkmdl$fitres$w.pred.kcf-1)/3,
                          sidx=bdata$sidx, N=N, ns=ns))
vals = c("b.di.mu","b.icf.mu","a0.icf.mu","a1.icf.mu","b.di.sig","b.icf.sig","a0.icf.sig",
         "b0.di.p","b1.di.p","b2.di.p","b1.icf.p","b2.icf.p","a0.icf.p")

fphdm.dimdl = run.jags(model=paste0(mdl.folder,"fpHDM_dimdl.txt"), monitor=vals, data=dat.di, n.chains=3, inits=c(inits1,inits2,inits3),
                       plots=FALSE, method="parallel", burnin=10000, sample=1000, adapt=1000, thin=20)

summary.fphdm.dimdl = add.summary(fphdm.dimdl)

ud = 3.4*(1-bdata$udc) + 2.5*bdata$udc
modadd = apply(sweep(res.rea.dkmdl$fitres$pmat.dk, 2, c(4,6), "*"), 1, sum)
tmpbdata = data.frame(open=bdata$open, ud=ud, add=modadd, kcf=(res.rea.dkmdl$fitres$w.pred.kcf-1)/3, deci=bdata$deci, icf=bdata$icf, sidx=bdata$sidx)
ud = 3.4*(1-mdata$udc) + 2.5*mdata$udc
modadd = apply(sweep(res.rea.dkmdl$valres$pmat.dk, 2, c(4,6), "*"), 1, sum)
tmpmdata = data.frame(open=mdata$open, ud=ud, add=modadd, kcf=(res.rea.dkmdl$valres$w.pred.kcf-1)/3, deci=mdata$deci, icf=mdata$icf, sidx=mdata$sidx)

res.fphdm.dimdl = model.validate(tmpbdata, tmpmdata, fphdm.dimdl$mcmc, rownames(summary.fphdm.dimdl$summaries), "fphdm.dimdl")


## HDM  model ############
## deck inference phase
# same as cw-HDM model

## decision phase
ud = 3.4*(1-bdata$udc) + 2.5*bdata$udc
cdk = 4*(1-bdata$dkest) + 6*bdata$dkest
dat.di = dump.format(list(open=bdata$open, ud=ud, add=cdk, deci=bdata$deci, icf=bdata$icf, sidx=bdata$sidx, N=N, ns=ns))
vals = c("b.di.mu","b.icf.mu","a0.icf.mu","a1.icf.mu","b.di.sig","b.icf.sig","a0.icf.sig",
         "b0.di.p","b1.di.p","b2.di.p","b1.icf.p","a0.icf.p")

hdm.dimdl = run.jags(model=paste0(mdl.folder,"basic_dimdl.txt"), monitor=vals, data=dat3, n.chains=3, inits=c(inits1,inits2,inits3),
                     plots=FALSE, method="parallel", burnin=10000, sample=1000, adapt=1000, thin=20)

summary.hdm.dimdl = add.summary(hdm.dimdl)

ud = 3.4*(1-bdata$udc) + 2.5*bdata$udc
cdk = 4*(1-bdata$dkest) + 6*bdata$dkest
tmpbdata = data.frame(open=bdata$open, ud=ud, add=cdk, deci=bdata$deci, icf=bdata$icf, sidx=bdata$sidx)
ud = 3.4*(1-mdata$udc) + 2.5*mdata$udc
cdk = 4*(1-mdata$dkest) + 6*mdata$dkest
tmpmdata = data.frame(open=mdata$open, ud=ud, add=cdk, deci=mdata$deci, icf=mdata$icf, sidx=mdata$sidx)

res.hdm.dimdl = model.validate(tmpbdata, tmpmdata, hdm.dimdl$mcmc, rownames(summary.hdm.dimdl$summaries), "hdm.dimdl")

ud = 3.4*(1-bdata$udc) + 2.5*bdata$udc
cdk = 4*(1-bdata$dkest) + 6*bdata$dkest
tmpbdata = data.frame(open=bdata$open, ud=ud, add=cdk, deci=bdata$deci, icf=bdata$icf, sidx=bdata$sidx)
ud = 3.4*(1-mdata$udc) + 2.5*mdata$udc
cdk = 4*(1-mdata$dkest) + 6*mdata$dkest
tmpmdata = data.frame(open=mdata$open, ud=ud, add=cdk, deci=mdata$deci, icf=mdata$icf, sidx=mdata$sidx)

res.hdm.dimdl = model.validate(tmpbdata, tmpmdata, hdm.dimdl$mcmc, rownames(summary.hdm.dimdl$summaries), "hdm.dimdl")


## PDM  model ############
## deck inference phase
dat.dk = dump.format(list(dkest=bdata$dkest, de=bdata$evd.mean, vr=bdata$vr, kcf=bdata$kcf, sidx=bdata$sidx, N=N, ns=ns))
vals   = c("b.dk.mu","b.kcf.mu","a0.kcf.mu","a1.kcf.mu","b.dk.sig","b.kcf.sig","a0.kcf.sig",
           "b0.dk.p","b1.dk.p","b1.kcf.p","b2.kcf.p","a0.kcf.p")

avg.dkmdl = run.jags(model=paste0(mdl.folder,"basic_dkmdl.txt"), monitor=vals, data=dat.dk, n.chains=3, inits=c(inits1,inits2,inits3),
                     plots=FALSE, method="parallel", burnin=10000, sample=1000, adapt=1000, thin=20)

summary.avg.dkmdl = add.summary(avg.dkmdl)

tmpbdata = data.frame(dkest=bdata$dkest, de=bdata$evd.mean, vr=bdata$vr, kcf=bdata$kcf, sidx=bdata$sidx)
tmpmdata = data.frame(dkest=mdata$dkest, de=mdata$evd.mean, vr=mdata$vr, kcf=mdata$kcf, sidx=mdata$sidx)

res.avg.dkmdl = model.validate(tmpbdata, tmpmdata, avg.dkmdl$mcmc, rownames(summary.avg.dkmdl$summaries), "avg.dkmdl")

## decision phase
ud = 3.4*(1-bdata$udc) + 2.5*bdata$udc
dat.di = dump.format(list(open=bdata$open, ud=ud, add=bdata$m.add, deci=bdata$deci, icf=bdata$icf, sidx=bdata$sidx, N=N, ns=ns))
vals = c("b.di.mu","b.icf.mu","a0.icf.mu","a1.icf.mu","b.di.sig","b.icf.sig","a0.icf.sig",
         "b0.di.p","b1.di.p","b2.di.p","b1.icf.p","a0.icf.p")

pdm.dimdl = run.jags(model=paste0(mdl.folder,"basic_dimdl.txt"), monitor=vals, data=dat3, n.chains=3, inits=c(inits1,inits2,inits3),
                     plots=FALSE, method="parallel", burnin=10000, sample=1000, adapt=1000, thin=20)

summary.pdm.dimdl = add.summary(pdm.dimdl)

ud = 3.4*(1-bdata$udc) + 2.5*bdata$udc
tmpbdata = data.frame(open=bdata$open, ud=ud, add=bdata$m.add, deci=bdata$deci, icf=bdata$icf, sidx=bdata$sidx)
ud = 3.4*(1-mdata$udc) + 2.5*mdata$udc
tmpmdata = data.frame(open=mdata$open, ud=ud, add=mdata$m.add, deci=mdata$deci, icf=mdata$icf, sidx=mdata$sidx)

res.pdm.dimdl = model.validate(tmpbdata, tmpmdata, pdm.dimdl$mcmc, rownames(summary.pdm.dimdl$summaries), "pdm.dimdl")


## deck inference alternative  model ############
## Bayesian ideal oberver
dat.dk = dump.format(list(dkest=bdata$dkest, ent=bdata$ent.bio, vr=bdata$vr, kcf=bdata$kcf, sidx=bdata$sidx, N=N, ns=ns))
vals = c("b.kcf.mu","a0.kcf.mu","a1.kcf.mu","b.kcf.sig","a0.kcf.sig","b1.kcf.p","b2.kcf.p","a0.kcf.p")

bio.dkmdl = run.jags(model=model=paste0(mdl.folder,"bio_dkmdl.txt"), monitor=vals, data=dat.dk, n.chains=3, inits=c(inits1,inits2,inits3),
                     plots=FALSE, method="parallel", burnin=10000, sample=1000, adapt=1000, thin=20)

summary.bio.dkmdl = add.summary(bio.dkmdl)

tmpbdata = data.frame(dkest=bdata$dkest, pd6=bdata$p6.bio, ent=bdata$ent.bio, vr=bdata$vr, kcf=bdata$kcf, sidx=bdata$sidx)
tmpmdata = data.frame(dkest=mdata$dkest, pd6=mdata$p6.bio, ent=mdata$ent.bio, vr=mdata$vr, kcf=mdata$kcf, sidx=mdata$sidx)

res.bio.dkmdl = model.validate(tmpbdata, tmpmdata, bio.dkmdl$mcmc, rownames(summary.bio.dkmdl$summaries), "bio.dkmdl")