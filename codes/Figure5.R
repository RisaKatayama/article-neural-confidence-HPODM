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

library(coin)
library(pracma)
library(dplyr)
library(ARTool)
library(extrafont)
library(runjags)
library(loo)
library(corrplot)
font_import()
fonts()
library(ggplot2)
library(ggsignif)
library(ggpubr)
library(ggh4x)
library(latex2exp)


.mod.env = new.env()
sys.source("modules.R", envir=.mod.env)
attach(.mod.env)

## data loading ####
dat.folder = "/Users/katayama-r/Documents/BlackJackgame/scripts/data"
mdlres.folder = "/Users/katayama-r/Documents/BlackJackgame/scripts/res_model"
ppi.folder = "/Users/katayama-r/Documents/BlackJackgame/scripts/PPI"
load(paste0(c(dat.folder,"bhvdata_mri.RData"),collapse="/"))

load(paste0(c(mdlres.folder,"cwHDM_dimdl.RData"),collapse="/"))

exsbj = c(7,22)
ns.org = max(Mridata$sidx)
insbj = c(1:ns.org)
insbj = insbj[!is.element(insbj,exsbj)]

mdata = Mridata[!is.element(Mridata$sidx,exsbj),]
mdata$sidx.org = mdata$sidx
for (i in 1:length(insbj)) {
    mdata$sidx[mdata$sidx.org==insbj[i]] = i
}

load.data.basic = function(data){
    ns = max(data$sidx)
    nses = max(data$ses)
    nt = 16
    N = dim(data)[1]

    dkest = (data$deckest==6)*1  # deck 6 = 1
    dktrue = (data$decktrue==6)*1  # deck 6 = 1
    kcf = data$cf.deck
    kcf.hl = rep(NaN,N)
    for (i in 1:ns) {
        tmpkcf=kcf[data$sidx==i]
        tmpkcfhl = (tmpkcf>nanmean(tmpkcf))*1
        kcf.hl[data$sidx==i] = tmpkcfhl
    }

    ev4.rel = c()
    ev6.rel = c()
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
    evd.rea = ev6.rel - ev4.rel

    vr = c()
    for (i in 1:ns) {
        for (s in 1:nses) {
        tmpdkest = dkest[data$sidx==i&data$ses==s]
        tmpdkest0 = (tmpdkest==0)*1
        tmpdkest1 = (tmpdkest==1)*1
        vr = c(vr, (sapply(1:16,function(j){nansum(tmpdkest0[1:j])})*tmpdkest0+sapply(1:16,function(j){nansum(tmpdkest1[1:j])})*tmpdkest1)/c(1:16))
        }
    }

    add = data$add # actual score of the additional card
    
    deci = data$decision # hit = 1, stay = 0
    open = data$open # face-up card score
    icf = data$cf.deci
    icf.hl = rep(NaN,N)
    for (i in 1:ns) {
        tmpicf=icf[data$sidx==i]
        tmpicfhl = (tmpicf>nanmean(tmpicf))*1
        icf.hl[data$sidx==i] = tmpicfhl
    }
    udc = (data$p.us==0.5)*1

    trl = data$trl
    ses = data$ses
    sidx = data$sidx

    seqid = rep(NaN, N)  # sequence id (predetermined displaying patterns of face-up/-down/additional scorecards)
    seqid[data$sidx==1] = data$ses[data$sidx==1]
    for (i in 2:ns) {
        tmppat = array(NaN, dim=c(nses,nt))
        tmp = t(array(data$add[data$sidx==i], dim=c(nt,nses)))
        for (s in 1:nses) {
            tmporg = data$add[data$sidx==1&data$ses==s]
            xi = which(apply(sweep(tmp,2,tmporg,"-")^2,1,sum)==0)
            tmppat[xi,] = s
        }
        seqid[data$sidx==i] = array(t(tmppat))
    }
    
    include = !is.na(data$deckest*data$cf.deck*data$decision*data$cf.deci)
    res = data.frame(dkest=dkest[include], dktrue=dktrue[include], kcf=kcf[include], kcf.hl=kcf.hl[include], add=add[include], 
                     evd.rea=evd.rea[include], vr=vr[include],
                     deci=deci[include], open=open[include], icf=icf[include], icf.hl=icf.hl[include], udc=udc[include],
                     seqid=seqid[include], trl=trl[include], ses=ses[include], sidx=sidx[include])
    return(res)
}
data = load.data.basic(mdata)
incidx = c(1:20,23:24) # sbj 21/22 (original id=23,24) are excluded due to insufficient confidence rating
incidx.org = c(1:6,8:21,25:26)
ns = length(incidx)

colors = c("#009900", "#ff6600")

tmpsummary = add.summary(res$dk.model)
tmpdata = data.frame(dkest=data$dkest, de=data$evd.rea, vr=data$vr, kcf=data$kcf, sidx=data$sidx)
tmpres = model.validate(NULL, tmpdata, res$dk.model$mcmc, rownames(tmpsummary$summaries), "rea.dkmdl")
res$dk.summaries = tmpsummary
res$dk.res = tmpres

tmpsummary = add.summary(res$di.model)
ud = 3.4*(1-data$udc) + 2.5*data$udc
cdk = 4*(1-data$dkest) + 6*data$dkest
udk = 6*(1-data$dkest) + 4*data$dkest
tmpdata = dump.format(list(open=data$open, ud=ud, cdk=cdk, udk=udk, kcf=(res$dk.res$fitres$w.pred.kcf-1)/3, deci=data$deci, icf=data$icf, sidx=data$sidx, N=N, ns=ns))
tmpres = model.validate(NULL, tmpdata, res$di.morel$mcmc, rownames(tmpsummary$summaries), "cwhdm.dimdl")
res$di.summaries = tmpsummary
res$di.res = tmpres

#####################
## Fig 5c (seed = vAIC)
tgroi = c("dACC","l_dAIC","pgACC","r_aPFC","dmPFC")

fc.kcfH = array(NaN, dim=c(ns,length(tgroi)))
fc.kcfL = array(NaN, dim=c(ns,length(tgroi)))

for (i in 1:ns) {
      sdts = readMat(paste0(ppi.folder,"/s",incidx.org[i],"/l_vAIC_8mm.mat"))
      sdts.H = c()
      sdts.L = c()

      tmpkcf = res$dk.res$valres$w.pred.kcf[data$sidx==incidx[i]]
      tmpses = data$ses[data$sidx==incidx[i]]

      fruns = c() # session including only high- (or low-) deck conf trial are excluded
      for (s in 1:8) {
            if (sum(tmpkcf[tmpses==s]>mean(tmpkcf))!=0 & sum(tmpkcf[tmpses==s]<=mean(tmpkcf))!=0){ fruns = c(fruns,s) }
      }

      for (s in fruns) {
            sdts.H = c(sdts.H, sdts$kcfH[[s]][[1]])
            sdts.L = c(sdts.L, sdts$kcfL[[s]][[1]])
      }

      for (r in 1:length(tgroi)) {
            tgts = readMat(paste0(ppi.folder,"/s",incidx.org[i],"/",tgroi[r],"_8mm.mat"))
            tgts.H = c()
            tgts.L = c()
            for (s in fruns) {
                  tgts.H = c(tgts.H, tgts$kcfH[[s]][[1]])
                  tgts.L = c(tgts.L, tgts$kcfL[[s]][[1]])
            }

            fH = lm(tgts.H~sdts.H, data=data.frame(t(rbind(sdts.H,tgts.H))))
            fL = lm(tgts.L~sdts.L, data=data.frame(t(rbind(sdts.L,tgts.L))))
            fc.kcfH[i,r] = coef(fH)[2]
            fc.kcfL[i,r] = coef(fL)[2]
      }
}

dat.fig5c1 = data.frame(fc=c(array(t(fc.kcfH)),array(t(fc.kcfL))), cf=factor(c(rep(1,ns*length(tgroi)),rep(0,ns*length(tgroi)))),
                        tg=factor(rep(1:length(tgroi),ns*2)))

fig5c1 = ggplot(data=dat.fig5c1, aes(x=tg, y=fc, group=interaction(tg,cf), fill=cf)) +
        stat_boxplot(geom='errorbar', linewidth=0.6, width=0.3, position=position_dodge(width=0.75)) +
        geom_boxplot(outlier.shape=3, outlier.size=1.2) +
        geom_hline(yintercept=0, linetype="dotted") +
        geom_signif(xmin=c(1,3,4,5)-0.2, xmax=c(1,3,4,5)+0.2, y_position=c(1.9,1.2,1.2,1.9), 
                    annotations=c('***','**','*','***'), tip_length=0, vjust=0.2, textsize=4.5) +
        scale_fill_manual(values=colors, labels=c("Low","High")) +
        scale_y_continuous(limits=c(-1,2),breaks=seq(-1,2,0.5),name="Coupling strength\n(vAIC ~ target)") +
        scale_x_discrete(labels=c("dACC","dAIC","pACC","aPFC","dmPFC"), name="Target regions") +
        theme_custom +
        theme(legend.position=c(0.05,0.05), legend.justification=c(0,0), legend.background=element_rect(colour="transparent"),
              legend.title=element_text(size=10),
              legend.text=element_text(size=10), legend.spacing.x=unit(0.1,'cm'), legend.spacing.y=unit(0.1,'cm'),
              legend.key.width = unit(0.4,"cm"), axis.text.x=element_text(angle=30, margin=margin(t=0.05,b=0.15,unit="cm"), hjust=1),
              axis.title.x=element_text(margin=margin(t=0.05,unit="cm"))) +
        guides(fill=guide_legend(title="Deck conf.", title.hjust=0.5, ncol=2))
fig5c1


## Fig 5c (seed = dACC)
tgroi = c("l_vAIC","l_dAIC","pgACC","r_aPFC","dmPFC")

fc.kcfH = array(NaN, dim=c(ns,length(tgroi)))
fc.kcfL = array(NaN, dim=c(ns,length(tgroi)))

for (i in 1:ns) {
      sdts = readMat(paste0(ppi.folder,"/s",incidx.org[i],"/dACC_8mm.mat"))
      sdts.H = c()
      sdts.L = c()

      tmpkcf = res$dk.res$valres$w.pred.kcf[data$sidx==incidx[i]]
      tmpses = data$ses[data$sidx==incidx[i]]

      fruns = c() # session including only high- (or low-) deck conf trial are excluded
      for (s in 1:8) {
            if (sum(tmpkcf[tmpses==s]>mean(tmpkcf))!=0 & sum(tmpkcf[tmpses==s]<=mean(tmpkcf))!=0){ fruns = c(fruns,s) }
      }

      for (s in fruns) {
            sdts.H = c(sdts.H, sdts$kcfH[[s]][[1]])
            sdts.L = c(sdts.L, sdts$kcfL[[s]][[1]])
      }

      for (r in 1:length(tgroi)) {
            tgts = readMat(paste0(ppi.folder,"/s",incidx.org[i],"/",tgroi[r],"_8mm.mat"))
            tgts.H = c()
            tgts.L = c()
            for (s in fruns) {
                  tgts.H = c(tgts.H, tgts$kcfH[[s]][[1]])
                  tgts.L = c(tgts.L, tgts$kcfL[[s]][[1]])
            }

            fH = lm(tgts.H~sdts.H, data=data.frame(t(rbind(sdts.H,tgts.H))))
            fL = lm(tgts.L~sdts.L, data=data.frame(t(rbind(sdts.L,tgts.L))))
            fc.kcfH[i,r] = coef(fH)[2]
            fc.kcfL[i,r] = coef(fL)[2]
      }
}

dat.fig5c2 = data.frame(fc=c(array(t(fc.kcfH)),array(t(fc.kcfL))), cf=factor(c(rep(1,ns*length(tgroi)),rep(0,ns*length(tgroi)))),
                        tg=factor(rep(1:length(tgroi),ns*2)))

fig5c2 = ggplot(data=dat.fig5c2, aes(x=tg, y=fc, group=interaction(tg,cf), fill=cf)) +
        stat_boxplot(geom='errorbar', linewidth=0.6, width=0.4, position=position_dodge(width=0.75)) +
        geom_boxplot(outlier.shape=3, outlier.size=1.2) +
        geom_hline(yintercept=0, linetype="dotted") +
        geom_signif(xmin=c(2,5)-0.2, xmax=c(2,5)+0.2, y_position=c(1.2,2), 
                    annotations=c('*','**'), tip_length=0, vjust=0.2, textsize=4.5) +
        scale_fill_manual(values=colors, labels=c("Low","High")) +
        scale_y_continuous(limits=c(-1,2.5),breaks=seq(-1,2.5,0.5),name="Coupling strength\n(dACC ~ target)") +
        scale_x_discrete(labels=c("vAIC","dAIC","pACC","aPFC","dmPFC"), name="Target regions") +
        theme_custom +
        theme(legend.position=c(0.05,0.05), legend.justification=c(0,0), legend.background=element_rect(colour="transparent"),
              legend.title=element_text(size=10),
              legend.text=element_text(size=10), legend.spacing.x=unit(0.1,'cm'), legend.spacing.y=unit(0.1,'cm'),
              legend.key.width = unit(0.4,"cm"), axis.text.x=element_text(angle=30, margin=margin(t=0.05,b=0.15,unit="cm"), hjust=1),
              axis.title.x=element_text(margin=margin(t=0.05,unit="cm"))) +
        guides(fill=guide_legend(title="Deck conf.", title.hjust=0.5, ncol=2))
fig5c2


## Fig 5d
g.kcf = sapply(incidx, function(i){res$di.summaries$summaries[paste0("g.kcf.p[",i,"]"),"Median"]})
g.kcf.hl = factor((g.kcf>mean(g.kcf))*1)

dfc.kcfHL.vAIC.dACC = rep(NaN, ns)
for (i in 1:ns) {
      sdts = readMat(paste0(ppi.folder,"/s",incidx.org[i],"/l_vAIC_8mm.mat"))
      sdts.H = c()
      sdts.L = c()

      tmpkcf = res$dk.res$valres$w.pred.kcf[data$sidx==incidx[i]]
      tmpses = data$ses[data$sidx==incidx[i]]

      fruns = c() # session including only high- (or low-) deck conf trial are excluded
      for (s in 1:8) {
            if (sum(tmpkcf[tmpses==s]>mean(tmpkcf))!=0 & sum(tmpkcf[tmpses==s]<=mean(tmpkcf))!=0){ fruns = c(fruns,s) }
      }

      for (s in fruns) {
            sdts.H = c(sdts.H, sdts$kcfH[[s]][[1]])
            sdts.L = c(sdts.L, sdts$kcfL[[s]][[1]])
      }

      tgts = readMat(paste0(ppi.folder,"/s",incidx.org[i],"/dACC_8mm.mat"))
      tgts.H = c()
      tgts.L = c()
      for (s in fruns) {
        tgts.H = c(tgts.H, tgts$kcfH[[s]][[1]])
        tgts.L = c(tgts.L, tgts$kcfL[[s]][[1]])
      }
      
      fH = lm(tgts.H~sdts.H, data=data.frame(t(rbind(sdts.H,tgts.H))))
      fL = lm(tgts.L~sdts.L, data=data.frame(t(rbind(sdts.L,tgts.L))))
      dfc.kcfHL.vAIC.dACC[i] = coef(fH)[2] - coef(fL)[2]      
}

wilcox_test(dfc.kcfHL.vAIC.dACC ~ g.kcf.hl, alternative="greater")

dat.fig5d = data.frame(gkcf=g.kcf.hl, dfc=dfc.kcfHL.vAIC.dACC)
fig5d = ggplot(data=dat.fig5d, aes(x=gkcf, y=dfc, group=gkcf)) +
        stat_boxplot(geom='errorbar', linewidth=0.6, width=0.4, position=position_dodge(width=0.75)) +
        geom_boxplot(outlier.shape=3, outlier.size=1.2, fill="gray") +
        geom_signif(comparisons=list(c(1,2)), annotations=c('*'), y_position=0.2, tip_length=0, vjust=0.2, textsize=4.5) +
        scale_y_continuous(limits=c(-1,0.51),breaks=seq(-1,0.5,0.5), name="Difference in coupling strength\nbtw high and low deck conf.\n(vAIC ~ dACC)") +
        scale_x_discrete(labels=c("Low","High"), name="Deck conf. weighting of\nvalue belief modulation") +
        theme_custom
fig5d


## Fig 5e
e.kcf = sapply(incidx, function(i){res$di.summaries$summaries[paste0("e.kcf.p[",i,"]"),"Median"]})
e.kcf.hl = factor((e.kcf>mean(e.kcf))*1)

fc = rep(NaN, ns)
for (i in 1:ns) {
      ps = readMat(paste0(ppi.folder,"/s",incidx.org[i],"/dACC_8mm.mat"))
      sdts = c()
      for (s in 1:8) {
            sdts = c(sdts,t(ps$dec[[s]][[1]]))
      }
      
      pt = readMat(paste0(ppi.folder,"/s",incidx.org[i],"/MCC_8mm.mat"))
      tgts = c()
      for (s in 1:8) {
        tgts = c(tgts,t(pt$dec[[s]][[1]]))
      }
      f = lm(tgts~sdts, data=data.frame(t(rbind(sdts,tgts))))
      fc[i] = coef(f)[2]
}

wilcox_test(fc ~ e.kcf.hl, alternative="less")

dat.fig5e = data.frame(ekcf=e.kcf.hl, fc=fc)

fig5e = ggplot(data=dat.fig5e, aes(x=ekcf, y=fc, group=ekcf)) +
        stat_boxplot(geom='errorbar', linewidth=0.6, width=0.4, position=position_dodge(width=0.75)) +
        geom_boxplot(outlier.shape=3, outlier.size=1.2, fill="gray") +
        geom_signif(comparisons=list(c(1,2)), annotations=c('*'), y_position=1.7, tip_length=0, vjust=0.2, textsize=4.5) +
        scale_y_continuous(limits=c(0,2), breaks=seq(0,2,0.5), name="Coupling strength (dACC ~ MCC)") +
        scale_x_discrete(labels=c("Low","High"), name="Deck inference\ndependency") +
        theme_custom
fig5e