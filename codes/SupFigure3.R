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
dat.folder = ""
mdl.folder = "res_models"
load(paste0(c(dat.folder,"bhvdata_bhv.RData"),collapse="/"))

load(paste0(c(dat.folder,mdl.folder,"BIO_dkmdl.RData"),collapse="/"))
res.bio = res 
load(paste0(c(dat.folder,mdl.folder,"AVE_dkmdl.RData"),collapse="/"))
res.ave = res 
load(paste0(c(dat.folder,mdl.folder,"REA_dkmdl.RData"),collapse="/"))
res.rea = res

inits1 = dump.format(list(.RNG.name="base::Super-Duper", .RNG.seed=111))
inits2 = dump.format(list(.RNG.name="base::Wichmann-Hill", .RNG.seed=222))
inits3 = dump.format(list(.RNG.name="base::Mersenne-Twister", .RNG.seed=333))

colors = c("#7ad151","#bebebe","#fde725")

exsbj = c(7,22) # subjects with low deck inference accuracy
ns.org = max(Bhvdata$sidx)
insbj = c(1:ns.org)
insbj = insbj[!is.element(insbj,exsbj)]
bdata = Bhvdata[!is.element(Bhvdata$sidx,exsbj),]
bdata$sidx.org = bdata$sidx
ns = length(insbj)
for (i in 1:ns) {
   bdata$sidx[bdata$sidx.org==insbj[i]] = i
}

load.data.basic = function(data){
    ns = max(data$sidx)
    nses = max(data$ses)
    nt = 16
    N = dim(data)[1]

    dkest = (data$deckest==6)*1  # deck 6 = 1
    dktrue = (data$decktrue==6)*1  # deck 6 = 1
    kcf = data$cf.deck
    add = data$add

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

    vr = c() # deck inference stability
    for (i in 1:ns) {
        for (s in 1:nses) {
            tmpdkest = dkest[data$sidx==i&data$ses==s]
            tmpdkest0 = (tmpdkest==0)*1
            tmpdkest1 = (tmpdkest==1)*1
            vr = c(vr, (sapply(1:16,function(j){nansum(tmpdkest0[1:j])})*tmpdkest0+sapply(1:16,function(j){nansum(tmpdkest1[1:j])})*tmpdkest1)/c(1:16))
        }
    }

    trl = data$trl
    ses = data$ses
    sidx = data$sidx
    
    include = !is.na(data$deckest*data$cf.deck*data$decision*data$cf.deci)
    res = data.frame(dkest=dkest[include], dktrue=dktrue[include], kcf=kcf[include], add=add[include],
               evd.rea=evd.rea[include], vr=vr[include], trl=trl[include], ses=ses[include], sidx=sidx[include])
    return(res)
}
data = load.data.basic(bdata)
ns = max(data$sidx)
ng = max(data$ses)
N = dim(data)[1]

calc.pred2 = function(p){
    m = length(p)
    if (is.element(p[m],which(p[1:(m-1)]==max(p[1:(m-1)])))){return(p[m])}
    else {return(which.max(p[1:(m-1)])-1)}
}

calc.pred4 = function(p){
    m = length(p)
    if (is.element(p[m],which(p[1:(m-1)]==max(p[1:(m-1)])))){return(p[m])}
    else {return(which.max(p[1:(m-1)]))}
}

calc.nmap = function(p){ return(nansum(p==max(p))) }

#####################
## Fig S3a
dkacc.trl = array(NaN,dim=c(ns,16))
for (i in 1:ns) {
    for (t in 1:16) {
        dkacc.trl[i,t] = sum(data$dkest==data$dktrue&data$trl==t&data$sidx==i)/sum(data$trl==t&data$sidx==i)
    }
}

corrdk.rea = (data$dktrue==res.rea$res$fitres$pred.dk)*1
corrdk.ave = (data$dktrue==res.ave$res$fitres$pred.dk)*1
corrdk.bio = (data$dktrue==res.bio$res$fitres$pred.dk)*1

simdkacc.trl.rea = array(NaN,dim=c(ns,16))
simdkacc.trl.ave = array(NaN,dim=c(ns,16))
simdkacc.trl.bio = array(NaN,dim=c(ns,16))
for (i in 1:ns) {
    for (t in 1:16) {
        simdkacc.trl.rea[i,t] = mean(corrdk.rea[data$sidx==i&data$trl==t])
        simdkacc.trl.ave[i,t] = mean(corrdk.ave[data$sidx==i&data$trl==t])
        simdkacc.trl.bio[i,t] = mean(corrdk.bio[data$sidx==i&data$trl==t])
    }
}

dat.figs3a = data.frame(acc=array(t(dkacc.trl)), trl=factor(rep(1:16,ns)))
dat.figs3a = dat.figs3a[!is.na(dat.figs3a$acc),]

simdat.figs3a = data.frame(acc=c(apply(simdkacc.trl.rea,2,nanmean),apply(simdkacc.trl.bio,2,nanmean),apply(simdkacc.trl.ave,2,nanmean)),
                           s=c(apply(simdkacc.trl.rea,2,sd)/sqrt(apply(!is.na(simdkacc.trl.rea),2,sum)),apply(simdkacc.trl.bio,2,nansd)/sqrt(apply(!is.na(simdkacc.trl.bio),2,sum)),
                               apply(simdkacc.trl.ave,2,nansd)/sqrt(apply(!is.na(simdkacc.trl.ave),2,sum))),
                           trl=rep(1:16,3), mdl=factor(c(rep('REA',16),rep('BIO',16),rep('AVG',16)),levels=c('REA','BIO','AVG')))
simdat.figs3a = simdat.figs3a[!is.na(simdat.figs3a$acc),]

figs3a = ggplot(data=dat.figs3a, aes(x=trl, y=acc)) +
        stat_boxplot(geom='errorbar', linewidth=0.6, width=0.3) +
        geom_boxplot(outlier.shape=3, outlier.size=1.2, fill="white") +
        geom_ribbon(data=simdat.figs3a, aes(x=trl, ymin=acc-s, ymax=acc+s, group=mdl, fill=mdl), colour=NA, alpha=0.5, show.legend=FALSE) +
        geom_line(data=simdat.figs3a, aes(x=trl, y=acc, group=mdl, color=mdl), linewidth=0.6, alpha=1) +
        geom_hline(yintercept=0.5, linetype="dotted") +
        scale_fill_manual(name=NULL, values=colors) +
        scale_color_manual(name=NULL, values=colors, labels=c(' REA',' BIO',' AVG')) +
        guides(color=guide_legend(override.aes=list(linewidth=1))) +
        scale_x_discrete(name="Trial index in a single game") +
        scale_y_continuous(limits=c(0,1), name="Deck inference accuracy") + 
        theme_custom +
        theme(legend.position=c(1,0.02),legend.justification=c(1,0),legend.background=element_rect(fill = "transparent"),
              legend.title=element_text(size=8),
              legend.text=element_text(size=8, margin=margin(r=0,l=0.1,unit='cm')),legend.spacing.x=unit(-0.05,'cm'),legend.spacing.y=unit(0.2,'cm'),
              legend.key.width=unit(0.5,"cm"), legend.key.height=unit(0.4,"cm"))
figs3a


## Fig S3b
qt.evd = quantile(abs(data$evd.rea),seq(0,1,0.2))
kcf.evd = array(NaN, dim=c(ns,length(qt.evd)-1))
for (i in 1:ns) {
    tmpkcf = data$kcf[data$sidx==i]
    tmpevd = data$evd.rea[data$sidx==i]
    kcf.evd[i,] = sapply(1:(length(qt.evd)-1), function(x){mean(tmpkcf[tmpevd>qt.evd[x]&tmpevd<=qt.evd[x+1]])})
}

dat.figs3b = data.frame(kcf=array(t(kcf.evd)), evd=factor(rep(1:(length(qt.evd)-1),ns)))
dat.figs3b = dat.figs3b[!is.na(dat.figs3b$kcf),]

figs3b = ggplot(data=dat.figs3b, aes(x=evd, y=kcf)) +
        stat_boxplot(geom='errorbar', linewidth=0.6, width=0.3) +
        geom_boxplot(outlier.shape=3, outlier.size=1.2, fill="gray") +
        scale_x_discrete(name=expression(paste("|",{E[d6]-E[d4]}, "| quantile", sep=""))) +
        scale_y_continuous(limits=c(0.99,4.01), breaks=c(1:4), name="Averaged deck confidence level", labels=scales::number_format(accuracy=0.1)) + 
        theme_custom
figs3b


## Fig S3c
qt.evd = quantile(abs(data$evd.rea),seq(0,1,0.2))
kcf.evd.s = array(NaN,dim=c(ns,length(qt.evd)-1,2))
for (i in 1:ns) {
    tmpkcf = data$kcf[data$sidx==i]
    tmpevd = data$evd.rea[data$sidx==i]
    tmps = data$vr[data$sidx==i]
    for (x in 1:length(qt.evd)-1) {
        kcf.evd.s[i,x,1] = nanmean(tmpkcf[tmps<=0.5&tmpevd>qt.evd[x]&tmpevd<=qt.evd[x+1]])
        kcf.evd.s[i,x,2] = nanmean(tmpkcf[tmps>0.5&tmpevd>qt.evd[x]&tmpevd<=qt.evd[x+1]])
    }
}

dat.figs3c = data.frame(kcf=array(aperm(kcf.evd.s,c(2,3,1))), evd=factor(rep(1:(length(qt.evd)-1),2*ns)), s=factor(rep(c(rep(0,length(qt.evd)-1),rep(1,length(qt.evd)-1)),ns)))
dat.figs3c = dat.figs3c[!is.na(dat.figs3c$kcf),]

figs3c = ggplot(data=dat.figs3c, aes(x=evd, y=kcf, group=interaction(evd,s), fill=s)) +
        stat_boxplot(geom='errorbar', linewidth=0.6, width=0.4, position=position_dodge(width=0.75)) +
        geom_boxplot(outlier.shape=3, outlier.size=1.2) +
        scale_fill_manual(values=c("lightgray","dimgray"), labels=c("Low","High")) +
        scale_y_continuous(limits=c(0.99,4.2), breaks=c(1:4), name="Averaged deck confidence level", labels=scales::number_format(accuracy=0.1)) + 
        scale_x_discrete(name=expression(paste("|",{E[d6]-E[d4]}, "| quantile", sep=""))) +
        geom_signif(xmin=c(1.8,2.8,3.8,4.8), xmax=c(2.2,3.2,4.2,5.2), annotations="***", y_position=4.12, tip_length=0, vjust=0.2, textsize=4.5) +
        theme_custom +
        theme(legend.position="top",legend.background=element_rect(colour="transparent"),
              legend.title=element_text(size=10),legend.title.align=0.5,
              legend.text=element_text(size=10), legend.spacing.x=unit(0.1,'cm'), legend.spacing.y=unit(0.1,'cm'),
              legend.key.width = unit(0.4,"cm"), axis.text.x=element_text(),
              axis.title.x=element_text(margin=margin(t=0.05,unit="cm"))) +
        guides(fill=guide_legend(title="Deck inference\nstability", ncol=2, title.position="left"))
figs3c


## Fig S3e
kcf.trl = array(NaN, dim=c(ns,16))
for (i in 1:ns) {
    for (t in 1:16) {
        kcf.trl[i,t] = mean(data$kcf[data$sidx==i&data$trl==t])
    }
}

simkcf.trl.rea = array(NaN, dim=c(ns,16))
simkcf.trl.ave = array(NaN, dim=c(ns,16))
simkcf.trl.bio = array(NaN, dim=c(ns,16))
for (i in 1:ns) {
    for (t in 1:16) {
       simkcf.trl.rea[i,t] = mean(res.rea$res$fitres$w.pred.kcf[data$sidx==i&data$trl==t])
       simkcf.trl.ave[i,t] = mean(res.ave$res$fitres$w.pred.kcf[data$sidx==i&data$trl==t])
       simkcf.trl.bio[i,t] = mean(res.bio$res$fitres$w.pred.kcf[data$sidx==i&data$trl==t])
    }
}

dat.figs3e = data.frame(cf=array(t(kcf.trl)), trl=factor(rep(1:16,ns)))
dat.figs3e = dat.figs3e[!is.na(dat.figs3e$cf),]

simdat.figs3e = data.frame(cf=c(apply(simkcf.trl.rea,2,nanmean),apply(simkcf.trl.bio,2,nanmean),apply(simkcf.trl.ave,2,nanmean)),
                           s=c(apply(simkcf.trl.rea,2,nansd)/sqrt(apply(!is.na(simkcf.trl.rea),2,sum)),apply(simkcf.trl.bio,2,nansd)/sqrt(apply(!is.na(simkcf.trl.bio),2,sum)),
                               apply(simkcf.trl.ave,2,nansd)/sqrt(apply(!is.na(simkcf.trl.ave),2,sum))),
                           trl=rep(1:16,3), mdl=factor(c(rep('REA',16),rep('BIO',16),rep('MEAN',16)),levels=c('REA','BIO','MEAN')))
simdat.figs3e = simdat.figs3e[!is.na(simdat.figs3e$cf),]

figs3e = ggplot(data=dat.figs3e, aes(x=trl, y=cf)) +
        stat_boxplot(geom='errorbar', linewidth=0.6, width=0.3) +
        geom_boxplot(outlier.shape=3, outlier.size=1.2, fill="white") +
        geom_ribbon(data=simdat.figs3e, aes(x=trl, ymin=cf-s, ymax=cf+s, group=mdl, fill=mdl), colour=NA, alpha=0.5, show.legend=FALSE) +
        geom_line(data=simdat.figs3e, aes(x=trl, y=cf, group=mdl, color=mdl), linewidth=0.6, alpha=1) +
        scale_fill_manual(name=NULL, values=colors) +
        scale_color_manual(name=NULL, values=colors, labels=c(' REA',' BIO',' AVG')) +
        guides(color=guide_legend(override.aes=list(linewidth=1))) +
        scale_x_discrete(name="Trial index in a single game") +
        scale_y_continuous(limits=c(1,4), name="Averaged Deck confidence level", labels=scales::number_format(accuracy=0.1)) + 
        theme_custom +
        theme(legend.position=c(1,0.02),legend.justification=c(1,0),legend.background=element_rect(fill = "transparent"),
              legend.title=element_text(size=8),
              legend.text=element_text(size=8, margin=margin(r=0,l=0.1,unit='cm')),legend.spacing.x=unit(-0.05,'cm'),legend.spacing.y=unit(0.1,'cm'),
              legend.key.width=unit(0.3,"cm"), legend.key.height=unit(0.3,"cm"))
figs3e