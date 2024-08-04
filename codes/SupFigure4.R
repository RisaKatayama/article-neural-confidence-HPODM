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
mdl.folder = ""
mdlres.folder = ""
load(paste0(c(dat.folder,"bhvdata_bhv.RData"),collapse="/"))

load(paste0(c(mdlres.folder,"HDM_dimdl.RData"),collapse="/"))
res.h = res
load(paste0(c(mdlres.folder,"PDM_dimdl.RData"),collapse="/"))
res.p = res

inits1 = dump.format(list(.RNG.name="base::Super-Duper", .RNG.seed=111))
inits2 = dump.format(list(.RNG.name="base::Wichmann-Hill", .RNG.seed=222))
inits3 = dump.format(list(.RNG.name="base::Mersenne-Twister", .RNG.seed=333))

colors = c("#ff9999", "#cc0033", "#9999ff", "#3300cc", "#009900", "#99ee00", "#ffbb00", "#ff6600")

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
    kcf.hl = rep(NaN,N)
    for (i in 1:ns) {
        tmpkcf=kcf[data$sidx==i]
        tmpkcfhl = (tmpkcf>nanmean(tmpkcf))*1
        kcf.hl[data$sidx==i] = tmpkcfhl
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
                     deci=deci[include], open=open[include], icf=icf[include], icf.hl=icf.hl[include],
                     seqid=seqid[include], trl=trl[include], ses=ses[include], sidx=sidx[include])
    return(res)
}
data = load.data.basic(bdata)
ns = max(data$sidx)
ng = max(data$ses)
N = dim(data)[1]


############################
## Fig S4bc/ij
phit.icfdk = array(NaN,dim=c(ns,11,2,2))
for (i in 1:ns) {
    for (o in 8:18) {
        for (c in 0:1) {
            for (d in 0:1) {
                phit.icfdk[i,o-7,c+1,d+1] = mean(data$deci[data$open==o&data$dkest==d&data$icf.hl==c&data$sidx==i])
            }
        }
    }
}

sim.phit.icfdk.h = array(NaN,dim=c(ns,11,2,2))
sim.phit.icfdk.p = array(NaN,dim=c(ns,11,2,2))
for (i in 1:ns) {
    for (d in 0:1) {
        for (o in 8:18) {
            tmpicf = res.h$di.res$fitres$w.pred.icf[data$sidx==i]
            sim.phit.icfdk.h[i,o-7,1,d+1] = mean(res.h$di.res$fitres$pmat.di[data$sidx==i&data$dkest==d&data$open==o&res.h$di.res$fitres$w.pred.icf<=mean(tmpicf),2])
            sim.phit.icfdk.h[i,o-7,2,d+1] = mean(res.h$di.res$fitres$pmat.di[data$sidx==i&data$dkest==d&data$open==o&res.h$di.res$fitres$w.pred.icf>mean(tmpicf),2])

            tmpicf = res.p$di.res$fitres$w.pred.icf[data$sidx==i]
            sim.phit.icfdk.p[i,o-7,1,d+1] = mean(res.p$di.res$fitres$pmat.di[data$sidx==i&data$dkest==d&data$open==o&res.p$di.res$fitres$w.pred.icf<=mean(tmpicf),2])
            sim.phit.icfdk.p[i,o-7,2,d+1] = mean(res.p$di.res$fitres$pmat.di[data$sidx==i&data$dkest==d&data$open==o&res.p$di.res$fitres$w.pred.icf>mean(tmpicf),2])
        }
    }
}

dat.figs4b = data.frame(prp=array(apply(phit.icfdk,c(2,3,4),nanmean)),
                       s=array(apply(phit.icfdk,c(2,3,4),nansd)/sqrt(apply(!is.na(phit.icfdk),c(2,3,4),nansum))),
                       opn=rep(8:18,2*2), cf=factor(rep(c(rep(0,11),rep(1,11)),2)), dk=factor(c(rep(0,11*2),rep(1,11*2))))
dat.figs4b = dat.figs4b[!is.na(dat.figs4b$prp),]

simdat.figs4b = data.frame(prp=array(apply(sim.phit.icfdk.h,c(2,3,4),nanmean)),
                          s=array(apply(sim.phit.icfdk.h,c(2,3,4),nansd)/sqrt(apply(!is.na(sim.phit.icfdk.h),c(2,3,4),nansum))),
                          opn=rep(8:18,2*2), cf=factor(rep(c(rep(0,11),rep(1,11)),2)), dk=factor(c(rep(0,11*2),rep(1,11*2))))
simdat.figs4b = simdat.figs4b[!is.na(simdat.figs4b$prp),]

figs4b = ggplot(data=dat.figs4b, aes(x=opn, y=prp, group=interaction(cf,dk), fill=interaction(cf,dk), color=interaction(cf,dk))) +
        geom_errorbar(aes(x=opn, ymin=prp-s, ymax=prp+s, group=interaction(cf,dk), color=interaction(cf,dk)), width=0.15, alpha=0.5) +
        geom_point(aes(fill=interaction(cf,dk)), shape=21, color="white", size=2.5, alpha=0.75) +
        geom_ribbon(data=simdat.figs4b, aes(x=opn, ymin=prp-s, ymax=prp+s, group=interaction(cf,dk), fill=interaction(cf,dk)), colour=NA, alpha=0.3, show.legend=FALSE) +
        geom_line(data=simdat.figs4b, aes(x=opn, y=prp, group=interaction(cf,dk), color=interaction(cf,dk)), linewidth=0.7, alpha=0.9) +
        scale_fill_manual(name="Decision conf.\n\n", values=colors[1:4], labels=c("",""," Deck 4 "," Deck 6"),
                          guide=guide_legend(override.aes=list(color=colors[c(1,3,2,4)], shape=c(15,15,15,15), size=4.5, alpha=1),ncol=2)) +
        scale_color_manual(name=NULL, values=colors[1:4], guide="none") +
        scale_x_continuous(breaks=seq(8,18,2), name="Face-up card score") +
        scale_y_continuous(limits=c(-0.01,1.01), breaks=seq(0,1,0.2), name="Proportion of Hit") +
        theme_custom +
        theme(legend.position=c(1,1),legend.justification=c(1,1),legend.background=element_rect(colour="transparent"),
              legend.title=element_text(size=10),
              legend.text=element_text(size=10, margin=margin(r=0,unit='cm')),legend.spacing.x=unit(-0.02,'cm'),legend.spacing.y=unit(0,'cm'),
              legend.key.width=unit(0.5,"cm"), legend.key.height=unit(0.3,"cm"))
figs4b

dat = data.frame(p.hit=array(aperm(sim.phit.icfdk.h,c(2,3,4,1))), open=rep(8:18,ns*2*2), cf=rep(c(rep(0,11),rep(1,11)),ns*2),
                    dk=rep(c(rep(0,11*2),rep(1,11*2)),ns), sidx=array(t(array(rep(1:ns,11*2*2),dim=c(ns,11*2*2)))))
dat = dat[!is.na(dat$p.hit),]
dat = dump.format(list(open=dat$open, cf=dat$cf, dk=dat$dk, p.hit=dat$p.hit, sidx=dat$sidx, ns=ns, N=dim(dat)[1]))
vals = c("a.mu","b.mu","a.sig","b.sig","a.p","b.p")

rgr.simicfdk.h = run.jags(model=paste0(mdl.folder,"SimpleRegress_pseudoLogistic.txt"), monitor=vals, data=dat, n.chains=3, inits=c(inits1,inits2,inits3),
                     plots=FALSE, method="parallel", burnin=5000, sample=1000, adapt=1000, thin=10)
summary.rgr.simicfdk.h = add.summary(rgr.simicfdk.h)

s.simicfdk.h = array(sapply(1:4,function(c){sapply(1:ns,function(i){summary.rgr.simicfdk.h$summaries[paste0("a.p[",i,",",c,"]"),"Median"]})}),dim=c(ns,2,2))
b.simicfdk.h = array(sapply(1:4,function(c){sapply(1:ns,function(i){summary.rgr.simicfdk.h$summaries[paste0("b.p[",i,",",c,"]"),"Median"]})}),dim=c(ns,2,2))

dat.figs4d = data.frame(slp=array(aperm(s.simicfdk.h,c(2,3,1))), sft=array(aperm(b.simicfdk.h,c(2,3,1))), cf=factor(rep(c('Low','High','Low','High'),ns),levels=c('Low','High')),
                       dk=factor(rep(c('Deck 4','Deck 4','Deck 6','Deck 6'),ns)), s=factor(array(t(array(rep(1:ns,2*2),dim=c(ns,2*2))))))

figs4dl = ggplot(data=dat.figs4d, aes(x=interaction(cf,dk), y=slp, fill=interaction(cf,dk))) +
        stat_boxplot(geom='errorbar',size=0.5, width=0.25) +
        geom_boxplot(aes(fill=interaction(cf,dk)), outlier.shape=3, outlier.size=1.2) +
        scale_fill_manual(values=colors[1:4]) +
        scale_x_discrete(guide="axis_nested") +
        scale_y_continuous(limits=c(-3,0), breaks=-5:1, labels=scales::number_format(accuracy=0.1), name="Slope") + 
        geom_signif(comparisons=list(c(1,2),c(3,4)), annotations=c('***'), y_position=-0.3, tip_length=0, vjust=0.2, textsize=4.5) +
        theme_custom +
        theme(axis.title.x=element_blank(), axis.text.x=element_text(angle=30, hjust=1, vjust=1, margin=margin(t=0.05, b=0, unit="cm")),
              ggh4x.axis.nestline.x=element_line(linewidth=0.25), ggh4x.axis.nesttext.x=element_text(angle=30, hjust=0.5, vjust=0.5, margin=margin(t=0.09, unit="cm")))
figs4dl

figs4dr = ggplot(data=dat.figs4d, aes(x=interaction(cf,dk), y=sft, fill=interaction(cf,dk))) +
        stat_boxplot(geom='errorbar',size=0.5, width=0.25) +
        geom_boxplot(aes(fill=interaction(cf,dk)), outlier.shape=3, outlier.size=1.2) +
        scale_fill_manual(values=colors[1:4]) +
        scale_x_discrete(guide="axis_nested") +
        scale_y_continuous(limits=c(9,15.1), breaks=seq(9,15,1.5), name="Shift") + 
        geom_signif(comparisons=list(c(1,3),c(2,4)), annotations=c('***'), y_position=c(9.6,9.1), tip_length=0, vjust=1.8, textsize=4.5) +
        theme_custom +
        theme(axis.title.x=element_blank(), axis.text.x=element_text(angle=30, hjust=1, vjust=1, margin=margin(t=0.05, b=0, unit="cm")),
              ggh4x.axis.nestline.x=element_line(linewidth=0.25), ggh4x.axis.nesttext.x=element_text(angle=30, hjust=0.5, vjust=0.5, margin=margin(t=0.09, unit="cm")))
figs4dr

##
simdat.figs4i = data.frame(prp=array(apply(sim.phit.icfdk.p,c(2,3,4),nanmean)),
                          s=array(apply(sim.phit.icfdk.p,c(2,3,4),nansd)/sqrt(apply(!is.na(sim.phit.icfdk.p),c(2,3,4),nansum))),
                          opn=rep(8:18,2*2), cf=factor(rep(c(rep(0,11),rep(1,11)),2)), dk=factor(c(rep(0,11*2),rep(1,11*2))))
simdat.figs4i = simdat.figs4i[!is.na(simdat.figs4i$prp),]

figs4i = ggplot(data=dat.figs4b, aes(x=opn, y=prp, group=interaction(cf,dk), fill=interaction(cf,dk), color=interaction(cf,dk))) +
        geom_errorbar(aes(x=opn, ymin=prp-s, ymax=prp+s, group=interaction(cf,dk), color=interaction(cf,dk)), width=0.15, alpha=0.5) +
        geom_point(aes(fill=interaction(cf,dk)), shape=21, color="white", size=2.5, alpha=0.75) +
        geom_ribbon(data=simdat.figs4i, aes(x=opn, ymin=prp-s, ymax=prp+s, group=interaction(cf,dk), fill=interaction(cf,dk)), colour=NA, alpha=0.3, show.legend=FALSE) +
        geom_line(data=simdat.figs4i, aes(x=opn, y=prp, group=interaction(cf,dk), color=interaction(cf,dk)), linewidth=0.7, alpha=0.9) +
        scale_fill_manual(name="Decision conf.\n\n", values=colors[1:4], labels=c("",""," Deck 4 "," Deck 6"),
                          guide=guide_legend(override.aes=list(color=colors[c(1,3,2,4)], shape=c(15,15,15,15), size=4.5, alpha=1),ncol=2)) +
        scale_color_manual(name=NULL, values=colors[1:4], guide="none") +
        scale_x_continuous(breaks=seq(8,18,2), name="Face-up card score") +
        scale_y_continuous(limits=c(-0.01,1.01), breaks=seq(0,1,0.2), name="Proportion of Hit") +
        theme_custom +
        theme(legend.position=c(1,1),legend.justification=c(1,1),legend.background=element_rect(colour="transparent"),
              legend.title=element_text(size=10),
              legend.text=element_text(size=10, margin=margin(r=0,unit='cm')),legend.spacing.x=unit(-0.02,'cm'),legend.spacing.y=unit(0,'cm'),
              legend.key.width=unit(0.5,"cm"), legend.key.height=unit(0.3,"cm"))
figs4i

dat = data.frame(p.hit=array(aperm(sim.phit.icfdk.p,c(2,3,4,1))), open=rep(8:18,ns*2*2), cf=rep(c(rep(0,11),rep(1,11)),ns*2),
                    dk=rep(c(rep(0,11*2),rep(1,11*2)),ns), sidx=array(t(array(rep(1:ns,11*2*2),dim=c(ns,11*2*2)))))
dat = dat[!is.na(dat$p.hit),]
dat = dump.format(list(open=dat$open, cf=dat$cf, dk=dat$dk, p.hit=dat$p.hit, sidx=dat$sidx, ns=ns, N=dim(dat)[1]))
vals = c("a.mu","b.mu","a.sig","b.sig","a.p","b.p")

rgr.simicfdk.p = run.jags(model=paste0(mdl.folder,"SimpleRegress_pseudoLogistic.txt"), monitor=vals, data=dat, n.chains=3, inits=c(inits1,inits2,inits3),
                     plots=FALSE, method="parallel", burnin=5000, sample=1000, adapt=1000, thin=10)
summary.rgr.simicfdk.p = add.summary(rgr.simicfdk.p)

s.simicfdk.p = array(sapply(1:4,function(c){sapply(1:ns,function(i){summary.rgr.simicfdk.p$summaries[paste0("a.p[",i,",",c,"]"),"Median"]})}),dim=c(ns,2,2))
b.simicfdk.p = array(sapply(1:4,function(c){sapply(1:ns,function(i){summary.rgr.simicfdk.p$summaries[paste0("b.p[",i,",",c,"]"),"Median"]})}),dim=c(ns,2,2))

dat.figs4j = data.frame(slp=array(aperm(s.simicfdk.p,c(2,3,1))), sft=array(aperm(b.simicfdk.p,c(2,3,1))), cf=factor(rep(c('Low','High','Low','High'),ns),levels=c('Low','High')),
                       dk=factor(rep(c('Deck 4','Deck 4','Deck 6','Deck 6'),ns)), s=factor(array(t(array(rep(1:ns,2*2),dim=c(ns,2*2))))))

figs4jl = ggplot(data=dat.figs4j, aes(x=interaction(cf,dk), y=slp, fill=interaction(cf,dk))) +
        stat_boxplot(geom='errorbar',size=0.5, width=0.25) +
        geom_boxplot(aes(fill=interaction(cf,dk)), outlier.shape=3, outlier.size=1.2) +
        scale_fill_manual(values=colors[1:4]) +
        scale_x_discrete(guide="axis_nested") +
        scale_y_continuous(limits=c(-3,0), breaks=-5:1, labels=scales::number_format(accuracy=0.1), name="Slope") + 
        geom_signif(comparisons=list(c(3,4)), annotations=c('***'), y_position=-0.4, tip_length=0, vjust=-0.2, textsize=4.5) +
        geom_signif(comparisons=list(c(2,4)), annotations=c('*'), y_position=-2.1, tip_length=0, vjust=2, textsize=4.5) +
        theme_custom +
        theme(axis.title.x=element_blank(), axis.text.x=element_text(angle=30, hjust=1, vjust=1, margin=margin(t=0.05, b=0, unit="cm")),
              ggh4x.axis.nestline.x=element_line(linewidth=0.25), ggh4x.axis.nesttext.x=element_text(angle=30, hjust=0.5, vjust=0.5, margin=margin(t=0.09, unit="cm")))
figs4jl

figs4jr = ggplot(data=dat.figs4j, aes(x=interaction(cf,dk), y=sft, fill=interaction(cf,dk))) +
        stat_boxplot(geom='errorbar',size=0.5, width=0.25) +
        geom_boxplot(aes(fill=interaction(cf,dk)), outlier.shape=3, outlier.size=1.2) +
        scale_fill_manual(values=colors[1:4]) +
        scale_x_discrete(guide="axis_nested") +
        scale_y_continuous(limits=c(9,15.2), breaks=seq(9,15,1.5), name="Shift") + 
        geom_signif(comparisons=list(c(1,2),c(3,4)), annotations=c('***'), y_position=14.5, tip_length=0, vjust=-0.2, textsize=4.5) +
        geom_signif(comparisons=list(c(1,3),c(2,4)), annotations=c('***'), y_position=c(11,10.3), tip_length=0, vjust=1.8, textsize=4.5) +
        theme_custom +
        theme(axis.title.x=element_blank(), axis.text.x=element_text(angle=30, hjust=1, vjust=1, margin=margin(t=0.05, b=0, unit="cm")),
              ggh4x.axis.nestline.x=element_line(linewidth=0.25), ggh4x.axis.nesttext.x=element_text(angle=30, hjust=0.5, vjust=0.5, margin=margin(t=0.09, unit="cm")))
figs4jr


## Fig S4de/kl
phit.kcfdk = array(NaN,dim=c(ns,11,2,2))
for (i in 1:ns) {
    for (o in 8:18) {
        for (c in 0:1) {
            for (d in 0:1) {
                phit.kcfdk[i,o-7,c+1,d+1] = mean(data$deci[data$open==o&data$dkest==d&data$kcf.hl==c&data$sidx==i])
            }
        }
    }
}

sim.phit.kcfdk.h = array(NaN,dim=c(ns,11,2,2))
sim.phit.kcfdk.p = array(NaN,dim=c(ns,11,2,2))
for (i in 1:ns) {
    for (d in 0:1) {
        for (o in 8:18) {
            tmpkcf = res.h$dk.res$fitres$w.pred.kcf[data$sidx==i]
            sim.phit.kcfdk.h[i,o-7,1,d+1] = mean(res.h$di.res$fitres$pmat.di[data$sidx==i&data$dkest==d&data$open==o&res.h$dk.res$fitres$w.pred.kcf<=mean(tmpkcf),2])
            sim.phit.kcfdk.h[i,o-7,2,d+1] = mean(res.h$di.res$fitres$pmat.di[data$sidx==i&data$dkest==d&data$open==o&res.h$dk.res$fitres$w.pred.kcf>mean(tmpkcf),2])

            tmpkcf = res.p$dk.res$fitres$w.pred.kcf[data$sidx==i]
            sim.phit.kcfdk.p[i,o-7,1,d+1] = mean(res.p$di.res$fitres$pmat.di[data$sidx==i&data$dkest==d&data$open==o&res.p$dk.res$fitres$w.pred.kcf<=mean(tmpkcf),2])
            sim.phit.kcfdk.p[i,o-7,2,d+1] = mean(res.p$di.res$fitres$pmat.di[data$sidx==i&data$dkest==d&data$open==o&res.p$dk.res$fitres$w.pred.kcf>mean(tmpkcf),2])
        }
    }
}

dat.figs4d = data.frame(prp=array(apply(phit.kcfdk,c(2,3,4),nanmean)),
                       s=array(apply(phit.kcfdk,c(2,3,4),nansd)/sqrt(apply(!is.na(phit.kcfdk),c(2,3,4),nansum))),
                       opn=rep(8:18,2*2), cf=factor(rep(c(rep(0,11),rep(1,11)),2)), dk=factor(c(rep(0,11*2),rep(1,11*2))))
dat.figs4d = dat.figs4d[!is.na(dat.figs4d$prp),]

simdat.figs4d = data.frame(prp=array(apply(sim.phit.kcfdk.h,c(2,3,4),nanmean)),
                          s=array(apply(sim.phit.kcfdk.h,c(2,3,4),nansd)/sqrt(apply(!is.na(sim.phit.kcfdk.h),c(2,3,4),nansum))),
                          opn=rep(8:18,2*2), cf=factor(rep(c(rep(0,11),rep(1,11)),2)), dk=factor(c(rep(0,11*2),rep(1,11*2))))
simdat.figs4d = simdat.figs4d[!is.na(simdat.figs4d$prp),]

figs4d = ggplot(data=dat.figs4d, aes(x=opn, y=prp, group=interaction(cf,dk), fill=interaction(cf,dk), color=interaction(cf,dk))) +
        geom_errorbar(aes(x=opn, ymin=prp-s, ymax=prp+s, group=interaction(cf,dk), color=interaction(cf,dk)), width=0.15, alpha=0.5) +
        geom_point(aes(fill=interaction(cf,dk)), shape=21, color="white", size=2.5, alpha=0.75) +
        geom_ribbon(data=simdat.figs4d, aes(x=opn, ymin=prp-s, ymax=prp+s, group=interaction(cf,dk), fill=interaction(cf,dk)), colour=NA, alpha=0.3, show.legend=FALSE) +
        geom_line(data=simdat.figs4d, aes(x=opn, y=prp, group=interaction(cf,dk), color=interaction(cf,dk)), linewidth=0.7, alpha=0.9) +
        scale_fill_manual(name="Decision conf.\n\n", values=colors[1:4], labels=c("",""," Deck 4 "," Deck 6"),
                          guide=guide_legend(override.aes=list(color=colors[c(1,3,2,4)], shape=c(15,15,15,15), size=4.5, alpha=1),ncol=2)) +
        scale_color_manual(name=NULL, values=colors[1:4], guide="none") +
        scale_x_continuous(breaks=seq(8,18,2), name="Face-up card score") +
        scale_y_continuous(limits=c(-.01,1.01), breaks=seq(0,1,0.2), name="Proportion of Hit") +
        theme_custom +
        theme(legend.position=c(1,1),legend.justification=c(1,1),legend.background=element_rect(colour="transparent"),
              legend.title=element_text(size=10),
              legend.text=element_text(size=10, margin=margin(r=0,unit='cm')),legend.spacing.x=unit(-0.02,'cm'),legend.spacing.y=unit(0,'cm'),
              legend.key.width=unit(0.5,"cm"), legend.key.height=unit(0.3,"cm"))
figs4d

dat = data.frame(p.hit=array(aperm(sim.phit.kcfdk.h,c(2,3,4,1))), open=rep(8:18,ns*2*2), cf=rep(c(rep(0,11),rep(1,11)),ns*2),
                 dk=rep(c(rep(0,11*2),rep(1,11*2)),ns), sidx=array(t(array(rep(1:ns,11*2*2),dim=c(ns,11*2*2)))))
dat = dat[!is.na(dat$p.hit),]
dat = dump.format(list(open=dat$open, cf=dat$cf, dk=dat$dk, p.hit=dat$p.hit, sidx=dat$sidx, ns=ns, N=dim(dat)[1]))
vals = c("a.mu","b.mu","a.sig","b.sig","a.p","b.p")

rgr.simkcfdk.h = run.jags(model=paste0(mdl.folder,"SimpleRegress_pseudoLogistic.txt"), monitor=vals, data=dat, n.chains=3, inits=c(inits1,inits2,inits3),
                     plots=FALSE, method="parallel", burnin=5000, sample=1000, adapt=1000, thin=10)
summary.rgr.simkcfdk.h = add.summary(rgr.simkcfdk.h)

s.simkcfdk.h = array(sapply(1:4,function(c){sapply(1:ns,function(i){summary.rgr.simkcfdk.h$summaries[paste0("a.p[",i,",",c,"]"),"Median"]})}),dim=c(ns,2,2))
b.simkcfdk.h = array(sapply(1:4,function(c){sapply(1:ns,function(i){summary.rgr.simkcfdk.h$summaries[paste0("b.p[",i,",",c,"]"),"Median"]})}),dim=c(ns,2,2))

dat.figs4e = data.frame(slp=array(aperm(s.simkcfdk.h,c(2,3,1))), sft=array(aperm(b.simkcfdk.h,c(2,3,1))), cf=factor(rep(c('Low','High','Low','High'),ns),levels=c('Low','High')),
                       dk=factor(rep(c('Deck 4','Deck 4','Deck 6','Deck 6'),ns)), s=factor(array(t(array(rep(1:ns,2*2),dim=c(ns,2*2))))))

figs4el = ggplot(data=dat.figs4e, aes(x=interaction(cf,dk), y=slp, fill=interaction(cf,dk))) +
        stat_boxplot(geom='errorbar',size=0.5, width=0.25) +
        geom_boxplot(aes(fill=interaction(cf,dk)), outlier.shape=3, outlier.size=1.2) +
        scale_fill_manual(values=colors[1:4]) +
        scale_x_discrete(guide="axis_nested") +
        scale_y_continuous(limits=c(-3,0), breaks=-5:1, labels=scales::number_format(accuracy=0.1), name="Slope") + 
        theme_custom +
        theme(axis.title.x=element_blank(), axis.text.x=element_text(angle=30, hjust=1, vjust=1, margin=margin(t=0.05, b=0, unit="cm")),
              ggh4x.axis.nestline.x=element_line(linewidth=0.25), ggh4x.axis.nesttext.x=element_text(angle=30, hjust=0.5, vjust=0.5, margin=margin(t=0.09, unit="cm")))
figs4el

figs4er = ggplot(data=dat.figs4e, aes(x=interaction(cf,dk), y=sft, fill=interaction(cf,dk))) +
        stat_boxplot(geom='errorbar',size=0.5, width=0.25) +
        geom_boxplot(aes(fill=interaction(cf,dk)), outlier.shape=3, outlier.size=1.2) +
        scale_fill_manual(values=colors[1:4]) +
        scale_x_discrete(guide="axis_nested") +
        scale_y_continuous(limits=c(9,15.2), breaks=seq(9,15,1.5), name="Shift") + 
        geom_signif(comparisons=list(c(1,3),c(2,4)), annotations=c('***'), y_position=c(9.6,9.1), tip_length=0, vjust=1.8, textsize=4.5) +
        theme_custom +
        theme(axis.title.x=element_blank(), axis.text.x=element_text(angle=30, hjust=1, vjust=1, margin=margin(t=0.05, b=0, unit="cm")),
              ggh4x.axis.nestline.x=element_line(linewidth=0.25), ggh4x.axis.nesttext.x=element_text(angle=30, hjust=0.5, vjust=0.5, margin=margin(t=0.09, unit="cm")))
figs4er


##
simdat.figs4k = data.frame(prp=array(apply(sim.phit.kcfdk.p,c(2,3,4),nanmean)),
                          s=array(apply(sim.phit.kcfdk.p,c(2,3,4),nansd)/sqrt(apply(!is.na(sim.phit.kcfdk.p),c(2,3,4),nansum))),
                          opn=rep(8:18,2*2), cf=factor(rep(c(rep(0,11),rep(1,11)),2)), dk=factor(c(rep(0,11*2),rep(1,11*2))))
simdat.figs4k = simdat.figs4k[!is.na(simdat.figs4k$prp),]

figs4k = ggplot(data=dat.figs4d, aes(x=opn, y=prp, group=interaction(cf,dk), fill=interaction(cf,dk), color=interaction(cf,dk))) +
        geom_errorbar(aes(x=opn, ymin=prp-s, ymax=prp+s, group=interaction(cf,dk), color=interaction(cf,dk)), width=0.15, alpha=0.5) +
        geom_point(aes(fill=interaction(cf,dk)), shape=21, color="white", size=2.5, alpha=0.75) +
        geom_ribbon(data=simdat.figs4k, aes(x=opn, ymin=prp-s, ymax=prp+s, group=interaction(cf,dk), fill=interaction(cf,dk)), colour=NA, alpha=0.3, show.legend=FALSE) +
        geom_line(data=simdat.figs4k, aes(x=opn, y=prp, group=interaction(cf,dk), color=interaction(cf,dk)), linewidth=0.7, alpha=0.9) +
        scale_fill_manual(name="Decision conf.\n\n", values=colors[1:4], labels=c("",""," Deck 4 "," Deck 6"),
                          guide=guide_legend(override.aes=list(color=colors[c(1,3,2,4)], shape=c(15,15,15,15), size=4.5, alpha=1),ncol=2)) +
        scale_color_manual(name=NULL, values=colors[1:4], guide="none") +
        scale_x_continuous(breaks=seq(8,18,2), name="Open card score") +
        scale_y_continuous(limits=c(-.01,1.01), breaks=seq(0,1,0.2), name="Proportion of Hit") +
        theme_custom +
        theme(legend.position=c(1,1),legend.justification=c(1,1),legend.background=element_rect(colour="transparent"),
              legend.title=element_text(size=10),
              legend.text=element_text(size=10, margin=margin(r=0,unit='cm')),legend.spacing.x=unit(-0.02,'cm'),legend.spacing.y=unit(0,'cm'),
              legend.key.width=unit(0.5,"cm"), legend.key.height=unit(0.3,"cm"))
figs4k

dat = data.frame(p.hit=array(aperm(sim.phit.kcfdk.p,c(2,3,4,1))), open=rep(8:18,ns*2*2), cf=rep(c(rep(0,11),rep(1,11)),ns*2),
                    dk=rep(c(rep(0,11*2),rep(1,11*2)),ns), sidx=array(t(array(rep(1:ns,11*2*2),dim=c(ns,11*2*2)))))
dat = dat[!is.na(dat$p.hit),]
dat = dump.format(list(open=dat$open, cf=dat$cf, dk=dat$dk, p.hit=dat$p.hit, sidx=dat$sidx, ns=ns, N=dim(dat)[1]))
vals = c("a.mu","b.mu","a.sig","b.sig","a.p","b.p")

rgr.simkcfdk.p = run.jags(model=paste0(mdl.folder,"SimpleRegress_pseudoLogistic.txt"), monitor=vals, data=dat, n.chains=3, inits=c(inits1,inits2,inits3),
                     plots=FALSE, method="parallel", burnin=5000, sample=1000, adapt=1000, thin=10)
summary.rgr.simkcfdk.p = add.summary(rgr.simkcfdk.p)

s.simkcfdk.p = array(sapply(1:4,function(c){sapply(1:ns,function(i){summary.rgr.simkcfdk.p$summaries[paste0("a.p[",i,",",c,"]"),"Median"]})}),dim=c(ns,2,2))
b.simkcfdk.p = array(sapply(1:4,function(c){sapply(1:ns,function(i){summary.rgr.simkcfdk.p$summaries[paste0("b.p[",i,",",c,"]"),"Median"]})}),dim=c(ns,2,2))

dat.figs4l = data.frame(slp=array(aperm(s.simkcfdk.p,c(2,3,1))), sft=array(aperm(b.simkcfdk.p,c(2,3,1))), cf=factor(rep(c('Low','High','Low','High'),ns),levels=c('Low','High')),
                       dk=factor(rep(c('Deck 4','Deck 4','Deck 6','Deck 6'),ns)), s=factor(array(t(array(rep(1:ns,2*2),dim=c(ns,2*2))))))

figs4ll = ggplot(data=dat.figs4l, aes(x=interaction(cf,dk), y=slp, fill=interaction(cf,dk))) +
        stat_boxplot(geom='errorbar',size=0.5, width=0.25) +
        geom_boxplot(aes(fill=interaction(cf,dk)), outlier.shape=3, outlier.size=1.2) +
        scale_fill_manual(values=colors[1:4]) +
        scale_x_discrete(guide="axis_nested") +
        scale_y_continuous(limits=c(-3,0), breaks=-5:1, labels=scales::number_format(accuracy=0.1), name="Slope") + 
        geom_signif(comparisons=list(c(3,4)), annotations=c('***'), y_position=-0.4, tip_length=0, vjust=-0.2, textsize=4.5) +
        geom_signif(comparisons=list(c(2,4)), annotations=c('*'), y_position=-2.1, tip_length=0, vjust=2, textsize=4.5) +
        theme_custom +
        theme(axis.title.x=element_blank(), axis.text.x=element_text(angle=30, hjust=1, vjust=1, margin=margin(t=0.05, b=0, unit="cm")),
              ggh4x.axis.nestline.x=element_line(linewidth=0.25), ggh4x.axis.nesttext.x=element_text(angle=30, hjust=0.5, vjust=0.5, margin=margin(t=0.09, unit="cm")))
figs4ll

figs4lr = ggplot(data=dat.figs4l, aes(x=interaction(cf,dk), y=sft, fill=interaction(cf,dk))) +
        stat_boxplot(geom='errorbar',size=0.5, width=0.25) +
        geom_boxplot(aes(fill=interaction(cf,dk)), outlier.shape=3, outlier.size=1.2) +
        scale_fill_manual(values=colors[1:4]) +
        scale_x_discrete(guide="axis_nested") +
        scale_y_continuous(limits=c(9,15.2), breaks=seq(9,15,1.5), name="Bias") + 
        geom_signif(comparisons=list(c(1,2),c(3,4)), annotations=c('***'), y_position=14.5, tip_length=0, vjust=-0.2, textsize=4.5) +
        geom_signif(comparisons=list(c(1,3),c(2,4)), annotations=c('***'), y_position=c(11,10.3), tip_length=0, vjust=1.8, textsize=4.5) +
        theme_custom +
        theme(axis.title.x=element_blank(), axis.text.x=element_text(angle=30, hjust=1, vjust=1, margin=margin(t=0.05, b=0, unit="cm")),
              ggh4x.axis.nestline.x=element_line(linewidth=0.25), ggh4x.axis.nesttext.x=element_text(angle=30, hjust=0.5, vjust=0.5, margin=margin(t=0.09, unit="cm")))
figs4lr


## Fig S4fg/mn
icf.kcf = array(NaN, dim=c(ns,11,4))
for (i in 1:ns) {
    for (o in 8:18) {
        for (c in 1:4) {
            icf.kcf[i,o-7,c] = mean(data$icf[data$sidx==i&data$open==o&data$kcf==c])
        }
    }
}

sim.icf.kcf.h = array(NaN,dim=c(ns,11,4))
sim.icf.kcf.p = array(NaN,dim=c(ns,11,4))
for (i in 1:ns){
    for (c in 1:4) {
        for (o in 8:18) {
           sim.icf.kcf.h[i,o-7,c] = mean(res.h$di.res$fitres$w.pred.icf[data$sidx==i&data$open==o&res.h$dk.res$fitres$pred.kcf==c])
           sim.icf.kcf.p[i,o-7,c] = mean(res.p$di.res$fitres$w.pred.icf[data$sidx==i&data$open==o&res.p$dk.res$fitres$pred.kcf==c])
        }
    }
}

dat.figs4f = data.frame(icf=array(apply(icf.kcf,c(2,3),nanmean)),
                       s=array(apply(icf.kcf,c(2,3),nansd)/sqrt(apply(!is.na(icf.kcf),c(2,3),nansum))),
                       opn=rep(8:18,4), cf=factor(c(rep(1,11),rep(2,11),rep(3,11),rep(4,11))))
dat.figs4f = dat.figs4f[!is.na(dat.figs4f$icf),]

simdat.figs4f = data.frame(icf=array(apply(sim.icf.kcf.h,c(2,3),nanmean)),
                          s=array(apply(sim.icf.kcf.h,c(2,3),nansd)/sqrt(apply(!is.na(sim.icf.kcf.h),c(2,3),nansum))),
                          opn=rep(8:18,4), cf=factor(c(rep(1,11),rep(2,11),rep(3,11),rep(4,11))))
simdat.figs4f = simdat.figs4f[!is.na(simdat.figs4f$icf),]

figs4f = ggplot(data=dat.figs4f, aes(x=opn, y=icf, group=cf, fill=cf, color=cf)) +
        geom_errorbar(aes(x=opn, ymin=icf-s, ymax=icf+s, group=cf, color=cf), width=0.15, alpha=0.5) +
        geom_point(aes(fill=cf), shape=21, color="white", size=2.5, alpha=0.75) +
        geom_ribbon(data=simdat.figs4f, aes(x=opn, ymin=icf-s, ymax=icf+s, group=cf, fill=cf), colour=NA, alpha=0.3, show.legend=FALSE) +
        geom_line(data=simdat.figs4f, aes(x=opn, y=icf, group=cf, color=cf), linewidth=0.7, alpha=0.9) +
        scale_fill_manual(name="Deck conf.", values=colors[5:8], labels=c("1","2","3","4"),
                          guide=guide_legend(override.aes=list(color=NA, shape=22, size=4, alpha=1),ncol=4,label.position="bottom",title.hjust=0.5)) +
        scale_color_manual(name=NULL, values=colors[5:8], guide="none") +
        scale_x_continuous(breaks=seq(8,18,2), name="Face-up card score") +
        scale_y_continuous(limits=c(1,4.05), breaks=1:4, labels=scales::number_format(accuracy=0.1), name="Averaged Decision\nconfidence level") +
        theme_custom +
        theme(legend.position=c(1,0.05),legend.justification=c(1,0),,legend.background=element_rect(colour="transparent"),
              legend.title=element_text(size=10),
              legend.text=element_text(size=10, margin=margin(r=0,unit='cm')),legend.spacing.x=unit(-0.05,'cm'),legend.spacing.y=unit(0.1,'cm'),
              legend.key.width=unit(0.5,"cm"), legend.key.height=unit(0.3,"cm"))
figs4f

sim.phit.kcf.h = array(NaN,dim=c(ns,11,4))
sim.phit.kcf.p = array(NaN,dim=c(ns,11,4))
for (i in 1:ns){
    for (o in 8:18) {
        for (c in 1:4) {
            sim.phit.kcf.h[i,o-7,c] = mean(res.h$di.res$fitres$pmat.di[data$sidx==i&data$open==o&res.h$dk.res$fitres$pred.kcf==c,2])
            sim.phit.kcf.p[i,o-7,c] = mean(res.p$di.res$fitres$pmat.di[data$sidx==i&data$open==o&res.p$dk.res$fitres$pred.kcf==c,2])
        }
    }
}

dat = data.frame(p.hit=array(aperm(sim.phit.kcf.h,c(2,3,1))), icf=array(aperm(sim.icf.kcf.h,c(2,3,1))), open=rep(8:18,ns*4), 
                 kcf=rep(c(rep(1,11),rep(2,11),rep(3,11),rep(4,11)),ns), sidx=array(t(array(rep(1:ns,11*2*2),dim=c(ns,11*2*2)))))
dat = dat[!is.na(dat$p.hit)&!is.na(dat$icf),]
dat$sv = abs(logit(dat$p.hit))
dat = dump.format(list(sv=dat$sv, kcf=dat$kcf, icf=dat$icf, sidx=dat$sidx, ns=ns, N=dim(dat)[1]))
vals = c("a.mu","b.mu","a.sig","b.sig","a.p","b.p","sig.mu")

rgr.simkcf.h = run.jags(model=paste0(mdl.folder,"SimpleRegress.txt"), monitor=vals, data=dat, n.chains=3, inits=c(inits1,inits2,inits3),
                     plots=FALSE, method="parallel", burnin=5000, sample=1000, adapt=1000, thin=10)
summary.rgr.simkcf.h = add.summary(rgr.simkcf.h)

s.simicf.kcf.h = sapply(1:4,function(c){sapply(1:ns,function(i){summary.rgr.simkcf.h$summaries[paste0("a.p[",i,",",c,"]"),"Median"]})})
b.simicf.kcf.h = sapply(1:4,function(c){sapply(1:ns,function(i){summary.rgr.simkcf.h$summaries[paste0("b.p[",i,",",c,"]"),"Median"]})})

dat.figs3g = data.frame(slp=array(t(s.simicf.kcf.h)), bis=array(t(b.simicf.kcf.h)), cf=factor(rep(1:4,ns)), s=factor(array(t(array(rep(1:ns,4),dim=c(ns,4))))))

figs3gl = ggplot(data=dat.figs3g, aes(x=cf, y=slp, fill=cf)) +
        stat_boxplot(geom='errorbar',size=0.5, width=0.25) +
        geom_boxplot(aes(fill=cf), outlier.shape=3, outlier.size=1.2) +
        scale_fill_manual(values=colors[5:8]) +
        scale_x_discrete(name="Deck confidence\nlevel") +
        scale_y_continuous(limits=c(-0,0.4), breaks=seq(0,0.4,0.1), name="Slope") + 
        theme_custom
figs3gl

figs3gr = ggplot(data=dat.figs3g, aes(x=cf, y=bis, fill=cf)) +
        stat_boxplot(geom='errorbar',size=0.5, width=0.25) +
        geom_boxplot(aes(fill=cf), outlier.shape=3, outlier.size=1.2) +
        scale_fill_manual(values=colors[5:8]) +
        scale_x_discrete(name="Deck confidence\nlevel") +
        scale_y_continuous(limits=c(1,4), labels=scales::number_format(accuracy=0.1), name="Bias") + 
        theme_custom
figs3gr


##
simdat.figs4m = data.frame(icf=array(apply(sim.icf.kcf.p,c(2,3),nanmean)),
                          s=array(apply(sim.icf.kcf.p,c(2,3),nansd)/sqrt(apply(!is.na(sim.icf.kcf.p),c(2,3),nansum))),
                          opn=rep(8:18,4), cf=factor(c(rep(1,11),rep(2,11),rep(3,11),rep(4,11))))
simdat.figs4m = simdat.figs4m[!is.na(simdat.figs4m$icf),]

figs4m = ggplot(data=dat.figs4f, aes(x=opn, y=icf, group=cf, fill=cf, color=cf)) +
        geom_errorbar(aes(x=opn, ymin=icf-s, ymax=icf+s, group=cf, color=cf), width=0.15, alpha=0.5) +
        geom_point(aes(fill=cf), shape=21, color="white", size=2.5, alpha=0.75) +
        geom_ribbon(data=simdat.figs4m, aes(x=opn, ymin=icf-s, ymax=icf+s, group=cf, fill=cf), colour=NA, alpha=0.3, show.legend=FALSE) +
        geom_line(data=simdat.figs4m, aes(x=opn, y=icf, group=cf, color=cf), linewidth=0.7, alpha=0.9) +
        scale_fill_manual(name="Deck conf.", values=colors[5:8], labels=c("1","2","3","4"),
                          guide=guide_legend(override.aes=list(color=NA, shape=22, size=4, alpha=1),ncol=4,label.position="bottom",title.hjust=0.5)) +
        scale_color_manual(name=NULL, values=colors[5:8], guide="none") +
        scale_x_continuous(breaks=seq(8,18,2), name="Face-up card score") +
        scale_y_continuous(limits=c(1,4.05), breaks=1:4, labels=scales::number_format(accuracy=0.1), name="Averaged Decision\nconfidence level") +
        theme_custom +
        theme(legend.position=c(1,0.05),legend.justification=c(1,0),,legend.background=element_rect(colour="transparent"),
              legend.title=element_text(size=10),
              legend.text=element_text(size=10, margin=margin(r=0,unit='cm')),legend.spacing.x=unit(-0.05,'cm'),legend.spacing.y=unit(0.1,'cm'),
              legend.key.width=unit(0.5,"cm"), legend.key.height=unit(0.3,"cm"))
figs4m

dat = data.frame(p.hit=array(aperm(sim.phit.kcf.p,c(2,3,1))), icf=array(aperm(sim.icf.kcf.p,c(2,3,1))), open=rep(8:18,ns*4), 
                 kcf=rep(c(rep(1,11),rep(2,11),rep(3,11),rep(4,11)),ns), sidx=array(t(array(rep(1:ns,11*2*2),dim=c(ns,11*2*2)))))
dat = dat[!is.na(dat$p.hit)&!is.na(dat$icf),]
dat$sv = abs(logit(dat$p.hit))
dat = dump.format(list(sv=dat$sv, kcf=dat$kcf, icf=dat$icf, sidx=dat$sidx, ns=ns, N=dim(dat)[1]))
vals = c("a.mu","b.mu","a.sig","b.sig","a.p","b.p","sig.mu")

rgr.simkcf.p = run.jags(model=paste0(mdl.folder,"SimpleRegress.txt"), monitor=vals, data=dat, n.chains=3, inits=c(inits1,inits2,inits3),
                     plots=FALSE, method="parallel", burnin=5000, sample=1000, adapt=1000, thin=10)
summary.rgr.simkcf.p = add.summary(rgr.simkcf.p)

s.simicf.kcf.p = sapply(1:4,function(c){sapply(1:ns,function(i){summary.rgr.simkcf.p$summaries[paste0("a.p[",i,",",c,"]"),"Median"]})})
b.simicf.kcf.p = sapply(1:4,function(c){sapply(1:ns,function(i){summary.rgr.simkcf.p$summaries[paste0("b.p[",i,",",c,"]"),"Median"]})})

dat.figs4n = data.frame(slp=array(t(s.simicf.kcf.p)), bis=array(t(b.simicf.kcf.p)), cf=factor(rep(1:4,ns)), s=factor(array(t(array(rep(1:ns,4),dim=c(ns,4))))))

figs4nl = ggplot(data=dat.figs4n, aes(x=cf, y=slp, fill=cf)) +
        stat_boxplot(geom='errorbar',size=0.5, width=0.25) +
        geom_boxplot(aes(fill=cf), outlier.shape=3, outlier.size=1.2) +
        scale_fill_manual(values=colors[5:8]) +
        scale_x_discrete(name="Deck confidence\nlevel") +
        scale_y_continuous(limits=c(-0,0.4), breaks=seq(0,0.4,0.1), name="Slope") + 
        theme_custom
figs4nl

figs4nr = ggplot(data=dat.figs4n, aes(x=cf, y=bis, fill=cf)) +
        stat_boxplot(geom='errorbar',size=0.5, width=0.25) +
        geom_boxplot(aes(fill=cf), outlier.shape=3, outlier.size=1.2) +
        scale_fill_manual(values=colors[5:8]) +
        scale_x_discrete(name="Deck confidence\nlevel") +
        scale_y_continuous(limits=c(1,4), labels=scales::number_format(accuracy=0.1), name="Bias") + 
        theme_custom
figs4nr