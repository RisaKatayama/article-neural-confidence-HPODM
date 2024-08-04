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
load(paste0(c(dat.folder,"bhvdata_mri.RData"),collapse="/"))

load(paste0(c(mdlres.folder,"cwHDM_dimdl.RData"),collapse="/"))

inits1 = dump.format(list(.RNG.name="base::Super-Duper", .RNG.seed=111))
inits2 = dump.format(list(.RNG.name="base::Wichmann-Hill", .RNG.seed=222))
inits3 = dump.format(list(.RNG.name="base::Mersenne-Twister", .RNG.seed=333))

colors = c("#ff9999", "#cc0033", "#9999ff", "#3300cc", "#009900", "#99ee00", "#ffbb00", "#ff6600", 
           "#7ad151", "#aadd32", "#336600", "#27ad81")

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
mdata = Mridata[!is.element(Mridata$sidx,exsbj),]
mdata$sidx.org = mdata$sidx
for (i in 1:ns) {
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
acc.dkpred = rep(NaN,ns)
acc.kcfpred = rep(NaN,ns)
acc.dipred = rep(NaN,ns)
acc.icfpred = rep(NaN,ns)

corr.dk = (data$dkest==apply(cbind(res$dk.res$fitres$pmat.dk,data$dkest), 1, calc.pred2)) / apply(res$dk.res$fitres$pmat.dk, 1, calc.nmap)
corr.kcf = (data$kcf==apply(cbind(res$dk.res$fitres$pmat.kcf,data$kcf), 1, calc.pred4)) / apply(res$dk.res$fitres$pmat.kcf, 1, calc.nmap)
corr.di = (data$deci==apply(cbind(res$di.res$fitres$pmat.di,data$deci), 1, calc.pred2)) / apply(res$di.res$fitres$pmat.di, 1, calc.nmap)
corr.icf = (data$icf==apply(cbind(res$di.res$fitres$pmat.icf,data$icf), 1, calc.pred4)) / apply(res$di.res$fitres$pmat.icf, 1, calc.nmap)

for (i in 1:ns){
    acc.dkpred[i] = mean(corr.dk[data$sidx==i])
    acc.kcfpred[i] = mean(corr.kcf[data$sidx==i])
    acc.dipred[i] = mean(corr.di[data$sidx==i])
    acc.icfpred[i] = mean(corr.icf[data$sidx==i])
}


## Figure 3bc
sim.dkest = apply(cbind(res$dk.res$fitres$pmat.dk,data$dkest), 1, calc.pred2)
sim.dkest[apply(res$dk.res$fitres$pmat.dk, 1, calc.nmap)==2] = 0.5

pdk6.seq = array(NaN,dim=c(ng,16))
simpdk6.seq = array(NaN,dim=c(ng,16))
kcf.seq = array(NaN,dim=c(ng,16))
simkcf.seq = array(NaN,dim=c(ng,16))
for (p in 1:ng) {
    for (t in 1:16) {
        pdk6.seq[p,t] = mean(data$dkest[data$seqid==p&data$trl==t])
        simpdk6.seq[p,t] = mean(sim.dkest[data$seqid==p&data$trl==t])
        kcf.seq[p,t] = mean(data$kcf[data$seqid==p&data$trl==t])
        simkcf.seq[p,t] = mean(res$dk.res$fitres$w.pred.kcf[data$seqid==p&data$trl==t])
    }
}

exseq = c(20,7) # example sequences
dat.fig3b = data.frame(prp=c(pdk6.seq[exseq[1],],NaN,pdk6.seq[exseq[2],],simpdk6.seq[exseq[1],],NaN,simpdk6.seq[exseq[2],]),
                       add=c(data$add[data$sidx==i&data$seqid==exseq[1]],NaN,data$add[data$sidx==i&data$seqid==exseq[2]]),
                       dktrue=factor(c(data$dktrue[data$sidx==i&data$seqid==exseq[1]],0,data$dktrue[data$sidx==i&data$seqid==exseq[2]])),
                       trl=rep(c(1:16,16.5,17:32),2), typ=factor(c(rep(0,33),rep(1,33))))

fig3b = ggplot(data=dat.fig3b, aes(x=trl, y=prp, color=typ, linetype=typ)) +
        geom_line(linewidth=0.75, alpha=1) +
        geom_vline(xintercept=16.5, linetype="dotted", color="dimgray", linewidth=0.5) +
        scale_color_manual(name=NULL, values=colors[c(9,10)],labels=c(" Data "," Model ")) +
        geom_text(data=dat.fig3b[dat.fig3b$dktrue==0,], aes(x=trl, label=add), y=-0.06, color="red", size=2.5, show.legend=FALSE) +
        geom_text(data=dat.fig3b[dat.fig3b$dktrue==1,], aes(x=trl, label=add), y=1.06, color="blue", size=2.5, show.legend=FALSE) +
        scale_linetype_manual(name=NULL, values=c('solid','dotdash'), labels=c(" Data "," Model "), 
                               guide=guide_legend(override.aes=list(color=colors[c(9,10)], linetype=c('solid','dotdash'), linewidth=c(0.75,0.75)))) +
        scale_x_continuous(breaks=1:32, labels=rep(c(1,"",3,"",5,"",7,"",9,"",11,"",13,"",15,""),2), expand=c(0.02,0.02),
                           name="Trial index within a single game") +
        scale_y_continuous(limits=c(-0.05,1.05), breaks=seq(0,1,0.2), name="Proportion of choosing\n\"deck 6\"") +
        theme_custom +
        theme(legend.position=c(1,0.05),legend.justification=c(1,0),legend.background=element_blank(),
              legend.box.background=element_rect(colour="gray"),legend.title=element_text(size=10),
              legend.text=element_text(size=10, margin=margin(r=0,t=0.1,unit='cm')),legend.spacing.x=unit(0,'cm'),legend.spacing.y=unit(0,'cm'),
              legend.key.width=unit(0.5,"cm"), legend.key.height=unit(0.3,"cm"))
fig3b

dat.fig3c = data.frame(kcf=c(kcf.seq[exseq[1],],NaN,kcf.seq[exseq[2],],simkcf.seq[exseq[1],],NaN,simkcf.seq[exseq[2],]),
                       add=c(data$add[data$sidx==i&data$seqid==exseq[1]],NaN,data$add[data$sidx==i&data$seqid==exseq[2]]),
                       dktrue=factor(c(data$dktrue[data$sidx==i&data$seqid==exseq[1]],0,data$dktrue[data$sidx==i&data$seqid==exseq[2]])),
                       trl=rep(c(1:16,16.5,17:32),2), typ=factor(c(rep(0,33),rep(1,33))))

fig3c = ggplot(data=dat.fig3c, aes(x=trl, y=kcf, color=typ, linetype=typ)) +
        geom_line(linewidth=0.75, alpha=1) +
        geom_vline(xintercept=16.5, linetype="dotted", color="dimgray", linewidth=0.5) +
        scale_color_manual(name=NULL, values=colors[c(11,12)],labels=c(" Data "," Model ")) +
        geom_text(data=dat.fig3c[dat.fig3c$dktrue==0,], aes(x=trl, label=add), y=0.96, color="red", size=2.5, show.legend=FALSE) +
        geom_text(data=dat.fig3c[dat.fig3c$dktrue==1,], aes(x=trl, label=add), y=4.06, color="blue", size=2.5, show.legend=FALSE) +
        scale_linetype_manual(name=NULL, values=c('solid','dotdash'), labels=c(" Data "," Model "), 
                               guide=guide_legend(override.aes=list(color=colors[c(11,12)], linetype=c('solid','dotdash'), linewidth=c(0.75,0.75)))) +
        scale_x_continuous(breaks=1:32, labels=rep(c(1,"",3,"",5,"",7,"",9,"",11,"",13,"",15,""),2), expand=c(0.02,0.02),
                           name="Trial index within a single game") +
        scale_y_continuous(limits=c(0.95,4.05), breaks=seq(1,4,1), name="Averaged deck\nconfidence level", labels=scales::number_format(accuracy=0.1)) +
        theme_custom +
        theme(legend.position=c(1,0.05),legend.justification=c(1,0),legend.background=element_blank(),
              legend.box.background=element_rect(colour="gray"),legend.title=element_text(size=10),
              legend.text=element_text(size=10, margin=margin(r=0,t=0.1,unit='cm')),legend.spacing.x=unit(0,'cm'),legend.spacing.y=unit(0,'cm'),
              legend.key.width=unit(0.5,"cm"), legend.key.height=unit(0.3,"cm"))
fig3c


## Fig 3de
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

sim.phit.icfdk.cwh = array(NaN,dim=c(ns,11,2,2))
for (i in 1:ns) {
    for (d in 0:1) {
        for (o in 8:18) {
            tmpicf = res$di.res$fitres$w.pred.icf[data$sidx==i]
            sim.phit.icfdk.cwh[i,o-7,1,d+1] = mean(res$di.res$fitres$pmat.di[data$sidx==i&data$dkest==d&data$open==o&res$di.res$fitres$w.pred.icf<=mean(tmpicf),2])
            sim.phit.icfdk.cwh[i,o-7,2,d+1] = mean(res$di.res$fitres$pmat.di[data$sidx==i&data$dkest==d&data$open==o&res$di.res$fitres$w.pred.icf>mean(tmpicf),2])
        }
    }
}

dat.fig3d = data.frame(prp=array(apply(phit.icfdk,c(2,3,4),nanmean)),
                       s=array(apply(phit.icfdk,c(2,3,4),nansd)/sqrt(apply(!is.na(phit.icfdk),c(2,3,4),nansum))),
                       opn=rep(8:18,2*2), cf=factor(rep(c(rep(0,11),rep(1,11)),2)), dk=factor(c(rep(0,11*2),rep(1,11*2))))

simdat.fig3d = data.frame(prp=array(apply(sim.phit.icfdk.cwh,c(2,3,4),nanmean)),
                          s=array(apply(sim.phit.icfdk.cwh,c(2,3,4),nansd)/sqrt(apply(!is.na(sim.phit.icfdk.cwh),c(2,3,4),nansum))),
                          opn=rep(8:18,2*2), cf=factor(rep(c(rep(0,11),rep(1,11)),2)), dk=factor(c(rep(0,11*2),rep(1,11*2))))

fig3d = ggplot(data=dat.fig3d, aes(x=opn, y=prp, group=interaction(cf,dk), fill=interaction(cf,dk), color=interaction(cf,dk))) +
        geom_errorbar(aes(x=opn, ymin=prp-s, ymax=prp+s, group=interaction(cf,dk), color=interaction(cf,dk)), width=0.15, alpha=0.5) +
        geom_point(aes(fill=interaction(cf,dk)), shape=21, color="white", size=2.5, alpha=0.75) +
        geom_ribbon(data=simdat.fig3d, aes(x=opn, ymin=prp-s, ymax=prp+s, group=interaction(cf,dk), fill=interaction(cf,dk)), colour=NA, alpha=0.3, show.legend=FALSE) +
        geom_line(data=simdat.fig3d, aes(x=opn, y=prp, group=interaction(cf,dk), color=interaction(cf,dk)), linewidth=0.7, alpha=0.9) +
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
fig3d

predat = data.frame(p.hit=array(aperm(sim.phit.icfdk.cwh,c(2,3,4,1))), open=rep(8:18,ns*2*2), cf=rep(c(rep(0,11),rep(1,11)),ns*2),
                 dk=rep(c(rep(0,11*2),rep(1,11*2)),ns), sidx=array(t(array(rep(1:ns,11*2*2),dim=c(ns,11*2*2)))))
predat = predat[!is.na(predat$p.hit),]

dat = dump.format(list(open=predat$open, cf=predat$cf, dk=predat$dk, p.hit=predat$p.hit, sidx=predat$sidx, ns=ns, N=dim(predat)[1]))
vals = c("a.mu","b.mu","a.sig","b.sig","a.p","b.p")

rgr.sim.icfdk = run.jags(model=paste0(mdl.folder,"SimpleRegress_pseudoLogistic.txt"), monitor=vals, data=dat, n.chains=3, inits=c(inits1,inits2,inits3),
                        plots=FALSE, method="parallel", burnin=5000, sample=1000, adapt=1000, thin=10)
summary.rgr.sim.icfdk = add.summary(rgr.sim.icfdk)

s.sim.icfdk = array(sapply(1:4,function(c){sapply(1:ns,function(i){summary.rgr.sim.icfdk$summaries[paste0("a.p[",i,",",c,"]"),"Median"]})}),dim=c(ns,2,2))
b.sim.icfdk = array(sapply(1:4,function(c){sapply(1:ns,function(i){summary.rgr.sim.icfdk$summaries[paste0("b.p[",i,",",c,"]"),"Median"]})}),dim=c(ns,2,2))

dat.fig3e = data.frame(slp=array(aperm(s.sim.icfdk,c(2,3,1))), sft=array(aperm(b.sim.icfdk,c(2,3,1))), cf=factor(rep(c('Low','High','Low','High'),ns),levels=c('Low','High')),
                       dk=factor(rep(c('Deck 4','Deck 4','Deck 6','Deck 6'),ns)), s=factor(array(t(array(rep(1:ns,2*2),dim=c(ns,2*2))))))

fig3el = ggplot(data=dat.fig3e, aes(x=interaction(cf,dk), y=slp, fill=interaction(cf,dk))) +
        stat_boxplot(geom='errorbar',size=0.5, width=0.25) +
        geom_boxplot(aes(fill=interaction(cf,dk)), outlier.shape=3, outlier.size=1.2) +
        scale_fill_manual(values=colors[1:4]) +
        scale_x_discrete(guide="axis_nested") +
        scale_y_continuous(limits=c(-5,0), breaks=-5:1, labels=scales::number_format(accuracy=0.1), name="Slope") + 
        geom_signif(comparisons=list(c(1,2),c(3,4)), annotations=c('***'), y_position=-0.3, tip_length=0, vjust=0.2, textsize=4.5) +
        geom_signif(comparisons=list(c(2,4)), annotations=c('**'), y_position=c(-4.2), tip_length=0, vjust=1.8, textsize=4.5) +
        theme_custom +
        theme(axis.title.x=element_blank(), axis.text.x=element_text(angle=30, hjust=1, vjust=1, margin=margin(t=0.05, b=0, unit="cm")),
              ggh4x.axis.nestline.x=element_line(linewidth=0.25), ggh4x.axis.nesttext.x=element_text(angle=30, hjust=0.5, vjust=0.5, margin=margin(t=0.09, unit="cm")))
fig3el

fig3er = ggplot(data=dat.fig3e, aes(x=interaction(cf,dk), y=sft, fill=interaction(cf,dk))) +
        stat_boxplot(geom='errorbar',size=0.5, width=0.25) +
        geom_boxplot(aes(fill=interaction(cf,dk)), outlier.shape=3, outlier.size=1.2) +
        scale_fill_manual(values=colors[1:4]) +
        scale_x_discrete(guide="axis_nested") +
        scale_y_continuous(limits=c(9,15), breaks=seq(9,15,1.5), name="Shift") + 
        geom_signif(comparisons=list(c(1,2)), annotations=c('*'), y_position=14.5, tip_length=0, vjust=-0.2, textsize=4.5) +
        geom_signif(comparisons=list(c(1,3),c(2,4)), annotations=c('***'), y_position=c(10.5,9.7), tip_length=0, vjust=1.8, textsize=4.5) +
        theme_custom +
        theme(axis.title.x=element_blank(), axis.text.x=element_text(angle=30, hjust=1, vjust=1, margin=margin(t=0.05, b=0, unit="cm")),
              ggh4x.axis.nestline.x=element_line(linewidth=0.25), ggh4x.axis.nesttext.x=element_text(angle=30, hjust=0.5, vjust=0.5, margin=margin(t=0.09, unit="cm")))
fig3er


## Fig 3fg
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

sim.phit.kcfdk.cwh = array(NaN,dim=c(ns,11,2,2))
for (i in 1:ns){
    tmpkcf = res$dk.res$fitres$w.pred.kcf[data$sidx==i]
    qt = mean(tmpkcf)
    for (d in 0:1) {
        for (o in 8:18) {
            sim.phit.kcfdk.cwh[i,o-7,1,d+1] = mean(res$di.res$fitres$pmat.di[data$sidx==i&data$dkest==d&data$open==o&res$dk.res$fitres$w.pred.kcf<=qt,2])
            sim.phit.kcfdk.cwh[i,o-7,2,d+1] = mean(res$di.res$fitres$pmat.di[data$sidx==i&data$dkest==d&data$open==o&res$dk.res$fitres$w.pred.kcf>qt,2])
        }
    }
}

dat.fig3f = data.frame(prp=array(apply(phit.kcfdk,c(2,3,4),nanmean)),
                       s=array(apply(phit.kcfdk,c(2,3,4),nansd)/sqrt(apply(!is.na(phit.kcfdk),c(2,3,4),nansum))),
                       opn=rep(8:18,2*2), cf=factor(rep(c(rep(0,11),rep(1,11)),2)), dk=factor(c(rep(0,11*2),rep(1,11*2))))

simdat.fig3f = data.frame(prp=array(apply(sim.phit.kcfdk.cwh,c(2,3,4),nanmean)),
                          s=array(apply(sim.phit.kcfdk.cwh,c(2,3,4),nansd)/sqrt(apply(!is.na(sim.phit.kcfdk.cwh),c(2,3,4),nansum))),
                          opn=rep(8:18,2*2), cf=factor(rep(c(rep(0,11),rep(1,11)),2)), dk=factor(c(rep(0,11*2),rep(1,11*2))))

fig3f = ggplot(data=dat.fig3f, aes(x=opn, y=prp, group=interaction(cf,dk), fill=interaction(cf,dk), color=interaction(cf,dk))) +
        geom_errorbar(aes(x=opn, ymin=prp-s, ymax=prp+s, group=interaction(cf,dk), color=interaction(cf,dk)), width=0.15, alpha=0.5) +
        geom_point(aes(fill=interaction(cf,dk)), shape=21, color="white", size=2.5, alpha=0.75) +
        geom_ribbon(data=simdat.fig3f, aes(x=opn, ymin=prp-s, ymax=prp+s, group=interaction(cf,dk), fill=interaction(cf,dk)), colour=NA, alpha=0.3, show.legend=FALSE) +
        geom_line(data=simdat.fig3f, aes(x=opn, y=prp, group=interaction(cf,dk), color=interaction(cf,dk)), linewidth=0.7, alpha=0.9) +
        scale_fill_manual(name="Deck conf.\n\n", values=colors[1:4], labels=c("",""," Deck 4 "," Deck 6"),
                          guide=guide_legend(override.aes=list(color=colors[c(1,3,2,4)], shape=c(15,15,15,15), size=4.5, alpha=1),ncol=2)) +
        scale_color_manual(name=NULL, values=colors[1:4], guide="none") +
        scale_x_continuous(breaks=seq(8,18,2), name="Face-up card score") +
        scale_y_continuous(limits=c(-0.01,1.01), breaks=seq(0,1,0.2), name="Proportion of Hit") +
        theme_custom +
        theme(legend.position=c(1,1),legend.justification=c(1,1),legend.background=element_rect(colour="transparent"),
              legend.title=element_text(size=10),
              legend.text=element_text(size=10, margin=margin(r=0,unit='cm')),legend.spacing.x=unit(-0.02,'cm'),legend.spacing.y=unit(0,'cm'),
              legend.key.width=unit(0.5,"cm"), legend.key.height=unit(0.3,"cm"))
fig3f

predat = data.frame(p.hit=array(aperm(sim.phit.kcfdk.cwh,c(2,3,4,1))), open=rep(8:18,ns*2*2), cf=rep(c(rep(0,11),rep(1,11)),ns*2),
                 dk=rep(c(rep(0,11*2),rep(1,11*2)),ns), sidx=array(t(array(rep(1:ns,11*2*2),dim=c(ns,11*2*2)))))
predat = predat[!is.na(predat$p.hit),]

dat = dump.format(list(open=predat$open, cf=predat$cf, dk=predat$dk, p.hit=predat$p.hit, sidx=predat$sidx, ns=ns, N=dim(predat)[1]))
vals = c("a.mu","b.mu","a.sig","b.sig","a.p","b.p")

rgr.sim.kcfdk = run.jags(model=paste0(mdl.folder,"SimpleRegress_pseudoLogistic.txt"), monitor=vals, data=dat, n.chains=3, inits=c(inits1,inits2,inits3),
                        plots=FALSE, method="parallel", burnin=5000, sample=1000, adapt=1000, thin=10)
summary.rgr.sim.kcfdk = add.summary(rgr.sim.kcfdk)

s.sim.kcfdk = array(sapply(1:4,function(c){sapply(1:ns,function(i){summary.rgr.sim.kcfdk$summaries[paste0("a.p[",i,",",c,"]"),"Median"]})}),dim=c(ns,2,2))
b.sim.kcfdk = array(sapply(1:4,function(c){sapply(1:ns,function(i){summary.rgr.sim.kcfdk$summaries[paste0("b.p[",i,",",c,"]"),"Median"]})}),dim=c(ns,2,2))

dat.fig3g = data.frame(slp=array(aperm(s.sim.kkcfdk,c(2,3,1))), sft=array(aperm(b.sim.kcfdk,c(2,3,1))), cf=factor(rep(c('Low','High','Low','High'),ns),levels=c('Low','High')),
                       dk=factor(rep(c('Deck 4','Deck 4','Deck 6','Deck 6'),ns)), s=factor(array(t(array(rep(1:ns,2*2),dim=c(ns,2*2))))))

art.dat.fig3g.slp = art(slp~cf*dk+Error(s+s:cf+s:dk+s:cf:dk), data=dat.fig3g)

fig3gl = ggplot(data=dat.fig3g, aes(x=interaction(cf,dk), y=slp, fill=interaction(cf,dk))) +
        stat_boxplot(geom='errorbar',linewidth=0.5, width=0.25) +
        geom_boxplot(aes(fill=interaction(cf,dk)), outlier.shape=3, outlier.size=1.2) +
        scale_fill_manual(values=colors[1:4]) +
        scale_x_discrete(guide="axis_nested") +
        scale_y_continuous(limits=c(-2.5,0.03), breaks=seq(-2.5,0,0.5), labels=scales::number_format(accuracy=0.1), name="Slope") + 
        theme_custom +
        theme(axis.title.x=element_blank(), axis.text.x=element_text(angle=30, hjust=1, vjust=1, margin=margin(t=0.05, b=0, unit="cm")),
              ggh4x.axis.nestline.x=element_line(linewidth=0.25), ggh4x.axis.nesttext.x=element_text(angle=30, hjust=0.5, vjust=0.5, margin=margin(t=0.09, unit="cm")))
fig3gl

fig3gr = ggplot(data=dat.fig3g, aes(x=interaction(cf,dk), y=bis, fill=interaction(cf,dk))) +
        stat_boxplot(geom='errorbar',size=0.5, width=0.25) +
        geom_boxplot(aes(fill=interaction(cf,dk)), outlier.shape=3, outlier.size=1.2) +
        scale_fill_manual(values=colors[1:4]) +
        scale_x_discrete(guide="axis_nested") +
        scale_y_continuous(limits=c(9,15), breaks=seq(9,15,1.5), name="Shift") + 
        geom_signif(comparisons=list(c(1,2),c(3,4)), annotations=c('***','***'), y_position=14.7, tip_length=0, vjust=0.2, textsize=4.5) +
        geom_signif(comparisons=list(c(1,3),c(2,4)), annotations=c('***','***'), y_position=c(10.5,9.7), tip_length=0, vjust=1.8, textsize=4.5) +
        theme_custom +
        theme(axis.title.x=element_blank(), axis.text.x=element_text(angle=30, hjust=1, vjust=1, margin=margin(t=0.05, b=0, unit="cm")),
              ggh4x.axis.nestline.x=element_line(linewidth=0.25), ggh4x.axis.nesttext.x=element_text(angle=30, hjust=0.5, vjust=0.5, margin=margin(t=0.09, unit="cm")))
fig3gr


## Fig 3hi
icf.kcf = array(NaN, dim=c(ns,11,4))
for (i in 1:ns) {
    for (o in 8:18) {
        for (c in 1:4) {
            icf.kcf[i,o-7,c] = mean(data$icf[data$sidx==i&data$open==o&data$kcf==c])
        }
    }
}

sim.icf.kcf.cwh = array(NaN,dim=c(ns,11,4))
for (i in 1:ns){
    for (c in 1:4) {
        for (o in 8:18) {
           sim.icf.kcf.cwh[i,o-7,c] = mean(res$di.res$fitres$w.pred.icf[data$sidx==i&data$open==o&res$dk.res$fitres$pred.kcf==c])
        }
    }
}

dat.fig3h = data.frame(icf=array(apply(icf.kcf,c(2,3),nanmean)),
                       s=array(apply(icf.kcf,c(2,3),nansd)/sqrt(apply(!is.na(icf.kcf),c(2,3),nansum))),
                       opn=rep(8:18,4), cf=factor(c(rep(1,11),rep(2,11),rep(3,11),rep(4,11))))

simdat.fig3h = data.frame(icf=array(apply(sim.icf.kcf.cwh,c(2,3),nanmean)),
                          s=array(apply(sim.icf.kcf.cwh,c(2,3),nansd)/sqrt(apply(!is.na(sim.icf.kcf.cwh),c(2,3),nansum))),
                          opn=rep(8:18,4), cf=factor(c(rep(1,11),rep(2,11),rep(3,11),rep(4,11))))

fig3h = ggplot(data=dat.fig3h, aes(x=opn, y=icf, group=cf, fill=cf, color=cf)) +
        geom_errorbar(aes(x=opn, ymin=icf-s, ymax=icf+s, group=cf, color=cf), width=0.15, alpha=0.5) +
        geom_point(aes(fill=cf), shape=21, color="white", size=2.5, alpha=0.75) +
        geom_ribbon(data=simdat.fig3h, aes(x=opn, ymin=icf-s, ymax=icf+s, group=cf, fill=cf), colour=NA, alpha=0.3, show.legend=FALSE) +
        geom_line(data=simdat.fig3h, aes(x=opn, y=icf, group=cf, color=cf), linewidth=0.7, alpha=0.9) +
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
fig3h

sim.phit.kcf.cwh = array(NaN,dim=c(ns,11,4))
for (i in 1:ns){
    for (o in 8:18) {
        for (c in 1:4) {
            sim.phit.kcf.cwh[i,o-7,c] = mean(res$di.res$fitres$pmat.di[data$sidx==i&data$open==o&res$dk.res$fitres$pred.kcf==c,2])
        }
    }
}

predat = data.frame(p.hit=array(aperm(sim.phit.kcf.cwh,c(2,3,1))), icf=array(aperm(sim.icf.kcf.cwh,c(2,3,1))),
                    open=rep(8:18,ns*4), kcf=rep(c(rep(1,11),rep(2,11),rep(3,11),rep(4,11)),ns),
                    sidx=array(t(array(rep(1:ns,11*2*2),dim=c(ns,11*2*2)))))
predat = predat[!is.na(predat$p.hit)&!is.na(predat$icf),]
predat$sv = abs(logit(predat$p.hit))

dat = dump.format(list(sv=predat$sv, kcf=predat$kcf, icf=predat$icf, sidx=predat$sidx, ns=ns, N=dim(predat)[1]))
vals = c("a.mu","b.mu","a.sig","b.sig","a.p","b.p","sig.mu")

rgr.sim.kcf = run.jags(model="SimpleRegress.txt", monitor=vals, data=dat, n.chains=3, inits=c(inits1,inits2,inits3),
                     plots=FALSE, method="parallel", burnin=5000, sample=1000, adapt=1000, thin=10)
summary.rgr.sim.kcf = add.summary(rgr.sim.kcf)

s.sim.icf.kcf = sapply(1:4,function(c){sapply(1:ns,function(i){summary.rgr.sim.kcf$summaries[paste0("a.p[",i,",",c,"]"),"Median"]})})
b.sim.icf.kcf = sapply(1:4,function(c){sapply(1:ns,function(i){summary.rgr.sim.kcf$summaries[paste0("b.p[",i,",",c,"]"),"Median"]})})

dat.fig3i = data.frame(slp=array(t(s.sim.icf.kcf)), bis=array(t(b.sim.icf.kcf)), cf=factor(rep(1:4,ns)), s=factor(array(t(array(rep(1:ns,4),dim=c(ns,4))))))

cor.test(as.numeric(dat.fig3i$cf),dat.fig3i$slp)
cor.test(as.numeric(dat.fig3i$cf),dat.fig3i$bis)

fig3il = ggplot(data=dat.fig3i, aes(x=cf, y=slp, fill=cf)) +
        stat_boxplot(geom='errorbar',size=0.5, width=0.25) +
        geom_boxplot(aes(fill=cf), outlier.shape=3, outlier.size=1.2) +
        scale_fill_manual(values=colors[5:8]) +
        scale_x_discrete(name="Deck confidence\nlevel") +
        scale_y_continuous(limits=c(-0,0.4), breaks=seq(0,0.4,0.1), name="Slope") + 
        theme_custom
fig3il

fig3ir = ggplot(data=dat.fig3i, aes(x=cf, y=bis, fill=cf)) +
        stat_boxplot(geom='errorbar',size=0.5, width=0.25) +
        geom_boxplot(aes(fill=cf), outlier.shape=3, outlier.size=1.2) +
        scale_fill_manual(values=colors[5:8]) +
        scale_x_discrete(name="Deck confidence\nlevel") +
        scale_y_continuous(limits=c(1,4), labels=scales::number_format(accuracy=0.1), name="Bias") + 
        theme_custom
fig3ir