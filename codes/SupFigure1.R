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
mdl.folder = "model/"
load(paste0(c(dat.folder,"bhvdata_mri.RData"),collapse="/"))

inits1 = dump.format(list(.RNG.name="base::Super-Duper", .RNG.seed=111))
inits2 = dump.format(list(.RNG.name="base::Wichmann-Hill", .RNG.seed=222))
inits3 = dump.format(list(.RNG.name="base::Mersenne-Twister", .RNG.seed=333))

colors = c("#ff9999", "#cc0033", "#9999ff", "#3300cc", "#009900", "#99ee00", "#ffbb00", "#ff6600")

exsbj = c(7,22) # subjects with low deck inference accuracy
ns.org = max(Mridata$sidx)
insbj = c(1:ns.org)
insbj = insbj[!is.element(insbj,exsbj)]
mdata = Mridata[!is.element(Mridata$sidx,exsbj),]
mdata$sidx.org = mdata$sidx
ns = length(insbj)
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

    undisc = data$undisc # actual score of the face-down card
    rw = data$rw

    trl = data$trl
    ses = data$ses
    sidx = data$sidx
    
    include = !is.na(data$deckest*data$cf.deck*data$decision*data$cf.deci)
    res = data.frame(dkest=dkest[include], dktrue=dktrue[include], kcf=kcf[include], kcf.hl=kcf.hl[include], add=add[include],
               deci=deci[include], open=open[include], undisc=undisc[include], icf=icf[include], icf.hl=icf.hl[include], 
               rw=rw[include], trl=trl[include], ses=ses[include], sidx=sidx[include])
    return(res)
}
data = load.data.basic(mdata)
ns = max(data$sidx)
ng = max(data$ses)
N = dim(data)[1]


#####################
## Fig S1a
dkacc.trl = array(NaN,dim=c(ns,16))
for (i in 1:ns) {
    for (t in 1:16) {
        dkacc.trl[i,t] = sum(data$dkest==data$dktrue&data$trl==t&data$sidx==i)/sum(data$trl==t&data$sidx==i)
    }
}

dat.figs1a = data.frame(acc=array(t(dkacc.trl)), trl=factor(rep(1:16,ns)))
dat.figs1a = dat.figs1a[!is.na(dat.figs1a$acc),]

figs2a = ggplot(data=dat.figs1a, aes(x=trl, y=acc)) +
        stat_boxplot(geom='errorbar', linewidth=0.6, width=0.3) +
        geom_boxplot(outlier.shape=3, outlier.size=1.2, fill="gray") +
        geom_hline(yintercept=0.5, linetype="dotted") +
        scale_x_discrete(name="Trial index in a single game") +
        scale_y_continuous(limits=c(0,1), name="Deck inference accuracy") + 
        theme_custom
figs2a


## Fig S1b
dkacc.kcf = array(NaN,dim=c(ns,4))
for (i in 1:ns) {
    for (c in 1:4) {
        dkacc.kcf[i,c] = sum(data$dkest==data$dktrue&data$kcf==c&data$sidx==i)/sum(data$kcf==c&data$sidx==i)
    }
}

dat.figs1b = data.frame(acc=array(t(dkacc.kcf)), cf=factor(rep(1:4,ns)))
dat.figs1b = dat.figs1b[!is.na(dat.figs1b$acc),]

fig2b = ggplot(data=dat.figs1b, aes(x=cf, y=acc)) +
        stat_boxplot(geom='errorbar', linewidth=0.6, width=0.3) +
        geom_boxplot(outlier.shape=3, outlier.size=1.2, fill="gray") +
        geom_hline(yintercept=0.5, linetype="dotted") +
        geom_text(data=data.frame(x=2:4, acc=rep(1.02,3), label=c('***','***','***')), aes(x=x, y=acc, label=label), size=4) +
        scale_x_discrete(name="Deck confidence\nlevel") +
        scale_y_continuous(limits=c(0,1.02), name="Deck inference accuracy") + 
        theme_custom
fig2b

cor.test(as.numeric(dat.figs1b$cf), dat.figs1b$acc)


## Fig S1c
rw.chosen = data$rw  # reward based on chosen action
rw.unchosen = data$open + data$undisc + data$add*(1-data$deci)
rw.unchosen[rw.unchosen>21] = 0  # (potential) reward based on unchosen action
difrw = rw.chosen - rw.unchosen
difrw.icf = array(NaN, dim=c(ns,4))
for (i in 1:ns) {
    for (c in 1:4) {
        difrw.icf[i,c] = mean(difrw[data$icf==c&data$sidx==i])
    }
}

dat.figs1c = data.frame(dg=array(t(difrw.icf)), cf=factor(rep(1:4,2*ns)))
dat.figs1c = dat.figs1c[!is.na(dat.figs1c$dg),]

fig2c = ggplot(data=dat.figs1c, aes(x=cf, y=dg)) +
        stat_boxplot(geom='errorbar', linewidth=0.6, width=0.3) +
        geom_boxplot(outlier.shape=3, outlier.size=1.2, fill="gray") +
        geom_hline(yintercept=0, linetype="dotted") +
        geom_text(data=data.frame(cf=1:4, dg=rep(15,4), label=c('**','***','***','***')), aes(x=cf, y=dg, label=label), size=4) +
        scale_x_discrete(name="Decision\nconfidence level") +
        scale_y_continuous(limits=c(-5,15), name="Relative gain of chosen action") + 
        theme_custom
fig2c

cor.test(as.numeric(dat.figs1c$cf), dat.figs1c$dg)


## Fig S1d
prp.kcf.icf =array(NaN, dim=c(4,4,ns))
for (i in 1:ns) {
    for (ck in 1:4) {
        for (ci in 1:4) {
            prp.kcf.icf[ck,ci,i] = sum(data$kcf==ck&data$icf==ci&data$sidx==i)/sum(data$sidx==i)
        }
    }
}

dat.figs1d = data.frame(ck=factor(rep(1:4,4)), ci=factor(array(t(array(rep(1:4,4),dim=c(4,4))))), prp=array(apply(prp.kcf.icf,c(1,2),nanmean)))

fig2d = ggplot(data=dat.figs1d, aes(x=ck, y=ci, fill=prp)) +
        geom_tile(color='black') +
        scale_fill_gradient(low='#ffffff', high="#333333", limit = c(0,0.15), breaks=seq(0,0.15,0.05), name="Porportion",
                            guide=guide_colourbar(title=NULL, frame.colour="black", frame.linewidth=0.8/.pt, ticks.colour='black')) +
        scale_x_discrete(expand=expansion(c(0.25,0.25)), name="Deck confidence level") +
        scale_y_discrete(expand=expansion(c(0.25,0.25)), name="Decision confidence level") +
        theme_custom +
        theme(legend.position="right", legend.title=element_text(size=10, angle=90), legend.title.align=0.5, legend.direction="vertical") +
        guides(fill=guide_colorbar(barwidth=1, barheight=7, title.position="left", frame.colour="black", ticks.colour="black"))
fig2d


## Fig S1ef
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

dat = dump.format(list(open=data$open, cf=data$icf.hl, dk=data$dkest, di=data$deci, sidx=data$sidx, ns=ns, N=N))
vals = c("a.mu","b.mu","a.sig","b.sig","a.p","b.p")

rgr.icfdk = run.jags(model=paste0(mdl.folder,"SimpleLogisticRegress.txt"), monitor=vals, data=dat, n.chains=3, inits=c(inits1,inits2,inits3),
                     plots=FALSE, method="parallel", burnin=5000, sample=1000, adapt=1000, thin=10)
summary.rgr.icfdk = add.summary(rgr.icfdk)

s.icfdk = array(sapply(1:4,function(c){sapply(1:ns,function(i){summary.rgr.icfdk$summaries[paste0("a.p[",i,",",c,"]"),"Median"]})}),dim=c(ns,2,2))
b.icfdk = array(sapply(1:4,function(c){sapply(1:ns,function(i){summary.rgr.icfdk$summaries[paste0("b.p[",i,",",c,"]"),"Median"]})}),dim=c(ns,2,2))

x = seq(8,18,0.1)
simphit.icfdk = array(NaN,dim=c(ns,length(x),2,2))
for (i in 1:ns) {
    for (c in 0:1) {
        for (d in 0:1) {
            simphit.icfdk[i,,c+1,d+1] = ilogit(s.icfdk[i,c+1,d+1]*(x-b.icfdk[i,c+1,d+1]))
        }
    }
}

dat.figs1e = data.frame(prp=array(aperm(phit.icfdk,c(2,3,4,1))), opn=rep(8:18,2*2*ns),
                       cf=factor(rep(c(rep(0,11),rep(1,11)),2*ns)), dk=factor(rep(c(rep(0,11*2),rep(1,11*2)),ns)))
dat.figs1e = dat.figs1e[!is.na(dat.figs1e$prp),]

sumdat.figs1e = data.frame(prp=array(apply(phit.icfdk,c(2,3,4),nanmean)),
                          s=array(apply(phit.icfdk,c(2,3,4),nansd)/sqrt(apply(!is.na(phit.icfdk),c(2,3,4),nansum))),
                          opn=rep(8:18,2*2), cf=factor(rep(c(rep(0,11),rep(1,11)),2)), dk=factor(c(rep(0,11*2),rep(1,11*2))))

simdat.figs1e = data.frame(prp=array(apply(simphit.icfdk,c(2,3,4),nanmean)),
                          s=array(apply(simphit.icfdk,c(2,3,4),nansd)/sqrt(apply(!is.na(simphit.icfdk),c(2,3,4),nansum))),
                          opn=rep(x,2*2), cf=factor(rep(c(rep(0,length(x)),rep(1,length(x))),2)), dk=factor(c(rep(0,length(x)*2),rep(1,length(x)*2))))

fig2e = ggplot(data=sumdat.figs1e, aes(x=opn, y=prp, group=interaction(cf,dk), fill=interaction(cf,dk), color=interaction(cf,dk))) +
        geom_line(data=simdat.figs1e, aes(x=opn, y=prp, group=interaction(cf,dk), color=interaction(cf,dk)), linewidth=0.6, alpha=0.7) +
        stat_boxplot(data=dat.figs1e, aes(x=opn, y=prp, group=interaction(opn,cf,dk), color=interaction(cf,dk)), geom='errorbar', linewidth=0.1, width=0.1, alpha=0.4, position=position_dodge(width=0.75)) +
        geom_boxplot(data=dat.figs1e, aes(x=opn, y=prp, group=interaction(opn,cf,dk), fill=interaction(cf,dk), color=interaction(cf,dk)), alpha=0.4, linewidth=0.15, width=0.25, outlier.shape=3, outlier.size=0.5, position=position_dodge(width=0.75)) +
        geom_point(aes(fill=interaction(cf,dk)), shape=21, position=position_dodge(0.75), color="white", size=2.5) +
        scale_fill_manual(name="Decision conf.\n\n", values=colors[1:4], labels=c("",""," Deck 4 "," Deck 6"),
                          guide=guide_legend(override.aes=list(color=colors[c(1,3,2,4)], shape=c(15,15,15,15), size=5),ncol=2)) +
        scale_color_manual(name=NULL, values=colors[1:4], guide="none") +
        scale_x_continuous(breaks=seq(8,18,2), name="Face-up card score") +
        scale_y_continuous(limits=c(0,1), breaks=seq(0,1,0.2), name="Proportion of Hit") +
        theme_custom +
        theme(legend.position=c(1,1),legend.justification=c(1,1),legend.background=element_rect(colour="transparent"),
              legend.title=element_text(size=10),
              legend.text=element_text(size=10, margin=margin(r=0,unit='cm')),legend.spacing.x=unit(-0.02,'cm'),legend.spacing.y=unit(0,'cm'),
              legend.key.width=unit(0.5,"cm"), legend.key.height=unit(0.3,"cm"))
fig2e

dat.figs1f = data.frame(slp=array(aperm(s.icfdk,c(2,3,1))), sft=array(aperm(b.icfdk,c(2,3,1))), cf=factor(rep(c('Low','High','Low','High'),ns),levels=c('Low','High')),
                       dk=factor(rep(c('Deck 4','Deck 4','Deck 6','Deck 6'),ns)), s=factor(array(t(array(rep(1:ns,2*2),dim=c(ns,2*2))))))

art.dat.figs1f.slp = art(slp~cf*dk+Error(s+s:cf+s:dk+s:cf:dk), data=dat.figs1f)
anova(art.dat.figs1f.slp)
wilcoxsign_test(s.icfdk[,1,1]~s.icfdk[,2,1])
wilcoxsign_test(s.icfdk[,1,2]~s.icfdk[,2,2])

art.dat.figs1f.sft = art(sft~cf*dk+Error(s+s:cf+s:dk+s:cf:dk), data=dat.figs1f)
anova(art.dat.figs1f.sft)
wilcoxsign_test(b.icfdk[,1,1]~b.icfdk[,1,2])
wilcoxsign_test(b.icfdk[,2,1]~b.icfdk[,2,2])

fig2fl = ggplot(data=dat.figs1f, aes(x=interaction(cf,dk), y=slp, fill=interaction(cf,dk))) +
        stat_boxplot(geom='errorbar',size=0.5, width=0.25) +
        geom_boxplot(aes(fill=interaction(cf,dk)), outlier.shape=3, outlier.size=1.2) +
        scale_fill_manual(values=colors[1:4]) +
        scale_x_discrete(guide="axis_nested") +
        scale_y_continuous(limits=c(-5,0), labels=scales::number_format(accuracy=0.1), , name="Slope") + 
        geom_signif(comparisons=list(c(1,2),c(3,4)), annotations=c('***'), y_position=-0.4, tip_length=0, vjust=0.2, textsize=4.5) +
        theme_custom +
        theme(axis.title.x=element_blank(), axis.text.x=element_text(angle=30, hjust=1, vjust=1, margin=margin(t=0.05, b=0, unit="cm")),
              ggh4x.axis.nestline.x=element_line(linewidth=0.25), ggh4x.axis.nesttext.x=element_text(angle=30, hjust=0.5, vjust=0.5, margin=margin(t=0.09, unit="cm")))
fig2fl

fig2fr = ggplot(data=dat.figs1f, aes(x=interaction(cf,dk), y=sft, fill=interaction(cf,dk))) +
        stat_boxplot(geom='errorbar',size=0.5, width=0.25) +
        geom_boxplot(aes(fill=interaction(cf,dk)), outlier.shape=3, outlier.size=1.2) +
        scale_fill_manual(values=colors[1:4]) +
        scale_x_discrete(guide="axis_nested") +
        scale_y_continuous(limits=c(10.5,15), breaks=seq(10.5,15,1.5), name="Shift") + 
        geom_signif(comparisons=list(c(1,3),c(2,4)), annotations=c('***'), y_position=c(10.9,10.5), tip_length=0, vjust=1.7, textsize=4.5) +
        theme_custom +
        theme(axis.title.x=element_blank(), axis.text.x=element_text(angle=30, hjust=1, vjust=1, margin=margin(t=0.05, b=0, unit="cm")),
              ggh4x.axis.nestline.x=element_line(linewidth=0.25), ggh4x.axis.nesttext.x=element_text(angle=30, hjust=0.5, vjust=0.5, margin=margin(t=0.09, unit="cm")))
fig2fr


## Fig S1gh
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

dat = dump.format(list(open=data$open, cf=data$kcf.hl, dk=data$dkest, di=data$deci, sidx=data$sidx, ns=ns, N=N))
vals = c("a.mu","b.mu","a.sig","b.sig","a.p","b.p")

rgr.kcfdk = run.jags(model=paste0(mdl.folder,"SimpleLogisticRegress.txt"), monitor=vals, data=dat, n.chains=3, inits=c(inits1,inits2,inits3),
                     plots=FALSE, method="parallel", burnin=5000, sample=1000, adapt=1000, thin=10)
summary.rgr.kcfdk = add.summary(rgr.kcfdk)

s.kcfdk = array(sapply(1:4,function(c){sapply(1:ns,function(i){summary.rgr.kcfdk$summaries[paste0("a.p[",i,",",c,"]"),"Median"]})}),dim=c(ns,2,2))
b.kcfdk = array(sapply(1:4,function(c){sapply(1:ns,function(i){summary.rgr.kcfdk$summaries[paste0("b.p[",i,",",c,"]"),"Median"]})}),dim=c(ns,2,2))

x = seq(8,18,0.1)
simphit.kcfdk = array(NaN,dim=c(ns,length(x),2,2))
for (i in 1:ns) {
    for (c in 0:1) {
        for (d in 0:1) {
            simphit.kcfdk[i,,c+1,d+1] = ilogit(s.kcfdk[i,c+1,d+1]*(x-b.kcfdk[i,c+1,d+1]))
        }
    }
}

dat.figs1g = data.frame(prp=array(aperm(phit.kcfdk,c(2,3,4,1))), opn=rep(8:18,2*2*ns),
                       cf=factor(rep(c(rep(0,11),rep(1,11)),2*ns)), dk=factor(rep(c(rep(0,11*2),rep(1,11*2)),ns)))
dat.figs1g = dat.figs1g[!is.na(dat.figs1g$prp),]

sumdat.figs1g = data.frame(prp=array(apply(phit.kcfdk,c(2,3,4),nanmean)),
                          s=array(apply(phit.kcfdk,c(2,3,4),nansd)/sqrt(apply(!is.na(phit.kcfdk),c(2,3,4),nansum))),
                          opn=rep(8:18,2*2), cf=factor(rep(c(rep(0,11),rep(1,11)),2)), dk=factor(c(rep(0,11*2),rep(1,11*2))))

simdat.figs1g = data.frame(prp=array(apply(simphit.kcfdk,c(2,3,4),nanmean)),
                          s=array(apply(simphit.kcfdk,c(2,3,4),nansd)/sqrt(apply(!is.na(simphit.kcfdk),c(2,3,4),nansum))),
                          opn=rep(x,2*2), cf=factor(rep(c(rep(0,length(x)),rep(1,length(x))),2)), dk=factor(c(rep(0,length(x)*2),rep(1,length(x)*2))))

fig2g = ggplot(data=sumdat.figs1g, aes(x=opn, y=prp, group=interaction(cf,dk), fill=interaction(cf,dk), color=interaction(cf,dk))) +
        geom_line(data=simdat.figs1g, aes(x=opn, y=prp, group=interaction(cf,dk), color=interaction(cf,dk)), linewidth=0.6, alpha=0.7) +
        stat_boxplot(data=dat.figs1g, aes(x=opn, y=prp, group=interaction(opn,cf,dk), color=interaction(cf,dk)), geom='errorbar', linewidth=0.1, width=0.1, alpha=0.4, position=position_dodge(width=0.75)) +
        geom_boxplot(data=dat.figs1g, aes(x=opn, y=prp, group=interaction(opn,cf,dk), fill=interaction(cf,dk), color=interaction(cf,dk)), alpha=0.4, linewidth=0.15, width=0.25, outlier.shape=3, outlier.size=0.5, position=position_dodge(width=0.75)) +
        geom_point(aes(fill=interaction(cf,dk)), shape=21, position=position_dodge(0.75), color="white", size=2.5) +
        scale_fill_manual(name="Deck conf.\n\n", values=colors[1:4], labels=c("",""," Deck 4 "," Deck 6"),
                          guide=guide_legend(override.aes=list(color=colors[c(1,3,2,4)], shape=c(15,15,15,15), size=5),ncol=2)) +
        scale_color_manual(name=NULL, values=colors[1:4], guide="none") +
        scale_x_continuous(breaks=seq(8,18,2), name="Face-up card score") +
        scale_y_continuous(limits=c(0,1), breaks=seq(0,1,0.2), name="Proportion of Hit") +
        theme_custom +
        theme(legend.position=c(1,1),legend.justification=c(1,1),legend.background=element_rect(colour="transparent"),
              legend.title=element_text(size=10),
              legend.text=element_text(size=10, margin=margin(r=0,unit='cm')),legend.spacing.x=unit(-0.02,'cm'),legend.spacing.y=unit(0,'cm'),
              legend.key.width=unit(0.5,"cm"), legend.key.height=unit(0.3,"cm"))
fig2g

dat.figs1h = data.frame(slp=array(aperm(s.kcfdk,c(2,3,1))), sft=array(aperm(b.kcfdk,c(2,3,1))), cf=factor(rep(c('Low','High','Low','High'),ns),levels=c('Low','High')),
                       dk=factor(rep(c('Deck 4','Deck 4','Deck 6','Deck 6'),ns)), s=factor(array(t(array(rep(1:ns,2*2),dim=c(ns,2*2))))))
art.dat.figs1h.slp = art(slp~cf*dk+Error(s+s:cf+s:dk+s:cf:dk), data=dat.figs1h)
anova(art.dat.figs1h.slp)
wilcoxsign_test(s.kcfdk[,1,1]~s.kcfdk[,2,1])
wilcoxsign_test(s.kcfdk[,1,2]~s.kcfdk[,2,2])
wilcoxsign_test(s.kcfdk[,1,1]~s.kcfdk[,1,2])
wilcoxsign_test(s.kcfdk[,2,1]~s.kcfdk[,2,2])

art.dat.figs1h.sft = art(sft~cf*dk+Error(s+s:cf+s:dk+s:cf:dk), data=dat.figs1h)
anova(art.dat.figs1h.sft)
wilcoxsign_test(b.kcfdk[,1,1]~b.kcfdk[,2,1])
wilcoxsign_test(b.kcfdk[,1,2]~b.kcfdk[,2,2])
wilcoxsign_test(b.kcfdk[,1,1]~b.kcfdk[,1,2])
wilcoxsign_test(b.kcfdk[,2,1]~b.kcfdk[,2,2])

fig2hl = ggplot(data=dat.figs1h, aes(x=interaction(cf,dk), y=slp, fill=interaction(cf,dk))) +
        stat_boxplot(geom='errorbar',size=0.5, width=0.25) +
        geom_boxplot(aes(fill=interaction(cf,dk)), outlier.shape=3, outlier.size=1.2) +
        scale_fill_manual(values=colors[1:4]) +
        scale_x_discrete(guide="axis_nested") +
        scale_y_continuous(limits=c(-3.5,-0.5),breaks=seq(-3.5,-0.5,0.5), name="Slope") + 
        geom_signif(comparisons=list(c(1,3)), annotations=c('**'), y_position=-3, tip_length=0, vjust=1.7, textsize=4.5) +
        theme_custom +
        theme(axis.title.x=element_blank(), axis.text.x=element_text(angle=30, hjust=1, vjust=1, margin=margin(t=0.05, b=0, unit="cm")),
              ggh4x.axis.nestline.x=element_line(linewidth=0.25), ggh4x.axis.nesttext.x=element_text(angle=30, hjust=0.5, vjust=0.5, margin=margin(t=0.09, unit="cm")))
fig2hl

fig2hr = ggplot(data=dat.figs1h, aes(x=interaction(cf,dk), y=sft, fill=interaction(cf,dk))) +
        stat_boxplot(geom='errorbar',size=0.5, width=0.25) +
        geom_boxplot(aes(fill=interaction(cf,dk)), outlier.shape=3, outlier.size=1.2) +
        scale_fill_manual(values=colors[1:4]) +
        scale_x_discrete(guide="axis_nested") +
        scale_y_continuous(limits=c(9.5,15.5), breaks=seq(9.5,15.5,1.5), name="Shift") + 
        geom_signif(comparisons=list(c(1,2),c(3,4)), annotations=c('**','**'), y_position=14.8, tip_length=0, vjust=0.2, textsize=4.5) +
        geom_signif(comparisons=list(c(1,3),c(2,4)), annotations=c('***','***'), y_position=c(10.7,9.8), tip_length=0, vjust=1.7, textsize=4.5) +
        theme_custom +
        theme(axis.title.x=element_blank(), axis.text.x=element_text(angle=30, hjust=1, vjust=1, margin=margin(t=0.05, b=0, unit="cm")),
              ggh4x.axis.nestline.x=element_line(linewidth=0.25), ggh4x.axis.nesttext.x=element_text(angle=30, hjust=0.5, vjust=0.5, margin=margin(t=0.09, unit="cm")))
fig2hr


## Fig S1ij
icf.kcf = array(NaN, dim=c(ns,11,4))
for (i in 1:ns) {
    for (o in 8:18) {
        for (c in 1:4) {
            icf.kcf[i,o-7,c] = mean(data$icf[data$sidx==i&data$open==o&data$kcf==c])
        }
    }
}

dat = dump.format(list(open=data$open, cf=data$kcf, dk=data$dkest, di=data$deci, sidx=data$sidx, ns=ns, N=N))
vals = c("a.mu","b.mu","a.sig","b.sig","a.p","b.p")

rgr.phit.kcf = run.jags(model=paste0(mdl.folder,"SimpleLogisticRegress_cf4.txt"), monitor=vals, data=dat, n.chains=3, inits=c(inits1,inits2,inits3),
                        plots=FALSE, method="parallel", burnin=5000, sample=1000, adapt=1000, thin=10)
summary.rgr.phit.kcf = add.summary(rgr.phit.kcf)

sv = rep(NaN,N)
for (i in 1:ns) {
    for (c in 1:4) {
        tmpa = summary.rgr.phit.kcf$summaries[paste0("a.p[",i,",",c,"]"),"Median"]
        tmpb = summary.rgr.phit.kcf$summaries[paste0("b.p[",i,",",c,"]"),"Median"]
        sv[data$sidx==i&data$kcf==c] = abs(tmpa*(data$open[data$sidx==i&data$kcf==c]-tmpb))
    }
}

dat = dump.format(list(sv=sv, icf=data$icf, kcf=data$kcf, sidx=data$sidx, ns=ns, N=N))
vals = c("a.mu","b.mu","a.sig","b.sig","a.p","b.p","sig.mu")

rgr.kcf = run.jags(model=paste0(mdl.folder,"SimpleRegress.txt"), monitor=vals, data=dat, n.chains=3, inits=c(inits1,inits2,inits3),
                     plots=FALSE, method="parallel", burnin=5000, sample=1000, adapt=1000, thin=10)
summary.rgr.kcf = add.summary(rgr.kcf)

s.icf.kcf = sapply(1:4,function(c){sapply(1:ns,function(i){summary.rgr.kcf$summaries[paste0("a.p[",i,",",c,"]"),"Median"]})})
b.icf.kcf = sapply(1:4,function(c){sapply(1:ns,function(i){summary.rgr.kcf$summaries[paste0("b.p[",i,",",c,"]"),"Median"]})})

sim.icf.kcf = array(NaN, dim=c(ns,length(x),4))
for (i in 1:ns) {
    for (c in 1:4) {
        tmpa = summary.rgr.phit.kcf$summaries[paste0("a.p[",i,",",c,"]"),"Median"]
        tmpb = summary.rgr.phit.kcf$summaries[paste0("b.p[",i,",",c,"]"),"Median"]
        tmpsv = abs(tmpa*(x-tmpb))
        sim.icf.kcf[i,,c] = s.icf.kcf[i,c]*tmpsv + b.icf.kcf[i,c]
    }
}

dat.figs1i = data.frame(icf=array(aperm(icf.kcf,c(2,3,1))), opn=rep(8:18,4*ns),
                       cf=factor(rep(c(rep(1,11),rep(2,11),rep(3,11),rep(4,11)),ns)))
dat.figs1i = dat.figs1i[!is.na(dat.figs1i$icf),]

sumdat.figs1i = data.frame(icf=array(apply(icf.kcf,c(2,3),nanmean)), s=array(apply(icf.kcf,c(2,3),nansd)/sqrt(apply(!is.na(icf.kcf),c(2,3),nansum))),
                          opn=rep(8:18,4), cf=factor(c(rep(1,11),rep(2,11),rep(3,11),rep(4,11))))

simdat.figs1i = data.frame(icf=array(apply(sim.icf.kcf,c(2,3),nanmean)), s=array(apply(sim.icf.kcf,c(2,3),nansd)/sqrt(apply(!is.na(sim.icf.kcf),c(2,3),nansum))),
                          opn=rep(x,4), cf=factor(c(rep(1,length(x)),rep(2,length(x)),rep(3,length(x)),rep(4,length(x)))))

fig2i = ggplot(data=sumdat.figs1i, aes(x=opn, y=icf, group=cf, fill=cf, color=cf)) +
        geom_line(data=simdat.figs1i, aes(x=opn, y=icf, group=cf, color=cf), linewidth=0.6, alpha=0.5) +
        stat_boxplot(data=dat.figs1i, aes(x=opn, y=icf, group=interaction(opn,cf), color=interaction(cf)), geom='errorbar', linewidth=0.1, width=0.1, alpha=0.4, position=position_dodge(width=0.75)) +
        geom_boxplot(data=dat.figs1i, aes(x=opn, y=icf, group=interaction(opn,cf), fill=interaction(cf), color=interaction(cf)), alpha=0.4, linewidth=0.15, width=0.25, outlier.shape=3, outlier.size=0.5, position=position_dodge(width=0.75)) +
        geom_point(aes(fill=cf), shape=21, position=position_dodge(0.75), color="white", size=2.5) +
        scale_fill_manual(name="Deck conf lv", values=colors[5:8], labels=c("1","2","3","4"),
                          guide=guide_legend(override.aes=list(color=NA, shape=22, size=4),ncol=4,label.position="bottom",title.hjust=0.5)) +
        scale_color_manual(name=NULL, values=colors[5:8], guide="none") +
        scale_x_continuous(breaks=seq(8,18,2), name="Face-up card score") +
        scale_y_continuous(limits=c(1,4.05), breaks=1:4, labels=scales::number_format(accuracy=0.1), name="Averaged Decision\nconfidence level") +
        theme_custom +
        theme(legend.position=c(0.02,0.05),legend.justification=c(0,0),legend.background=element_rect(fill = "transparent"),
              legend.title=element_text(size=8),
              legend.text=element_text(size=8, margin=margin(r=0,unit='cm')),legend.spacing.x=unit(-0.05,'cm'),legend.spacing.y=unit(0.1,'cm'),
              legend.key.width=unit(0.3,"cm"), legend.key.height=unit(0.3,"cm"))
fig2i

dat.figs1j = data.frame(slp=array(t(s.icf.kcf)), bis=array(t(b.icf.kcf)), cf=factor(rep(1:4,ns)), s=factor(array(t(array(rep(1:ns,4),dim=c(ns,4))))))
art.dat.figs1j.slp = art(slp~cf+Error(s+s:cf), data=dat.figs1j)
anova(art.dat.figs1j.slp)
cor.test(as.numeric(dat.figs1j$cf),dat.figs1j$slp)

art.dat.figs1j.bis = art(bis~cf+Error(s+s:cf), data=dat.figs1j)
anova(art.dat.figs1j.bis)
cor.test(as.numeric(dat.figs1j$cf),dat.figs1j$bis)