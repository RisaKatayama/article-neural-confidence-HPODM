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
load(paste0(c(dat.folder,"bhvdata_bhv.RData"),collapse="/"))

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

    udc = (data$p.us==0.5)*1 # face-down card deck H = 1, deck L = 0
    undisc = data$undisc # actual score of the face-down card
    rw = data$rw

    trl = data$trl
    ses = data$ses
    sidx = data$sidx
    
    include = !is.na(data$deckest*data$cf.deck*data$decision*data$cf.deci)
    res = data.frame(dkest=dkest[include], dktrue=dktrue[include], kcf=kcf[include], kcf.hl=kcf.hl[include], add=add[include],
               deci=deci[include], open=open[include], undisc=undisc[include], udc=udc[include], icf=icf[include], icf.hl=icf.hl[include], 
               rw=rw[include], trl=trl[include], ses=ses[include], sidx=sidx[include])
    return(res)
}
data = load.data.basic(bdata)
ns = max(data$sidx)
ng = max(data$ses)
N = dim(data)[1]


#####################
## Fig S2ab
phit.uddk = array(NaN,dim=c(ns,11,2,2))
for (i in 1:ns) {
    for (o in 8:18) {
        for (u in 0:1) {
            for (d in 0:1) {
                phit.uddk[i,o-7,u+1,d+1] = mean(data$deci[data$open==o&data$dkest==d&data$udc==u&data$sidx==i])
            }
        }
    }
}

dat = dump.format(list(open=data$open, cf=data$udc, dk=data$dkest, di=data$deci, sidx=data$sidx, ns=ns, N=N))
vals = c("a.mu","b.mu","a.sig","b.sig","a.p","b.p")

rgr.uddk = run.jags(model="SimpleLogisticRegress.txt", monitor=vals, data=dat, n.chains=3, inits=c(inits1,inits2,inits3),
                     plots=FALSE, method="parallel", burnin=5000, sample=1000, adapt=1000, thin=10)
summary.rgr.uddk = add.summary(rgr.uddk)

s.uddk = array(sapply(1:4,function(c){sapply(1:ns,function(i){summary.rgr.uddk$summaries[paste0("a.p[",i,",",c,"]"),"Median"]})}),dim=c(ns,2,2))
b.uddk = array(sapply(1:4,function(c){sapply(1:ns,function(i){summary.rgr.uddk$summaries[paste0("b.p[",i,",",c,"]"),"Median"]})}),dim=c(ns,2,2))

x = seq(8,18,0.1)
simphit.uddk = array(NaN,dim=c(ns,length(x),2,2))
for (i in 1:ns) {
    for (u in 0:1) {
        for (d in 0:1) {
            simphit.uddk[i,,u+1,d+1] = ilogit(s.uddk[i,u+1,d+1]*(x-b.uddk[i,u+1,d+1]))
        }
    }
}

dat.figs2a = data.frame(prp=array(aperm(phit.uddk,c(2,3,4,1))), opn=rep(8:18,2*2*ns),
                       ud=factor(rep(c(rep(1,11),rep(0,11)),2*ns)), dk=factor(rep(c(rep(0,11*2),rep(1,11*2)),ns)))
dat.figs2a = dat.figs2a[!is.na(dat.figs2a$prp),]

sumdat.figs2a = data.frame(prp=array(apply(phit.uddk,c(2,3,4),nanmean)),
                       s=array(apply(phit.uddk,c(2,3,4),nansd)/sqrt(apply(!is.na(phit.uddk),c(2,3,4),nansum))),
                       opn=rep(8:18,2*2), ud=factor(rep(c(rep(1,11),rep(0,11)),2)), dk=factor(c(rep(0,11*2),rep(1,11*2))))

simdat.figs2a = data.frame(prp=array(apply(simphit.uddk,c(2,3,4),nanmean)),
                          s=array(apply(simphit.uddk,c(2,3,4),nansd)/sqrt(apply(!is.na(simphit.uddk),c(2,3,4),nansum))),
                          opn=rep(x,2*2), ud=factor(rep(c(rep(1,length(x)),rep(0,length(x))),2)), dk=factor(c(rep(0,length(x)*2),rep(1,length(x)*2))))

figs2a = ggplot(data=sumdat.figs2a, aes(x=opn, y=prp, group=interaction(ud,dk), fill=interaction(ud,dk), color=interaction(ud,dk))) +
        geom_line(data=simdat.figs2a, aes(x=opn, y=prp, group=interaction(ud,dk), color=interaction(ud,dk)), linewidth=0.6, alpha=0.7) +
        stat_boxplot(data=dat.figs2a, aes(x=opn, y=prp, group=interaction(opn,ud,dk), color=interaction(ud,dk)), geom='errorbar', linewidth=0.1, width=0.1, alpha=0.4, position=position_dodge(width=0.75)) +
        geom_boxplot(data=dat.figs2a, aes(x=opn, y=prp, group=interaction(opn,ud,dk), fill=interaction(ud,dk), color=interaction(ud,dk)), alpha=0.4, linewidth=0.15, width=0.25, outlier.shape=3, outlier.size=0.5, position=position_dodge(width=0.75)) +
        geom_point(aes(fill=interaction(ud,dk)), shape=21, position=position_dodge(0.75), color="white", size=2.5) +
        scale_fill_manual(name="Face-down\ncard deck\n\n", values=colors[1:4], labels=c("",""," Deck 4 "," Deck 6"),
                          guide=guide_legend(override.aes=list(color=colors[c(1,3,2,4)], shape=c(15,15,15,15), size=5),ncol=2)) +
        scale_color_manual(name=NULL, values=colors[1:4], guide="none") +
        scale_x_continuous(breaks=seq(8,18,2), name="Face-up card score") +
        scale_y_continuous(limits=c(0,1), breaks=seq(0,1,0.2), name="Proportion of Hit") +
        theme_custom +
        theme(legend.position=c(1,1),legend.justification=c(1,1),legend.background=element_rect(colour="transparent"),
              legend.title=element_text(size=10),
              legend.text=element_text(size=10, margin=margin(r=0,unit='cm')),legend.spacing.x=unit(-0.02,'cm'),legend.spacing.y=unit(0,'cm'),
              legend.key.width=unit(0.5,"cm"), legend.key.height=unit(0.3,"cm"))
figs2a

dat.figs2b = data.frame(slp=array(aperm(s.uddk,c(2,3,1))), sft=array(aperm(b.uddk,c(2,3,1))), ud=factor(rep(c('deck L','deck H','deck L','deck H'),ns),levels=c('deck H','deck L')),
                        dk=factor(rep(c('Deck 4','Deck 4','Deck 6','Deck 6'),ns)), s=factor(array(t(array(rep(1:ns,2*2),dim=c(ns,2*2))))))
art.dat.figs2bslp = art(slp~ud*dk+Error(s+s:ud+s:dk+s:ud:dk), data=dat.figs2b)
anova(art.dat.figs2b.slp)

art.dat.figs2b.sft = art(sft~ud*dk+Error(s+s:ud+s:dk+s:ud:dk), data=dat.figs2b)
anova(art.dat.figs2b.sft)

figs2bl = ggplot(data=dat.figs2b, aes(x=interaction(ud,dk), y=slp, fill=interaction(ud,dk))) +
        stat_boxplot(geom='errorbar',size=0.5, width=0.25) +
        geom_boxplot(aes(fill=interaction(ud,dk)), outlier.shape=3, outlier.size=1.2) +
        scale_fill_manual(values=colors[1:4]) +
        scale_x_discrete(guide="axis_nested") +
        scale_y_continuous(limits=c(-3.5,-0.4),breaks=seq(-3.5,-0.5,0.5), name="Slope") + 
        geom_signif(comparisons=list(c(1,3)), annotations=c('*'), y_position=-3.2, tip_length=0, vjust=2, textsize=4.5) +
        geom_signif(comparisons=list(c(3,4)), annotations=c('*'), y_position=-0.7, tip_length=0, vjust=-0.2, textsize=4.5) +
        theme_custom +
        theme(axis.title.x=element_blank(), axis.text.x=element_text(angle=30, hjust=1, vjust=1, margin=margin(t=0.05, b=0, unit="cm")),
              ggh4x.axis.nestline.x=element_line(linewidth=0.25), ggh4x.axis.nesttext.x=element_text(angle=30, hjust=0.5, vjust=0.5, margin=margin(t=0.09, unit="cm")))
figs2bl

figs2br = ggplot(data=dat.figs2b, aes(x=interaction(ud,dk), y=bis, fill=interaction(ud,dk))) +
        stat_boxplot(geom='errorbar',size=0.5, width=0.25) +
        geom_boxplot(aes(fill=interaction(ud,dk)), outlier.shape=3, outlier.size=1.2) +
        scale_fill_manual(values=colors[1:4]) +
        scale_x_discrete(guide="axis_nested") +
        scale_y_continuous(limits=c(10,15), breaks=seq(10,15,1), name="Bias") + 
        geom_signif(comparisons=list(c(1,2),c(3,4)), annotations=c('***','***'), y_position=14.7, tip_length=0, vjust=0.2, textsize=4.5) +
        theme_custom +
        theme(axis.title.x=element_blank(), axis.text.x=element_text(angle=30, hjust=1, vjust=1, margin=margin(t=0.05, b=0, unit="cm")),
              ggh4x.axis.nestline.x=element_line(linewidth=0.25), ggh4x.axis.nesttext.x=element_text(angle=30, hjust=0.5, vjust=0.5, margin=margin(t=0.09, unit="cm")))
figs2br


## Fig S2cd
icf.ud = array(NaN, dim=c(ns,11,2))
for (i in 1:ns) {
    for (o in 8:18) {
        for (u in 0:1) {
            icf.ud[i,o-7,u+1] = mean(data$icf[data$sidx==i&data$open==o&data$udc==u])
        }
    }
}

dat = dump.format(list(open=data$open, dk=data$udc, di=data$deci, sidx=data$sidx, ns=ns, N=N))
vals = c("a.mu","b.mu","a.sig","b.sig","a.p","b.p")

rgr.phit.ud = run.jags(model="SimpleLogisticRegress2.txt", monitor=vals, data=dat, n.chains=3, inits=c(inits1,inits2,inits3),
                     plots=FALSE, method="parallel", burnin=5000, sample=1000, adapt=1000, thin=10)
summary.rgr.phit.ud = add.summary(rgr.phit.ud)

sv = rep(NaN,N)
for (i in 1:ns) {
    for (u in 0:1) {
        tmpa = summary.rgr.phit.ud$summaries[paste0("a.p[",i,",",u+1,"]"),"Median"]
        tmpb = summary.rgr.phit.ud$summaries[paste0("b.p[",i,",",u+1,"]"),"Median"]
        sv[data$sidx==i&data$udc==u] = abs(tmpa*(data$open[data$sidx==i&data$udc==u]-tmpb))
    }
}

dat = dump.format(list(sv=sv, icf=data$icf, udc=data$udc, sidx=data$sidx, ns=ns, N=N))
vals = c("a.mu","b.mu","a.sig","b.sig","a.p","b.p")

rgr.ud = run.jags(model="SimpleRegress4.txt", monitor=vals, data=dat, n.chains=3, inits=c(inits1,inits2,inits3),
                     plots=FALSE, method="parallel", burnin=5000, sample=1000, adapt=1000, thin=10)
summary.rgr.ud = add.summary(rgr.ud)

s.icf.ud = sapply(1:2,function(c){sapply(1:ns,function(i){summary.rgr.ud$summaries[paste0("a.p[",i,",",c,"]"),"Median"]})})
b.icf.ud = sapply(1:2,function(c){sapply(1:ns,function(i){summary.rgr.ud$summaries[paste0("b.p[",i,",",c,"]"),"Median"]})})

sim.icf.ud = array(NaN, dim=c(ns,length(x),2))
for (i in 1:ns) {
    for (u in 1:2) {
        tmpa = summary.rgr.phit.ud$summaries[paste0("a.p[",i,",",u,"]"),"Median"]
        tmpb = summary.rgr.phit.ud$summaries[paste0("b.p[",i,",",u,"]"),"Median"]
        tmpsv = abs(tmpa*(x-tmpb))
        sim.icf.ud[i,,u] = s.icf.ud[i,u]*tmpsv + b.icf.ud[i,u]
    }
}

dat.figs2c = data.frame(icf=array(aperm(icf.ud,c(2,3,1))), opn=rep(8:18,ns*2), ud=factor(rep(c(rep('deck L',11*2),rep('deck H',11*2)),ns),levels=c('deck H','deck L')))
dat.figs2c = dat.figs2c[!is.na(dat.figs2c$icf),]

sumdat.figs2c = data.frame(icf=array(apply(icf.ud,c(2,3),nanmean)), s=array(apply(icf.ud,c(2,3),nansd)/sqrt(apply(!is.na(icf.ud),c(2,3),nansum))),
                           opn=rep(8:18,2), ud=factor(c(rep('deck L',11),rep('deck H',11)),levels=c('deck H','deck L')))

simdat.figs2c = data.frame(icf=array(apply(sim.icf.ud,c(2,3),nanmean)), s=array(apply(sim.icf.ud,c(2,3),nansd)/sqrt(apply(!is.na(sim.icf.ud),c(2,3),nansum))),
                           opn=rep(x,2), ud=factor(c(rep('deck L',length(x)),rep('deck H',length(x))),levels=c('deck H','deck L')))

figs2c = ggplot(data=sumdat.figs2c, aes(x=opn, y=icf, group=ud, fill=ud, color=ud)) +
        geom_line(data=simdat.figs2c, aes(x=opn, y=icf, group=ud, color=ud), linewidth=0.6, alpha=0.5) +
        stat_boxplot(data=dat.figs2c, aes(x=opn, y=icf, group=interaction(opn,ud), color=interaction(ud)), geom='errorbar', linewidth=0.1, width=0.1, alpha=0.4, position=position_dodge(width=0.75)) +
        geom_boxplot(data=dat.figs2c, aes(x=opn, y=icf, group=interaction(opn,ud), fill=interaction(ud), color=interaction(ud)), alpha=0.4, linewidth=0.15, width=0.25, outlier.shape=3, outlier.size=0.5, position=position_dodge(width=0.75)) +
        geom_point(aes(fill=ud), shape=21, position=position_dodge(0.75), color="white", size=2.5) +
        scale_fill_manual(name="Face-down\ncard deck", values=colors[c(5,8)], labels=c(" deck H"," deck L"),
                          guide=guide_legend(override.aes=list(color=NA, shape=22, size=4),label.position="right",title.hjust=0.5)) +
        scale_color_manual(name=NULL, values=colors[c(5,8)], guide="none") +
        scale_x_continuous(breaks=seq(8,18,2), name="Face-up card score") +
        scale_y_continuous(limits=c(1,4.05), breaks=1:4, labels=scales::number_format(accuracy=0.1), name="Averaged Decision\nconfidence level") +
        theme_custom +
        theme(legend.position=c(0.02,0.05),legend.justification=c(0,0),legend.background=element_rect(fill = "transparent"),
              legend.title=element_text(size=8), legend.title.align=0,
              legend.text=element_text(size=8, margin=margin(r=0,unit='cm')),legend.spacing.x=unit(-0.05,'cm'),legend.spacing.y=unit(0.1,'cm'),
              legend.key.width=unit(0.3,"cm"), legend.key.height=unit(0.3,"cm"))
figs2c

dat.figs2d = data.frame(slp=array(aperm(s.icf.ud,c(2,1))), bis=array(aperm(b.icf.ud,c(2,1))), ud=factor(rep(c('deck L','deck H'),ns),levels=c('deck H','deck L')),
                        s=factor(array(t(array(rep(1:ns,2),dim=c(ns,2))))))

figs2dl = ggplot(data=dat.figs2d, aes(x=ud, y=slp, fill=ud)) +
        stat_boxplot(geom='errorbar',size=0.5, width=0.25) +
        geom_boxplot(aes(fill=ud), outlier.shape=3, outlier.size=1.2) +
        scale_fill_manual(values=colors[c(5,8)]) +
        scale_x_discrete(name="Face-down\ncard deck") +
        scale_y_continuous(limits=c(0,0.6), name="Slope", labels=scales::number_format(accuracy=0.1)) + 
        theme_custom +
        theme(axis.text.x=element_text(angle=30, hjust=1, vjust=1, margin=margin(t=0.05, b=0, unit="cm")))
figs2dl

figs2dr = ggplot(data=dat.figs2d, aes(x=ud, y=bis, fill=ud)) +
        stat_boxplot(geom='errorbar',size=0.5, width=0.25) +
        geom_boxplot(aes(fill=ud), outlier.shape=3, outlier.size=1.2) +
        scale_fill_manual(values=colors[c(5,8)]) +
        scale_x_discrete(name="Face-down\ncard deck") +
        scale_y_continuous(limits=c(1,4), labels=scales::number_format(accuracy=0.1), name="Bias") + 
        theme_custom +
        theme(axis.text.x=element_text(angle=30, hjust=1, vjust=1, margin=margin(t=0.05, b=0, unit="cm")))
figs2dr