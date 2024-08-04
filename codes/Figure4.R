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
mdlres.folder = ""
ts.folder = ""
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

include = !is.na(mdata$deckest*mdata$cf.deck*mdata$decision*mdata$cf.deci)
data = mdata[include,]

incidx = c(1:20,23:24) # sbj 21/22 (original id=23,24) are excluded due to insufficient confidence rating
incidx.org = c(1:6,8:21,25:26)
ns = length(incidx)

event.ons = c(0,5,6.844,9.775,12.275,15.561,17.469) 
# averaged event onset; step 5-9 (in trial t) + step 2-3 (in trial t+1), step 5 in trial t as a reference onset


#####################
## Fig 4d (vAIC)
b.kcf.vaic = array(NaN, dim=c(ns,201))
scl.kcf = scaling.sbj(res$dk.res$valres$w.pred.kcf, data$sidx)

for (i in 1:ns){
    tmpkcf = scl.kcf[data$sidx==incidx[i]]
    f = readMat(paste0(ts.folder,"/s",incidx.org[i],"_l_vAIC_8mm.mat")) 
    # BOLD timeseries (0-1 rescaling & centering), Ntrl x 201 timepoints (0.1TR, [-2, 18]s from the decision onset [step 5])

    for (t in 1:201){
        tmpy = f$ts[,t]
        tmpfit = lm(y~x, data=data.frame(x=tmpkcf[!is.na(tmpy)],y=tmpy[!is.na(tmpy)]))
        b.kcf.vaic[i,t] = coef(tmpfit)[2]
    }
}

p.b.kcf.vaic = rep(NaN,201)
for (t in 1:201) {
   w = wilcoxsign_test(b.kcf.vaic[,t]~rep(0,ns),lower.tail=FALSE)
   p.b.kcf.vaic[t] = 1 - pnorm(w@statistic@teststatistic)
}
(which(p.adjust(p.b.kcf.vaic,"fdr")<0.05)-1)*0.1-2

t = seq(-2,18,0.1)
dat.fig4d1 = data.frame(t=t, m=apply(b.kcf.vaic,2,nanmean), s=apply(b.kcf.vaic,2,nansd)/sqrt(ns))
dat.fig4d1.p = data.frame(ti=c(-2,4.6,6.4,8.4), te=c(-1.8,5.9,7.2,10), t=c(0,0,0,0), m=c(0,0,0,0))

fig4d1 = ggplot(data=dat.fig4d1, aes(x=t, y=m)) +
        geom_rect(data=dat.fig4d1.p, aes(xmin=ti, xmax=te, ymin=-Inf, ymax=Inf), alpha=0.25) +
        geom_hline(yintercept=0, linewidth=0.2) +
        geom_vline(xintercept=event.ons[1], linewidth=0.2) +
        sapply(2:length(event.ons),function(x){geom_vline(xintercept=event.ons[x], linewidth=0.2, linetype="dotted")}) +
        geom_ribbon(aes(x=t, ymin=m-s, ymax=m+s), colour=NA, fill="#336600", alpha=0.5, show.legend=FALSE) +
        geom_line(aes(x=t, y=m), color="#336600", linewidth=0.6, alpha=0.9) +
        scale_x_continuous(limits=c(-2,18), breaks=seq(-2,18,2), name="Time from Decision onset (s)") +
        scale_y_continuous(limits=c(-0.02,0.06), name="Effect on BOLD (a.u.)") +
        theme_custom
fig4d1


## Fig 4d (dACC)
b.kcf.dacc = array(NaN, dim=c(ns,201))

for (i in 1:ns){
    tmpkcf = scl.kcf[data$sidx==incidx[i]]
    f = readMat(paste0(ts.folder,"/s",incidx.org[i],"_dACC_8mm.mat")) 
    # BOLD timeseries (0-1 rescaling & centering), Ntrl x 201 timepoints (0.1TR, [-2, 18]s from the decision onset [step 5])

    for (t in 1:201){
        tmpy = f$ts[,t]
        tmpfit = lm(y~x, data=data.frame(x=tmpkcf[!is.na(tmpy)],y=tmpy[!is.na(tmpy)]))
        b.kcf.dacc[i,t] = coef(tmpfit)[2]
    }
}

p.b.kcf.dacc = rep(NaN,201)
for (t in 1:201) {
   w = wilcoxsign_test(b.kcf.dacc[,t]~rep(0,ns),lower.tail=FALSE)
   p.b.kcf.dacc[t] = 1 - pnorm(w@statistic@teststatistic)
}
(which(p.adjust(p.b.kcf.dacc,"fdr")<0.05)-1)*0.1-2

t = seq(-2,18,0.1)
dat.fig4d2 = data.frame(t=t, m=apply(b.kcf.dacc,2,nanmean), s=apply(b.kcf.dacc,2,nansd)/sqrt(ns))
dat.fig4d2.p = data.frame(ti=c(1.3,5.1), te=c(4.8,10.6), t=c(0,0), m=c(0,0))

fig4d2 = ggplot(data=dat.fig4d2, aes(x=t, y=m)) +
        geom_rect(data=dat.fig4d2.p, aes(xmin=ti, xmax=te, ymin=-Inf, ymax=Inf), alpha=0.25) +
        geom_hline(yintercept=0, linewidth=0.2) +
        geom_vline(xintercept=event.ons[1], linewidth=0.2) +
        sapply(2:length(event.ons),function(x){geom_vline(xintercept=event.ons[x], linewidth=0.2, linetype="dotted")}) +
        geom_ribbon(aes(x=t, ymin=m-s, ymax=m+s), colour=NA, fill="#336600", alpha=0.5, show.legend=FALSE) +
        geom_line(aes(x=t, y=m), color="#336600", linewidth=0.6, alpha=0.9) +
        scale_x_continuous(limits=c(-2,18), breaks=seq(-2,18,2), name="Time from Decision onset (s)") +
        scale_y_continuous(limits=c(-0.02,0.06), name="Effect on BOLD (a.u.)") +
        theme_custom
fig4d2