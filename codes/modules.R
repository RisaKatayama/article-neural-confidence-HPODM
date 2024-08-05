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

nanmin = function(x){ min(x,na.rm=TRUE) }
nanmax = function(x){ max(x,na.rm=TRUE) }
nanmean = function(x){ mean(x,na.rm=TRUE) }
nansd = function(x){ sd(x,na.rm=TRUE) }
nanmedian = function(x){ median(x,na.rm=TRUE) }
nansum = function(x){ sum(x,na.rm=TRUE) }

ilogit = function(x) { 1/(1+exp(-x))}


calc.weighted.conf = function(pmat){

    w.pmat = sweep(pmat, 2, c(1:4), "*")
    estcf = apply(w.pmat, 1, sum)
    
    return(estcf)
}


model.validate = function(data, tsdata, model, pnames, modeltype){

    if (modeltype == "rea.dkmdl"){
        if (!is.null(data)){
            message("Training results calculating : DK Model : REA ...")  
            fitres = basic.dkmdl(data, model, pnames)
            if (!is.null(tsdata)) {
                message("Testing results calculating  : DK Model : REA ...")
                valres = basic.dkmdl(tsdata, model, pnames)
                valres = valres[names(valres)!="model"]
            }
        } else {
            message("Testing results calculating  : DK Model : REA ...")
            valres = basic.dkmdl(tsdata, model, pnames)
        }
    
    }else if (modeltype == "avg.dkmdl"){
        if (!is.null(data)){
            message("Training results calculating : DK Model : AVG ...")  
            fitres = basic.dkmdl(data, model, pnames)
            if (!is.null(tsdata)) {
                message("Testing results calculating  : DK Model : AVG ...")
                valres = basic.dkmdl(tsdata, model, pnames)
                valres = valres[names(valres)!="model"]
            }
        } else {
            message("Testing results calculating  : DK Model : AVG ...")
            valres = basic.dkmdl(tsdata, model, pnames)
        }
    
    }else if (modeltype == "bio.dkmdl"){
        if (!is.null(data)){
            message("Training results calculating : DK Model : BIO ...")  
            fitres = bio.dkmdl(data, model, pnames)
            if (!is.null(tsdata)) {
                message("Testing results calculating  : DK Model : BIO ...")
                valres = bio.dkmdl(tsdata, model, pnames)
                valres = valres[names(valres)!="model"]
            }
        } else {
            message("Testing results calculating  : DK Model : BIO ...")
            valres = bio.dkmdl(tsdata, model, pnames)
        }

    } else if (modeltype == "cwhdm.dimdl") {
        if (!is.null(data)){
            message("Training results calculating : DI Model : cwHDM ...")  
            fitres = cwhdm.dimdl(data, model, pnames)
            if (!is.null(tsdata)) {
                message("Testing results calculating  : DI Model : cwHDM ...")
                valres = cwhdm.dimdl(tsdata, model, pnames)
                valres = valres[names(valres)!="model"]
            }
        } else {
            message("Testing results calculating  : DI Model : cwHDM ...")
            valres = cwhdm.dimdl(tsdata, model, pnames)
        }

    } else if (modeltype == "fphdm.dimdl") {
        if (!is.null(data)){
            message("Training results calculating : DI Model : fpHDM ...")  
            fitres = fphdm.dimdl(data, model, pnames)
            if (!is.null(tsdata)) {
                message("Testing results calculating  : DI Model : fpHDM ...")
                valres = fphdm.dimdl(tsdata, model, pnames)
                valres = valres[names(valres)!="model"]
            }
        } else {
            message("Testing results calculating  : DI Model : fpHDM ...")
            valres = fphdm.dimdl(tsdata, model, pnames)
        }

    } else if (modeltype == "hdm.dimdl") {
        if (!is.null(data)){
            message("Training results calculating : DI Model : HDM ...")  
            fitres = basic.dimdl(data, model, pnames)
            if (!is.null(tsdata)) {
                message("Testing results calculating  : DI Model : HDM ...")
                valres = basic.dimdl(tsdata, model, pnames)
                valres = valres[names(valres)!="model"]
            }
        } else {
            message("Testing results calculating  : DI Model : HDM ...")
            valres = basic.dimdl(tsdata, model, pnames)
        }

    } else if (modeltype == "pdm.dimdl") {
        if (!is.null(data)){
            message("Training results calculating : DI Model : PDM ...")  
            fitres = basic.dimdl(data, model, pnames)
            if (!is.null(tsdata)) {
                message("Testing results calculating  : DI Model : PDM ...")
                valres = basic.dimdl(tsdata, model, pnames)
                valres = valres[names(valres)!="model"]
            }
        } else {
            message("Testing results calculating  : DI Model : PDM ...")
            valres = basic.dimdl(tsdata, model, pnames)
        }
    }

    if (!is.null(data)&!is.null(tsdata)){
        res = list(fitres=fitres, valres=valres)
    } else if (!is.null(data)) {
        res = list(fitres=fitres)
    } else if (!is.null(tsdata)) {
        res = list(valres=valres)
    }

    return(res)
}


basic.dkmdl = function(data, model, pnames){

    ns = max(data$sidx)
    N = dim(data)[1]
    nsmp = dim(model[[1]])[1]
    b.dk = array(NaN, dim=c(3*nsmp,N,2))
    b.kcf = array(NaN, dim=c(3*nsmp,N,2))
    theta.kcf = array(NaN, dim=c(3*nsmp,N,3))

    mdl.mcmc = rbind(model[[1]],model[[2]],model[[3]])
    for (i in 1:ns) {
        n = sum(data$sidx==i)
        b.dk[,data$sidx==i,1] = array(rep(mdl.mcmc[,pnames==paste0("b0.dk.p[",i,"]")],n), dim=c(3*nsmp,n))
        b.dk[,data$sidx==i,2] = array(rep(mdl.mcmc[,pnames==paste0("b1.dk.p[",i,"]")],n), dim=c(3*nsmp,n))

        b.kcf[,data$sidx==i,1] = array(rep(mdl.mcmc[,pnames==paste0("b1.kcf.p[",i,"]")],n), dim=c(3*nsmp,n))
        b.kcf[,data$sidx==i,2] = array(rep(mdl.mcmc[,pnames==paste0("b2.kcf.p[",i,"]")],n), dim=c(3*nsmp,n))

        a0 = array(rep(mdl.mcmc[,pnames==paste0("a0.kcf.p[",i,"]")],n), dim=c(3*nsmp,n))
        a1 = array(rep(mdl.mcmc[,"a1.kcf.mu"],n), dim=c(3*nsmp,n))
        theta.kcf[,data$sidx==i,1] = a0-a1
        theta.kcf[,data$sidx==i,2] = a0
        theta.kcf[,data$sidx==i,3] = a0+a1
    }
    m.b.dk = apply(b.dk,c(2,3),median)
    m.b.kcf = apply(b.kcf,c(2,3),median)
    m.theta.kcf = apply(theta.kcf,c(2,3),median)
  
    # deck estim
    message("deck inference phase...")

    X = cbind(array(1,dim=c(N,1)),data$de)
    Z.dk = apply(sweep(b.dk,c(2,3),X,"*"), c(1,2), sum)
    p.d6 = ilogit(Z.dk)
    p.dk = sweep((1-ilogit(Z.dk)),2,(1-data$dkest),"*") + sweep(ilogit(Z.dk),2,data$dkest,"*")
    
    z.dk = apply(X*m.b.dk, 1, sum)
    p.d6 = ilogit(z.dk)
    pred.dk = (p.d6>0.5)*1
    pred.dk[p.d6==0.5] = data$dkest[p.d6==0.5]
    
    # deck conf
    message("deck confidence phase...")

    X = cbind(abs(data$de), data$vr) ##
    Z.kcf = apply(sweep(b.kcf,c(2,3),X,"*"), c(1,2), sum)
    Q = ilogit(sweep(theta.kcf,c(1,2),Z.kcf,"+"))
    p.kcf.mat = array(NaN,dim=c(3*nsmp,N,4))
    p.kcf.mat[,,1] = 1 - Q[,,3]
    p.kcf.mat[,,2] = Q[,,3] - Q[,,2]
    p.kcf.mat[,,3] = Q[,,2] - Q[,,1]
    p.kcf.mat[,,4] = Q[,,1]
    cmat = sapply(1:4, function(j){ data$kcf==j })
    p.kcf = apply(sweep(p.kcf.mat, c(2,3), cmat, "*"), c(1,2), sum)

    z.kcf = apply(X*m.b.kcf, 1, sum)
    Q = ilogit(sweep(m.theta.kcf,1,z.kcf,"+"))
    p.kcf.mat = array(NaN,dim=c(N,4))
    p.kcf.mat[,1] = 1 - Q[,3]
    p.kcf.mat[,2] = Q[,3] - Q[,2]
    p.kcf.mat[,3] = Q[,2] - Q[,1]
    p.kcf.mat[,4] = Q[,1]
    pred.kcf = mapply(prediction, 1:N, MoreArgs=list(pmat=p.kcf.mat,y=data$kcf))
    w.pred.kcf = calc.weighted.conf(p.kcf.mat)

    res = list(pmat.dk=array(c(1-p.d6,p.d6),dim=c(N,2)), pred.dk=pred.dk, pmat.kcf=p.kcf.mat, pred.kcf=pred.kcf, w.pred.kcf=w.pred.kcf,
               loglik=list(dk=log(p.dk), kcf=log(p.kcf)))
    
    return(res)
}


bio.dkmdl = function(data, model, pnames){

    ns = max(data$sidx)
    N = dim(data)[1]
    nsmp = dim(model[[1]])[1]
    b.kcf = array(NaN, dim=c(3*nsmp,N,2))
    theta.kcf = array(NaN, dim=c(3*nsmp,N,3))
    
    mdl.mcmc = rbind(model[[1]],model[[2]],model[[3]])
    for (i in 1:ns) {
        n = sum(data$sidx==i)
        
        b.kcf[,data$sidx==i,1] = array(rep(mdl.mcmc[,pnames==paste0("b1.kcf.p[",i,"]")],n), dim=c(3*nsmp,n))
        b.kcf[,data$sidx==i,2] = array(rep(mdl.mcmc[,pnames==paste0("b2.kcf.p[",i,"]")],n), dim=c(3*nsmp,n))

        a0 = array(rep(mdl.mcmc[,pnames==paste0("a0.kcf.p[",i,"]")],n), dim=c(3*nsmp,n))
        a1 = array(rep(mdl.mcmc[,"a1.kcf.mu"],n), dim=c(3*nsmp,n))
        theta.kcf[,data$sidx==i,1] = a0-a1
        theta.kcf[,data$sidx==i,2] = a0
        theta.kcf[,data$sidx==i,3] = a0+a1
    }
    m.b.kcf = apply(b.kcf,c(2,3),median)
    m.theta.kcf = apply(theta.kcf,c(2,3),median)
    
    # deck estim
    message("deck estimation phase...")

    p.d6 = ilogit(t(array(rep(data$pd6,3*nsmp),dim=c(N,3*nsmp))))
    p.dk = sweep((1-p.d6),2,(1-data$dkest),"*") + sweep(p.d6,2,data$dkest,"*")
    
    pred.dk = (data$pd6>0.5)*1
    pred.dk[data$p6.bio==0.5] = data$dkest[data$p6.bio==0.5]
    
    # deck conf
    message("deck confidence phase...")

    X = cbind(data$ent, data$vr) ##
    Z.kcf = apply(sweep(b.kcf,c(2,3),X,"*"), c(1,2), sum)
    Q = ilogit(sweep(theta.kcf,c(1,2),Z.kcf,"+"))
    p.kcf.mat = array(NaN,dim=c(3*nsmp,N,4))
    p.kcf.mat[,,1] = 1 - Q[,,3]
    p.kcf.mat[,,2] = Q[,,3] - Q[,,2]
    p.kcf.mat[,,3] = Q[,,2] - Q[,,1]
    p.kcf.mat[,,4] = Q[,,1]
    cmat = sapply(1:4, function(j){ data$kcf==j })
    p.kcf = apply(sweep(p.kcf.mat,c(2,3),cmat,"*"), c(1,2), sum)

    z.kcf = apply(X*m.b.kcf, 1, sum)
    Q = ilogit(sweep(m.theta.kcf,1,z.kcf,"+"))
    p.kcf.mat = array(NaN,dim=c(N,4))
    p.kcf.mat[,1] = 1 - Q[,3]
    p.kcf.mat[,2] = Q[,3] - Q[,2]
    p.kcf.mat[,3] = Q[,2] - Q[,1]
    p.kcf.mat[,4] = Q[,1]
    pred.kcf = mapply(prediction, 1:N, MoreArgs=list(pmat=p.kcf.mat,y=data$kcf))
    w.pred.kcf = calc.weighted.conf(p.kcf.mat)

    res = list(pmat.dk=array(c(1-data$pd6,data$pd6),dim=c(N,2)), pred.dk=pred.dk, pmat.kcf=p.kcf.mat, pred.kcf=pred.kcf, w.pred.kcf=w.pred.kcf,
                loglik=list(dk=log(p.dk), kcf=log(p.kcf)))
   
    return(res)
}


cwhdm.dimdl = function(data, model, pnames){

    ns = max(data$sidx)
    N = dim(data)[1]
    nsmp = dim(model[[1]])[1]
    b.di = array(NaN, dim=c(3*nsmp,N,3))
    g.kcf = array(NaN, dim=c(3*nsmp,N))
    e.kcf = array(NaN, dim=c(3*nsmp,N))
    b.icf = array(NaN, dim=c(3*nsmp,N,2))
    theta.icf = array(NaN, dim=c(3*nsmp,N,3))
  
    mdl.mcmc = rbind(model[[1]],model[[2]],model[[3]])
    for (i in 1:ns) {
        n = sum(data$sidx==i)
        b.di[,data$sidx==i,1] = array(rep(mdl.mcmc[,pnames==paste0("b0.di.p[",i,"]")],n), dim=c(3*nsmp,n))
        b.di[,data$sidx==i,2] = array(rep(mdl.mcmc[,pnames==paste0("b1.di.p[",i,"]")],n), dim=c(3*nsmp,n))
        b.di[,data$sidx==i,3] = array(rep(mdl.mcmc[,pnames==paste0("b2.di.p[",i,"]")],n), dim=c(3*nsmp,n))
    
        g.kcf[,data$sidx==i] = array(rep(mdl.mcmc[,pnames==paste0("g.kcf.p[",i,"]")],n), dim=c(3*nsmp,n))
        e.kcf[,data$sidx==i] = array(rep(mdl.mcmc[,pnames==paste0("e.kcf.p[",i,"]")],n), dim=c(3*nsmp,n))
    
        b.icf[,data$sidx==i,1] = array(rep(mdl.mcmc[,pnames==paste0("b1.icf.p[",i,"]")],n), dim=c(3*nsmp,n))
        b.icf[,data$sidx==i,2] = array(rep(mdl.mcmc[,pnames==paste0("b2.icf.p[",i,"]")],n), dim=c(3*nsmp,n))
    
        a0 = array(rep(mdl.mcmc[,pnames==paste0("a0.icf.p[",i,"]")],n), dim=c(3*nsmp,n))
        a1 = array(rep(mdl.mcmc[,"a1.icf.mu"],n), dim=c(3*nsmp,n))
        theta.icf[,data$sidx==i,1] = a0-a1
        theta.icf[,data$sidx==i,2] = a0
        theta.icf[,data$sidx==i,3] = a0+a1
    }
    m.b.di = apply(b.di,c(2,3),median)
    m.g.kcf = apply(g.kcf,2,median)
    m.e.kcf = apply(e.kcf,2,median)
    m.b.icf = apply(b.icf,c(2,3),median)
    m.theta.icf = apply(theta.icf,c(2,3),median)
  
    # decision
    message("decision phase...")

    X = array(1,dim=c(3*nsmp,N,3))
    w.cdk = e.kcf - sweep(g.kcf,2,(1-data$kcf),"*")
    modadd = sweep(w.cdk,2,data$cdk,"*") + sweep(1-w.cdk,2,data$udk,"*")
    Vh = 21 - sweep(modadd,2,data$open+data$ud,"+")
    X[,,2] = Vh
    X[,,3] = Vh * abs(Vh)
    Z.di = apply(X*b.di, c(1,2), sum)
    p.hit.mat = ilogit(Z.di)
    p.di = sweep((1-p.hit.mat),2,(1-data$deci),"*") + sweep(p.hit.mat,2,data$deci,"*")

    w.cdk = m.e.kcf - m.g.kcf*(1-data$kcf)
    modadd = w.cdk*data$cdk + (1-w.cdk)*data$udk
    vh = 21 - (data$open + data$ud + modadd)
    x = cbind(array(1,dim=c(N,1)), vh, vh*abs(vh))
    z.di = apply(x*m.b.di, 1, sum)
    p.hit = ilogit(z.di)
    pred.di = (p.hit>0.5)*1
    pred.di[p.hit==0.5] = data$deci[p.hit==0.5]
    p.ch = (1-p.hit)*(1-data$deci) + p.hit*data$deci
  
    # decision confidence
    message("decision confidence phase...")
  
    SV = abs(Z.di)
    Z.icf = b.icf[,,1]*SV + sweep(b.icf[,,2],2,data$kcf,"*")
    Q = ilogit(sweep(theta.icf,c(1,2),Z.icf,"+"))
    p.icf.mat = array(NaN,dim=c(3*nsmp,N,4))
    p.icf.mat[,,1] = 1 - Q[,,3]
    p.icf.mat[,,2] = Q[,,3] - Q[,,2]
    p.icf.mat[,,3] = Q[,,2] - Q[,,1]
    p.icf.mat[,,4] = Q[,,1]
    cmat = sapply(1:4, function(j){ data$icf==j })
    p.icf =apply(sweep(p.icf.mat,c(2,3),cmat,"*"), c(1,2), sum)
  
    sv = abs(z.di)
    z.icf = m.b.icf[,1]*sv + m.b.icf[,2]*data$kcf
    Q = ilogit(sweep(m.theta.icf,1,z.icf,"+"))
    p.icf.mat = array(NaN,dim=c(N,4))
    p.icf.mat[,1] = 1 - Q[,3]
    p.icf.mat[,2] = Q[,3] - Q[,2]
    p.icf.mat[,3] = Q[,2] - Q[,1]
    p.icf.mat[,4] = Q[,1]
    pred.icf = mapply(prediction, 1:N, MoreArgs=list(pmat=p.icf.mat,y=data$icf))
    w.pred.icf = calc.weighted.conf(p.icf.mat)
  
    res = list(pmat.di=array(c(1-p.hit,p.hit),dim=c(N,2)), pred.di=pred.di, pmat.icf=p.icf.mat, pred.icf=pred.icf, w.pred.icf=w.pred.icf, 
               pred.vh.kcf=vh, pred.sv=sv, pred.util.icf=z.icf, modadd=modadd, w.cdk=w.cdk, g.kcf=m.g.kcf, e.kcf=m.e.kcf,
               loglik=list(di=log(p.di), icf=log(p.icf)))
  
    return(res)
}


fphdm.dimdl = function(data, model, pnames){
    
    ns = max(data$sidx)
    N = dim(data)[1]
    nsmp = dim(model[[1]])[1]
    b.di = array(NaN, dim=c(3*nsmp,N,3))
    b.icf = array(NaN, dim=c(3*nsmp,N,2))
    theta.icf = array(NaN, dim=c(3*nsmp,N,3))
  
    mdl.mcmc = rbind(model[[1]],model[[2]],model[[3]])
    for (i in 1:ns) {
        n = sum(data$sidx==i)
        b.di[,data$sidx==i,1] = array(rep(mdl.mcmc[,pnames==paste0("b0.di.p[",i,"]")],n), dim=c(3*nsmp,n))
        b.di[,data$sidx==i,2] = array(rep(mdl.mcmc[,pnames==paste0("b1.di.p[",i,"]")],n), dim=c(3*nsmp,n))
        b.di[,data$sidx==i,3] = array(rep(mdl.mcmc[,pnames==paste0("b2.di.p[",i,"]")],n), dim=c(3*nsmp,n))
        
        b.icf[,data$sidx==i,1] = array(rep(mdl.mcmc[,pnames==paste0("b1.icf.p[",i,"]")],n), dim=c(3*nsmp,n))
        b.icf[,data$sidx==i,2] = array(rep(mdl.mcmc[,pnames==paste0("b2.icf.p[",i,"]")],n), dim=c(3*nsmp,n))
        
        a0 = array(rep(mdl.mcmc[,pnames==paste0("a0.icf.p[",i,"]")],n), dim=c(3*nsmp,n))
        a1 = array(rep(mdl.mcmc[,"a1.icf.mu"],n), dim=c(3*nsmp,n))
        theta.icf[,data$sidx==i,1] = a0-a1
        theta.icf[,data$sidx==i,2] = a0
        theta.icf[,data$sidx==i,3] = a0+a1
    }
    m.b.di = apply(b.di,c(2,3),median)
    m.b.icf = apply(b.icf,c(2,3),median)
    m.theta.icf = apply(theta.icf,c(2,3),median)
    
    # decision
    message("decision phase...")

    X = array(1,dim=c(3*nsmp,N,3))
    X[,,2] = rep(1,3*nsmp) %o% (21 - (data$open + data$ud + data$add))
    X[,,3] = X[,,2] * abs(X[,,2])
    Z.di = apply(X*b.di, c(1,2), sum)
    p.hit.mat = ilogit(Z.di)
    p.di = sweep((1-p.hit.mat),2,(1-data$deci),"*") + sweep(p.hit.mat,2,data$deci,"*")

    vh = 21 - (data$open + data$ud + data$add)
    x = cbind(array(1,dim=c(N,1)), vh, vh*abs(vh))
    z.di = apply(x*m.b.di, 1, sum)
    p.hit = ilogit(z.di)
    pred.di = (p.hit>0.5)*1
    pred.di[p.hit==0.5] = data$deci[p.hit==0.5]
    p.ch = (1-p.hit)*(1-data$deci) + p.hit*data$deci
    
    # decision confidence
    message("decision confidence phase...")
    
    SV = abs(Z.di)
    Z.icf = b.icf[,,1]*SV + sweep(b.icf[,,2],2,data$kcf,"*")
    Q = ilogit(sweep(theta.icf,c(1,2),Z.icf,"+"))
    p.icf.mat = array(NaN,dim=c(3*nsmp,N,4))
    p.icf.mat[,,1] = 1 - Q[,,3]
    p.icf.mat[,,2] = Q[,,3] - Q[,,2]
    p.icf.mat[,,3] = Q[,,2] - Q[,,1]
    p.icf.mat[,,4] = Q[,,1]
    cmat = sapply(1:4, function(j){ data$icf==j })
    p.icf =apply(sweep(p.icf.mat,c(2,3),cmat,"*"), c(1,2), sum)
    
    sv = abs(z.di)
    z.icf = m.b.icf[,1]*sv + m.b.icf[,2]*data$kcf
    Q = ilogit(sweep(m.theta.icf,1,z.icf,"+"))
    p.icf.mat = array(NaN,dim=c(N,4))
    p.icf.mat[,1] = 1 - Q[,3]
    p.icf.mat[,2] = Q[,3] - Q[,2]
    p.icf.mat[,3] = Q[,2] - Q[,1]
    p.icf.mat[,4] = Q[,1]
    pred.icf = mapply(prediction, 1:N, MoreArgs=list(pmat=p.icf.mat,y=data$icf))
    w.pred.icf = calc.weighted.conf(p.icf.mat)
    
    res = list(pmat.di=array(c(1-p.hit,p.hit),dim=c(N,2)), pred.di=pred.di, pmat.icf=p.icf.mat, pred.icf=pred.icf, w.pred.icf=w.pred.icf, 
               loglik=list(di=log(p.di), icf=log(p.icf)))
    
    return(res)
}


basic.dimdl = function(data, model, pnames){

    ns = max(data$sidx)
    N = dim(data)[1]
    nsmp = dim(model[[1]])[1]
    b.di = array(NaN, dim=c(3*nsmp,N,3))
    b.icf = array(NaN, dim=c(3*nsmp,N,1))
    theta.icf = array(NaN, dim=c(3*nsmp,N,3))
    
    mdl.mcmc = rbind(model[[1]],model[[2]],model[[3]])
    for (i in 1:ns) {
        n = sum(data$sidx==i)
        b.di[,data$sidx==i,1] = array(rep(mdl.mcmc[,pnames==paste0("b0.di.p[",i,"]")],n), dim=c(3*nsmp,n))
        b.di[,data$sidx==i,2] = array(rep(mdl.mcmc[,pnames==paste0("b1.di.p[",i,"]")],n), dim=c(3*nsmp,n))
        b.di[,data$sidx==i,3] = array(rep(mdl.mcmc[,pnames==paste0("b2.di.p[",i,"]")],n), dim=c(3*nsmp,n))
        
        b.icf[,data$sidx==i,1] = array(rep(mdl.mcmc[,pnames==paste0("b1.icf.p[",i,"]")],n), dim=c(3*nsmp,n))

        a0 = array(rep(mdl.mcmc[,pnames==paste0("a0.icf.p[",i,"]")],n), dim=c(3*nsmp,n))
        a1 = array(rep(mdl.mcmc[,"a1.icf.mu"],n), dim=c(3*nsmp,n))
        theta.icf[,data$sidx==i,1] = a0-a1
        theta.icf[,data$sidx==i,2] = a0
        theta.icf[,data$sidx==i,3] = a0+a1
    }
    m.b.di = apply(b.di,c(2,3),median)
    m.b.icf = apply(b.icf,c(2,3),median)
    m.theta.icf = apply(theta.icf,c(2,3),median)
    
    # decision
    message("decision phase...")

    X = array(1,dim=c(3*nsmp,N,3))
    X[,,2] = rep(1,3*nsmp) %o% (21 - (data$open + data$ud + data$add))
    X[,,3] = sweep(X[,,2],2,abs(21 - (data$open + data$ud + data$add)),"*")
    Z.di = apply(X*b.di, c(1,2), sum)
    p.hit.mat = ilogit(Z.di)
    p.di = sweep((1-p.hit.mat),2,(1-data$deci),"*") + sweep(p.hit.mat,2,data$deci,"*")

    x = cbind(array(1,dim=c(N,1)), (21 - (data$open + data$ud + data$add)), (21 - (data$open + data$ud + data$add))*abs(21 - (data$open + data$ud + data$add)))
    z.di = apply(x*m.b.di, 1, sum)
    p.hit = ilogit(z.di)
    pred.di = (p.hit>0.5)*1
    pred.di[p.hit==0.5] = data$deci[p.hit==0.5]
    p.ch = (1-p.hit)*(1-data$deci) + p.hit*data$deci
    
    # decision confidence
    message("decision confidence phase...")
    
    SV = abs(Z.di)
    Z.icf = b.icf[,,1] * SV
    Q = ilogit(sweep(theta.icf,c(1,2),Z.icf,"+"))
    p.icf.mat = array(NaN,dim=c(3*nsmp,N,4))
    p.icf.mat[,,1] = 1 - Q[,,3]
    p.icf.mat[,,2] = Q[,,3] - Q[,,2]
    p.icf.mat[,,3] = Q[,,2] - Q[,,1]
    p.icf.mat[,,4] = Q[,,1]
    cmat = sapply(1:4, function(j){ data$icf==j })
    p.icf =apply(sweep(p.icf.mat,c(2,3),cmat,"*"), c(1,2), sum)
    
    sv = abs(z.di)
    z.icf = m.b.icf[,1] * sv
    Q = ilogit(sweep(m.theta.icf,1,z.icf,"+"))
    p.icf.mat = array(NaN,dim=c(N,4))
    p.icf.mat[,1] = 1 - Q[,3]
    p.icf.mat[,2] = Q[,3] - Q[,2]
    p.icf.mat[,3] = Q[,2] - Q[,1]
    p.icf.mat[,4] = Q[,1]
    pred.icf = mapply(prediction, 1:N, MoreArgs=list(pmat=p.icf.mat,y=data$icf))
    w.pred.icf = calc.weighted.conf(p.icf.mat)
    
    res = list(pmat.di=array(c(1-p.hit,p.hit),dim=c(N,2)), pred.di=pred.di, pmat.icf=p.icf.mat, pred.icf=pred.icf, w.pred.icf=w.pred.icf, 
                loglik=list(di=log(p.di), icf=log(p.icf)))

    return(res)
}


theme_custom = theme_classic()+
               theme(text=element_text(size=13,color="black",family="Arial"), axis.text=element_text(size=13,color="black"),
                     axis.title.x = element_text(margin=margin(t=5)), axis.title.y = element_text(margin=margin(r=6)),
                     legend.position = "none")


simulate.bhv = function(mdlres, simns, dat){
    set.seed(444)
    n = dim(dat)[1]
    nses = 20
    ntrl = 16

    ### parameter setting
    m.b0.dk = mdlres$dk.summaries$summaries["b.dk.mu[1]","Mean"]
    m.b1.dk = mdlres$dk.summaries$summaries["b.dk.mu[2]","Mean"]
    s.b0.dk = mdlres$dk.summaries$summaries["b.dk.sig[1]","Mean"]
    s.b1.dk = mdlres$dk.summaries$summaries["b.dk.sig[2]","Mean"]

    m.b1.kcf = mdlres$dk.summaries$summaries["b.kcf.mu[1]","Mean"]
    m.b2.kcf = mdlres$dk.summaries$summaries["b.kcf.mu[2]","Mean"]
    s.b1.kcf = mdlres$dk.summaries$summaries["b.kcf.sig[1]","Mean"]
    s.b2.kcf = mdlres$dk.summaries$summaries["b.kcf.sig[2]","Mean"]

    m.a0.kcf = mdlres$dk.summaries$summaries["a0.kcf.mu","Mean"]
    s.a0.kcf = mdlres$dk.summaries$summaries["a0.kcf.sig","Mean"]
    m.a1.kcf = mdlres$dk.summaries$summaries["a1.kcf.mu","Mean"]

    m.b0.di = mdlres$di.summaries$summaries["b.di.mu[1]","Mean"]
    m.b1.di = mdlres$di.summaries$summaries["b.di.mu[2]","Mean"]
    m.b2.di = mdlres$di.summaries$summaries["b.di.mu[3]","Mean"]
    s.b0.di = mdlres$di.summaries$summaries["b.di.sig[1]","Mean"]
    s.b1.di = mdlres$di.summaries$summaries["b.di.sig[2]","Mean"]
    s.b2.di = mdlres$di.summaries$summaries["b.di.sig[3]","Mean"]

    m.b1.icf = mdlres$di.summaries$summaries["b.icf.mu[1]","Mean"]
    m.b2.icf = mdlres$di.summaries$summaries["b.icf.mu[2]","Mean"]
    s.b1.icf = mdlres$di.summaries$summaries["b.icf.sig[1]","Mean"]
    s.b2.icf = mdlres$di.summaries$summaries["b.icf.sig[2]","Mean"]
    
    m.a0.icf = mdlres$di.summaries$summaries["a0.icf.mu","Mean"]
    s.a0.icf = mdlres$di.summaries$summaries["a0.icf.sig","Mean"]
    m.a1.icf = mdlres$di.summaries$summaries["a1.icf.mu","Mean"]

    b.dk = array(c(rnorm(simns, mean=m.b0.dk, sd=s.b0.dk), rnorm(simns, mean=m.b1.dk, sd=s.b1.dk)), dim=c(simns,2))
    b.kcf = array(c(rnorm(simns, mean=m.b1.kcf, sd=s.b1.kcf), rnorm(simns, mean=m.b2.kcf, sd=s.b2.kcf)), dim=c(simns,2))
    a0.kcf = rnorm(simns, mean=m.a0.kcf, sd=s.a0.kcf)

    b.di = array(c(rnorm(simns, mean=m.b0.di, sd=s.b0.di), rnorm(simns, mean=m.b1.di, sd=s.b1.di), rnorm(simns, mean=m.b2.di, sd=s.b2.di)), 
                dim=c(simns,3))
    b.icf = array(c(rnorm(simns, mean=m.b1.icf, sd=s.b1.icf), rnorm(simns, mean=m.b2.icf, sd=s.b2.icf)), dim=c(simns,2))
    e.kcf = runif(simns, min=0.5, max=1)
    g.kcf = sapply(1:simns,function(i){runif(1, min=0, max=e.kcf[i])})
    a0.icf = rnorm(simns, mean=m.a0.icf, sd=s.a0.icf)

    sidx = array(t(array(rep(1:simns,n),dim=c(simns,n))))
    ses = array(t(array(rep(1:nses,simns*ntrl),dim=c(nses*simns,ntrl))))
    trl = rep(1:16, simns*nses)

    ### deck estim
    dkest = c()
    for (i in 1:simns) {
        X = cbind(array(1,dim=c(n,1)),dat$evd.rea)
        z.dk = apply(sweep(X, 2, b.dk[i,], "*"),1,sum)

        dkest = c(dkest, sapply(1:n, function(l){rbinom(1,size=1,prob=ilogit(z.dk[l]))}))
    }

    ### deck conf
    vr = c()
    for (i in 1:simns) {
        for (s in 1:nses) {
        tmpdkest = dkest[sidx==i&ses==s]
        tmpdkest0 = (tmpdkest==0)*1
        tmpdkest1 = (tmpdkest==1)*1
        vr = c(vr, (sapply(1:16,function(j){nansum(tmpdkest0[1:j])})*tmpdkest0+sapply(1:16,function(j){nansum(tmpdkest1[1:j])})*tmpdkest1)/c(1:16))
        }
    }

    kcf = c()
    for (i in 1:simns) {
        X = cbind(abs(dat$evd.rea), vr[sidx==i])
        z.kcf = apply(sweep(X, 2, b.kcf[i,], "*"), 1, sum)
        p.kcf.mat = array(NaN,dim=c(n,4))
        p.kcf.mat[,1] = 1 - ilogit(z.kcf+a0.kcf[i]+m.a1.kcf)
        p.kcf.mat[,2] = ilogit(z.kcf+a0.kcf[i]+m.a1.kcf) - ilogit(z.kcf+a0.kcf[i])
        p.kcf.mat[,3] = ilogit(z.kcf+a0.kcf[i]) - ilogit(z.kcf+a0.kcf[i]-m.a1.kcf)
        p.kcf.mat[,4] = ilogit(z.kcf+a0.kcf[i]-m.a1.kcf)
        kcf = c(kcf, rcat(n,p.kcf.mat))
    }

    ### deci, deci conf
    sclkcf = (kcf-1)/3
    deci = c()
    icf = c()
    ud = rep(3.4, n)
    ud[dat$udc==1] = 2.5
    for (i in 1:simns) {
        w.cdk = e.kcf[i] - g.kcf[i]*(1-sclkcf[sidx==i])
        cdk = 4*(1-dkest[sidx==i]) + 6*dkest[sidx==i]
        udk = 4*dkest[sidx==i] + 6*(1-dkest[sidx==i])
        modadd = w.cdk*cdk + (1-w.cdk)*udk
        vh = 21 - (dat$open + ud + modadd)
        X = cbind(array(1,dim=c(n,1)), vh, vh*abs(vh))
        z.di = apply(sweep(X, 2, b.di[i,], "*"), 1, sum)

        deci = c(deci, sapply(1:n, function(l){rbinom(1,size=1,prob=ilogit(z.di[l]))}))

        sv = abs(z.di)
        z.icf = b.icf[i,1]*sv + b.icf[i,2]*sclkcf[sidx==i]
        p.icf.mat = array(NaN,dim=c(n,4))
        p.icf.mat[,1] = 1 - ilogit(z.icf+a0.icf[i]+m.a1.icf)
        p.icf.mat[,2] = ilogit(z.icf+a0.icf[i]+m.a1.icf) - ilogit(z.icf+a0.icf[i])
        p.icf.mat[,3] = ilogit(z.icf+a0.icf[i]) - ilogit(z.icf+a0.icf[i]-m.a1.icf)
        p.icf.mat[,4] = ilogit(z.icf+a0.icf[i]-m.a1.icf)
        icf = c(icf, rcat(n,p.icf.mat))
    }

    simdat = data.frame(sidx=sidx, ses=ses, dkest=dkest, dktrue=rep(dat$dktrue, simns), evd.rea=rep(dat$evd.rea, simns), vr=vr, kcf=kcf,
                        open=dat$open, udc=rep(dat$udc, simns), deci=deci, icf=icf, trl=trl)
    simparam = list(m.b0.dk=m.b0.dk, m.b1.dk=m.b1.dk, m.b1.kcf=m.b1.kcf, m.b2.kcf=m.b2.kcf, m.a0.kcf=m.a0.kcf, m.a1.kcf=m.a1.kcf,  
                    m.b0.di=m.b0.di, m.b1.di=m.b1.di, m.b2.di=m.b2.di, m.b1.icf=m.b1.icf, m.b2.icf=m.b2.icf, m.a0.icf=m.a0.icf, m.a1.icf=m.a1.icf,   
                    b.dk=b.dk, b.kcf=b.kcf, a0.kcf=a0.kcf,
                    b.di=b.di, g.kcf=g.kcf, e.kcf=e.kcf, b.icf=b.icf, a0.icf=a0.kcf)
    simdat = list(dat=simdat, prm=simparam)

    return(simdat)
}