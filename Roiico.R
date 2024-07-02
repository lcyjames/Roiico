list.of.packages <- c("MASS","BPST","splines2")
new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
if(length(new.packages)) install.packages(new.packages)
library(BPST)
library(MASS)
library(splines2)

Bernstein<-function(Y, V.est, Tr.est, d.est, r=1, Z, lambda = 10^seq(-6,6,0.5)){
  n <- nrow(Y)
  Bfull.est <- basis(V.est,Tr.est,d.est,r,Z)
  B <- Bfull.est$B
  ind.inside <- Bfull.est$Ind.inside
  Q2 <- Bfull.est$Q2
  K <- Bfull.est$K
  Y <- matrix(Y[,ind.inside],nrow=n)
  lambda <- as.matrix(lambda)
  t.area = Bfull.est$tria.all 
  
  this.call <- match.call()
  n <- nrow(Y)
  npix <- ncol(Y)
  J <- ncol(Q2)
  
  W <- as.matrix(B%*%Q2)
  WW <- crossprod(W,W)
  rhs <- crossprod(W,t(Y))
  D <- crossprod(t(crossprod(Q2,as.matrix(K))),Q2)
  D <- as.matrix(D)
  
  flag <- (rankMatrix(WW)<J)
  if(!flag){
    Ainv <- chol(WW,pivot=TRUE)
    A <- solve(t(Ainv))
    ADA <- A%*%D%*%t(A)
    eigs <- eigen(ADA)
    Cval <- eigs$values
  }
  
  nl <- length(lambda)
  
  gcv_all <- sapply(lambda,FUN=function(Lam){  
    Dlam <- Lam*D
    lhs <- WW+Dlam
    lhs.inv <- chol2inv(chol(lhs));
    theta <- crossprod(t(lhs.inv),rhs)
    gamma <- crossprod(t(Q2),theta)
    Yhat <- crossprod(t(W),theta)
    res <- t(Y)-Yhat
    sse <- apply(res^2,2,sum)
    if(!flag){
      df <- sum(1/(1+Cval*Lam))
    }
    if(flag){
      Hmtx <- crossprod(t(crossprod(t(W),lhs.inv)),t(W))
      df <- sum(diag(Hmtx))
    }
    gcv <- npix*sse/(npix-df)^2
  })
  gcv_all <- matrix(gcv_all,nrow=n)
  lam.ind <- apply(gcv_all,1,which.min)
  lambdac <- lambda[lam.ind]
  
  theta <- c()
  gamma <- c()
  Yhat <- c()
  df <- c()
  for (i in 1:n){
    lamc.tmp <- lambdac[i]
    Dlam <- lamc.tmp*D
    lhs <- WW+Dlam
    lhs.inv <- chol2inv(chol(lhs));
    rhs.tmp <- as.matrix(rhs[,i],ncol=1)
    theta.tmp <- crossprod(t(lhs.inv),rhs.tmp)
    theta <- cbind(theta,theta.tmp)
    gamma.tmp <- crossprod(t(Q2),theta.tmp) 
    gamma <- cbind(gamma,gamma.tmp)
    Yhat.tmp <- crossprod(t(W),theta.tmp)
    Yhat <- cbind(Yhat,Yhat.tmp)
    if(!flag){
      df.tmp <- sum(1/(1+Cval*lamc.tmp))
    }
    if(flag){
      Hmtx <- crossprod(t(crossprod(t(W),lhs.inv)),t(W))
      df.tmp <- sum(diag(Hmtx))
    }
    df <- c(df,df.tmp)
  }
  
  return(list(B, gamma, lambdac, t.area))
}
lik_gaml <-function(par, ind, gam, xi.mat, Z.beta.d, Exi, wLR, pen, EABij, J){
  gaml        <-gam
  gaml[ind]   <-par
  N.gaml      <-norm(gaml[ind],type = "2")
  xi.gam      <-c(xi.mat%*%gaml)
  xibar_den   <-c(c(exp(Z.beta.d+xi.gam)*Exi)%*%wLR)
  return(-sum(EABij*(replicate(J,c(Z.beta.d+xi.gam))-t(replicate(nrow(xi.mat),log(xibar_den)))))+pen*N.gaml)
}
RoiicoSIM<-function(seed = NA, n, deg = 2, a, beta, gamma, rho, VT, pattern){
  if(is.numeric(seed)){set.seed(seed)}
  pw     <-0.4
  pb     <-0.1
  Sigma  <-0.1
  
  block_size <-factorial(deg+2)/factorial(deg)/factorial(2)
  num_blocks <-nrow(VT$Tr)
  total_size <-num_blocks*block_size
  num.select<-floor(num_blocks*a)   #number of selected triangles
  select    <-pattern[1:num.select] #the index of the selected triangles
  
  Bfull  <-basis(V = VT$V,Tr = VT$Tr, d = deg, r = 1, Z = VT$Z_grid)
  Basis  <-as.matrix(Bfull$B)
  
  matrix_blockdiag <- matrix(pb, nrow = total_size, ncol = total_size)
  for (i in 1:num_blocks) {
    start_row <- (i-1)*block_size+1
    end_row   <- i*block_size
    matrix_blockdiag[start_row:end_row, start_row:end_row]<-pw
  }
  diag(matrix_blockdiag)<-1
  aijk <- mvrnorm(n = n,mu = rep(0,total_size), Sigma = matrix_blockdiag)
  
  Y_grid <- mu_grid <- matrix(NA,nrow=n,ncol=nrow(VT$Z_grid))
  if(n==1){
    mu_grid<-as.numeric(aijk%*%t(Basis))
    Y_grid <-mu_grid+rnorm(n=nrow(VT$Z_grid),0,Sigma)
  }else{
    mu_grid<- as.matrix(aijk%*%t(Basis))
    Y_grid <- mu_grid+matrix(rnorm(n=n*nrow(VT$Z_grid),0,Sigma),nrow = n,ncol = nrow(VT$Z_grid))
  }
  
  gijk         <-matrix(0, nrow=num_blocks,ncol = block_size, byrow=T)
  gijk[select,]<-t(replicate(length(select),c(0.1, 0.2, 0.3, 0.5, 0.6, 0.4)))
  gamma        <-as.matrix(c(t(gijk)))
  xi           <-aijk%*%t(Basis)%*%Basis
  xi.norm      <-as.matrix(scale(xi))
  
  Z1       <-rbinom(n = n,size = 1,prob = 0.5)
  Z2       <-rnorm(n = n, mean = 0, sd = 1)
  regressor<-exp(cbind(Z1,Z2)%*%beta+c(xi.norm%*%gamma))
  Ui       <-runif(n = n,min = 0,max = 1)
  if(rho==0)  {ti<-(-log(Ui)/regressor/0.25)^(1/2)} 
  if(rho==0.5){ti<-((Ui^{-0.5}-1)/regressor/0.125)^(1/2)}
  if(rho==1)  {ti<-((Ui^{-1}-1)/regressor/0.25)^(1/2)}
  
  checkup.times<-lapply(1:n,function(w){x<-cumsum(runif(100,min = 0,max = 1))
  x<-c(0,x[which(x<3+2)]); return(x)})
  which.interval<-sapply(1:n,function(x){findInterval(ti[x],checkup.times[[x]])})
  
  Li<-ifelse(ti>=3, 3, sapply(1:n,function(y){checkup.times[[y]][which.interval[y]]}))
  Ri<-ifelse(ti>=3, Inf, pmin(3,sapply(1:n,function(z){checkup.times[[z]][which.interval[z]+1]})))
  DeltaL<-ifelse(Li==0,1,0)
  DeltaI<-ifelse(!(Li==0)&!is.infinite(Ri),1,0)
  
  store<-as.data.frame(cbind(c(1:n),Li,Ri,DeltaL,DeltaI,cbind(Z1,Z2)))
  names(store)<-c("id","Li","Ri","DL","DI","Z1","Z2")
  
  return(list(Y.data=Y_grid, surv.data=store, selected.true=sort(select)))
}
RoiicoEST<-function(Yi, Zi, Li, Ri, DL, DI, rho, VT, deg1 = 2, deg2 = 3, J = 7, tolerance = 10^{-4}, lambda.grid = 10^seq(6,-6,-0.1), TRACE = FALSE){
  n     <- nrow(Yi)
  numTr <- nrow(VT$Tr)
  Z.mat <- as.matrix(Zi)
  P     <- ncol(Z.mat)
  est   <- Bernstein(Y=Yi, V.est=VT$V, Tr.est=VT$Tr, d.est=deg1, Z = VT$Z_grid)
  basis <- as.matrix(est[[1]])
  scores<- as.matrix(t(est[[2]]))
  xi.mat<- scale(scores%*%t(basis)%*%basis)
  Q     <- ncol(xi.mat)
  
  grp.id<- rep(1:numTr,each=choose(deg1+2,2))
  sKl   <- rep(sqrt(choose(deg1+2,2)),numTr)
  
  part   <-J-deg2+1
  LRi    <-sort(unique(c(Li,Ri)[is.finite(c(Li,Ri))]))
  DR     <-1-DL-DI
  if(J<deg2){print("ERROR: J cannot be less than the degree of polynomial splines (deg).");break}
  if(J==deg2){knots<-NULL}else{knots <- quantile(LRi,probs = seq(1/part,(part-1)/part,1/part),na.rm = T)}  #adaptive
  ibsMat   <-iSpline(LRi, knots = knots, degree = deg2, intercept = FALSE)
  basis_Li <-ibsMat[match(Li,LRi),]
  basis_Ri <-ibsMat[match(Ri,LRi),]; basis_Ri[is.na(basis_Ri)]<-1
  
  if(rho<0){print("ERROR: rho cannot be negative.");break}
  
  #search for lambda
  Store.H<-c()
  tri.selected.lambda <- vector(mode = "list", length = length(lambda.grid))
  KS<-1; tot<-0
  for(H in 1:length(lambda.grid)){
    conv<-1
    while(conv!=0){
      lambda <-lambda.grid[H]
      if(KS==1){omega.d <- rep(0.1,J); beta.d  <- rep(0,P); gam.d   <- rep(0,Q); KS<-0}
      
      I<-0
      tick<-0
      obs.like2<-0
      complete_cycle<-1
      active_update <-0  #whether there is a newly added member in the set
      active_set_old<-c()
      repeat{
        I<-I+1
        omega.c<- omega.d; beta.c <- beta.d; gam.c  <- gam.d
        
        Lambda.Li<-c(basis_Li%*%omega.d)
        Lambda.Ri<-c(basis_Ri%*%omega.d)
        Z.beta.d <-c(Z.mat%*%beta.d)
        xi.gam.d <-c(xi.mat%*%gam.d)
        
        SLi    <- Lambda.Li*exp(Z.beta.d+xi.gam.d)
        SRi    <- Lambda.Ri*exp(Z.beta.d+xi.gam.d)
        mu_ij  <- as.matrix(t(replicate(n,omega.d))*(replicate(J,DL)*basis_Ri+replicate(J,DI)*basis_Li)*(replicate(J,exp(Z.beta.d+xi.gam.d))))
        eta_ij <- as.matrix(t(replicate(n,omega.d))*(replicate(J,DI)*(basis_Ri-basis_Li)+replicate(J,DR)*basis_Li)*(replicate(J,exp(Z.beta.d+xi.gam.d))))
        
        if(rho==0){
          e_GLi<- exp(-SLi)
          e_GRi<- ifelse(Ri==Inf, 0, exp(-SRi))
          GdSLi<- rep(1,n)
          GdSRi<- rep(1,n)
        }else{
          e_GLi<- (1+rho*SLi)^(-1/rho)
          e_GRi<- ifelse(Ri==Inf, 0,(1+rho*SRi)^(-1/rho))
          GdSLi<- (1+rho*SLi)^(-1)
          GdSRi<- ifelse(Ri==Inf, 0, (1+rho*SRi)^(-1))
        }
        e_G_diff<-e_GLi-e_GRi
        
        Exi <- numeric(n)
        EAij <- EBij <- matrix(0, n, J) 
        for(i in c(1:n)){
          if(DL[i]==1){EAij[i,]<-mu_ij[i,]/e_G_diff[i]}
          if(DI[i]==1){EBij[i,]<-eta_ij[i,]*e_GLi[i]*GdSLi[i]/e_G_diff[i]}
        }
        Exi<-(e_GLi*GdSLi-e_GRi*GdSRi)/e_G_diff
        
        obs.like1<-sum(log(ifelse(DL, 1-e_GRi, ifelse(DI, e_G_diff, e_GLi))))
        like.diff<-obs.like1-obs.like2
        
        wLR     <-replicate(J,DR)*basis_Li+replicate(J,(1-DR))*basis_Ri
        EABij   <-EAij+EBij
        
        Zbar_den<-c(t(as.matrix(exp(Z.beta.d+xi.gam.d)*Exi))%*%wLR)
        Zbar    <-t((Z.mat)*replicate(P,exp(Z.beta.d+xi.gam.d)*Exi))%*%wLR/t(replicate(P,Zbar_den))
        
        score.beta<-numeric(P)
        for(o in 1:P){score.beta[o]<-sum(EABij*(replicate(J,Z.mat[,o])-t(replicate(n,Zbar[o,]))))}
        ZZmat<-array(0, dim = c(P,P,J))
        for(j in 1:J){
          ZZmat[,,j]<- t(Z.mat)%*%(replicate(P,c(exp(Z.beta.d+xi.gam.d)*Exi*(wLR[,j])))*Z.mat)/matrix(replicate(P*P,Zbar_den[j]),nrow=P)
        }
        hess.beta<-matrix(0,nrow = P,ncol = P)
        for(o1 in 1:P){for(o2 in 1:P){
          hess.beta[o1,o2]<-sum(-EABij*t(replicate(n,(ZZmat[o1,o2,]-Zbar[o1,]*Zbar[o2,]))))
        }}
        beta.d   <-beta.d-as.numeric(ginv(hess.beta)%*%score.beta)
        Z.beta.d <-c(Z.mat%*%beta.d)
        
        if(complete_cycle==1){
          active_set<-c()
          for(r in 1:numTr){
            ind.r     <-which(grp.id==r)
            idlen     <-length(ind.r)
            sub.xi.mat<-xi.mat[,ind.r]
            
            gam_0r       <-gam.d
            gam_0r[ind.r]<-0
            xi.gam_0r    <-c(xi.mat%*%gam_0r)
            expbzxi_0r   <-c(exp(Z.beta.d+xi.gam_0r)*Exi)
            xi.den_0r    <-as.numeric(expbzxi_0r%*%wLR)
            xibar_0r     <-t((sub.xi.mat)*replicate(idlen,expbzxi_0r))%*%wLR/t(replicate(idlen,xi.den_0r))
            grad.lc_0r   <-numeric(idlen)
            for(rr in 1:idlen){grad.lc_0r[rr]<-sum(EABij*(replicate(J,sub.xi.mat[,rr])-t(replicate(n,xibar_0r[rr,]))))}
            
            if(norm(grad.lc_0r,type="2")<=lambda*sKl[r]){gam.d[ind.r]<-0}else{
              active_set<-c(active_set,r)
              fit<-optim(par = gam.d[ind.r],fn = lik_gaml,gr = NULL,method = "BFGS", 
                         ind=ind.r,gam=gam.d,xi.mat=xi.mat,Z.beta.d=Z.beta.d,Exi=Exi,wLR=wLR,pen=lambda*sKl[r],EABij=EABij,J=J,
                         control=list(reltol=10^{-6}))
              gam.d[ind.r]<-fit$par
            }
            if(length(active_set)==length(active_set_old)){active_update<-0}else{active_update<-1}
          }
          complete_cycle<-0
          active_set_old<-active_set
        }else if(complete_cycle==0){
          for(r in active_set){
            ind.r     <-which(grp.id==r)
            fit<-optim(par = gam.d[ind.r],fn = lik_gaml,gr = NULL,method = "BFGS", 
                       ind=ind.r,gam=gam.d,xi.mat=xi.mat,Z.beta.d=Z.beta.d,Exi=Exi,wLR=wLR,pen=lambda*sKl[r],EABij=EABij,J=J,
                       control=list(reltol=10^{-6}))
            gam.d[ind.r]<-fit$par
          }
        }
        
        Z.beta.d <-c(Z.mat%*%beta.d)
        xi.gam.d <-c(xi.mat%*%gam.d)
        ome1     <- apply(EABij,2,sum)
        ome2     <- apply(replicate(J,Exi*exp(Z.beta.d+xi.gam.d))*wLR,2,sum)
        omega.d  <- ome1/ome2
        
        obs.like2<-obs.like1
        
        dist<-max(abs(c(omega.c-omega.d, beta.c-beta.d, gam.c-gam.d)))
        DFg <-sum(abs(gam.d)>0)
        AIC <- -2*obs.like1+2*(P+DFg)
        
        if((dist<tolerance)&(active_update==1)){complete_cycle<-1}
        if((dist<tolerance)&(active_update==0)){
          tot<-tot+1
          conv<-0
          selected<-apply(abs(matrix(gam.d,nrow=numTr,byrow=T))>0,1,prod)
          tri.selected.lambda[[H]] <- c(lambda[H], selected)
          if(TRACE==TRUE){print(round(c(H,lambda,conv,DFg,obs.like2,AIC),3))}
          Store.H <-rbind(Store.H,c(selected, lambda, obs.like1, AIC, DFg))
          break
        }
      }
    }
  }
  
  Store.H.df       <-as.data.frame(Store.H)
  names(Store.H.df)<-c(paste0("tri",1:numTr),"lambda","obslike","AIC","DF")
  opt.Store.AIC    <-Store.H.df[which.min(Store.H.df$AIC),]
  opt.select.AIC   <-which(opt.Store.AIC[1:numTr]==1)
  
  return(list(selected=opt.select.AIC, AIC.lambda=Store.H.df$AIC, opt.lambda=Store.H.df$lambda, tri.selected.lambda=tri.selected.lambda, lambda=lambda.grid))
}





