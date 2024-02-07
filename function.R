#### back transformation from L^2 to B^2 ####
clr2density <- function(z, z_step, clr)
{
  if(is.fd(clr))
    return(exp(eval.fd(z,clr))/trapzc(z_step,exp(eval.fd(z,clr))))
  if(!is.fd(clr))
    return(exp(clr)/trapzc(z_step,exp(clr)))
}

#### basis with zero integral ####
ZsplineBasis = function(knots,k)
{
  library(fda)
  r = length(knots)
  lambda_index = c(0:(r-1)) 
  g = lambda_index[length(lambda_index) - 1]
  
  N = g+(k-1)+1
  
  lambda = c(rep(min(knots),k-1),knots,rep(max(knots),k-1))
  division = seq(min(lambda), max(lambda), length = 1000) 
  
  # standard B-spline basis; collocation matrix := C
  splajn.basis = create.bspline.basis(range(knots),nbasis = N , norder = k, breaks = knots)
  C = eval.basis(division, splajn.basis)
  
  # Matrix D
  diffrence = lambda[(1+k):(r+2*(k-1))] - lambda[(1:(r+k-2))]
  D = (k)*diag(1/diffrence)
  
  # Matrix L
  L = array(0, c(N,N-1))
  L[1,1]=1
  L[N,N-1]=-1
  
  for (j in (2:(N-1))){
    L[j,j-1] = (-1)
    L[j,j] = 1
  }
  
  # Spline0 basis: collocation matrix C0
  C0 = C%*%D%*%L
  
  # Matrix M - function for computing integral
  SLP=function(step, c){
    integral = step*(0.5*c[1]+sum(c[2:(length(c)-1)]) +0.5*c[length(c)])
    return (integral)
  }
  
  Division = seq(min(lambda), max(lambda),  length = 10000);step=diff(Division[1:2])
  CC = eval.basis(Division, splajn.basis)
  
  CC0 = CC%*%D%*%L
  
  M=array(0, c(N-1,N-1))
  for (i in 1:(N-1)){
    for (j in 1:(N-1)){
      nozero = c()
      multi = CC0[,i]*CC0[,j]
      for (m in 1:length(Division)){
        if (multi[m] != 0) {nozero[m] = multi[m]}
      }
      M[i,j]=SLP(step, multi)
    }
  }
  return(list(C0 = C0, M0 = M, K = L, D = D))
}

#### smoothing spline ####
SmoothingSpline01 = function(knots,t,f,w,k,der,alfa,ch){
  # INPUT:  knots= sequence of knots of a spline,
  #         t = points of approximation,
  #         f = function values at points of approximation,
  #         w = weighting coefficients,
  #         k = order of spline, degree = k-1
  #         der = derivation
  #         alfa = smoothing parametr
  #         ch ={1,2};  ch=1: functional with (1-alfa) and alfa
  #                     ch=2: funkcional with alfa   
  
  r = length(knots) 
  library(splines)
  
  Final_length = 2*(k-1) + r 
  y = c() 
  for (i in 1:(Final_length)){
    if (i <= k-1){y[i] = knots[1]}
    if ((i > k-1) && (i <= r + k-1)){y[i] = knots[i-(k-1)]}
    if (i > r+ k-1){y[i] = knots[r]}
  }
  
  # Collocation matrix K
  K = splineDesign(y, t_mid, k, outer.ok = TRUE)
  
  # Diag matrix with weights
  W = diag(w)
  
  # Collocation matrix C 
  division = seq(min(y), max(y),   length = 1000)    
  lambda = c(0:(r-1)) 
  g = lambda[length(lambda) - 1]
  # Dimension(space of splines)
  N = g+(k-1)+1
  C = array(0, c(length(division),N))
  l = c()
  for(i in (1:N)){
    for (j in 1:(k+1)){
      l[j] = y[i+j-1]
    }
    C[ ,i] = splineDesign(l, division, k, outer.ok = TRUE) 
  }
  # Verification of full column rank of collocation matrix K
  if (length(t_mid) <= N) stop ('length(t) must be higher then Dimension(space of splines)')
  if (qr(K)$rank != N) stop ('Collocaton matrix does not have full column rank.')
  
  # Matrix S
  S = array(0)
  if (der == 0){
    S = diag(1, c(N,N))
  }
  if (der > 0){
    i=der
    while (i>0){
      D = array(0)
      diffrence = y[(1+k):(N+k-i)] - y[(1+i):(N)]
      D = (k-i)*diag(1/diffrence)
      L = array(0, c(N-i,N-i+1))
      for (j in (1:(N-i))){
        L[j,j] = (-1)
        L[j,j+1] = 1
      }
      if (i==der){
        S = D%*%L
      } else{
        S = S%*%D%*%L
      }
      i=i-1
    }
  }
  
  # Matrix M -  order of spline = k-der
  kk = k-der 
  
  # Matrix M - augment knot sequence
  final_length = 2*(kk-1) + r                                         
  Y = c()
  for (i in 1:final_length){
    if (i <= (kk-1)){Y[i] = knots[1]}
    if ((i > kk-1) && (i <= r + kk-1)){Y[i] = knots[i-(kk-1)]}
    if (i > (r+(kk-1))){Y[i] = knots[r]}
  }  
  
  Division = seq(min(Y), max(Y),  length = 10000)    
  Lambda = c(0:(r-1)) 
  G = Lambda[length(Lambda) - 1]
  
  # Matrix M - spline space dimension
  NN = G+(kk-1)+1
  
  # Matrix M - collocation matrix KK
  CC = splineDesign(Y, Division, kk, outer.ok=TRUE)
  # Matrix M - function for computing integral
  SLP=function(step, c){
    integral = step*(0.5*c[1]+sum(c[2:(length(c)-1)]) +0.5*c[length(c)])
    return (integral)
  }
  
  step=diff(Division[1:2])
  
  # Matrix M
  M=array(0, c(NN,NN))
  for (i in 1:NN){
    for (j in 1:NN){
      nozero = c()
      multi = CC[,i]*CC[,j]
      for (m in 1:length(Division)){
        if (multi[m] != 0) {nozero[m] = multi[m]}
      }
      M[i,j]=SLP(step, multi)
    }
  }
  
  # Matrix D
  diffrence = y[(1+k):(r+2*(k-1))] - y[(1:(r+k-2))]
  D = (k)*diag(1/diffrence)
  
  # Matrix K
  KK = array(0, c(N,N-1))
  KK[1,1]=1
  KK[N,N-1]=-1
  
  for (j in (2:(N-1))){
    KK[j,j-1] = (-1)
    KK[j,j] = 1
  }
  KK
  
  # Matrix U
  U = S%*%D%*%KK
  
  # Matrix G, vector g
  if (ch==1){GG = t(U)%*%((1-alfa)*M)%*%U + alfa * t(KK)%*%t(D)%*%t(K)%*%W%*%K%*%D%*%KK}
  if (ch==2){GG = t(U)%*%M%*%U +            alfa * t(KK)%*%t(D)%*%t(K)%*%W%*%K%*%D%*%KK}
  
  gg = alfa*t(KK)%*%t(D)%*%t(K)%*%W%*%f
  
  # vector of B-spline coefficients := z
  z = solve(GG)%*%gg
  
  # B-spline basis in L20
  Bbasis = C%*%D%*%KK
  
  # Resulting spline
  spline0 = (C%*%D%*%KK)%*%z
  
  matplot(division,spline0, type="l",las=1,xlab="t",ylab="",col="darkblue",lwd=2,
          ylim = c(min(c(min(f),min(spline0))),max(c(max(f),max(spline0)))),
          main = paste("Spline of order k =",k))
  matpoints(t_mid,f, pch = 8,col=1)
  abline(h=0,col="red",lty=2,lwd=1)
  
  integral = SLP(step,spline0)
  
  if (ch==1) {J = (1-alfa)*t(z)%*%t(U)%*%M%*%U%*%z + alfa*t(f-K%*%D%*%KK%*%z)%*%W%*%(f-K%*%D%*%KK%*%z)}
  if (ch==2) {J =          t(z)%*%t(U)%*%M%*%U%*%z + alfa*t(f-K%*%D%*%KK%*%z)%*%W%*%(f-K%*%D%*%KK%*%z)}
  return(list(J=J,z=z,spline0=spline0))
}

#### backwards pivot coordinates ####
bpc1 <- function(x){
  x.bpc <- matrix(NA, nrow = nrow(x), ncol = ncol(x) - 1)
  D <- ncol(x)
  for (i in 1:ncol(x.bpc)) {
    x.bpc[, i] <- sqrt(i/(i+1))*log((x[, i+1])/apply(as.matrix(x[,1:i]), 1, gm))
  }
  return(x.bpc)
}

#### isotemporal substitution ####
is=function(from,to,bp,area){
  mat=matrix(nrow = 100,ncol = 1000)
  y_mat=matrix(nrow=10,ncol=100)
  t.fine=seq(from=from,to=to,l=1000)
  breakpoint=t.fine[min(which(t.fine>=bp))]
  dist1=breakpoint-from
  dist2=to-breakpoint
  for (j in 1:10) {   #pre vymenu intervalov
    for (i in 1:length(area)){   
      a1=1-area[i]
      w1=a1/dist2
      w2=area[i]/dist1
      mat[i,]=c(rep(w1,(j-1)*100),rep(w2,100),rep(w1,(10-j)*100))
    }
    list(mat=mat)
    W=mat   
    W_clr_trans=cenLR((W))$x.clr
    list(W=W)
    est_mean_spline=matrix(splines_mean,nrow=1000,ncol=100)
    f_star=est_mean_spline+t(W_clr_trans)
    
    s_scalar=t(estimate.l)%*%f_star
    scalar=s_scalar*(1/1000)
    scalar2=(t(estimate.l)%*%mean_spline)*(1/1000)
    BETA0_new=BETA0-scalar2
    y_final=as.vector(BETA0_new)+scalar
    matplot(t.fine,t(W),xlab="log(PA intensity)",ylab="PDF g",
            type="l",col=magma(100),
            lty=1,ylim = c(0,1.2))
    y_mat[j,]=y_final
  }
  return(list(y_mat=y_mat))
}
