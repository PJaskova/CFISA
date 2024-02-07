setwd("~/Desktop/github")

# packages
install.packages("zCompositions")
pack=c("zCompositions","robCompositions","compositions","fda","scatterplot3d","chemometrics",
       "carData","car","depthTools","boot","Epi","viridis","boot")
lapply(pack, require, character.only = TRUE,quietly =TRUE)
rm(pack)

source("function.R")

# loading data
density=read.csv("hustoty.csv")
clr_density = cenLR((density))$x.clr #clr transformation

# response
x=read.csv2("response.csv",header = TRUE,dec = ",",sep = ";")
resp=x$FM_PERCENT
l.response=car::logit(resp)  #logit-transformed body fat percentage

# covariates
covar=read.csv("covariates.csv")
cov2=cbind(covar$SLEEP_ADJUSTED24,covar$SB_ADJUSTED24)
colnames(cov2)=c("Sleep","SB")

age_parnition=seq(log10(36/1000),0.5,length=23)  #boundaries

age_cut = c()    #intervals
for(i in 1:(length(age_parnition)-1)) age_cut[i] = paste("[",age_parnition[i],",",age_parnition[i+1], ")",sep="")
n_class = length(age_parnition)-1

t_mid = c() #centres of intervals
for (i in 1:n_class){t_mid[i]=(age_parnition[i+1]+age_parnition[i])/2}
t_mid  

# setting of smoothing
knots=c(-1.4436975,-0.75,-0.5,  -0.3,0,   0.45582506)  # chosen sequence of knots
w = rep(1,22)                  # weights
k = 3                          # degree of spline, k=3 cubic spline
der = 1                        # first derivation continuous (derivation in functional)
alfa = 0.9                     # smoothing parameter
ch = 1                         # functional 1-alfa
f=t(clr_density)
n=ncol(f)

#z-coefficients
z_coef = NULL
spline_coef=matrix(nrow = 1000, ncol=n)
J = c()
for (i in 1:n){
  spline = SmoothingSpline01(knots=knots, t=t_mid, f=f[,i], w=w, k=k, der=der, alfa=alfa, ch=ch)
  abline(v=knots,col="gray",lty=2)
  z_coef = cbind(z_coef,spline[[2]])     #z-coefficients
  J[i] = spline[[1]]   
  spline_coef[,i]=SmoothingSpline01(knots=knots, t=t_mid, f=f[,i], w=w, k=k, der=der, alfa=alfa, ch=ch)[[3]]
}
sum(J)

# ZB-spline basis
Z = ZsplineBasis(knots=knots,k) 
Zco= Z$C0

data.l=Zco%*%(z_coef)

splines_mean=apply(spline_coef,1,mean)
t.fine = seq(knots[1],max(knots),l=1000) 
t.step=diff(t.fine[1:2])  

data1 = NULL
for (i in 1:n)
{
  data1 = cbind(data1, clr2density(t.fine, t.step, data.l[,i]))
}

par(mfrow=c(1,2))
matplot(t.fine,spline_coef,
        lty=1, type="l",las=1,
        ylab="clr(PDFs)",xlab="log(PA intensity)",col="darkblue")
abline(h=0,col="red",lty=1)

matplot(t.fine,data1, 
        lty=1, type="l",
        ylab="PDFs",xlab="log(PA intensity)", col="darkblue")
dev.off()

####  B-spline representation ####
# expression of compositional splines using B-spline basis: b=DKz
basisB = create.bspline.basis(range(knots),nbasis = dim(z_coef)[1]+1,norder=k,breaks=knots)
b_coef =Z$D%*%Z$K%*%(z_coef)    # b-coefficients
data.fd = fd(b_coef, basisB)    # defined functional data object from b-spline coefficients

#### SFPCA ####
ncom=2     # setting the number of components
pcafd = pca.fd(data.fd,ncom, centerfns = T)   #SFPCA
summary(pcafd)
attach(pcafd)
scores=pcafd$scores    #scores

#variation of estimated function
var_coef=var(t(b_coef))
bbasis=fd(diag(1,nrow=nrow(var_coef),ncol=ncol(var_coef)),basisB)
B_spline=eval.fd(t.fine,bbasis)
var_beta1=B_spline%*%var_coef%*%t(B_spline)
var_beta_1_new=var_beta1/n   

#### regression ####
X = as.matrix(cbind(rep(1,n),scores))

# (ordinary) least squares
LS = solve(t(X)%*%X)%*%t(X)%*%as.vector(l.response)
LS

# Estimate beta0 and beta1
BETA0=LS[1,]
BETA0

BETA1=LS[-1,]
BETA1

estimate.l = 0
estimate.b = c()
for (j in 1:ncom)
{
  estimate.l  =  estimate.l + BETA1[j]*eval.fd(t.fine, pcafd$harmonics[j,])
}
estimate.b = clr2density(t.fine, t.step, estimate.l)

matplot(t.fine, (estimate.l),type="l",col="darkblue",lty=1,las=1, ylab="clr(estimate of beta1)",xlab="log(PA intensity)",
        ylim = c(-0.45,0.45))


# confidence bands
LB=estimate.l-2*sqrt(diag(var_beta_1_new))
UB=estimate.l+2*sqrt(diag(var_beta_1_new))
y_mat=cbind(estimate.l,LB,UB)
matshade(x=t.fine,y=y_mat)

#### ISOTEMPORAL SUBSTITUTION ####    
mean_spline=matrix(splines_mean,nrow=1000,ncol=1)

## INPUT
from=knots[1]
to=max(knots)
bp=t.fine[100]
area2=seq(from=1/10,to=2/10,length=100)

par(mfrow=c(1,2),mar=c(4,4,2,2))
y_mat_final=is(knots[1],max(knots),bp,area2)$y_mat
y_mat_final_inv=inv.logit(y_mat_final)*100
dev.off()

plot(area2,y_mat_final_inv[1,],type="l",xlim=c(0.1,0.25),ylim=range(y_mat_final_inv),
     xlab = "Values of coefficient K (divided by 10)",ylab = "Body fat percentage")
lines(area2,y_mat_final_inv[2,],col="red")
lines(area2,y_mat_final_inv[3,],col="green")
lines(area2,y_mat_final_inv[4,],col="blue")
lines(area2,y_mat_final_inv[5,],col="lightsalmon")
lines(area2,y_mat_final_inv[6,],col="#660066")
lines(area2,y_mat_final_inv[7,],col="brown")
lines(area2,y_mat_final_inv[8,],col="grey")
lines(area2,y_mat_final_inv[9,],col="violet")
lines(area2,y_mat_final_inv[10,],col="gold")
legend("right", legend=c("36 - 56", "56 - 86","86 - 133","133 - 207","207 - 320",
                         "320 - 496","496 - 768","768 - 1190","1190 - 1844","1844 - 2856"),title="Intervals (mg)",
       col=c("black","red","green", "blue","lightsalmon","#660066","brown","grey","violet","gold"),lty=1,cex=0.8)


