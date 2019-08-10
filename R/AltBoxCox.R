# Get M spline and I spline matrix with order
#
# @param x is a row vector
# @param k is the order of I spline
# @param knots are a sequence of increasing points
# The number of free parameters in M spline is the length of knots plus 1.
# This function returns a list which consists of I splines and M splines basis
MIspline<-function(x,order,knots){
  ### get Mspline bases ###
  k1=order
  m=length(knots)
  n1=m-2+k1 # number of parameters
  t1=c(rep(1,k1)*knots[1], knots[2:(m-1)], rep(1,k1)*knots[m]) # newknots

  tem1=array(rep(0,(n1+k1-1)*length(x)),dim=c(n1+k1-1, length(x)))
  for (l in k1:n1){
    tem1[l,]=(x>=t1[l] & x<t1[l+1])/(t1[l+1]-t1[l])
  }

  if (order==1){
    mbases=tem1
  }else{
    mbases=tem1
    for (ii in 1:(order-1)){
      tem=array(rep(0,(n1+k1-1-ii)*length(x)),dim=c(n1+k1-1-ii, length(x)))
      for (i in (k1-ii):n1){
        tem[i,]=(ii+1)*((x-t1[i])*mbases[i,]+(t1[i+ii+1]-x)*mbases[i+1,])/(t1[i+ii+1]-t1[i])/ii
      }
      mbases=tem
    }
  }

  ### get Ispline bases ###
  k=order+1
  n=m-2+k # number of parameters
  t=c(rep(1,k)*knots[1], knots[2:(m-1)], rep(1,k)*knots[m]) # newknots

  yy1=array(rep(0,(n+k-1)*length(x)),dim=c(n+k-1, length(x)))
  for (l in k:n){
    yy1[l,]=(x>=t[l] & x<t[l+1])/(t[l+1]-t[l])
  }

  yytem1=yy1
  for (ii in 1:order){
    yytem2=array(rep(0,(n+k-1-ii)*length(x)),dim=c(n+k-1-ii, length(x)))
    for (i in (k-ii):n){
      yytem2[i,]=(ii+1)*((x-t[i])*yytem1[i,]+(t[i+ii+1]-x)*yytem1[i+1,])/(t[i+ii+1]-t[i])/ii
    }
    yytem1=yytem2
  }


  index=rep(0,length(x))
  for (i in 1:length(x)){
    index[i]=sum(t<=x[i])
  }

  ibases=array(rep(0,(n-1)*length(x)),dim=c(n-1,length(x)))

  if (order==1){
    for (i in 2:n){
      ibases[i-1,]=(i<index-order+1)+(i==index)*(t[i+order+1]-t[i])*yytem2[i,]/(order+1)
    }
  }else{
    for (j in 1:length(x)){
      for (i in 2:n){
        if (i<(index[j]-order+1)){
          ibases[i-1,j]=1
        }else if ((i<=index[j]) && (i>=(index[j]-order+1))){
          ibases[i-1,j]=(t[(i+order+1):(index[j]+order+1)]-t[i:index[j]])%*%yytem2[i:index[j],j]/(order+1)
        }else{
          ibases[i-1,j]=0
        }
      }
    }
  }
  return(list(mbases,ibases))
}



# Estimate the trasnformation G on the response y such that we have linear relationship between x and G(y)
# We use monotone splines to estimate G
#
# @param x the design matrix which does not contain the intercept
# @param y the response vector
# @param eps the stopping criteria for the internal EM algorithm
# This function returns a list which contains estimated coefficients of monoton splines and of the covarites
NonparaBoxCox<-function(x,y,eps=0.000001){
  y.modified = y
  y.modified[which.max(y.modified)] = y.modified[which.max(y.modified)] + .01
  y.modified[which.min(y.modified)] = y.modified[which.min(y.modified)] - .01
  knots = quantile(y.modified,seq(0,1,length=15))
  K = length(knots)
  MI=MIspline(t(y),2,knots)

  Isp=MI[[2]]
  Msp=MI[[1]]
  all.knots[it,] = knots

  gammahat.cd = c(rep(0.01,nrow(Isp)),-5)
  betahat = rep(.1,ncol(x))
  new.gammahat.cd = gammahat.cd

  repeat{ # EM
    post = t(Msp)*matrix(rep(gammahat.cd[-(K+1)],n),nrow=n,byrow = T)
    post = post/rowSums(post)
    for(i in 1:(K+1)){
      if(i<K+1) {
        Wpost = t(Isp)[,-i]*matrix(rep(new.gammahat.cd[-c(i,K+1)],n),nrow=n,byrow = T)*matrix(rep(t(Isp)[,i],K-1),nrow=n)
        c=sum(post[,i])
        b=-( sum(Wpost) +sum(t(Isp)[,i])*new.gammahat.cd[K+1] - sum(x%*%betahat*t(Isp)[,i]))
        a=-sum((t(Isp)[,i])^2)
        if(a==0){new.gammahat.cd[i] = -c/b}
        if(!a==0){new.gammahat.cd[i] = (-b-sqrt(b^2-4*a*c))/(2*a)}
      }
      if(i==K+1){new.gammahat.cd[i] = -mean(t(Isp)%*%new.gammahat.cd[-(K+1)] - x%*%betahat)}
    }
    new.betahat = solve(t(x)%*%x)%*%t(x)%*%(t(Isp)%*%new.gammahat.cd[1:K]+new.gammahat.cd[K+1])

    if(max(abs(new.gammahat.cd-gammahat.cd),abs(new.betahat-betahat))<eps){break}
    gammahat.cd = new.gammahat.cd
    betahat=new.betahat
  }
  return(list(gammahat.cd,betahat))
}
