gammaMisrepEM<-function(y,v_star,lambda = c(0.6,0.4), epsilon = 1e-08, maxit = 10000, maxrestarts=20, verb = FALSE) {
  library(MASS)
  # obtain initial values from gamma glm with log link
  naive=glm(y~v_star, family="Gamma"(link='log'))
  coef.reg=as.numeric(naive$coefficients)    # coef.reg[1]=estimate of intercept; coef.reg[2]=estimate of v_star
  alpha=as.numeric(gamma.shape(naive))[1]
  theta.reg <- c(coef.reg,alpha)           # parameters under regression setting

  mu=c(exp(coef.reg[1]),exp(sum(coef.reg))) # gamma mean
  beta = alpha/mu
  theta = c(alpha,beta)

  iter <- 0
  mr <- 0
  diff <- epsilon+1
  n <- length(y)
  
  # observed loglikelihood
  obs.ll <- function(lambda,theta){
    sum(v_star*log(dgamma(y,theta[1],theta[3])))+
    sum((1-v_star)*log(lambda[2]*dgamma(y,theta[1],theta[3])+lambda[1]*dgamma(y,theta[1],theta[2])))
  }
 
  # M step loglikelihood
  mstep.ll <- function(theta,z){
    -sum(log(dgamma(y[v_star==1],theta[1],theta[3])))-
    sum((1-z[v_star==0])*log(dgamma(y[v_star==0],theta[1],theta[3]))+
         z[v_star==0]*log(dgamma(y[v_star==0],theta[1],theta[2])))
  }
    
  old.obs.ll <- obs.ll(lambda,theta)
  ll <- old.obs.ll

  while(diff > epsilon && iter < maxit){

      # E-step
	dens1=lambda[1]*dgamma(y,alpha,beta[1])
	dens2=lambda[2]*dgamma(y,alpha,beta[2])
	z=dens1/(dens1+dens2)
	lambda.hat=c(mean(z[v_star==0]),(1-mean(z[v_star==0])))

      # pseudo observations for M-step based on glm() 
      #res=y[v_star==1]
      #cov=v_star[v_star==1]  
      #wt=rep(1,length(res))
      #res=c(res,y[v_star==0],y[v_star==0])
      #cov=c(cov,v_star[v_star==0],1-v_star[v_star==0])
      #wt=c(wt,z[v_star==0],1-z[v_star==0])

      # M-step, cannot use glm() since the weights change the dispersion parameter
      m=try(suppressWarnings(nlm(mstep.ll,p=theta,z=z)),silent=TRUE)

	#mstep=glm(res~cov, weights=wt, family="Gamma"(link='log'))
      #coef.reg.hat=as.numeric(mstep$coefficients)
      #alpha.hat=as.numeric(gamma.shape(mstep))[1]          # from non linear maximization
      #theta.reg.hat <- c(coef.reg.hat,alpha.hat)           # parameters under regression setting

      #mu.hat=c(exp(coef.reg.hat[1]),exp(sum(coef.reg.hat))) # gamma mean
      #beta.hat = alpha.hat[1]/mu.hat
      #theta.hat = c(alpha.hat,beta.hat)
      alpha.hat=m$estimate[1]
      beta.hat=m$estimate[2:3]
      theta.hat=m$estimate

	new.obs.ll <- obs.ll(lambda.hat,theta.hat)
	diff <- new.obs.ll-old.obs.ll
	old.obs.ll <- new.obs.ll
	ll <- c(ll,old.obs.ll)
	lambda=lambda.hat
	theta=theta.hat
	alpha=alpha.hat
	beta=beta.hat
	iter=iter+1
         if (verb) {
            cat("iteration =", iter, " log-lik diff =", diff, " log-lik =", 
                new.obs.ll, "\n")
          }
   }

    if (iter == maxit) {
        cat("WARNING! NOT CONVERGENT!", "\n")
    }
    cat("number of iterations=", iter, "\n")
    names(theta)=c("alpha","beta1","beta2")

    # transform the parameters back to regression setting
    m1=alpha/beta[1]
    m2=alpha/beta[2]
    b0=log(m1)
    b1=log(m2/m1)
    theta.reg=c(b0,b1,alpha)
    names(theta.reg)=c("intercept","v_star","alpha")

    a=list(y=y, lambda = lambda, reg.pars = theta.reg, gamma.pars = theta, loglik = new.obs.ll, 
        posterior = z, all.loglik=ll, ft="gammaMisrepEM")
    class(a) = "misrepEM"
    a
}	

