#Work.Dir <- "E:/Dropbox/Xia/08 EM algorithm for misrepresentation/CAS paper/code"
#setwd(Work.Dir)

##############
n=1000
theta0=0.5
p0=0.25
V = rbinom(n,1,theta0)
V_star = V
V_star[V==1] = rbinom(sum(V==1),1,1-p0)

### Generate Y ###
aaa0 = 1.2
aaa1 = 1
#aaa2 = 0.5
#mu = exp(aaa0+aaa1*V+aaa2*X)
mu = exp(aaa0+aaa1*V)
phi = 0.2
alpha0 = 1/phi
beta = 1/mu/phi
Y = rgamma(n, alpha0, rate=beta)
y=Y
v_star=V_star
plot(density(y[v_star==0]),main="",lty=2)
lines(density(y[v_star==1]))
m=rep(NA,2)
v=m
m[1]<-mean(y[v_star==0])
m[2]<-mean(y[v_star==1])
v<-var(y[v_star==1])
beta=m/v
alpha=m^2/v

# EM-algorithm
source("gammaMisrepEM 09-01-2018 Xia.R")

out=gammaMisrepEM(y,v_star, verb = TRUE)
out$reg.pars
out$gamma.pars
out$lambda

# Compare with naive and true estimates
glm.naive=glm(y~v_star, family="Gamma"(link='log'))
glm.naive$coefficients
as.numeric(gamma.shape(glm.naive))[1]

glm.true=glm(y~V, family="Gamma"(link='log'))
glm.true$coefficients
as.numeric(gamma.shape(glm.true))[1]

