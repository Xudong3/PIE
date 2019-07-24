
install.packages(c("logitnorm","pracma","lme4"))
library("logitnorm","pracma","lme4")

n.sample=2
n.sim=150

#mu.x,sig.x
true.x.par=c(0,1)
#beta.0, beta.1
true.beta.par=c(-1,1)
#mu.0,mu.1,sig.0,sig.1
true.random.par=c(2,2,1.5,1.5)

#generate data
generate.data<-function(x.par,  beta.par, random.par){
  
  mu.x=x.par[1]
  sig.x=x.par[2]
  
  beta.0=beta.par[1]
  beta.1=beta.par[2]
  
  mu.0=random.par[1]
  mu.1=random.par[2]
  sig.0=random.par[3]
  sig.1=random.par[4]
  
  x=rnorm(n.sample,mu.x, sqrt(sig.x))
  
  # calculate P(y=1)
    p.y=exp(beta.0+beta.1*x)/(1+exp(beta.0+beta.1*x))
    y=rbinom(n.sample, 1, p.y)
  
  #generate random effect
    random.0=rlogitnorm(n.sample, mu.0, sqrt(sig.0))
    random.1=rlogitnorm(n.sample, mu.1, sqrt(sig.1))
  
  #calculate P(s=1)
    p.s=(1-random.0)+(random.0+random.1-1)*p.y
    s=rbinom(n.sample, 1, p.s)

  #return data
    data.frame(x,s,y, random.0, random.1, index=1:n.sample )
}

complete.data=generate.data(true.x.par, true.beta.par, true.random.par)
observed.data=complete.data[c("x","s")]
missing.data=complete.data[c("y", "random.0", "random.1")]

summary(missing.data["y"])

summary(missing.data[,"random.0"])
summary(missing.data[,"random.1"])

## naive fit using  missing data 
glm(s~ x, data=complete.data,family="binomial")
naive.est=coef(glm(s~ x, data=complete.data,family="binomial"))
bias.naive=naive.est-true.beta.par


#generate random effect alpha_i0 and alpha_i1 conditional on observed data using metropolis hasting algorithm i.e., conditional distribution f(alpha_i0\vert s_i) and 
#f(alpha_i1 \vert s_i )
#observed.data: first column is X, second column is S

mh<-function(beta.par, random.par, start=c(0.5,0.5)){
  
  beta.0=beta.par[1]
  beta.1=beta.par[2]
  
  mu.0=random.par[1]
  mu.1=random.par[2]
  sig.0=random.par[3]
  sig.1=random.par[4]
 
  
#set up matrix random[, ] conditional on s. column is for i (people) and row is for t (simulation)
  random.0=matrix(NA, nrow=n.sim, ncol=n.sample)
  random.1=matrix(NA, nrow=n.sim, ncol=n.sample)
  random.0[1, ]=rep(start[1], n.sample)
  random.1[1, ]=rep(start[2], n.sample)


#count the number of accepted proposals
  accept=rep(0, n.sample)
  
#loop
  for (t in 2:n.sim){
    
    #generate a candidate value from logisticnormal distribution
      random.0.star=rlogitnorm(n.sample,mu.0,sqrt(sig.0))
      random.1.star=rlogitnorm(n.sample,mu.1,sqrt(sig.1))
      
    #generate a random number from uniform distribution
      u=runif(n.sample,min=0, max=1)
   
      
    #calculate acceptance probability
      #p(s=1) conditional on random.0.star and random.1.star
      p.star=(1-random.0.star)+(random.0.star+random.1.star-1)*exp(beta.0+beta.1*observed.data[,"x"])/(1+exp(beta.0+beta.1*observed.data[,"x"]))
      
      
      #p(s=1) conditional on previous random.0 and random.1
      p.pre=(1-random.0[t-1, ])+(random.0[t-1, ]+random.1[t-1,]-1)*exp(beta.0+beta.1*observed.data[,"x"])/(1+exp(beta.0+beta.1*observed.data[,"x"]))
      
      #calculuate the ratio
      r=dbinom(observed.data[,"s"],size=1,p.star)/dbinom(observed.data[,"s"],size=1,p.pre)
      
      #a= min(r,1)
      a=sapply(r, function(x) min(x,1))
      
      
      #loop over i
        for (i in 1:n.sample){
          if (u[i]<a[i]){
            accept[i]=accept[i]+1
            random.0[t,i]=random.0.star[i]
            random.1[t,i]=random.1.star[i]
          }else{
            random.0[t,i]=random.0[t-1,i]
            random.1[t,i]=random.1[t-1,i]
          }
        }
      }
   list(random.0,random.1, accept/n.sim)
  
    }



#test mh at true value
mh(true.beta.par, true.random.par)


#EM algorithm
#observed.data: first column is X, second column is S
em<- function(start){
  
  #number of iterations
  n.iters=0
  
  #beta.0, beta.1
  beta.par.cur=start[c(1,2)]
  #mu.0, mu.1, sig.0, sig.1
  random.par.cur=start[3:6]
  
  beta.par.next=start[c(1,2)]
  random.par.next=start[3:6]
  
  #discrepany=1
 
  #Loop 
  for (i in 1:400){
  
    
    #update
    beta.par.cur=beta.par.next
    random.par.cur=random.par.next
    
    #apply EM algorithm to get the random effect
        random.data.cur=mh(beta.par.cur, random.par.cur)
  
    #E step: define the empirical q function
      q<-function(beta.par.next,random.par.next=true.random.par ){
        
        #p(y=1) conditional on beta.par.next
        p.y.next= exp(beta.par.next[1]+ beta.par.next[2]*observed.data[,"x"])/(1+exp(beta.par.next[1]+beta.par.next[2]*observed.data[,'x']))
        
        
        #P(s=1) conditional on current random effect
        p.s.cur=matrix(NA, nrow=n.sim, ncol=n.sample)
        for (i in 1:n.sample){
          p.s.cur[ , i]=(1-random.data.cur[[1]][, i])+(random.data.cur[[1]][,i]+random.data.cur[[1]][ ,i]-1)*p.y.next[i]
          
        }
        
        #logdensity function for s conditional on current random effect and sum over all possible i and t 
        log.density.s.cur=sum(log(dbinom(t(replicate(n.sim, observed.data[,"s"])),size=1,p.s.cur)))        
        
        
        #logdensity function for random effect conditional on current random effect and sum over all possible i and t
        log.density.random.cur=sum(log(dlogitnorm(random.data.cur[[1]],random.par.next[1], sqrt(random.par.next[3]))
                                       *dlogitnorm(random.data.cur[[1]],random.par.next[2], sqrt(random.par.next[4]))))
        
        
        result =log.density.s.cur/n.sim+ log.density.random.cur/n.sim
        result
      }
    
      
      q(beta.par.next=true.beta.par)
    #M step 
    
     pars=rep(0.5,2)
     par.next= optim(pars,q, method="BFGS",control=list(fnscale=-1,parscale=c(1/n.sample, 1/n.sample))) #1/n.sim,1/n.sim,1/n.sim,1/n.sim)))
     par.next=unname(unlist(par.next["par"]))
     par.next
     
     nlminb
     
     print(c(par.next,i))
     
     
  }

}

em(rep(0.2,6))

##log likelihood for the whole data 
log.full.likelihood.data<-function(beta.par, random.par=true.random.par){

    ##log likelihood  for a single observation
    log.full.likelihood.point<-function(beta.par, random.par, observed.point){
  
        #integrand function
        integrand<-function(random.0, random.1){
  
                  #p(y=1)
                  p.y=exp(beta.par[1]+beta.par[2]*observed.point)/(1+exp(beta.par[1]+beta.par[2]*observed.point))
  
                  #p(s=1)
                  p.s=(1-random.0)+(random.0+random.1-1)*p.y
  
                  #integrand
                  dbinom(observed.point,size=1, prob=p.s)*dlogitnorm(random.0, true.random.par[1],sqrt(true.random.par[3]))*
                  dlogitnorm(random.1, true.random.par[2],sqrt(true.random.par[4]))
                }

    #quadrature function
    return(log(pracma::dblquad(integrand, xa=0,xb=1,ya=0,yb=1)))
  }

-sum(sapply(observed.data[,"s"], function(x) log.full.likelihood.point(beta.par, random.par, x)))
}


pars=rep(0.8,2)
optim(pars, log.full.likelihood.data, method="BFGS",control=list(fnscale=-1,parscale=c(1/n.sample, 1/n.sample)))
nlminb(start=rep(0.7,2), full.likelihood.data)   

full.likelihood.point<-function( beta.par, random.par, observed.point=1){
  
  #integrand function
  integrand<-function(random.0, random.1){
    
    #p(y=1)
    p.y=exp(beta.par[1]+beta.par[2]*observed.point)/(1+exp(beta.par[1]+beta.par[2]*observed.point))
    
    #p(s=1)
    p.s=(1-random.0)+(random.0+random.1-1)*p.y
    
    
    #integrand
    dbinom(observed.point,size=1, prob=p.s)*dlogitnorm(random.0, true.random.par[1],sqrt(true.random.par[3]))*dlogitnorm(random.1, true.random.par[2],sqrt(true.random.par[4]))
  }
  
  #quadrature function
  return(pracma::dblquad(integrand, xa=0,xb=1,ya=0,yb=1,dim=2))
}

full.likelihood.point(beta.par=c(1,1))


