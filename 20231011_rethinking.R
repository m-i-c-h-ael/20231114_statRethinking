#install.packages("cmdstanr", repos = c("https://mc-stan.org/r-packages/", getOption("repos")))
#library('devtools')
#library('githubinstall')
#githubinstall('rethinking')

library('rethinking')

#Chapter 2
globe= quap(
  alist( W ~ dbinom(W+L, p),
         p ~ dunif(0,1)),
  data= list( W=6, L=3)
)

precis(globe)

curve(dnorm(x, mean= precis(globe)$mean,sd= precis(globe)$sd),type='l',from=0,to=1)
curve(dbinom(x=6,size=6+3,prob=x)*dunif(x,min=0,max=1)/
        integrate(function(x) {dbinom(x=6,size=6+3,prob=x)*dunif(x,min=0,max=1)},lower=0,upper=1)$value,
      add=TRUE, col='red',type='l')


# 2M1: Grid approx
theta= seq(0,1,.1)
len=length(theta)

P_theta= rep(1,len) #uniform prior

plot(theta,P_theta,main='Prior')

L1= theta^3
post1= L1*P_theta / sum(L1[1:len-1]*P_theta[1:len-1]*diff(theta))
plot(theta,post1,main='Post1')

L2= theta^3 * (1-theta)
post2= L2*P_theta / sum(L2[1:len-1]*P_theta[1:len-1]*diff(theta))
plot(theta,post2,main='Post2')

L3= theta^5 * (1-theta)^2
post3= L3*P_theta / sum(L3[1:len-1]*P_theta[1:len-1]*diff(theta))
plot(theta,post3,main='Post3')

Prior2= rep(0,length(theta))
Prior2[theta>= 0.5]= 2
post1a= L1*Prior2 / sum(L1[1:len-1]*Prior2[1:len-1]*diff(theta))
post2a= L2*Prior2 / sum(L2[1:len-1]*Prior2[1:len-1]*diff(theta))
post3a= L3*Prior2 / sum(L3[1:len-1]*Prior2[1:len-1]*diff(theta))
plot(theta,Prior2,ylim=c(0,4))
lines(theta,post1a,col='blue')
lines(theta,post2a,col='green')
lines(theta,post3a,col='red')


########
# 3H1
rm(list=ls())
cat('\014')
graphics.off()

library('rethinking')
data(homeworkch3)

theta_vec= seq(0,1,0.001)
prior= rep(1,length(theta_vec))
birth1Tab= table(birth1)
birth2Tab= table(birth2)
allBirths= table(c(birth1,birth2))

l= dbinom(x=allBirths[2],size=sum(allBirths),p=theta_vec)
post= prior*l/sum(prior*l)
#sum(post)
plot(theta_vec,post,type='b')

#posterior maximum
(theta_postMax= theta_vec[ which(post==max(post))])

#boy next
l_boy= dbinom(x=1,size=1,prob=theta_vec)
(P_boy_next= sum(l_boy * post))

# 3H2: sample from posterior
postSamp= sample(theta_vec,size=10000,replace=TRUE,prob=post)

(p50= rethinking::HPDI(postSamp,prob=.5))
(p89= rethinking::HPDI(postSamp,prob=.89))
(p97= rethinking::HPDI(postSamp,prob=.97))

# 3H3: simulate number boys
simFromPost= rbinom(length(postSamp),size=200,p=postSamp)
plot(density(simFromPost))
abline(v=allBirths[2],col='red')


# Chapter 4: Distribution from product
ranProd= replicate(10000, prod(runif(n=12, min=1, max=1.1)))
plot(density(ranProd))
1.1^12
plot(density(log(ranProd)))

big= replicate(10000, prod(1+ runif(n=12,min=0,max=0.5)))
plot(density(big))

### Howell1
data(Howell1)
d= Howell1
precis(d)

# Prior predictive distribution of heights
mu_prior= rnorm(1e4, mean=178, sd=20)
sig_prior= runif(1e4, min=0, max=50)
pp= rnorm(1e4,mean= mu_prior,sd=sig_prior)

plot(density(pp))


d2= d[d$age>18,]
mu_vec= 150:160
sig_vec= seq(7,9,.1)
grd= expand.grid(mu= mu_vec,sig= sig_vec)
grd$LL= apply(grd,1,function(x){ sum(log(dnorm(x=d2$height,mean=x[1],sd=x[2]))) })
 #apply(grd,1,function(x){ exp( sum(log(dnorm(x=d2$height,mean=x[1],sd=x[2]))) )})
 
#LL <- sapply( 1:nrow(grd), function(i) 
  #sum(  dnorm( d2$height, grd$mu[i], grd$sig[i],log=TRUE)) )
                                                           
grd$muPri= dnorm(grd$mu,mean=178,sd=20)
grd$sdPri= dunif(grd$sig,min=0,max=50)
grd$logPost= grd$LL+log(grd$muPri)+log(grd$sdPri)
grd$Post= exp( grd$logPost - max(grd$logPost) )  #why?

library('plotly')
plotly::plot_ly() %>%
  add_markers(data= grd,
              x= ~mu,y= ~sig,
              z= ~Post, type='scatter3d',
              #color=DF_3D_melt$prob,mode='markers',
              size=0.1) %>%
  layout(showlegend= FALSE) #%>%
  #hide_colorbar() %>%
  #layout(scene= scene, title= '...') 

unique(grd$mu)
unique(grd$sig)

contour_xyz <- function( x , y , z , ... ) {
  ux <- unique(x)
  uy <- unique(y)
  n <- length(ux)
  ny <- length(uy)
  m <- matrix( z , nrow=n , ncol=ny )
  contour( ux , uy , m , ... )
}

contour_xyz(x=grd$mu, y=grd$sig,z=grd$Post)
