#install.packages("cmdstanr", repos = c("https://mc-stan.org/r-packages/", getOption("repos")))
#library('devtools')
#install.packages('githubinstall')
#library('githubinstall')
#githubinstall('rethinking')

library('rethinking')
library('ellipse')
library('MASS')

data(Howell1)
d = Howell1
d_adult = d[d$age > 18, ]
precis(d_adult)

hm = quap(alist(height ~ dnorm(mu, sig),
                mu ~ dnorm(178, 20),
                sig ~ dunif(0, 50)),
          data = d_adult)

(hm_coef = precis(hm))
(hm_covMtx= vcov(hm))

cov2cor(hm_covMtx)  #correlation

#simplistic: ignoring co-variance
mu_sim = rnorm(5000, mean = hm_coef[1, 1], sd = hm_coef[1, 2])
sd_sim = rnorm(5000, mean = hm_coef[2, 1], sd = hm_coef[2, 2])
par(mfrow = c(3, 1))
for (i in 1:3) {
  curve(dnorm(x, mu_sim[i], sd_sim[i]), from = 120, to = 230)
}
x_sim = rnorm(1000, mu_sim, sd_sim)
par(mfrow = c(1, 1))
hist(x_sim, col = 'skyblue')

plot(mu_sim, sd_sim)
elp = ellipse(hm_covMtx, centre = hm_coef[, 1])
lines(elp, col = 'red', lwd = 2)


# More correct sampling (considering covariance)
mvRnd = mvrnorm(n = 5000, mu = hm_coef[, 1], hm_covMtx)
plot(mvRnd)
lines(elp, col = 'red', lwd = 2)

# log-normal distribution
rnd = seq(0.1, 10, .1)
prods = sapply(rep(10, 10000), function(x) {
  prod(sample(rnd, x, replace = TRUE))
})
plot(density(prods), xlim = c(0, 10000000))

log_mean = mean(log(prods))
log_sd = sd(log(prods))

curve(dlnorm(x, meanlog = log_mean, sdlog = log_sd),
      add = TRUE,
      col = 'red')

# chapter 5




##########
# B-Splines
library('rethinking')
library('splines')

data(cherry_blossoms)
head(cherry_blossoms)
cher = cherry_blossoms[, 1:2]
cher = cher[apply(cher, 1, function(x) {
  !any(is.na(x))
}), ]

knots = quantile(cher$year, c(.25, .5, .75))

# B_mtx= matrix(0,ncol=length(cher$year),nrow= length(knots))
# B_mtx[1, cher$year<=knots[1]]= 1/knots[1] * cher$year[cher$year<=knots[1]]
# B_mtx[2, cher$year>knots[1] & cher$year<=knots[2] ]=
#           1/(knots[2]-knots[1]) * (cher$year[cher$year>knots[1]&cher$year<=knots[2]] -knots[1])
# B_mtx[2, cher$year>knots[2] & cher$year<=knots[3] ]=
#   1/(knots[2]-knots[3]) * (cher$year[cher$year>knots[2]&cher$year<=knots[3]] -knots[2])

B = bs(
  x = cher$year,
  knots = knots,
  degree = 1,
  intercept = TRUE
)
for (i in 1:dim(B)[2]) {
  plot(cher$year, B[, i], main = paste(i))
}

m4.7 = quap(
  alist(
    doy ~ dnorm(mu, sig),
    mu <- a + B %*% w,
    a ~ dnorm(100, 10),
    w ~ dnorm(0, 10),
    sig ~ dexp(1)
  ),
  data = list(B = B, doy = cher$doy),
  start = list(w = rep(0, dim(B)[2]))
)


# Plot the cher data with the posterior mean curve and HPDI
plot(cher$year, cher$doy, pch = 16, cex = 0.6)
post = extract.samples(m4.7)

a_mean = mean(post$a)
w_means = apply(post$w, 2, mean)
mu_post = a_mean + B %*% w_means
lines(cher$year, mu_post, lwd = 2, col = 'red')

mu_samp = matrix(NA, nrow = dim(post$w)[1], ncol = dim(B)[1])
for (i in 1:dim(post$w)[1]) {
  mu_samp[i, ] = post$a[i] + B %*% post$w[i, ]
}
mu_hpdi= apply(mu_samp,2,HPDI, prob=0.95)
shade( mu_hpdi, cher$year)
lines(cher$year, mu_hpdi[1,], lwd = 2, col ='red', lty = 2)
lines(cher$year, mu_hpdi[2,], lwd = 2, col = 'red', lty = 2)

# Visualizing multivariate distributions
library('rethinking')
library('scatterplot3d')
library('vrmlgen')

#list all datasets in the rethinking-package
data(package = 'rethinking')
data("WaffleDivorce")
wd = WaffleDivorce
head(wd)

# Plot the data
plot(wd$MedianAgeMarriage, wd$Divorce, pch=16, cex=0.6)

scatterplot3d(x=wd$MedianAgeMarriage, y=wd$Marriage, z=wd$Divorce) #, pch=16, cex=0.6)
mod3D= lm(Divorce ~ MedianAgeMarriage + Marriage, data=wd)
summary(mod3D)

# Add regression lines to 3D scatterplot
scatterplot3d(x=wd$MedianAgeMarriage, y=wd$Marriage, z=wd$Divorce, pch=16)
x = 23:30
y = 0:35
z_pred = outer(
  X = x,
  Y = y,
  FUN = function(x, y) {
    coef(mod3D)[1] + coef(mod3D)[2] * x + coef(mod3D)[3] * y
  }
)

#######
coef_x = 3
coef_y = 0
x = seq(0, 10, .5)
y = seq(0, 10, .5)
#df= expand.grid(x=x,y=y)
z = outer(
  X = x,
  Y = y,
  FUN = function(x, y) {
    coef_x * x + coef_y * y
  }
)
persp(x,y,z, col = 'white', theta = -40, phi = 30)

################
### Milk calories vs female body mass and neocortex size
rm(list = ls())
cat('\014')
graphics.off()

library('rethinking')
library('tidyverse')
library('GGally')
data(milk)
head(milk)
milk = milk[!is.na(milk$neocortex), ]

# Standardize
milk$kcal_s = (milk$kcal.per.g - mean(milk$kcal.per.g)) / sd(milk$kcal.per.g)
milk$nec_s = (milk$neocortex - mean(milk$neocortex)) / sd(milk$neocortex)
milk$mass_s = (milk$mass - mean(milk$mass)) / sd(milk$mass)

# Neocortex-only model
mN = quap(alist(
  kcal_s ~ dnorm(mu, sig),
  mu <- a + bN * nec_s,
  
  #Priors
  a ~ dnorm(0, .2),
  bN ~ dnorm(0, .5),
  sig ~ dexp(1)
),
data = milk)

# Prior predictive
mN_priorSamp = extract.prior(mN)
mN_priorPred = with(mN_priorSamp, sapply(1:length(a), function(i)
  a[i] + bN[i] * milk$nec))
#mN_priorPred_mu= apply(mN_priorPred,2,mean)
#mN_priorPred_hpdi= apply(mN_priorPred,2,HPDI, prob=0.89)
plot(type = 'n', x = 0, y = 0, xlim = c(-2,2), ylim = c(-2,2))
for(i in 1:100) {
  lines(milk$nec_s, mN_priorPred[, i], col = rgb(0, 0, 0, .1))
}

# Get posterior samples
mN_post_coef = precis(mN)
mN_postSamp = extract.samples(mN)
nec_vec = seq(-2, 2, .1)
mN_postPred = with(mN_postSamp, sapply(1:length(a), function(i)
  a[i] + bN[i] * nec_vec))
mN_postPred_mu = apply(mN_postPred, 1, mean)
mN_postPred_hpdi = apply(mN_postPred, 1, HPDI, prob = 0.89)

# Plot posterior predictive
plot(milk$nec_s, milk$kcal_s, pch = 16, cex = 0.6,col = 'blue')
lines(nec_vec, mN_postPred_mu, lwd = 2, col = 'red')
shade(mN_postPred_hpdi, nec_vec)

# Pairs-plot
GGally::ggpairs(milk[, c('kcal_s', 'nec_s', 'mass_s')])

# Plot in 3D with projections to the sides
 #http://www.sthda.com/english/wiki/impressive-package-for-3d-and-4d-graph-r-software-and-data-visualization
library('plot3D')
scatter3D_fancy <- function(x, y, z,...,
          xlim,
          mainCol = 'black',XYwallCol = 'blue', XZwallCol = 'blue', YZwallCol = 'blue')
{
  panelfirst <- function(pmat) {
    #XY <- trans3D(x, y, z = rep(min(z), length(z)), pmat = pmat)
    #scatter2D(XY$x, XY$y, col= XYwallCol, pch = ".", 
    #          cex = 4, add = TRUE, colkey = FALSE)
    
    XY <- trans3D(x = rep(min(x), length(x)), y, z, pmat = pmat)
    scatter2D(XY$x, XY$y, col=YZwallCol, pch = ".", 
              cex = 4, add = TRUE, colkey = FALSE)
    
    XY <- trans3D(x , y = rep(max(y), length(y)), z, pmat = pmat)
    scatter2D(XY$x, XY$y, col=XZwallCol, pch = ".", 
              cex = 4, add = TRUE, colkey = FALSE)
  }
  scatter3D(x, y, z, ..., col = mainCol, panel.first = panelfirst)
}

# scatter3D_fancy(x=milk$nec_s,y=milk$mass_s,z=milk$kcal_s, bty='b2')

#######
data(WaffleDivorce)
dv = WaffleDivorce
head(dv)
scatter3D_fancy(x = dv$MedianAgeMarriage, y = dv$Marriage, z = dv$Divorce, 
                bty ='b2', theta = 30, phi = 30,
                xlim = c(20,30), ylim = c(13,30) )

A = min(dv$MedianAgeMarriage):max(dv$MedianAgeMarriage)
M = min(dv$Marriage):max(dv$Marriage)
AMgrid = expand.grid(A = A, M = M)

# Age as full predictor
lmA = lm(Divorce ~ MedianAgeMarriage, data = dv)
AMgrid$D_lmA = coef(lmA)[1] + coef(lmA)[2] * AMgrid$A
scatter3D_fancy(x = dv$MedianAgeMarriage, y = dv$Marriage, z = dv$Divorce, 
                bty = 'b2', theta = 30, phi = 30,
                xlim = c(20,30), ylim = c(13,30))
scatter3D_fancy(x = AMgrid$A, y = AMgrid$M, z = AMgrid$D_lmA,
                bty = 'b2', theta = 30, phi = 30, add = TRUE,
                mainCol = 'orange', type = 'l', XZwallCol = 'green', YZwallCol = NA)

# Marriage rate as full predictor
lmM = lm(Divorce ~ Marriage, data = dv)
AMgrid$D_lmM = coef(lmM)[1] + coef(lmM)[2] * AMgrid$M
scatter3D_fancy(x = dv$MedianAgeMarriage, y = dv$Marriage, z = dv$Divorce,
                bty = 'b2', theta = 30, phi = 30,
                xlim = c(20,30), ylim = c(13,30))
scatter3D_fancy(x = AMgrid$A, y = AMgrid$M, z = AMgrid$D_lmM, bty ='b2',
                theta = 30, phi = 30, add = TRUE,
                mainCol = 'orange', type = 'l', XZwallCol = NA, YZwallCol = 'green')

# Age and marriage rate as predictors
lmAM = lm(Divorce ~ MedianAgeMarriage + Marriage , data = dv)
AMgrid$D_lmAM = coef(lmAM)[1] + coef(lmAM)[2] * AMgrid$A + coef(lmAM)[3] *
  AMgrid$M
scatter3D_fancy(x = dv$MedianAgeMarriage, y = dv$Marriage, z = dv$Divorce,
                bty = 'b2', theta = 30, phi = 30,
                xlim = c(20,30), ylim = c(13,30))
scatter3D_fancy(x = AMgrid$A, y = AMgrid$M, z = AMgrid$D_lmAM,
                bty = 'b2', theta = 30, phi = 30, add = TRUE,
                mainCol = 'orange', type = 'l', XZwallCol = NA, YZwallCol = NA)

library('dagitty')

# Make causal predictions based on model (counterfactual plot)
g <- dagitty("dag{
  M -> D ;
  M <- A -> D
 }")
plot(g)

# Causal effect of marriage rate
 #A can be simulated independent of M (as done for AMgrid) as M is not caused by A
scatter3D_fancy(x = AMgrid$A, y = AMgrid$M, z = AMgrid$D_lmAM,
                bty ='b2', theta = 30, phi = 30,
                xlim = c(20,30), ylim = c(13,30))
scatter3D_fancy(x = AMgrid$A, y = AMgrid$M, z = AMgrid$D_lmAM,
                bty = 'b2', theta = 30, phi = 30, add = TRUE,
                mainCol = 'orange', type = 'l', XZwallCol = NA, YZwallCol = NA)

# Causal effect of age
lm_MfromA = lm(Marriage ~ MedianAgeMarriage, data = dv)
AMgrid$MfromA = coef(lm_MfromA)[1] + coef(lm_MfromA)[2] * AMgrid$A
AMgrid$D_causA = coef(lmAM)[1] + coef(lmAM)[2] * AMgrid$A + coef(lmAM)[3] *
  AMgrid$MfromA
scatter3D_fancy(x = AMgrid$A, y = AMgrid$MfromA, z = AMgrid$D_causA,
                bty = 'b2', theta = 30, phi = 30,
                xlim = c(20,30), ylim = c(13,30))
scatter3D_fancy(x = AMgrid$A, y = AMgrid$M, z = AMgrid$D_lmAM,
                bty = 'b2', theta = 30, phi = 30, add = TRUE,
                mainCol = 'orange', type ='l', XZwallCol = NA, YZwallCol = NA)

########
#Fungus

library('dagitty')

fungusDAG= dagitty('dag {H0 -> H1 <- M -> F <- T}')
coordinates(fungusDAG)= list(x=c(H0=0,H1=1,M=2,F=3,T=4),
                             y=c(H0=0,H1=0,M=0,F=0,T=0))
drawdag(fungusDAG)

#simulate data
set.seed(71)
N = 1000
H0 = rnorm(n = N, mean = 10, sd = 2)
M = rbinom(n = N, size = 1, prob = 0.5)
T = rep(c(0, 1), each = N / 2)
H1 = H0 + rnorm(n = N, mean = 5 + 3 * M, sd = 1)
F = rbinom(n = N,
           size = 1,
           prob = 0.5 - 0.4 * T + 0.4 * M)
d2= data.frame(H0=H0,H1=H1,M=M,F=F,T=T)

#plot data
scatter3D_fancy(x=d2$H0,y=d2$H1,z=d2$F, bty='b2',theta=30,phi=30,
                xlim=c(0,20),ylim=c(0,20))

#fit model
mFung = quap(
  alist(
    H1 ~ dnorm(mu, sigma),
    mu <- H0 * p,
    p <- a + bF * F + bT * T,
    sigma ~ dexp(1),
    a ~ dlnorm(0, .2),
    bF ~ dnorm(0, 5),
    bT ~ dnorm(0, .5)
  ),
  data = d2
)
precis(mFung)

# simulate from posterior according to fitting model
post = extract.samples(mFung, n = 1000)
p = post$a + post$bF * d2$F + post$bT * d2$T
mu = d2$H0 * p
H1sim = rnorm(n = 1000, mean = mu, sd = post$sigma)
d2$H1sim = H1sim

ggplot(data = d2, aes(x = as.factor(T), y = H1sim)) +
  geom_boxplot()

#fit treatment-only model
mFungT = quap(alist(
  H1 ~ dnorm(mu, sigma),
  mu <- H0 * p,
  p <- a +  bT * T,
  sigma ~ dexp(1),
  a ~ dlnorm(0, .2),
  bT ~ dnorm(0, .5)
),
data = d2)
precis(mFungT)

################
# Chapter 6

# 6H1
library('rethinking')
data(WaffleDivorce)
d6H1 = WaffleDivorce
d6H1$South = d6H1$South + 1
head(d6H1)

#scale
d6H1$D_scale = (d6H1$Divorce - mean(d6H1$Divorce)) / sd(d6H1$Divorce)
d6H1$W_scale = (d6H1$WaffleHouses - mean(d6H1$WaffleHouses)) / sd(d6H1$WaffleHouses)

WH_mod = quap(
  alist(
    D_scale ~ dnorm(mu, sigma),
    mu <- a + bW[South] * W_scale + bS[South],
    
    #priors
    a ~ dnorm(0, 1),
    bW[South] ~ dnorm(0, .5),
    bS[South] ~ dnorm(0, .5),
    sigma ~ dexp(1)
  ),
  data = d6H1
)

priorSamp = extract.prior(WH_mod, n = 50)  #size matched to data
d6H1$mu_prior_sim = rep(NA, dim(d6H1)[1])
for (i in 1:dim(d6H1)[1]) {
  d6H1$mu_prior_sim[i] =
    priorSamp$a[i] + priorSamp$bW[i, d6H1$South[i]] * d6H1$W_scale[i] + priorSamp$bS[i, d6H1$South[i]]
}

colPal = c('red', 'blue')
plot(d6H1$W_scale,d6H1$mu_prior_sim, pch=16, cex=.5, col=colPal[d6H1$South],
     xlab='Scaled Waffle Houses', ylab='Scaled Divorce Rate')

plot(type='n',x=0,y=0,xlim=c(-2,5),ylim=c(-4,4),
     xlab='Scaled Waffle Houses', ylab='Scaled Divorce Rate', main='Prior Means')
for(i in 1:50){
  abline(a=priorSamp$a[i]+priorSamp$bS[i,d6H1$South[i]],
         b=priorSamp$bW[i,d6H1$South[i]], col=colPal[d6H1$South[i]], lty=1)
}


precis(WH_mod, depth = 2)
plot(precis(WH_mod,depth=2))
postSamp= extract.samples(WH_mod, n= 1000)
hist(postSamp$bW[,1],main='bW (South=0)',xlab='bW',col= rgb(1,0,0,.5))
hist(postSamp$bW[,2],main='bW (South=1)',xlab='bW',col=rgb(0,0,1,.5),add=TRUE)

hist( postSamp$bS[,2] - postSamp$bS[,1] )
rethinking::HPDI(postSamp$bS[, 2] - postSamp$bS[, 1])

###-----------------########
# Foxes
library('rethinking')
library('dagitty')

data(foxes)
head(foxes)

dag = dagitty('dag{A -> F -> G -> W;
             F -> W}')
drawdag(dag)

#rescale
foxes$F_scaled = (foxes$avgfood - mean(foxes$avgfood)) / sd(foxes$avgfood)
foxes$A_scaled = (foxes$area - mean(foxes$area)) / sd(foxes$area)
foxes$G_scaled = (foxes$group - mean(foxes$group)) / sd(foxes$group)
foxes$W_scaled = (foxes$weight - mean(foxes$weight)) / sd(foxes$weight)

#a Total causal influence of A on F
foxM1 = quap(alist(
  F_scaled ~ dnorm(mu, sigma),
  mu <- a + bA * A_scaled,
  a ~ dnorm(0, 1),
  bA ~ dnorm(0, 1),
  sigma ~ dexp(1)
),
data = foxes)

precis(foxM1)
plot(precis(foxM1))

#transform to original scale
a_postMean= precis(foxM1)$mean[1]
bA_postMean= precis(foxM1)$mean[2]
( a_postMean_orig= a_postMean * sd(foxes$avgfood) -
  bA_postMean*sd(foxes$avgfood)/sd(foxes$area)*mean(foxes$area)+
  mean(foxes$avgfood) )
( bA_postMean_orig= bA_postMean*sd(foxes$avgfood)/sd(foxes$area) )

# b, Total causal influence of F on W
foxM2 = quap(
  alist(
    W_scaled ~ dnorm(mu, sigma),
    mu <- a + bFW * F_scaled,
    a ~ dnorm(0, 5),
    bFW ~ dnorm(0, 0.2),
    sigma ~ dexp(1)
  ),
  data = foxes
)

precis(foxM2)

#simulate from posterior
postSamp = extract.samples(foxM2, n = 50)
F_scaled_vec = seq(-2, 2, .1)
W_sim = matrix(NA, nrow = length(F_scaled_vec), ncol = 50)
for (i in 1:50) {
  W_sim[, i] = postSamp$a[i] + postSamp$bFW[i] * F_scaled_vec
}
plot(NA, xlim=c(-2,2), ylim=c(-2,2), xlab='F_scaled', ylab='W_scaled')
for(i in 1:50) {
  lines(F_scaled_vec, W_sim[, i], col = rgb(0, 0, 0, .1))
}

#via link-function
W_link = link(foxM2, data = list(F_scaled = F_scaled_vec), n = 50)
plot(NA, xlim=c(-2,2), ylim=c(-2,2), xlab='F_scaled', ylab='W_scaled')
lines(F_scaled_vec, apply(W_link,2,mean), col=rgb(1,0,0,1),lwd=3)
shade(apply(W_link, 2, HPDI), F_scaled_vec, col = rgb(1, 0, 0, .5))

# c, Direct causal influence of F on W
foxM3 = quap(
  alist(
    W_scaled ~ dnorm(mu, sigma),
    mu <- a + bFW * F_scaled + bGW * G_scaled,
    a ~ dnorm(0, 5),
    c(bFW, bGW) ~ dnorm(0, 0.2),
    sigma ~ dexp(1)
  ),
  data = foxes
)

precis(foxM3)

####################
# Chapter 7
rm(list = ls())
cat('\014')
graphics.off()

library('rethinking')
library('reshape2')
library('tidyverse')

# simulate in- and out-of-sample accuracy
# simulate training and test data
set.seed(09112023)
a1 = 0.15
a2 = -0.4
sd = 0.01
plot(curve(a1 * x + a2 * x ^ 2, 0, 1))
x20 = runif(20, 0, 1)
x100 = runif(100, 0, 1)

y20_train = rnorm(20, mean = a1 * x20 + a2 * x20 ^ 2, sd = sd)
y20_test = rnorm(20, mean = a1 * x20 + a2 * x20 ^ 2, sd = sd)
y100_train = rnorm(100, mean = a1 * x100 + a2 * x100 ^ 2, sd = sd)
y100_test = rnorm(100, mean = a1 * x100 + a2 * x100 ^ 2, sd = sd)

# plot simulated data
#par(mfrow=c(1,2))
plot(x20, y20_train, pch=16, col=rgb(0,0,0,.5), ylim=c(-.2,.2),main='N=20')
points(x20, y20_test, pch=16, col=rgb(1,0,0,.5))
curve(a1*x + a2*x^2, add=TRUE, col=rgb(0,1,0,.5), lwd=3)
# plot(x100, y100_train, pch=16, col=rgb(0,0,0,.5), ylim=c(-2,2),main='N=100')
# points(x100, y100_test, pch=16, col=rgb(1,0,0,.5))
# curve(a1*x + a2*x^2, add=TRUE, col=rgb(0,1,0,.5), lwd=3)

# fit models
trainFit20_deg2 = quap(
  alist(
    y ~ dnorm(mu, sigma),
    mu <- a1 * x + a2 * x ^ 2,
    a1 ~ dnorm(0, 0.5),
    a2 ~ dnorm(0, 0.5),
    sigma ~ dexp(1)
  ),
  data = data.frame(x = x20, y = y20_train)
)

trainFit20_deg3 = quap(
  alist(
    y ~ dnorm(mu, sigma),
    mu <- a1 * x + a2 * x ^ 2 + a3 * x ^ 3,
    a1 ~ dnorm(0, 0.5),
    a2 ~ dnorm(0, 0.5),
    a3 ~ dnorm(0, 0.5),
    sigma ~ dexp(1)
  ),
  data = data.frame(x = x20, y = y20_train)
)

trainFit20_deg4 = quap(
  alist(
    y ~ dnorm(mu, sigma),
    mu <- a1 * x + a2 * x ^ 2 + a3 * x ^ 3 + a4 *
      x ^ 4,
    a1 ~ dnorm(0, 0.5),
    a2 ~ dnorm(0, 0.5),
    a3 ~ dnorm(0, 0.5),
    a4 ~ dnorm(0, 0.5),
    sigma ~ dexp(1)
  ),
  data = data.frame(x = x20, y = y20_train)
)

# Plot fits
meanPar_deg2 = precis(trainFit20_deg2)[, 1]
meanPar_deg3 = precis(trainFit20_deg3)[, 1]
meanPar_deg4 = precis(trainFit20_deg4)[, 1]
curve(meanPar_deg2[1]*x + meanPar_deg2[2]*x^2, col=rgb(0,0,0,.5), lwd=3,add=TRUE)
curve(meanPar_deg3[1]*x + meanPar_deg3[2]*x^2 + meanPar_deg3[3]*x^3, 
      col=rgb(0,0,0.3,.5), lwd=3,add=TRUE)
curve(meanPar_deg4[1]*x + meanPar_deg4[2]*x^2+ meanPar_deg4[3]*x^3 +
        meanPar_deg4[4]*x^4, col=rgb(0,0.3,0.3,.5), lwd=3,add=TRUE)

#evaluate posterior predictive distribution
getPostPred_deg2= function(fit, x_pred, y_pred){
  #get posterior predictive probabilities for given x_pred and y_pred
  postSamp= extract.samples(fit, n=100)
  log_postPred= apply(postSamp,1,FUN=function(coefs){
    sum(log(dnorm(y20_train, coefs[1]*x_pred + coefs[2]*x_pred^2, coefs[3])))})
  return(log_postPred)
}
postPred_train20_deg2 = getPostPred_deg2(trainFit20_deg2, x20, y20_train)

getPostPred_deg3= function(fit, x_pred, y_pred){
  #get posterior predictive probabilities for given x_pred and y_pred
  postSamp= extract.samples(fit, n=100)
  log_postPred= apply(postSamp,1,FUN=function(coefs){
    sum(log(dnorm(y20_train, 
          coefs[1]*x_pred + coefs[2]*x_pred^2 + coefs[3]*x_pred^3, coefs[4])))})
  return(log_postPred)
}
postPred_train20_deg3 = getPostPred_deg3(trainFit20_deg3, x20, y20_train)

getPostPred_deg4= function(fit, x_pred, y_pred){
  #get posterior predictive probabilities for given x_pred and y_pred
  postSamp= extract.samples(fit, n=100)
  log_postPred= apply(postSamp,1,FUN=function(coefs){
    sum(log(dnorm(y20_train, 
        coefs[1]*x_pred + coefs[2]*x_pred^2 + coefs[3]*x_pred^3 + coefs[4]*x_pred^4, coefs[5])))})
  return(log_postPred)
}
postPred_train20_deg4 = getPostPred_deg4(trainFit20_deg4, x20, y20_train)

#log probability DF
postPred_train20_DF = cbind.data.frame(deg2 = postPred_train20_deg2,
                                       deg3 = postPred_train20_deg3,
                                       deg4 = postPred_train20_deg4)
train20_means = data.frame(
  name = colnames(postPred_train20_DF),
  mean = apply(postPred_train20_DF, 2, mean),
  sd = apply(postPred_train20_DF, 2, sd)
)

# plot posterior predictive distribution
train20_melt = melt(postPred_train20_DF)
ggplot(data=train20_means)+
  geom_point(aes(x=name, y= mean), size=5)+
  geom_linerange(aes(x=name, ymin=mean-sd, ymax=mean+sd))+
  scale_y_continuous(name='log(probability)')


# probability DF
postPred_train20_DF = cbind.data.frame(
  deg2 = exp(postPred_train20_deg2),
  deg3 = exp(postPred_train20_deg3),
  deg4 = exp(postPred_train20_deg4)
)
train20_means = data.frame(
  name = colnames(postPred_train20_DF),
  mean = apply(postPred_train20_DF, 2, mean),
  sd = apply(postPred_train20_DF, 2, sd)
)

# plot posterior predictive distribution
train20_melt = melt(postPred_train20_DF)
ggplot(data=train20_means)+
  geom_point(aes(x=name, y= mean), size=5)+
  geom_linerange(aes(x=name, ymin=mean-sd, ymax=mean+sd))+
  scale_y_continuous(name='Probability')


#Deviance DF (-2*log(probability))
postPred_train20_DF = cbind.data.frame(
  deg2 = -2 * postPred_train20_deg2,
  deg3 = -2 * postPred_train20_deg3,
  deg4 = -2 * postPred_train20_deg4
)
train20_means = data.frame(
  name = colnames(postPred_train20_DF),
  mean = apply(postPred_train20_DF, 2, mean),
  sd = apply(postPred_train20_DF, 2, sd)
)

# plot posterior predictive distribution
train20_melt = melt(postPred_train20_DF)
ggplot(data=train20_means)+
  geom_point(aes(x=name, y= mean), size=5)+
  geom_linerange(aes(x=name, ymin=mean-sd, ymax=mean+sd))+
  scale_y_continuous(name='Deviance [-2*log(probability)]')

## More cross-validation, WAIC
library('rethinking')
data(cars)
head(cars)
#cars$speed= standardize(cars$speed)
#cars$dist= standardize(cars$dist)

carMod = quap(alist(
  dist ~ dnorm(mu, sigma),
  mu <- a + b * speed,
  a ~ dnorm(0, 100),
  b ~ dnorm(0, 10),
  sigma ~ dexp(1)
),
data = cars)

carPost = extract.samples(carMod, n = 1000)
#postPred= sapply(1:1000, function(i){ carPost$a[i] + carPost$b[i]*cars$speed})
#postLogProb= apply(postPred, 2, function(x){sum( log(dnorm(cars$dist, x, carPost$sigma))) })
logprob = sapply(1:1000, function(i) {
  mu = carPost$a[i] + carPost$b[i] * cars$speed
  return(log(dnorm(cars$dist, mu, carPost$sigma[i])))
})
ncases = dim(cars)[1]
#lppi= sapply(1:ncases, function(i){sum(logprob[i,])/1000})
lppi =  sapply(1:ncases, function(i) {
  log_sum_exp(logprob[i, ]) - log(1000)
})
pWAIC = sapply(1:ncases, function(i) {
  var(logprob[i, ])
})
(WAIC = -2 * (sum(lppi) - sum(pWAIC)))

WAIC(carMod)

########
# Week 3 - ex. 4
rm(list = ls())
graphics.off()
cat('\014')

library('dagitty')
library('rethinking')

data(foxes)
foxW = foxes
colnames(foxW) = c('group', 'F', 'G', 'A', 'W')
foxW$F = standardize(foxW$F)
foxW$G = standardize(foxW$G)
foxW$W = standardize(foxW$W)
foxW$A = standardize(foxW$A)

# DAG
dag= dagitty('dag{
             A -> F
             F <- U -> G
             G <- F -> W
             G -> W}')
coordinates(dag)= list(x= c(A=0, F=0, U=1, G=1, W=0.5), 
                       y= c(A=0, F=0.5, U=0, G=0.5, W=1))
plot(dag)

# total effect - cannot be estimated: backdoor path through U
totFox = quap(alist(
  W ~ dnorm(mu, sigma),
  mu <- a + b * F,
  a ~ dnorm(0, .2),
  b ~ dnorm(0, .5),
  sigma ~ dexp(1)
),
data = foxW)

precis(totFox)


# direct effect
dirFox = quap(
  alist(
    W ~ dnorm(mu, sigma),
    mu <- a + b_F * F + b_G * G,
    a ~ dnorm(0, .2),
    b_F ~ dnorm(0, .5),
    b_G ~ dnorm(0, .5),
    sigma ~ dexp(1)
  ),
  data = foxW
)

precis(dirFox)

# Week 4 - ex. 1
happy = sim_happiness(seed = 1977, N_years = 1000)
# happy$age= standardize(happy$age)
happy$married = happy$married + 1
# happy$happiness= standardize(happy$happiness)
happy = happy[happy$age >= 18, ] #only adults
happy$age = (happy$age - 18) / (65 - 18)  # rescale age to 0-1

# conditioned on collider
m6.9 = quap(
  alist(
    happiness ~ dnorm(mu, sigma),
    mu <- a[married] + b_AH * age,
    a[married] ~ dnorm(0, 1),
    b_AH ~ dnorm(0, 2),
    sigma ~ dexp(1)
  ),
  data = happy
)
precis(m6.9, depth = 2)

# not conditioned on collider
m6.10 = quap(alist(
  happiness <- dnorm(mu, sigma),
  mu <- a + b_AH * age,
  a ~ dnorm(0, 1),
  b_AH ~ dnorm(0, 2),
  sigma ~ dexp(1)
),
data = happy)
precis(m6.10)

WAIC(m6.9)
WAIC(m6.10)

PSIS(m6.9)
PSIS(m6.10)

compare(m6.9, m6.10, func = WAIC)
plot(compare(m6.9, m6.10, func = WAIC))

compare(m6.9, m6.10, func = PSIS)
plot(compare(m6.9, m6.10, func = PSIS))

# result: model that conditions on collider gives the lower WAIC / PSIS
# WAIC / PSIS focus on prediction, but do not help with causal structure

# Week 4 - ex. 2
#run again fox-models from last week

# total effect of area
areaFox = quap(alist(
  W ~ dnorm(mu, sigma),
  mu <- a + b_A * A,
  a ~ dnorm(0, .2),
  b_A ~ dnorm(0, .5),
  sigma ~ dexp(1)
),
data = foxW)

compare(totFox, dirFox, areaFox, func = WAIC)
plot(compare(totFox, dirFox, areaFox, func = WAIC))
#the direct effect model (conditioning on group size) is better (lower WAIC), but SEs overlap


compare(totFox, dirFox, areaFox, func = PSIS)
plot(compare(totFox, dirFox, areaFox, func = PSIS))
PSIS(areaFox, pointwise = TRUE)


precis(dirFox, depth = 2)
#causal interpretation: food availability increases fox weight; group size, which is correlated
#with food availability counters this effect through negative effect on weight

# Chapter 8 - terrain ruggedness vs GDP
data("rugged")
head(rugged)

rugFilt= rugged[,colnames(rugged) %in% c('cont_africa','cont_asia','cont_europe',
                                         'rugged', 'rgdppc_2000')]
rugFilt = rugFilt[apply(rugFilt, 1, function(x) {
  all(!is.na(x))
}), ]
rugFilt$logGDP = log(rugFilt$rgdppc_2000)
rugFilt$logGDP_std = rugFilt$logGDP / mean(rugFilt$logGDP)
rugFilt$rug_std = rugFilt$rugged / max(rugFilt$rugged)

m_rugInterac = quap(
  alist(
    logGDP_std ~ dnorm(mu, sigma),
    mu <-
      a + bR * (rug_std - 0.215) + bAf * cont_africa + bEu * cont_asia + bAs *
      cont_asia +
      bRaf * (rug_std - 0.215) * cont_africa +
      bReu * (rug_std - 0.215) * cont_europe +
      bRas * (rug_std - 0.215) * cont_asia,
    #using hardcoded mean of standardized ruggedness
    a ~ dnorm(1, .1),
    bR ~ dnorm(0, 1),
    c(bAf, bEu, bAs) ~ dnorm(0, 1),
    c(bRaf, bReu, bRas) ~ dnorm(0, .3),
    sigma <- dexp(1)
  ),
  data = rugFilt
)

precis(m_rugInterac, depth = 2)

head(rugFilt)
# add Africa-index (1= Africa, 2= not Africa)
rugFilt$CID = ifelse(rugFilt$cont_africa == 1, 1, 2) # CID = continent ID

# model that does not distinguish between Africa or not
m8.1 = quap(
  alist(
    logGDP_std ~ dnorm(mu, sigma),
    mu <- a + bR * (rug_std - 0.215),
    #using hardcoded mean of standardized ruggedness
    a ~ dnorm(1, .1),
    bR ~ dnorm(0, 1),
    sigma <- dexp(1)
  ),
  data = rugFilt
)


m8.2 = quap(
  alist(
    logGDP_std ~ dnorm(mu, sigma),
    mu <- a[CID] + bR * (rug_std - 0.215),
    #using hardcoded mean of standardized ruggedness
    a[CID] ~ dnorm(1, .1),
    bR ~ dnorm(0, 1),
    sigma <- dexp(1)
  ),
  data = rugFilt
)


m8.3 = quap(
  alist(
    logGDP_std ~ dnorm(mu, sigma),
    mu <- a[CID] + bR[CID] * (rug_std - 0.215),
    #using hardcoded mean of standardized ruggedness
    a[CID] ~ dnorm(1, .1),
    bR[CID] ~ dnorm(0, 1),
    sigma <- dexp(1)
  ),
  data = rugFilt
)

# compare models; check for influential datapoints
compare(m8.1, m8.2, m8.3, func = PSIS)
plot( PSIS(m8.1, pointwise = TRUE)$k )
plot( PSIS(m8.2, pointwise = TRUE)$k )
plot( PSIS(m8.3, pointwise = TRUE)$k )

rugVec= rep( seq(0,1,0.01), 2)
#CID_vec= rep( c(1,2), each= length(rugVec) )
postPred8.3_Af= link(m8.3, data= cbind.data.frame(rug_std= rugVec, CID= rep(1, length(rugVec))))
postPred8.3_nonAf= link(m8.3, data= cbind.data.frame(rug_std= rugVec, CID= rep(2, length(rugVec))))
Af_mu= apply(postPred8.3_Af, 2, mean)
AF_ci= apply(postPred8.3_Af, 2, PI, prob= 0.95)
nonAf_mu= apply(postPred8.3_nonAf, 2, mean)
nonAF_ci= apply(postPred8.3_nonAf, 2, PI, prob= 0.95)

cols= c('black', 'blue')
plot(rugFilt$rug_std, rugFilt$logGDP_std, pch= 16, col= cols[rugFilt$CID])
shade(AF_ci, rugVec)
lines(rugVec, Af_mu, col= 'black')
shade(nonAF_ci, rugVec, col= col.alpha(rangi2, .3))
lines(rugVec, nonAf_mu, col= 'blue')
