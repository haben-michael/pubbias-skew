## Figure. Visualize effect of selection on residuals.

set.seed(1)
alpha <- .3
n <- 7e2
v.bw <- 0
v <- runif(n,1,5)
y <- rnorm(n,sd=sqrt(v+v.bw))
keep.idx <- y/sqrt(v) > qnorm(1-alpha/2) ## p-val thresholding
data.sets <- list(preselect=list(y=y,v=v),postselect=list(y=y[keep.idx],v=v[keep.idx]))
file.name <- NULL
file.name <- './figs/heuristic.tex'
if(!is.null(file.name))tikzDevice::tikz(file.name, width = 6.5, height = 3)        
op <- par(mfrow=c(1,2))
invisible(lapply(data.sets, function(data.set) {
    with(data.set, {
        lm0 <- lm(I((y-mean(y))/sqrt(v+v.bw)) ~ I(1/sqrt(v+v.bw)))
        theta.fe <- sum(y*v)/sum(v)
        lm0 <- lm(I((y-theta.fe)/sqrt(v)) ~ I(1/sqrt(v)))
        hist(resid(lm0),breaks=30,xlab="residuals",freq=FALSE,main=NULL)
    })
}))
par(op)
if(!is.null(file.name))dev.off()





## Simulation. Gaussian model.
set.seed(1)
source('utils.R')
rZ <- rnorm #1
rVar <- function(n)runif(n,1,4)^2
start <- Sys.time()
ns <- 5:100
p.vals <- parallel::mclapply(ns, mc.cores=4, FUN=function(n) {
    print(n)
    B <- round(n^2)
    replicate(B, {
        z <- rZ(n)
        s <- 1/sqrt(rVar(n))
        y <- z/s
        v <- 1/s^2
        c(Egger=egger.test(y=y,v=v)$p.value,
          skew=skew.test(y=y,v=v)$p.value)
    })
})
print(Sys.time() - start)
alpha <- .1
p.vals.egger <- sapply(p.vals,function(x)x['Egger',])
p.vals.skew <- sapply(p.vals,function(x)x['skew',])
by.test <- list(Egger=p.vals.egger,skew=p.vals.skew)

file.name <- NULL
file.name <- './figs/gaussian.tex'
if(!is.null(file.name))tikzDevice::tikz(file.name, width = 6.5, height = 3)        
edgeworth <- list('Egger'=function(q,n)pnorm(q) - dnorm(q)*1/4*(q^3+q)/n, 'skew'=function(q,n)pnorm(q) - dnorm(q)*(-9*q/2/n +5/6/n*(q^3-3*q) )  )
op <- par(mfrow=c(1,2))
for(test.name in names(by.test)) {
    p.vals.test <- by.test[[test.name]]
    fpr <- sapply(p.vals.test, function(x)mean(x<alpha))
    cdf.ew <- edgeworth[[test.name]]
    fpr.error.ew <- function(n)abs(cdf.ew(qnorm(alpha/2),n) + 1-cdf.ew(qnorm(1-alpha/2),n) - alpha)
    plot(ns,abs(fpr-alpha),xlab='number of studies',ylab=paste0('$|$observed - nominal FPR$|$'),main=paste0(test.name,' test'),cex=1/2)
    curve(fpr.error.ew,add=TRUE)
}
par(op)
if(!is.null(file.name))dev.off()




## Common simulation and plotting code for remaining simulations.

fpr.sim <- function(dgp,ns,n.to.B,alpha) {
    fprs <- lapply(ns, FUN=function(n) {
        print(n)
        B <- n.to.B(n)
        p.vals <- replicate(B, {
            study.data <- dgp(n)
            y <- study.data$y; v <- study.data$v
            c(Egger=egger.test(y=y,v=v)$p.value,
              skew=skew.test(y=y,v=v)$p.value)
        })
        rowMeans(p.vals<=alpha)
    })
    fprs <- simplify2array(fprs)
    return(fprs)
}
plot.fprs <- function(ns,fprs,param.name,param.levels,file.name=NULL,smooth=FALSE) {
    if(!is.null(file.name)) {
        tikzDevice::tikz(file.name, width = 6.5, height = 3)        
    }
    op <- par(mfrow=c(1,2),  mar = c(4, 4, 3, 1), mgp = c(2.2, 1, 0))
    for(test.name in c("Egger","skew")) {
        fpr.test <- sapply(fprs,function(mat)mat[test.name,])
        vals <- abs(fpr.test-alpha)
        if(smooth)vals <- apply(vals,2,function(y)predict(loess(y~ns)))
        matplot(ns,vals,pch=1,col=1,type='l',xlab='number of studies',ylab="$|$observed - nominal FPR$|$", main=paste(test.name,'test'))
    }
    legend('topright',lty=1:length(param.levels),legend=param.levels,title=param.name,cex=.8)
    par(op)
    if(!is.null(file.name)) dev.off()
}







## Simulation. Moment model.

set.seed(1)
ns <- 5:100
alpha <- .1
n.to.B <- function(n)n^2
skew.levels <- c(a=1,b=2,c=3)
mu <- 0
sigma.to.skew <- function(sigma)(exp(sigma^2)+2)*sqrt(exp(sigma^2)-1)
params <- sigmas <- sapply(skew.levels, function(skew.level) uniroot(function(sigma)sigma.to.skew(sigma)-skew.level,c(0,10))$root)
fprs <- lapply(sigmas, function(sigma) {
    dgp <- function(n) {
        ms <- rpois(n,lambda=20)
        y <- sapply(ms, function(m) mean(exp(rnorm(m,mean=0,sd=sigma))))
        v <- sqrt((exp(sigma^2)-1)*exp(2*mu+sigma^2)) / ms
        return(list(y=y,v=v))
    }
    fpr.sim(dgp,ns,n.to.B,alpha)
})
file.name <- './figs/moment.tex'
plot.fprs(ns,fprs,param.name="skew",param.levels=skew.levels,file.name=file.name,smooth=TRUE)



## Simulation. Misspecified model--SMDs.
set.seed(1)
source('utils.R')
ns <- 5:100
alpha <- .1
n.to.B <- function(n)1e3
skew.level <- 1
sigma.to.skew <- function(sigma)(exp(sigma^2)+2)*sqrt(exp(sigma^2)-1)
sigma <- uniroot(function(sigma)sigma.to.skew(sigma)-skew.level,c(0,10))$root
rW <- function(m)exp(rnorm(m,0,sd=sigma))
n.to.ms <- list(function(n)rpois(n,lambda=rep(20,n)),function(n)2+rpois(n,lambda=rep(n,n)),function(n)2+rpois(n,rep(n^2,n)))
n.to.ms <- rev(n.to.ms)
start <- Sys.time()
fprs <- lapply(n.to.ms, function(n.to.m) {
    dgp <- function(n) {
        ms <- n.to.m(n)+2
        w <- lapply(ms,rW)
        y <- sapply(w,mean)
        v <- sapply(w,var)/ms
        return(list(y=y,v=v))
    }
    fpr.sim(dgp,ns,n.to.B,alpha)
})
save.image('250817c.RData')
file.name <- './figs/misspecified-smd.tex'
print(Sys.time()-start)
plot.fprs(ns,fprs,param.name="order of $m_i$",param.levels=c("$O(1)$","$O(n)$","$O(n^2)$"),file.name=file.name,smooth=TRUE)






## Simulation. Misspecified model--log-odds--effect of mu = average
## log odds.
set.seed(1)
source('utils.R')
ns <- 20:100
alpha <- .1
n.to.B <- function(n)1e3
n.to.m <- function(n)rep(n^2,n)
delta <- .05
mus <- seq(0,.1,len=3)
mus <- rev(mus)
start <- Sys.time()
fprs <- lapply(mus, function(mu) {
    pi1 <- plogis(mu-delta/2)
    pi2 <- plogis(mu+delta/2)
    dgp <- function(n) {
        ms <- n.to.m(n)
        pi1.hats <- rbinom(n,ms,pi1)/ms
        pi2.hats <- rbinom(n,ms,pi2)/ms
        y <- log(pi1.hats/(1-pi1.hats)) - log(pi2.hats/(1-pi2.hats))
        v <- ( 1/pi1.hats + 1/(1-pi1.hats) + 1/(pi2.hats) + 1/(1-pi2.hats) )/ms
        if(any(is.infinite(y)))browser()
        return(list(y=y,v=v))
    }
    fpr.sim(dgp,ns,n.to.B,alpha)
})
file.name <- './figs/misspecified-logodds-mu.tex'
## file.name <- NULL
print(Sys.time()-start)
plot.fprs(ns,fprs,param.name="$\\mu$",param.levels=round(mus,2),file.name=file.name,smooth=TRUE)








