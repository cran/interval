`scoresFromFit` <- function(icFIT, scores){
    A<-icFIT$A
    k<-dim(A)[[2]]
    n<-dim(A)[[1]]
    phat<-c(0,icFIT$pf)
    Shat<- 1-cumsum(phat)
    Stilde<- exp( - c(0,cumsum(phat[2:(k+1)]/Shat[1:k])) )
    lr2<- (1/phat[2:(k+1)])*c( 
        Shat[1:(k-1)]* log( Shat[1:(k-1)]) - Shat[2:(k)]* log( Shat[2:(k)]),
        Shat[k]* log( Shat[k]))
    lr1<- (1/phat[2:(k+1)])*( Shat[1:k]* log(Stilde[1:k]) - Shat[2:(k+1)]* log(Stilde[2:(k+1)]) )
    wmw<-Shat[1:k] + Shat[2:(k+1)] - 1
    ckstar<-switch(scores,
        logrank2=lr2,
        logrank1=lr1,
        wmw=wmw)
    p<-phat[2:(k+1)]
    tempfunc<-function(Arow){
        sum(Arow*p*ckstar)/sum(Arow*p)
    }     
    cc<- apply(A,1,tempfunc)
    cc
}

