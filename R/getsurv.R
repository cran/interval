`getsurv` <-
function(times,icfit,nonUMLE.method="interpolation"){
    p<-icfit$pf
    S<-c(1-cumsum(p))
    L<-icfit$intmap[1,]
    R<-icfit$intmap[2,]
    k<-length(p)
    #print(cbind(L,R,p,S))
    
    ntimes<-length(times)
    Sout<-rep(NA,ntimes)
    mle<-rep(TRUE,ntimes)

    method<-match.arg(nonUMLE.method,c("interpolation","left","right"))
    nonUMLE.func<-switch(method,
                    interpolation=function(i,Time){
                        if (R[i]==Inf) stop("cannot interpolate when R[i]=Inf")
                        S[i-1] + ((Time - L[i])/(R[i]-L[i]))*(S[i]-S[i-1])},
                    left=function(i,Time){ S[i-1]},
                    right=function(i,Time){ S[i]})


    for (i in 1:ntimes){
        if (any(times[i]==R)) Sout[i]<-S[times[i]==R]
        else if (times[i]<=L[1]) Sout[i]<-1
        else if (times[i]>=R[k]) Sout[i]<-0
        else {
            if  (times[i]>L[k]){ 
                Sout[i]<-nonUMLE.func(k,times[i])
                mle[i]<-FALSE
            } else {
                iLup<-min((1:k)[L>=times[i]])
                if (R[iLup-1]<=times[i] | L[iLup]==times[i]) Sout[i]<-S[iLup-1]
                else {
                    Sout[i]<- nonUMLE.func(iLup-1,times[i])
                    mle[i]<-FALSE
                }
            }
       }
    }

    out<-list(S=Sout,times=times,unique.mle=mle,nonUMLE.method=method)
    out
}

