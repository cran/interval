`print.icfit`<-function(x,...){
    temp<-x
    if (!is.null(temp$A)){
        n<-dim(temp$A)[[1]]
        k<-dim(temp$A)[[2]]
        temp$A<-paste("A matrix (",n," by ",k,") not printed but part of the list")
    }
    class(temp)<-"list"
    print(temp)
}



`summary.icfit`<-function(object,digits=4,...){
    
    if (!all(object$converge)){
        if (length(object$converge)>1) warning("fit did not converge for at least one strata")
        else warning("fit did not converge")
    }

    intmap<-object$intmap
    k<-dim(intmap)[[2]]

    LRin<- attr(intmap,"LRin")
    if (is.null(LRin)) LRin<-matrix(TRUE,2,k)
    Lbracket<-rep("(",k)
    Lbracket[LRin[1,]]<-"["
    Rbracket<-rep(")",k)
    Rbracket[LRin[2,]]<-"]"
    intname<-paste(Lbracket,round(intmap[1,],digits),",",round(intmap[2,],digits),Rbracket,sep="")

    tab<-data.frame(Interval=intname,Probability=round(object$pf,digits))
    if (!is.null(object$strata) & length(object$strata)>1){
        cnt<-1
        for (i in 1:length(object$strata)){
            cat(paste(names(object$strata)[i],":",sep=""))
            endCnt<-cnt+object$strata[i]-1
            cat("\n")
            tabi<-tab[cnt:endCnt,]
            dimnames(tabi)[[1]]<-1:object$strata[i]
            print(tabi,digits=digits)
            cnt<-endCnt+1
        }
    } else {
        dimnames(tab)[[1]]<- 1:k 
        print(tab,digits=digits) 
    }
}

`plot.icfit` <-
function(x,XLAB="time",YLAB="Survival",COL=gray((8:1)*.1),LTY=1:9,LEGEND=NULL,
    XLEG=NULL,YLEG=.1,...){
   time<-c(0,as.vector(x$intmap)) 
   plot(range(time[time!=Inf]),c(0,1),type="n",xlab=XLAB,ylab=YLAB,...)
    lines.icfit<-function(x,LTY){
        time<-c(0,as.vector(x$intmap))
        S<-c(1,1-cumsum(x$pf))
        S<-rep(S,each=2)[-2*length(S)]
        time[time==Inf]<-max(time)
        lines(time,S,lty=LTY)
    }
    polygon.icfit<-function(x,COL){
        S<-c(1,1-cumsum(x$pf))
        S<-rep(S,each=2)[-2*length(S)]
        time<-c(0,as.vector(x$intmap)) 
        if (any(time==Inf)){
            maxtime<-max(time[time<Inf])
            time[time==Inf]<-maxtime
            Satmax<-S[time==maxtime][1]
            polygon(c(maxtime,maxtime,2*maxtime,2*maxtime,maxtime),
                c(Satmax,0,0,Satmax,Satmax),col=COL,border=NA)
        }
        tt<-rep(time,each=2)
        tt<-c(tt[-1],tt[(length(tt)-1):1])
        SS<-rep(S,each=2)
        SS<-c(SS[-length(SS)],SS[length(SS):2]) 
        polygon(tt,SS,col=COL,border=NA)
    }
    
    nstrata<-length(x$strata)
    if (nstrata==0) nstrata<-1
    if (nstrata>1){
        if (length(COL)<nstrata) COL<-rep(COL[1],nstrata)
        if (length(LTY)<nstrata) LTY<-rep(LTY[1],nstrata)
        
        for (i in 1:nstrata){
            polygon.icfit(x[i],COL[i])
        }
       for (i in 1:nstrata){
            lines.icfit(x[i],LTY[i])
        }
    } else {
        polygon.icfit(x,COL[1])
        lines.icfit(x,LTY[1])
    }
    if (is.null(XLEG)) XLEG<- max(0,min(x$intmap[1,]))
    if (is.null(YLEG)) YLEG<- .1
    legend.list<-list(x=XLEG,y=YLEG,legend=names(x$strata),
        fill=COL[1:nstrata],lty=LTY[1:nstrata],bty="n")
    if (is.null(LEGEND)){
        if (nstrata>1) do.call("legend",legend.list)
    } else if (LEGEND) do.call("legend",legend.list)
    invisible(legend.list)   
}


`[.icfit` <-
function(x,i){
    if (is.null(x$strata)){
        warning("no strata element in icfit object")
        out<-x
    } else{
        nstrata<-length(x$strata)
        if (!any((1:nstrata)==i)){
            warning(paste("number of strata=",nstrata,"but index=",i))
        }
        pick.strata.part<-function(X,i,strata=x$strata){
            if (is.vector(X) && length(X)==length(strata)) part<-X[i]
            else if (is.vector(X) && length(X)==sum(strata)){
                part<-X[(sum(strata[0:(i-1)])+1):sum(strata[0:i])]
            } else if (is.matrix(X)){
                part<-X[,(sum(strata[0:(i-1)])+1):sum(strata[0:i])]
            } else part<-X
            
            return(part)
        }
        out<-mapply(pick.strata.part,x,i)
    }
    class(out)<-c("icfit")
    return(out)
}



