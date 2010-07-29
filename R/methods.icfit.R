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
    if (!is.null(object$converse) && !all(object$converge)){
        if (length(object$converge)>1) warning("fit did not converge for at least one strata")
        else warning("fit did not converge")
    }

    intmap<-object$intmap
    k<-dim(intmap)[[2]]

    # make interval description match LRin attributes
    # ( and ) denote LRin value is FALSE
    # [ and ] denote LRin value is TRUE
    LRin<- attr(intmap,"LRin")
    # default to [ and ] if LRin is not given
    if (is.null(LRin)) LRin<-matrix(TRUE,2,k)
    Lbracket<-rep("(",k)
    Lbracket[LRin[1,]]<-"["
    Rbracket<-rep(")",k)
    Rbracket[LRin[2,]]<-"]"
    intname<-paste(Lbracket,round(intmap[1,],digits),",",round(intmap[2,],digits),Rbracket,sep="")

    tab<-data.frame(Interval=intname,Probability=round(object$pf,digits))
    # if more than one strata, then print each strata in order
    if (!is.null(object$strata) && length(object$strata)>1){
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
function(x,XLAB="time",YLAB=NULL,COL=gray((8:1)*.1),LTY=1:9,LEGEND=NULL,
    XLEG=NULL,YLEG=NULL,shade=TRUE,dtype="survival",
    dlink=function(x){log(-log(1-x))},...){

   ## time is union of 0 and all left and right endpoints
   time<-c(0,as.vector(x$intmap)) 
   
   ## dtype determines what y-axis represents, so make labels accordingly
   if (is.null(YLAB)){
       YLAB<-switch(dtype,
           survival="Survival",
           cdf="Cumulative Distribution",
           link="Transformed Distribution")
   }

   ## shading not supported for dtype='link'   
   if (dtype=="link") shade<-FALSE

   ## for dtype='link' need to calculate range of y for plot
   calc.ylim<-function(x){
        if (dtype=="survival" | dtype=="cdf"){ YLIM<-c(0,1)
        } else {
            if (dtype!="link") stop("dtype must be 'survival', 'cdf' or 'link' ")
            nstrata<-length(x$strata)
            ymin<-ymax<-NA
            ## calculate range of y for each strata and combine ranges
            for (i in 1:nstrata){
                H<- cumsum(x[i]$pf)
                H<-H[H>0 & H<1]
                ylim<- range( dlink(H) )
                if (i==1){
                    YLIM<-ylim
                } else {
                    YLIM<-c(min(ylim[1],YLIM[1]),max(ylim[2],YLIM[2]))
                }
            }
        }
        YLIM
    }

    YLIM<-calc.ylim(x)
    plot(range(time[time!=Inf]),YLIM,type="n",xlab=XLAB,ylab=YLAB,...)


   ## lines.icfit puts in lines 
    lines.icfit<-function(x,LTY){
        time<-c(0,as.vector(x$intmap))
        S<-c(1,1-cumsum(x$pf))
        S<-rep(S,each=2)[-2*length(S)]
        time[time==Inf]<-max(time)
        ## change y values depending on dtype
        if (dtype=="survival"){
            lines(time,S,lty=LTY)
        } else if (dtype=="cdf"){
            lines(time,1-S,lty=LTY)
        } else if (dtype=="link"){
            ## links will be inverse distributions, 
            ## so 0 and 1 give -Inf and +Inf, do not plot those
            I<-S>0 & S<1
            H<-1-S
            lines(time[I],dlink(H[I]),lty=LTY)
        }
    }

    ## function to do shading for dtype="survival" or dtype="cdf"
    polygon.icfit<-function(x,COL){
        S<-c(1,1-cumsum(x$pf))
        S<-rep(S,each=2)[-2*length(S)]
        time<-c(0,as.vector(x$intmap)) 
        if (any(time==Inf)){
            maxtime<-max(time[time<Inf])
            time[time==Inf]<-maxtime
            Satmax<-S[time==maxtime][1]
            ## change y values depending on dtype
            ## do not do shading for dtype="link"
            if (dtype=="survival"){
                polygon(c(maxtime,maxtime,2*maxtime,2*maxtime,maxtime),
                    c(Satmax,0,0,Satmax,Satmax),col=COL,border=NA)
            } else if (dtype=="cdf"){
                polygon(c(maxtime,maxtime,2*maxtime,2*maxtime,maxtime),
                    c(1-Satmax,1,1,1-Satmax,1-Satmax),col=COL,border=NA)
            } 
        }
        tt<-rep(time,each=2)
        tt<-c(tt[-1],tt[(length(tt)-1):1])
        SS<-rep(S,each=2)
        SS<-c(SS[-length(SS)],SS[length(SS):2]) 
        ## change y values depending on dtype
        ## do not do shading for dtype='link'
        if (dtype=="survival"){ polygon(tt,SS,col=COL,border=NA)
        } else if (dtype=="cdf"){ polygon(tt,1-SS,col=COL,border=NA) }
    }
    
    nstrata<-length(x$strata)
    if (nstrata==0) nstrata<-1

    if (nstrata>1){
        if (length(COL)<nstrata) COL<-rep(COL[1],nstrata)
        if (length(LTY)<nstrata) LTY<-rep(LTY[1],nstrata)
    
        ## remember shading not supported for dtype='link'    
        if (shade){
            for (i in 1:nstrata){polygon.icfit(x[i],COL[i])}
        }
       for (i in 1:nstrata){
            lines.icfit(x[i],LTY[i])
        }
    } else {
        ## remember shading not supported for dtype='link'    
        if (shade){ polygon.icfit(x,COL[1]) }
        lines.icfit(x,LTY[1])
    }

    ## xylegfunc is a function to find good place (maybe) to put legend
    ## if these are bad guesses then user should input XLEG and YLEG directly 
    xylegfunc<-function(x){
        if (dtype=="survival"){
            XLEG<-max(0,min(x$intmap[1,]))
            YLEG<-.1
        } else if (dtype=="cdf"){
            XLEG<-max(0,min(x$intmap[1,]))
            YLEG<-.9
        } else if (dtype=="link"){
            XLEG<-max(0,min(x$intmap[1,]))
            YLEG<-dlink(.9)
        }
        out<-list(x=XLEG,y=YLEG)
        out
    }
    xyleg<-xylegfunc(x)

    if (is.null(XLEG)) XLEG<- xyleg$x
    if (is.null(YLEG)) YLEG<- xyleg$y
    legend.list<-list(x=XLEG,y=YLEG,legend=names(x$strata),
        lty=LTY[1:nstrata],bty="n")
    if (shade) legend.list<-c(legend.list,list(fill=COL[1:nstrata]))

    if (is.null(LEGEND)){
        if (nstrata>1) do.call("legend",legend.list)
    } else if (LEGEND) do.call("legend",legend.list)
    ## return legend.list in case you want to first plot with LEGEND=FALSE, then put 
    ## the legend on manually using do.call such as above  except 
    ## after perhaps changing legend.list$x and legend.list$y
    invisible(legend.list)   
}



`[.icfit` <-
function(x,i){
    if (is.null(x$strata)){
        warning("no strata element in icfit object")
        out<-x
    } else{
        nstrata<-length(x$strata)
        # if you try to pick out the ith strata, but there are n<i strata, give error
        if (!any((1:nstrata)==i)){
            #warning(paste("number of strata=",nstrata,"but index=",i))
            stop("must use an integer to select a stratum")
        }
        pick.strata.part<-function(X,i,strata=x$strata){
            if (is.vector(X) && length(X)==length(strata)) part<-X[i]
            else if (is.vector(X) && length(X)==sum(strata)){
                part<-X[(sum(strata[0:(i-1)])+1):sum(strata[0:i])]
            } else if (is.matrix(X)){
                part<-X[,(sum(strata[0:(i-1)])+1):sum(strata[0:i])]
                if (!is.null(attr(X,"LRin"))){
                    attr(part,"LRin")<-attr(X,"LRin")[,(sum(strata[0:(i-1)])+1):sum(strata[0:i])]
                }
            } else part<-X
            
            return(part)
        }
        out<-mapply(pick.strata.part,x,i)
    }
    class(out)<-c("icfit")
    return(out)
}



