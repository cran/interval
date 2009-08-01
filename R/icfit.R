`icfit` <-
function(L,...){
    UseMethod("icfit")
}

`icfit.formula` <-
function (formula, data,...) 
{
    ## Most of this function is copied or slightly modified from survfit
    ## Copied starting from here:
    call <- match.call()
    if ((mode(call[[2]]) == "call" && call[[2]][[1]] == as.name("Surv")) || 
        inherits(formula, "Surv")) {
        formula <- eval(parse(text = paste(deparse(call[[2]]), 
            1, sep = "~")))
        environment(formula) <- parent.frame()
    }
    ## change survfit code to UseMethod "icfit" 
    if (!inherits(formula, "formula")) 
        temp <- UseMethod("icfit")
    else {
        m <- match.call(expand.dots = FALSE)
        m$... <- NULL
        Terms <- terms(formula, "strata")
        ord <- attr(Terms, "order")
        if (length(ord) & any(ord != 1)) 
            stop("Interaction terms are not valid for this function")
        m$formula <- Terms
        m[[1]] <- as.name("model.frame")
        m <- eval(m, parent.frame())
        n <- nrow(m)
        Y <- model.extract(m, "response")
        ## change survfit code next few lines
        if (!is.Surv(Y)){
            if (is.numeric(Y) & is.vector(Y)) Y<-Surv(Y,rep(1,length(Y))) 
            else stop("Response must be a survival object or numeric vector")
        }
        casewt <- model.extract(m, "weights")
        if (is.null(casewt)) 
            casewt <- rep(1, n)
        if (!is.null(attr(Terms, "offset"))) 
            warning("Offset term ignored")
        ll <- attr(Terms, "term.labels")
        if (length(ll) == 0) 
            X <- factor(rep(1, n))
        else X <- strata(m[ll])
        ### end of copied code from survfit.

        ### do a separate fit for each level of the factor
        group<-levels(X)
        nstrata<-length(group)
        sbind<-function(x,y){
            if (is.vector(x) & is.vector(y)) out<-c(x,y)
            if (is.matrix(x) & is.matrix(y)) out<-cbind(x,y) 
            if (!is.null(attr(x,"LRin")) & !is.null(attr(y,"LRin"))){
                xLRin<-attr(x,"LRin")
                yLRin<-attr(y,"LRin")
                attr(out,"LRin")<-cbind(xLRin,yLRin)
            } 
            return(out)
        }
        Y<- SurvLR(Y)
        for (i in 1:nstrata){
             tempout<-icfit.default(Y$L[X==group[i]],Y$R[X==group[i]],...)
             tempout$strata<-length(tempout$pf)
             names(tempout$strata)<-group[i]
             if (i==1){ icout<-tempout
             } else{
                 icout$A<-NULL
                 tempout$A<-NULL
                 icout<-mapply(sbind,icout,tempout)
             }
        }
        class(icout) <- c("icfit")
        if (!is.null(attr(m, "na.action"))) 
            icout$na.action <- attr(m, "na.action")
    }
    icout$call <- call
    icout
}




`Aintmap`<-function(L,R,Lin=NULL,Rin=NULL){
    n<-length(L)
    if (is.null(Lin) & is.null(Rin)){
        Lin<-rep(FALSE,n)
        Rin<-rep(TRUE,n)
        Lin[L==R]<-TRUE
        Rin[R==Inf]<-FALSE
    } else if (length(Lin)==1 & length(Rin)==1 & is.logical(Lin) & is.logical(Rin) ){
        Lin<-rep(Lin,n)
        Rin<-rep(Rin,n)
    } else if (length(Lin)!=n | length(Rin)!=n | !all(is.logical(Lin)) | !all(is.logical(Rin)) ){
        stop("Lin and Rin should be either NULL, logical length 1 or length same as L,R")
    } 
    if(n != length(R))
        stop("length of L and R must be the same")
    # calculate a small number, eps, to differentiate between e.g.,  [L,R] and (L,R]
    # we will treat, (L,R] as [L+eps,R], and [L,R) as [L,R-eps] 
    # since eps is only 
    # used in ranking we do not need to make it super small
    # just smaller than the smallest difference
    LRvalues<-sort(unique(c(0,L,R,Inf)))
    eps<- min(diff(LRvalues))/2
    Le<-L
    Re<-R
    Le[!Lin]<-L[!Lin]+eps
    Re[!Rin]<-R[!Rin]-eps
    # let s be the vector of ordered L and R values with 
    # R values later when there are ties
    # then intmap are values s[i] and s[i+1] where s[i] is 
    # associated with L and s[i+1] is associated with R
    oLR<-order(c(Le,Re+eps/2) )
    Leq1.Req2<-c(rep(1,n),rep(2,n))
    flag<- c(0,diff( Leq1.Req2[oLR] ))
    R.right.of.L<- (1:(2*n))[flag==1]
    intmapR<- c(L,R)[oLR][R.right.of.L]
    intmapL<- c(L,R)[oLR][R.right.of.L - 1]
    intmapRin<- c(Lin,Rin)[oLR][R.right.of.L]
    intmapLin<- c(Lin,Rin)[oLR][R.right.of.L - 1]
    intmap<-matrix(c(intmapL,intmapR),byrow=TRUE,nrow=2)
    attr(intmap,"LRin")<-matrix(c(intmapLin,intmapRin),byrow=TRUE,nrow=2)
    k<-dim(intmap)[[2]]
    Lbracket<-rep("(",k)
    Lbracket[intmapLin]<-"["
    Rbracket<-rep(")",k)
    Rbracket[intmapRin]<-"]"
    intname<-paste(Lbracket,intmapL,",",intmapR,Rbracket,sep="")
    A<-matrix(0,n,k,dimnames=list(1:n,intname))
    intmapLe<-intmapL
    intmapLe[!intmapLin]<-intmapL[!intmapLin]+eps
    intmapRe<-intmapR
    intmapRe[!intmapRin]<-intmapR[!intmapRin]-eps
    for (i in 1:n){
        tempint<- Le[i]<=intmapRe & Re[i]>=intmapLe
        A[i,tempint]<-1
    }
    out<-list(A=A,intmap=intmap)
    out
}




`icfit.default` <-
function(L, R, initfit = NULL, control=icfitControl(), Lin=NULL, Rin=NULL,...)
{
    epsilon<-control$epsilon
    maxit<-control$maxit
    AI<-Aintmap(L,R,Lin,Rin)
    A<-AI$A
    if (any(apply(A,1,sum)==0)) stop("A row all zeros. Appears that there are some R<L")
    n<-dim(A)[[1]]
    k<-dim(A)[[2]]
    intmap<-AI$intmap
    if (k>1){
##################################################################
### perform primary reductions on A  
### see Aragon and Eberly (1992) J of Computational and Graphical
###     Statistics 1:129-140 for discussion of primary reduction
##################################################################
	colsums <- apply(A, 2, sum)
	pairmult <- rep(0, k - 1)
	mark.to.keep <- rep(TRUE, k)
	for(i in 1:(k - 1)) {
		pairmult[i] <- sum(A[, i] * A[, i + 1])
		if(pairmult[i] == colsums[i]) {
			if(colsums[i] < colsums[i + 1]) {
				mark.to.keep[i] <- FALSE
			}
		}
		if(pairmult[i] == colsums[i + 1]) {
			if(colsums[i] >= colsums[i + 1]) {
				mark.to.keep[i + 1] <- FALSE
			}
		}
	}
	A <- A[, mark.to.keep]

      LRin<-attr(intmap,"LRin")[,mark.to.keep]
      intmap<-intmap[,mark.to.keep]
      attr(intmap,"LRin")<-LRin
     } # end primary reductions
### come up with the initial estimates
	if(is.null(initfit)) {
		pbar <- apply(A/apply(A, 1, sum), 2, mean)
	} else {
            if (is.null(initfit$pf) | is.null(initfit$intmap)) stop("initfit should be list with elements pf and intmap")
            nkeep<- dim(intmap)[[2]]
            pbar<-rep(0,nkeep)
            for (i in 1:nkeep){
                index<-initfit$intmap[1,]==intmap[1,i] & initfit$intmap[2,]==intmap[2,i]
                if (any(index)){
                    if (length(initfit$pf[index])>1) stop("initfit has non-unique intmap columns")
                    pbar[i]<- initfit$pf[index]
                }               
            }
            if (sum(pbar)==0) stop("initfit has no matching intmap elements with L and R")
            pbar<-pbar/sum(pbar)
	}
    em<-function(A,pbar,lower.bound=.01,startcount=1){
        converge<-FALSE 
        message<-"normal convergence"
        A.pbar<-as.vector(A %*% pbar)
        if (any(A.pbar==0)) stop("initfit$pf does not have nonzero values associated with needed intmap spaces")
        for (i in startcount:maxit){
 	      A.div.A.pbar <- A/A.pbar
            newpbar<- apply(t(A.div.A.pbar)*pbar,1,mean)
	      d <- apply(A.div.A.pbar, 2, sum)
	      u <-  - d + n
	      u[newpbar > 0] <- 0
	      error <- max(d + u - n)
            if (error<epsilon){
                pbar<-newpbar
                converge<-TRUE
                if(any(u < 0)){
		        message<-"Kuhn-Tucker conditions not met, self-consistent estimator not MLE"
                    pbar[u<0]<-min(pbar[pbar>0])
                    pbar<-pbar/sum(pbar)
                } 
               break()
            }
            newpbar[newpbar<lower.bound]<-0
            A.pbar<-as.vector(A %*% newpbar)
            if (any(A.pbar==0)){ 
                message<-"lower bound too high"
                break()}
            pbar<-newpbar
            if (i==maxit) message<-"maxit reached"
        }
        out<-list(A=A,pbar=pbar,count=i,message=message,converge=converge,error=error)
        out
    }

    ## try different values to "polish" the estimates, where if 
    ## the jump in the distribution is less than the lower.bound then 
    ## set it equal to zero. If it turns out that a jump should not 
    ## have been set to zero, the em function checks that and spits out 
    ## the last pbar before the jumps were set to zero 
    ## (see Gentleman and Geyer, 1994)
    lower.bounds<- 10^(0:ceiling(log10(epsilon)))
    lower.bounds<-lower.bounds[lower.bounds<=max(1/n,epsilon)]
    emout<-list(pbar=pbar,count=0)
    for (lb in lower.bounds){
        emout<-em(A,pbar=emout$pbar,lower.bound=lb,startcount=emout$count+1)

        if (emout$message=="normal convergence") break()
        if (emout$count==maxit) {
		emout$message<-"problem with convergence, increase maxit"
            break()
        }
        #print(emout$message)
    }
    keep<- !emout$pbar==0
    if (!all(keep)){
        LRin<-attr(intmap,"LRin")[,keep]
        intmap<-intmap[,keep]
        attr(intmap,"LRin")<-LRin
        A<-A[,keep]
    }
    pf<-emout$pbar[keep]
    strata<-length(pf)
    # if there is only one strata, title describes output:  NPMLE
    names(strata)<-"NPMLE"
    out <- list(A=A, strata=strata, error = emout$error, numit = emout$count, pf = pf, intmap = 
        intmap, converge= emout$converge, message= emout$message)
    class(out)<-c("icfit","list")
    return(out)
}

