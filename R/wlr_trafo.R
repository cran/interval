`wlr_trafo` <- function(x,...){
    UseMethod("wlr_trafo")
}

`wlr_trafo.Surv`<-function(x,...){
    LR<-SurvLR(x)
    wlr_trafo.default(LR$L,R=LR$R,...)
}

`wlr_trafo.default` <-
function(x, R=NULL, 
    scores =c("logrank1","logrank2","wmw"), 
    icFIT=NULL,
    initfit=NULL, 
    control=icfitControl(),
    Lin=NULL,
    Rin=NULL,...){
    L<-x
    scores<-match.arg(scores)
    if (scores!="logrank1" & scores!="logrank2" & scores!="wmw") stop("scores must equal 'logrank1' or 'logrank2' or 'wmw' ")


    if (is.null(R)) R<-L
    if (is.null(icFIT)){ 
        icFIT<-icfit(L,R,initfit,control,Lin,Rin)
        if (icFIT$message!="normal convergence") warning("icFIT does not have normal convergence")   
    }  

    A<-icFIT$A
    k<-dim(A)[[2]]
    n<-dim(A)[[1]]
    if (length(icFIT$pf)!=k) stop("icFIT$pf not proper length")

    cc<-scoresFromFit(icFIT,scores)   
    cc
}

