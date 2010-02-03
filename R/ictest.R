
`ictest` <-
function (L, ...){
    UseMethod("ictest")
}

`ictest.formula`<-function (formula, data, subset, na.action,...) 
{
    call <- match.call()
    m <- match.call(expand.dots = FALSE)
    m[[1]] <- as.name("model.frame")
    m$... <- NULL
    m <- eval(m, parent.frame())
    Y <- model.extract(m, "response")
    nmf<-length(m)
    if (nmf!=2) stop("formula should be in the form y~x, where y may be a Surv object")
    if (!is.Surv(Y)){
        if (is.numeric(Y) & is.vector(Y)) Y<-Surv(Y,rep(1,length(Y))) 
        else stop("Response must be a survival object or numeric vector")
    }
    LR<-SurvLR(Y)
    group.name<-attr(attr(m,"terms"),"term.labels")
    if (!is.numeric(m[[2]])){
        group<-paste(group.name,"=",as.character(m[[2]]),sep="")
    } else { 
        group<-m[[2]] 
    }
    icout<-do.call("ictest",c(list(L=LR$L,R=LR$R,group=group),list(...)))
    icout$call <- call

    icout$data.name <- paste(names(m), collapse = " by ")
    icout
}

`ictest.default` <-
function(L, R, group,  
    scores = c("logrank1","logrank2","wmw","normal","general"),
    rho=NULL,
    alternative= c("two.sided", "less", "greater","two.sidedAbs"),   
    icFIT=NULL,
    initfit="initcomputeMLE", 
    icontrol=icfitControl(),
    exact=NULL,
    method=NULL,
    methodRule=methodRuleIC1,
    mcontrol=mControl(),
    Lin=NULL,
    Rin=NULL,
    dqfunc=NULL,...){
    
    ## translate rho into scores
    if (!is.null(rho)){
        if (length(scores)!=5){ warning("scores ignored because rho is non-NULL") }
        if (rho==0) scores<-"logrank1"
        else if (rho==1) scores<-"wmw"
        else stop("rho must be either NULL, 0 or 1")
    }
    
    alternative <- match.arg(alternative)
    scores <- match.arg(scores)
    ## find NPMLE based on all the data (ignoring group membership)
    ## unless icFIT is not null
    if (is.null(icFIT)){ 
        icFIT<-icfit(L,R, initfit, control=icontrol, Lin=Lin, Rin=Rin)
        if (icFIT$message!="normal convergence") warning("icFIT does not have normal convergence")   
    }  
    ## calculate scores
    cc<-wlr_trafo(L,R, scores=scores, icFIT=icFIT, 
        initfit=initfit, control=icontrol, Lin=Lin, Rin=Rin, dqfunc=dqfunc)


    ## if group is numeric but only two unique levels then treat as two sample case
    if (is.factor(group) | length(unique(group))==2) group<-as.character(group)

    ## for 2- or k-sample: calculate efficient score statistics, U, and sample size per group, N
    if (is.character(group)){
        ug<-unique(group)
        ng<-length(ug)
        if (ng==1) stop("group variable same for every individual")
        U<-rep(NA,ng)
        names(U)<-ug
        N<-U
        for (j in 1:ng){
            U[j]<-sum(cc[group==ug[j]])
            N[j]<-length(cc[group==ug[j]])
        }
    } else if (is.numeric(group)){
    ## for trend tests: U is still efficient score statistic, N is total sample size
        U<-sum(cc*group) 
        N<-length(cc)
        ng<-0
    } else { stop("group should be a factor, character, or numeric vector") }

    ## get method and check it is allowed
    if (is.null(method))    method<-methodRule(cc,group,exact)
    method.OK<-(method=="pclt" | method=="exact.mc" | method=="exact.network" | method=="exact.ce" | method=="scoretest" | method=="wsr.mc" | method=="wsr.HLY" | method=="wsr.pclt")
    if (!method.OK) stop("method not one of: 'pclt', 'exact.mc'. 'exact.network', 'exact.ce', 'scoretest', 'wsr.mc', 'wsr.HLY', 'wsr.pclt'")

    ## TEST describes results
    if (method=="exact.network" || method=="exact.ce" || method=="exact.mc" || method=="wsr.mc") TEST<-"Exact"
    else TEST<-"Asymptotic"
    if (scores=="logrank1" || scores=="logrank2") TEST<-paste(TEST,"Logrank")
    else if (scores=="wmw")  TEST<-paste(TEST,"Wilcoxon")
    else if (scores=="normal")  TEST<-paste(TEST,"Normal Scores")
    else if (scores=="general")  TEST<-paste(TEST,"General Scores")
  
    ## normal and general scores only work for permutation tests
    if ((scores=="normal" | scores=="general") & (method=="scoretest" | method=="wsr.mc" | 
          method=="wsr.HLY" | method=="wsr.pclt")) stop("normal or general scores only programmed for permutation methods")
 

    mcontrol<-c(list(exact=exact,method=method),mcontrol)
    
    ## if all scores equal 0, then no need to do any calculations, all p-values=1
    if (all(cc==0)){
        p.values<-c(p.twosided=1,p.lte=1,p.gte=1,p.twosidedAbs=1)
        pout<-list(p.value=1,p.values=p.values)
        TEST<-paste(TEST,"(all scores=0, score option irrelevant)")
        alt.phrase<-"(p-values for all alternatives equal 1)"
    } else {
    ## Next is main calculation section, ng=0 is for trend tests
    ## ng=2 is 2-sample tests
    ## ng>2 is k-sample tests
    if (ng==0){ 
        if (method=="scoretest"){
            pout<-icScoreTest(icFIT,group,scores,alternative,tol.svd=mcontrol$tol.svd)
            TEST<-paste(TEST,"trend test(score form)")
       } else if (method=="wsr.mc" | method=="wsr.HLY" | method=="wsr.pclt"){
           pout<-icWSR(icFIT,group,scores,alternative,type=method,control=mcontrol)
            if (method=="wsr.mc"){
                TEST<-paste(TEST,"trend test(WSR Monte Carlo)")
            } else if (method=="wsr.pclt"){
                TEST<-paste(TEST,"trend test(WSR pclt)")
            }
       } else {
            ## all other methods are permutation test methods using method options in permTREND
            pout <- do.call("permTREND", list(x=cc,y=group, alternative=alternative, 
                exact=exact,method=method,control=mcontrol))
            TEST<-paste(TEST,"trend test(permutation form)")
        }
        ## alt.phrase is used to describe alternative hypothesis
        alt.phrase<-switch(alternative,
            two.sided="survival distributions not equal",
            two.sidedAbs="survival distributions not equal",
            less=paste("higher independent variable implies earlier failure times than expected"),
            greater=paste("higher independent variable implies later failure times than expected"))
    } else if (ng==2){
        if (method=="scoretest"){
            pout<-icScoreTest(icFIT,group,scores,alternative,tol.svd=mcontrol$tol.svd)
            TEST<-paste(TEST,"two-sample test (score form)")
        } else if (method=="wsr.mc" | method=="wsr.HLY" | method=="wsr.pclt"){
            pout<-icWSR(icFIT,group,scores,alternative,type=method,control=mcontrol)
            if (method=="wsr.mc"){
                TEST<-paste(TEST,"2-sample test(WSR Monte Carlo)")
            } else if (method=="wsr.pclt"){
                TEST<-paste(TEST,"2-sample test(WSR pclt)")
            } else if (method=="wsr.HLY"){
                TEST<-paste(TEST,"2-sample test(WSR HLY)")
            }
       } else {
            ## all other methods are permutation test methods using method options in permTS
            X<-cc[group==ug[1]]
            Y<-cc[group==ug[2]]
            pout <- do.call("permTS", list(x=X,y=Y, alternative=alternative, 
                exact=exact,method=method, control=mcontrol))       
            TEST<-paste(TEST,"two-sample test (permutation form)")
        }
        ## alt.phrase is used to describe alternative hypothesis
        alt.phrase<-switch(alternative,
            two.sided="survival distributions not equal",
            two.sidedAbs="survival distributions not equal",
            less=paste(ug[2]," has earlier failure times than expected"),
            greater=paste(ug[2]," has later failure times than expected"))
    } else if (ng>2){
        if (method=="scoretest"){
            pout<-icScoreTest(icFIT,group,scores,alternative,tol.svd=mcontrol$tol.svd)
            TEST<-paste(TEST,"k-sample test (score form)")
        }  else if (method=="wsr.mc" | method=="wsr.HLY" | method=="wsr.pclt"){
            pout<-icWSR(icFIT,group,scores,alternative,type=method,control=mcontrol)

            if (method=="wsr.mc"){
                TEST<-paste(TEST,"2-sample test(WSR Monte Carlo)")
            } else if (method=="wsr.pclt"){
                TEST<-paste(TEST,"2-sample test(WSR pclt)")
            } else if (method=="wsr.HLY"){
                TEST<-paste(TEST,"2-sample test(WSR HLY)")
            }
       } else {
            ## all other methods are permutation test methods using method options in permKS
            pout <- do.call("permKS", list(x=cc,g=group,
                exact=exact,method=method, control=mcontrol))      
            TEST<-paste(TEST,"k-sample test (permutation form)")
        }
        alt.phrase<-"survival distributions not equal"
        if (alternative=="less" | alternative=="greater") warning("alternative ignored, group is factor with more than 2 groups")
    }
    if (scores=="logrank1") TEST<-paste(TEST,", Sun's scores",sep="")
    else if (scores=="logrank2") TEST<-paste(TEST,", Finkelstein's scores",sep="")

    } # end else associated with "if (all(cc==0))"

    call <- match.call()
    getName<-function(CALL,defaultName=""){
        Name<-as.character(CALL)
        if (length(Name)>1) Name<-defaultName
        Name
    }

    L.name<-getName(call[[2]],"L")
    R.name<-getName(call[[3]],"R")
    group.name<-getName(call[[4]],"group")
    pout$data.name<-paste("{",L.name,",",R.name,"}"," by ",group.name,sep="")

    pout$method<-TEST    
    pout$scores <- cc
    pout$U<-U
    pout$N<-N
    pout$algorithm<-method
    pout$alt.phrase<-alt.phrase
    pout$fit<-icFIT

    if (method=="wsr.HLY" | method=="wsr.pclt" | method=="wsr.mc") pout$nwsr<-mcontrol$nwsr
    if (method=="wsr.mc") pout$np<-mcontrol$np

    class(pout)<-"ictest"
    pout
}



