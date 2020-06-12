oglmx<-function(formulaMEAN,formulaSD=NULL,data,start=NULL,weights=NULL,link="probit",constantMEAN=TRUE,constantSD=TRUE,beta=NULL,delta=NULL,threshparam=NULL,analhessian=TRUE,sdmodel=expression(exp(z)),SameModelMEANSD=FALSE,na.action,savemodelframe=FALSE,Force=FALSE,robust=FALSE){
  call<-match.call()
  names(call)[match("formulaMEAN",names(call),0)]<-"formula"
  m<-match(c("formula","data","subset","weights", "na.action", "offset"),names(call),0)
  mf<-call[c(1L,m)]
  mf$drop.unused.levels <- TRUE
  mf[[1L]]<- quote(stats::model.frame)
  mf<-eval(mf,parent.frame())
  weights<-as.vector(model.weights(mf))
  termsMEAN<-terms(mf)
  X<-model.matrix(formulaMEAN,mf)
  Y<-model.response(mf,"numeric")
  Keep<-!is.na(Y) & !apply(X,1,.IsNARow)
  if (!is.null(formulaSD) & !SameModelMEANSD){
    dataframeSD<-model.frame(formulaSD,data=data,na.action=na.pass)
    Z<-model.matrix(formulaSD,dataframeSD)
    Keep<-Keep & !apply(Z,1,.IsNARow)
    Z<-Z[Keep, ,drop=FALSE]
    termsSD<-terms(dataframeSD)
    # If specified that there is no constant in the model remove it from the data frame.
    if (!constantSD & ncol(Z)>1){
      savecolnamesZ<-colnames(Z)[colnames(Z)!="(Intercept)"]
      Z<-Z[,colnames(Z)!="(Intercept)",drop=FALSE]
      colnames(Z)<-savecolnamesZ
    }
  } else {
    termsSD<-NULL
  }
  X<-X[Keep, ,drop=FALSE]
  Y<-Y[Keep]
  if (!is.null(weights)){
    weights<-weights[Keep]
    robust<-TRUE
    X<-X[!is.na(weights), ,drop=FALSE]
    Y<-Y[!is.na(weights)]
    if (!is.null(formulaSD) & !SameModelMEANSD){
      Z<-Z[!is.na(weights), ,drop=FALSE]
    }
    weights<-weights[!is.na(weights)]
    weighted<-TRUE
    weights<-weights*length(weights)/sum(weights)
  } else {
    weights<-rep_len(1,length(Y))
    weighted<-FALSE
  }
  NoVarModData<-data.frame(Y,weights)

  No.Obs<-length(Y)
  termsMODEL<-list(termsMEAN,termsSD)
  formulaMODEL<-list(formulaMEAN,formulaSD)
  # If specified that there is no constant in the model remove it from the data frame.
  if (!constantMEAN){
    savecolnames<-colnames(X)[colnames(X)!="(Intercept)"]
    X<-X[,colnames(X)!="(Intercept)",drop=FALSE]
    colnames(X)<-savecolnames
  }

  # collect variable means and check which variables are binary.
  XVarMeans<-apply(X,2,mean)
  XVarBinary<-apply(X,2,.checkbinary)

  Heteroskedastic<-TRUE # set model to heteroskedastic, check call and then switch if not.
  if (!is.null(formulaSD)){
    # will check if the formula for the mean is identical to the
    # the right hand side variables of a formula are accessible via terms in the attribute term.labels
    meaneqnames<-attr(terms(formulaMEAN),"term.labels")
    sdeqnames<-attr(terms(formulaSD),"term.labels")
    if (sum(is.na(match(meaneqnames,sdeqnames)))==sum(is.na(match(sdeqnames,meaneqnames))) & sum(is.na(match(sdeqnames,meaneqnames)))==0){
      if (constantSD==constantMEAN){
        SameModelMEANSD<-TRUE
        # collect the names and column numbers of variables that are in both the mean and variance equation
        meanandvarNAME<-colnames(X)
        meanandvarLOC<-c(1:ncol(X))
        meanandvarLOCZ<-meanandvarLOC<-meanandvarLOC[meanandvarNAME!="(Intercept)"]
        meanandvarNAME<-meanandvarNAME[meanandvarNAME!="(Intercept)"]
        BothMeanVar<-data.frame(meanandvarNAME,meanandvarLOC,meanandvarLOCZ,stringsAsFactors=FALSE)
        ZVarMeans<-XVarMeans
        ZVarBinary<-XVarBinary
      }
    }
  }

  if (!is.null(formulaSD) & !SameModelMEANSD){
    if (!constantSD){Z<-Z[,colnames(Z)!="(Intercept)",drop=FALSE]}
    ZVarMeans<-apply(Z,2,mean)
    ZVarBinary<-apply(Z,2,.checkbinary)
    # collect the names and column numbers of variables that are in both the mean and variance equation
    # find the colnames of Z that are the same as the colnames of X
    meanandvarLOC<-c(1:ncol(X))[!is.na(match(colnames(X),colnames(Z)))]
    # find the columns of Z that are in Z and X
    meanandvarLOCZ<-match(colnames(X),colnames(Z))[!is.na(match(colnames(X),colnames(Z)))]
    meanandvarNAME<-colnames(X)[meanandvarLOC]
    meanandvarLOC<-meanandvarLOC[meanandvarNAME!="(Intercept)"]
    meanandvarLOCZ<-meanandvarLOCZ[meanandvarNAME!="(Intercept)"]
    meanandvarNAME<-meanandvarNAME[meanandvarNAME!="(Intercept)"]
    BothMeanVar<-data.frame(meanandvarNAME,meanandvarLOC,meanandvarLOCZ,stringsAsFactors=FALSE)
  } else if (is.null(formulaSD) & !SameModelMEANSD){
    Z<-as.matrix(rep(1,nrow(X)),ncol=1)
    ZVarMeans<-NULL
    ZVarBinary<-NULL
    Heteroskedastic<-FALSE
    if (is.null(delta) & is.null(threshparam)){
      # if no formula for the standard deviation is given, the threshold parameters are not specified or the standard deviation is not specified then use the unit variance assumption.
      calcdelta<-function(x){eval({z<-x;sdmodel})-1}
      delta<-uniroot(calcdelta,c(-10,10),extendInt="yes",tol=.Machine$double.eps)$root # solve for delta to get the unit variance assumption
    }
    BothMeanVar=NULL
  }


  checkoutcomes<-.checkoutcomes(Y,Force=Force)
  listoutcomes<-checkoutcomes[[1]]
  no.outcomes<-checkoutcomes[[2]]
  #return(list(Y,X,X,weights,link,sdmodel,beta,delta,threshparam,analhessian,robust,start))
  output<-oglmx.fit(Y,X,Z,w=weights,link = link,sdmodel = sdmodel,beta=beta,delta=delta,threshparam=threshparam,analhessian=analhessian,robustmatrix=robust,start=start,savemodelframe=savemodelframe)
  output<-append(output,list(call=call,terms=termsMODEL,formula=formulaMODEL,NoVarModData=NoVarModData,Hetero=Heteroskedastic,BothEq=BothMeanVar,varMeans=list(XVarMeans,ZVarMeans),varBinary=list(XVarBinary,ZVarBinary)))
  class(output)<-"oglmx"
  output
}
