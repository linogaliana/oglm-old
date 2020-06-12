oglmx.fit<-function(Y,X,Z=NULL,w,link="probit",sdmodel=expression(exp(z)),SameModelMEANSD=FALSE,beta=NULL,delta=NULL,threshparam=NULL,analhessian=TRUE,robustmatrix=FALSE,start=NULL,savemodelframe=FALSE){
  Xr<-split.data.frame(X,Y,drop=FALSE)
  no.outcomes<-length(Xr)
  No.Obs<-length(Y)
  listoutcomes<- as.numeric(levels(as.factor(Y)))[order(as.numeric(levels(as.factor(Y))))]
  weightsr<-split(w,Y,drop=FALSE)

  if (is.null(Z) & SameModelMEANSD){
    Zr<-Xr
  } else if (is.null(Z) & !SameModelMEANSD){
    Z<-matrix(rep(1,No.Obs),ncol=1)
    Zr<-split.data.frame(Z,Y,drop = FALSE)
  } else {
    Zr<-split.data.frame(Z,Y,drop=FALSE)
  }

  ProbFunc<-.cdf.func(link)
  ProbFuncD<-.pdf.func(link)
  ProbFuncDD<-.Dpdf.func(link)

  no.Xvar<-ncol(X)
  if (!SameModelMEANSD){no.Zvar<-ncol(Z)} else {no.Zvar<-no.Xvar}

  # set the prespecified and not prespecified parts.
  # beta
  if (!is.null(beta) & length(beta)==1){
    storebeta<-beta
    beta<-rep(NA,no.Xvar)
    beta[1]<-storebeta
  } else if (!is.null(beta) & length(beta)>1){
    # check that the specified vector is of correct length
    if (length(beta)!=no.Xvar){stop("Specified beta vector of incorrect length.")}
  } else if (is.null(beta)){
    beta<-rep(NA,no.Xvar)
  }
  # delta
  if (!is.null(delta) & length(delta)==1){
    storedelta<-delta
    delta<-rep(NA,no.Zvar)
    delta[1]<-storedelta
  } else if (!is.null(delta) & length(delta)>1){
    # check that the specified vector is of correct length
    if (length(delta)!=no.Zvar){stop("Specified delta vector of incorrect length.")}
  } else if (is.null(delta)){
    if (is.null(Z)){
      no.Zvar<-1
    }
    delta<-rep(NA,no.Zvar)
  }
  # threshparam
  if (!is.null(threshparam) & length(threshparam)==1){
    storethreshparam<-threshparam
    threshparam<-rep(NA,no.outcomes-1)
    threshparam[1]<-storethreshparam
  } else if (!is.null(threshparam) & length(threshparam)>1){
    # check that the specified vector is of correct length
    if (length(threshparam)!=no.outcomes-1){stop("Specified vector of threshold parameters of incorrect length.")}
  } else if (is.null(threshparam)){
    threshparam<-rep(NA,no.outcomes-1)
  }

  if (!is.null(beta)){
    collectbeta<-is.na(beta)
  } else {
    collectbeta<-!vector("logical",no.Xvar)
  }
  if (!is.null(delta)){
    collectdelta<-is.na(delta)
  } else {
    collectdelta<-!vector("logical",no.Zvar)
  }
  if (!is.null(threshparam)){
    collectthreshparam<-is.na(threshparam)
  } else {
    collectthreshparam<-!vector("logical",no.outcomes-1)
  }

  no.betaparams<-sum(collectbeta)
  no.deltaparams<-sum(collectdelta)
  no.threshparams<-sum(collectthreshparam)
  no.parameters<-no.betaparams+no.deltaparams+no.threshparams
  Est.Parameters<-list(beta=collectbeta,delta=collectdelta,alpha=collectthreshparam)

  # specify the start vector in the case that it is given as null, give error message if not of correct length.
  if (is.null(start)){
    # start with a vector of zeros for the betas
    start<-vector("numeric",no.parameters)
    # if none of the delta parameters are set, set the first element so that the initial standard deviation is 0.5
    calcstartdelta<-function(x){eval({z<-x;sdmodel})-0.5}
    startdelta<-uniroot(calcstartdelta,c(-10,10),extendInt="yes")$root
    if (no.deltaparams==length(delta)){start[no.betaparams+1]<-startdelta}
    # more complicated for threshparam, should respect the order of the prespecified values
    if (no.threshparams>0){
      cutoff<-1
      savethreshparam<-threshparam
      for (i in 1:length(threshparam)){
        if (collectthreshparam[i]){
          if (i==1){
            if (length(threshparam)>1){
              if (!all(collectthreshparam[2:length(threshparam)])){
                start[no.betaparams+no.deltaparams+cutoff]<-min(savethreshparam[2:length(threshparam)],na.rm = TRUE)-1
              } else {
                start[no.betaparams+no.deltaparams+cutoff]<-(listoutcomes[i]+listoutcomes[i+1])/2
              }
            } else {start[no.betaparams+no.deltaparams+cutoff]<-(listoutcomes[i]+listoutcomes[i+1])/2}
          } else {
            if (i<length(threshparam)){
              if (!all(collectthreshparam[(i+1):length(threshparam)])){
                right<-min(savethreshparam[(i+1):length(threshparam)],na.rm = TRUE)
                start[no.betaparams+no.deltaparams+cutoff]<-(savethreshparam[i-1]+right)/2
              } else {start[no.betaparams+no.deltaparams+cutoff]<-savethreshparam[i-1]+0.5}
            } else {start[no.betaparams+no.deltaparams+cutoff]<-savethreshparam[i-1]+0.5}
          }
          savethreshparam[i]<-start[no.betaparams+no.deltaparams+cutoff]
          cutoff<-cutoff+1
        }
      }
    }
  } else if (length(start)!=no.parameters){
    stop("Specified vector of start values for parameters of incorrect length.")
  }
  # create the vectors of names matched to parameter estimates
  collectXnames<-colnames(X)[collectbeta]
  if (!SameModelMEANSD){
    collectZnames<-colnames(Z)[collectdelta]
  } else {
    collectZnames<-colnames(X)[collectdelta]
  }

  outputnames<-vector("character",0)
  if (no.betaparams>0){
    outputnames<-c(outputnames,collectXnames)
  }
  if (no.deltaparams>0){
    outputnames<-c(outputnames,collectZnames)
  }
  if (no.threshparams>0){
    collectnumbers<-c(1:length(threshparam))[collectthreshparam]
    outputnames<-c(outputnames,sapply(collectnumbers,function(x){paste("Threshold (", listoutcomes[x],"->",listoutcomes[x+1],")",sep="")}))
  }

  # save model frame if requested
  if (savemodelframe){
    modelframes<-list(X)
    if (!SameModelMEANSD){
      modelframes[[2]]<-Z
    } else {
      modelframes[[2]]<-X
    }
  } else {
    modelframes<-NULL
  }


  # function that calculates log likelihood, gradient and hessian given a set of parameter values
  LLoglmx<-function(param,beta=NULL,delta=NULL,threshparam=NULL,analhessian=FALSE,robustmatrix=FALSE){
    if (is.null(beta)){ # if elements of the vector beta are not prespecified then all are included in param
      beta<-param[1:no.Xvar]
      param<-param[(no.Xvar+1):length(param)] # remove from param the elements allocated to the beta vector
    } else { # if not then fill NAs in beta with the first elements in param
      countNA<-sum(is.na(beta))
      if (countNA>0){beta[is.na(beta)]<-param[1:countNA]}
      param<-param[(countNA+1):length(param)]
    }
    # repeat to fill the delta vector
    if (is.null(delta)){ # if elements of the vector delta are not prespecified then all are included in param
      if (SameModelMEANSD){
        delta<-param[1:no.Xvar]
        param<-param[(no.Xvar+1):length(param)] # remove from param the elements allocated to the beta vector
      } else {
        delta<-param[1:no.Zvar]
        param<-param[(no.Zvar+1):length(param)] # remove from param the elements allocated to the beta vector
      }
    } else { # if not then fill NAs in delta with the first elements in param
      countNA<-sum(is.na(delta))
      if (countNA>0){delta[is.na(delta)]<-param[1:countNA]}
      param<-param[(countNA+1):length(param)]
    }
    # repeat to fill threshparam vector.
    if (is.null(threshparam)){
      threshparam<-param
    } else {
      if (length(param)>0){threshparam[is.na(threshparam)]<-param} else {stop("Insufficient number of parameters specified.")}
    }
    threshparam<-c(-Inf,threshparam,Inf)

    calcprobs<-function(outcome){
      # function that calculates relevant probabilities relevant to outcome
      # only to be called inside LLoglmx
      Xb<-Xr[[outcome]]%*%beta
      if (!SameModelMEANSD){
        Zdinv<-1/eval({z<-Zr[[outcome]]%*%delta;sdmodel})
      } else {
        Zdinv<-1/eval({z<-Xr[[outcome]]%*%delta;sdmodel})
      }
      Probs<-ProbFunc((threshparam[outcome+1]-Xb)*Zdinv)-ProbFunc((threshparam[outcome]-Xb)*Zdinv)
    }

    sdmodfirstderiv<-D(sdmodel,"z")
    sdmodsecondderiv<-D(sdmodfirstderiv,"z")

    delta<-as.matrix(delta)
    beta<-as.matrix(beta)

    vectorsprobs<-lapply(c(1:no.outcomes),calcprobs)
    # if weights are used then the standard log likelihood is no longer a relevant measure of model suitability.
    # need to calculate a pseudo-log likelihood, also for the baseline log-likelihood should take account of weights

    #loglikelihood<-sum(sapply(suppressWarnings(lapply(vectorsprobs,log)),sum))
    wloglikelihoodvecs<-list()
    for (i in 1:no.outcomes){
      wloglikelihoodvecs[[i]]<-suppressWarnings(log(vectorsprobs[[i]]))*weightsr[[i]]
    }
    loglikelihood<-sum(sapply(wloglikelihoodvecs,sum))

    # write function to work with the matrices for each outcome separately
    # produce the relevant gradient and hessian.
    # afterwards can sum. Allows the use of the Map function.

    getLLgradhess<-function(X,Z,w,index){
      j<-index
      Xb<-X%*%beta
      Zd<-Z%*%delta
      Zdinv<-1/eval({z<-Zd;sdmodel})
      if (j==1){
        frac0<-rep(-Inf,length(Zd))
        frac1<-(threshparam[j+1]-Xb)*Zdinv
      } else if (j==no.outcomes){
        frac1<-rep(Inf,length(Zd))
        frac0<-(threshparam[j]-Xb)*Zdinv
      } else {
        frac1<-(threshparam[j+1]-Xb)*Zdinv
        frac0<-(threshparam[j]-Xb)*Zdinv
      }

      # functions used to calculate score and hessian.
      calcscorebeta<-function(beta){
        # only to be called inside LLoglmx
        if (j==1){
          ProbderivBeta<-X[,beta]*Zdinv*(-ProbFuncD(frac1))/vectorsprobs[[j]]
        } else if (j==no.outcomes){
          ProbderivBeta<-X[,beta]*Zdinv*(ProbFuncD(frac0))/vectorsprobs[[j]]
        } else {
          ProbderivBeta<-X[,beta]*Zdinv*(ProbFuncD(frac0)-ProbFuncD(frac1))/vectorsprobs[[j]]
        }
        ProbderivBeta
      }

      calcscoredelta<-function(delta){
        # only to be called inside LLoglmx
        if (j==1){
          ProbderivDelta<- -Z[,delta]*eval({z<-Zd;sdmodfirstderiv})*Zdinv*frac1*ProbFuncD(frac1)/vectorsprobs[[j]]
        } else if (j==no.outcomes){
          ProbderivDelta<- Z[,delta]*eval({z<-Zd;sdmodfirstderiv})*Zdinv*frac0*ProbFuncD(frac0)/vectorsprobs[[j]]
        } else {
          ProbderivDelta<- -Z[,delta]*eval({z<-Zd;sdmodfirstderiv})*Zdinv*(frac1*ProbFuncD(frac1)-frac0*ProbFuncD(frac0))/vectorsprobs[[j]]
        }
        ProbderivDelta
      }

      calcscorethreshparam<-function(alpha){
        # only to be called inside LLoglmx
        if (alpha==j){
          Probderivalpha<- -Zdinv*ProbFuncD(frac0)/vectorsprobs[[j]]
        } else if (alpha==j+1){
          Probderivalpha<- Zdinv*ProbFuncD(frac1)/vectorsprobs[[j]]
        } else {
          Probderivalpha<-vector("numeric",nrow(X))
        }
        Probderivalpha
      }

      scorevector<-vector("numeric",no.parameters)

      if (no.betaparams>0){
        vectorprobderivbeta<-sapply(c(1:length(beta))[collectbeta],calcscorebeta)
        scorevector[1:sum(collectbeta)]<-scorevector[1:sum(collectbeta)]+apply(vectorprobderivbeta*w,2,sum)
      }

      if (no.deltaparams>0){
        vectorprobderivdelta<-sapply(c(1:length(delta))[collectdelta],calcscoredelta)
        scorevector[(1+no.betaparams):(no.betaparams+no.deltaparams)]<-scorevector[(1+no.betaparams):(no.betaparams+no.deltaparams)]+apply(vectorprobderivdelta*w,2,sum)
      }

      if (no.threshparams>0){
        vectorprobderivthreshparam<-sapply(c(2:(length(threshparam)-1))[collectthreshparam],calcscorethreshparam)
        scorevector[(1+no.betaparams+no.deltaparams):no.parameters]<-scorevector[(1+no.betaparams+no.deltaparams):no.parameters]+apply(vectorprobderivthreshparam*w,2,sum)
      }

      if (robustmatrix){
        if (no.betaparams>0 & no.deltaparams>0){
          scorevecs<-cbind(vectorprobderivbeta,vectorprobderivdelta)
          if (no.threshparams>0){
            scorevecs<-cbind(scorevecs,vectorprobderivthreshparam)
          }
        } else if (no.betaparams>0){
          scorevecs<-vectorprobderivbeta
          if (no.threshparams>0){
            scorevecs<-cbind(scorevecs,vectorprobderivthreshparam)
          }
        } else if (no.deltaparams>0){
          scorevecs<-vectorprobderivdelta
          if (no.threshparams>0){
            scorevecs<-cbind(scorevecs,vectorprobderivthreshparam)
          }
        }
        BHHHmatrix<-matrix(vector("numeric",ncol(scorevecs)^2),nrow=ncol(scorevecs),ncol=ncol(scorevecs))
        for (v in 1:nrow(scorevecs)){
          obmat<-scorevecs[v,]%*%t(scorevecs[v,])
          BHHHmatrix<-BHHHmatrix+(w[v]^2)*obmat
        }
      } else {
        BHHHmatrix<-NULL
      }


      if (analhessian){
        calc2ndderivprobbetabeta<-function(x){ # x is a two element vector specifying location of each coefficient
          # only to be called inside LLoglmx
          if (j==1){
            probderiv2beta2<-sum((X[,x[1]]*X[,x[2]]*(Zdinv^2)*(ProbFuncDD(frac1))/vectorsprobs[[j]])*w)
          } else if (j==no.outcomes){
            probderiv2beta2<-sum((X[,x[1]]*X[,x[2]]*(Zdinv^2)*(-ProbFuncDD(frac0))/vectorsprobs[[j]])*w)
          } else {
            probderiv2beta2<-sum((X[,x[1]]*X[,x[2]]*(Zdinv^2)*(ProbFuncDD(frac1)-ProbFuncDD(frac0))/vectorsprobs[[j]])*w)
          }
          probderiv2beta2
        }

        calc2ndderivprobbetadelta<-function(x){ # x is a two element vector specifying location of each coefficient, beta 1st, delta 2nd
          # only to be called inside LLoglmx
          if (j==1){
            probderiv2betadelta<- X[,x[1]]*Z[,x[2]]*(Zdinv^2)*eval({z<-Zd;sdmodfirstderiv})*(ProbFuncD(frac1)+frac1*ProbFuncDD(frac1))/vectorsprobs[[j]]
          } else if (j==no.outcomes){
            probderiv2betadelta<- X[,x[1]]*Z[,x[2]]*(Zdinv^2)*eval({z<-Zd;sdmodfirstderiv})*(-ProbFuncD(frac0)-frac0*ProbFuncDD(frac0))/vectorsprobs[[j]]
          } else {
            probderiv2betadelta<- X[,x[1]]*Z[,x[2]]*(Zdinv^2)*eval({z<-Zd;sdmodfirstderiv})*(ProbFuncD(frac1)-ProbFuncD(frac0)+frac1*ProbFuncDD(frac1)-frac0*ProbFuncDD(frac0))/vectorsprobs[[j]]
          }
          sum(probderiv2betadelta*w)
        }

        calc2ndderivprobbetaalpha<-function(x){ # x is a two element vector specifying location of each coefficient, beta 1st, alpha 2nd
          # only to be called inside LLoglmx
          if ((x[2]==j) & (j>1)){
            probderiv2betaalpha<- X[,x[1]]*(Zdinv^2)*ProbFuncDD(frac0)/vectorsprobs[[j]]
          } else if ((x[2]==j+1) & (j<no.outcomes)){
            probderiv2betaalpha<- -X[,x[1]]*(Zdinv^2)*ProbFuncDD(frac1)/vectorsprobs[[j]]
          } else {
            probderiv2betaalpha<-vector("numeric",nrow(Xr[[j]]))
          }
          sum(probderiv2betaalpha*w)
        }

        calc2ndderivprobdeltadelta<-function(x){ # x is a two element vector specifying location of each coefficient
          # only to be called inside LLoglmx
          if (j==1){
            probderiv2deltadelta<- Z[,x[1]]*Z[,x[2]]*Zdinv*((2*(eval({z<-Zd;sdmodfirstderiv})^2)*Zdinv-eval({z<-Zd;sdmodsecondderiv}))*(frac1*ProbFuncD(frac1))+(eval({z<-Zd;sdmodfirstderiv})^2)*Zdinv*((frac1^2)*ProbFuncDD(frac1)))/vectorsprobs[[j]]
          } else if (j==no.outcomes){
            probderiv2deltadelta<- Z[,x[1]]*Z[,x[2]]*Zdinv*((2*(eval({z<-Zd;sdmodfirstderiv})^2)*Zdinv-eval({z<-Zd;sdmodsecondderiv}))*(-frac0*ProbFuncD(frac0))+(eval({z<-Zd;sdmodfirstderiv})^2)*Zdinv*(-(frac0^2)*ProbFuncDD(frac0)))/vectorsprobs[[j]]
          } else {
            probderiv2deltadelta<- Z[,x[1]]*Z[,x[2]]*Zdinv*((2*(eval({z<-Zd;sdmodfirstderiv})^2)*Zdinv-eval({z<-Zd;sdmodsecondderiv}))*(frac1*ProbFuncD(frac1)-frac0*ProbFuncD(frac0))+(eval({z<-Zd;sdmodfirstderiv})^2)*Zdinv*((frac1^2)*ProbFuncDD(frac1)-(frac0^2)*ProbFuncDD(frac0)))/vectorsprobs[[j]]
          }
          sum(probderiv2deltadelta*w)
        }

        calc2ndderivprobdeltaalpha<-function(x){
          # only to be called inside LLoglmx
          if (x[2]==j & j>1){
            probderiv2deltaalpha<- Z[,x[1]]*(Zdinv^2)*eval({z<-Zd;sdmodfirstderiv})*(ProbFuncD(frac0)+frac0*ProbFuncDD(frac0))/vectorsprobs[[j]]
          } else if (x[2]==j+1 & j<no.outcomes){
            probderiv2deltaalpha<-  -Z[,x[1]]*(Zdinv^2)*eval({z<-Zd;sdmodfirstderiv})*(ProbFuncD(frac1)+frac1*ProbFuncDD(frac1))/vectorsprobs[[j]]
          } else {
            probderiv2deltaalpha<-vector("numeric",nrow(Z))
          }
          sum(probderiv2deltaalpha*w)
        }

        calc2ndderivprobalphaalpha<-function(x){
          # only to be called inside LLoglmx
          if (x[1]==x[2] & x[1]==j & j>1){
            probderiv2alphaalpha<- -(Zdinv^2)*ProbFuncDD(frac0)/vectorsprobs[[j]]
          } else if (x[1]==x[2] & x[1]==j+1 & j<no.outcomes){
            probderiv2alphaalpha<- (Zdinv^2)*ProbFuncDD(frac1)/vectorsprobs[[j]]
          } else {
            probderiv2alphaalpha<-vector("numeric",nrow(X))
          }
          sum(probderiv2alphaalpha*w)
        }

        hessian<-matrix(vector("numeric",no.parameters^2),nrow=no.parameters,ncol=no.parameters)
        # first add the term that is derived from the cross product of gradients
        if (no.betaparams>0 & no.deltaparams>0){
          scorevecs<-cbind(vectorprobderivbeta,vectorprobderivdelta)
          if (no.threshparams>0){
            scorevecs<-cbind(scorevecs,vectorprobderivthreshparam)
          }
          crossprodterms<- -t(scorevecs)%*%(scorevecs*w)
          crossprodterms[upper.tri(crossprodterms)]<-0
          hessian<-hessian+crossprodterms
        } else if (no.betaparams>0){
          scorevecs<-vectorprobderivbeta
          if (no.threshparams>0){
            scorevecs<-cbind(scorevecs,vectorprobderivthreshparam)
          }
          crossprodterms<- -t(scorevecs)%*%(scorevecs*w)
          crossprodterms[upper.tri(crossprodterms)]<-0
          hessian<-hessian+crossprodterms
        } else if (no.deltaparams>0){
          scorevecs<-vectorprobderivdelta
          if (no.threshparams>0){
            scorevecs<-cbind(scorevecs,vectorprobderivthreshparam)
          }
          crossprodterms<- -t(scorevecs)%*%(scorevecs*w)
          crossprodterms[upper.tri(crossprodterms)]<-0
          hessian<-hessian+crossprodterms
        }
        # then do the same with the cross derivatives
        # clear space in memory, remove first derivatives
        if (no.deltaparams>0){rm(vectorprobderivdelta)}
        if (no.betaparams>0){rm(vectorprobderivbeta)}
        if (no.threshparams>0){rm(vectorprobderivthreshparam)}

        if (no.betaparams>0){
          listderivs<-.paircombn(c(1:length(beta))[collectbeta])
          secondderivterms<-sapply(listderivs,calc2ndderivprobbetabeta)
          startrow<-1
          endrow<-no.betaparams
          startcol_1<-0
          counter<-1
          while (counter<=no.betaparams){
            hessian[(startcol_1+counter):endrow,counter+startcol_1]<-hessian[(startcol_1+counter):endrow,counter+startcol_1]+secondderivterms[(no.betaparams*(counter-1)-(counter-2)*(counter-1)/2+1):(no.betaparams*counter-(counter*(counter-1)/2))]
            counter<-counter+1
          }
        }

        if (no.deltaparams>0){
          listderivs<-.paircombn(c(1:length(delta))[collectdelta])
          secondderivterms<-sapply(listderivs,calc2ndderivprobdeltadelta)
          startrow<-no.betaparams+1
          endrow<-no.betaparams+no.deltaparams
          startcol_1<-no.betaparams
          counter<-1
          while (counter<=no.deltaparams){
            hessian[(startcol_1+counter):endrow,counter+startcol_1]<-hessian[(startcol_1+counter):endrow,counter+startcol_1]+secondderivterms[(no.deltaparams*(counter-1)-(counter-2)*(counter-1)/2+1):(no.deltaparams*counter-(counter*(counter-1)/2))]
            counter<-counter+1
          }
        }

        if (no.threshparams>0){
          listderivs<-.paircombn(c(2:(length(threshparam)-1))[collectthreshparam])
          secondderivterms<-sapply(listderivs,calc2ndderivprobalphaalpha)
          startrow<-no.betaparams+no.deltaparams+1
          endrow<-no.betaparams+no.deltaparams+no.threshparams
          startcol_1<-no.betaparams+no.deltaparams
          counter<-1
          while (counter<=no.threshparams){
            hessian[(startcol_1+counter):endrow,counter+startcol_1]<-hessian[(startcol_1+counter):endrow,counter+startcol_1]+secondderivterms[((no.threshparams*(counter-1)-(counter-2)*(counter-1)/2+1)):(no.threshparams*counter-(counter*(counter-1)/2))]
            counter<-counter+1
          }
        }

        if (no.betaparams>0 & no.deltaparams>0){
          listderivs<-.paircombn(c(1:length(beta))[collectbeta],c(1:length(delta))[collectdelta],same=FALSE)
          secondderivterms<-sapply(listderivs,calc2ndderivprobbetadelta)
          startrow<-no.betaparams+1
          startcol_1<-0
          endrow<-no.betaparams+no.deltaparams
          counter<-1
          while (counter<=no.betaparams){
            hessian[startrow:endrow,counter]<-hessian[startrow:endrow,counter]+secondderivterms[(1+(counter-1)*no.deltaparams):(counter*no.deltaparams)]
            counter<-counter+1
          }
        }

        if (no.betaparams>0 & no.threshparams>0){
          listderivs<-.paircombn(c(1:length(beta))[collectbeta],c(2:(length(threshparam)-1))[collectthreshparam],same=FALSE)
          secondderivterms<-sapply(listderivs,calc2ndderivprobbetaalpha)
          startrow<-no.betaparams+no.deltaparams+1
          startcol_1<-0
          endrow<-no.betaparams+no.deltaparams+no.threshparams
          counter<-1
          while (counter<=no.betaparams){
            hessian[startrow:endrow,counter]<-hessian[startrow:endrow,counter]+secondderivterms[(1+(counter-1)*no.threshparams):(counter*no.threshparams)]
            counter<-counter+1
          }
        }

        if (no.deltaparams>0 & no.threshparams>0){
          listderivs<-.paircombn(c(1:length(delta))[collectdelta],c(2:(length(threshparam)-1))[collectthreshparam],same=FALSE)
          secondderivterms<-sapply(listderivs,calc2ndderivprobdeltaalpha)
          startrow<-no.betaparams+no.deltaparams+1
          startcol_1<-no.betaparams
          endrow<-no.betaparams+no.deltaparams+no.threshparams
          counter<-1
          while (counter<=no.deltaparams){
            hessian[startrow:endrow,counter+startcol_1]<-hessian[startrow:endrow,counter+startcol_1]+secondderivterms[(1+(counter-1)*no.threshparams):(counter*no.threshparams)]
            counter<-counter+1
          }
        }
      }
      output<-list(scorevector,hessian,BHHHmatrix)
    }
    if (SameModelMEANSD){
      collectresults<-Map(getLLgradhess,Xr,Xr,weightsr,c(1:no.outcomes))
    } else {
      collectresults<-Map(getLLgradhess,Xr,Zr,weightsr,c(1:no.outcomes))
    }
    scorevector<-Reduce("+",lapply(collectresults,function(x){x[[1]]}))

    attr(loglikelihood,"gradient")<-scorevector
    if (analhessian){
      # if coded correctly the hessian calculated up to now is lower triangular, need to make it symmetric
      hessian<-Reduce("+",lapply(collectresults,function(x){x[[2]]}))
      hessian<-hessian+t(hessian)-diag(diag(hessian))
      attr(loglikelihood,"hessian")<-hessian
    }
    if (robustmatrix){
      BHHHmatrix<-Reduce("+",lapply(collectresults, function(x){x[[3]]}))
      attr(loglikelihood,"BHHHhessian")<-BHHHmatrix
    }
    loglikelihood
  }

  # function that calculates the log-likelihood, as a function of the parameters that are not prespecified.
  # the maxLik function used to maximise the log-likelihood requires a function of only the estimated parameters.
  LLoglmxTOP<-function(param){LLoglmx(param,beta=beta,delta=delta,threshparam=threshparam,analhessian=analhessian)}
  # call maxLik to estimate parameters
  maxLikRes<-maxLik(LLoglmxTOP,start=start,iterlim=300,finalHessian=TRUE,method="NR")
  # just remains to extract results and store.
  #stop("getshere")
  coefficients<-maxLikRes$estimate
  if (robustmatrix){
    LLoutput<-LLoglmx(coefficients,beta=beta,delta=delta,threshparam=threshparam,analhessian=analhessian,robustmatrix = TRUE)
    BHHHmatrix<-attr(LLoutput,"BHHHhessian")
  } else {
    BHHHmatrix<-NULL
  }

  # separate coefficients into beta, delta and alpha, with prespecified values
  betacoeffs<-deltacoeffs<-threshparamcoeffs<-vector("logical",length(coefficients))
  if (no.betaparams>0){
    beta[is.na(beta)]<-coefficients[1:no.betaparams]
    betacoeffs[1:no.betaparams]<-TRUE
  }
  if (no.deltaparams>0){
    delta[is.na(delta)]<-coefficients[(no.betaparams+1):(no.betaparams+no.deltaparams)]
    deltacoeffs[(no.betaparams+1):(no.betaparams+no.deltaparams)]<-TRUE
  }
  if (no.threshparams>0){
    threshparam[is.na(threshparam)]<-coefficients[(no.betaparams+no.deltaparams+1):(no.betaparams+no.deltaparams+no.threshparams)]
    threshparamcoeffs[(no.betaparams+no.deltaparams+1):(no.betaparams+no.deltaparams+no.threshparams)]<-TRUE
  }
  allparameters<-list(beta,delta,threshparam)

  names(allparameters)<-c("beta","delta","threshparam")
  coeff_type<-list(betacoeffs,deltacoeffs,threshparamcoeffs)
  names(coefficients)<-outputnames
  attr(coefficients,"coefftypes")<-coeff_type
  loglikelihood<-maxLikRes$maximum
  #attr(loglikelihood,"BaselineLL")<-BaselineLL
  attr(loglikelihood,"No.Obs")<-No.Obs
  if (sum(w==rep(1,length(w)))!=length(w)){
    weights<-weights
  } else {
    weights<-NULL
  }

  results<-list(loglikelihood=loglikelihood,link=link,no.iterations=maxLikRes$iterations,coefficients=coefficients,returnCode=maxLikRes$code,gradient=maxLikRes$gradient
                ,hessian=maxLikRes$hessian,BHHHhessian=BHHHmatrix,NOutcomes=no.outcomes,Outcomes=listoutcomes,sdmodel=sdmodel,allparams=allparameters,Est.Parameters=Est.Parameters,modelframes=modelframes)
  class(results)<-c("oglmx.fit")
  invisible(results)
}
