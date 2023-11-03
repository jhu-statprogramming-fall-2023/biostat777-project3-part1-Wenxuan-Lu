mu_func<-function(x,family){
    if(family=='gaussian'){return(x)}
    else if(family=='binomial'){return(expit(x))}
}
mu_prime_func<-function(x,family){
    if(family=='gaussian'){return(1)}
    else if(family=='binomial'){return(expit(x)*(1-expit(x)))}
}
Delta_opt<-function(y,Z,W,family,
                    study_info,A=NULL,hat_thetaA=NULL,
                    beta=NULL){
    if(is.null(A)){
        pA=0
    }else if(is.null(dim(A)[1])){
        pA=1
    }else{
        pA=ncol(A)
    }
    X=cbind(A,Z,W)
    XR=cbind(A,Z)
    n_main=nrow(Z)
    tilde_thetaZ=study_info[[1]]$Coeff
    tilde_theta=c(hat_thetaA,tilde_thetaZ)
    mu_X_beta=mu_func(X%*%beta,family)
    mu_XR_theta=mu_func(XR%*%tilde_theta,family)
    V_U1=(1/n_main)*crossprod(X*c(mu_X_beta-y))
    V_U2=(1/n_main)*crossprod(Z*c(mu_X_beta-mu_XR_theta))
    Cov_U1U2=(1/n_main)*crossprod(X*c(mu_X_beta-y),Z*c(mu_X_beta-mu_XR_theta))
    Gamma2Z=(1/n_main)*crossprod(Z*mu_prime_func(mu_XR_theta,family),Z)
    V_thetaZ=study_info[[1]]$Covariance
    Delta22 = V_U2
    +Gamma2Z%*%V_thetaZ%*%t(Gamma2Z)
    Delta12 = Cov_U1U2
    if(pA!=0){
        Gamma2A=(1/n_main)*crossprod(Z*mu_prime_func(mu_XR_theta,family),A)
        inv_Gamma3=ginv((1/n_main)*crossprod(XR*mu_prime_func(mu_XR_theta,family),XR))
        V_thetaA=inv_Gamma3[1:pA,]%*%(1/n_main)*crossprod(XR*c(mu_XR_theta-y))%*%inv_Gamma3[,1:pA]

        Cov_U1theta=(1/n_main)*crossprod(X*c(mu_X_beta-y),XR*c(mu_XR_theta-y))%*%inv_Gamma3[,1:pA]%*%t(Gamma2A)
        Cov_U2theta=(1/n_main)*crossprod(Z*c(mu_X_beta-mu_XR_theta),XR*c(mu_XR_theta-y))%*%inv_Gamma3[,1:pA]%*%t(Gamma2A)

        Delta22 = Delta22 + Gamma2A%*%V_thetaA%*%t(Gamma2A)
        + Cov_U2theta+t(Cov_U2theta)
        Delta12 = Delta12 + Cov_U1theta
    }

    Delta = rbind(cbind(V_U1,Delta12),cbind(t(Delta12),Delta22))
    Delta
}

pseudo_Xy_gaussian<-function(
        C_half,Z,W,A,y,beta=NULL,hat_thetaA=NULL,study_info=NULL){
    tilde_thetaZ=study_info[[1]]$Coeff
    tilde_theta=c(hat_thetaA,tilde_thetaZ)
    X=cbind(A,Z,W)
    XR=cbind(A,Z)
    pseudo_X=C_half%*%crossprod(cbind(X,Z),X)
    pseudo_y1=c(crossprod(y,X))
    pseudo_y2=c(crossprod(Z,XR)%*%tilde_theta)
    pseudo_y=c(c(pseudo_y1,pseudo_y2)%*%C_half)
    list("pseudo_X"=pseudo_X,"pseudo_y"=pseudo_y)
}

pseudo_Xy_binomial<-function(
        C_half,Z,W,A,y,beta=NULL,hat_thetaA=NULL,study_info=NULL){
    tilde_thetaZ=study_info[[1]]$Coeff
    tilde_theta=c(hat_thetaA,tilde_thetaZ)
    X=cbind(A,Z,W)
    XR=cbind(A,Z)
    expit_beta=c(expit(X%*%beta))
    dexpit_beta=c(expit_beta*(1-expit_beta))
    pseudo_X=C_half%*%crossprod(cbind(X,Z),X*dexpit_beta)
    u1=crossprod((expit_beta-y),X)
    u2=c(expit_beta-c(expit(XR%*%tilde_theta)))%*%Z
    pseudo_y= -C_half%*%c(u1,u2) + c(pseudo_X%*%beta)
    list("pseudo_X"=pseudo_X,"pseudo_y"=pseudo_y)
}







## cross validation function for continuous y with lambda only
cv_mse_lambda_func<-function(index_fold,Z,W,A,y,
                             C_half,beta_initial,hat_thetaA,
                             study_info,lambda_list,
                             w_adaptive,final_alpha,
                             pseudo_Xy){
    fold_mse_lambda<-sapply(1:length(index_fold), function(cur_fold){
        index_test<-index_fold[[cur_fold]]
        Ztrain<-Z[-index_test,]
        Ztest<-Z[index_test,]
        if(!is.null(W)){
            Wtrain<-W[-index_test,]
            Wtest<-W[index_test,]
        }else{
            Wtrain<-NULL
            Wtest<-NULL}
        if(!is.null(A)){
            Atrain<-A[-index_test,]
            Atest<-A[index_test,]
        }else{
            Atrain<-NULL
            Atest<-NULL}
        ytrain<-y[-index_test]
        ytest<-y[index_test]
        pseudo_Xy_list_train<-pseudo_Xy(C_half,Ztrain,Wtrain,Atrain,
                                        ytrain,beta = beta_initial,hat_thetaA = hat_thetaA,
                                        study_info=study_info)
        initial_sf_train<-nrow(Ztrain)/sqrt(nrow(pseudo_Xy_list_train$pseudo_X))
        pseudo_X_train<-pseudo_Xy_list_train$pseudo_X/initial_sf_train
        pseudo_y_train<-pseudo_Xy_list_train$pseudo_y/initial_sf_train
        mse_lam<-sapply(lambda_list,function(cur_lam){
            cv_fit<-glmnet(x= (pseudo_X_train),y= (pseudo_y_train),
                           standardize=F,intercept=F,
                           alpha = final_alpha,penalty.factor = w_adaptive,lambda = cur_lam)
            cur_beta<-coef.glmnet(cv_fit)[-1]
            mean((cbind(Atest,Ztest,Wtest)%*%cur_beta - ytest)^2)
        })
        mse_lam
    })
    cv_mse_lambda<-rowSums(fold_mse_lambda)
    cv_mse_lambda
}

## cross validation function for continuous y with lambda and ratio
cv_mse_lambda_ratio_func<-function(index_fold,Z,W,A,y,
                                   C_half,beta_initial,hat_thetaA,
                                   study_info,lambda_list,
                                   ratio_range,pZ,pW,pA,
                                   w_adaptive,final_alpha,
                                   pseudo_Xy){
    fold_mse_lambda_ratio<-lapply(1:length(index_fold), function(cur_fold){
        index_test<-index_fold[[cur_fold]]
        Ztrain<-Z[-index_test,]
        Ztest<-Z[index_test,]
        if(!is.null(W)){
            Wtrain<-W[-index_test,]
            Wtest<-W[index_test,]
        }else{
            Wtrain<-NULL
            Wtest<-NULL}
        if(!is.null(A)){
            Atrain<-A[-index_test,]
            Atest<-A[index_test,]
        }else{
            Atrain<-NULL
            Atest<-NULL}
        ytrain<-y[-index_test]
        ytest<-y[index_test]
        pseudo_Xy_list_train<-pseudo_Xy(C_half,Ztrain,Wtrain,Atrain,
                                        ytrain,beta = beta_initial,hat_thetaA = hat_thetaA,
                                        study_info=study_info)
        initial_sf_train<-nrow(Ztrain)/sqrt(nrow(pseudo_Xy_list_train$pseudo_X))
        pseudo_X_train<-pseudo_Xy_list_train$pseudo_X/initial_sf_train
        pseudo_y_train<-pseudo_Xy_list_train$pseudo_y/initial_sf_train
        mse_lam_ratio_fold<-sapply(lambda_list,function(cur_lam){
            sapply(ratio_range,function(cur_ratio){
                ratio_vec<-c(rep(1,pA),rep(cur_ratio,pZ),rep(1,pW))
                w_adaptive_ratio<-w_adaptive*ratio_vec
                cv_fit<-glmnet(x=pseudo_X_train,y=pseudo_y_train,
                               standardize=F,intercept=F,alpha = final_alpha,
                               penalty.factor = w_adaptive_ratio,
                               lambda = cur_lam)
                cur_beta<-coef.glmnet(cv_fit)[-1]
                mean((cbind(Atest,Ztest,Wtest)%*%cur_beta - ytest)^2)
            })
        }) # row is ratio_range & col is lambda_list
        mse_lam_ratio_fold
    })
    cv_mse_lambda_ratio<-Reduce(`+`, fold_mse_lambda_ratio)
    cv_mse_lambda_ratio
}


## cross validation function for binary y with lambda only
cv_dev_lambda_func<-function(index_fold,Z,W,A,y,
                             C_half,beta_initial,hat_thetaA,
                             study_info,lambda_list,
                             w_adaptive,final_alpha,
                             pseudo_Xy){
    dev_fold<-lapply(1:length(index_fold), function(cur_fold){
        index_test<-index_fold[[cur_fold]]
        Ztrain<-Z[-index_test,]
        Ztest<-Z[index_test,]
        if(!is.null(W)){
            Wtrain<-W[-index_test,]
            Wtest<-W[index_test,]
        }else{
            Wtrain<-NULL
            Wtest<-NULL}
        if(!is.null(A)){
            Atrain<-A[-index_test,]
            Atest<-A[index_test,]
        }else{
            Atrain<-NULL
            Atest<-NULL}
        ytrain<-y[-index_test]
        ytest<-y[index_test]
        pseudo_Xy_list_train<-pseudo_Xy(C_half,Ztrain,Wtrain,Atrain,
                                        ytrain,beta = beta_initial,hat_thetaA = hat_thetaA,
                                        study_info=study_info)
        initial_sf_train<-nrow(Ztrain)/sqrt(nrow(pseudo_Xy_list_train$pseudo_X))
        pseudo_X_train<-pseudo_Xy_list_train$pseudo_X/initial_sf_train
        pseudo_y_train<-pseudo_Xy_list_train$pseudo_y/initial_sf_train
        dev_lam<-sapply(lambda_list,function(cur_lam){
            cv_fit<-glmnet(x= (pseudo_X_train),y= (pseudo_y_train),
                           standardize=F,intercept=F,
                           alpha = final_alpha,penalty.factor = w_adaptive,lambda = cur_lam)
            cur_beta<-coef.glmnet(cv_fit)[-1]
            probtest <- c(expit(cbind(Atest,Ztest,Wtest)%*%cur_beta))
            cur_dev <- -2*sum( ytest * log(probtest) + (1 - ytest) * log(1 - probtest) )
            suppressMessages(cur_auc<-c(auc(ytest,c(expit(cbind(Atest,Ztest,Wtest)%*%cur_beta)),direction = "<")))
            #sum((cbind(Ztest,Wtest,Atest)%*%cur_beta - ytest)^2)
            c(cur_dev,cur_auc)
        })
        dev_lam
    })
    sum_dev_lam<-Reduce(`+`, dev_fold)
    list("deviance"=sum_dev_lam[1,],"auc"=sum_dev_lam[2,])
}

## cross validation function for continuous y with lambda and ratio
cv_dev_lambda_ratio_func<-function(index_fold,Z,W,A,y,
                                   C_half,beta_initial,hat_thetaA,
                                   study_info,lambda_list,
                                   ratio_range,pZ,pW,pA,
                                   w_adaptive,final_alpha,
                                   pseudo_Xy){
    dev_lam_ratio<-lapply(1:length(index_fold), function(cur_fold){
        index_test<-index_fold[[cur_fold]]
        Ztrain<-Z[-index_test,]
        Ztest<-Z[index_test,]
        if(!is.null(W)){
            Wtrain<-W[-index_test,]
            Wtest<-W[index_test,]
        }else{
            Wtrain<-NULL
            Wtest<-NULL}
        if(!is.null(A)){
            Atrain<-A[-index_test,]
            Atest<-A[index_test,]
        }else{
            Atrain<-NULL
            Atest<-NULL}
        ytrain<-y[-index_test]
        ytest<-y[index_test]
        pseudo_Xy_list_train<-pseudo_Xy(C_half,Ztrain,Wtrain,Atrain,
                                        ytrain,beta = beta_initial,hat_thetaA = hat_thetaA,
                                        study_info=study_info)
        initial_sf_train<-nrow(Ztrain)/sqrt(nrow(pseudo_Xy_list_train$pseudo_X))
        pseudo_X_train<-pseudo_Xy_list_train$pseudo_X/initial_sf_train
        pseudo_y_train<-pseudo_Xy_list_train$pseudo_y/initial_sf_train
        dev_lam_ratio_fold<-lapply(ratio_range,function(cur_ratio){
            sapply(lambda_list,function(cur_lam){
                ratio_vec<-c(rep(1,pA),rep(cur_ratio,pZ),rep(1,pW))
                w_adaptive_ratio<-w_adaptive*ratio_vec
                cv_fit<-glmnet(x=pseudo_X_train,y=pseudo_y_train,
                               standardize=F,intercept=F,alpha = final_alpha,
                               penalty.factor = w_adaptive_ratio,
                               lambda = cur_lam)
                cur_beta<-coef.glmnet(cv_fit)[-1]
                probtest <- c(expit(cbind(Atest,Ztest,Wtest)%*%cur_beta))
                cur_dev <- -2*sum( ytest * log(probtest) + (1 - ytest) * log(1 - probtest) )
                suppressMessages(cur_auc<-c(auc(ytest,c(expit(cbind(Atest,Ztest,Wtest)%*%cur_beta)),direction = "<")))
                c(cur_dev,cur_auc)
            })
        }) # row is ratio_range & col is lambda_list

        list("deviance"=do.call(rbind, lapply(dev_lam_ratio_fold, function(m) m[1,])),
             "auc"=do.call(rbind, lapply(dev_lam_ratio_fold, function(m) m[2,])))
    })
    dev_lam_ratio1<-lapply(1:length(index_fold), function(cur_fold){
        dev_lam_ratio[[cur_fold]]$deviance
    })
    dev_lam_ratio2<-lapply(1:length(index_fold), function(cur_fold){
        dev_lam_ratio[[cur_fold]]$auc
    })
    list("deviance"=Reduce(`+`, dev_lam_ratio1),
         "auc"=Reduce(`+`, dev_lam_ratio2))
}




htlgmm.default<-function(
        y,Z,W=NULL,
        study_info=NULL,
        A=1,
        penalty_type = "lasso",
        family = "gaussian",
        initial_with_type = "ridge",
        beta_initial = NULL,
        hat_thetaA = NULL,
        remove_penalty_Z = FALSE,
        remove_penalty_W = FALSE,
        inference = TRUE,
        use_cv = TRUE,
        type_measure = "default",
        nfolds = 10,
        fix_lambda = NULL,
        lambda_list = NULL,
        tune_ratio = FALSE,
        fix_ratio = NULL,
        ratio_list = NULL,
        gamma_adaptivelasso = 1/2,
        use_sparseC = TRUE,
        seed.use = 97
){
    set.seed(seed.use)
    if (is.null(study_info)){stop("Please input study_info as trained model")}
    if(!penalty_type %in% c("ols","adaptivelasso","lasso","ridge")){
        stop("Select penalty type from c('ols','adaptivelasso','lasso','ridge').")
    }
    if(!type_measure%in% c("default", "mse", "deviance", "auc")){
        stop("Select type_measure from c('default','deviance','auc')")
    }
    final_alpha = 1
    if(penalty_type == "ridge"){final_alpha = 0}
    if(penalty_type == "ols"){
        initial_with_type = "ols"
        use_cv = FALSE
        tune_ratio = FALSE
    }

    if(!is.null(fix_ratio)){
        if(tune_ratio){
            stop("If ratio is fixed, please set tune_ratio as FALSE")
        }else if(remove_penalty_Z | remove_penalty_W){
            stop("If ratio is fixed, please set remove_penalty's as FALSE")
        }
    }

    nZ=nrow(Z)
    nZext=study_info[[1]]$Sample_size
    if(A==1){A=matrix(1,nrow=nZ,ncol=1)}
    pZ=ncol(Z)
    if(is.null(W)){pW=0}else{pW=ncol(W)}
    if(is.null(A)){pA=0}else{pA=ncol(A)}
    if(nZ<2*pZ+pW+pA){use_sparseC=TRUE}

    if(family == "gaussian"){pseudo_Xy=pseudo_Xy_gaussian}
    else if(family == "binomial"){pseudo_Xy=pseudo_Xy_binomial}

    if(pA!=0 & is.null(hat_thetaA)){
        hat_thetaA_glm=glm(y~0+.,data = data.frame(y,A,Z),family = family)
        hat_thetaA=hat_thetaA_glm$coefficients[1:pA]
    }

    X<-cbind(A,Z,W)
    Aid<-1:pA
    Zid<-(pA+1):(pA+pZ)
    Wid<-(pA+pZ+1):(pA+pZ+pW)
    fix_penalty<-c(rep(0,pA),rep(1,pZ+pW))
    if(remove_penalty_Z){fix_penalty[Zid]<-0}else{
        if(!is.null(fix_ratio)){fix_penalty[Zid]<-fix_ratio}
    }
    if(remove_penalty_W){fix_penalty[Wid]<-0}
    if(remove_penalty_Z & remove_penalty_W){
        penalty_type = "ols"
        warning("All penalties are removed, turn to ols without penalties!")
    }

    if(!is.null(beta_initial) & length(beta_initial)!=pA+pZ+pW){
        warning("beta_initial should be from A,Z,W.\n Length not match, compute default initial instead.")
        beta_initial=NULL
    }
    if(is.null(beta_initial) & initial_with_type %in% c("ols","ridge","lasso")){
        if(initial_with_type == "ridge"){initial_alpha=0}else{initial_alpha=1}
        if(initial_with_type == "ols"){
            fit_initial=glm(y~0+.,data = data.frame(y,X),family = family)
            beta_initial=fit_initial$coefficients
        }else if(pA == 0){
            fit_initial=cv.glmnet(x=X,y= y,
                                  alpha = initial_alpha,
                                  penalty.factor = fix_penalty,
                                  family=family)
            beta_initial=c(coef.glmnet(fit_initial,s="lambda.min")[-1])
        }else if(unique(A[,1])==1){
            fit_initial=cv.glmnet(x=X[,-1],y=y,
                                  alpha = initial_alpha,
                                  penalty.factor = fix_penalty[-1],
                                  family=family)
            beta_initial=c(coef.glmnet(fit_initial,s="lambda.min"))
        }else{
            stop("With A, the first column of A should be 1 for intercept.")
        }
    }else{stop("Select Initial Type from c('ridge','lasso')")}
    if (penalty_type == "adaptivelasso"){
        w_adaptive<-1/abs(beta_initial)^gamma_adaptivelasso
        w_adaptive[is.infinite(w_adaptive)]<-max(w_adaptive[!is.infinite(w_adaptive)])*100
        w_adaptive<-w_adaptive*fix_penalty
    }else{w_adaptive<-fix_penalty}
    # Estimation of C
    inv_C = Delta_opt(y=y,Z=Z,W=W,
                      family=family,
                      study_info=study_info,
                      A=A,hat_thetaA=hat_thetaA,
                      beta=beta_initial)
    if(use_sparseC){
        C_half<-diag(1/sqrt(diag(inv_C)))
    }else{
        inv_C_svd<-fast.svd(inv_C+diag(1e-15,nrow(inv_C)))
        C_half<-inv_C_svd$v%*%diag(1/sqrt(inv_C_svd$d))%*%t(inv_C_svd$u)
    }

    # Prepare for final model

    pseudo_Xy_list<-pseudo_Xy(C_half,Z,W,A,y,beta_initial,hat_thetaA,study_info)
    initial_sf<-nZ/sqrt(nrow(pseudo_Xy_list$pseudo_X))
    pseudo_X<-pseudo_Xy_list$pseudo_X/initial_sf
    pseudo_y<-pseudo_Xy_list$pseudo_y/initial_sf

    # Fit final model

    if(penalty_type != "ols"|is.null(fix_lambda)|is.null(lambda_list)){
        fit_final<-glmnet(x= pseudo_X,y= pseudo_y,standardize=F,
                          intercept=F,alpha = final_alpha,penalty.factor = w_adaptive)
        lambda_list<-fit_final$lambda
        lambda_list<-lambda_list[!is.na(lambda_list)]
    }
    if(!is.null(fix_lambda)){
        use_cv = FALSE
        if(fix_lambda<0){stop("The fixed lambda should be nonnegative.")}
    }

    if(tune_ratio & !remove_penalty_Z & !remove_penalty_W){
        if(is.null(ratio_list)){
            ratio_lower<-sqrt(nZ/(nZ+nZext))/2
            ratio_upper<-(nZ)^(1/3)/2
            ratio_count<-10
            ratio_list<-(seq(sqrt(ratio_lower),sqrt(ratio_upper),(sqrt(ratio_upper)-sqrt(ratio_lower))/ratio_count)^2)
            ratio_list<-c(1,ratio_list)
        }
    }else{tune_ratio<-FALSE}
    if(!use_cv){
        if(penalty_type == "ols"){
            fit_final_ols=lm(y~0+.,data = data.frame(y= pseudo_y,pseudo_X))
            beta=fit_final_ols$coefficients
            return_list<-list("beta"=beta)
        }else{
            if(!is.null(fix_lambda)){
                fit_final_fixed_lambda<-glmnet(x= pseudo_X,y= pseudo_y,standardize=F,
                                               intercept=F,alpha = final_alpha,
                                               penalty.factor = w_adaptive,
                                               lambda = fix_lambda)
                beta<-coef.glmnet(fit_final_fixed_lambda)[-1]
                return_list<-list("beta"=beta,
                                  "fix_lambda"=fix_lambda,
                                  "fix_ratio"=fix_ratio)
            }else{
                fit_final_lambda_list<-glmnet(x= pseudo_X,y= pseudo_y,standardize=F,
                                              intercept=F,alpha = final_alpha,
                                              penalty.factor = w_adaptive,
                                              lambda = lambda_list)
                return_list<-list("beta"=fit_final_lambda_list$beta,
                                  "lambda_list"=fit_final_lambda_list$lambda)
                inference=FALSE
            }
        }
    }else{
        if(length(unique(y)) <= 2){
            index_fold<-createFolds(as.numeric(y>0),k = nfolds)
        }else{index_fold<-createFolds(y,k = nfolds)}

        if(tune_ratio){
            if(family == "gaussian"){
                cv_mse<-cv_mse_lambda_ratio_func(index_fold,Z,W,A,y,
                                                 C_half,beta_initial,hat_thetaA,
                                                 study_info,lambda_list,
                                                 ratio_list,pZ,pW,pA,
                                                 w_adaptive,final_alpha,pseudo_Xy)
                ids<-which(cv_mse==min(cv_mse),arr.ind = TRUE)
            }else if(family == "binomial"){
                cv_dev<-cv_dev_lambda_ratio_func(index_fold,Z,W,A,y,
                                                 C_half,beta_initial,hat_thetaA,
                                                 study_info,lambda_list,
                                                 ratio_list,pZ,pW,pA,
                                                 w_adaptive,final_alpha,pseudo_Xy)
                if(type_measure == "auc"){
                    cv_dev1<-cv_dev$auc
                    ids<-which(cv_dev1==max(cv_dev1),arr.ind = TRUE)
                }else{
                    cv_dev1<-cv_dev$deviance
                    ids<-which(cv_dev1==min(cv_dev1),arr.ind = TRUE)
                }
            }
            final.ratio.min<-ratio_list[ids[1]]
            final.lambda.min<-lambda_list[ids[2]]
        }else{
            if(family == "gaussian"){
                cv_mse<-cv_mse_lambda_func(index_fold,Z,W,A,y,
                                           C_half,beta_initial,hat_thetaA,
                                           study_info,lambda_list,
                                           w_adaptive,final_alpha,pseudo_Xy)
                final.lambda.min<-lambda_list[which.min(cv_mse)]
            }else if(family == "binomial"){
                cv_dev<-cv_dev_lambda_func(index_fold,Z,W,A,y,
                                           C_half,beta_initial,hat_thetaA,
                                           study_info,lambda_list,
                                           w_adaptive,final_alpha,pseudo_Xy)
                if(type_measure == "auc"){
                    cv_dev1<-cv_dev$auc
                    final.lambda.min<-lambda_list[which.max(cv_dev1)]
                }else{
                    cv_dev1<-cv_dev$deviance
                    final.lambda.min<-lambda_list[which.min(cv_dev1)]
                }
            }
            final.ratio.min<-1
        }
        ratio_vec<-c(rep(final.ratio.min,pZ),rep(1,pW+pA))
        w_adaptive_ratio<-w_adaptive*ratio_vec
        fit_final_lam_ratio<-glmnet(x= pseudo_X,y= pseudo_y,standardize=F,
                                    intercept=F,alpha = final_alpha,
                                    penalty.factor = w_adaptive_ratio,
                                    lambda = final.lambda.min)

        beta<-coef.glmnet(fit_final_lam_ratio)[-1]
        return_list<-list("beta"=beta,
                          "lambda_list"=lambda_list,
                          "ratio_list"=ratio_list,
                          "lambda_min"=final.lambda.min,
                          "ratio_min"=final.ratio.min)
        if(family == "gaussian"){
            return_list<-c(return_list,
                           list("cv_mse"=cv_mse/nfolds))
        }else if(family == "binomial"){
            return_list<-c(return_list,
                           list("cv_dev"=cv_dev$deviance/nfolds,
                                "cv_auc"=cv_dev$auc/nfolds))
        }
    }

    index_nonzero<-which(beta!=0)
    if(inference & length(index_nonzero) > 1){
        if(penalty_type == "lasso"){
            warning("Current penalty is lasso, please turn to adaptivelasso for inference")
        }
        # refine C
        inv_C = Delta_opt(y=y,Z=Z,W=W,
                          family=family,
                          study_info=study_info,
                          A=A,hat_thetaA=hat_thetaA,
                          beta=beta)
        inv_C_svd<-fast.svd(inv_C+diag(1e-15,nrow(inv_C)))
        C_half<-inv_C_svd$v%*%diag(1/sqrt(inv_C_svd$d))%*%t(inv_C_svd$u)

        pseudo_Xy_list<-pseudo_Xy(C_half,Z,W,A,y,beta,hat_thetaA,study_info)
        Sigsum_half<-pseudo_Xy_list$pseudo_X/nZ

        #expit_beta<-c(expit(X%*%beta))
        #dexpit_beta<-expit_beta*(1-expit_beta)
        #pseudo_X<-C_half%*%rbind(t(X),t(Z))%*%(X*c(dexpit_beta))/nZ

        ## gaussian: Sigsum_half<-cbind(ZWtZW/nZ,crossprod(ZW,Z)/nZ)%*%C_half

        Sigsum_scaled<-crossprod(Sigsum_half)
        Sigsum_scaled_nonzero<-Sigsum_scaled[index_nonzero,index_nonzero]
        inv_Sigsum_scaled_nonzero<-solve(Sigsum_scaled_nonzero)
        final_v<-diag(inv_Sigsum_scaled_nonzero)/nZ

        pval_final<-pchisq(beta[index_nonzero]^2/final_v,1,lower.tail = F)
        pval_final1<-p.adjust(pval_final,method = "BH")
        selected_pos<-index_nonzero[which(pval_final1<0.05)]
        return_list<-c(return_list,
                       list("selected_vars"=
                                list("position"=index_nonzero,
                                     "coef"=beta[index_nonzero],
                                     "pval"=pval_final,
                                     "variance"=final_v,
                                     "FDR_adjust_position"=selected_pos)
                       ))
    }
    return(return_list)
}