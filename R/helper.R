mu_func<-function(x,family){
    if(family=='gaussian'){return(x)}
    else if(family=='binomial'){return(expit(x))}
}
Delta_opt<-function(y,Z,W,family,
                    study_info,A=NULL,pA=NULL,
                    beta=NULL,hat_thetaA=NULL,
                    V_thetaA=NULL,use_offset=TRUE){
    X=cbind(A,Z,W)
    XR=cbind(A,Z)
    n_main=length(y)
    tilde_thetaZ=study_info[[1]]$Coeff
    tilde_theta=c(hat_thetaA,tilde_thetaZ)
    mu_X_beta=mu_func(X%*%beta,family) # check
    mu_XR_theta=mu_func(XR%*%tilde_theta,family)
    if(family == "binomial"){mu_prime_XR_theta=mu_XR_theta*(1-mu_XR_theta)}
    else{mu_prime_XR_theta=1}
    V_U1=(1/n_main)*crossprod(X*c(mu_X_beta-y))
    V_U2=(1/n_main)*crossprod(Z*c(mu_X_beta-mu_XR_theta))
    Cov_U1U2=(1/n_main)*crossprod(X*c(mu_X_beta-y),Z*c(mu_X_beta-mu_XR_theta))
    GammaZZ=(1/n_main)*crossprod(Z*c(mu_prime_XR_theta),Z)
    V_thetaZ=study_info[[1]]$Covariance
    Delta22=V_U2 +GammaZZ%*%(n_main*V_thetaZ)%*%t(GammaZZ)
    Delta12=Cov_U1U2
    if(pA!=0){
        GammaZA=(1/n_main)*crossprod(Z*c(mu_prime_XR_theta),A)
        if(use_offset){
          inv_GammaAA=ginv((1/n_main)*crossprod(A*c(mu_prime_XR_theta),A))
          #V_thetaA = (1/n_main)*inv_GammaAA%*%((1/n_main)* crossprod(A*c(mu_XR_theta-y)) )%*%inv_GammaAA
          Cov_U1theta=(1/n_main)*crossprod(X*c(mu_X_beta-y),A*c(mu_XR_theta-y))%*%(inv_GammaAA%*%t(GammaZA))
          Cov_U2theta=(1/n_main)*crossprod(Z*c(mu_X_beta-mu_XR_theta),A*c(mu_XR_theta-y))%*%(inv_GammaAA%*%t(GammaZA))
        }else{
          inv_GammaXRXR=ginv((1/n_main)*crossprod(XR*c(mu_prime_XR_theta),XR))
          #V_thetaA = (1/n_main)*inv_GammaXRXR[1:pA,]%*%((1/n_main)* crossprod(XR*c(mu_XR_theta-y)) )%*%inv_GammaXRXR[,1:pA]
          Cov_U1theta=(1/n_main)*crossprod(X*c(mu_X_beta-y),XR*c(mu_XR_theta-y))%*%(inv_GammaXRXR[,1:pA]%*%t(GammaZA))
          Cov_U2theta=(1/n_main)*crossprod(Z*c(mu_X_beta-mu_XR_theta),XR*c(mu_XR_theta-y))%*%(inv_GammaXRXR[,1:pA]%*%t(GammaZA))
        }

        Delta22 = Delta22 + GammaZA%*%(n_main*V_thetaA)%*%t(GammaZA)
        + Cov_U2theta+t(Cov_U2theta)
        Delta12 = Delta12 + Cov_U1theta
    }

    Delta = rbind(cbind(V_U1,Delta12),cbind(t(Delta12),Delta22))
    Delta
}

vcov_sandwich<-function(y,A,Z,family,study_info,pA,
                        hat_thetaA,use_offset){
    n_main = length(y)
    XR = cbind(A,Z)
    tilde_thetaZ=study_info[[1]]$Coeff
    tilde_theta=c(hat_thetaA,tilde_thetaZ)
    mu_XR_theta=mu_func(cbind(A,Z)%*%tilde_theta,family)
    if(family == "binomial"){mu_prime_XR_theta=mu_XR_theta*(1-mu_XR_theta)}
    else{mu_prime_XR_theta=1}
    if(use_offset){
        inv_GammaAA=ginv((1/n_main)*crossprod(A*c(mu_prime_XR_theta),A))
        V_thetaA = (1/n_main)*inv_GammaAA%*%((1/n_main)* crossprod(A*c(mu_XR_theta-y)) )%*%inv_GammaAA
    }else{
        inv_GammaXRXR=ginv((1/n_main)*crossprod(XR*c(mu_prime_XR_theta),XR))
        V_thetaA = (1/n_main)*inv_GammaXRXR[1:pA,]%*%((1/n_main)* crossprod(XR*c(mu_XR_theta-y)) )%*%inv_GammaXRXR[,1:pA]
    }
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
        Ztrain<-Z[-index_test,,drop=FALSE]
        Ztest<-Z[index_test,,drop=FALSE]
        if(!is.null(W)){
            Wtrain<-W[-index_test,,drop=FALSE]
            Wtest<-W[index_test,,drop=FALSE]
        }else{
            Wtrain<-NULL
            Wtest<-NULL}
        if(!is.null(A)){
            Atrain<-A[-index_test,,drop=FALSE]
            Atest<-A[index_test,,drop=FALSE]
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
        Ztrain<-Z[-index_test,,drop=FALSE]
        Ztest<-Z[index_test,,drop=FALSE]
        if(!is.null(W)){
            Wtrain<-W[-index_test,,drop=FALSE]
            Wtest<-W[index_test,,drop=FALSE]
        }else{
            Wtrain<-NULL
            Wtest<-NULL}
        if(!is.null(A)){
            Atrain<-A[-index_test,,drop=FALSE]
            Atest<-A[index_test,,drop=FALSE]
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
        Ztrain<-Z[-index_test,,drop=FALSE]
        Ztest<-Z[index_test,,drop=FALSE]
        if(!is.null(W)){
            Wtrain<-W[-index_test,,drop=FALSE]
            Wtest<-W[index_test,,drop=FALSE]
        }else{
            Wtrain<-NULL
            Wtest<-NULL}
        if(!is.null(A)){
            Atrain<-A[-index_test,,drop=FALSE]
            Atest<-A[index_test,,drop=FALSE]
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
        Ztrain<-Z[-index_test,,drop=FALSE]
        Ztest<-Z[index_test,,drop=FALSE]
        if(!is.null(W)){
            Wtrain<-W[-index_test,,drop=FALSE]
            Wtest<-W[index_test,,drop=FALSE]
        }else{
            Wtrain<-NULL
            Wtest<-NULL}
        if(!is.null(A)){
            Atrain<-A[-index_test,,drop=FALSE]
            Atest<-A[index_test,,drop=FALSE]
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
        V_thetaA = NULL,
        use_offset = TRUE,
        V_thetaA_sandwich = TRUE,
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
    if(!penalty_type %in% c("none","adaptivelasso","lasso","ridge")){
        stop("Select penalty type from c('none','adaptivelasso','lasso','ridge').")
    }
    if(!type_measure%in% c("default", "mse", "deviance", "auc")){
        stop("Select type_measure from c('default','deviance','auc')")
    }
    if(is.null(dim(Z)[1])){
        warning("Z is input as a vector, convert Z into matrix with size nZ*1")
        Z=matrix(Z,ncol=1)
    }
    final_alpha = 1
    if(penalty_type == "ridge"){final_alpha = 0}

    if(!is.null(fix_ratio)){
        if(tune_ratio){
            stop("If ratio is fixed, please set tune_ratio as FALSE")
        }else if(remove_penalty_Z | remove_penalty_W){
            stop("If ratio is fixed, please set remove_penalty's as FALSE")
        }
    }

    nZ=nrow(Z)
    nZext=study_info[[1]]$Sample_size
    pZ=ncol(Z)
    if(is.null(W)){pW=0}else{pW=ncol(W)}
    if(is.null(A)){pA=0}else{
        if(is.null(dim(A)[1])){
            if(length(A)==1){
                if(A==1){A=matrix(1,nrow=nZ,ncol=1)}
            }
        }
        pA=ncol(A)
    }
    if(nZ<2*pZ+pW+pA){use_sparseC=TRUE}

    if(family == "gaussian"){pseudo_Xy=pseudo_Xy_gaussian
    }else if(family == "binomial"){pseudo_Xy=pseudo_Xy_binomial}


    Aid<-1:pA
    Zid<-(pA+1):(pA+pZ)
    Wid<-(pA+pZ+1):(pA+pZ+pW)
    Acolnames=NULL
    if(pA>0){
        Acolnames=colnames(A)
        if(is.null(Acolnames[1])){
            Acolnames=paste0('A',1:pA)
        }
        if(length(unique(A[,1])) == 1){
            if(unique(A[,1]) == 1){
                Acolnames[1]='intercept'
            }
        }
        colnames(A)=Acolnames
    }
    Zcolnames=colnames(Z)
    Wcolnames=colnames(W)
    if(is.null(Zcolnames[1])){
        Zcolnames=paste0('Z',1:pZ)
        colnames(Z)=Zcolnames
    }
    if(is.null(Wcolnames[1])){
        Wcolnames=paste0('W',1:pW)
        colnames(W)=Wcolnames
    }
    Xcolnames<-c(Acolnames,Zcolnames,Wcolnames)
    X<-cbind(A,Z,W)

    if(pA!=0){
        if(is.null(hat_thetaA)){
            if(!is.null(V_thetaA)){
                stop("With customized hat_thetaA input, V_thetaA is also needed")
            }
            if(use_offset){
                offset_term = c(Z%*%study_info[[1]]$Coeff)
                df=data.frame(y,A)
                if(family=="binomial"){
                    hat_thetaA_glm=speedglm(y~0+.,data = df,offset = offset_term,family = binomial())
                }else if(family=="gaussian"){
                    hat_thetaA_glm=speedlm(y~0+.,data = df,offset = offset_term)
                    hat_thetaA_glm=lm(y~0+.,data = df,offset = offset_term)
                }
                hat_thetaA=hat_thetaA_glm$coefficients
                if(V_thetaA_sandwich){
                    V_thetaA=vcov_sandwich(y,A,Z,family,study_info,pA,
                                           hat_thetaA,use_offset)
                }else{V_thetaA=vcov(hat_thetaA_glm)}
            }else{
                df=data.frame(y,A,Z)
                if(family=="binomial"){
                    hat_thetaA_glm=speedglm(y~0+.,data = df,family = binomial())
                }else if(family=="gaussian"){
                    hat_thetaA_glm=speedlm(y~0+.,data = df)
                }
                hat_thetaA=hat_thetaA_glm$coefficients[1:pA]
                if(V_thetaA_sandwich){
                    V_thetaA=vcov_sandwich(y,A,Z,family,study_info,pA,
                                           hat_thetaA,use_offset)
                }else{V_thetaA=vcov(hat_thetaA_glm)[1:pA,1:pA]}
            }
        }
    }


    fix_penalty<-c(rep(0,pA),rep(1,pZ+pW))
    if(remove_penalty_Z){fix_penalty[Zid]<-0}else{
        if(!is.null(fix_ratio)){fix_penalty[Zid]<-fix_ratio}
    }
    if(remove_penalty_W){fix_penalty[Wid]<-0}
    if((remove_penalty_Z & remove_penalty_W)|(length(unique(fix_penalty))==1 & unique(fix_penalty)[1] == 0) ){
        penalty_type = "none"
        warning("All penalties are removed, turn to no penalties!")
    }
    if(penalty_type == "none"){
        initial_with_type = "glm"
        use_cv = FALSE
        tune_ratio = FALSE
    }
    if(!is.null(beta_initial) & length(beta_initial)!=pA+pZ+pW){
        warning("beta_initial should be from A,Z,W.\n Length not match, compute default initial instead.")
        beta_initial=NULL
    }
    if(is.null(beta_initial)){
        if(initial_with_type %in% c("glm","ridge","lasso")){
            if(initial_with_type == "ridge"){initial_alpha=0}else{initial_alpha=1}
            if(initial_with_type == "glm"){
                df=data.frame(y,X)
                if(family=="binomial"){
                    fit_initial=speedglm(y~0+.,data = df,family = binomial())
                }else if(family=="gaussian"){
                    fit_initial=speedlm(y~0+.,data = df)
                }
                beta_initial=fit_initial$coefficients
            }else if(pA == 0){
                fit_initial=cv.glmnet(x=X,y= y,
                                      alpha = initial_alpha,
                                      penalty.factor = fix_penalty,
                                      family=family)
                beta_initial=c(coef.glmnet(fit_initial,s="lambda.min")[-1])
            }else if(length(unique(A[,1]))==1){
                if(unique(A[,1])==1){
                    fit_initial=cv.glmnet(x=X[,-1],y=y,
                                          alpha = initial_alpha,
                                          penalty.factor = fix_penalty[-1],
                                          family=family)
                    beta_initial=as.vector(coef.glmnet(fit_initial,s="lambda.min"))
                }else{
                    stop("The first column of A is constant, then it should be 1 for intercept.")
                }
            }else{
                fit_initial=cv.glmnet(x=X,y= y,
                                      alpha = initial_alpha,
                                      penalty.factor = fix_penalty,
                                      family=family)
                beta_initial=c(coef.glmnet(fit_initial,s="lambda.min")[-1])
            }
        }else{stop("Select Initial Type from c('glm','ridge','lasso')")}
    }
    if (penalty_type == "adaptivelasso"){
        w_adaptive<-1/abs(beta_initial)^gamma_adaptivelasso
        w_adaptive[is.infinite(w_adaptive)]<-max(w_adaptive[!is.infinite(w_adaptive)])*100
        w_adaptive<-w_adaptive*fix_penalty
    }else{w_adaptive<-fix_penalty}
    # Estimation of C
    inv_C = Delta_opt(y=y,Z=Z,W=W,
                      family=family,
                      study_info=study_info,
                      A=A,pA=pA,beta=beta_initial,
                      hat_thetaA=hat_thetaA,
                      V_thetaA = V_thetaA)
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

    # generate lambda list from glmnet
    if(penalty_type != "none"){
        if(is.null(fix_lambda)&is.null(lambda_list)){
            fit_final<-glmnet(x= pseudo_X,y= pseudo_y,standardize=F,
                              intercept=F,alpha = final_alpha,penalty.factor = w_adaptive)
            lambda_list<-fit_final$lambda
            lambda_list<-lambda_list[!is.na(lambda_list)]
        }
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
        if(penalty_type == "none"){
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
                                  "fix_lambda"=fix_lambda)
            }else{
                fit_final_lambda_list<-glmnet(x= pseudo_X,y= pseudo_y,standardize=F,
                                              intercept=F,alpha = final_alpha,
                                              penalty.factor = w_adaptive,
                                              lambda = lambda_list)
                return_list<-list("beta"=fit_final_lambda_list$beta,
                                  "lambda_list"=fit_final_lambda_list$lambda)
                if(inference){warnings("When use_cv=F,fix_lambda is NULL, no inference will be done")}
                inference=FALSE
            }
            if(!is.null(fix_ratio)){
                return_list<-c(return_list,
                               list("fix_ratio"=fix_ratio))
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
    if(inference){
        index_nonzero<-which(beta!=0)
    if(length(index_nonzero) > 1){
        if(penalty_type == "lasso"){
            warning("Current penalty is lasso, please turn to adaptivelasso for inference")
        }
        # refine C
        inv_C = Delta_opt(y=y,Z=Z,W=W,
                          family=family,
                          study_info=study_info,
                          A=A,pA=pA,beta=beta_initial,
                          hat_thetaA=hat_thetaA,
                          V_thetaA = V_thetaA)
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
                                     "name"=Xcolnames[index_nonzero],
                                     "coef"=beta[index_nonzero],
                                     "variance"=final_v,
                                     "pval"=pval_final,
                                     "FDR_adjust_position"=selected_pos,
                                     "FDR_adjust_name"=Xcolnames[selected_pos])
                       ))
    }}
    return(return_list)
}
