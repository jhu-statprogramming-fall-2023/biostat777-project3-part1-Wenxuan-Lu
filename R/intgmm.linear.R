#pacman::p_load(expm,magic,glmnet,cluster,MASS,locfit,corpcor,caret,mvtnorm)
pseudo_Xy_multiv_addition<-function(
        C_half,X,A,y,study_info){
    XA<-cbind(X,A)
    xatxa<-crossprod(XA)
    xtxa<-crossprod(X,XA)
    pseudo_X<-C_half%*%rbind(xatxa,xtxa)
    pseudo_y1<-crossprod(y,XA)
    pseudo_y2<-crossprod(X)%*%study_info[[1]]$Coeff
    pseudo_y<-c(c(pseudo_y1,pseudo_y2)%*%C_half)
    list("pseudo_X"=pseudo_X,"pseudo_y"=pseudo_y)
}
# work for computation of weighting matrix C
var_U_beta_theta_func_multiv_addition<-function(
        X,A,y,beta,
        study_info){
    XA<-cbind(X,A)
    xatbeta<-(XA%*%beta)
    xttheta<-(X%*%study_info[[1]]$Coeff)
    var_11<-crossprod(XA*c(xatbeta-y))
    var_22<-crossprod(X*c(xatbeta-xttheta))
    var_12<-crossprod(XA*c(xatbeta-y),X*c(xatbeta-xttheta))
    (1/nrow(X))*rbind(cbind(var_11,var_12),cbind(t(var_12),var_22))
}

pseudo_Xy_univ_addition<-function(
        C_half,X,A,y,study_info){
    XA<-cbind(X,A)
    xatxa<-crossprod(XA)
    xtxa<-crossprod(X,XA)
    pseudo_X<-C_half%*%rbind(xatxa,xtxa)
    pseudo_y1<-crossprod(y,XA)
    pseudo_y2<-sapply(1:ncol(X), function(id){
        u2_id<-c(X[,id]*study_info[[id]]$Coeff)
        c(u2_id%*%X[,id])
    })
    pseudo_y<-c(c(pseudo_y1,pseudo_y2)%*%C_half)
    list("pseudo_X"=pseudo_X,"pseudo_y"=pseudo_y)
}
# work for computation of weighting matrix C
var_U_beta_theta_func_univ_addition<-function(
        X,A,y,beta,
        study_info){
    XA<-cbind(X,A)
    xatbeta<-(XA%*%beta)
    var_11<-crossprod(XA*c(xatbeta-y))
    u2_theta_coef<-sapply(1:ncol(X), function(id){
        xtgamma_id<-c(X[,id]*study_info[[id]]$Coeff)
        xatbeta-xtgamma_id
    }) #row is sample #col is SNP
    var_22<-crossprod(u2_theta_coef*X)
    var_12<-crossprod(XA*c(xatbeta-y),u2_theta_coef*X)
    (1/nrow(X))*rbind(cbind(var_11,var_12),cbind(t(var_12),var_22))
}




## inside helper function
cv_mse_lambda_ratio_func<-function(index_fold,X,A,y,C_half,
                                   study_info,lambda_list,
                                   ratio_range,pX,pA,w_adaptive,
                                   final_alpha,pseudo_Xy){
    fold_mse_lambda_ratio<-lapply(1:length(index_fold), function(cur_fold){
        index_test<-index_fold[[cur_fold]]
        Xtrain<-X[-index_test,]
        Xtest<-X[index_test,]
        if(!is.null(A)){
            Atrain<-A[-index_test,]
            Atest<-A[index_test,]
        }else{
            Atrain<-NULL
            Atest<-NULL}
        ytrain<-y[-index_test]
        ytest<-y[index_test]
        pseudo_Xy_list_train<-pseudo_Xy(C_half,Xtrain,Atrain,
                                        ytrain,study_info)
        initial_sf_train<-nrow(Xtrain)/sqrt(nrow(pseudo_Xy_list_train$pseudo_X))
        pseudo_X_train<-pseudo_Xy_list_train$pseudo_X/initial_sf_train
        pseudo_y_train<-pseudo_Xy_list_train$pseudo_y/initial_sf_train
        mse_lam_ratio_fold<-sapply(lambda_list,function(cur_lam){
            sapply(ratio_range,function(cur_ratio){
                ratio_vec<-c(rep(cur_ratio,pX),rep(1,pA))
                w_adaptive_ratio<-w_adaptive*ratio_vec
                cv_fit<-glmnet(x=pseudo_X_train,y=pseudo_y_train,
                               standardize=F,intercept=F,alpha = final_alpha,
                               penalty.factor = w_adaptive_ratio,
                               lambda = cur_lam)
                cur_beta<-coef.glmnet(cv_fit)[-1]
                mean((cbind(Xtest,Atest)%*%cur_beta - ytest)^2)
            })
        }) # row is ratio_range & col is lambda_list
        mse_lam_ratio_fold
    })
    cv_mse_lambda_ratio<-Reduce(`+`, fold_mse_lambda_ratio)
    cv_mse_lambda_ratio
}

# helper function 2
cv_mse_lambda_func<-function(index_fold,X,A,y,
                             C_half,study_info,
                             lambda_list,final_alpha,
                             w_adaptive,pseudo_Xy){
    fold_mse_lambda<-sapply(1:length(index_fold), function(cur_fold){
        index_test<-index_fold[[cur_fold]]
        Xtrain<-X[-index_test,]
        Xtest<-X[index_test,]
        if(!is.null(A)){
            Atrain<-A[-index_test,]
            Atest<-A[index_test,]
        }else{
            Atrain<-NULL
            Atest<-NULL}
        ytrain<-y[-index_test]
        ytest<-y[index_test]
        pseudo_Xy_list_train<-pseudo_Xy(C_half,Xtrain,Atrain,
                                        ytrain,study_info)
        initial_sf_train<-nrow(Xtrain)/sqrt(nrow(pseudo_Xy_list_train$pseudo_X))
        pseudo_X_train<-pseudo_Xy_list_train$pseudo_X/initial_sf_train
        pseudo_y_train<-pseudo_Xy_list_train$pseudo_y/initial_sf_train
        mse_lam<-sapply(lambda_list,function(cur_lam){
            cv_fit<-glmnet(x= (pseudo_X_train),y= (pseudo_y_train),
                           standardize=F,intercept=F,
                           alpha = final_alpha,penalty.factor = w_adaptive,lambda = cur_lam)
            cur_beta<-coef.glmnet(cv_fit)[-1]
            mean((cbind(Xtest,Atest)%*%cur_beta - ytest)^2)
        })
        mse_lam
    })
    cv_mse_lambda<-rowMeans(fold_mse_lambda)
    cv_mse_lambda
}

# helper function 3
holdout_mse_lambda_func<-function(lambda_list,pseudo_X_train,pseudo_y_train,
                                  final_alpha,w_adaptive,Xtest,Atest,ytest){
    mse_lam<-sapply(lambda_list,function(cur_lam){
        cv_fit<-glmnet(x= (pseudo_X_train),y= (pseudo_y_train),
                       standardize=F,intercept=F,
                       alpha = final_alpha,penalty.factor = w_adaptive,lambda = cur_lam)
        cur_beta<-coef.glmnet(cv_fit)[-1]
        mean((cbind(Xtest,Atest)%*%cur_beta - ytest)^2)
    })
    mse_lam
}

#helper function 4
holdout_mse_lambda_ratio_func<-function(lambda_list,ratio_range,pX,pA,
                                        w_adaptive,pseudo_X_train,pseudo_y_train,
                                        final_alpha,Xtest,Atest,ytest){
    mse_lam_ratio<-sapply(lambda_list,function(cur_lam){
        sapply(ratio_range,function(cur_ratio){
            ratio_vec<-c(rep(cur_ratio,pX),rep(1,pA))
            w_adaptive_ratio<-w_adaptive*ratio_vec
            cv_fit<-glmnet(x=pseudo_X_train,y=pseudo_y_train,
                           standardize=F,intercept=F,alpha = final_alpha,
                           penalty.factor = w_adaptive_ratio,
                           lambda = cur_lam)
            cur_beta<-coef.glmnet(cv_fit)[-1]
            mean((cbind(Xtest,Atest)%*%cur_beta - ytest)^2)
        })
    }) # row is ratio_range & col is lambda_list
    mse_lam_ratio
}

intgmm.linear<-function(
        y,X,A=NULL,
        study_info=NULL,
        summary_type = "multi",
        penalty_type = "lasso",
        initial_with_type = "ridge",
        beta_initial = NULL,
        remove_penalty_X = FALSE,
        remove_penalty_A = FALSE,
        tune_ratio = TRUE,
        fix_lambda = NULL,
        lambda_list = NULL,
        fix_ratio = NULL,
        ratio_lower = NULL,
        ratio_upper = NULL,
        ratio_count = 10,
        ratio_range = NULL,
        gamma_adaptivelasso = 1/2,
        inference = FALSE,
        validation_type = "cv",
        nfolds = 10,
        holdout_p = 0.2,
        use_sparseC = FALSE,
        shrink = FALSE,
        seed.use = 97
){
    set.seed(seed.use)
    if (is.null(study_info)){stop("Please input study_info as summary statistics")}
    if (!summary_type %in% c("uni","multi")){
        stop("Select summary_type for input summary statistics(study_info).
       Use 'uni' for univariate summary statistics.
       Use 'multi' for multivariate summary statistics.")
    }
    if(!penalty_type %in% c("adaptivelasso","lasso","ridge")){
        stop("Select penalty type from c('adaptivelasso','lasso','ridge').")
    }
    final_alpha = 1
    if(penalty_type == "ridge"){final_alpha = 0}
    if (summary_type == "uni"){
        pseudo_Xy <- pseudo_Xy_univ_addition
        var_U_beta_theta_func <- var_U_beta_theta_func_univ_addition
    }else if (summary_type == "multi"){
        pseudo_Xy <- pseudo_Xy_multiv_addition
        var_U_beta_theta_func <- var_U_beta_theta_func_multiv_addition
    }

    if(!is.null(fix_ratio)){
        if(tune_ratio){
            stop("If ratio is fixed, please set tune_ratio as FALSE")
        }else if(remove_penalty_X | remove_penalty_A){
            stop("If ratio is fixed, please set remove_penalty's as FALSE")
        }
    }
    # X and y should be scaled to avoid potential mismatch of intercept term.
    # n for sample size; p for feature size
    nX<-nrow(X)
    nXext<-study_info[[1]]$Sample_size
    pX<-ncol(X)
    if(is.null(A)){pA<-0}else{pA<-ncol(A)}
    if(nX<2*pX+pA){use_sparseC<-TRUE}
    XA<-cbind(X,A)
    Xid<-1:pX
    Aid<-(pX+1):(pX+pA)
    fix_penalty<-rep(1,pX+pA)
    if(remove_penalty_X){fix_penalty[Xid]<-0}
    if(remove_penalty_A){fix_penalty[Aid]<-0}
    if(!is.null(fix_ratio)){fix_penalty[Xid]<-fix_ratio}
    xatxa<-crossprod(XA)
    xtx<-crossprod(X)
    if (summary_type == "uni"){
        diag_theta_sd<-sapply(1:pX,function(id){
            sqrt(study_info[[id]]$Covariance)*c(X[,id]%*%X[,id])})
        var_U2<-t(diag_theta_sd*cor(X))*diag_theta_sd/nX
    }else if (summary_type == "multi"){
        var_U2<-xtx%*%study_info[[1]]$Covariance%*%xtx/nX
    }
    if(!is.null(beta_initial)){
        if(length(beta_initial)!=pX+pA){
            warning("beta_initial should be from X,A,G.\n Length not match, compute default initial instead.")
            beta_initial<-NULL
        }
    }

    if(is.null(beta_initial) & initial_with_type %in% c("ridge","lasso")){
        if(initial_with_type == "ridge"){initial_alpha=0}else{initial_alpha=1}
        fit_initial<-cv.glmnet(x=XA,y= y,alpha = initial_alpha,penalty.factor = fix_penalty)
        beta_initial<-c(coef.glmnet(fit_initial,s="lambda.min")[-1])
    }else if(is.null(beta_initial)){
        # Initial estimation of C by GMM
        var_U1<-xatxa/nX*c(var(y))
        C_11_half<-diag(1/sqrt(diag(var_U1)))
        C_22_half<-diag(1/sqrt(diag(var_U2)))
        C_half<-adiag(C_11_half,C_22_half)

        # Prepare for initial model
        pseudo_Xy_list<-pseudo_Xy(C_half,X,A,y,study_info)
        initial_sf<-nX/sqrt(nrow(pseudo_Xy_list$pseudo_X))
        pseudo_X<-pseudo_Xy_list$pseudo_X/initial_sf
        pseudo_y<-pseudo_Xy_list$pseudo_y/initial_sf
        fit_initial<-cv.glmnet(x= pseudo_X,y= pseudo_y,standardize=F,intercept=F,alpha = 0,penalty.factor = fix_penalty)
        beta_initial<-c(coef.glmnet(fit_initial,s ='lambda.min')[-1])
    }
    if (penalty_type == "adaptivelasso"){
        w_adaptive<-1/abs(beta_initial)^gamma_adaptivelasso
        w_adaptive[is.infinite(w_adaptive)]<-max(w_adaptive[!is.infinite(w_adaptive)])*100
        w_adaptive<-w_adaptive*fix_penalty
    }else{w_adaptive<-fix_penalty}
    # Refined estimation of C
    var_1st_U_beta_theta<-var_U_beta_theta_func(X = X,A = A, y = y,beta = beta_initial,
                                                study_info = study_info)
    var_2nd_grad_times_theta_hat = adiag(matrix(0,nrow = pX+pA,ncol = pX+pA),var_U2)
    inv_C = var_1st_U_beta_theta + var_2nd_grad_times_theta_hat

    if(use_sparseC){
        C_half<-diag(1/sqrt(diag(inv_C)))
    }else{
        inv_C_svd<-fast.svd(inv_C+diag(1e-15,nrow(inv_C)))
        C_half<-inv_C_svd$v%*%diag(1/sqrt(inv_C_svd$d))%*%t(inv_C_svd$u)
    }

    pseudo_Xy_list<-pseudo_Xy(C_half,X,A,y,study_info)
    initial_sf<-nX/sqrt(nrow(pseudo_Xy_list$pseudo_X))
    pseudo_X<-pseudo_Xy_list$pseudo_X/initial_sf
    pseudo_y<-pseudo_Xy_list$pseudo_y/initial_sf

    # Fit final model
    fit_final<-cv.glmnet(x= pseudo_X,y= pseudo_y,standardize=F,
                         intercept=F,alpha = final_alpha,penalty.factor = w_adaptive)
    if(is.null(lambda_list)){lambda_list<-fit_final$lambda
    if(shrink){lambda_list<-lambda_list[1:50]}}
    if(!is.null(fix_lambda)){
        validation_type<-"None"
        if(fix_lambda<0){stop("The fixed lambda should be nonnegative.")}
    }

    if(tune_ratio & !remove_penalty_X & !remove_penalty_A){
        if(is.null(ratio_range)){
            if(is.null(ratio_lower)){ratio_lower<-sqrt(nX/(nX+nXext))/2}
            if(is.null(ratio_upper)){ratio_upper<-(nX)^(1/3)/2}
            #ratio_range<-exp(seq(log(ratio_lower),log(ratio_upper),(log(ratio_upper)-log(ratio_lower))/ratio_count))
            ratio_range<-(seq(sqrt(ratio_lower),sqrt(ratio_upper),(sqrt(ratio_upper)-sqrt(ratio_lower))/ratio_count)^2)
            ratio_range<-c(1,ratio_range)
        }
    }else{tune_ratio<-F}
    if(validation_type == "None"){
        fit_final_fixed_lambda<-glmnet(x= pseudo_X,y= pseudo_y,standardize=F,
                                    intercept=F,alpha = final_alpha,
                                    penalty.factor = w_adaptive,
                                    lambda = fix_lambda)
        beta<-coef.glmnet(fit_final_fixed_lambda)[-1]
        return_list<-list("beta"=beta,
                          "fix_lambda"=fix_lambda,
                          "fix_ratio"=fix_ratio)
    }else if(validation_type == "holdout"){
        if(length(unique(y)) <= 2){
            index_valid<-createDataPartition(as.numeric(y>0),p = holdout_p)$Resample1
        }else{index_valid<-createDataPartition(y,p = holdout_p)$Resample1}
        Xtrain<-X[-index_valid,]
        Xtest<-X[index_valid,]
        if(!is.null(A)){
            Atrain<-A[-index_valid,]
            Atest<-A[index_valid,]
        }else{
            Atrain<-NULL
            Atest<-NULL}
        ytrain<-y[-index_valid]
        ytest<-y[index_valid]

        pseudo_Xy_list_train<-pseudo_Xy(C_half,Xtrain,Atrain,
                                        ytrain,study_info)
        initial_sf_train<-nrow(Xtrain)/sqrt(nrow(pseudo_Xy_list_train$pseudo_X))
        pseudo_X_train<-pseudo_Xy_list_train$pseudo_X/initial_sf_train
        pseudo_y_train<-pseudo_Xy_list_train$pseudo_y/initial_sf_train
        if(tune_ratio){
            holdout_mse<-holdout_mse_lambda_ratio_func(lambda_list,ratio_range,pX,pA,
                                                       w_adaptive,pseudo_X_train,pseudo_y_train,
                                                       final_alpha,Xtest,Atest,ytest)

            ids<-which(holdout_mse==min(holdout_mse),arr.ind = TRUE)
            final.ratio.min<-ratio_range[ids[1]]
            final.lambda.min<-lambda_list[ids[2]]
        }else{
            holdout_mse<-holdout_mse_lambda_func(lambda_list,pseudo_X_train,pseudo_y_train,
                                                 final_alpha,w_adaptive,Xtest,Atest,ytest)
            final.lambda.min<-lambda_list[which.min(holdout_mse)]
            final.ratio.min<-1
        }
        ratio_vec<-c(rep(final.ratio.min,pX),rep(1,pA))
        w_adaptive_ratio<-w_adaptive*ratio_vec
        fit_final_lam_ratio<-glmnet(x= pseudo_X_train,y= pseudo_y_train,standardize=F,
                                    intercept=F,alpha = final_alpha,
                                    penalty.factor = w_adaptive_ratio,
                                    lambda = final.lambda.min)
        beta<-coef.glmnet(fit_final_lam_ratio)[-1]
        return_list<-list("beta"=beta,
                          "lambda_list"=lambda_list,
                          "ratio_list"=ratio_range,
                          "holdout_mse"=holdout_mse,
                          "lambda_min"=final.lambda.min,
                          "ratio_min"=final.ratio.min)
    }else if(validation_type == "cv"){
        ## detailed cross validation
        if(length(unique(y)) <= 2){
            index_fold<-createFolds(as.numeric(y>0),k = nfolds)
        }else{index_fold<-createFolds(y,k = nfolds)}
        if(tune_ratio){
            cv_mse<-cv_mse_lambda_ratio_func(index_fold,X,A,y,C_half,
                                             study_info,lambda_list,
                                             ratio_range,pX,pA,w_adaptive,
                                             final_alpha,pseudo_Xy)
            ids<-which(cv_mse==min(cv_mse),arr.ind = TRUE)
            final.ratio.min<-ratio_range[ids[1]]
            final.lambda.min<-lambda_list[ids[2]]
        }else{
            cv_mse<-cv_mse_lambda_func(index_fold,X,A,y,
                                       C_half,study_info,
                                       lambda_list,final_alpha,
                                       w_adaptive,pseudo_Xy)
            final.lambda.min<-lambda_list[which.min(cv_mse)]
            final.ratio.min<-1
        }
        ratio_vec<-c(rep(final.ratio.min,pX),rep(1,pA))
        w_adaptive_ratio<-w_adaptive*ratio_vec
        fit_final_lam_ratio<-glmnet(x= pseudo_X,y= pseudo_y,standardize=F,
                                    intercept=F,alpha = final_alpha,
                                    penalty.factor = w_adaptive_ratio,
                                    lambda = final.lambda.min)
        beta<-coef.glmnet(fit_final_lam_ratio)[-1]
        return_list<-list("beta"=beta,
                          "lambda_list"=lambda_list,
                          "ratio_list"=ratio_range,
                          "cv_mse"=cv_mse,
                          "lambda_min"=final.lambda.min,
                          "ratio_min"=final.ratio.min)
    }else{stop("Validation Type should be chosen from c('cv','holdout'). Or use fix_lambda without validation.")}

    index_nonzero<-which(beta!=0)
    if(inference & length(index_nonzero) > 1){
        # Use final estimated beta to refine the C
        var_1st_U_beta_theta<-var_U_beta_theta_func(X = X,A = A, y = y,beta = beta,
                                                    study_info = study_info)
        var_2nd_grad_times_theta_hat = adiag(matrix(0,nrow = pX+pA,ncol = pX+pA),var_U2)
        inv_C = var_1st_U_beta_theta + var_2nd_grad_times_theta_hat

        inv_C_svd<-fast.svd(inv_C+diag(1e-15,nrow(inv_C)))
        C_half<-inv_C_svd$v%*%diag(1/sqrt(inv_C_svd$d))%*%t(inv_C_svd$u)


        Sigsum_half<-cbind(xatxa/nX,crossprod(XA,X)/nX)%*%C_half
        Sigsum_scaled<-Sigsum_half%*%t(Sigsum_half)
        Sigsum_scaled_nonzero<-Sigsum_scaled[index_nonzero,index_nonzero]
        inv_Sigsum_scaled_nonzero<-solve(Sigsum_scaled_nonzero)
        final_v<-diag(inv_Sigsum_scaled_nonzero)

        pval_final<-pchisq(nX*beta[index_nonzero]^2/final_v,1,lower.tail = F)
        pval_final1<-p.adjust(pval_final,method = "BH")
        corrected_pos<-index_nonzero[which(pval_final1<0.05)]
        return_list<-c(return_list,
                       list("corrected_pos"=corrected_pos,
                            "nonzero_pos"=index_nonzero,
                            "pval"=pval_final,
                            "nonzero_var"=final_v))
    }
    return(return_list)
}
