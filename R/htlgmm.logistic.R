#pacman::p_load(expm,magic,glmnet,cluster,MASS,locfit,corpcor,caret,pROC,mvtnorm)

U_func_binary_univ<-function(Z,W,A,y,beta,hat_thetaA,study_info){
    ZWA<-cbind(Z,W,A)
    expit_beta<-c(expit(ZWA%*%beta))
    u1<-crossprod((expit_beta-y),ZWA)
    u2_beta<-expit_beta%*%Z
    u2_theta<-sapply(1:ncol(Z), function(id){
        u2_id<-c(expit(cbind(A,Z[,id])%*%c(hat_thetaA,study_info[[id]]$Coeff)))
        c(u2_id%*%Z[,id])})
    u2<-u2_beta - u2_theta
    u<-c(u1,u2)*(1/nrow(Z))
    u
}

grad_U2_wrt_theta_func_binary_univ<-function(
        Z,A,y,hat_thetaA,study_info){
    u2_theta<-sapply(1:ncol(Z), function(id){
        dexpit_id<-c(Z[,id])*c(dexpit(cbind(A,Z[,id])%*%c(hat_thetaA,study_info[[id]]$Coeff)))
        u2id_thetaA_gradient<-dexpit_id%*%A
        u2id_thetaZ_gradient<-dexpit_id%*%Z[,id]
        c(u2id_thetaA_gradient,u2id_thetaZ_gradient)
    })
    #u2_theta is (1+len_GPC+1)*N_SNP
    # Gradient is N_equation*N_variables
    -(1/nrow(Z))*t(rbind(u2_theta[-nrow(u2_theta),],diag(u2_theta[nrow(u2_theta),])))
}

var_U_beta_theta_func_binary_univ<-function(
        Z,W,A,y,beta,hat_thetaA,study_info
){
    ZWA<-cbind(Z,W,A)
    expit_beta<-expit(ZWA%*%beta)
    var_11<-crossprod(ZWA*c(expit_beta-y))
    u2_theta<-sapply(1:ncol(Z), function(id){
        expit_id<-c(expit(cbind(A,Z[,id])%*%c(hat_thetaA,study_info[[id]]$Coeff)))
        expit_beta-expit_id
    }) #col is SNP #row is sample
    var_22<-crossprod(u2_theta*Z)
    var_12<-crossprod(ZWA*c(expit_beta-y),u2_theta*Z)
    (1/nrow(Z))*rbind(cbind(var_11,var_12),cbind(t(var_12),var_22))
}


var_theta_hat_vec_func_binary_univ<-function(study_info){
    var_vec<-sapply(1:length(study_info), function(i){
        study_info[[i]]$Covariance
    })
}

cov_U_with_theta_hat_func_binary_univ<-function(
        Z,W,A,y,beta,hat_thetaA,study_info){
    ZWA<-cbind(Z,W,A)
    expit_beta<-c(expit(ZWA%*%beta))
    expit_thetaA<-c(expit(A%*%hat_thetaA))
    cov_U1_theta_hat<-(1/nrow(Z))*crossprod(ZWA*c(expit_beta-y),A*c(expit_thetaA-y))
    u2_theta<-sapply(1:ncol(Z), function(id){
        expit_id<-c(expit(cbind(A,Z[,id])%*%c(hat_thetaA,study_info[[id]]$Coeff)))
        expit_beta-expit_id
    }) #col is SNP #row is sample
    cov_U2_theta_hat<-(1/nrow(Z))*crossprod(u2_theta*Z,A*c(expit_thetaA-y))
    rbind(cov_U1_theta_hat,cov_U2_theta_hat)
}

final_var_U_beta_theta_hat_func_binary_univ<-function(
        Z,W,A,y,beta,hat_thetaA,study_info){
    if(is.null(Z)){pZ<-0}else{pZ<-ncol(Z)}
    if(is.null(W)){pW<-0}else{pW<-ncol(W)}
    if(is.null(A)){pA<-0}else{pA<-ncol(A)}
    var_1st_U_beta_theta<-var_U_beta_theta_func_binary_univ(Z=Z,W=W,A=A,y=y,beta=beta,hat_thetaA=hat_thetaA,study_info=study_info)
    var_theta_vec<-var_theta_hat_vec_func_binary_univ(study_info = study_info)
    U_theta_gradient<-rbind(grad_U1_wrt_theta_func(pZ=pZ,pW=pW,pA=pA),
                            grad_U2_wrt_theta_func_binary_univ(Z=Z,A=A,y=y,hat_thetaA=hat_thetaA,study_info=study_info))
    var_thetaA<-var_thetaA_hat_func(A=A,y=y,hat_thetaA=hat_thetaA,study_info=study_info)
    var_grad_times_thetaA_hat<-U_theta_gradient[,1:ncol(A)]%*%var_thetaA%*%t(U_theta_gradient[,1:ncol(A)])
    var_grad_times_thetaZ_hat<-U_theta_gradient[,-c(1:ncol(A))]%*%(var_theta_vec*t(U_theta_gradient[,-c(1:ncol(A))]))*nrow(Z)
    var_2nd_grad_times_theta_hat<-var_grad_times_thetaA_hat+var_grad_times_thetaZ_hat
    mat_outside<-inv_grad_U3_wrt_thetaA_func(A=A,y=y,hat_thetaA=hat_thetaA)
    cov_U<-cov_U_with_theta_hat_func_binary_univ(Z=Z,W=W,A=A,y=y,beta=beta,hat_thetaA=hat_thetaA,study_info=study_info)
    cov_3rd_between_1st_2nd<-cov_U%*%mat_outside%*%t(U_theta_gradient[,1:ncol(A)])

    res<-(var_1st_U_beta_theta+var_2nd_grad_times_theta_hat+cov_3rd_between_1st_2nd+t(cov_3rd_between_1st_2nd))
    res
}

###### multivariate case

# theta1 = thetaA theta2=theta
# thetaA from A, theta from Z, from summary statistics
U_func_binary_multiv<-function(Z,W,A,y,beta,hat_thetaA,study_info){
    ZWA<-cbind(Z,W,A)
    expit_beta<-c(expit(ZWA%*%beta))
    u1<-crossprod((expit_beta-y),ZWA)
    u2<-c(expit_beta%*%Z-c(expit(cbind(A,Z)%*%c(hat_thetaA,study_info[[1]]$Coeff)))%*%Z)
    u<-c(u1,u2)*(1/nrow(Z))
    u
}
dexpit<-function(x){expit(x)*(1-expit(x))}

#grad_U_wrt_beta_func<-function(Z,W,A,y,beta){
#    ZWA<-cbind(Z,W,A)
#dexpit_beta<-dexpit(UKBB_pop[,-1]%*%beta)
#    dexpit_beta<-dexpit(ZWA%*%beta)
#U1_beta_gradient<-crossprod(UKBB_pop[,-1]*c(dexpit_beta),UKBB_pop[,-1])*(-1)
#    U1_beta_gradient<-crossprod(ZWA*c(dexpit_beta),ZWA)
#U2_beta_gradient<-crossprod(UKBB_pop[,var_SNP]*c(dexpit_beta),UKBB_pop[,-1])
#    U2_beta_gradient<-crossprod(Z*c(dexpit_beta),ZWA)
#rbind(U1_beta_gradient,U2_beta_gradient)*(1/N_Pop)
#    rbind(U1_beta_gradient,U2_beta_gradient)*(1/nrow(Z))}

#U3_func<-function(A,y,hat_thetaA){
#    expit_thetaA<-expit(A%*%hat_thetaA)
#    crossprod(c(expit_thetaA-y),A)}

### Generate coefficients of GPCs.

inv_grad_U3_wrt_thetaA_func<-function(A,y,hat_thetaA){
    dexpit_thetaA<-dexpit(A%*%hat_thetaA)
    mat<-(-1/nrow(A))*crossprod(A*c(dexpit_thetaA),A)
    -ginv(mat)
}

grad_U1_wrt_theta_func<-function(pZ,pW,pA){matrix(0,nrow=pZ+pW+pA,ncol=pZ+pA)}

grad_U2_wrt_theta_func_binary_multiv<-function(
        Z,A,y,hat_thetaA,study_info){
    dexpit_theta<-dexpit(cbind(Z,A)%*%c(study_info[[1]]$Coeff,hat_thetaA))
    U2_theta_gradient<-crossprod(Z*c(dexpit_theta),cbind(A,Z))
    -(1/nrow(Z))*U2_theta_gradient
}

var_U_beta_theta_func_binary_multiv<-function(
        Z,W,A,y,beta,hat_thetaA,study_info){
    ZWA<-cbind(Z,W,A)
    expit_beta<-c(expit(ZWA%*%beta))
    var_11<-crossprod(ZWA*c(expit_beta-y))
    u2_theta<-expit_beta-c(expit(cbind(Z,A)%*%c(study_info[[1]]$Coeff,hat_thetaA)))
    var_22<-crossprod(u2_theta*Z)
    var_12<-crossprod(ZWA*c(expit_beta-y),u2_theta*Z)
    (1/nrow(Z))*rbind(cbind(var_11,var_12),cbind(t(var_12),var_22))
}

var_thetaA_hat_func<-function(
        A,y,hat_thetaA,study_info){
    expit_thetaA<-expit(A%*%hat_thetaA)
    mat_inside<-(1/nrow(A))*crossprod(c(expit_thetaA-y)*A)
    mat_outside<-inv_grad_U3_wrt_thetaA_func(A=A,y=y,hat_thetaA=hat_thetaA)
    mat_outside%*%mat_inside%*%t(mat_outside)
}
var_theta_hat_vec_func_binary_multiv<-function(study_info){
    study_info[[1]]$Covariance
}

cov_U_with_theta_hat_func_binary_multiv<-function(
        Z,W,A,y,beta,hat_thetaA,study_info){
    ZWA<-cbind(Z,W,A)
    expit_beta<-c(expit(ZWA%*%beta))
    expit_thetaA<-c(expit(A%*%hat_thetaA))
    cov_U1_theta_hat<-(1/nrow(Z))*crossprod(ZWA*c(expit_beta-y),A*c(expit_thetaA-y))
    u2_theta<-expit_beta-c(expit(cbind(Z,A)%*%c(study_info[[1]]$Coeff,hat_thetaA)))
    cov_U2_theta_hat<-(1/nrow(Z))*crossprod(u2_theta*Z,A*c(expit_thetaA-y))
    rbind(cov_U1_theta_hat,cov_U2_theta_hat)
}

final_var_U_beta_theta_hat_func_binary_multiv<-function(
        Z,W,A,y,beta,hat_thetaA,study_info){
    if(is.null(Z)){pZ<-0}else{pZ<-ncol(Z)}
    if(is.null(W)){pW<-0}else{pW<-ncol(W)}
    if(is.null(A)){pA<-0}else{pA<-ncol(A)}
    var_1st_U_beta_theta<-var_U_beta_theta_func_binary_multiv(Z=Z,W=W,A=A,y=y,beta=beta,hat_thetaA=hat_thetaA,study_info=study_info)
    var_theta_mat<-var_theta_hat_vec_func_binary_multiv(study_info = study_info)
    U_theta_gradient<-rbind(grad_U1_wrt_theta_func(pZ=pZ,pW=pW,pA=pA),
                            grad_U2_wrt_theta_func_binary_multiv(Z=Z,A=A,y=y,hat_thetaA=hat_thetaA,study_info=study_info))
    var_thetaA<-var_thetaA_hat_func(A=A,y=y,hat_thetaA=hat_thetaA,study_info=study_info)

    var_grad_times_thetaA_hat<-U_theta_gradient[,1:ncol(A)]%*%var_thetaA%*%t(U_theta_gradient[,1:ncol(A)])
    var_grad_times_thetaZ_hat<-U_theta_gradient[,-c(1:ncol(A))]%*%var_theta_mat%*%t(U_theta_gradient[,-c(1:ncol(A))])*nrow(Z)
    var_2nd_grad_times_theta_hat<-var_grad_times_thetaA_hat+var_grad_times_thetaZ_hat
    mat_outside<-inv_grad_U3_wrt_thetaA_func(A=A,y=y,hat_thetaA=hat_thetaA)
    cov_U<-cov_U_with_theta_hat_func_binary_multiv(Z=Z,W=W,A=A,y=y,beta=beta,hat_thetaA=hat_thetaA,study_info=study_info)
    cov_3rd_between_1st_2nd<-cov_U%*%mat_outside%*%t(U_theta_gradient[,1:ncol(A)])
    res<-(var_1st_U_beta_theta+var_2nd_grad_times_theta_hat+cov_3rd_between_1st_2nd+t(cov_3rd_between_1st_2nd))
    res
}

############## Beta
pseudo_Xy_binary_multiv<-function(C_half,
                                  Z,W,A,y,beta,hat_thetaA,study_info){
    ZWA<-cbind(Z,W,A)
    expit_beta<-c(expit(ZWA%*%beta))
    dexpit_beta<-expit_beta*(1-expit_beta)
    pseudo_X<-C_half%*%rbind(t(ZWA),t(Z))%*%(ZWA*c(dexpit_beta))
    u<-c(C_half%*%U_func_binary_multiv(Z=Z,W=W,A=A,y=y,beta=beta,hat_thetaA=hat_thetaA,study_info=study_info))
    pseudo_y<- -u*nrow(Z) + c(pseudo_X%*%beta)
    newList<-list("pseudo_X"=pseudo_X,"pseudo_y"=pseudo_y)
    newList
}

pseudo_Xy_binary_univ<-function(C_half,
                                Z,W,A,y,beta,hat_thetaA,study_info){
    ZWA<-cbind(Z,W,A)
    expit_beta<-c(expit(ZWA%*%beta))
    dexpit_beta<-expit_beta*(1-expit_beta)
    pseudo_X<-C_half%*%rbind(t(ZWA),t(Z))%*%(ZWA*c(dexpit_beta))
    u<-c(C_half%*%U_func_binary_univ(Z=Z,W=W,A=A,y=y,beta=beta,hat_thetaA=hat_thetaA,study_info=study_info))
    pseudo_y<- -u*nrow(Z) + c(pseudo_X%*%beta)
    newList<-list("pseudo_X"=pseudo_X,"pseudo_y"=pseudo_y)
    newList
}

## cross validation helper function 1

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
        dev_lam_ratio_fold<-sapply(lambda_list,function(cur_lam){
            sapply(ratio_range,function(cur_ratio){
                ratio_vec<-c(rep(cur_ratio,pZ),rep(1,pW+pA))
                w_adaptive_ratio<-w_adaptive*ratio_vec
                cv_fit<-glmnet(x=pseudo_X_train,y=pseudo_y_train,
                               standardize=F,intercept=F,alpha = final_alpha,
                               penalty.factor = w_adaptive_ratio,
                               lambda = cur_lam)
                cur_beta<-coef.glmnet(cv_fit)[-1]
                probtest <- c(expit(cbind(Ztest,Wtest,Atest)%*%cur_beta))
                cur_dev <- -2*sum( ytest * log(probtest) + (1 - ytest) * log(1 - probtest) )
                #suppressMessages(cur_auc<-c(auc(ytest,c(expit(cbind(Ztest,Wtest,Atest)%*%cur_beta)),direction = "<")))
                cur_dev
            })
        }) # row is ratio_range & col is lambda_list
        dev_lam_ratio_fold
    })

    sum_dev_lam_ratio<-Reduce(`+`, dev_lam_ratio)
}


## cross validation helper function 2

cv_dev_lambda_func<-function(index_fold,Z,W,A,y,
                             C_half,beta_initial,hat_thetaA,
                             study_info,lambda_list,
                             w_adaptive,final_alpha,
                             pseudo_Xy){
    dev_fold<-sapply(1:length(index_fold), function(cur_fold){
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
            probtest <- c(expit(cbind(Ztest,Wtest,Atest)%*%cur_beta))
            cur_dev <- -2*sum( ytest * log(probtest) + (1 - ytest) * log(1 - probtest) )
            #suppressMessages(cur_auc<-c(auc(ytest,c(expit(cbind(Ztest,Wtest,Atest)%*%cur_beta)),direction = "<")))
            #sum((cbind(Ztest,Wtest,Atest)%*%cur_beta - ytest)^2)
            cur_dev
        })
        dev_lam
    })
    rowMeans(dev_fold)
}

## cross validation helper function 3

holdout_dev_lambda_ratio_func<-function(lambda_list,ratio_range,pZ,pW,pA,
                                        w_adaptive,pseudo_X_train,pseudo_y_train,
                                        final_alpha,Ztest,Wtest,Atest,ytest){
    dev_lam_ratio<-sapply(lambda_list,function(cur_lam){
        sapply(ratio_range,function(cur_ratio){
            ratio_vec<-c(rep(cur_ratio,pZ),rep(1,pW+pA))
            w_adaptive_ratio<-w_adaptive*ratio_vec
            cv_fit<-glmnet(x=pseudo_X_train,y=pseudo_y_train,
                           standardize=F,intercept=F,alpha = final_alpha,
                           penalty.factor = w_adaptive_ratio,
                           lambda = cur_lam)
            cur_beta<-coef.glmnet(cv_fit)[-1]
            probtest <- c(expit(cbind(Ztest,Wtest,Atest)%*%cur_beta))
            cur_dev <- -2*sum( ytest * log(probtest) + (1 - ytest) * log(1 - probtest) )
            #suppressMessages(cur_auc<-c(auc(ytest,c(expit(cbind(Ztest,Wtest,Atest)%*%cur_beta)),direction = "<")))
            cur_dev
        })
    }) # row is ratio_range & col is lambda_list
    sum_dev_lam_ratio<-Reduce(`+`, dev_lam_ratio)
}


## cross validation helper function 4

holdout_dev_lambda_func<-function(lambda_list,pseudo_X_train,pseudo_y_train,
                                  final_alpha,w_adaptive,Ztest,Wtest,Atest,ytest){
    dev_lam<-sapply(lambda_list,function(cur_lam){
        cv_fit<-glmnet(x= (pseudo_X_train),y= (pseudo_y_train),
                       standardize=F,intercept=F,
                       alpha = final_alpha,penalty.factor = w_adaptive,
                       lambda = cur_lam)
        cur_beta<-coef.glmnet(cv_fit)[-1]
        probtest <- c(expit(cbind(Ztest,Wtest,Atest)%*%cur_beta))
        cur_dev <- -2*sum( ytest * log(probtest) + (1 - ytest) * log(1 - probtest) )
        #suppressMessages(cur_auc<-c(auc(ytest,c(expit(cbind(Ztest,Wtest,Atest)%*%cur_beta)),direction = "<")))
        cur_dev
    })
    dev_lam
}

#' @import glmnet
#' @import stats
#' @importFrom caret createFolds createDataPartition
#' @importFrom locfit expit
#' @importFrom corpcor fast.svd
#' @importFrom magic adiag
#' @importFrom MASS ginv
#'
htlgmm.binary<-function(
        y,Z,W=NULL,A=1,
        study_info=NULL,
        summary_type = "multi",
        penalty_type = "lasso",
        initial_with_type = "ridge",
        beta_initial = NULL,
        remove_penalty_Z = FALSE,
        remove_penalty_W = FALSE,
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

    if(!is.null(fix_ratio)){
        if(tune_ratio){
            stop("If ratio is fixed, please set tune_ratio as FALSE")
        }else if(remove_penalty_Z | remove_penalty_W){
            stop("If ratio is fixed, please set remove_penalty's as FALSE")
        }
    }

    nZ<-nrow(Z)
    if(A==1){A<-matrix(1,nrow=nZ,ncol=1)}
    nZext<-study_info[[1]]$Sample_size
    pZ<-ncol(Z)
    if(is.null(W)){pW<-0}else{pW<-ncol(W)}
    if(is.null(A)){pA<-0}else{pA<-ncol(A)}
    if(nZ<2*pZ+pW+pA){use_sparseC<-TRUE}
    if (summary_type == "uni"){
        pseudo_Xy <- pseudo_Xy_binary_univ
        inv_C_func <- final_var_U_beta_theta_hat_func_binary_univ

        hat_thetaA_glm<-glm(y~0+.,data = data.frame(y,A),family = "binomial")
        hat_thetaA<-hat_thetaA_glm$coefficients[1:pA]
    }else if (summary_type == "multi"){
        pseudo_Xy <- pseudo_Xy_binary_multiv
        inv_C_func <- final_var_U_beta_theta_hat_func_binary_multiv

        hat_thetaA_glm<-glm(y~0+.,data = data.frame(y,A,Z),family = "binomial")
        hat_thetaA<-hat_thetaA_glm$coefficients[1:pA]
    }
    ZWA<-cbind(Z,W,A)
    Zid<-1:pZ
    Wid<-(pZ+1):(pZ+pW)
    Aid<-(pZ+pW+1):(pZ+pW+pA)
    fix_penalty<-c(rep(1,pZ+pW),rep(0,pA))
    if(remove_penalty_Z){fix_penalty[Zid]<-0}
    if(remove_penalty_W){fix_penalty[Wid]<-0}
    if(!is.null(fix_ratio)){fix_penalty[Zid]<-fix_ratio}
    if(!is.null(beta_initial)){
        if(length(beta_initial)!=pZ+pW+pA){
            warning("beta_initial should be from Z,W,A.\n Length not match, compute default initial instead.")
            beta_initial<-NULL
        }
    }
    if(is.null(beta_initial) & initial_with_type %in% c("ridge","lasso")){
        if(initial_with_type == "ridge"){initial_alpha=0}else{initial_alpha=1}
        if(pA == 1 & unique(A) == 1){
            fit_initial<-cv.glmnet(x=cbind(Z,W),y=y,alpha = initial_alpha,penalty.factor = fix_penalty[c(Zid,Wid)],family="binomial")
            beta_initial<-c(coef.glmnet(fit_initial,s="lambda.min")[-1])
            beta_initial<-c(beta_initial,coef.glmnet(fit_initial,s="lambda.min")[1])
        }else if(pA > 1 & unique(A[,1])==1){
            fit_initial<-cv.glmnet(x=cbind(Z,W,A[,-1]),y=y,alpha = initial_alpha,penalty.factor = fix_penalty[c(Zid,Wid,Aid[-1])],family="binomial")
            beta_initial<-c(coef.glmnet(fit_initial,s="lambda.min")[-1])
            beta_initial<-c(beta_initial[c(Zid,Wid)],coef.glmnet(fit_initial,s="lambda.min")[1],beta_initial[c(Aid[-1])])
        }else{
            stop("The first column of A should be 1 for intercept.")
        }
    }else{stop("Select Initial Type from c('ridge','lasso')")}
    if (penalty_type == "adaptivelasso"){
        w_adaptive<-1/abs(beta_initial)^gamma_adaptivelasso
        w_adaptive[is.infinite(w_adaptive)]<-max(w_adaptive[!is.infinite(w_adaptive)])*100
        w_adaptive<-w_adaptive*fix_penalty
    }else{w_adaptive<-fix_penalty}
    # Estimation of C
    inv_C = inv_C_func(Z=Z,W=W,A=A,y=y,
                       beta=beta_initial,hat_thetaA=hat_thetaA,
                       study_info=study_info)
    if(use_sparseC){
        C_half<-diag(1/sqrt(diag(inv_C)))
    }else{
        inv_C_svd<-fast.svd(inv_C+diag(1e-15,nrow(inv_C)))
        C_half<-inv_C_svd$v%*%diag(1/sqrt(inv_C_svd$d))%*%t(inv_C_svd$u)
    }

    # Prepare for final model
    pseudo_Xy_list<-pseudo_Xy(C_half,Z,W,A,y,beta = beta_initial,hat_thetaA = hat_thetaA,study_info=study_info)
    initial_sf<-nZ/sqrt(nrow(pseudo_Xy_list$pseudo_X))
    pseudo_X<-pseudo_Xy_list$pseudo_X/initial_sf
    pseudo_y<-pseudo_Xy_list$pseudo_y/initial_sf

    # Fit final model
    fit_final<-cv.glmnet(x= pseudo_X,y= pseudo_y,standardize=F,
                         intercept=F,alpha = final_alpha,penalty.factor = w_adaptive)
    if(is.null(lambda_list)){
        lambda_list<-fit_final$lambda
        lambda_list<-lambda_list[!is.na(lambda_list)]
    }
    if(!is.null(fix_lambda)){
        validation_type<-"None"
        if(fix_lambda<0){stop("The fixed lambda should be nonnegative.")}
    }

    if(tune_ratio & !remove_penalty_Z & !remove_penalty_W){
        if(is.null(ratio_range)){
            if(is.null(ratio_lower)){ratio_lower<-sqrt(nZ/(nZ+nZext))/2}
            if(is.null(ratio_upper)){ratio_upper<-(nZ)^(1/3)/2}
            #ratio_range<-exp(seq(log(ratio_lower),log(ratio_upper),(log(ratio_upper)-log(ratio_lower))/ratio_count))
            ratio_range<-(seq(sqrt(ratio_lower),sqrt(ratio_upper),(sqrt(ratio_upper)-sqrt(ratio_lower))/ratio_count)^2)
            ratio_range<-c(1,ratio_range)
        }
    }else{tune_ratio<-FALSE}
    if(validation_type == "None"){
        fit_final_fixed_lambda<-glmnet(x= pseudo_X,y= pseudo_y,standardize=F,
                                       intercept=F,alpha = final_alpha,
                                       penalty.factor = w_adaptive,
                                       lambda = fix_lambda)
        beta<-coef.glmnet(fit_final_fixed_lambda)[-1]
        return_list<-list("beta"=beta,
                          "fix_lambda"=fix_lambda,
                          "fix_ratio"=fix_ratio)
    }else if(validation_type == 'holdout'){
        index_valid<-createDatapartition(y,p = holdout_p)$Resample1
        Ztrain<-Z[-index_valid,]
        Ztest<-Z[index_valid,]
        if(!is.null(W)){
            Wtrain<-W[-index_valid,]
            Wtest<-W[index_valid,]
        }else{
            Wtrain<-NULL
            Wtest<-NULL}
        if(!is.null(A)){
            Atrain<-A[-index_valid,]
            Atest<-A[index_valid,]
        }else{
            Atrain<-NULL
            Atest<-NULL}
        ytrain<-y[-index_valid]
        ytest<-y[index_valid]

        pseudo_Xy_list_train<-pseudo_Xy(C_half,Ztrain,Wtrain,Atrain,
                                        ytrain,beta = beta_initial,hat_thetaA = hat_thetaA,
                                        study_info=study_info)
        initial_sf_train<-nrow(Ztrain)/sqrt(nrow(pseudo_Xy_list_train$pseudo_X))
        pseudo_X_train<-pseudo_Xy_list_train$pseudo_X/initial_sf_train
        pseudo_y_train<-pseudo_Xy_list_train$pseudo_y/initial_sf_train


        if(tune_ratio){
            holdout_dev<-holdout_dev_lambda_ratio_func(lambda_list,ratio_range,pZ,pW,pA,
                                                       w_adaptive,pseudo_X_train,pseudo_y_train,
                                                       final_alpha,Ztest,Wtest,Atest,ytest)
            ids<-which(holdout_dev==min(holdout_dev),arr.ind = TRUE)
            final.ratio.min<-ratio_range[ids[1]]
            final.lambda.min<-lambda_list[ids[2]]
        }else{
            holdout_dev<-holdout_dev_lambda_func(lambda_list,pseudo_X_train,pseudo_y_train,
                                                 final_alpha,w_adaptive,Ztest,Wtest,Atest,ytest)
            final.lambda.min<-lambda_list[which.min(holdout_dev)]
            final.ratio.min<-1
        }
        ratio_vec<-c(rep(final.ratio.min,pZ),rep(1,pW+pA))
        w_adaptive_ratio<-w_adaptive*ratio_vec
        fit_final_lam_ratio<-glmnet(x= pseudo_X_train,y= pseudo_y_train,standardize=F,
                                    intercept=F,alpha = final_alpha,
                                    penalty.factor = w_adaptive_ratio,
                                    lambda = final.lambda.min)
        beta<-coef.glmnet(fit_final_lam_ratio)[-1]
        return_list<-list("beta"=beta,#[-c((length(beta)-pA+1):length(beta)) ],
                          "lambda_list"=lambda_list,
                          "ratio_list"=ratio_range,
                          "holdout_dev"= holdout_dev,
                          "lambda_min"=final.lambda.min,
                          "ratio_min"=final.ratio.min)

    }else if(validation_type == 'cv'){
        index_fold<-createFolds(y,k = nfolds)

        if(tune_ratio){
            cv_dev<-cv_dev_lambda_ratio_func(index_fold,Z,W,A,y,
                                             C_half,beta_initial,hat_thetaA,
                                             study_info,lambda_list,
                                             ratio_range,pZ,pW,pA,
                                             w_adaptive,final_alpha,pseudo_Xy)
            ids<-which(cv_dev==min(cv_dev),arr.ind = TRUE)
            final.ratio.min<-ratio_range[ids[1]]
            final.lambda.min<-lambda_list[ids[2]]
        }else{
            cv_dev<-cv_dev_lambda_func(index_fold,Z,W,A,y,
                                       C_half,beta_initial,hat_thetaA,
                                       study_info,lambda_list,
                                       w_adaptive,final_alpha,pseudo_Xy)
            final.lambda.min<-lambda_list[which.min(cv_dev)]
            final.ratio.min<-1
        }
        ratio_vec<-c(rep(final.ratio.min,pZ),rep(1,pW+pA))
        w_adaptive_ratio<-w_adaptive*ratio_vec
        fit_final_lam_ratio<-glmnet(x= pseudo_X,y= pseudo_y,standardize=F,
                                    intercept=F,alpha = final_alpha,
                                    penalty.factor = w_adaptive_ratio,
                                    lambda = final.lambda.min)

        beta<-coef.glmnet(fit_final_lam_ratio)[-1]
        return_list<-list("beta"=beta,#[-c((length(beta)-pA+1):length(beta)) ],
                          "lambda_list"=lambda_list,
                          "ratio_list"=ratio_range,
                          "cv_dev"=cv_dev,
                          "lambda_min"=final.lambda.min,
                          "ratio_min"=final.ratio.min)
    }

    index_nonzero<-which(beta!=0)
    if(inference & length(index_nonzero) > 1){
        # refine C

        inv_C = inv_C_func(Z=Z,W=W,A=A,y=y,
                           beta=beta,hat_thetaA=hat_thetaA,
                           study_info=study_info)
        inv_C_svd<-fast.svd(inv_C+diag(1e-15,nrow(inv_C)))
        C_half<-inv_C_svd$v%*%diag(1/sqrt(inv_C_svd$d))%*%t(inv_C_svd$u)
        ZWA<-cbind(Z,W,A)
        expit_beta<-c(expit(ZWA%*%beta))
        dexpit_beta<-expit_beta*(1-expit_beta)
        pseudo_X<-C_half%*%rbind(t(ZWA),t(Z))%*%(ZWA*c(dexpit_beta))/nZ
        #Sigsum_half<-cbind(ZWAtZWA/nZ,crossprod(ZWA,Z)/nZ)%*%C_half
        #Sigsum_scaled<-Sigsum_half%*%t(Sigsum_half)
        Sigsum_scaled<-crossprod(pseudo_X)
        Sigsum_scaled_nonzero<-Sigsum_scaled[index_nonzero,index_nonzero]
        inv_Sigsum_scaled_nonzero<-solve(Sigsum_scaled_nonzero)
        final_v<-diag(inv_Sigsum_scaled_nonzero)

        pval_final<-pchisq(nZ*beta[index_nonzero]^2/final_v,1,lower.tail = F)
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
