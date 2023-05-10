#pacman::p_load(expm,magic,glmnet,cluster,MASS,locfit,corpcor,caret,pROC,mvtnorm)
# theta1 = thetaG theta2=theta
# thetaG from G, theta from X, from summary statistics
U_func_binary_univ<-function(X,A,G,y,beta,hat_thetaG,study_info){
    XAG<-cbind(X,A,G)
    expit_beta<-c(expit(XAG%*%beta))
    u1<-crossprod((expit_beta-y),XAG)
    u2_beta<-expit_beta%*%X
    u2_theta<-sapply(1:ncol(X), function(id){
        u2_id<-c(expit(cbind(G,X[,id])%*%c(hat_thetaG,study_info[[id]]$Coeff)))
        c(u2_id%*%X[,id])})
    u2<-u2_beta - u2_theta
    u<-c(u1,u2)*(1/nrow(X))
    u
}

grad_U2_wrt_theta_func_binary_univ<-function(
        X,G,y,hat_thetaG,study_info){
    u2_theta<-sapply(1:ncol(X), function(id){
        dexpit_id<-c(X[,id])*c(dexpit(cbind(G,X[,id])%*%c(hat_thetaG,study_info[[id]]$Coeff)))
        u2id_thetaG_gradient<-dexpit_id%*%G
        u2id_thetaX_gradient<-dexpit_id%*%X[,id]
        c(u2id_thetaG_gradient,u2id_thetaX_gradient)
    })
    #u2_theta is (1+len_GPC+1)*N_SNP
    # Gradient is N_equation*N_variables
    -(1/nrow(X))*t(rbind(u2_theta[-nrow(u2_theta),],diag(u2_theta[nrow(u2_theta),])))
}

var_U_beta_theta_func_binary_univ<-function(
        X,A,G,y,beta,hat_thetaG,study_info
){
    XAG<-cbind(X,A,G)
    expit_beta<-expit(XAG%*%beta)
    var_11<-crossprod(XAG*c(expit_beta-y))
    u2_theta<-sapply(1:ncol(X), function(id){
        expit_id<-c(expit(cbind(G,X[,id])%*%c(hat_thetaG,study_info[[id]]$Coeff)))
        expit_beta-expit_id
    }) #col is SNP #row is sample
    var_22<-crossprod(u2_theta*X)
    var_12<-crossprod(XAG*c(expit_beta-y),u2_theta*X)
    (1/nrow(X))*rbind(cbind(var_11,var_12),cbind(t(var_12),var_22))
}


var_theta_hat_vec_func_binary_univ<-function(study_info){
    var_vec<-sapply(1:length(study_info), function(i){
        study_info[[i]]$Covariance
    })
}

cov_U_with_theta_hat_func_binary_univ<-function(
        X,A,G,y,beta,hat_thetaG,study_info){
    XAG<-cbind(X,A,G)
    expit_beta<-c(expit(XAG%*%beta))
    expit_thetaG<-c(expit(G%*%hat_thetaG))
    cov_U1_theta_hat<-(1/nrow(X))*crossprod(XAG*c(expit_beta-y),G*c(expit_thetaG-y))
    u2_theta<-sapply(1:ncol(X), function(id){
        expit_id<-c(expit(cbind(G,X[,id])%*%c(hat_thetaG,study_info[[id]]$Coeff)))
        expit_beta-expit_id
    }) #col is SNP #row is sample
    cov_U2_theta_hat<-(1/nrow(X))*crossprod(u2_theta*X,G*c(expit_thetaG-y))
    rbind(cov_U1_theta_hat,cov_U2_theta_hat)
}

final_var_U_beta_theta_hat_func_binary_univ<-function(
        X,A,G,y,beta,hat_thetaG,study_info){
    if(is.null(X)){pX<-0}else{pX<-ncol(X)}
    if(is.null(A)){pA<-0}else{pA<-ncol(A)}
    if(is.null(G)){pG<-0}else{pG<-ncol(G)}
    var_1st_U_beta_theta<-var_U_beta_theta_func_binary_univ(X=X,A=A,G=G,y=y,beta=beta,hat_thetaG=hat_thetaG,study_info=study_info)
    var_theta_vec<-var_theta_hat_vec_func_binary_univ(study_info = study_info)
    U_theta_gradient<-rbind(grad_U1_wrt_theta_func(pX=pX,pA=pA,pG=pG),
                            grad_U2_wrt_theta_func_binary_univ(X=X,G=G,y=y,hat_thetaG=hat_thetaG,study_info=study_info))
    var_thetaG<-var_thetaG_hat_func(G=G,y=y,hat_thetaG=hat_thetaG,study_info=study_info)
    var_grad_times_thetaG_hat<-U_theta_gradient[,1:ncol(G)]%*%var_thetaG%*%t(U_theta_gradient[,1:ncol(G)])
    var_grad_times_thetaX_hat<-U_theta_gradient[,-c(1:ncol(G))]%*%(var_theta_vec*t(U_theta_gradient[,-c(1:ncol(G))]))*nrow(X)
    var_2nd_grad_times_theta_hat<-var_grad_times_thetaG_hat+var_grad_times_thetaX_hat
    mat_outside<-inv_grad_U3_wrt_thetaG_func(G=G,y=y,hat_thetaG=hat_thetaG)
    cov_U<-cov_U_with_theta_hat_func_binary_univ(X=X,A=A,G=G,y=y,beta=beta,hat_thetaG=hat_thetaG,study_info=study_info)
    cov_3rd_between_1st_2nd<-cov_U%*%mat_outside%*%t(U_theta_gradient[,1:ncol(G)])

    res<-(var_1st_U_beta_theta+var_2nd_grad_times_theta_hat+cov_3rd_between_1st_2nd+t(cov_3rd_between_1st_2nd))
    res
}

###### multivariate case

# theta1 = thetaG theta2=theta
# thetaG from G, theta from X, from summary statistics
U_func_binary_multiv<-function(X,A,G,y,beta,hat_thetaG,study_info){
    XAG<-cbind(X,A,G)
    expit_beta<-c(expit(XAG%*%beta))
    u1<-crossprod((expit_beta-y),XAG)
    u2<-c(expit_beta%*%X-c(expit(cbind(G,X)%*%c(hat_thetaG,study_info[[1]]$Coeff)))%*%X)
    u<-c(u1,u2)*(1/nrow(X))
    u
}
dexpit<-function(x){expit(x)*(1-expit(x))}

#grad_U_wrt_beta_func<-function(X,A,G,y,beta){
#    XAG<-cbind(X,A,G)
#dexpit_beta<-dexpit(UKBB_pop[,-1]%*%beta)
#    dexpit_beta<-dexpit(XAG%*%beta)
#U1_beta_gradient<-crossprod(UKBB_pop[,-1]*c(dexpit_beta),UKBB_pop[,-1])*(-1)
#    U1_beta_gradient<-crossprod(XAG*c(dexpit_beta),XAG)
#U2_beta_gradient<-crossprod(UKBB_pop[,var_SNP]*c(dexpit_beta),UKBB_pop[,-1])
#    U2_beta_gradient<-crossprod(X*c(dexpit_beta),XAG)
#rbind(U1_beta_gradient,U2_beta_gradient)*(1/N_Pop)
#    rbind(U1_beta_gradient,U2_beta_gradient)*(1/nrow(X))}

#U3_func<-function(G,y,hat_thetaG){
#    expit_thetaG<-expit(G%*%hat_thetaG)
#    crossprod(c(expit_thetaG-y),G)}

### Generate coefficients of GPCs.

inv_grad_U3_wrt_thetaG_func<-function(G,y,hat_thetaG){
    dexpit_thetaG<-dexpit(G%*%hat_thetaG)
    mat<-(-1/nrow(G))*crossprod(G*c(dexpit_thetaG),G)
    -ginv(mat)
}

grad_U1_wrt_theta_func<-function(pX,pA,pG){matrix(0,nrow=pX+pA+pG,ncol=pX+pG)}

grad_U2_wrt_theta_func_binary_multiv<-function(
        X,G,y,hat_thetaG,study_info){
    dexpit_theta<-dexpit(cbind(X,G)%*%c(study_info[[1]]$Coeff,hat_thetaG))
    U2_theta_gradient<-crossprod(X*c(dexpit_theta),cbind(G,X))
    -(1/nrow(X))*U2_theta_gradient
}

var_U_beta_theta_func_binary_multiv<-function(
        X,A,G,y,beta,hat_thetaG,study_info){
    XAG<-cbind(X,A,G)
    expit_beta<-c(expit(XAG%*%beta))
    var_11<-crossprod(XAG*c(expit_beta-y))
    u2_theta<-expit_beta-c(expit(cbind(X,G)%*%c(study_info[[1]]$Coeff,hat_thetaG)))
    var_22<-crossprod(u2_theta*X)
    var_12<-crossprod(XAG*c(expit_beta-y),u2_theta*X)
    (1/nrow(X))*rbind(cbind(var_11,var_12),cbind(t(var_12),var_22))
}

var_thetaG_hat_func<-function(
        G,y,hat_thetaG,study_info){
    expit_thetaG<-expit(G%*%hat_thetaG)
    mat_inside<-(1/nrow(G))*crossprod(c(expit_thetaG-y)*G)
    mat_outside<-inv_grad_U3_wrt_thetaG_func(G=G,y=y,hat_thetaG=hat_thetaG)
    mat_outside%*%mat_inside%*%t(mat_outside)
}
var_theta_hat_vec_func_binary_multiv<-function(study_info){
    study_info[[1]]$Covariance
}

cov_U_with_theta_hat_func_binary_multiv<-function(
        X,A,G,y,beta,hat_thetaG,study_info){
    XAG<-cbind(X,A,G)
    expit_beta<-c(expit(XAG%*%beta))
    expit_thetaG<-c(expit(G%*%hat_thetaG))
    cov_U1_theta_hat<-(1/nrow(X))*crossprod(XAG*c(expit_beta-y),G*c(expit_thetaG-y))
    u2_theta<-expit_beta-c(expit(cbind(X,G)%*%c(study_info[[1]]$Coeff,hat_thetaG)))
    cov_U2_theta_hat<-(1/nrow(X))*crossprod(u2_theta*X,G*c(expit_thetaG-y))
    rbind(cov_U1_theta_hat,cov_U2_theta_hat)
}

final_var_U_beta_theta_hat_func_binary_multiv<-function(
        X,A,G,y,beta,hat_thetaG,study_info){
    if(is.null(X)){pX<-0}else{pX<-ncol(X)}
    if(is.null(A)){pA<-0}else{pA<-ncol(A)}
    if(is.null(G)){pG<-0}else{pG<-ncol(G)}
    var_1st_U_beta_theta<-var_U_beta_theta_func_binary_multiv(X=X,A=A,G=G,y=y,beta=beta,hat_thetaG=hat_thetaG,study_info=study_info)
    var_theta_mat<-var_theta_hat_vec_func_binary_multiv(study_info = study_info)
    U_theta_gradient<-rbind(grad_U1_wrt_theta_func(pX=pX,pA=pA,pG=pG),
                            grad_U2_wrt_theta_func_binary_multiv(X=X,G=G,y=y,hat_thetaG=hat_thetaG,study_info=study_info))
    var_thetaG<-var_thetaG_hat_func(G=G,y=y,hat_thetaG=hat_thetaG,study_info=study_info)

    var_grad_times_thetaG_hat<-U_theta_gradient[,1:ncol(G)]%*%var_thetaG%*%t(U_theta_gradient[,1:ncol(G)])
    var_grad_times_thetaX_hat<-U_theta_gradient[,-c(1:ncol(G))]%*%var_theta_mat%*%t(U_theta_gradient[,-c(1:ncol(G))])*nrow(X)
    var_2nd_grad_times_theta_hat<-var_grad_times_thetaG_hat+var_grad_times_thetaX_hat
    mat_outside<-inv_grad_U3_wrt_thetaG_func(G=G,y=y,hat_thetaG=hat_thetaG)
    cov_U<-cov_U_with_theta_hat_func_binary_multiv(X=X,A=A,G=G,y=y,beta=beta,hat_thetaG=hat_thetaG,study_info=study_info)
    cov_3rd_between_1st_2nd<-cov_U%*%mat_outside%*%t(U_theta_gradient[,1:ncol(G)])
    res<-(var_1st_U_beta_theta+var_2nd_grad_times_theta_hat+cov_3rd_between_1st_2nd+t(cov_3rd_between_1st_2nd))
    res
}

############## Beta
pseudo_Xy_binary_multiv<-function(C_half,
        X,A,G,y,beta,hat_thetaG,study_info){
    XAG<-cbind(X,A,G)
    expit_beta<-c(expit(XAG%*%beta))
    dexpit_beta<-expit_beta*(1-expit_beta)
    pseudo_X<-C_half%*%rbind(t(XAG),t(X))%*%(XAG*c(dexpit_beta))
    u<-c(C_half%*%U_func_binary_multiv(X=X,A=A,G=G,y=y,beta=beta,hat_thetaG=hat_thetaG,study_info=study_info))
    pseudo_y<- -u*nrow(X) + c(pseudo_X%*%beta)
    newList<-list("pseudo_X"=pseudo_X,"pseudo_y"=pseudo_y)
    newList
}

pseudo_Xy_binary_univ<-function(C_half,
        X,A,G,y,beta,hat_thetaG,study_info){
    XAG<-cbind(X,A,G)
    expit_beta<-c(expit(XAG%*%beta))
    dexpit_beta<-expit_beta*(1-expit_beta)
    pseudo_X<-C_half%*%rbind(t(XAG),t(X))%*%(XAG*c(dexpit_beta))
    u<-c(C_half%*%U_func_binary_univ(X=X,A=A,G=G,y=y,beta=beta,hat_thetaG=hat_thetaG,study_info=study_info))
    pseudo_y<- -u*nrow(X) + c(pseudo_X%*%beta)
    newList<-list("pseudo_X"=pseudo_X,"pseudo_y"=pseudo_y)
    newList
}

## cross validation helper function 1

cv_auc_lambda_ratio_func<-function(index_fold,X,A,G,y,
                                   C_half,beta_initial,hat_thetaG,
                                   study_info,lambda_list,
                                   ratio_range,pX,pA,pG,
                                   w_adaptive,final_alpha,
                                   pseudo_Xy){
    auc_lam_ratio<-lapply(1:length(index_fold), function(cur_fold){
        index_test<-index_fold[[cur_fold]]
        Xtrain<-X[-index_test,]
        Xtest<-X[index_test,]
        if(!is.null(A)){
            Atrain<-A[-index_test,]
            Atest<-A[index_test,]
        }else{
            Atrain<-NULL
            Atest<-NULL}
        if(!is.null(G)){
            Gtrain<-G[-index_test,]
            Gtest<-G[index_test,]
        }else{
            Gtrain<-NULL
            Gtest<-NULL}
        ytrain<-y[-index_test]
        ytest<-y[index_test]
        pseudo_Xy_list_train<-pseudo_Xy(C_half,Xtrain,Atrain,Gtrain,
                                        ytrain,beta = beta_initial,hat_thetaG = hat_thetaG,
                                        study_info=study_info)
        initial_sf_train<-nrow(Xtrain)/sqrt(nrow(pseudo_Xy_list_train$pseudo_X))
        pseudo_X_train<-pseudo_Xy_list_train$pseudo_X/initial_sf_train
        pseudo_y_train<-pseudo_Xy_list_train$pseudo_y/initial_sf_train
        auc_lam_ratio_fold<-sapply(lambda_list,function(cur_lam){
            sapply(ratio_range,function(cur_ratio){
                ratio_vec<-c(rep(cur_ratio,pX),rep(1,pA+pG))
                w_adaptive_ratio<-w_adaptive*ratio_vec
                cv_fit<-glmnet(x=pseudo_X_train,y=pseudo_y_train,
                               standardize=F,intercept=F,alpha = final_alpha,
                               penalty.factor = w_adaptive_ratio,
                               lambda = cur_lam)
                cur_beta<-coef.glmnet(cv_fit)[-1]
                suppressMessages(cur_auc<-c(auc(ytest,c(expit(cbind(Xtest,Atest,Gtest)%*%cur_beta)),direction = "<")))
                cur_auc
            })
        }) # row is ratio_range & col is lambda_list
        auc_lam_ratio_fold
    })

    sum_auc_lam_ratio<-Reduce(`+`, auc_lam_ratio)
}


## cross validation helper function 2

cv_auc_lambda_func<-function(index_fold,X,A,G,y,
                             C_half,beta_initial,hat_thetaG,
                             study_info,lambda_list,
                             w_adaptive,final_alpha,
                             pseudo_Xy){
    auc_fold<-sapply(1:length(index_fold), function(cur_fold){
        index_test<-index_fold[[cur_fold]]
        Xtrain<-X[-index_test,]
        Xtest<-X[index_test,]
        if(!is.null(A)){
            Atrain<-A[-index_test,]
            Atest<-A[index_test,]
        }else{
            Atrain<-NULL
            Atest<-NULL}
        if(!is.null(G)){
            Gtrain<-G[-index_test,]
            Gtest<-G[index_test,]
        }else{
            Gtrain<-NULL
            Gtest<-NULL}
        ytrain<-y[-index_test]
        ytest<-y[index_test]
        pseudo_Xy_list_train<-pseudo_Xy(C_half,Xtrain,Atrain,Gtrain,
                                        ytrain,beta = beta_initial,hat_thetaG = hat_thetaG,
                                        study_info=study_info)
        initial_sf_train<-nrow(Xtrain)/sqrt(nrow(pseudo_Xy_list_train$pseudo_X))
        pseudo_X_train<-pseudo_Xy_list_train$pseudo_X/initial_sf_train
        pseudo_y_train<-pseudo_Xy_list_train$pseudo_y/initial_sf_train
        auc_lam<-sapply(lambda_list,function(cur_lam){
            cv_fit<-glmnet(x= (pseudo_X_train),y= (pseudo_y_train),
                           standardize=F,intercept=F,
                           alpha = final_alpha,penalty.factor = w_adaptive,lambda = cur_lam)
            cur_beta<-coef.glmnet(cv_fit)[-1]
            suppressMessages(cur_auc<-c(auc(ytest,c(expit(cbind(Xtest,Atest,Gtest)%*%cur_beta)),direction = "<")))
            #sum((cbind(Xtest,Atest,Gtest)%*%cur_beta - ytest)^2)
            cur_auc
        })
        auc_lam
    })
    rowMeans(auc_fold)
}

## cross validation helper function 3

holdout_auc_lambda_ratio_func<-function(lambda_list,ratio_range,pX,pA,pG,
                                        w_adaptive,pseudo_X_train,pseudo_y_train,
                                        final_alpha,Xtest,Atest,Gtest,ytest){
    auc_lam_ratio<-sapply(lambda_list,function(cur_lam){
        sapply(ratio_range,function(cur_ratio){
            ratio_vec<-c(rep(cur_ratio,pX),rep(1,pA+pG))
            w_adaptive_ratio<-w_adaptive*ratio_vec
            cv_fit<-glmnet(x=pseudo_X_train,y=pseudo_y_train,
                           standardize=F,intercept=F,alpha = final_alpha,
                           penalty.factor = w_adaptive_ratio,
                           lambda = cur_lam)
            cur_beta<-coef.glmnet(cv_fit)[-1]
            suppressMessages(cur_auc<-c(auc(ytest,c(expit(cbind(Xtest,Atest,Gtest)%*%cur_beta)),direction = "<")))
            cur_auc
        })
    }) # row is ratio_range & col is lambda_list
    sum_auc_lam_ratio<-Reduce(`+`, auc_lam_ratio)
}


## cross validation helper function 4

holdout_auc_lambda_func<-function(lambda_list,pseudo_X_train,pseudo_y_train,
                                  final_alpha,w_adaptive,Xtest,Atest,Gtest,ytest){
    auc_lam<-sapply(lambda_list,function(cur_lam){
        cv_fit<-glmnet(x= (pseudo_X_train),y= (pseudo_y_train),
                       standardize=F,intercept=F,
                       alpha = final_alpha,penalty.factor = w_adaptive,
                       lambda = cur_lam)
        cur_beta<-coef.glmnet(cv_fit)[-1]
        suppressMessages(cur_auc<-c(auc(ytest,c(expit(cbind(Xtest,Atest,Gtest)%*%cur_beta)),direction = "<")))
        cur_auc
    })
    auc_lam
}

intgmm.binary<-function(
        y,X,A=NULL,G=1,
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
        }else if(remove_penalty_X | remove_penalty_A){
            stop("If ratio is fixed, please set remove_penalty's as FALSE")
        }
    }

    nX<-nrow(X)
    if(G==1){G<-matrix(1,nrow=nX,ncol=1)}
    nXext<-study_info[[1]]$Sample_size
    pX<-ncol(X)
    if(is.null(A)){pA<-0}else{pA<-ncol(A)}
    if(is.null(G)){pG<-0}else{pG<-ncol(G)}
    if(nX<2*pX+pA+pG){use_sparseC<-TRUE}
    if (summary_type == "uni"){
        pseudo_Xy <- pseudo_Xy_binary_univ
        inv_C_func <- final_var_U_beta_theta_hat_func_binary_univ

        hat_thetaG_glm<-glm(y~0+.,data = data.frame(y,G),family = "binomial")
        hat_thetaG<-hat_thetaG_glm$coefficients[1:pG]
    }else if (summary_type == "multi"){
        pseudo_Xy <- pseudo_Xy_binary_multiv
        inv_C_func <- final_var_U_beta_theta_hat_func_binary_multiv

        hat_thetaG_glm<-glm(y~0+.,data = data.frame(y,G,X),family = "binomial")
        hat_thetaG<-hat_thetaG_glm$coefficients[1:pG]
    }
    XAG<-cbind(X,A,G)
    Xid<-1:pX
    Aid<-(pX+1):(pX+pA)
    Gid<-(pX+pA+1):(pX+pA+pG)
    fix_penalty<-c(rep(1,pX+pA),rep(0,pG))
    if(remove_penalty_X){fix_penalty[Xid]<-0}
    if(remove_penalty_A){fix_penalty[Aid]<-0}
    if(!is.null(fix_ratio)){fix_penalty[Xid]<-fix_ratio}
    if(!is.null(beta_initial)){
        if(length(beta_initial)!=pX+pA+pG){
            warning("beta_initial should be from X,A,G.\n Length not match, compute default initial instead.")
            beta_initial<-NULL
        }
    }
    if(is.null(beta_initial) & initial_with_type %in% c("ridge","lasso")){
        if(initial_with_type == "ridge"){initial_alpha=0}else{initial_alpha=1}
        if(pG == 1 & unique(G) == 1){
            fit_initial<-cv.glmnet(x=cbind(X,A),y=y,alpha = initial_alpha,penalty.factor = fix_penalty[c(Xid,Aid)],family="binomial")
            beta_initial<-c(coef.glmnet(fit_initial,s="lambda.min")[-1])
            beta_initial<-c(beta_initial,coef.glmnet(fit_initial,s="lambda.min")[1])
        }else if(pG > 1 & unique(G[,1])==1){
            fit_initial<-cv.glmnet(x=cbind(X,A,G[,-1]),y=y,alpha = initial_alpha,penalty.factor = fix_penalty[c(Xid,Aid,Gid[-1])],family="binomial")
            beta_initial<-c(coef.glmnet(fit_initial,s="lambda.min")[-1])
            beta_initial<-c(beta_initial[c(Xid,Aid)],coef.glmnet(fit_initial,s="lambda.min")[1],beta_initial[c(Gid[-1])])
        }else{
            stop("The first column of G should be 1 for intercept.")
        }
    }else{stop("Select Initial Type from c('ridge','lasso')")}
    if (penalty_type == "adaptivelasso"){
        w_adaptive<-1/abs(beta_initial)^gamma_adaptivelasso
        w_adaptive[is.infinite(w_adaptive)]<-max(w_adaptive[!is.infinite(w_adaptive)])*100
        w_adaptive<-w_adaptive*fix_penalty
    }else{w_adaptive<-fix_penalty}
    # Estimation of C
    inv_C = inv_C_func(X=X,A=A,G=G,y=y,
                       beta=beta_initial,hat_thetaG=hat_thetaG,
                       study_info=study_info)
    if(use_sparseC){
        C_half<-diag(1/sqrt(diag(inv_C)))
    }else{
        inv_C_svd<-fast.svd(inv_C+diag(1e-15,nrow(inv_C)))
        C_half<-inv_C_svd$v%*%diag(1/sqrt(inv_C_svd$d))%*%t(inv_C_svd$u)
    }

    # Prepare for final model
    pseudo_Xy_list<-pseudo_Xy(C_half,X,A,G,y,beta = beta_initial,hat_thetaG = hat_thetaG,study_info=study_info)
    initial_sf<-nX/sqrt(nrow(pseudo_Xy_list$pseudo_X))
    pseudo_X<-pseudo_Xy_list$pseudo_X/initial_sf
    pseudo_y<-pseudo_Xy_list$pseudo_y/initial_sf

    # Fit final model
    fit_final<-cv.glmnet(x= pseudo_X,y= pseudo_y,standardize=F,
                         intercept=F,alpha = final_alpha,penalty.factor = w_adaptive)
    if(is.null(lambda_list)){lambda_list<-fit_final$lambda}
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
        index_valid<-createDataPartition(y,p = holdout_p)$Resample1
        Xtrain<-X[-index_valid,]
        Xtest<-X[index_valid,]
        if(!is.null(A)){
            Atrain<-A[-index_valid,]
            Atest<-A[index_valid,]
        }else{
            Atrain<-NULL
            Atest<-NULL}
        if(!is.null(G)){
            Gtrain<-G[-index_valid,]
            Gtest<-G[index_valid,]
        }else{
            Gtrain<-NULL
            Gtest<-NULL}
        ytrain<-y[-index_valid]
        ytest<-y[index_valid]

        pseudo_Xy_list_train<-pseudo_Xy(C_half,Xtrain,Atrain,Gtrain,
                                        ytrain,beta = beta_initial,hat_thetaG = hat_thetaG,
                                        study_info=study_info)
        initial_sf_train<-nrow(Xtrain)/sqrt(nrow(pseudo_Xy_list_train$pseudo_X))
        pseudo_X_train<-pseudo_Xy_list_train$pseudo_X/initial_sf_train
        pseudo_y_train<-pseudo_Xy_list_train$pseudo_y/initial_sf_train


        if(tune_ratio){
            holdout_auc<-holdout_auc_lambda_ratio_func(lambda_list,ratio_range,pX,pA,pG,
                                   w_adaptive,pseudo_X_train,pseudo_y_train,
                                   final_alpha,Xtest,Atest,Gtest,ytest)
            ids<-which(holdout_auc==max(holdout_auc),arr.ind = TRUE)
            final.ratio.min<-ratio_range[ids[1]]
            final.lambda.min<-lambda_list[ids[2]]
        }else{
            holdout_auc<-holdout_auc_lambda_func(lambda_list,pseudo_X_train,pseudo_y_train,
                                 final_alpha,w_adaptive,Xtest,Atest,Gtest,ytest)
            final.lambda.min<-lambda_list[which.max(holdout_auc)]
            final.ratio.min<-1
        }
        ratio_vec<-c(rep(final.ratio.min,pX),rep(1,pA+pG))
        w_adaptive_ratio<-w_adaptive*ratio_vec
        fit_final_lam_ratio<-glmnet(x= pseudo_X_train,y= pseudo_y_train,standardize=F,
                                    intercept=F,alpha = final_alpha,
                                    penalty.factor = w_adaptive_ratio,
                                    lambda = final.lambda.min)
        beta<-coef.glmnet(fit_final_lam_ratio)[-1]
        return_list<-list("beta"=beta,#[-c((length(beta)-pG+1):length(beta)) ],
                          "lambda_list"=lambda_list,
                          "ratio_list"=ratio_range,
                          "holdout_auc"= holdout_auc,
                          "lambda_min"=final.lambda.min,
                          "ratio_min"=final.ratio.min)

    }else if(validation_type == 'cv'){
        index_fold<-createFolds(y,k = nfolds)

        if(tune_ratio){
            cv_auc<-cv_auc_lambda_ratio_func(index_fold,X,A,G,y,
                                             C_half,beta_initial,hat_thetaG,
                                             study_info,lambda_list,
                                             ratio_range,pX,pA,pG,
                                             w_adaptive,final_alpha,pseudo_Xy)
            ids<-which(cv_auc==max(cv_auc),arr.ind = TRUE)
            final.ratio.min<-ratio_range[ids[1]]
            final.lambda.min<-lambda_list[ids[2]]
        }else{
            cv_auc<-cv_auc_lambda_func(index_fold,X,A,G,y,
                                       C_half,beta_initial,hat_thetaG,
                                       study_info,lambda_list,
                                       w_adaptive,final_alpha,pseudo_Xy)
            final.lambda.min<-lambda_list[which.max(cv_auc)]
            final.ratio.min<-1
        }
        ratio_vec<-c(rep(final.ratio.min,pX),rep(1,pA+pG))
        w_adaptive_ratio<-w_adaptive*ratio_vec
        fit_final_lam_ratio<-glmnet(x= pseudo_X,y= pseudo_y,standardize=F,
                                    intercept=F,alpha = final_alpha,
                                    penalty.factor = w_adaptive_ratio,
                                    lambda = final.lambda.min)

        beta<-coef.glmnet(fit_final_lam_ratio)[-1]
        return_list<-list("beta"=beta,#[-c((length(beta)-pG+1):length(beta)) ],
                          "lambda_list"=lambda_list,
                          "ratio_list"=ratio_range,
                          "cv_auc"=cv_auc,
                          "lambda_min"=final.lambda.min,
                          "ratio_min"=final.ratio.min)
    }

    index_nonzero<-which(beta!=0)
    if(inference & length(index_nonzero) > 1){
        # refine C

        inv_C = inv_C_func(X=X,A=A,G=G,y=y,
                           beta=beta,hat_thetaG=hat_thetaG,
                           study_info=study_info)
        inv_C_svd<-fast.svd(inv_C+diag(1e-15,nrow(inv_C)))
        C_half<-inv_C_svd$v%*%diag(1/sqrt(inv_C_svd$d))%*%t(inv_C_svd$u)
        XAG<-cbind(X,A,G)
        expit_beta<-c(expit(XAG%*%beta))
        dexpit_beta<-expit_beta*(1-expit_beta)
        pseudo_X<-C_half%*%rbind(t(XAG),t(X))%*%(XAG*c(dexpit_beta))/nX
        #Sigsum_half<-cbind(xagtxag/nX,crossprod(XAG,X)/nX)%*%C_half
        #Sigsum_scaled<-Sigsum_half%*%t(Sigsum_half)
        Sigsum_scaled<-crossprod(pseudo_X)
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
