#' htlgmm
#'
#' @details htlgmm: Transfer Learning via generalized method of moments(GMM).
#'
#' @param y The y for response variable, which can be continous or binary.
#' @param X The matched features for internal and external data.
#' @param Z The mismatched features only in internal data, the default is NULL.
#' @param study_info The summary statistics for X only from external data,
#' which can be summarized in the type of multivariate version or univariate version.
#' @param summary_type The summary statistics type, chosen from c("multi","uni"), where the default is "multi".
#' @param family The family is chosen from c("gaussian","binomial"). Linear regression for "gaussian" and logistic regression for "binomial".
#' @param A Usually used in "binomial" family, e.g. the intercept term for logistic regression, where the default is 1.
#' Usually not used in "gaussian" family.
#' G are the features working for adjustment in reduced model, but G is not summarized in summary statistics(input:study_info).
#' @param penalty_type The penalty type for htlgmm, chosen from c("lasso","adaptivelasso","ridge"). The default is "lasso".
#' @param initial_with_type Get initial estimation for beta using internal data only
#'  by cross validation using penalty regression, chosen from c("ridge","lasso"). The default is "ridge".
#' @param beta_initial The initial estimation for beta if a consistent estimator is available.
#' E.g., one may input htlgmm result as beta_initial for more rounds to refine the final estimation.
#' The default is NULL, and internal data is used for initial estimation.
#' @param remove_penalty_X Not penalize X if it is TRUE. The default is FALSE.
#' @param remove_penalty_Z Not penalize Z if it is TRUE. The default is FALSE.
#' @param lambda Without cross validation, fix the lambda. The default is NULL.
#' @param ratio The fixed ratio of X for bi-lambda strategy. The default is NULL. If it is NULL, select the best ratio via cross validation or holdout validation.
#' @param gamma_adaptivelasso The gamma for adaptive lasso. Select from c(1/2,1,2). The default is 1/2.
#' @param inference Whether to do post-selection inference, only work for adaptivelasso. The default is FALSE.
#' @param use_sparseC Whether to use approximate version of weighting matrix C.
#' If approximation, use the diagonal of inverse of C(inv_C) to approximate the inv_C. The default is FALSE.
#' When internal data sample size is limited, use_sparseC = TRUE is recommended.
#' When internal data sample size is large enough, use_sparseC = FALSE is recommended.
#'
#' @return beta The target coefficient estimation.
#' @return fix_lambda The output lambda.
#' @return fix_ratio The output ratio.
#' @return corrected_pos For post-selection inference, they are the corrected position passing significant level 0.05 after BH adjustment (Benjamini & Hochberg).
#' @return nonzero_pos For estimated beta, the nonzero positions.
#' @return pval For nonzero_pos, the calculated p values.
#' @return nonzero_var For nonzero_pos, the calculated variances.
#' @return \itemize{
#'  \item{beta:} The target coefficient estimation.
#'  \item{fix_lambda:} The output lambda.
#'  \item{fix_ratio:} The output ratio.
#'  \item{corrected_pos:} For post-selection inference, they are the corrected position passing significant level 0.05 after BH adjustment (Benjamini & Hochberg).
#'  \item{nonzero_pos:} For estimated beta, the nonzero positions.
#'  \item{pval:} For nonzero_pos, the calculated p values.
#'  \item{nonzero_var:} For nonzero_pos, the calculated variances.
#' }
#'
#' @import glmnet
#' @import stats
#' @importFrom caret createFolds createDataPartition
#' @importFrom locfit expit
#' @importFrom corpcor fast.svd
#' @importFrom magic adiag
#' @importFrom MASS ginv
#' @export
#'
#' @examples
#' set.seed(1)
#' X<-matrix(rnorm(18000),900,20)
#' Z<-matrix(rnorm(2700),900,3)
#' X<-scale(X)
#' Z<-scale(Z)
#' coefXZ<-c(rep(0,23))
#' coefXZ[1:3]<-0.5
#' coefXZ[21:22]<-0.5
#' internal_index<-1:100
#' external_index<-101:900
#' y<-cbind(X,Z)%*%coefXZ+rnorm(900,0,1)
#' #y_binary<-rbinom(n=900,size=1,prob=locfit::expit(cbind(X,Z)%*%coefXZ))
#' study_info_multi<-list()
#' reslm<-lm(y~.,data = data.frame(y=y[external_index],X[external_index,]))
#' study.m = list(Coeff=reslm$coefficients[-1],
#'                Covariance=vcov(reslm)[-1,-1],Sample_size=800)
#' study_info_multi[[1]] <- study.m
#' study_info_uni<-list()
#' for(i in 1:20){
#' reslm<-lm(y~.,data = data.frame(y=y[external_index],X[external_index,i]))
#' study.m = list(Coeff=reslm$coefficients[-1],
#'                Covariance=vcov(reslm)[-1,-1],Sample_size=800)
#' study_info_uni[[i]] <- study.m}
#'
#' y<-scale(y,scale = FALSE)
#' library(glmnet)
#' res_glm<-glmnet(x=cbind(X[internal_index,],A[internal_index,]),y=y[internal_index],lambda=0)
#' res_htlgmm_multi<-htlgmm(y[internal_index],X[internal_index,],A[internal_index,],
#'     summary_type = "multi",study_info = study_info_multi,lambda=0,use_sparseC = TRUE)
#' res_htlgmm_uni<-htlgmm(y[internal_index],X[internal_index,],A[internal_index,],
#'     summary_type = "uni",study_info = study_info_uni,lambda=0,use_sparseC = TRUE)
#'     ee_lasso<-round(sum((coefXZ-coef.glmnet(res_glm)[-1])^2),4)
#' ee_htlgmm_lasso_multi<-round(sum((coefXZ-res_htlgmm_multi$beta)^2),4)
#' ee_htlgmm_lasso_uni<-round(sum((coefXZ-res_htlgmm_uni$beta)^2),4)
#' print(paste0("Estimation Error: ","lasso(",ee_lasso,"); htlgmm_lasso_multi(",
#'              ee_htlgmm_lasso_multi,"); htlgmm_lasso_uni(",ee_htlgmm_lasso_uni,")"))
#'
#'

htlgmm<-function(
        y,X,Z=NULL,
        study_info=NULL,
        summary_type = "multi",
        family = "gaussian",
        A=1,
        penalty_type = "lasso",
        initial_with_type = "ridge",
        beta_initial = NULL,
        remove_penalty_X = FALSE,
        remove_penalty_Z = FALSE,
        #tune_ratio = TRUE,
        lambda = 0,
        #lambda_list = NULL,
        ratio = NULL,
        #ratio_lower = NULL,
        #ratio_upper = NULL,
        #ratio_count = 10,
        #ratio_range = NULL,
        gamma_adaptivelasso = 1/2,
        inference = FALSE,
        #validation_type = "cv",
        #nfolds = 10,
        #holdout_p = 0.2,
        use_sparseC = FALSE
){
    tune_ratio = FALSE
    lambda_list = NULL
    ratio_lower = NULL
    ratio_upper = NULL
    ratio_count = 10
    ratio_range = NULL
    validation_type = "None"
    nfolds = 10
    holdout_p = 0.2
    fix_lambda = lambda
    fix_ratio = ratio
    seed.use = 97
    if(!family %in% c("gaussian","binomial")){
        stop("Select family from c('gaussian','binomial')")
    }
    if(family == "gaussian"){
        #warning("For gaussian family, no A is used. \n Just assume A is cancelled for full and reduced model.")
        res<-htlgmm.linear(y,X,Z,study_info,summary_type,penalty_type,
                           initial_with_type,beta_initial,
                           remove_penalty_X,remove_penalty_Z,
                           tune_ratio,fix_lambda,lambda_list,
                           fix_ratio,ratio_lower,ratio_upper,
                           ratio_count,ratio_range,gamma_adaptivelasso,
                           inference,validation_type,
                           nfolds,holdout_p,use_sparseC,seed.use)
    }else{
        res<-htlgmm.binary(y,X,Z,A,study_info,summary_type,penalty_type,
                           initial_with_type,beta_initial,
                           remove_penalty_X,remove_penalty_Z,
                           tune_ratio,fix_lambda,lambda_list,
                           fix_ratio,ratio_lower,ratio_upper,
                           ratio_count,ratio_range,gamma_adaptivelasso,
                           inference,validation_type,nfolds,holdout_p,
                           use_sparseC,seed.use)
    }
    return(res)
}



#' cv.htlgmm: cross validtaion for htlgmm(integration with generalized method of moments).
#'
#' @details cross validtaion for htlgmm: integration with generalized method of moments(gmm).
#'
#' @param y The y for response variable, which can be continouse or binary.
#' @param X The matched features for internal and external data.
#' @param Z The mismatched features only in internal data, the default is NULL.
#' @param study_info The summary statistics for X only from external data,
#' which can be summarized in the type of multivariate version or univariate version.
#' @param summary_type The summary statistics type, chosen from c("multi","uni"), where the default is "multi".
#' @param family The family is chosen from c("gaussian","binomial"). Linear regression for "gaussian" and logistic regression for "binomial".
#' @param A Usually used in "binomial" family, e.g. the intercept term for logistic regression, where the default is 1.
#' Usually not used in "gaussian" family.
#' G are the features working for adjustment in reduced model, but G is not summarized in summary statistics(input:study_info).
#' @param penalty_type The penalty type for htlgmm, chosen from c("lasso","adaptivelasso","ridge"). The default is "lasso".
#' @param initial_with_type Get initial estimation for beta using internal data only
#'  by cross validation using penalty regression, chosen from c("ridge","lasso"). The default is "ridge".
#' @param beta_initial The initial estimation for beta if a consistent estimator is available.
#' E.g., one may input htlgmm result as beta_initial for more rounds to refine the final estimation.
#' The default is NULL, and internal data is used for initial estimation.
#' @param remove_penalty_X Not penalize X if it is TRUE. The default is FALSE.
#' @param remove_penalty_Z Not penalize Z if it is TRUE. The default is FALSE.
#' @param tune_ratio Whether to use bi-lambda stratgey. The default is TRUE.
#' @param fix_lambda Without cross validation, fix the lambda. The default is NULL.
#' @param lambda_list Customize the input lambda list for validation. The default is NULL.
#' @param fix_ratio The fixed ratio of X for bi-lambda strategy. The default is NULL. If it is NULL, select the best ratio via cross validation or holdout validation.
#' @param ratio_lower The lower bound for ratio range. The default is NULL.
#' @param ratio_upper The upper bound for ratio range. The default is NULL.
#' @param ratio_count The lengte of ratio list. The default is 10.
#' @param ratio_range The ratio range if it is preset. The default is NULL.
#' @param gamma_adaptivelasso The gamma for adaptive lasso. Select from c(1/2,1,2). The default is 1/2.
#' @param inference Whether to do post-selection inference, only work for adaptivelasso. The default is FALSE.
#' @param validation_type How to perform validation to find the best lamdba or ratio.
#' Select from c("cv","holdout"). The default is "cv".
#' @param nfolds The fold number for cross validation. Only work for validation_type = "cv".The default is 10.
#' @param holdout_p The holdout validation data proportion. Only work for validation_type = "holdout". The default is 0.2.
#' @param use_sparseC Whether to use approximate version of weighting matrix C.
#' If approximation, use the diagonal of inverse of C(inv_C) to approximate the inv_C. The default is FALSE.
#' When internal data sample size is limited, use_sparseC = TRUE is recommended.
#' When internal data sample size is large enough, use_sparseC = FALSE is recommended.
#' @param seed.use The seed for  97.
#'
#' @return \itemize{
#'  \item{beta:} The target coefficient estimation.
#'  \item{lambda_list:} The lambda list for validation (cross validation or holdout validation).
#'  \item{ratio_list:} The ratio list for validation (cross validation or holdout validation).
#'  \item{holdout_mse:} The mean square error(mse) when family = "gaussian", and validation_type = "holdout".
#'  \item{cv_mse:} The mean square error(mse) when family = "gaussian", and validation_type = "cv".
#'  \item{holdout_dev:} The deviance(dev) when family = "binomial", and validation_type = "holdout".
#'  \item{cv_dev:} The deviance(dev) when family = "binomial", and validation_type = "cv".
#'  \item{lambda_min:} The selected best lambda.
#'  \item{ratio_min:} The selected best ratio.
#'  \item{corrected_pos:} For post-selection inference, they are the corrected position passing significant level 0.05 after BH adjustment (Benjamini & Hochberg).
#'  \item{nonzero_pos:} For estimated beta, the nonzero positions.
#'  \item{pval:} For nonzero_pos, the calculated p values.
#'  \item{nonzero_var:} For nonzero_pos, the calculated variances.
#' }
#'
#'
#' @import glmnet
#' @import stats
#' @importFrom caret createFolds createDataPartition
#' @importFrom locfit expit
#' @importFrom corpcor fast.svd
#' @importFrom magic adiag
#' @importFrom MASS ginv
#' @export
#'
#' @examples
#' set.seed(1)
#' X<-matrix(rnorm(18000),900,20)
#' Z<-matrix(rnorm(2700),900,3)
#' X<-scale(X)
#' Z<-scale(Z)
#' coefXZ<-c(rep(0,23))
#' coefXZ[1:3]<-0.5
#' coefXZ[21:22]<-0.5
#' internal_index<-1:100
#' external_index<-101:900
#' y<-cbind(X,Z)%*%coefXZ+rnorm(900,0,1)
#' #y_binary<-rbinom(n=900,size=1,prob=locfit::expit(cbind(X,Z)%*%coefXZ))
#' study_info_multi<-list()
#' reslm<-lm(y~.,data = data.frame(y=y[external_index],X[external_index,]))
#' study.m = list(Coeff=reslm$coefficients[-1],
#'                Covariance=vcov(reslm)[-1,-1],Sample_size=800)
#' study_info_multi[[1]] <- study.m
#' study_info_uni<-list()
#' for(i in 1:20){
#' reslm<-lm(y~.,data = data.frame(y=y[external_index],X[external_index,i]))
#' study.m = list(Coeff=reslm$coefficients[-1],
#'                Covariance=vcov(reslm)[-1,-1],Sample_size=800)
#' study_info_uni[[i]] <- study.m}
#'
#' y<-scale(y,scale = FALSE)
#' library(glmnet)
#' res_glm<-cv.glmnet(x=cbind(X[internal_index,],A[internal_index,]),y=y[internal_index])
#' res_htlgmm_multi<-cv.htlgmm(y[internal_index],X[internal_index,],A[internal_index,],
#'     summary_type = "multi",study_info = study_info_multi,tune_ratio = FALSE,use_sparseC = TRUE)
#' res_htlgmm_uni<-cv.htlgmm(y[internal_index],X[internal_index,],A[internal_index,],
#'     summary_type = "uni",study_info = study_info_uni,tune_ratio = FALSE,use_sparseC = TRUE)
#' ee_lasso<-round(sum((coefXZ-coef.glmnet(res_glm,s="lambda.min")[-1])^2),4)
#' ee_htlgmm_lasso_multi<-round(sum((coefXZ-res_htlgmm_multi$beta)^2),4)
#' ee_htlgmm_lasso_uni<-round(sum((coefXZ-res_htlgmm_uni$beta)^2),4)
#' print(paste0("Estimation Error: ","lasso(",ee_lasso,"); htlgmm_lasso_multi(",
#'              ee_htlgmm_lasso_multi,"); htlgmm_lasso_uni(",ee_htlgmm_lasso_uni,")"))
#'
#'
cv.htlgmm<-function(
        y,X,Z=NULL,
        study_info=NULL,
        summary_type = "multi",
        family = "gaussian",
        A=1,
        penalty_type = "lasso",
        initial_with_type = "ridge",
        beta_initial = NULL,
        remove_penalty_X = FALSE,
        remove_penalty_Z = FALSE,
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
    if(!family %in% c("gaussian","binomial")){
        stop("Select family from c('gaussian','binomial')")
    }
    if(family == "gaussian"){
        #warning("For gaussian family, no A is used. \n Just assume A is cancelled for full and reduced model.")
        res<-htlgmm.linear(y,X,Z,study_info,summary_type,penalty_type,
                           initial_with_type,beta_initial,
                           remove_penalty_X,remove_penalty_Z,
                           tune_ratio,fix_lambda,lambda_list,
                           fix_ratio,ratio_lower,ratio_upper,
                           ratio_count,ratio_range,gamma_adaptivelasso,
                           inference,validation_type,
                           nfolds,holdout_p,use_sparseC,seed.use)
    }else{
        res<-htlgmm.binary(y,X,Z,A,study_info,summary_type,penalty_type,
                           initial_with_type,beta_initial,
                           remove_penalty_X,remove_penalty_Z,
                           tune_ratio,fix_lambda,lambda_list,
                           fix_ratio,ratio_lower,ratio_upper,
                           ratio_count,ratio_range,gamma_adaptivelasso,
                           inference,validation_type,nfolds,holdout_p,
                           use_sparseC,seed.use)
    }
    return(res)
}
