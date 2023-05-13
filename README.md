# README



This works for README file for `intgmm` function. 

### $p_{\bf X}=20,p_{\bf A}=3$

```R
# Simulation Example
set.seed(1)
X<-matrix(rnorm(18000),900,20)
A<-matrix(rnorm(2700),900,3)
X<-scale(X)
A<-scale(A)
coefXA<-c(rep(0,23))
coefXA[1:3]<-0.5
coefXA[21:22]<-0.5
internal_index<-1:100
external_index<-101:900
y<-cbind(X,A)%*%coefXA+rnorm(900,0,1)
#y_binary<-rbinom(n=900,size=1,prob=locfit::expit(cbind(X,A)%*%coefXA))
study_info_multi<-list()
reslm<-lm(y~.,data = data.frame(y=y[external_index],X[external_index,]))
study.m = list(Coeff=reslm$coefficients[-1],
               Covariance=vcov(reslm)[-1,-1],Sample_size=800)
study_info_multi[[1]] <- study.m
study_info_uni<-list()
for(i in 1:20){
reslm<-lm(y~.,data = data.frame(y=y[external_index],X[external_index,i]))
study.m = list(Coeff=reslm$coefficients[-1],
               Covariance=vcov(reslm)[-1,-1],Sample_size=800)
study_info_uni[[i]] <- study.m}

y<-scale(y,scale = FALSE)
library(glmnet)
res_glm<-cv.glmnet(x=cbind(X[internal_index,],A[internal_index,]),y=y[internal_index])
res_intgmm_multi<-cv.intgmm(y[internal_index],X[internal_index,],A[internal_index,],
    summary_type = "multi",study_info = study_info_multi,tune_ratio = FALSE,use_sparseC = TRUE)
res_intgmm_uni<-cv.intgmm(y[internal_index],X[internal_index,],A[internal_index,],
    summary_type = "uni",study_info = study_info_uni,tune_ratio = FALSE,use_sparseC = TRUE)
ee_lasso<-round(sum((coefXA-coef.glmnet(res_glm,s="lambda.min")[-1])^2),4)
ee_intgmm_lasso_multi<-round(sum((coefXA-res_intgmm_multi$beta)^2),4)
ee_intgmm_lasso_uni<-round(sum((coefXA-res_intgmm_uni$beta)^2),4)
print(paste0("Estimation Error: ","lasso(",ee_lasso,"); intgmm_lasso_multi(",
             ee_intgmm_lasso_multi,"); intgmm_lasso_uni(",ee_intgmm_lasso_uni,")"))

```

```R
"Estimation Error: lasso(0.2699); intgmm_lasso_multi(0.0866); intgmm_lasso_uni(0.1394)"
```



### $p_{\bf X}=500,p_{\bf A}=30$

```R
set.seed(1)
X<-matrix(rnorm(900*500),900,500)
A<-matrix(rnorm(900*30),900,30)
X<-scale(X)
A<-scale(A)
coefXA<-c(rep(0,530))
coefXA[1:5]<-0.5
coefXA[521:522]<-0.5
internal_index<-1:100
external_index<-101:900
y<-cbind(X,A)%*%coefXA+rnorm(900,0,1)
#y_binary<-rbinom(n=900,size=1,prob=locfit::expit(cbind(X,A)%*%coefXA))
study_info_multi<-list()
reslm<-lm(y~.,data = data.frame(y=y[external_index],X[external_index,]))
study.m = list(Coeff=reslm$coefficients[-1],
               Covariance=vcov(reslm)[-1,-1],Sample_size=800)
study_info_multi[[1]] <- study.m
study_info_uni<-list()
for(i in 1:500){
  reslm<-lm(y~.,data = data.frame(y=y[external_index],X[external_index,i]))
  study.m = list(Coeff=reslm$coefficients[-1],
                 Covariance=vcov(reslm)[-1,-1],Sample_size=800)
  study_info_uni[[i]] <- study.m}

y<-scale(y,scale = FALSE)
library(glmnet)
res_glm<-glmnet(x=cbind(X[internal_index,],A[internal_index,]),y=y[internal_index],lambda=0)
res_intgmm_multi<-intgmm(y[internal_index],X[internal_index,],A[internal_index,],
                         summary_type = "multi",study_info = study_info_multi,lambda=0,use_sparseC = F)
res_intgmm_uni<-intgmm(y[internal_index],X[internal_index,],A[internal_index,],
                       summary_type = "uni",study_info = study_info_uni,lambda=0,use_sparseC = F)
ee_lasso<-round(sum((coefXA-coef.glmnet(res_glm)[-1])^2),4)
ee_intgmm_lasso_multi<-round(sum((coefXA-res_intgmm_multi$beta)^2),4)
ee_intgmm_lasso_uni<-round(sum((coefXA-res_intgmm_uni$beta)^2),4)
print(paste0("Estimation Error: ","lasso(",ee_lasso,"); intgmm_lasso_multi(",
             ee_intgmm_lasso_multi,"); intgmm_lasso_uni(",ee_intgmm_lasso_uni,")"))
```

```R
"Estimation Error: lasso(1.9504); intgmm_lasso_multi(1.5171); intgmm_lasso_uni(1.2041)"
```

