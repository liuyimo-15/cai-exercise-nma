options(stringsAsFactors=FALSE)
.libPaths(c(normalizePath('work/cinema_20261112/Rlib'),.libPaths()))
suppressPackageStartupMessages(library(metafor))
b<-'work/nmr_20261113';models<-readRDS(file.path(b,'models.rds'));rows<-list();checks<-list()
for(id in names(models)){
 m<-models[[id]];k<-length(unique(m$input$ct$id));p<-ncol(m$fit$X);df<-k-p
 row<-data.frame(model_id=id,cluster_k=k,parameters=p,cluster_residual_df=df,p_global_CR1=NA_real_,shared_lower_CR1=NA_real_,shared_upper_CR1=NA_real_,shared_p_CR1=NA_real_,status='Cluster residual df <5: CR1 not reported')
 if(df>=5){
  dat<-m$input$ct;dat$inner<-ave(seq_len(nrow(dat)),dat$id,FUN=seq_along)
  f<-rma.mv(yi,V=m$input$V,mods=m$fit$X,intercept=FALSE,random=~inner|id,struct='CS',rho=.5,tau2=m$fit$tau2,data=dat,method='REML',test='t',dfs='residual')
  rr<-tryCatch(robust(f,cluster=dat$id,adjust=TRUE),error=function(e)e)
  if(inherits(rr,'error'))row$status<-conditionMessage(rr)else{
   # Some packages mark exactly zero-variance saturated basis directions NA.
   # Construct CR1 directly before transforming to identifiable slope contrasts.
   X<-m$fit$X;W<-solve(m$input$V+m$fit$tau2*m$input$R);e<-m$input$ct$yi-drop(X%*%m$fit$beta);we<-drop(W%*%e);meat<-matrix(0,p,p)
   for(cl in unique(dat$id)){ii<-which(dat$id==cl);score<-drop(crossprod(X[ii,,drop=FALSE],we[ii]));meat<-meat+tcrossprod(score)}
   CR<-m$fit$vcov%*%meat%*%m$fit$vcov*k/(k-p)
   ok<-is.finite(rr$vb);err<-max(abs(CR[ok]-rr$vb[ok]));stopifnot(err<1e-5);checks[[id]]<-data.frame(model_id=id,max_finite_covariance_error=err,package_NA_covariance_cells=sum(!ok))
   C<-m$G%*%CR%*%t(m$G);v<-drop(m$L%*%m$beta);CV<-m$L%*%C%*%t(m$L)
   if(any(!is.finite(CV)))row$status<-'Robust interaction covariance not finite: no global CR1 test'else if(qr(CV,tol=1e-9)$rank<nrow(CV))row$status<-'Robust interaction covariance singular: no global CR1 test'else{
    F<-drop(crossprod(v,solve(CV,v)))/m$interaction_rank;row$p_global_CR1<-pf(F,m$interaction_rank,df,lower.tail=FALSE);row$status<-'CR1 sensitivity at report-ID cluster level'
    if(m$kind=='Common_active'){
     j<-length(m$beta);se<-sqrt(C[j,j]);row$shared_lower_CR1<-m$beta[j]-qt(.975,df)*se;row$shared_upper_CR1<-m$beta[j]+qt(.975,df)*se;row$shared_p_CR1<-2*pt(-abs(m$beta[j]/se),df)
    }
   }
  }
 }
 rows[[id]]<-row
}
write.csv(do.call(rbind,rows),file.path(b,'robust.csv'),row.names=FALSE,na='')
write.csv(do.call(rbind,checks),file.path(b,'robust_verification.csv'),row.names=FALSE,na='')
cat('CR1 sensitivity evaluated for',length(rows),'models\n')
