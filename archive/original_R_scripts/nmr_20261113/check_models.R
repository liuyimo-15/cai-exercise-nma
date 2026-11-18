options(stringsAsFactors=FALSE)
b<-'work/nmr_20261113';m<-readRDS(file.path(b,'models.rds'));checks<-list()
for(id in names(m)){
 a<-m[[id]];ip<-a$input;W<-solve(ip$V+a$fit$tau2*ip$R);X<-ip$X;y<-ip$ct$yi
 beta0<-solve(crossprod(X,W%*%X),crossprod(X,W%*%y));e0<-y-drop(X%*%beta0);e1<-y-drop(a$fit$X%*%a$fit$beta)
 F_extra<-as.numeric((crossprod(e0,W%*%e0)-crossprod(e1,W%*%e1))/a$interaction_rank)
 err<-abs(F_extra-a$F);stopifnot(err<1e-7)
 # The fitted reduced model spans the full treatment-by-covariate design.
 pp<-a$p;X0<-ip$X;xc<-(a$cov[[a$v$variable]][match(ip$ct$id,a$cov$id)]-a$x0)/a$v$increment
 H<-if(a$kind=='Treatment_specific')cbind(X0,X0*xc)else cbind(X0,(as.integer(ip$ct$treatment!='NAT')-as.integer(ip$ct$baseline!='NAT'))*xc)
 pred_error<-max(abs(drop(H%*%a$beta)-drop(a$fit$X%*%a$fit$beta)));stopifnot(pred_error<1e-8)
 checks[[id]]<-data.frame(model_id=id,F_extra_error=err,prediction_parameterisation_error=pred_error)
}
write.csv(do.call(rbind,checks),file.path(b,'algebra_verification.csv'),row.names=FALSE)
cat('Independent extra-weighted-Q test and full/reduced prediction checks:',length(checks),'models\n')
