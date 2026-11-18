options(stringsAsFactors=FALSE)
ll<-readLines('work/sensitivity_20261112/analyse.R');eval(parse(text=paste(ll[seq_len(grep('^all_runs<-',ll)[1]-1L)],collapse='\n')))
b<-'work/nmr_20261113';aaall<-read.csv(file.path(b,'arms.csv'));ccall<-read.csv(file.path(b,'covariates.csv'));variables<-read.csv(file.path(b,'variables.csv'))
rank_basis<-function(X,tol=1e-9){
 z<-svd(X,nu=0,nv=min(dim(X)));keep<-z$d>max(z$d)*tol;z$v[,keep,drop=FALSE]
}
estimable<-function(l,G)sum((l-drop(G%*%crossprod(G,l)))^2)<1e-10*max(1,sum(l^2))
infer<-function(l,beta,C,df,G=NULL){
 if(!is.null(G)&&!estimable(l,G))return(c(estimate=NA,se=NA,lower=NA,upper=NA,p=NA))
 v<-sum(l*beta);se<-sqrt(max(0,drop(l%*%C%*%l)));q<-qt(.975,df)
 c(estimate=v,se=se,lower=v-q*se,upper=v+q*se,p=2*pt(-abs(v/se),df))
}
buildfit<-function(a,c,v,kind='Treatment_specific',min_k=10){
 c<-c[!is.na(c[[v$variable]]),];a<-a[a$id%in%c$id,];inp<-tryCatch(input(a,spec$sm[match(unique(a$outcome),spec$key)]),error=function(e)e)
 if(inherits(inp,'error'))return(list(status=conditionMessage(inp)))
 c<-c[match(unique(inp$d$id),c$id),];x0<-median(c[[v$variable]]);xx<-(c[[v$variable]]-x0)/v$increment;names(xx)<-c$id
 if(nrow(c)<min_k)return(list(status='Fewer than 10 complete NAT-connected report IDs',input=inp))
 if(length(unique(xx))<2)return(list(status='No covariate variation',input=inp))
 xc<-xx[inp$ct$id];X<-inp$X;p<-ncol(X)
 if(kind=='Treatment_specific'){
  W<-X*as.numeric(xc);H<-cbind(X,W);colnames(H)<-c(paste0('d_',colnames(X)),paste0('beta_',colnames(X)))
  residualW<-qr.resid(qr(X),W);A<-rank_basis(residualW);r<-ncol(A)
  # rank_basis treats a numerically zero matrix specially below.
  if(max(abs(residualW))<1e-9)r<-0L
  if(r==0)return(list(status='No interaction function identifiable separately from treatment effects',input=inp))
  L<-cbind(matrix(0,r,p),t(A));G<-rank_basis(H);design<-H%*%G
 }else{
  z<-(as.integer(inp$ct$treatment!='NAT')-as.integer(inp$ct$baseline!='NAT'))*as.numeric(xc)
  H<-cbind(X,beta_common=z);colnames(H)<-c(paste0('d_',colnames(X)),'beta_common');G<-diag(ncol(H));design<-H;r<-1L;L<-matrix(c(rep(0,p),1),1)
  if(qr(design,tol=1e-9)$rank<ncol(design))return(list(status='Common slope not identifiable',input=inp))
 }
 rdf<-nrow(design)-ncol(design)
 if(rdf<5)return(list(status='Fewer than 5 residual contrast degrees of freedom',input=inp))
 f<-tryCatch(fit_gls(inp,design),error=function(e)e);if(inherits(f,'error'))return(list(status=conditionMessage(f),input=inp))
 beta<-drop(G%*%f$beta);C<-G%*%f$vcov%*%t(G);z<-drop(L%*%beta);CV<-L%*%C%*%t(L);F<-drop(crossprod(z,solve(CV,z)))/r
 un<-fit_gls(inp);pglobal<-pf(F,r,rdf,lower.tail=FALSE)
 list(status='Estimated',input=inp,cov=c,x0=x0,v=v,kind=kind,X=H,G=G,fit=f,beta=beta,vcov=C,p=p,interaction_rank=r,L=L,df=rdf,F=F,pglobal=pglobal,unadjusted=un)
}
summaries<-slopes<-predictions<-nodecover<-loo<-list();objects<-list()
record<-function(m,key,v,kind,scenario,allc,removed=character()){
 id<-paste(key,v$variable,kind,scenario,sep='|');ip<-m$input
 row<-data.frame(model_id=id,outcome=key,variable=v$variable,kind=kind,scenario=scenario,k_total=nrow(allc),k_complete=sum(!is.na(allc[[v$variable]])),k_connected=if(is.null(ip))0 else length(unique(ip$d$id)),n=if(is.null(ip))0 else sum(ip$d$n),nodes=if(is.null(ip))0 else length(ip$tt),missing_ids=paste(allc$id[is.na(allc[[v$variable]])],collapse=';'),removed_ids=paste(removed,collapse=';'),disconnected_ids=if(is.null(ip))''else paste(ip$excluded_ids,collapse=';'),x0=NA_real_,x_min=NA_real_,x_max=NA_real_,unique_x=NA_integer_,interaction_df=NA_integer_,residual_df=NA_integer_,F=NA_real_,p_global=NA_real_,tau2_unadjusted=NA_real_,tau2_adjusted=NA_real_,tau2_change=NA_real_,tau2_reduction_fraction=NA_real_,shared_beta=NA_real_,shared_lower=NA_real_,shared_upper=NA_real_,shared_p=NA_real_,status=m$status)
 if(m$status=='Estimated'){
  x<-m$cov[[v$variable]];row$x0<-m$x0;row$x_min<-min(x);row$x_max<-max(x);row$unique_x<-length(unique(x));row$interaction_df<-m$interaction_rank;row$residual_df<-m$df;row$F<-m$F;row$p_global<-m$pglobal;row$tau2_unadjusted<-m$unadjusted$tau2;row$tau2_adjusted<-m$fit$tau2;row$tau2_change<-m$fit$tau2-m$unadjusted$tau2
  if(m$unadjusted$tau2>1e-10)row$tau2_reduction_fraction<-(m$unadjusted$tau2-m$fit$tau2)/m$unadjusted$tau2
  tt<-colnames(m$input$X);pp<-m$p
  if(kind=='Common_active'){
   z<-infer(c(rep(0,pp),1),m$beta,m$vcov,m$df);row$shared_beta<-z['estimate'];row$shared_lower<-z['lower'];row$shared_upper<-z['upper'];row$shared_p<-z['p']
  }else{
   so<-lapply(seq_along(tt),function(j){
    l<-rep(0,2*pp);l[pp+j]<-1;z<-infer(l,m$beta,m$vcov,m$df,m$G);ci<-m$cov[[v$variable]][m$cov$id%in%m$input$d$id[m$input$d$t==tt[j]]]
    data.frame(model_id=id,outcome=key,variable=v$variable,scenario=scenario,treatment=tt[j],reference='NAT',increment=v$increment,estimate=z['estimate'],se=z['se'],lower=z['lower'],upper=z['upper'],p=z['p'],estimable=!is.na(z['estimate']),node_k=length(ci),node_unique_x=length(unique(ci)),node_x_min=min(ci),node_x_max=max(ci),p_global=m$pglobal,residual_df=m$df)
   });slopes[[id]]<<-do.call(rbind,so)
   pr<-list();q<-unique(as.numeric(quantile(x,c(.25,.5,.75))))
   for(j in seq_along(tt))for(xvalue in q){
    ci<-m$cov[[v$variable]][m$cov$id%in%m$input$d$id[m$input$d$t==tt[j]]];inside<-xvalue>=min(ci)-1e-9&&xvalue<=max(ci)+1e-9
    l<-rep(0,2*pp);l[j]<-1;l[pp+j]<-(xvalue-m$x0)/v$increment;z<-if(inside)infer(l,m$beta,m$vcov,m$df,m$G)else c(estimate=NA,se=NA,lower=NA,upper=NA,p=NA)
    bu<-infer(replace(rep(0,pp),j,1),m$unadjusted$beta,m$unadjusted$vcov,nrow(m$input$X)-pp)
    pr[[length(pr)+1L]]<-data.frame(model_id=id,outcome=key,variable=v$variable,scenario=scenario,treatment=tt[j],reference='NAT',x_value=xvalue,inside_node_range=inside,estimate=z['estimate'],lower=z['lower'],upper=z['upper'],unadjusted_estimate=bu['estimate'],unadjusted_lower=bu['lower'],unadjusted_upper=bu['upper'],status=if(!inside)'Outside treatment covariate support; not predicted'else if(is.na(z['estimate']))'Not identifiable at this covariate value'else'Estimated')
   };predictions[[id]]<<-do.call(rbind,pr)
  }
  objects[[id]]<<-m
 }
 summaries[[id]]<<-row
}
for(key in spec$key){
 a<-aaall[aaall$outcome==key,];c<-ccall[ccall$outcome==key,];bf<-fit_gls(input(a,spec$sm[match(key,spec$key)]));ref<-selected[[key]]
 stopifnot(max(abs(bf$beta-ref$fit$beta))<1e-7,abs(bf$tau2-ref$fit$tau2)<1e-7)
 for(i in seq_len(nrow(variables))){
  v<-variables[i,];complete_ids<-c$id[!is.na(c[[v$variable]])];rem<-character()
  if(all(c('RPT023','RPT184')%in%complete_ids))rem<-c(rem,'RPT184')
  if(all(c('RPT092_306','RPT141')%in%complete_ids))rem<-c(rem,'RPT092_306')
  scenarios<-list(R050=character());if(length(rem))scenarios[['Cohort_cautious']]<-rem
  for(sc in names(scenarios))for(kind in c('Treatment_specific','Common_active')){
   rr<-scenarios[[sc]];ac<-a[!a$id%in%rr,];cc<-c[!c$id%in%rr,];m<-buildfit(ac,cc,v,kind);record(m,key,v,kind,sc,c,rr)
  }
  # Node support before fitting, including missingness and all possible levels.
  for(t in sort(unique(a$t))){
   ids<-unique(a$id[a$t==t]);x<-c[[v$variable]][match(ids,c$id)];nodecover[[paste(key,v$variable,t)]]<-data.frame(outcome=key,variable=v$variable,treatment=t,node_k=length(ids),complete_k=sum(!is.na(x)),unique_x=length(unique(x[!is.na(x)])),x_min=if(all(is.na(x)))NA else min(x,na.rm=TRUE),x_max=if(all(is.na(x)))NA else max(x,na.rm=TRUE),ids=paste(ids,collapse=';'))
  }
 }
 cat(key,'primary and cohort models done\n');flush.console()
}
summaries_df<-do.call(rbind,summaries);summaries_df$q_BH<-NA_real_
for(sc in unique(summaries_df$scenario))for(kind in unique(summaries_df$kind))for(core in c(0,1)){
 vars<-variables$variable[variables$core==core];ii<-which(summaries_df$scenario==sc&summaries_df$kind==kind&summaries_df$variable%in%vars&!is.na(summaries_df$p_global));summaries_df$q_BH[ii]<-p.adjust(summaries_df$p_global[ii],method='BH')
}
# Leave-one-report-ID-out fits for the two core variables, with no selection by P value.
for(id in names(objects)){
 mo<-objects[[id]];bits<-strsplit(id,'|',fixed=TRUE)[[1]];if(bits[4]!='R050'||!mo$v$core)next
 key<-bits[1];a<-aaall[aaall$outcome==key&aaall$id%in%mo$cov$id,];c<-mo$cov
 for(removed in mo$cov$id){
  fit<-buildfit(a[a$id!=removed,],c[c$id!=removed,],mo$v,mo$kind);z<-rep(NA_real_,5);names(z)<-c('estimate','se','lower','upper','p');pag<-z
  if(fit$status=='Estimated'){
   if(mo$kind=='Common_active')l<-c(rep(0,fit$p),1)else{l<-rep(0,2*fit$p);j<-match('BAL',colnames(fit$input$X));if(!is.na(j))l[fit$p+j]<-1 else l<-NULL}
   if(!is.null(l))z<-infer(l,fit$beta,fit$vcov,fit$df,if(mo$kind=='Treatment_specific')fit$G else NULL)
   if(mo$kind=='Treatment_specific'&&'PAG'%in%colnames(fit$input$X)){
    l<-rep(0,2*fit$p);l[fit$p+match('PAG',colnames(fit$input$X))]<-1;pag<-infer(l,fit$beta,fit$vcov,fit$df,fit$G)
   }
  }
  lo<-data.frame(model_id=id,outcome=key,variable=mo$v$variable,kind=mo$kind,removed_id=removed,status=fit$status,k=if(is.null(fit$input))0 else length(unique(fit$input$d$id)),p_global=if(fit$status=='Estimated')fit$pglobal else NA,tau2=if(fit$status=='Estimated')fit$fit$tau2 else NA,interaction_df=if(fit$status=='Estimated')fit$interaction_rank else NA,tracked_slope=if(mo$kind=='Common_active')'All active vs NAT common slope'else'BAL vs NAT slope',estimate=z['estimate'],lower=z['lower'],upper=z['upper'],p=z['p'],lost_nodes=if(is.null(fit$input))paste(mo$input$tt,collapse=';')else paste(setdiff(mo$input$tt,fit$input$tt),collapse=';'))
  lo$PAG_beta<-pag['estimate'];lo$PAG_lower<-pag['lower'];lo$PAG_upper<-pag['upper'];lo$PAG_p<-pag['p'];loo[[paste(id,removed)]]<-lo
 }
 cat('LOO',id,'done\n');flush.console()
}
for(pair in list(c('summary','summaries_df'),c('slopes','slopes'),c('predictions','predictions'),c('node_support','nodecover'),c('leave_one_out','loo'))){x<-get(pair[2]);if(is.list(x)&&!is.data.frame(x))x<-do.call(rbind,x);write.csv(x,file.path(b,paste0(pair[1],'.csv')),row.names=FALSE,na='')}
saveRDS(objects,file.path(b,'models.rds'))
# Independent package checks: re-estimate REML with metafor for every R050 fit.
suppressPackageStartupMessages(library(metafor));checks<-list()
for(id in names(objects)){
 if(!endsWith(id,'|R050'))next
 m<-objects[[id]];dat<-m$input$ct;dat$inner<-ave(seq_len(nrow(dat)),dat$id,FUN=seq_along)
 fit<-rma.mv(yi,V=m$input$V,mods=m$fit$X,intercept=FALSE,random=~inner|id,struct='CS',rho=.5,data=dat,method='REML',test='t',dfs='residual',control=list(rel.tol=1e-10))
 err_beta<-max(abs(drop(fit$beta)-m$fit$beta));err_se<-max(abs(sqrt(diag(fit$vb))-sqrt(diag(m$fit$vcov))));err_tau<-abs(fit$tau2-m$fit$tau2)
 stopifnot(err_beta<1e-4,err_se<1e-4,err_tau<1e-4)
 checks[[id]]<-data.frame(model_id=id,beta_error=err_beta,se_error=err_se,tau2_error=err_tau)
}
write.csv(do.call(rbind,checks),file.path(b,'independent_verification.csv'),row.names=FALSE);writeLines(capture.output(sessionInfo()),file.path(b,'sessionInfo.txt'))
cat('Done:',length(objects),'fitted models;',length(loo),'LOO records;',length(checks),'independent REML package checks\n')
