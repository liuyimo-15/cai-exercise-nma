options(stringsAsFactors=FALSE)
base<-'work/sensitivity_20261112';old<-'work/cinema_20261112'
source("R/bootstrap.R")
spec<-read.csv(file.path(base,'spec.csv'));ratings<-read.csv(file.path(base,'ratings.csv'))
fits<-readRDS(file.path(old,'candidate_fits.rds'));balance<-readRDS(file.path(old,'matched_balance.rds'))
selected<-c(list(CAIT=fits[['CAIT_raw_Cohen_N']],FAAM_ADL=fits[['FAAM_ADL_raw_Cohen_N']],Sports=fits[['Sports_raw_Hedges_N']]),balance)
q95<-qnorm(.975)
bdiag_base<-function(blocks){
 n<-sum(sapply(blocks,nrow));out<-matrix(0,n,n);i<-1L
 for(x in blocks){ii<-seq.int(i,length.out=nrow(x));out[ii,ii]<-x;i<-max(ii)+1L};out
}
input<-function(raw,sm){
 if(!nrow(raw)||!'NAT'%in%raw$t)stop('NAT absent')
 connected<-'NAT';repeat{before<-connected;ss<-unique(raw$id[raw$t%in%connected]);connected<-unique(raw$t[raw$id%in%ss]);if(setequal(before,connected))break}
 d<-raw[raw$t%in%connected,];if(length(unique(d$t))<2)stop('No treatment connected to NAT')
 vb<-rb<-ct<-list()
 for(id in unique(d$id)){
  a<-d[d$id==id,];a<-a[order(a$t),];k<-nrow(a);m<-k-1L;N<-sum(a$n);df<-N-k
  stopifnot(k>=2L)
  if(sm=='MD'){
   y<-a$mean[-1]-a$mean[1];va<-a$sd^2/a$n;V<-matrix(va[1],m,m);diag(V)<-va[-1]+va[1]
  }else{
   pooled<-sqrt(sum((a$n-1)*a$sd^2)/df);J<-exp(lgamma(df/2)-log(sqrt(df/2))-lgamma((df-1)/2))
   y<-J*(a$mean[-1]-a$mean[1])/pooled
   V<-matrix(J^2/a$n[1],m,m)+tcrossprod(y)/(2*N);diag(V)<-J^2*(1/a$n[-1]+1/a$n[1])+y^2/(2*N)
  }
  R<-matrix(.5,m,m);diag(R)<-1;vb[[id]]<-V;rb[[id]]<-R
  ct[[id]]<-data.frame(id=id,treatment=a$t[-1],baseline=a$t[1],yi=y)
 }
 ct<-do.call(rbind,ct);tt<-sort(unique(d$t));Z<-diag(length(tt));dimnames(Z)<-list(tt,tt);Z<-Z[,setdiff(tt,'NAT'),drop=FALSE]
 list(d=d,ct=ct,V=bdiag_base(vb),R=bdiag_base(rb),X=Z[ct$treatment,,drop=FALSE]-Z[ct$baseline,,drop=FALSE],Z=Z,tt=tt,sm=sm,excluded_ids=setdiff(unique(raw$id),unique(d$id)))
}
fit_gls<-function(inp,X=inp$X,method='REML'){
 y<-inp$ct$yi;V<-inp$V;R<-inp$R;n<-length(y);p<-ncol(X)
 if(qr(X,tol=1e-9)$rank<p||n<=p)stop('Rank deficient or no residual degrees of freedom')
 eval<-function(tau2,details=FALSE){
  U<-chol(V+tau2*R);W<-chol2inv(U);I<-crossprod(X,W%*%X);C<-solve(I);b<-drop(C%*%crossprod(X,W%*%y));names(b)<-colnames(X);e<-y-drop(X%*%b);Q<-drop(crossprod(e,W%*%e))
  obj<-.5*(2*sum(log(diag(U)))+as.numeric(determinant(I,logarithm=TRUE)$modulus)+Q)
  if(details)list(beta=b,vcov=C,tau2=tau2,Q=Q,n=n,p=p,X=X,objective=obj)else obj
 }
 if(method=='COMMON')tau2<-0 else{
  upper<-max(1,var(y),max(diag(V)));repeat{opt<-optimize(eval,c(0,upper),tol=1e-10);if(opt$minimum<.95*upper)break;upper<-upper*4;if(upper>1e10)stop('Tau2 optimisation bound')}
  tau2<-if(eval(0)<=opt$objective+1e-10)0 else opt$minimum
 }
 eval(tau2,TRUE)
}
component_fit<-function(inp,terms){
 Z<-sapply(terms,function(term){bits<-strsplit(term,':',fixed=TRUE)[[1]];vapply(strsplit(inp$tt,'+',fixed=TRUE),function(reg)as.numeric(all(bits%in%reg)),numeric(1))})
 rownames(Z)<-inp$tt;colnames(Z)<-terms
 X<-Z[inp$ct$treatment,,drop=FALSE]-Z[inp$ct$baseline,,drop=FALSE]
 ev<-eigen(crossprod(X),symmetric=TRUE);use<-ev$values>max(ev$values)*1e-9;A<-ev$vectors[,use,drop=FALSE]
 if(!ncol(A))stop('No estimable component functions')
 f<-fit_gls(inp,X%*%A);b<-drop(A%*%f$beta);names(b)<-terms;C<-A%*%f$vcov%*%t(A);dimnames(C)<-list(terms,terms)
 estimable<-colSums((diag(length(terms))-tcrossprod(A))^2)<1e-12
 list(fit=f,beta=b,vcov=C,estimable=estimable,terms=terms,X=X,Z=Z,rank=ncol(A))
}
effects<-function(f,Z,pairs){
 out<-lapply(seq_len(nrow(pairs)),function(i){
  a<-pairs$treatment[i];b<-pairs$reference[i];ok<-a%in%rownames(Z)&&b%in%rownames(Z)
  if(ok){z<-Z[a,]-Z[b,];est<-sum(z*f$beta);se<-sqrt(max(0,drop(z%*%f$vcov%*%z)));lo<-est-q95*se;hi<-est+q95*se;p<-2*pnorm(-abs(est/se))}else{est<-se<-lo<-hi<-p<-NA_real_}
  data.frame(treatment=a,reference=b,primary=b=='NAT',estimate=est,se=se,lower=lo,upper=hi,p=p,estimable=ok,status=if(ok)'Estimated'else'Node absent or disconnected')
 });do.call(rbind,out)
}
rank_scores<-function(f,Z){
 tt<-rownames(Z);p<-sapply(tt,function(a)mean(vapply(setdiff(tt,a),function(b){z<-Z[a,]-Z[b,];pnorm(sum(z*f$beta)/sqrt(max(0,drop(z%*%f$vcov%*%z))))},numeric(1))))
 data.frame(treatment=tt,Pscore=unname(p),rank=rank(-p,ties.method='average'))
}
direction<-function(x)ifelse(is.na(x),NA_character_,ifelse(abs(x)<1e-10,'Zero',ifelse(x>0,'Positive','Negative')))
conclusion<-function(lo,hi)ifelse(is.na(lo),NA_character_,ifelse(lo>0,'Benefit',ifelse(hi<0,'Harm','Not significant')))

all_runs<-all_eff<-all_rank<-all_coef<-all_input<-list();objects<-list();checks<-list()
for(k in seq_len(nrow(spec))){
 key<-spec$key[k];obj<-selected[[key]];main_d<-obj$input$d;tt0<-sort(unique(main_d$t));terms<-strsplit(spec$terms[k],';',fixed=TRUE)[[1]]
 pairs<-do.call(rbind,lapply(combn(tt0,2,simplify=FALSE),function(ab){if('NAT'%in%ab)data.frame(treatment=setdiff(ab,'NAT'),reference='NAT')else data.frame(treatment=ab[1],reference=ab[2])}))
 rr<-ratings[ratings$outcome==key,];high<-intersect(unique(main_d$id),rr$id[rr$Overall=='High risk']);unknown<-setdiff(unique(main_d$id),rr$id)
 input_risk<-data.frame(outcome=key,id=unique(main_d$id));m<-match(input_risk$id,rr$id)
 input_risk$author<-rr$author[m];input_risk$pdf_label<-rr$pdf_label[m];input_risk$Overall<-rr$Overall[m];input_risk$Original_Overall<-rr$Original_Overall[m]
 input_risk$Overall[is.na(input_risk$Overall)]<-'Not present in provided RoB2 PDF';input_risk$author[is.na(input_risk$author)]<-'Combined report ID (Reyes / Foreman 2024)';input_risk$pdf_label[is.na(input_risk$pdf_label)]<-'Reyes / Foreman (2024)'
 input_risk$arms<-vapply(input_risk$id,function(id)sum(main_d$id==id),integer(1));input_risk$n<-vapply(input_risk$id,function(id)sum(main_d$n[main_d$id==id]),numeric(1));all_input[[key]]<-input_risk
 scenarios<-list(list(id='Primary',type='Primary',remove=character()),list(id='Exclude_High',type='Exclude_High',remove=high))
 for(id in unique(main_d$id))scenarios[[length(scenarios)+1L]]<-list(id=paste0('LOO_',id),type='LOO',remove=id)
 if(length(unknown))scenarios[[length(scenarios)+1L]]<-list(id='Exclude_High_Unrated',type='Supplement_unrated',remove=union(high,unknown))
 for(g in list(c('RPT023','RPT184'),c('RPT092_306','RPT141'))){
  gg<-intersect(unique(main_d$id),g)
  if(length(gg)>1)scenarios[[length(scenarios)+1L]]<-list(id=paste0('Joint_',paste(gg,collapse='_')),type='Supplement_linked',remove=gg)
 }
 main_fit<-main_inp<-main_eff<-main_coef<-main_rank<-NULL
 for(s in scenarios){
  raw<-main_d[!main_d$id%in%s$remove,];inp<-tryCatch(input(raw,spec$sm[k]),error=function(e)e)
  status<-'Estimated';fit<-cf<-NULL
  if(inherits(inp,'error'))status<-conditionMessage(inp) else{
   fit<-tryCatch(fit_gls(inp),error=function(e)e);if(inherits(fit,'error')){status<-conditionMessage(fit);fit<-NULL}
  }
  if(is.null(fit))stop(paste('Cannot estimate sensitivity model:',key,s$id,status))
  ef<-effects(fit,inp$Z,pairs);rk<-rank_scores(fit,inp$Z)
  cf<-tryCatch(component_fit(inp,terms),error=function(e)e)
  if(inherits(cf,'error')){
   co<-data.frame(term=terms,estimate=NA_real_,se=NA_real_,lower=NA_real_,upper=NA_real_,p=NA_real_,estimable=FALSE,status=conditionMessage(cf));ctau<-NA_real_;crank<-0L
  }else{
   se<-sqrt(pmax(0,diag(cf$vcov)));est<-cf$beta;est[!cf$estimable]<-NA_real_;se[!cf$estimable]<-NA_real_
   co<-data.frame(term=terms,estimate=unname(est),se=unname(se),lower=unname(est-q95*se),upper=unname(est+q95*se),p=unname(2*pnorm(-abs(est/se))),estimable=cf$estimable,status=ifelse(cf$estimable,'Estimated','Coefficient not separately identifiable'))
   ctau<-cf$fit$tau2;crank<-cf$rank
  }
  if(s$type=='Primary'){
   main_fit<-fit;main_inp<-inp;main_eff<-ef;main_rank<-rk;main_coef<-co
   stopifnot(max(abs(fit$beta-obj$fit$beta))<1e-7,abs(fit$tau2-obj$fit$tau2)<1e-7)
  }
  for(tab in c('ef','co')){
   x<-get(tab);b<-if(tab=='ef')main_eff else main_coef
   x$baseline_estimate<-b$estimate;x$baseline_se<-b$se;x$baseline_lower<-b$lower;x$baseline_upper<-b$upper
   x$delta<-x$estimate-b$estimate;x$delta_over_main_se<-x$delta/b$se
   x$baseline_direction<-direction(b$estimate);x$direction<-direction(x$estimate)
   x$baseline_conclusion<-conclusion(b$lower,b$upper);x$conclusion<-conclusion(x$lower,x$upper)
   x$point_reversal<-ifelse(x$estimable,x$direction!=x$baseline_direction,NA)
   x$significance_changed<-ifelse(x$estimable,(x$conclusion=='Not significant')!=(x$baseline_conclusion=='Not significant'),NA)
   x$conclusion_changed<-ifelse(x$estimable,x$conclusion!=x$baseline_conclusion,NA)
   x$significant_reversal<-ifelse(x$estimable,(x$conclusion=='Benefit'&x$baseline_conclusion=='Harm')|(x$conclusion=='Harm'&x$baseline_conclusion=='Benefit'),NA)
   x$outcome<-key;x$scenario<-s$id;x$type<-s$type;x$removed<-paste(s$remove,collapse=';');x$tau2<-if(tab=='ef')fit$tau2 else ctau
   assign(tab,x)
  }
  common<-main_inp$Z[match(rk$treatment,rownames(main_inp$Z)),,drop=FALSE];bc<-rank_scores(main_fit,common)
  rfull<-data.frame(treatment=tt0);i<-match(tt0,rk$treatment);j<-match(tt0,main_rank$treatment);z<-match(tt0,bc$treatment)
  rfull$Pscore<-rk$Pscore[i];rfull$rank<-rk$rank[i];rfull$baseline_Pscore_full<-main_rank$Pscore[j];rfull$baseline_rank_full<-main_rank$rank[j]
  rfull$baseline_Pscore_common<-bc$Pscore[z];rfull$baseline_rank_common<-bc$rank[z];rfull$delta_Pscore_common<-rfull$Pscore-rfull$baseline_Pscore_common;rfull$delta_rank_common<-rfull$rank-rfull$baseline_rank_common
  rfull$estimable<-!is.na(rfull$Pscore);rfull$outcome<-key;rfull$scenario<-s$id;rfull$type<-s$type;rfull$removed<-paste(s$remove,collapse=';')
  top<-paste(rk$treatment[rk$rank==min(rk$rank)],collapse=';');base_top<-paste(main_rank$treatment[main_rank$rank==min(main_rank$rank)],collapse=';');base_common_top<-paste(bc$treatment[bc$rank==min(bc$rank)],collapse=';')
  top_retained<-all(strsplit(base_top,';',fixed=TRUE)[[1]]%in%rk$treatment)
  stats<-data.frame(outcome=key,scenario=s$id,type=s$type,removed=paste(s$remove,collapse=';'),k_removed=length(s$remove),k_remaining=length(unique(inp$d$id)),n_remaining=sum(inp$d$n),nodes=length(inp$tt),lost_nodes=paste(setdiff(tt0,inp$tt),collapse=';'),disconnected_ids=paste(inp$excluded_ids,collapse=';'),tau2=fit$tau2,tau2_primary=main_fit$tau2,CNMA_tau2=ctau,CNMA_rank=crank,CNMA_terms=length(terms),top=top,primary_top=base_top,primary_top_retained=top_retained,primary_top_changed=if(top_retained)base_top!=top else NA,baseline_common_top=base_common_top,top_changed_common=base_common_top!=top,spearman_common=cor(bc$rank[match(rk$treatment,bc$treatment)],rk$rank,method='spearman'),max_abs_rank_change_common=max(abs(rfull$delta_rank_common),na.rm=TRUE),max_abs_Pscore_change_common=max(abs(rfull$delta_Pscore_common),na.rm=TRUE),primary_estimable=sum(ef$primary&ef$estimable),primary_unavailable=sum(ef$primary&!ef$estimable),primary_point_reversal=sum(ef$primary&ef$point_reversal,na.rm=TRUE),primary_sig_changed=sum(ef$primary&ef$significance_changed,na.rm=TRUE),primary_sig_reversal=sum(ef$primary&ef$significant_reversal,na.rm=TRUE),all_point_reversal=sum(ef$point_reversal,na.rm=TRUE),all_sig_changed=sum(ef$significance_changed,na.rm=TRUE),CNMA_unavailable=sum(!co$estimable),CNMA_point_reversal=sum(co$point_reversal,na.rm=TRUE),CNMA_sig_changed=sum(co$significance_changed,na.rm=TRUE),status=status)
  id<-paste(key,s$id,sep='|');all_runs[[id]]<-stats;all_eff[[id]]<-ef;all_rank[[id]]<-rfull;all_coef[[id]]<-co
  objects[[id]]<-list(input=inp,fit=fit,component=if(inherits(cf,'error'))NULL else cf,scenario=s)
 }
 cat(key,':',length(scenarios),'fits; high-risk removed',paste(high,collapse=','),'; primary top',main_rank$treatment[which.min(main_rank$rank)],'\n');flush.console()
}
for(item in list(c('scenarios','all_runs'),c('effects','all_eff'),c('rankings','all_rank'),c('components','all_coef'),c('risk_input','all_input'))){
 write.csv(do.call(rbind,get(item[2])),file.path(base,paste0(item[1],'.csv')),row.names=FALSE,na='')
}
saveRDS(objects,file.path(base,'sensitivity_models.rds'))
# Independent package verification of all primary/exclude-high and one LOO per outcome.
suppressPackageStartupMessages(library(netmeta));verification<-list()
for(key in spec$key){
 ids<-c(paste(key,c('Primary','Exclude_High'),sep='|'),grep(paste0('^',key,'\\|LOO_'),names(objects),value=TRUE)[1])
 for(id in ids){
  o<-objects[[id]];inp<-o$input;ct<-inp$ct;pairs<-list()
  for(study in unique(ct$id)){
   ii<-which(ct$id==study);tt<-c(ct$baseline[ii[1]],ct$treatment[ii]);y<-c(0,ct$yi[ii]);V<-matrix(0,length(tt),length(tt));V[-1,-1]<-inp$V[ii,ii,drop=FALSE]
   for(ab in combn(seq_along(tt),2,simplify=FALSE)){a<-ab[2];b<-ab[1];pairs[[length(pairs)+1L]]<-data.frame(studlab=study,treat1=tt[a],treat2=tt[b],TE=y[a]-y[b],seTE=sqrt(V[a,a]+V[b,b]-2*V[a,b]),correlated=inp$sm=='SMD')}
  }
  pd<-do.call(rbind,pairs);net<-netmeta(TE,seTE,treat1,treat2,studlab,data=pd,correlated=correlated,sm=inp$sm,common=FALSE,random=TRUE,tau.preset=sqrt(o$fit$tau2),reference.group='NAT')
  tt<-setdiff(inp$tt,'NAT');delta<-max(abs(net$TE.random[tt,'NAT']-o$fit$beta[tt]),abs(net$seTE.random[tt,'NAT']-sqrt(diag(o$fit$vcov))[tt]))
  ps<-netrank(net,small.values='bad')$ranking.random;calc<-rank_scores(o$fit,inp$Z);dp<-max(abs(ps[calc$treatment]-calc$Pscore))
  stopifnot(delta<1e-6,dp<1e-6);verification[[id]]<-data.frame(model=id,effect_se_error=delta,Pscore_error=dp)
 }
}
write.csv(do.call(rbind,verification),file.path(base,'independent_verification.csv'),row.names=FALSE)
ref<-read.csv(file.path(base,'cnma_reference.csv'));out<-do.call(rbind,all_coef);out<-out[out$type=='Primary',];err<-numeric()
for(i in seq_len(nrow(ref))){a<-out[out$outcome==ref$outcome[i]&out$term==ref$term[i],c('estimate','lower','upper')];stopifnot(nrow(a)==1);err<-c(err,max(abs(as.numeric(a)-as.numeric(ref[i,c('estimate','lower','upper')]))))}
stopifnot(max(err)<.00501);write.csv(data.frame(reference_rows=nrow(ref),max_difference=max(err)),file.path(base,'CNMA_pdf_validation.csv'),row.names=FALSE)
writeLines(capture.output(sessionInfo()),file.path(base,'sessionInfo.txt'))
cat('Validated CNMA coefficients:',nrow(ref),'max display error',max(err),'Total models:',length(objects),'\n')
