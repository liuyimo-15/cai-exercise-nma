options(stringsAsFactors=FALSE)
lines<-readLines('work/sensitivity_20261112/analyse.R');end<-grep('^all_runs<-',lines)[1]-1L
eval(parse(text=paste(lines[seq_len(end)],collapse='\n')))
b<-'work/subgroups_20261113';rawall<-read.csv(file.path(b,'arms.csv'),check.names=FALSE)
summary<-efout<-intout<-adultout<-verifications<-list();models<-list();idx<-0L
bdiag_rect<-function(xs){
 out<-matrix(0,sum(vapply(xs,nrow,integer(1))),sum(vapply(xs,ncol,integer(1))));r<-c<-0L
 for(x in xs){out[r+seq_len(nrow(x)),c+seq_len(ncol(x))]<-x;r<-r+nrow(x);c<-c+ncol(x)};out
}
# Joint subgroup model: subgroup-specific treatment effects, common REML tau^2,
# exact existing within-study V blocks, and multi-arm random-effect correlation .5.
runjoint<-function(aa,key,var,groups,scenario,removed){
 dd<-aa[aa[[var]]%in%groups,];gg<-lapply(groups,function(g)dd[dd[[var]]==g,]);names(gg)<-groups
 ins<-lapply(gg,function(x)tryCatch(input(x,spec$sm[match(key,spec$key)]),error=function(e)e))
 ks<-sapply(ins,function(x)if(inherits(x,'error'))0L else length(unique(x$d$id)))
 ns<-sapply(gg,function(x)length(unique(x$id)))
 common<-if(all(!sapply(ins,inherits,'error')))Reduce(intersect,lapply(ins,function(x)x$tt))else character()
 meta<-data.frame(outcome=key,variable=var,scenario=scenario,group0=groups[1],group1=groups[2],k0_available=ns[1],k1_available=ns[2],k0=ks[1],k1=ks[2],common_nodes=paste(common,collapse=';'),removed=paste(removed,collapse=';'),missing_ids=paste(setdiff(unique(aa$id),unique(dd$id)),collapse=';'),disconnected_ids=paste(unique(unlist(lapply(ins,function(x)if(inherits(x,'error'))character()else x$excluded_ids))),collapse=';'),tau2=NA_real_,wald=NA_real_,df=NA_integer_,p_interaction=NA_real_,residual_df=NA_integer_,status='')
 if(any(ks<5L)){meta$status<-'Not fitted: fewer than 5 NAT-connected report IDs in a subgroup';return(list(summary=meta))}
 tt<-setdiff(common,'NAT');if(!length(tt)){meta$status<-'Not fitted: no shared active treatment connected to NAT';return(list(summary=meta))}
 ji<-list(ct=do.call(rbind,lapply(ins,function(x)x$ct)),V=bdiag_base(lapply(ins,function(x)x$V)),R=bdiag_base(lapply(ins,function(x)x$R)),X=bdiag_rect(lapply(ins,function(x)x$X)))
 colnames(ji$X)<-unlist(lapply(seq_along(groups),function(i)paste(groups[i],colnames(ins[[i]]$X),sep='|')))
 meta$residual_df<-nrow(ji$X)-ncol(ji$X)
 if(meta$residual_df<5L){meta$status<-'Not fitted: fewer than 5 residual contrast degrees of freedom';return(list(summary=meta))}
 f<-tryCatch(fit_gls(ji),error=function(e)e);if(inherits(f,'error')){meta$status<-conditionMessage(f);return(list(summary=meta))}
 off<-c(0L,ncol(ins[[1]]$X));C<-matrix(0,length(tt),ncol(ji$X));rownames(C)<-tt
 for(j in seq_along(tt)){C[j,match(paste(groups[1],tt[j],sep='|'),colnames(ji$X))]<--1;C[j,match(paste(groups[2],tt[j],sep='|'),colnames(ji$X))]<-1}
 dif<-drop(C%*%f$beta);cv<-C%*%f$vcov%*%t(C);se<-sqrt(diag(cv));Q<-drop(crossprod(dif,solve(cv,dif)))
 meta$tau2<-f$tau2;meta$wald<-Q;meta$df<-nrow(C);meta$p_interaction<-pchisq(Q,nrow(C),lower.tail=FALSE);meta$status<-'Estimated: common tau2; shared-node Wald interaction'
 pair<-data.frame(treatment=tt,reference='NAT');efs<-list();ints<-data.frame(outcome=key,variable=var,scenario=scenario,treatment=tt,reference='NAT',difference_group1_minus0=dif,se=se,lower=dif-q95*se,upper=dif+q95*se,p=2*pnorm(-abs(dif/se)),group0=groups[1],group1=groups[2])
 for(i in seq_along(groups)){
  pos<-off[i]+seq_len(ncol(ins[[i]]$X));ff<-f;ff$beta<-f$beta[pos];names(ff$beta)<-colnames(ins[[i]]$X);ff$vcov<-f$vcov[pos,pos,drop=FALSE]
  pairs<-data.frame(treatment=setdiff(ins[[i]]$tt,'NAT'),reference='NAT');ef<-effects(ff,ins[[i]]$Z,pairs);ef$outcome<-key;ef$variable<-var;ef$scenario<-scenario;ef$group<-groups[i];ef$k<-ks[i];ef$n<-sum(ins[[i]]$d$n);ef$common_node<-ef$treatment%in%tt;ef$tau2<-f$tau2
  rk<-rank_scores(ff,ins[[i]]$Z[common,,drop=FALSE]);ef$Pscore_common<-rk$Pscore[match(ef$treatment,rk$treatment)];ef$rank_common<-rk$rank[match(ef$treatment,rk$treatment)];ef$report_ids<-vapply(ef$treatment,function(t)paste(unique(ins[[i]]$d$id[ins[[i]]$d$t==t]),collapse=';'),character(1));ef$node_k<-vapply(ef$treatment,function(t)length(unique(ins[[i]]$d$id[ins[[i]]$d$t==t])),integer(1));efs[[i]]<-ef
 }
 list(summary=meta,effects=do.call(rbind,efs),interaction=ints,model=list(input=ji,fit=f,ins=ins,groups=groups,common=common))
}
for(key in spec$key){
 aa<-rawall[rawall$outcome==key,];baseobj<-selected[[key]];inp<-input(aa,spec$sm[match(key,spec$key)]);f<-fit_gls(inp)
 stopifnot(max(abs(f$beta-baseobj$fit$beta))<1e-7,abs(f$tau2-baseobj$fit$tau2)<1e-7)
 suspect<-character();ids<-unique(aa$id)
 if(all(c('RPT023','RPT184')%in%ids))suspect<-c(suspect,'RPT184')
 if(all(c('RPT092_306','RPT141')%in%ids))suspect<-c(suspect,'RPT092_306')
 scenarios<-list(R050=character());if(length(suspect))scenarios[['Cohort_cautious']]<-suspect
 for(sc in names(scenarios)){
  rem<-scenarios[[sc]];a<-aa[!aa$id%in%rem,]
  jobs<-list(duration=c('Short_lt6','Long_ge6'),population_group=c('Recreational_active','Athlete'),supervision_group=c('Onsite','Home_remote'))
  if(key%in%c('Anterior','Posteromedial','Posterolateral'))jobs[['tool']]<-c('SEBT','YBT')
  for(v in names(jobs)){
   res<-runjoint(a,key,v,jobs[[v]],sc,rem);id<-paste(key,v,sc,sep='|');summary[[id]]<-res$summary
   if(!is.null(res$effects)){efout[[id]]<-res$effects;intout[[id]]<-res$interaction;models[[id]]<-res$model}
  }
 }
 for(mode in c('Adult_verified','Exclude_known_under18')){
  keep<-if(mode=='Adult_verified')aa$age_scope=='Adult_verified'else aa$age_scope!='Includes_under18';a<-aa[keep,];ip<-tryCatch(input(a,spec$sm[match(key,spec$key)]),error=function(e)e)
  if(inherits(ip,'error'))next
  ff<-tryCatch(fit_gls(ip),error=function(e)e);if(inherits(ff,'error'))next
  pairs<-data.frame(treatment=setdiff(inp$tt,'NAT'),reference='NAT');ef<-effects(ff,ip$Z,pairs);bf<-effects(f,inp$Z,pairs)
  ef$outcome<-key;ef$restriction<-mode;ef$k<-length(unique(ip$d$id));ef$n<-sum(ip$d$n);ef$tau2<-ff$tau2;ef$baseline_estimate<-bf$estimate;ef$baseline_lower<-bf$lower;ef$baseline_upper<-bf$upper
  ef$point_reversal<-ef$estimate*bf$estimate<0;ef$significance_changed<-(ef$lower>0|ef$upper<0)!=(bf$lower>0|bf$upper<0);ef$removed_ids<-paste(setdiff(unique(aa$id),unique(a$id)),collapse=';');ef$disconnected_ids<-paste(ip$excluded_ids,collapse=';');adultout[[paste(key,mode)]]<-ef
 }
 cat(key,': subgroup scenarios processed\n');flush.console()
}
summary<-do.call(rbind,summary);summary$q_BH<-NA_real_
for(sc in unique(summary$scenario)){ii<-which(summary$scenario==sc&!is.na(summary$p_interaction));summary$q_BH[ii]<-p.adjust(summary$p_interaction[ii],method='BH')}
for(item in list(c('summary','summary'),c('effects','efout'),c('interactions','intout'),c('adult_restrictions','adultout'))){
 x<-get(item[2]);if(is.list(x)&&!is.data.frame(x))x<-do.call(rbind,x);write.csv(x,file.path(b,paste0(item[1],'.csv')),row.names=FALSE,na='')
}
saveRDS(models,file.path(b,'subgroup_models.rds'))
# Independent verification of every subgroup model effect/SE using netmeta at the same tau.
suppressPackageStartupMessages(library(netmeta))
for(id in names(models)){
 mo<-models[[id]]
 for(i in seq_along(mo$ins)){
  ip<-mo$ins[[i]];pd<-list()
  for(st in unique(ip$ct$id)){
   ii<-which(ip$ct$id==st);tt<-c(ip$ct$baseline[ii[1]],ip$ct$treatment[ii]);y<-c(0,ip$ct$yi[ii]);V<-matrix(0,length(tt),length(tt));V[-1,-1]<-ip$V[ii,ii,drop=FALSE]
   for(ab in combn(seq_along(tt),2,simplify=FALSE)){a<-ab[2];c<-ab[1];pd[[length(pd)+1L]]<-data.frame(studlab=st,treat1=tt[a],treat2=tt[c],TE=y[a]-y[c],seTE=sqrt(V[a,a]+V[c,c]-2*V[a,c]),correlated=ip$sm=='SMD')}
  }
  dat<-do.call(rbind,pd);nf<-netmeta(TE,seTE,treat1,treat2,studlab,data=dat,correlated=correlated,sm=ip$sm,common=FALSE,random=TRUE,tau.preset=sqrt(mo$fit$tau2),reference.group='NAT')
  ef<-efout[[id]];ef<-ef[ef$group==mo$groups[i],];tt<-ef$treatment;er<-max(abs(ef$estimate-nf$TE.random[tt,'NAT']),abs(ef$se-nf$seTE.random[tt,'NAT']))
  stopifnot(er<1e-6);verifications[[paste(id,i)]]<-data.frame(model=id,group=mo$groups[i],max_error=er)
 }
}
write.csv(do.call(rbind,verifications),file.path(b,'independent_verification.csv'),row.names=FALSE)
writeLines(capture.output(sessionInfo()),file.path(b,'sessionInfo.txt'))
cat('Fitted',sum(!is.na(summary$p_interaction)),'joint subgroup models; checked',length(verifications),'strata\n')
