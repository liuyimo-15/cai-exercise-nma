options(stringsAsFactors=FALSE);base<-'work/cinema_20261112'
.libPaths(c(normalizePath(file.path(base,'Rlib')),.libPaths()))
suppressPackageStartupMessages(library(netmeta))
fits<-readRDS(file.path(base,'candidate_fits.rds'));balanced<-readRDS(file.path(base,'matched_balance.rds'))
selected<-c(list(CAIT=fits[['CAIT_raw_Cohen_N']],FAAM_ADL=fits[['FAAM_ADL_raw_Cohen_N']],Sports=fits[['Sports_raw_Hedges_N']]),balanced)
outdir<-file.path(base,'diagnostics');dir.create(outdir,showWarnings=FALSE)
summary<-all_effects<-list();objects<-list()
for(key in names(selected)){
 obj<-selected[[key]];inp<-obj$input;fit<-obj$fit;ct<-inp$ct
 od<-file.path(outdir,key);dir.create(od,showWarnings=FALSE)
 pairs<-list()
 for(id in unique(ct$id)){
  ii<-which(ct$id==id);tt<-c(ct$baseline[ii[1]],ct$treatment[ii]);y<-c(0,ct$yi[ii]);V<-matrix(0,length(tt),length(tt));V[-1,-1]<-inp$V[ii,ii,drop=FALSE]
  for(ab in combn(seq_along(tt),2,simplify=FALSE)){
   a<-ab[2];b<-ab[1]
   pairs[[length(pairs)+1]]<-data.frame(studlab=id,treat1=tt[a],treat2=tt[b],TE=y[a]-y[b],seTE=sqrt(V[a,a]+V[b,b]-2*V[a,b]),n1=inp$d$n[match(paste(id,tt[a]),paste(inp$d$id,inp$d$t))],n2=inp$d$n[match(paste(id,tt[b]),paste(inp$d$id,inp$d$t))],correlated=inp$sm=='SMD')
  }
 }
 pd<-do.call(rbind,pairs)
 net<-netmeta(TE,seTE,treat1,treat2,studlab,data=pd,correlated=correlated,n1=n1,n2=n2,sm=inp$sm,common=FALSE,random=TRUE,tau.preset=sqrt(fit$tau2),reference.group='NAT',prediction=TRUE,method.tau='REML')
 rg<-obj$pred;match_delta<-max(abs(net$TE.random[rg$treatment,'NAT']-rg$estimate),abs(net$seTE.random[rg$treatment,'NAT']-(rg$upper-rg$lower)/(2*qnorm(.975))))
 stopifnot(match_delta<1e-6)
 split<-netsplit(net,method='Back-calculation',common=FALSE,random=TRUE)
 contrib<-netcontrib(net,method='shortestpath',common=FALSE,random=TRUE,study=TRUE)
 write.csv(pd,file.path(od,'pairwise_input.csv'),row.names=FALSE)
 write.csv(inp$d,file.path(od,'analysed_arms.csv'),row.names=FALSE)
 for(name in c('direct.random','indirect.random','compare.random'))write.csv(as.data.frame(split[[name]]),file.path(od,paste0(name,'.csv')),row.names=FALSE)
 write.csv(as.data.frame(contrib$study.random),file.path(od,'study_contributions.csv'),row.names=FALSE)
 write.csv(as.data.frame(contrib$random),file.path(od,'direct_contributions.csv'),row.names=TRUE)
 # Global design-by-treatment Wald with the same primary REML heterogeneity.
 design<-character(nrow(ct))
 for(id in unique(ct$id)){ii<-which(ct$id==id);design[ii]<-paste(sort(unique(c(ct$treatment[ii],ct$baseline[ii]))),collapse='|')}
 cols<-unique(paste(design,ct$treatment,sep='::'));Xd<-matrix(0,nrow(ct),length(cols));for(i in seq_len(nrow(ct)))Xd[i,match(paste(design[i],ct$treatment[i],sep='::'),cols)]<-1
 residQ<-function(X){b<-solve(crossprod(X,fit$W%*%X),crossprod(X,fit$W%*%ct$yi));e<-ct$yi-X%*%b;drop(crossprod(e,fit$W%*%e))}
 df_inc<-qr(Xd)$rank-ncol(inp$X);Q_inc<-max(0,residQ(inp$X)-residQ(Xd));p_inc<-if(df_inc>0)pchisq(Q_inc,df_inc,lower.tail=FALSE)else NA_real_
 global<-data.frame(outcome=key,k=length(unique(inp$d$id)),nodes=length(inp$tt),arms=nrow(inp$d),participants=sum(inp$d$n),tau2=fit$tau2,tau=sqrt(fit$tau2),Q_common=obj$common$Q,df_common=obj$common$n-obj$common$p,I2=100*max(0,1-(obj$common$n-obj$common$p)/obj$common$Q),global_Q=Q_inc,global_df=df_inc,global_p=p_inc,netmeta_agreement=match_delta)
 summary[[key]]<-global
 nms<-colnames(inp$X);Z<-matrix(0,length(inp$tt),length(nms),dimnames=list(inp$tt,nms));for(n in nms)Z[n,n]<-1
 effects<-list()
 for(ab in combn(inp$tt,2,simplify=FALSE)){
  # Positive always means the first intervention improves the outcome more.
  if('NAT'%in%ab){a<-setdiff(ab,'NAT');b<-'NAT'}else{a<-ab[1];b<-ab[2]}
  z<-Z[a,]-Z[b,];est<-sum(z*fit$beta);se<-sqrt(drop(z%*%fit$vcov%*%z));half<-qnorm(.975)*sqrt(se^2+fit$tau2)
  effects[[length(effects)+1]]<-data.frame(outcome=key,treatment=a,reference=b,sm=inp$sm,estimate=est,se=se,lower=est-qnorm(.975)*se,upper=est+qnorm(.975)*se,pi_lower=est-half,pi_upper=est+half,primary=b=='NAT')
 }
 ef<-do.call(rbind,effects);all_effects[[key]]<-ef;write.csv(ef,file.path(od,'all_estimates.csv'),row.names=FALSE)
 objects[[key]]<-list(primary=obj,net=net,split=split,contrib=contrib,global=global)
 cat(key,'k',global$k,'tau2',global$tau2,'global P',p_inc,'agreement',match_delta,'\n')
 if(key=='CAIT'){cat('CONTRIB COLUMNS\n');print(names(contrib$study.random));print(head(contrib$study.random,3));cat('SIDE\n');print(head(split$compare.random,2));print(head(split$direct.random,2));print(head(split$indirect.random,2))}
}
write.csv(do.call(rbind,summary),file.path(base,'network_summary.csv'),row.names=FALSE)
write.csv(do.call(rbind,all_effects),file.path(base,'all_estimates.csv'),row.names=FALSE)
saveRDS(objects,file.path(base,'diagnostic_objects.rds'));writeLines(capture.output(sessionInfo()),file.path(base,'analysis_sessionInfo.txt'))
