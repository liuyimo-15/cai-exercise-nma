quiet_reproduce<-TRUE;source('work/cinema_20261112/reproduce.R')
out<-list();matched<-list()
for(key in c('Anterior','Posteromedial','Posterolateral')){
 raw<-read.csv(file.path(base,paste0(key,'_arms.csv')));raw<-raw[raw$id!='RPT128',]
 connected<-make_input(raw,'SMD','Hedges','N')$d
 drops<-c('NONE',unique(connected$id))
 if(key=='Posteromedial')drops<-c(drops,paste('RPT289',setdiff(unique(connected$id),'RPT289'),sep=';'))
 for(dropid in drops){
  d<-if(dropid=='NONE')connected else connected[!connected$id%in%strsplit(dropid,';',fixed=TRUE)[[1]],]
  inp<-make_input(d,'SMD','Hedges','N');fit<-fit_gaussian(inp);cm<-fit_gaussian(inp,'COMMON')
  ref<-pdfref[pdfref$outcome==key,];if(any(!ref$treatment%in%names(fit$beta)))next
  b<-fit$beta[ref$treatment];se<-sqrt(diag(fit$vcov))[ref$treatment]
  pred<-data.frame(outcome=key,treatment=ref$treatment,estimate=b,lower=b-qnorm(.975)*se,upper=b+qnorm(.975)*se)
  delta<-max(abs(as.matrix(pred[,3:5])-as.matrix(ref[,3:5])))
  out[[length(out)+1]]<-data.frame(outcome=key,drop=dropid,k=length(unique(inp$d$id)),tau2=fit$tau2,Q=cm$Q,df=cm$n-cm$p,max_delta=delta)
  if(delta<=.00501)matched[[key]]<-list(input=inp,fit=fit,common=cm,pred=pred,excluded=c('RPT128','RPT008','RPT284',if(dropid!='NONE')strsplit(dropid,';',fixed=TRUE)[[1]]))
 }
}
results<-do.call(rbind,out);write.csv(results,file.path(base,'balance_subset_matching.csv'),row.names=FALSE);saveRDS(matched,file.path(base,'matched_balance.rds'))
print(do.call(rbind,lapply(split(results,results$outcome),function(x)head(x[order(x$max_delta),],3))),row.names=FALSE)
