options(stringsAsFactors=FALSE)
base <- 'work/cinema_20261112';source("R/bootstrap.R")
spec <- data.frame(key=c('CAIT','FAAM_ADL','Sports','Anterior','Posteromedial','Posterolateral'),sm=c('MD','MD',rep('SMD',4)))
pdfref <- read.csv(file.path(base,'pdf_reference.csv'))
bdiag_base <- function(blocks){
 n<-sum(sapply(blocks,nrow));out<-matrix(0,n,n);i<-1
 for(x in blocks){ii<-i:(i+nrow(x)-1);out[ii,ii]<-x;i<-max(ii)+1};out
}
make_input <- function(d,sm,type='Cohen',den='N'){
 connected<-'NAT';repeat{old<-connected;studies<-unique(d$id[d$t%in%connected]);connected<-unique(d$t[d$id%in%studies]);if(setequal(old,connected))break}
 raw<-d;d<-d[d$t%in%connected,];excluded<-setdiff(unique(raw$id),unique(d$id))
 vb<-rb<-ct<-list();studies<-unique(d$id)
 for(sid in studies){
  a<-d[d$id==sid,];a<-a[order(a$t),];k<-nrow(a);m<-k-1;N<-sum(a$n);df<-N-k
  if(sm=='MD'){
   y<-a$mean[-1]-a$mean[1];va<-a$sd^2/a$n;V<-matrix(va[1],m,m);diag(V)<-va[-1]+va[1]
  }else{
   sdpool<-sqrt(sum((a$n-1)*a$sd^2)/df);J<-if(type=='Hedges')exp(lgamma(df/2)-log(sqrt(df/2))-lgamma((df-1)/2))else 1
   y<-J*(a$mean[-1]-a$mean[1])/sdpool;denom<-if(den=='N')N else df
   V<-matrix(J^2/a$n[1],m,m)+tcrossprod(y)/(2*denom);diag(V)<-J^2*(1/a$n[-1]+1/a$n[1])+y^2/(2*denom)
  }
  R<-matrix(.5,m,m);diag(R)<-1;vb[[sid]]<-V;rb[[sid]]<-R
  ct[[sid]]<-data.frame(id=sid,treatment=a$t[-1],baseline=a$t[1],yi=y)
 }
 ct<-do.call(rbind,ct);tt<-sort(unique(d$t));Z<-diag(length(tt));dimnames(Z)<-list(tt,tt);Z<-Z[,setdiff(tt,'NAT'),drop=FALSE]
 X<-Z[ct$treatment,,drop=FALSE]-Z[ct$baseline,,drop=FALSE]
 list(d=d,ct=ct,V=bdiag_base(vb),R=bdiag_base(rb),X=X,tt=tt,excluded=excluded,sm=sm)
}
fit_gaussian <- function(inp,method='REML'){
 y<-inp$ct$yi;X<-inp$X;V<-inp$V;R<-inp$R;n<-length(y);p<-ncol(X)
 eval<-function(tau2,details=FALSE){
  S<-V+tau2*R;U<-chol(S);W<-chol2inv(U);I<-crossprod(X,W%*%X);C<-solve(I);b<-drop(C%*%crossprod(X,W%*%y));names(b)<-colnames(X);e<-y-drop(X%*%b);Q<-drop(crossprod(e,W%*%e))
  obj<-.5*(2*sum(log(diag(U)))+as.numeric(determinant(I,logarithm=TRUE)$modulus)+Q)
  if(details)list(beta=b,vcov=C,tau2=tau2,Q=Q,n=n,p=p,S=S,W=W)else obj
 }
 if(method=='COMMON')tau2<-0 else{
  upper<-max(1,var(y),max(diag(V)));repeat{opt<-optimize(eval,c(0,upper),tol=1e-10);if(opt$minimum<.95*upper)break;upper<-upper*4};tau2<-if(eval(0)<=opt$objective+1e-10)0 else opt$minimum
 }
 eval(tau2,TRUE)
}
all<-list();fits<-list()
for(i in seq_len(nrow(spec))){
 key<-spec$key[i];d<-read.csv(file.path(base,paste0(key,'_arms.csv')))
 variants<-if(i>=4)c('raw','regular_id','regular_id_no_FAM')else 'raw'
 for(variant in variants){
  dd<-d
  if(variant!='raw')dd<-dd[grepl('^RPT[0-9]{3}$',dd$id),]
  if(variant=='regular_id_no_FAM')dd<-dd[!dd$t%in%c('FAM','BAL+FAM'),]
  dd<-dd[dd$id%in%names(which(table(dd$id)>=2)),]
 for(type in if(spec$sm[i]=='MD')'Cohen' else c('Cohen','Hedges'))for(den in if(spec$sm[i]=='MD')'N' else c('N','df')){
  inp<-make_input(dd,spec$sm[i],type,den);fit<-fit_gaussian(inp);cm<-fit_gaussian(inp,'COMMON')
  ref<-pdfref[pdfref$outcome==key,];b<-fit$beta[ref$treatment];se<-sqrt(diag(fit$vcov))[ref$treatment];pred<-data.frame(outcome=key,treatment=ref$treatment,estimate=b,lower=b-qnorm(.975)*se,upper=b+qnorm(.975)*se)
  delta<-max(abs(as.matrix(pred[,3:5])-as.matrix(ref[,3:5])))
  all[[length(all)+1]]<-data.frame(outcome=key,variant=variant,type=type,den=den,k=length(unique(inp$d$id)),nodes=length(inp$tt),tau2=fit$tau2,Q=cm$Q,df=cm$n-cm$p,max_delta_pdf=delta,excluded=paste(inp$excluded,collapse=';'))
  fits[[paste(key,variant,type,den,sep='_')]]<-list(input=inp,fit=fit,common=cm,pred=pred)
 }
 }
}
result<-do.call(rbind,all);write.csv(result,file.path(base,'reproduction_candidates.csv'),row.names=FALSE);saveRDS(fits,file.path(base,'candidate_fits.rds'));if(!exists('quiet_reproduce'))print(result,row.names=FALSE)
