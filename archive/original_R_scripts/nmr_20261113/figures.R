options(stringsAsFactors=FALSE)
b<-'work/nmr_20261113';o<-'outputs/2026-11-13_网状Meta回归';models<-readRDS(file.path(b,'models.rds'));su<-read.csv(file.path(b,'summary.csv'))
display<-c(CAIT='CAIT (MD)',FAAM_ADL='FAAM-ADL (MD)',Sports='Sports function (Hedges g)',Anterior='Anterior balance (Hedges g)',Posteromedial='Posteromedial balance (Hedges g)',Posterolateral='Posterolateral balance (Hedges g)')
blue<-'#2D6484';gray<-'#52606D';navy<-'#243B53'
setup<-function(){par(mfrow=c(3,2),oma=c(1,1,1,1),family='Arial',ps=10,cex=1,mgp=c(1.9,.5,0),tcl=-.2,xpd=FALSE)}
check<-function(l,G)sum((l-drop(G%*%crossprod(G,l)))^2)<1e-10*max(1,sum(l^2))
ci<-function(l,m){
 if(!check(l,m$G))return(rep(NA_real_,3));e<-sum(l*m$beta);se<-sqrt(max(0,drop(l%*%m$vcov%*%l)));c(e,e-qt(.975,m$df)*se,e+qt(.975,m$df)*se)
}
forest<-function(v){
 setup()
 for(key in names(display)){
  id<-paste(key,v,'Treatment_specific','R050',sep='|');m<-models[[id]];s<-su[su$model_id==id,];tt<-colnames(m$input$X);a<-matrix(NA_real_,length(tt),3)
  for(j in seq_along(tt)){l<-rep(0,2*m$p);l[m$p+j]<-1;a[j,]<-ci(l,m)}
  keep<-is.finite(a[,1]);a<-a[keep,,drop=FALSE];tt<-tt[keep];ys<-rev(seq_along(tt));xr<-range(c(0,a));xr<-xr+c(-.08,.08)*max(.1,diff(xr))
  par(mar=c(3.5,7,3.4,1.2));plot(NA,xlim=xr,ylim=c(.4,length(tt)+.7),xlab=if(v=='weeks')'Beta per 2-week increase (95% CI)'else'Beta per 5-point lower baseline CAIT (95% CI)',ylab='',yaxt='n',bty='n',cex.lab=.82,cex.axis=.8)
  segments(xr[1],ys,xr[2],ys,col='#E9EDF1',lwd=.6);abline(v=0,lty=2,col='#8A929C',lwd=.8);axis(2,at=ys,labels=tt,las=1,tick=FALSE,cex.axis=.8);segments(a[,2],ys,a[,3],ys,col=blue,lwd=1.25);points(a[,1],ys,pch=16,col=blue,cex=.9)
  mtext(display[key],side=3,line=2.0,adj=0,font=2,cex=.98,col=navy);mtext(sprintf('Global P=%.4f; BH q=%.4f; k=%d',s$p_global,s$q_BH,s$k_connected),side=3,line=.65,adj=0,cex=.78,col=gray)
 }
}
bubbles<-list();curves<-list()
curveplot<-function(v){
 setup()
 for(key in names(display)){
  id<-paste(key,v,'Treatment_specific','R050',sep='|');m<-models[[id]];ip<-m$input;j<-match('BAL',colnames(ip$X));balids<-unique(ip$d$id[ip$d$t=='BAL']);bx<-m$cov[[v]][m$cov$id%in%balids];xx<-seq(min(bx),max(bx),length.out=101)
  vals<-t(vapply(xx,function(x){l<-rep(0,2*m$p);l[j]<-1;l[m$p+j]<-(x-m$x0)/m$v$increment;ci(l,m)},numeric(3)))
  pts<-list()
  for(st in unique(ip$d$id)){
   a<-ip$d[ip$d$id==st,];if(!all(c('BAL','NAT')%in%a$t))next;ii<-which(ip$ct$id==st);z<-as.numeric(ip$ct$treatment[ii]=='BAL')-as.numeric(ip$ct$treatment[ii]=='NAT');yi<-sum(z*ip$ct$yi[ii]);va<-drop(z%*%ip$V[ii,ii,drop=FALSE]%*%z)
   pts[[st]]<-data.frame(outcome=key,variable=v,id=st,x=m$cov[[v]][match(st,m$cov$id)],effect=yi,variance=va)
  }
  p<-if(length(pts))do.call(rbind,pts)else data.frame(x=numeric(),effect=numeric(),variance=numeric());yr<-range(c(0,vals,p$effect),na.rm=TRUE);yr<-yr+c(-.08,.12)*max(.1,diff(yr))
  par(mar=c(3.5,5.3,2.3,1.3));plot(NA,xlim=range(xx),ylim=yr,xlab=if(v=='weeks')'Treatment duration (weeks)'else'Baseline adult CAIT (points)',ylab=if(key%in%c('CAIT','FAAM_ADL'))'MD vs NAT'else'Hedges g vs NAT',bty='n',cex.lab=.84,cex.axis=.8)
  abline(h=0,lty=2,col='#A7B1BC',lwd=.7);polygon(c(xx,rev(xx)),c(vals[,2],rev(vals[,3])),col='#DDE8EE',border=NA);lines(xx,vals[,1],col=blue,lwd=1.6)
  if(nrow(p)){w<-1/p$variance;points(p$x,p$effect,pch=21,bg=adjustcolor(blue,.65),col=blue,cex=2.2*sqrt(w/max(w)))}
  mtext(paste(display[key],': BAL vs NAT',sep=''),side=3,line=.85,adj=0,font=2,cex=.94,col=navy)
  legend('topright',legend=c('Network fit','Direct comparison'),col=blue,lty=c(1,NA),pch=c(NA,21),pt.bg=blue,bty='n',cex=.69)
  bubbles[[paste(key,v)]]<<-p;curves[[paste(key,v)]]<<-data.frame(outcome=key,variable=v,x=xx,estimate=vals[,1],lower=vals[,2],upper=vals[,3])
 }
}
quartz(type='pdf',file=file.path(o,'六结局_Meta回归_核心变量森林图与曲线.pdf'),width=12,height=11.8,family='Arial',pointsize=10)
forest('weeks');forest('cait');curveplot('weeks');curveplot('cait');dev.off()
for(v in c('weeks','cait'))for(type in c('forest','curve')){
 png(file.path(b,paste0(type,'_',v,'.png')),width=12,height=11.8,units='in',res=180,type='quartz',pointsize=10)
 if(type=='forest')forest(v)else curveplot(v);dev.off()
}
write.csv(do.call(rbind,bubbles),file.path(b,'bubble_points.csv'),row.names=FALSE,na='');write.csv(do.call(rbind,curves),file.path(b,'curve_values.csv'),row.names=FALSE,na='')
cat('4-page vector PDF: core-variable slopes and BAL network regression curves\n')
