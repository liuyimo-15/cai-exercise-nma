options(stringsAsFactors=FALSE)
base<-'work/sensitivity_20261112';out<-'results/figures/sensitivity';figdir<-file.path(out,'Figures');dir.create(figdir,recursive=TRUE,showWarnings=FALSE)
source("R/bootstrap.R")
spec<-read.csv(file.path(base,'spec.csv'));sc<-read.csv(file.path(base,'scenarios.csv'));ef<-read.csv(file.path(base,'effects.csv'));rk<-read.csv(file.path(base,'rankings.csv'));co<-read.csv(file.path(base,'components.csv'));risk<-read.csv(file.path(base,'risk_input.csv'))
display<-c(CAIT='CAIT',FAAM_ADL='FAAM-ADL',Sports='Sports function',Anterior='Anterior balance',Posteromedial='Posteromedial balance',Posterolateral='Posterolateral balance')
blue<-'#0072B2';orange<-'#D55E00';teal<-'#009E73';purple<-'#8B5AA8';navy<-'#243B53';gray<-'#6B7280';light<-'#E8EDF2'
safe<-function(x)stringi::stri_trans_general(enc2utf8(x),'Latin-ASCII')
report_label<-function(key,id){a<-risk[risk$outcome==key&risk$id==id,];if(!nrow(a))id else paste(id,safe(a$pdf_label[1]))}
save_fig<-function(stem,fun,width=7.8,height=10.8){
 pdf(file=file.path(figdir,paste0(stem,'.pdf')),width=width,height=height,family='Helvetica',pointsize=9);fun();dev.off()
 png(file.path(base,paste0(stem,'_preview.png')),width=width,height=height,units='in',res=180,type=nma_png_type(),pointsize=9);fun();dev.off()
 cat('Figure',stem,'\n');flush.console()
}
overlay<-function(title,subtitle,legend_labels,legend_cols,legend_pch,note){
 par(fig=c(0,1,0,1),new=TRUE,mar=c(0,0,0,0),oma=c(0,0,0,0),xpd=NA);plot.new();plot.window(xlim=c(0,1),ylim=c(0,1),xaxs='i',yaxs='i')
 text(.035,.983,title,adj=c(0,1),font=2,cex=1.17,col=navy)
 text(.035,.958,subtitle,adj=c(0,1),cex=.76,col=gray)
 legend(.5,.043,legend=legend_labels,col=legend_cols,pch=legend_pch,bty='n',horiz=TRUE,xjust=.5,yjust=.5,cex=.70,pt.cex=.95,x.intersp=.7)
 text(.5,.012,note,cex=.65,col=gray)
}
setup6<-function(){par(mfrow=c(3,2),oma=c(4.8,0,5,0),family='Helvetica',ps=9,cex=.9,mgp=c(1.8,.45,0),tcl=-.22)}
yrules<-function(y,xr){segments(xr[1],y,xr[2],y,col='#F0F2F5',lwd=.55)}

forest6<-function(component=FALSE){
 setup6()
 for(i in seq_len(nrow(spec))){
  key<-spec$key[i];a<-if(component)co[co$outcome==key&co$type=='Primary',]else ef[ef$outcome==key&ef$type=='Primary'&ef$primary,]
  b<-if(component)co[co$outcome==key&co$type=='Exclude_High',]else ef[ef$outcome==key&ef$type=='Exclude_High'&ef$primary,]
  if(component){tr<-a$term;labels<-gsub(':',' x ',tr,fixed=TRUE)}else{r<-rk[rk$outcome==key&rk$type=='Primary'&rk$treatment!='NAT',];tr<-r$treatment[order(r$rank)];a<-a[match(tr,a$treatment),];b<-b[match(tr,b$treatment),];labels<-tr}
  n<-nrow(a);y<-rev(seq_len(n));xr<-range(c(0,a$lower,a$upper,b$lower,b$upper),na.rm=TRUE);span<-diff(xr);xr<-xr+c(-.04,.06)*span
  par(mar=c(3.0,8.8,2.6,1.1));plot(NA,xlim=xr,ylim=c(.4,n+.7),xlab=paste(spec$sm[i],'(95% CI)'),ylab='',yaxt='n',bty='n',cex.axis=.72,cex.lab=.82)
  yrules(y,xr);abline(v=0,lty=2,col='#8A929C',lwd=.7);axis(2,at=y,labels=labels,las=1,tick=FALSE,cex.axis=.70)
  segments(a$lower,y+.12,a$upper,y+.12,col=blue,lwd=.9);points(a$estimate,y+.12,pch=16,col=blue,cex=.65)
  ok<-b$estimable;segments(b$lower[ok],y[ok]-.12,b$upper[ok],y[ok]-.12,col=orange,lwd=.9);points(b$estimate[ok],y[ok]-.12,pch=15,col=orange,cex=.65)
  if(any(!ok))text(xr[2]-.02*diff(xr),y[!ok]-.12,'NE',adj=1,col=gray,cex=.62)
  ss<-sc[sc$outcome==key&sc$type=='Exclude_High',]
  mtext(paste0(LETTERS[i],'. ',display[key]),side=3,line=1.35,adj=0,font=2,cex=.85,col=navy)
  mtext(paste0('Excluded ',ss$k_removed,' high-risk reports; ',sum(!ok),' term(s) unavailable'),side=3,line=.25,adj=0,cex=.62,col=gray)
 }
 if(!component)overlay('Effect estimates after exclusion of high-risk studies','Complete-regimen NMA; outcome-specific Overall RoB2 judgments supplied for this analysis',c('Primary analysis','Exclude high risk'),c(blue,orange),c(16,15),'Positive values favour the named regimen over NAT. NE: absent or disconnected node. Nominal 95% CIs.')
 else overlay('Component estimates after exclusion of high-risk studies','Fixed primary CNMA specification; REML refitted without reselecting interactions',c('Primary CNMA','Exclude high risk'),c(blue,orange),c(16,15),'NE: coefficient not separately identifiable. Negative interactions indicate departure from additivity, not regimen-level harm.')
}
save_fig('Figure_1_Exclude_high_effects',function()forest6(FALSE))
save_fig('Figure_4_CNMA_Exclude_high',function()forest6(TRUE))

rank6<-function(){
 setup6()
 for(i in seq_len(nrow(spec))){
  key<-spec$key[i];d<-rk[rk$outcome==key&rk$type=='Exclude_High',];d<-d[order(d$baseline_rank_full),];n<-nrow(d);y<-rev(seq_len(n));ok<-d$estimable
  par(mar=c(3.0,8.8,2.6,1.1));plot(NA,xlim=c(0,1.29),ylim=c(.4,n+.7),xlab='P-score',ylab='',yaxt='n',xaxt='n',bty='n',cex.lab=.82)
  axis(1,at=seq(0,1,.25),labels=c('0','.25','.50','.75','1'),cex.axis=.72);axis(2,at=y,labels=d$treatment,las=1,tick=FALSE,cex.axis=.70);yrules(y,c(0,1))
  segments(d$baseline_Pscore_common[ok],y[ok],d$Pscore[ok],y[ok],col='#AAB4BF',lwd=1.3)
  points(d$baseline_Pscore_common[ok],y[ok]+.07,pch=16,col=blue,cex=.65);points(d$Pscore[ok],y[ok]-.07,pch=15,col=orange,cex=.65)
  text(1.04,y[ok],paste0(d$baseline_rank_common[ok],' -> ',d$rank[ok]),adj=0,cex=.65,col=navy)
  if(any(!ok))text(.57,y[!ok],'not estimable',cex=.63,col=gray)
  ss<-sc[sc$outcome==key&sc$type=='Exclude_High',]
  mtext(paste0(LETTERS[i],'. ',display[key]),side=3,line=1.35,adj=0,font=2,cex=.85,col=navy)
  mtext(sprintf('Same-node rank correlation: %.3f',ss$spearman_common),side=3,line=.25,adj=0,cex=.65,col=gray)
 }
 overlay('Treatment ranking after exclusion of high-risk studies','Primary P-scores recomputed over the same retained nodes; rank labels: primary -> sensitivity',c('Primary, same node set','Exclude high risk'),c(blue,orange),c(16,15),'Higher P-scores rank better. P-scores are not probabilities of being best. Missing nodes are not assigned a rank.')
}
save_fig('Figure_2_Exclude_high_rankings',rank6)

loo_ranges<-function(){
 setup6()
 for(i in seq_len(nrow(spec))){
  key<-spec$key[i];a<-ef[ef$outcome==key&ef$type=='Primary'&ef$primary,];d<-ef[ef$outcome==key&ef$type=='LOO'&ef$primary,];r<-rk[rk$outcome==key&rk$type=='Primary'&rk$treatment!='NAT',];tr<-r$treatment[order(r$rank)];a<-a[match(tr,a$treatment),];y<-rev(seq_along(tr));xr<-range(c(a$estimate,d$estimate),na.rm=TRUE);span<-max(.1,diff(xr));plotr<-xr+c(-.10,.38)*span
  par(mar=c(3,8.8,2.6,1.1));plot(NA,xlim=plotr,ylim=c(.4,length(tr)+.7),xlab=paste(spec$sm[i],'(point estimates)'),ylab='',yaxt='n',bty='n',cex.axis=.72,cex.lab=.82)
  yrules(y,plotr);abline(v=0,lty=2,col='#AAB4BF',lwd=.7);axis(2,at=y,labels=tr,las=1,tick=FALSE,cex.axis=.70)
  attempts<-length(unique(d$scenario))
  for(j in seq_along(tr)){
   z<-d[d$treatment==tr[j]&d$estimable,];yy<-y[j];segments(min(z$estimate),yy,max(z$estimate),yy,col=teal,lwd=2)
   points(z$estimate,yy+seq(-.075,.075,length.out=nrow(z)),pch=16,col=adjustcolor(gray,.50),cex=.38)
   s<-z[z$significance_changed,];if(nrow(s))points(s$estimate,rep(yy,nrow(s)),pch=15,col=orange,cex=.85)
   s<-z[z$point_reversal,];if(nrow(s))points(s$estimate,rep(yy,nrow(s)),pch=2,col=purple,cex=1.0,lwd=1.2)
   text(xr[2]+.09*span,yy,paste0(nrow(z),'/',attempts),cex=.65,adj=0,col=gray)
  }
  points(a$estimate,y+.15,pch=18,col=navy,cex=.82)
  mtext(paste0(LETTERS[i],'. ',display[key]),side=3,line=1.35,adj=0,font=2,cex=.85,col=navy)
  mtext(paste0(attempts,' leave-one-report-out refits; n/k shown at right'),side=3,line=.25,adj=0,cex=.65,col=gray)
 }
 overlay('Leave-one-out variation in regimen effect estimates','Ranges and points reflect repeated refits, not confidence intervals; full CIs are in the supplementary forests',c('Primary','LOO estimates','Significance change','Point-sign reversal'),c(navy,gray,orange,purple),c(18,16,15,2),'n/k: estimable / attempted LOO refits. Both observed point-sign reversals remain statistically nonsignificant.')
}
save_fig('Figure_3_LOO_effect_ranges',loo_ranges)

heatmap_outcome<-function(key,index){
 a<-ef[ef$outcome==key&ef$type=='Primary'&ef$primary,];r<-rk[rk$outcome==key&rk$type=='Primary'&rk$treatment!='NAT',];tr<-r$treatment[order(r$rank)];ids<-sort(unique(ef$removed[ef$outcome==key&ef$type=='LOO']));cases<-c('Primary',paste0('LOO_',ids));labels<-c('Primary analysis',vapply(ids,function(id)report_label(key,id),character(1)));n<-length(cases);m<-length(tr)
 effect<-matrix(NA_real_,n,m);changed<-revsign<-matrix(FALSE,n,m);rd<-matrix(NA_real_,n,m)
 for(i in seq_along(cases)){
  d<-ef[ef$outcome==key&ef$scenario==cases[i]&ef$primary,];d<-d[match(tr,d$treatment),]
  effect[i,]<-ifelse(!d$estimable,NA,ifelse(d$conclusion=='Benefit',1,ifelse(d$conclusion=='Harm',3,2)));changed[i,]<-!is.na(d$significance_changed)&d$significance_changed;revsign[i,]<-!is.na(d$point_reversal)&d$point_reversal
  z<-rk[rk$outcome==key&rk$scenario==cases[i],];rd[i,]<-z$delta_rank_common[match(tr,z$treatment)]
 }
 layout(matrix(c(1,2),1,2),widths=c(1.16,1));par(oma=c(5.0,0,4.3,0),family='Helvetica',ps=9,cex=.85)
 for(panel in 1:2){
  par(mar=c(7.8,if(panel==1)13.7 else 1.3,2.2,1.0),xpd=NA)
  plot(NA,xlim=c(.5,m+.5),ylim=c(.5,n+.5),xlab='',ylab='',axes=FALSE,xaxs='i',yaxs='i')
  if(panel==1){
   for(i in seq_len(n))for(j in seq_len(m)){
    y<-n-i+1;v<-effect[i,j];fill<-if(is.na(v))'#B4BBC3'else c(blue,'#E8EDF2',orange)[v]
    rect(j-.5,y-.5,j+.5,y+.5,col=fill,border='white',lwd=.5)
    if(is.na(v))text(j,y,'NE',cex=.52,col='white')
    if(changed[i,j])text(j,y,'*',cex=1.0,font=2,col=if(!is.na(v)&&v==1)'white'else'black')
    if(revsign[i,j])points(j,y,pch=2,cex=.85,col='black',lwd=1.2)
   }
   axis(2,at=rev(seq_len(n)),labels=labels,las=1,tick=FALSE,cex.axis=.64)
   title('A. Effect and significance vs NAT',line=.65,cex.main=.95,col.main=navy,adj=0)
  }else{
   lim<-max(1,abs(rd),na.rm=TRUE);pal<-colorRampPalette(c(blue,'#FFFFFF','#E69F00'))(2*lim+1)
   for(i in seq_len(n))for(j in seq_len(m)){
    y<-n-i+1;v<-rd[i,j];fill<-if(is.na(v))'#B4BBC3'else pal[round(v)+lim+1]
    rect(j-.5,y-.5,j+.5,y+.5,col=fill,border='#DFE5EA',lwd=.5)
    text(j,y,if(is.na(v))'NE'else if(v>0)paste0('+',v)else as.character(v),cex=.58,col=if(!is.na(v)&&v < -0.5*lim)'white'else '#374151')
   }
   title('B. Rank change on common nodes',line=.65,cex.main=.95,col.main=navy,adj=0)
  }
  text(seq_len(m),rep(.21,m),labels=tr,srt=55,adj=1,cex=.64)
 }
 par(fig=c(0,1,0,1),new=TRUE,mar=c(0,0,0,0),oma=c(0,0,0,0),xpd=NA);plot.new();plot.window(xlim=c(0,1),ylim=c(0,1),xaxs='i',yaxs='i')
 text(.04,.983,paste0(display[key],': leave-one-report-out sensitivity'),adj=c(0,1),font=2,cex=1.2,col=navy)
 text(.04,.953,'Each row omits all arms from one R050 report unit. Models and heterogeneity are refitted.',adj=c(0,1),cex=.78,col=gray)
 legend(.50,.070,legend=c('Significant benefit','Not significant','Significant harm','Not estimable'),fill=c(blue,'#E8EDF2',orange,'#B4BBC3'),border=NA,bty='n',horiz=TRUE,xjust=.5,cex=.73)
 text(.5,.035,'* CI significance changed; triangle: point-sign reversal. Rank change: negative = improved, positive = worsened.',cex=.72,col=gray)
 text(.5,.013,'Ranks are compared over the same retained nodes. NE denotes absence/disconnection, not a neutral effect.',cex=.68,col=gray)
}
for(i in seq_len(nrow(spec))){key<-spec$key[i];num<-sprintf('%02d',i);save_fig(paste0('Figure_S',num,'_',key,'_LOO_heatmap'),function()heatmap_outcome(key,i),width=10.8,height=max(7.3,2.7+.19*(1+length(unique(ef$scenario[ef$outcome==key&ef$type=='LOO'])))))}

full_forest<-function(component=FALSE){
 filename<-if(component)'Supplement_CNMA_LOO_forest_43pages.pdf'else'Supplement_NMA_LOO_forest_48pages.pdf'
 pdf(file=file.path(figdir,filename),width=9.6,height=8.3,family='Helvetica',pointsize=9)
 for(k in seq_len(nrow(spec))){
  key<-spec$key[k];data<-if(component)co[co$outcome==key,]else ef[ef$outcome==key&ef$primary,]
  items<-if(component)unique(data$term)else unique(data$treatment)
  for(item in items){
   z<-data[if(component)data$term==item else data$treatment==item,];z<-z[z$type%in%c('Primary','LOO','Exclude_High'),]
   z<-z[order(match(z$type,c('Primary','LOO','Exclude_High')),z$removed),];n<-nrow(z);y<-rev(seq_len(n));ok<-z$estimable
   labels<-vapply(seq_len(n),function(i)if(z$type[i]=='Primary')'Primary analysis'else if(z$type[i]=='Exclude_High')'Exclude all high-risk reports'else report_label(key,z$removed[i]),character(1))
   xr<-range(c(0,z$lower,z$upper),na.rm=TRUE);span<-max(.1,diff(xr));xr<-xr+c(-.06,.08)*span
   par(mar=c(5.5,14.0,5.5,14.0),family='Helvetica',ps=9,cex=.9,xpd=FALSE,mgp=c(1.8,.5,0),tcl=-.2)
   plot(NA,xlim=xr,ylim=c(.3,n+.9),yaxt='n',ylab='',xlab=paste(spec$sm[k],'(95% CI)'),bty='n',cex.axis=.8,cex.lab=.9)
   for(i in which(z$type!='LOO'))rect(xr[1],y[i]-.45,xr[2],y[i]+.45,col=if(z$type[i]=='Primary')'#E9F2F8'else'#FFF1E8',border=NA)
   yrules(y,xr);abline(v=0,col='#9099A5',lty=2,lwd=.7);abline(v=z$estimate[z$type=='Primary'],col=blue,lty=3,lwd=.8)
   cc<-ifelse(z$type=='Primary',blue,ifelse(z$type=='Exclude_High',orange,gray));cc[z$significance_changed%in%TRUE]<-orange;cc[z$point_reversal%in%TRUE]<-purple
   pch<-ifelse(z$type=='Primary',18,ifelse(z$type=='Exclude_High',15,16));pch[z$point_reversal%in%TRUE]<-2
   segments(z$lower[ok],y[ok],z$upper[ok],y[ok],col=cc[ok],lwd=ifelse(z$type[ok]=='Primary',1.4,.8));points(z$estimate[ok],y[ok],col=cc[ok],pch=pch[ok],cex=.68)
   axis(2,at=y,labels=labels,las=1,tick=FALSE,cex.axis=.70)
   txt<-ifelse(ok,sprintf('%.3f [%.3f, %.3f]',z$estimate,z$lower,z$upper),'Not estimable')
   text(xr[2]+.08*diff(xr),y,txt,adj=0,cex=.74,col=cc,xpd=NA)
   text(xr[2]+.08*diff(xr),n+.75,'Estimate [95% CI]',adj=0,font=2,cex=.78,col=navy,xpd=NA)
   title(paste0(display[key],': ',if(component)gsub(':',' x ',item,fixed=TRUE)else paste(item,'vs NAT')),line=3.3,cex.main=1.1,col.main=navy,adj=0)
   mtext(if(component)'Fixed primary CNMA; one-report deletions and high-risk exclusion'else'Complete-regimen NMA; one-report deletions and high-risk exclusion',side=3,line=2.05,adj=0,cex=.80,col=gray)
   mtext('Random-effects Gaussian REML, refitted heterogeneity; nominal normal 95% CIs',side=3,line=1.0,adj=0,cex=.75,col=gray)
   mtext('Orange: changed significance or high-risk exclusion. Purple triangle: point-sign reversal. Gray: other LOO refits.',side=1,line=3.1,cex=.64,col=gray)
   if(component)mtext('Interaction coefficients are conditional departures from additivity; model structure was held fixed.',side=1,line=4.1,cex=.64,col=gray)
  }
 }
 dev.off();cat(filename,'\n')
}
full_forest(FALSE);full_forest(TRUE)
