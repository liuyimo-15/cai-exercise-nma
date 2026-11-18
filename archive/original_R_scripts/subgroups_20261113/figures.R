options(stringsAsFactors=FALSE)
clean<-'--clean'%in%commandArgs(trailingOnly=TRUE)
b<-'work/subgroups_20261113';o<-'outputs/2026-11-13_亚组分析补充'
su<-read.csv(file.path(b,'summary.csv'));ef<-read.csv(file.path(b,'effects.csv'));ad<-read.csv(file.path(b,'adult_restrictions.csv'))
display<-c(CAIT='CAIT',FAAM_ADL='FAAM-ADL',Sports='Sports function',Anterior='Anterior balance',Posteromedial='Posteromedial balance',Posterolateral='Posterolateral balance')
cols<-c('#2D6484','#B96834');gray<-'#52606D';navy<-'#243B53'
groups<-list(duration=c('Short_lt6','Long_ge6'),tool=c('SEBT','YBT'),population_group=c('Recreational_active','Athlete'))
labels<-list(duration=c('<6 weeks','>=6 weeks'),tool=c('SEBT / mSEBT','YBT'),population_group=c('Recreational / active','Athlete samples'))
setup<-function(){par(mfrow=c(3,2),oma=if(clean)c(1,1,1,1)else c(6,1,5,1),family='Arial',ps=10,cex=1,mgp=c(1.9,.55,0),tcl=-.22,xpd=FALSE)}
overlay<-function(title,subtitle,foot){
 if(clean)return(invisible(NULL))
 par(fig=c(0,1,0,1),new=TRUE,mar=c(0,0,0,0),oma=c(0,0,0,0),xpd=NA);plot.new();plot.window(xlim=c(0,1),ylim=c(0,1),xaxs='i',yaxs='i')
 text(.045,.98,title,adj=c(0,1),font=2,cex=1.18,col=navy);text(.045,.952,subtitle,adj=c(0,1),cex=.82,col=gray)
 text(.045,.064,foot,adj=c(0,1),cex=.77,col=gray)
}
forest<-function(var){
 setup();ss<-su[su$variable==var&su$scenario=='R050'&!is.na(su$p_interaction),]
 for(i in 1:6){
  if(i>nrow(ss)){plot.new();next};s<-ss[i,];rr<-ef[ef$outcome==s$outcome&ef$variable==var&ef$scenario=='R050'&ef$common_node,];tt<-sort(unique(rr$treatment));ys<-rev(seq_along(tt));xr<-range(c(0,rr$lower,rr$upper));xr<-xr+c(-.05,.10)*max(.1,diff(xr))
  par(mar=c(3.5,7,3.8,1.2));plot(NA,xlim=xr,ylim=c(.4,length(tt)+1.1),xlab=paste(if(s$outcome%in%c('CAIT','FAAM_ADL'))'MD'else'Hedges g','vs NAT (95% CI)'),ylab='',yaxt='n',bty='n',cex.lab=.87,cex.axis=.82)
  segments(xr[1],ys,xr[2],ys,col='#E9EDF1',lwd=.6);abline(v=0,lty=2,col='#8A929C',lwd=.8);axis(2,at=ys,labels=tt,las=1,tick=FALSE,cex.axis=.82)
  for(j in 1:2){a<-rr[rr$group==groups[[var]][j],];a<-a[match(tt,a$treatment),];yy<-ys+if(j==1).13 else-.13;segments(a$lower,yy,a$upper,yy,col=cols[j],lwd=1.15);points(a$estimate,yy,pch=if(j==1)16 else 15,col=cols[j],cex=.82)}
  mtext(display[s$outcome],side=3,line=2.4,adj=0,font=2,cex=1.02,col=navy);mtext(sprintf('Global interaction P=%.3f; BH q=%.3f',s$p_interaction,s$q_BH),side=3,line=1.15,adj=0,cex=.79,col=gray)
  legend('top',legend=paste0(labels[[var]],' (k=',c(s$k0,s$k1),')'),col=cols,pch=c(16,15),bty='n',horiz=TRUE,cex=.69,x.intersp=.65,y.intersp=.7,inset=c(0,-.02))
 }
 title<-switch(var,duration='Treatment duration',tool='Balance measurement tool',population_group='Sport / activity population')
 overlay(paste('Exploratory network subgroup analysis:',title),'Subgroup-specific treatment effects; common REML tau-squared; R050 input preserved.',paste('Only shared NAT-connected nodes shown. k counts analysis/report IDs, not confirmed independent trials.',
 'Positive effects follow the original improvement direction. Post hoc analyses; nominal 95% confidence intervals.',
 'No global interaction P <0.05. Absence of interaction evidence does not establish equivalence.',sep='\n'))
}
adult<-function(){
 setup()
 for(key in names(display)){
  a<-ad[ad$outcome==key&ad$restriction=='Adult_verified'&ad$estimable,];a<-a[order(a$treatment),];ys<-rev(seq_len(nrow(a)));xr<-range(c(0,a$lower,a$upper,a$baseline_lower,a$baseline_upper));xr<-xr+c(-.05,.10)*max(.1,diff(xr))
  par(mar=c(3.5,7,3.8,1.2));plot(NA,xlim=xr,ylim=c(.4,nrow(a)+1.1),xlab=paste(if(key%in%c('CAIT','FAAM_ADL'))'MD'else'Hedges g','vs NAT (95% CI)'),ylab='',yaxt='n',bty='n',cex.lab=.87,cex.axis=.82)
  segments(xr[1],ys,xr[2],ys,col='#E9EDF1',lwd=.6);abline(v=0,lty=2,col='#8A929C',lwd=.8);axis(2,at=ys,labels=a$treatment,las=1,tick=FALSE,cex.axis=.78)
  segments(a$baseline_lower,ys+.13,a$baseline_upper,ys+.13,col=cols[1],lwd=1.15);points(a$baseline_estimate,ys+.13,pch=16,col=cols[1],cex=.75)
  segments(a$lower,ys-.13,a$upper,ys-.13,col=cols[2],lwd=1.15);points(a$estimate,ys-.13,pch=15,col=cols[2],cex=.75)
  mtext(display[key],side=3,line=2.4,adj=0,font=2,cex=1.02,col=navy)
  if(!clean)mtext(paste('Verified-adult NAT-connected report IDs: k=',a$k[1],sep=''),side=3,line=1.15,adj=0,cex=.79,col=gray)
  legend('top',legend=if(clean)c('Primary analysis',paste0('Adult-only (k=',a$k[1],')'))else c('Original R050','Verified-adult restriction'),col=cols,pch=c(16,15),bty='n',horiz=TRUE,cex=.69,x.intersp=.65,y.intersp=.7,inset=c(0,-.02))
 }
 overlay('Age restriction sensitivity: verified adult-only samples','Each restricted network re-estimates REML tau-squared; original R050 shown for reference.',paste('Age-unconfirmed reports are excluded from this strict adult restriction. Missing/disconnected nodes are in the workbook.',
 'Only retained NAT-connected nodes are plotted. The second, less restrictive age sensitivity is tabulated separately.',
 'This is a restriction sensitivity analysis, not an adult-versus-youth interaction test.',sep='\n'))
}
quartz(type='pdf',file=file.path(o,if(clean)'六结局_亚组分析森林图_精简版.pdf'else'六结局_亚组分析森林图.pdf'),width=12.0,height=11.8,family='Arial',pointsize=10)
for(v in names(groups))forest(v);adult();dev.off()
for(v in c(names(groups),'adult')){
 png(file.path(b,paste0('forest_',v,if(clean)'_clean'else'','.png')),width=12,height=11.8,units='in',res=180,type='quartz',pointsize=10)
 if(v=='adult')adult()else forest(v);dev.off()
}
cat('4-page vector forest PDF exported\n')
