setwd('C:/Users/Judson/Documents/QC metrics/Ahrne 2016/MaxQuant analysis/PSM,protein FDR = 0.01/combined/txt')

mqout = read.table('evidence.txt',sep='\t',header=TRUE)

#Confirm that FDR was NOT 1% on initial build
revcount = length(which(mqout$Reverse == '+'))
revcount/nrow(mqout)

#Make collated data frame
mqout$seq_charge = paste(mqout$Modified.sequence,mqout$Charge,sep='_') 
seq_charge = unique(mqout$seq_charge)

dat = data.frame(id=seq_charge, sequence=NA, charge=NA,
                 parea0.1_rep1 = NA, parea0.1_rep2 = NA, parea0.1_rep3 = NA, parea0.1_rep4 = NA, parea0.1_rep5 = NA,
                 parea0.02_rep1 = NA, parea0.02_rep2 = NA, parea0.02_rep3 = NA, parea0.02_rep4 = NA, parea0.02_rep5 = NA,
                 parea0.015_rep1 = NA, parea0.015_rep2 = NA, parea0.015_rep3 = NA, parea0.015_rep4 = NA, parea0.015_rep5 = NA,
                 parea0.01_rep1 = NA, parea0.01_rep2 = NA, parea0.01_rep3 = NA, parea0.01_rep4 = NA, parea0.01_rep5 = NA)

exptmap = list(parea0.1_rep1 = 'A14_01043_bhens_MQ', parea0.1_rep2 = 'A14_01044_bhens_MQ', parea0.1_rep3 = 'A14_01045_bhens_MQ', parea0.1_rep4 = 'A14_01046_bhens_MQ', parea0.1_rep5 = 'A14_01047_bhens_MQ',
               parea0.02_rep1 = 'A14_01049_bhens_MQ', parea0.02_rep2 = 'A14_01050_bhens_MQ', parea0.02_rep3 = 'A14_01051_bhens_MQ', parea0.02_rep4 = 'A14_01052_bhens_MQ', parea0.02_rep5 = 'A14_01053_bhens_MQ',
               parea0.015_rep1 = 'A14_01055_bhens_MQ', parea0.015_rep2 = 'A14_01056_bhens_MQ', parea0.015_rep3 = 'A14_01057_bhens_MQ', parea0.015_rep4 = 'A14_01058_bhens_MQ', parea0.015_rep5 = 'A14_01059_bhens_MQ',
               parea0.01_rep1 = 'A14_01061_bhens_MQ', parea0.01_rep2 = 'A14_01062_bhens_MQ', parea0.01_rep3 = 'A14_01063_bhens_MQ', parea0.01_rep4 = 'A14_01064_bhens_MQ', parea0.01_rep5 = 'A14_01065_bhens_MQ')

#Get peak areas from individual replicates to populate collated data frame
peakareas = c()
uniqueraws = as.character(unique(mqout$Raw.file))

for(r in 1:nrow(mqout)){
  rdat = which(dat$id == mqout$seq_charge[r])
  dat$sequence[rdat] = mqout$Sequence[r]
  dat$charge[rdat] = mqout$Charge[r]

  if(length(mqout$Intensity[r]) > 1){
    print(paste('ERROR: row ',r,' has >1 intensity value',sep=''))
  }
  
  if(mqout$seq_charge[r] %in% names(peakareas)){  
      peakareas[[mqout$seq_charge[r]]][,mqout$Raw.file[r]] = mqout$Intensity[r]
  } else {
    peakareas[[mqout$seq_charge[r]]] = matrix(nrow=1, ncol=length(uniqueraws))
    colnames(peakareas[[mqout$seq_charge[r]]]) = uniqueraws
    peakareas[[mqout$seq_charge[r]]][,mqout$Raw.file[r]] = mqout$Intensity[r]
  }
}

#Populate data frame
for(pep in names(peakareas)){
  for(rep in names(exptmap)){
    rdat = which(dat$id == pep)
#    print(paste('dat[rdat,rep] before: ',dat[rdat,rep],sep=''))
#    print(peakareas[[pep]][1,exptmap[[rep]]])
    dat[rdat,rep] = peakareas[[pep]][1,exptmap[[rep]]]
#    print(paste('dat[rdat,rep] after: ',dat[rdat,rep],sep=''))
  }
}
#Calculate peak area avgs

getMeanPA = function(x,amount){
  M=mean(c(as.numeric(x[[eval(paste('parea',amount,'_rep1',sep=''))]]),
    as.numeric(x[[eval(paste('parea',amount,'_rep2',sep=''))]]),
    as.numeric(x[[eval(paste('parea',amount,'_rep3',sep=''))]]),
    as.numeric(x[[eval(paste('parea',amount,'_rep4',sep=''))]]),
    as.numeric(x[[eval(paste('parea',amount,'_rep5',sep=''))]])),na.rm=TRUE)
}

dat$parea0.1 = apply(dat,1,getMeanPA,'0.1')
dat$parea0.02 = apply(dat,1,getMeanPA,'0.02')
dat$parea0.015 = apply(dat,1,getMeanPA,'0.015')
dat$parea0.01 = apply(dat,1,getMeanPA,'0.01')

#Build 1:5, 1:2, 1.5:1, and 10:1 ratios
dat$ratio1to5 = dat$parea0.02 / dat$parea0.1
dat$ratio1to2 = dat$parea0.01 / dat$parea0.02
dat$ratio1.5to1 = dat$parea0.015 / dat$parea0.01
dat$ratio10to1 = dat$parea0.1 / dat$parea0.01

dat$l2ratio1to5 = log(dat$ratio1to5,base=2)
dat$l2ratio1to2 = log(dat$ratio1to2, base=2)
dat$l2ratio1.5to1 = log(dat$ratio1.5to1, base=2)
dat$l2ratio10to1 = log(dat$ratio10to1, base=2)

######################################################################
#Replicate the accuracy boxplot from the Ahrne paper
######################################################################
library(ggplot2)

#order of expected ratios in ratio heatmap columns 1-4 is 0.02/0.1 (5:1); 0.01/0.02 (1:2); 0.015/0.01 (1.5:1); 0.1/0.01 (10:1)
df = data.frame(ratio = c(dat$ratio1to5,dat$ratio1to2,dat$l2ratio1.5to1,dat$l2ratio10to1),
                Dilution = c(rep('x1vs5',nrow(dat)),rep('x1vs2',nrow(dat)),rep('x1vs1.5',nrow(dat)),rep('x10vs1',nrow(dat)))
)

df$l2ratio = log(df$ratio,base=2)

#Reorder factors to match the ordering of Ahrne's plot
df$Dilution = factor(df$Dilution, levels = levels(factor(df$Dilution))[c(4,3,2,1)])

median1vs5 = median(df[df$label=='1:5','l2ratio'],na.rm=TRUE)
median1vs2 = median(df[df$label=='1:2','l2ratio'],na.rm=TRUE)
median1vs1.5 = median(df[df$label=='1:1.5','l2ratio'],na.rm=TRUE)
median10vs1 = median(df[df$label=='10:1','l2ratio'],na.rm=TRUE)

p=ggplot(aes(y = l2ratio, x = Dilution, fill = Dilution), data = df) + 
  geom_boxplot(linetype = 'dashed', color = 'black', outlier.shape=1, outlier.size=3) +
  geom_boxplot(aes(ymin=..lower.., ymax=..upper..), outlier.shape='none') +
  scale_fill_manual(values=c('darkgreen', 'blue', 'purple', 'orange')) +
  scale_x_discrete('Dilution', labels=c('1:5','1:2','1:1.5','10:1')) +
  ylab('LFQ ratio (log2)') +
  scale_y_continuous(limits=c(-4,4),breaks=seq(-4,4,by=2),labels=seq(-4,4,by=2)) + 
  geom_hline(yintercept=log(1/5,base=2),linetype='dashed',color='darkgreen', size=0.75) +
  geom_hline(yintercept=log(1/2,base=2),linetype='dashed',color='blue', size=0.75) +
  geom_hline(yintercept=log(1,base=2),linetype='dashed',color='black', size=0.75) +
  geom_hline(yintercept=log(1.5,base=2),linetype='dashed',color='purple', size=0.75) +
  geom_hline(yintercept=log(10,base=2),linetype='dashed',color='orange', size=0.75) +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.background = element_blank(),
    panel.border = element_rect(color='black', size=1.25, fill=NA, linetype='solid'),
    axis.text.x = element_text(colour='black', size=22),
    axis.text.y = element_text(colour='black', size=22),
    axis.title.x = element_blank(),
    axis.title.y = element_text(margin=margin(t=0, r=10.5, b=0, l=0), size=26),
    axis.ticks.length = unit(0.33, "cm"),
    legend.position = 'none'
  ) +
  labs(x = "Dilution")

ggsave('LFQ accuracy log2 ratio boxplot--Maxquant reanalysis.png',plot=p)

######################################################################
#Calculate q values
######################################################################
library(qvalue)

#Get columns to be used in pvalue calculations
parea0.1cols = grep('parea0.1_rep',colnames(dat))
parea0.02cols = grep('parea0.02_rep',colnames(dat))
parea0.015cols = grep('parea0.015_rep',colnames(dat))
parea0.01cols = grep('parea0.01_rep',colnames(dat))

doTtest = function(df,colrange1,colrange2){
  pareas1 = as.numeric(df[colrange1])
  pareas2 = as.numeric(df[colrange2])
  if(sum(!is.na(pareas1))<3 | sum(!is.na(pareas2))<3){
      return(NA)
  } else {
    return(t.test(pareas1,pareas2)$p.value)
  }
}

#Calculate p values
dat$pvals10to1 = apply(dat,1,doTtest,parea0.1cols,parea0.01cols)
dat$pvals1to5 = apply(dat,1,doTtest,parea0.1cols,parea0.02cols)
dat$pvals1to2 = apply(dat,1,doTtest,parea0.01cols,parea0.02cols)
dat$pvals1.5to1 = apply(dat,1,doTtest,parea0.015cols,parea0.01cols)

pvals10to1NotNA = which(!is.na(dat$pvals10to1))
pvals1to5NotNA = which(!is.na(dat$pvals1to5))
pvals1to2NotNA = which(!is.na(dat$pvals1to2))
pvals1.5to1NotNA = which(!is.na(dat$pvals1.5to1))

#Now do them q values
dat$qvals10to1 = dat$qvals10to1BH = NA
dat$qvals1to5 = dat$qvals1to5BH = NA
dat$qvals1to2 = dat$qvals1to2BH = NA
dat$qvals1.5to1 = dat$qvals1.5to1BH = NA

dat$qvals10to1[pvals10to1NotNA] = qvalue(p = dat$pvals10to1[pvals10to1NotNA])$qvalue
dat$qvals1to5[pvals1to5NotNA] = qvalue(p = dat$pvals1to5[pvals1to5NotNA])$qvalue
dat$qvals1to2[pvals1to2NotNA] = qvalue(p = dat$pvals1to2[pvals1to2NotNA])$qvalue
dat$qvals1.5to1[pvals1.5to1NotNA] = qvalue(p = dat$pvals1.5to1[pvals1.5to1NotNA])$qvalue

dat$qvals10to1BH[pvals10to1NotNA] = qvalue(p = dat$pvals10to1[pvals10to1NotNA],lambda=0)$qvalue
dat$qvals1to5BH[pvals1to5NotNA] = qvalue(p = dat$pvals1to5[pvals1to5NotNA],lambda=0)$qvalue
dat$qvals1to2BH[pvals1to2NotNA] = qvalue(p = dat$pvals1to2[pvals1to2NotNA],lambda=0)$qvalue
dat$qvals1.5to1BH[pvals1.5to1NotNA] = qvalue(p = dat$pvals1.5to1[pvals1.5to1NotNA],lambda=0)$qvalue

#Alternatively--what if we controlled FDR considering all pooled pvalues?