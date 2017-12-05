library(gdata)

setwd("C:/Users/Judson/Documents/QC metrics/Ahrne 2016")

dat = read.xls("Ahrne_bhens_1fdr_HsapMS2peakNorm_exportPt1.xls")
datmq = read.table("Ahrne MQ reanalysis (PSM,protein FDR = 0.01).txt",sep='\t',header=TRUE) #output file from 'make maxquant ratios for accuracy analysis.r' script

#have to replace tryptic sites with '_' in our pipeline's output
dat$sequence = gsub('^.\\.','_',dat$peptide.data.2..assigned.sequence)
dat$sequence = gsub('\\..','_',dat$sequence)
#also remove mod characters
dat$sequence = gsub('#','',dat$sequence)

dat$seq_charge = paste(dat$sequence,dat$charge.state.peptide,sep='_')

rownames(datmq) = datmq$seq_charge

#now make merged data frame of peptides observed in both pipelines.
#autofill output has redundant peptides, so take the one with most qvals<.01 and the best replicate support
# NOTE: Not finished
# intersectors = intersect(dat$seq_charge, datmq$id)
# merged = lapply(intersectors,function(pepid){
  # df = data.frame(id = pepids, mqratio1to5 = datmq[pepid,''], mqratio1to5 = datmq[pepid,''], 
                  # mqqval = NA, afratio = NA, afqval = NA)
  # mqrow
  # return(df)  
# })

#################################################################################################
#Get number of inverters in each dataset
#################################################################################################
comparison = '10to1'
coln = 4
ratiostocomparemq = na.omit(datmq[,paste('ratio',comparison,sep='')])
#order of expected ratios in ratio heatmap columns 1-4 is 0.02/0.1 (1:5); 0.01/0.02 (1:2); 0.015/0.01 (1.5:1); 0.1/0.01 (10:1)
ratiostocompare = na.omit(dat[,paste('SILAC.ratio.32.for.user.selected.SILAC.timepoint',coln,sep='')])
#how many negative fold changes (num < denom)?
nfcmq = sum(ratiostocomparemq <= 1)
reciprocal = length(ratiostocomparemq) - nfcmq
print(paste('MQ not inverted: ',nfcmq, ' (',nfcmq/length(ratiostocomparemq) * 100,'%)',sep=''))
print(paste('MQ inverted: ',reciprocal, ' (',reciprocal/length(ratiostocomparemq) * 100,'%)',sep=''))

nfc = sum(ratiostocompare <= 1)
reciprocal = length(ratiostocompare) - nfc
print(paste('AF not inverted: ',nfc, ' (',nfc/length(ratiostocompare) * 100,'%)',sep=''))
print(paste('AF inverted: ',reciprocal, ' (',reciprocal/length(ratiostocompare) * 100,'%)',sep=''))
#################################################################################################
#Same as above, but only consider ratios with qval < .01 when calculating proportion of inversions
#################################################################################################
comparison = '10to1'
qvalcolstring = paste('qvals',comparison,sep='')
mqindices = which(datmq[,qvalcolstring]<.01)
ratiosmq = datmq[mqindices,paste('ratio',comparison,sep='')]

nfcmq = sum(ratiosmq < 1)
pfcmq = length(ratiosmq) - nfcmq
#how many negative fold changes (num < denom)?
print(paste('MQ - fold change: ',nfcmq, ' (',nfcmq/length(ratiosmq) * 100,'%)',sep=''))
#how many positive fold changes (num < denom)?
print(paste('MQ + fold change: ',pfcmq, ' (',pfcmq/length(ratiosmq) * 100,'%)',sep=''))

coln = 4
#order of expected ratios in ratio heatmap columns 1-4 is 0.02/0.1 (1:5); 0.01/0.02 (1:2); 0.015/0.01 (1.5:1); 0.1/0.01 (10:1)
qvalcolstring = paste('qvalues.for.SILAC.timepoint',coln,sep='')
afindices = which(dat[,qvalcolstring]<.01)
ratios = dat[afindices,paste('SILAC.ratio.32.for.user.selected.SILAC.timepoint',coln,sep='')]

nfc = sum(ratios < 1)
pfc = length(ratios) - nfc
print(paste('AF - fold change: ',nfc, ' (',nfc/length(ratios) * 100,'%)',sep=''))
print(paste('AF + fold change: ',pfc, ' (',pfc/length(ratios) * 100,'%)',sep=''))