# 3 methods for deriving the global normalization factors applied to each LC/MS run in Ahrne's dataset. We do not know what these normalization
# factors are, so these methods represent 3 best guesses for how their normalization was performed.

setwd("C:/Users/Judson/Documents/QC metrics/Ahrne 2016")

### Method 1: sum up the peak areas from all detected (OR autofilled) peaks
library(gdata)

dat=read.xls("Ahrne_hsapiens_1fdr_allStds1_metricsPt1.xls")

A14.01043 = sum(as.numeric(dat$peakarea.manual.1.rep1.thresholded.timepoint4), na.rm=TRUE)
A14.01044 = sum(as.numeric(dat$peakarea.manual.1.rep2.thresholded.timepoint4), na.rm=TRUE)
A14.01045 = sum(as.numeric(dat$peakarea.manual.1.rep3.thresholded.timepoint4), na.rm=TRUE)
A14.01046 = sum(as.numeric(dat$peakarea.manual.1.rep4.thresholded.timepoint4), na.rm=TRUE)
A14.01047 = sum(as.numeric(dat$peakarea.manual.1.rep5.thresholded.timepoint4), na.rm=TRUE)
A14.01049 = sum(as.numeric(dat$peakarea.manual.1.rep1.thresholded.timepoint3), na.rm=TRUE)
A14.01050 = sum(as.numeric(dat$peakarea.manual.1.rep2.thresholded.timepoint3), na.rm=TRUE)
A14.01051 = sum(as.numeric(dat$peakarea.manual.1.rep3.thresholded.timepoint3), na.rm=TRUE)
A14.01052 = sum(as.numeric(dat$peakarea.manual.1.rep4.thresholded.timepoint3), na.rm=TRUE)
A14.01053 = sum(as.numeric(dat$peakarea.manual.1.rep5.thresholded.timepoint3), na.rm=TRUE)
A14.01055 = sum(as.numeric(dat$peakarea.manual.1.rep1.thresholded.timepoint2), na.rm=TRUE)
A14.01056 = sum(as.numeric(dat$peakarea.manual.1.rep2.thresholded.timepoint2), na.rm=TRUE)
A14.01057 = sum(as.numeric(dat$peakarea.manual.1.rep3.thresholded.timepoint2), na.rm=TRUE)
A14.01058 = sum(as.numeric(dat$peakarea.manual.1.rep4.thresholded.timepoint2), na.rm=TRUE)
A14.01059 = sum(as.numeric(dat$peakarea.manual.1.rep5.thresholded.timepoint2), na.rm=TRUE)
A14.01061 = sum(as.numeric(dat$peakarea.manual.1.rep1.thresholded.timepoint1), na.rm=TRUE)
A14.01062 = sum(as.numeric(dat$peakarea.manual.1.rep2.thresholded.timepoint1), na.rm=TRUE)
A14.01063 = sum(as.numeric(dat$peakarea.manual.1.rep3.thresholded.timepoint1), na.rm=TRUE)
A14.01064 = sum(as.numeric(dat$peakarea.manual.1.rep4.thresholded.timepoint1), na.rm=TRUE)
A14.01065 = sum(as.numeric(dat$peakarea.manual.1.rep5.thresholded.timepoint1), na.rm=TRUE)

### Method 2: get TIC area directly from raw data
library(xcms)

path = 'mzXMLs'
mzXMLlist = list.files(path)

TIClist = lapply(mzXMLlist, function(mzXML){
	mzXMLpath = paste(path,mzXML,sep='/')
	print(mzXMLpath)
	xcmsR = xcmsRaw(mzXMLpath, profstep=0)
	TICarea = sum(xcmsR@tic)
	return(data.frame(file=mzXML, TICarea = TICarea))
})

TICdf = do.call(rbind,TIClist)
write.table(TICdf, 'Ahrne LFQ data TICs.txt', sep='\t', quote=FALSE, row.names=FALSE, col.names=TRUE)

### Method 3: sum detected peak areas only in each experiment
library(readxl)

path = 'A14_hsapiens_expt_dumps'
dumplist = list.files(path)

exptlist = lapply(dumplist, function(exptname){
  exptpath = paste(path, exptname, sep='/')
  print(exptpath)
  expt = read_excel(exptpath)
  peakareaSum = sum(expt$'peak area final',na.rm=TRUE)
  return(data.frame(file=exptname, peakareaSum = peakareaSum))
})

exptdf = do.call(rbind,exptlist)
write.table(exptdf, 'Ahrne LFQ data ms2 detected peakarea sums.txt', sep='\t', quote=FALSE, row.names=FALSE, col.names=TRUE)
