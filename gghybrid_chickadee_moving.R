library(gghybrid)
library(dplyr)
library(tidyr)
library(splitstackshape)

args<-commandArgs(trailingOnly=TRUE)

FILE_IN = paste0(args[1], "_genos.txt")
FILE_OUT = paste0(args[1], "_genos_formatted.txt")
GC_OUT = paste0(args[1], "_clines_posterior.txt")

HI_OUT = paste0(args[1], "_HI.png")
CLINE_OUT = paste0(args[1], "_clines.png")

genos<-read.table(FILE_IN,header=F,stringsAsFactors=FALSE) #Should be chr, pos, and geno columns of vcf 
pos = paste0(genos$V1,"_",genos$V2) #Paste CHR and POS into a single chr_pos coordinate
pop = sort(rep(seq(1,(dim(genos)[2]-2)/15),15)) #15 individuals from each pop
inds = as.character(genos[1,3:length(genos[1,])])
genos_only = genos[2:length(genos[,1]),-c(1,2)]
genos_t = as.data.frame(t(genos_only))
data = cbind(inds,pop,genos_t)
colnames(data) = c("INDLABEL","POPID",seq(1,length(pos[2:length(pos)])))
data = as.data.frame(cSplit(data,as.character(seq(1,length(pos[2:length(pos)]))),"|")) #Split the genotype into alleles
al1<-data[,-c(seq(4,length(data[1,]),2))] #Allele 1 
colnames(al1) = c("INDLABEL","POPID",seq(1,length(pos[2:length(pos)])))
al2<-data[c(1,2,seq(4,length(data[1,]),2))] #Allele 2
colnames(al2) = c("INDLABEL","POPID",seq(1,length(pos[2:length(pos)])))
struc<-NULL #Empty object to merge the two allele dataframes
for (i in 1:length(al1[,1])){
  struc<-rbind(struc,al1[i,],al2[i,]) #Take one row at a time from each allele df, and merge into 1
}
write.table(struc,file=FILE_OUT, quote = F, row.names = F)

dat=read.data(FILE_OUT,nprecol=2,MISSINGVAL = -9,
              NUMINDS = dim(genos)[2]-2, NUMLOCI = dim(genos)[1]-1)

prepdata=data.prep(data=dat$data,
                   loci=dat$loci,
                   alleles=dat$alleles,
                   S0="1", #POPID names for the first parental reference set#
                   S1=as.character(pop[length(pop)]), #POPID names for the second parental reference set#
                   precols=dat$precols,
                   #max.S.MAF = 0.01,	#Filtering by parental minor allele frequency#
                   return.genotype.table=T,
                   return.locus.table=T)

hindlabel= esth(data.prep.object = prepdata$data.prep,
                read.data.precols = dat$precols,
                include.Source = TRUE,	#Set to TRUE if you want hybrid indices for the parental reference individuals#
                #plot.ind = c("P?08-141","PD11-255","PH08-442","PI07-243",
                #             "PI08-501","PX08-520"),
                #plot.col = c("blue","green","cyan","purple","magenta","red"),
                nitt=1000,burnin=100)

hindlabel

setkey(hindlabel$hi,POPID)	#function from data.table, for rapid sorting and subsetting#

#png(filename=HI_OUT,width=1200,height=500)
#abc = plot_h(data=hindlabel$hi,#Subset of POPIDs#
#             test.subject=hindlabel$test.subject,
#             mean.h.by="POPID",			#Calculate the mean hybrid index for each value of the "POPID" column#
#             sort.by=c("mean_h","POPID","h_posterior_mode"),	#Order test subjects along the x axis by the mean hybrid 
#             #index calculated above and also by individual hybrid index#
#             col.group="POPID",
#             group.sep="POPID",
#             fill.source=TRUE,
#             basic.lines=FALSE,
#             source.col=c("blue","red"),
#             source.limits=c("blue","red"),
#             cex=1,pch=16,
#             cex.lab=1.5,cex.main=1.5,ylim=c(0,1))
#dev.off()

hybrids = length(hindlabel$hi$h_posterior_mode[c(hindlabel$hi$h_posterior_mode > 0.35 & hindlabel$hi$h_posterior_mode < 0.65)])
if (hybrids == 0) {
  error_message = paste0(FILE_IN, " doesn't contain any F1-like hybrids.")
  write(error_message,file = "WARNINGS.txt", append = T)
} else {
loc_idxs = which(prepdata$data.prep$locus == "101")
gc=ggcline(
  data.prep.object=prepdata$data.prep[loc_idxs,],    
  esth.object=hindlabel,                  
  #esth.colname="h_posterior_mode",          #Default. Can be replaced with "beta_mean" from the esth object#
  test.subject = "locus",
  poolv = FALSE,
  poolcentre = FALSE,
  include.Source = TRUE,                  #Default is FALSE#
  return.likmeans = TRUE,                 #Default is FALSE#
  read.data.precols=dat$precols,          #Needs an entry#
  fix.subject.v = FALSE,
  #fix.value.v,
  fix.subject.centre = FALSE,
  #fix.value.centre,
  #plot.test.subject = c("1","100"), #Plots are just to check the MCMC is working#
  #plot.col = c("orange","cyan"),                #Plots are just to check the MCMC is working#
  #plot.ylim = c(-3, 5),                         #Plots are just to check the MCMC is working#
  #plot.pch.v.centre = c(1, 3),                  #Plots are just to check the MCMC is working#
  #prior.logitcentre = c(0, sqrt(50)),         #Default#
  #prior.logv = c(0, sqrt(10)),                #Default#
  nitt=3000,                                     #Needs an entry, 5000 should be sufficient for standard per-locus estimates#
  burnin=1000,                                   #Needs an entry, 2000 should be sufficient for standard per-locus estimates#
  start.v = NULL,
  start.centre = NULL,
  init.var.v = NULL,
  init.var.centre = NULL,
  init.cov.vcentre = NULL,
  print.k = 50
)
}
#png(filename=CLINE_OUT,width=1200,height=500)
#plot_clinecurve(
#  ggcline.object=gc$gc,
#  cline.locus=c("101"),
#  locus.column="locus",
#  cline.col=c("orange","blue","green","red","magenta"),
#  null.line.locus=c("20","60","101","120","190"),
#  null.line.col=c("orange","blue","green","red","magenta"),
#  cline.centre.line=c("20","60","101","120","190"),
#  cline.centre.col=c("orange","blue","green","red","magenta")
#)
#Add a title and axis labels:
#title(main = "20,60,101,120,190",xlab="Hybrid index",ylab="Locus allele frequency",cex.main=1.5,cex.lab=1.5)
#dev.off()

write.table(gc$gc,file=GC_OUT,quote = F,row.names = F)
