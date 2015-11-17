########################
## read stata data and correct one data entry error that is inconsistent with the source paper

## download the .dta file from http://journals.plos.org/plosone/article/asset?unique&id=info:doi/10.1371/journal.pone.0114266.s003 and rename it as zinc_iron.DTA
########################

library(foreign); zinDt = read.dta("zinc_iron.DTA"); 
zinDt$studyid = as.factor(zinDt$studyid)

zinDt$studyid = as.factor(as.character(zinDt$studyid))
length(levels(zinDt$studyid))

zinDt[is.na(zinDt$zinc_grp),c("studyid")]

#  zinDt[zinDt$studyid%in%c("188491","194808") ,c("studyid","zinc_grp","iron_grp")]

zinDt$zinc_grp[zinDt$studyid%in%c("188491")]=c(0,0,0)
zinDt$iron_grp[zinDt$studyid%in%c("188491")]=c(1,1,1)

zinDt=zinDt[!zinDt$studyid%in%c("194808"),]
zinDt$studyid=as.factor(as.character(zinDt$studyid)) #  length(levels(zinDt$studyid))

zinc_grpInd = tapply(zinDt$zinc_grp,zinDt$studyid,mean)
iron_grpInd = tapply(zinDt$iron_grp,zinDt$studyid,mean)

table(zinc_grpInd,iron_grpInd) # this table is consistent with the descriptions at the source paper http://journals.plos.org/plosone/article?id=10.1371/journal.pone.0114266

####################### Pure Zinc group
Dtemp = zinDt[zinDt$zinc_grp==1 & zinDt$iron_grp==0 ,c("studyid","visit","week","b_temp","sum_tem","sum_feed")]

## use "zinDt$zinc_grp==0 & zinDt$iron_grp==1" to similarly obtain the file for the pure iron group, "Irontem3121.txt"

## use "zinDt$zinc_grp==1 & zinDt$iron_grp==1" to similarly obtain the file for the zinc plus iron group, "ZincIrontem3121.txt"

## use "zinDt$zinc_grp==0 & zinDt$iron_grp==0" to similarly obtain the file for the placebo group, "Placebotem3121.txt"

DtempID = unique(as.character(Dtemp$studyid))

Sz = length(DtempID);

for(k in 1:Sz)
{ dt = Dtemp[Dtemp$studyid==DtempID[k],c("studyid","week","sum_tem","sum_feed","visit")]
  tem_21[k]=(dt[2,3]  -  dt[1,3])/ ( dt[2,2]-dt[1,2] )
  tem_31[k]=(dt[3,3]  -  dt[1,3])/ ( dt[3,2]-dt[1,2] )
}
 
dt = data.frame(x1=-tem_21*1,x2=-tem_31*1)

library(MissMech)
TestMCARNormality(data = dt,del.lesscases = 2)
# Hawkins test described in Jamshidian 2010

dt = data.frame(x1=-tem_21*50,x2=-tem_31*50)
DataCppProgram = dt;DataCppProgramInd=is.na(DataCppProgram); DataCppProgram[DataCppProgramInd]=1000; write.table(rbind(DataCppProgram,DataCppProgramInd),"/home/guohai/Dropbox/Academic/Thesis/CppEx/Test/Zinctem3121.txt",row.names=FALSE,col.names=FALSE,sep="     ")
