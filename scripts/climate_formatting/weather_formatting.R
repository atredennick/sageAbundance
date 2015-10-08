# format Daymet data and calculate seasonal totals

## STUDY AREA 1 FILES ##
# dFile <- "/Users/atredenn/Dropbox/sagebrush_class_2013/Wyoming/studyarea1/climate/DAYMET/center_weather_data_SA1.csv"
# outfile <- "/Users/atredenn/Dropbox/sagebrush_class_2013/Wyoming/studyarea1/climate/DAYMET/FormattedClimate_WY_SA1.csv"

## STUDY AREA 2 FILES ##
dFile <- "/Users/atredenn/Dropbox/sagebrush_class_2013/Wyoming/studyarea2/climate/DAYMET/center_weather_data_SA2.csv"
outfile <- "/Users/atredenn/Dropbox/sagebrush_class_2013/Wyoming/studyarea2/climate/DAYMET/FormattedClimate_WY_SA2.csv"

D=read.csv(dFile)

D$tmean=rowMeans(D[,3:4])
D=D[,c("year","yday","tmean","prcp..mm.day.")]
names(D)[4]="prcp"
D$CalendarYear=D$year
D$year=ifelse(D$yday>=274,D$year+1,D$year)
D$summer=ifelse(D$yday<274 & D$yday>181,"yes","no")
D$spring=ifelse(D$yday<=181 & D$yday>90,"yes","no")

tmeanD=aggregate(tmean~year,data=D[which(D$spring=="yes"),],FUN=mean)
names(tmeanD)[2]="TmeanSpr2"
tmeanD1=tmeanD
tmeanD1$year=tmeanD1$year+1
names(tmeanD1)[2]="TmeanSpr1"
tmeanD=merge(tmeanD,tmeanD1,all.x=T)

pptD=aggregate(prcp~year+summer,data=D,FUN=sum)
pptD=reshape(pptD,direction="wide",idvar=c("year"),timevar="summer")
names(pptD)[2:3]=c("ppt2","pptSummer2")
tmp=pptD
tmp$year=tmp$year+1
names(tmp)[2:3]=c("ppt1","pptSummer1")
pptD=merge(pptD,tmp,all.x=T)
tmp=pptD
tmp$year=tmp$year+2
tmp$pptLag=tmp$ppt2+tmp$pptSummer2
tmp=tmp[,c("year","pptLag")]
pptD=merge(pptD,tmp,all.x=T)

climD=merge(pptD,tmeanD,all.x=T)

write.csv(climD,outfile,row.names=F)
