library("grid")
library("ggplot2")
library("gplots")
library("ggcorrplot")
library("gridExtra")
library("viridis")
library("reshape2")
library("gridExtra")
library("RColorBrewer")
library("scales")
library("stringr")
library("rstatix")
library("dplyr")
library("tidyr")
source('./COVIDvaccines_PlottingFunctions.R')
sessionInfo()

knitr::opts_chunk$set(fig.width=12, fig.height=8) 


mergedData <- readRDS(file = "mergedData.Rds")
demog.melt <- readRDS(file = "demog.melt.Rds")
# colors:   COVID-exp B5B2F1 (purple)             COVID-naive FFDFB1 (orange)

#' ------------------ Cohort description --------------------------
#'
# ggplot(demog, aes(x = "", y=Vaccine.Type, fill=Vaccine.Type)) + geom_bar(stat='identity') + coord_polar("y",start=0) + theme_void() + scale_fill_brewer(palette = "Set1")
# ggsave(filename = "Images/VaccineTypes.pdf")


temp <- mergedData
temp$timeCategory <- factor(temp$timeCategory, levels = c("Baseline", "Post 1st dose", "two Weeks", "Pre 2nd dose", "Post 2nd dose", 
                                                          "One month post\n2nd dose", "Four months post\n2nd dose"))
temp <- temp[! is.na(temp$timeCategory), ]
temp <- temp[-which(temp$Alias == "PHI-077" | temp$Alias == "PHI-054" | temp$Alias =="PHI-075" | temp$Alias =="HV-095" ), ]
temp$Record.ID <- factor(temp$Record.ID)
a <- ggplot(data = temp, aes(x = DPV, y = Record.ID, group = Record.ID)) + geom_vline(xintercept = 0, linetype = "dashed", alpha=0.5) + geom_path() + 
  geom_point(aes(fill = timeCategory, shape = timeCategory), size=4) + theme_bw()  + ggtitle("Cohort overview") + xlab("Days relative to vaccine dose 1") + ylab("") + 
  scale_shape_manual(values=c(20:25,16)) + scale_fill_viridis_d() + 
  theme(axis.text = element_text(color="black",size=16, angle=45, hjust = 1, vjust=0), plot.title = element_text(color="black",size=24), axis.text.y = element_blank(), 
        legend.text = element_text(size=16), legend.title = element_blank(), axis.title = element_text(color="black", size=16)) + 
  scale_y_discrete(limits = rev(levels(temp$Record.ID)))+ 
  scale_x_continuous( breaks = seq(-100,200,20))
a
# ggsave(filename = "Images/FullCohort_overview.pdf", width=7)
# plotly::ggplotly(a)

#' ------------------ First x subjects and no longitudinal --------------------------
#'
# 
# if("PHI-046" %in% mergedData$Alias)  {mergedData <- mergedData[-which(mergedData$Alias == "PHI-046"),] }          # only Moderna recipient, will temporarily exclude
# if("HV-083" %in% mergedData$Alias)  {mergedData <- mergedData[-which(mergedData$Alias == "HV-083"),]   }          # exclude due to active COVID at time of vaccination so uncertain cohort    
# if("HV-089" %in% mergedData$Alias)  {mergedData <- mergedData[-which(mergedData$Alias == "HV-089"),]   }          # transplant patient who was already fully vaccinated at enrollment    
# if("PHI-077" %in% mergedData$Alias)  {mergedData <- mergedData[-which(mergedData$Alias == "PHI-077"),]   }          # already fully vaccinated at enrollment      

keepList <- paste0("CV-",sprintf("%03d", seq(1,43,1)))
keepList <- keepList[-which(keepList == "CV-013" |       # developed COVID shortly after 1st vaccination so unclear status
                              keepList == "CV-031" |       # received Moderna
                              keepList == "CV-038" |       # fully vax transplant patient at enrollment
                              keepList == "CV-040" |       # enrolled late so don't have recent baseline or post 1st dose   
                              keepList == "CV-041" |       # incomplete timecourse, vax date 2 unknown and lost to followup? 
                              keepList == "CV-042")]       # received Moderna
mergedData <- mergedData[ mergedData$Record.ID %in% keepList, ]
temp <- mergedData[which(mergedData$Record.ID %in% keepList),] %>% group_by(Prior.COVID.infection., timeCategory) %>% get_summary_stats(type = 'common')
# write.csv(temp, file = "summaryStatistics.csv")

subsetData <- subset(mergedData, timeCategory == "Baseline")
table(subsetData$Prior.COVID.infection., subsetData$Sex)
table(subsetData$Prior.COVID.infection., subsetData$Race)
median(subsetData$DPO.covid, na.rm=T)
range(subsetData$DPO.covid, na.rm=T)


temp <- mergedData[which(mergedData$Record.ID %in% keepList), ]  
temp <- temp[-which(temp$timeCategory == "2 wks Post 2nd dose" | temp$timeCategory == "" | temp$shortForm == "5W" | temp$shortForm == "4M"),]
temp$timeCategory <- factor(temp$timeCategory, levels = c("Baseline", "Post 1st dose", "two Weeks", "Pre 2nd dose", "Post 2nd dose", "One month post\n2nd dose")) #,"Four months post\n2nd dose"))
if( anyNA(temp$timeCategory))  {   temp <- temp[-is.na(temp$timeCategory),]  }

ggplot(data = temp, aes(x = DPV, y = Record.ID, group = Record.ID)) + geom_vline(xintercept = 0, linetype = "dashed", alpha=0.5) + geom_path() + 
  geom_point(aes(fill = timeCategory, shape = timeCategory), size=4) + theme_bw()  + xlab("Days relative to vaccine dose 1") + ylab("") + 
  scale_shape_manual(values=c(18:25)) + scale_fill_viridis_d() + 
  theme(axis.text = element_text(color="black",size=16), axis.title = element_text(color="black",size=16), axis.text.y = element_blank()) + 
  scale_x_continuous(breaks = seq(-100,200,20)) 
# ggsave(filename = "./Images/Subject_timecourse_overview.pdf")


#' -----------------------  FINAL DATA OBJECT  --------------------------------

# write.csv(mergedData, file = "../Data/mergedData_postExclusions.csv")
# saveRDS(mergedData, file = "mergedData_postExclusions.Rds")



#' -----------------------  OMIQ metaData export  --------------------------------
#' 

##  *****   REMEMBER EXCLUSIONS ABOVE !!!!    *******
# ##  *****   also deleted CV-034_PHI-058_bL_CPT because it had no CD4 stain    **********
# 
# metaData <- mergedData[,c("fcsFile","Alias","Label","timeCategory","Sex","Prior.COVID.infection.","DPO.covid","DPV")]
# metaData <- metaData[-which(metaData$Label == "PHI-021_V1"),]
# # fileList <- data.frame(original = list.files(path = "../../COVIDvax_export_CD16loCD14lo/"))   # all the fcs files that will go to OMIQ.ai
# fileList[! fileList %in% metaData$fcsFile]
# metaData$fcsFile[! metaData$fcsFile %in% fileList]
# 
# fileList <-  list.files(path = "../../COVIDvax_export_CD16loCD14lo/")
# metaData <- metaData[which(metaData$fcsFile %in% fileList),]
# metaData <- metaData[,c("fcsFile","Alias","Label","timeCategory","Sex","Prior.COVID.infection.","DPO.covid","DPV")]
# omiqID <- read.csv(file = "../../COVIDvax_export_CD16loCD14lo/20210316_metadata.csv")
# omiqID$Filename[! omiqID$Filename %in% metaData$fcsFile]
# 
# omiqID <- merge(x = omiqID, y = metaData, by.x = "Filename", by.y = "fcsFile")
# omiqID <- omiqID[, c(2,1, 3:ncol(omiqID))]
# # write.csv(omiqID, file = "../Flow cytometry/COVIDvax_metaData_OMIQ.ai.csv", row.names = F)


#' ------------------ Activated T cells analyses --------------------------
#'

subsetData <- subset(mergedData,  timeCategory != "two Weeks"& timeCategory != "2 wks post 2nd dose");
# subsetData <- subsetData[which(subsetData$Tube == "HEP"),] 
subsetData <- subset(subsetData, Record.ID != "CV-011" & Record.ID != "CV-012" )        # absence of Ki67 stain
subsetData$timeCategory <- factor(subsetData$timeCategory, levels = c("Baseline", "Post 1st dose", "Pre 2nd dose", "Post 2nd dose", "One month post\n2nd dose"))
prePostTime(subsetData, xData = "DPV", yData="CD4_.CD38.Ki67._FreqParent", fillParam = "Prior.COVID.infection.", groupby="Record.ID", title = "CD4 - Post 1st dose",
            xLabel = "Days", yLabel = "Ki67+CD38+ (% CD4)", repMeasures = F, exponential=F) +  coord_cartesian(xlim = c(-1,12)) + geom_vline(xintercept = 0,linetype="dashed", alpha=0.5)
# ggsave(filename = "./Images/ActivCD4_bothCohorts_Vax1_continuousTime.pdf")
prePostTime(subsetData, xData = "DPV", yData="CD8_.CD38.Ki67._FreqParent", fillParam = "Prior.COVID.infection.", groupby="Record.ID", title = "CD8 - Post 1st dose",
            xLabel = "Days", yLabel = "Ki67+CD38+ (% CD8)", repMeasures = F, exponential=F) +  coord_cartesian(xlim = c(-1,12)) + geom_vline(xintercept = 0,linetype="dashed" , alpha=0.5)
# ggsave(filename = "./Images/ActivCD8_bothCohorts_Vax1_continuousTime.pdf")
subsetData$DPV <- subsetData$DPV - as.numeric(difftime(subsetData$Vaccine.2.date, subsetData$Vaccine.1.date, units="days" ) )
prePostTime(subsetData, xData = "DPV", yData="CD4_.CD38.Ki67._FreqParent", fillParam = "Prior.COVID.infection.", groupby="Record.ID", title = "CD4 - Post 2nd dose",
            xLabel = "Days", yLabel = "Ki67+CD38+ (% CD4)", repMeasures = F, exponential=F) +  coord_cartesian(xlim = c(-5,12)) + geom_vline(xintercept = 0,linetype="dashed", alpha=0.5)
# ggsave(filename = "./Images/ActivCD4_bothCohorts_Vax2_continuousTime.pdf")
prePostTime(subsetData, xData = "DPV", yData="CD8_.CD38.Ki67._FreqParent", fillParam = "Prior.COVID.infection.", groupby="Record.ID", title = "CD8 - Post 2nd dose",
            xLabel = "Days", yLabel = "Ki67+CD38+ (% CD8)", repMeasures = F, exponential=F) +  coord_cartesian(xlim = c(-5,12))  + geom_vline(xintercept = 0,linetype="dashed" , alpha=0.5)
# ggsave(filename = "./Images/ActivCD8_bothCohorts_Vax2_continuousTime.pdf")

##' *************** activated CD4 ********************
prePostTime(data = subsetData, xData = "timeCategory", yData="CD4_.CD38.Ki67._FreqParent", fillParam = "Prior.COVID.infection.", groupby="Record.ID", title = "Activated CD4", 
            xLabel = " ", yLabel = "Ki67+CD38+ (% CD4)", repMeasures = F, exponential=F, newform = T)   
# ggsave(filename = "./Images/ActivCD4_bothCohorts_overTime.pdf")

bartlett.test(CD4_.CD38.Ki67._FreqParent ~ timeCategory, data=subset(subsetData, Prior.COVID.infection. == 'No'))
kruskal_test(formula = CD4_.CD38.Ki67._FreqParent ~ timeCategory, data=subset(subsetData, Prior.COVID.infection. == 'No')) 
dunn_test(CD4_.CD38.Ki67._FreqParent ~ timeCategory, data=subset(subsetData, Prior.COVID.infection. == 'No'))
bartlett.test(CD4_.CD38.Ki67._FreqParent ~ timeCategory, data=subset(subsetData, Prior.COVID.infection. == 'Yes'))
tukey_hsd( aov(CD4_.CD38.Ki67._FreqParent ~ timeCategory, data=subset(subsetData, Prior.COVID.infection. == 'Yes')) )

prePostTime(data = subsetData, xData = "timeCategory", yData="CD4_.CD38.Ki67._FreqParent", fillParam = "Prior.COVID.infection.", groupby="Record.ID", title = "Activated CD4", 
            xLabel = " ", yLabel = "Ki67+CD38+ (% CD4)", repMeasures = F, exponential=F, newform = T, recentCOVID = T)  
# ggsave(filename = "./Images/ActivCD4_bothCohorts_overTime_recentCOVID.pdf")


##' *************** activated CD8 ********************
prePostTime(data = subsetData, xData = "timeCategory", yData="CD8_.CD38.Ki67._FreqParent", fillParam = "Prior.COVID.infection.", groupby="Record.ID", title = "Activated CD8", 
            xLabel = " ", yLabel = "Ki67+CD38+ (% CD8)", repMeasures = F, exponential=F, newform = T)  + scale_y_continuous(breaks=seq(0,10,1), limits = c(0,5.5)) # +
  # ggrepel::geom_text_repel(data = subset(mergedData, Alias =="PHI-071"), aes(label = Label))
# ggsave(filename = "./Images/ActivCD8_bothCohorts_overTime.pdf")

bartlett.test(CD8_.CD38.Ki67._FreqParent ~ timeCategory, data=subset(subsetData, Prior.COVID.infection. == 'No'))
kruskal_test(formula = CD8_.CD38.Ki67._FreqParent ~ timeCategory, data=subset(subsetData, Prior.COVID.infection. == 'No')) 
dunn_test(CD8_.CD38.Ki67._FreqParent ~ timeCategory, data=subset(subsetData, Prior.COVID.infection. == 'No'))
bartlett.test(CD8_.CD38.Ki67._FreqParent ~ timeCategory, data=subset(subsetData, Prior.COVID.infection. == 'Yes'))
kruskal_test(formula = CD8_.CD38.Ki67._FreqParent ~ timeCategory, data=subset(subsetData, Prior.COVID.infection. == 'Yes')) 
dunn_test(CD8_.CD38.Ki67._FreqParent ~ timeCategory, data=subset(subsetData, Prior.COVID.infection. == 'Yes'))


prePostTime(data = subsetData, xData = "timeCategory", yData="CD8_.CD38.Ki67._FreqParent", fillParam = "Prior.COVID.infection.", groupby="Record.ID", title = "Activated CD8", 
            xLabel = " ", yLabel = "Ki67+CD38+ (% CD8)", repMeasures = F, exponential=F, newform = T, recentCOVID = T)  
# ggsave(filename = "./Images/ActivCD8_bothCohorts_overTime_recentCOVID.pdf")




subsetData <- subset(mergedData, Record.ID != "CV-011" & Record.ID != "CV-012" & Record.ID != "CV-005")        # absence of Ki67 stain
subsetData <- subset(subsetData, timeCategory == "Baseline")
twoSampleBar(data = subsetData, xData = "Prior.COVID.infection.", yData = "FCActivCD4_Vax1", fillParam = "Prior.COVID.infection.", title = "Activ CD4", 
             yLabel = "Fold-change at one week", nonparam = T)




#' ------------------ AIM analyses --------------------------
#' 

subsetData <- subset(mergedData,  timeCategory != "two Weeks"& timeCategory != "2 wks post 2nd dose")


##' ******************** CD4 analyses ********************
##'

twoSampleBar(data = subset(subsetData, timeCategory == "Baseline"), xData = "Prior.COVID.infection.", yData = "AIM_CD4.CD69.CD200._FreqParent", 
             fillParam = "Prior.COVID.infection.", title = "Baseline", yLabel = "CD69+ CD200+ (% CD4)", nonparam = T)
# ggsave(filename = "./Images/AIM_CD4_CD69CD200_baseline.pdf", width = 5)

prePostTime(subsetData, xData = "timeCategory", yData="AIM_CD4.CD69.CD200._FreqParent", fillParam = "Prior.COVID.infection.", groupby="Record.ID", title = "CD69+CD200+ CD4", 
            xLabel = " ", yLabel = "CD69+CD200+ (% CD4)", repMeasures = F, exponential=F, newform = T)  + scale_y_continuous(trans="log10", limits = c(0.002, 0.8))
# ggsave(filename = "./Images/AIM_CD4_CD69CD200_overTime.pdf")
bartlett.test(AIM_CD4.CD69.CD200._FreqParent ~ timeCategory, data=subsetData)
kruskal_test(formula = AIM_CD4.CD69.CD200._FreqParent ~ timeCategory, data=subset(subsetData, Prior.COVID.infection. == 'No'))
dunn_test(formula = AIM_CD4.CD69.CD200._FreqParent ~ timeCategory, data=subset(subsetData, Prior.COVID.infection. == 'No'))
tukey_hsd(aov(AIM_CD4.CD69.CD200._FreqParent ~ timeCategory, data=subset(subsetData, Prior.COVID.infection. == 'Yes')))
# kruskal_test(formula = AIM_CD4.CD69.CD200._FreqParent ~ timeCategory, data=subset(subsetData, Prior.COVID.infection. == 'Yes')) 
# dunn_test(formula = AIM_CD4.CD69.CD200._FreqParent ~ timeCategory, data=subset(subsetData, Prior.COVID.infection. == 'Yes'))

subsetData <- subset(mergedData, timeCategory == "Baseline" | timeCategory == "One month post\n2nd dose")
FC_response2 <- dcast( subsetData, `Record.ID`+`Alias` + `Prior.COVID.infection.` ~`timeCategory`, value.var = c("AIM_CD4.CD69.CD200._FreqParent")) 
FC_response2$FoldChange <- FC_response2$`One month post\n2nd dose`/FC_response2$`Baseline`; FC_response2$Cohort <- NULL
FC_response2 <- FC_response2[!is.infinite(FC_response2$FoldChange), ]
FC_response2 %>% group_by(Prior.COVID.infection.) %>% get_summary_stats(type = "common") 
wilcox_test( data = FC_response2, formula = FoldChange ~ Prior.COVID.infection.)

subsetData <- subset(mergedData,  timeCategory != "two Weeks" & timeCategory != "2 wks post 2nd dose")

prePostTime(subsetData, xData = "timeCategory", yData="AIM_CD4.CD69.CD200._FreqParent", fillParam = "Prior.COVID.infection.", groupby="Record.ID", title = "CD69+CD200+ CD4", 
            xLabel = " ", yLabel = "CD69+CD200+ (% CD4)", repMeasures = F, exponential=F, newform = T, recentCOVID = T)  + scale_y_continuous(trans="log10", limits = c(0.002, 0.8))
# ggsave(filename = "./Images/AIM_CD4_CD69CD200_overTime_recentCOVID.pdf")


prePostTime(subsetData, xData = "timeCategory", yData="AIM_CD4.Ox40.CD137._FreqParent", fillParam = "Prior.COVID.infection.", groupby="Record.ID", title = "OX40+CD137+ CD4", 
            xLabel = " ", yLabel = "CD137+OX40+ (% CD4)", repMeasures = F, exponential=F, newform = T)  + scale_y_continuous(trans = "log10")
# ggsave(filename = "./Images/AIM_CD4_CD137Ox40_overTime.pdf")

bartlett.test(AIM_CD4.Ox40.CD137._FreqParent ~ timeCategory, data=subsetData)
kruskal_test(formula = AIM_CD4.Ox40.CD137._FreqParent ~ timeCategory, data=subset(subsetData, Prior.COVID.infection. == 'No'))
dunn_test(formula = AIM_CD4.Ox40.CD137._FreqParent ~ timeCategory, data=subset(subsetData, Prior.COVID.infection. == 'No'))
kruskal_test(formula = AIM_CD4.Ox40.CD137._FreqParent ~ timeCategory, data=subset(subsetData, Prior.COVID.infection. == 'Yes')) 
dunn_test(formula = AIM_CD4.Ox40.CD137._FreqParent ~ timeCategory, data=subset(subsetData, Prior.COVID.infection. == 'Yes'))

##' ******************** CD8 analyses ********************
##' 

twoSampleBar(data = subset(subsetData, timeCategory == "Baseline"), xData = "Prior.COVID.infection.", yData = "AIM_CD8.CD137.IFNg._FreqParent", 
             fillParam = "Prior.COVID.infection.", title = "Baseline", yLabel = "CD137+ IFNg+ (% CD8)", nonparam = T)
# ggsave(filename = "./Images/AIM_CD8_CD137IFNg_baseline.pdf", width = 5)

prePostTime(subsetData, xData = "timeCategory", yData="AIM_CD8.CD137.IFNg._FreqParent", fillParam = "Prior.COVID.infection.", groupby="Record.ID", title = "CD137+IFNg+ CD8", 
            xLabel = " ", yLabel = "CD137+IFNg+ (% CD8)", repMeasures = F, exponential=F, newform = T) + scale_y_continuous(trans = "log10", limits = c(0.0003,0.5))
# ggsave(filename = "./Images/AIM_CD8_CD137IFNg_overTime.pdf")
bartlett.test(AIM_CD8.CD137.IFNg._FreqParent ~ timeCategory, data=subsetData)
kruskal_test(formula = AIM_CD8.CD137.IFNg._FreqParent ~ timeCategory, data=subset(subsetData, Prior.COVID.infection. == 'No'))
dunn_test(formula = AIM_CD8.CD137.IFNg._FreqParent ~ timeCategory, data=subset(subsetData, Prior.COVID.infection. == 'No'))
kruskal_test(formula = AIM_CD8.CD137.IFNg._FreqParent ~ timeCategory, data=subset(subsetData, Prior.COVID.infection. == 'Yes')) 
# dunn_test(formula = AIM_CD8.CD137.IFNg._FreqParent ~ timeCategory, data=subset(subsetData, Prior.COVID.infection. == 'Yes'))

subsetData <- subset(mergedData, timeCategory == "Baseline" | timeCategory == "One month post\n2nd dose")
FC_response2 <- dcast( subsetData, `Record.ID`+`Alias` + `Prior.COVID.infection.` ~`timeCategory`, value.var = c("AIM_CD8.CD137.IFNg._FreqParent")) 
FC_response2$FoldChange <- FC_response2$`One month post\n2nd dose`/FC_response2$`Baseline`; FC_response2$Cohort <- NULL
FC_response2 <- FC_response2[!is.infinite(FC_response2$FoldChange), ]
FC_response2 %>% group_by(Prior.COVID.infection.) %>% get_summary_stats(type = "common") 
wilcox_test( data = FC_response2, formula = FoldChange ~ Prior.COVID.infection.)


subsetData <- subset(mergedData,  timeCategory != "two Weeks" & timeCategory != "2 wks post 2nd dose")

prePostTime(subsetData, xData = "timeCategory", yData="AIM_CD8.CD137.IFNg._FreqParent", fillParam = "Prior.COVID.infection.", groupby="Record.ID", title = "CD137+IFNg+ CD8", 
            xLabel = " ", yLabel = "CD137+IFNg+ (% CD8)", repMeasures = F, exponential=F, newform = T, recentCOVID = T) + scale_y_continuous(trans = "log10", limits = c(0.0003,0.5))
# ggsave(filename = "./Images/AIM_CD8_CD137IFNg_overTime_recentCOVID.pdf")


prePostTime(subsetData, xData = "timeCategory", yData="AIM_CD8.CD69.CD200._FreqParent", fillParam = "Prior.COVID.infection.", groupby="Record.ID", title = "CD69+CD200+ CD8", 
            xLabel = " ", yLabel = "CD69+CD200+ (% CD8)", repMeasures = F, exponential=F, newform = T) 
# ggsave(filename = "./Images/AIM_CD8_CD69CD200_overTime.pdf")

prePostTime(subsetData, xData = "timeCategory", yData="AIM_CD8.Ox40.CD137._FreqParent", fillParam = "Prior.COVID.infection.", groupby="Record.ID", title = "OX40+CD137+ CD8", 
            xLabel = " ", yLabel = "CD137+OX40+  (% CD8)", repMeasures = F, exponential=F, newform = T) + scale_y_continuous(trans = "log10")
# ggsave(filename = "./Images/AIM_CD8_CD137Ox40_overTime.pdf")
bartlett.test(AIM_CD8.Ox40.CD137._FreqParent ~ timeCategory, data=subsetData)
kruskal_test(formula = AIM_CD8.Ox40.CD137._FreqParent ~ timeCategory, data=subset(subsetData, Prior.COVID.infection. == 'No'))
dunn_test(formula = AIM_CD8.Ox40.CD137._FreqParent ~ timeCategory, data=subset(subsetData, Prior.COVID.infection. == 'No'))
kruskal_test(formula = AIM_CD8.Ox40.CD137._FreqParent ~ timeCategory, data=subset(subsetData, Prior.COVID.infection. == 'Yes')) 
dunn_test(formula = AIM_CD8.Ox40.CD137._FreqParent ~ timeCategory, data=subset(subsetData, Prior.COVID.infection. == 'Yes'))


prePostTime(subsetData, xData = "timeCategory", yData="AIM_CD4.TNF._FreqParent", fillParam = "Prior.COVID.infection.", groupby="Record.ID", title = "TNF+ CD4", 
            xLabel = " ", yLabel = "TNF+ (% CD4)", repMeasures = F, exponential=F, newform = T)  + scale_y_continuous(trans="log10")
# ggsave(filename = "./Images/AIM_CD4_TNF_overTime.pdf")

bartlett.test(AIM_CD4.TNF._FreqParent ~ timeCategory, data=subsetData)
kruskal_test(formula = AIM_CD4.TNF._FreqParent ~ timeCategory, data=subset(subsetData, Prior.COVID.infection. == 'No'))
dunn_test(formula = AIM_CD4.TNF._FreqParent ~ timeCategory, data=subset(subsetData, Prior.COVID.infection. == 'No'))
kruskal_test(formula = AIM_CD4.TNF._FreqParent ~ timeCategory, data=subset(subsetData, Prior.COVID.infection. == 'Yes')) 
dunn_test(formula = AIM_CD4.TNF._FreqParent ~ timeCategory, data=subset(subsetData, Prior.COVID.infection. == 'Yes'))

prePostTime(subsetData, xData = "timeCategory", yData="AIM_CD8.TNF._FreqParent", fillParam = "Prior.COVID.infection.", groupby="Record.ID", title = "TNF+ CD8", 
            xLabel = " ", yLabel = "TNF+ (% CD8)", repMeasures = F, exponential=F, newform = T)  + scale_y_continuous(trans="log10")
# ggsave(filename = "./Images/AIM_CD8_TNF_overTime.pdf")

bartlett.test(AIM_CD8.TNF._FreqParent ~ timeCategory, data=subsetData)
kruskal_test(formula = AIM_CD8.TNF._FreqParent ~ timeCategory, data=subset(subsetData, Prior.COVID.infection. == 'No'))
dunn_test(formula = AIM_CD8.TNF._FreqParent ~ timeCategory, data=subset(subsetData, Prior.COVID.infection. == 'No'))
kruskal_test(formula = AIM_CD8.TNF._FreqParent ~ timeCategory, data=subset(subsetData, Prior.COVID.infection. == 'Yes')) 
dunn_test(formula = AIM_CD8.TNF._FreqParent ~ timeCategory, data=subset(subsetData, Prior.COVID.infection. == 'Yes'))


prePostTime(subsetData, xData = "timeCategory", yData="AIM_CD4.IFNg._FreqParent", fillParam = "Prior.COVID.infection.", groupby="Record.ID", title = "IFNg+ CD4", 
            xLabel = " ", yLabel = "IFNg+ (% CD4)", repMeasures = F, exponential=F, newform = T)  + scale_y_continuous(trans="log10")
# ggsave(filename = "./Images/AIM_CD4_IFNg_overTime.pdf")

bartlett.test(AIM_CD4.IFNg._FreqParent ~ timeCategory, data=subsetData)
kruskal_test(formula = AIM_CD4.IFNg._FreqParent ~ timeCategory, data=subset(subsetData, Prior.COVID.infection. == 'No'))
dunn_test(formula = AIM_CD4.IFNg._FreqParent ~ timeCategory, data=subset(subsetData, Prior.COVID.infection. == 'No'))
kruskal_test(formula = AIM_CD4.IFNg._FreqParent ~ timeCategory, data=subset(subsetData, Prior.COVID.infection. == 'Yes')) 
dunn_test(formula = AIM_CD4.IFNg._FreqParent ~ timeCategory, data=subset(subsetData, Prior.COVID.infection. == 'Yes'))


prePostTime(subsetData, xData = "timeCategory", yData="AIM_CD8.IFNg._FreqParent", fillParam = "Prior.COVID.infection.", groupby="Record.ID", title = "IFNg+ CD8", 
            xLabel = " ", yLabel = "IFNg+ (% CD8)", repMeasures = F, exponential=F, newform = T) + scale_y_continuous(trans = "log10")
# ggsave(filename = "./Images/AIM_CD8_IFNg_overTime.pdf")

bartlett.test(AIM_CD8.IFNg._FreqParent ~ timeCategory, data=subsetData)
kruskal_test(formula = AIM_CD8.IFNg._FreqParent ~ timeCategory, data=subset(subsetData, Prior.COVID.infection. == 'No'))
dunn_test(formula = AIM_CD8.IFNg._FreqParent ~ timeCategory, data=subset(subsetData, Prior.COVID.infection. == 'No'))
kruskal_test(formula = AIM_CD8.IFNg._FreqParent ~ timeCategory, data=subset(subsetData, Prior.COVID.infection. == 'Yes')) 
dunn_test(formula = AIM_CD8.IFNg._FreqParent ~ timeCategory, data=subset(subsetData, Prior.COVID.infection. == 'Yes'))

twoSampleBar(data = subset(subsetData, timeCategory == "Baseline"), xData = "Prior.COVID.infection.", yData = "AIM_CD8.IFNg._FreqParent", 
             fillParam = "Prior.COVID.infection.", title = "Baseline", yLabel = "CD69+ CD200+ (% CD4)", nonparam = T)
# ggsave(filename = "./Images/AIM_CD4_CD69CD200_baseline.pdf", width = 5)

subsetData <- subset(mergedData, timeCategory == "Baseline" | timeCategory == "One month post\n2nd dose")
FC_response2 <- dcast( subsetData, `Record.ID`+`Alias` + `Prior.COVID.infection.` ~`timeCategory`, value.var = c("AIM_CD8.IFNg._FreqParent")) 
FC_response2$FC_AIMCD8_IFNg <- FC_response2$`One month post\n2nd dose`/FC_response2$`Baseline`; FC_response2$Cohort <- NULL
FC_response2 %>% group_by(Prior.COVID.infection.) %>% get_summary_stats(type = "common") 



subsetData <- subset(mergedData,  timeCategory != "two Weeks" & timeCategory != "2 wks post 2nd dose")


dunn_test(AIM_CD4.CD69.CD137._FreqParent ~ timeCategory, data=subset(subsetData, Prior.COVID.infection. == 'No')) 
dunn_test(AIM_CD4.CD71.CD137._FreqParent ~ timeCategory, data=subset(subsetData, Prior.COVID.infection. == 'No')) 
dunn_test(AIM_CD4.CD40L._FreqParent ~ timeCategory, data=subset(subsetData, Prior.COVID.infection. == 'No'))  
dunn_test(AIM_CD8.CD40L._FreqParent ~ timeCategory, data=subset(subsetData, Prior.COVID.infection. == 'No'))  
dunn_test(AIM_CD8.TNF._FreqParent ~ timeCategory, data=subset(subsetData, Prior.COVID.infection. == 'No'))  
dunn_test(AIM_CD8.TNF._FreqParent ~ timeCategory, data=subset(subsetData, Prior.COVID.infection. == 'No'))  

bivScatter(data1 = subset(mergedData, Prior.COVID.infection. == 'No' & timeCategory == "Post 2nd dose" & Tube == "HEP"), name1 = "Naive", 
           data2 = subset(mergedData, Prior.COVID.infection. == 'Yes' & timeCategory == "Post 2nd dose" & Tube == "HEP"), name2 = "Experienced", 
           xData = "AIM_CD4.CD69.CD200._FreqParent", yData = "FCActivCD8_Vax1", fillParam = "Prior.COVID.infection.", title = "Vaccine dose 2", 
           xLabel = "FC_AIMCD4_69200", yLabel = "Fold-change ICOS+CD38+ cTfh", nonparam = T) #+   scale_y_continuous(breaks=seq(0,30,1), limits = c(0,6) ) + scale_x_continuous(limits = c(20,70))
# ggsave(filename = "./Images/Age_correl_FCtfh_Vax2.pdf", width=8)

  
  
subsetData <- subset(mergedData, Prior.COVID.infection. == 'No' & timeCategory == "Post 2nd dose")
subsetData <- subsetData[,grep(paste(c("^FC", "Age"), collapse = "|"), names(subsetData))] 
temp <- names(subsetData) ; temp <-  do.call(rbind.data.frame, strsplit(temp, split = "_"))
# subsetData <- subsetData[, -grep("SampleID",names(subsetData))]
subsetData <- subsetData[, -grep("Btet", names(subsetData))]
cor.AIM <- cor(subsetData, method="kendall" , use="pairwise.complete.obs" )
cor.AIM.pmat <- ggcorrplot::cor_pmat(subsetData, method="kendall", use="pairwise.complete.obs"  )
ggcorrplot::ggcorrplot(corr = cor.AIM, p.mat = cor.AIM.pmat, title = "AIM post 2nd dose", legend.title = "Kendall tau", tl.cex = 10,pch.cex = 1, insig = "blank")
# ggsave(filename = "./Images/AIM_ggcorrplot_Naive_post2nd.pdf")
ggcorrplot::ggcorrplot(corr = cor.AIM, p.mat = cor.AIM.pmat, title = "AIM post 2nd dose", legend.title = "Kendall tau", tl.cex = 18,pch.cex = 1) #, insig = "blank"
# ggsave(filename = "./Images/AIM_ggcorrplot_Naive_post2nd_full.pdf")

subsetData <- subset(mergedData, Prior.COVID.infection. == 'Yes' & timeCategory == "Post 2nd dose")
subsetData <- subsetData[,grep("^AIM_", names(subsetData))] 
temp <- names(subsetData) ; temp <-  do.call(rbind.data.frame, strsplit(temp, split = "_"))
subsetData <- subsetData[, -grep("SampleID",names(subsetData))]
cor.AIM <- cor(subsetData, method="kendall" , use="pairwise.complete.obs" )
cor.AIM.pmat <- ggcorrplot::cor_pmat(subsetData, method="kendall", use="pairwise.complete.obs"  )
ggcorrplot::ggcorrplot(corr = cor.AIM, p.mat = cor.AIM.pmat, title = "AIM post 2nd dose", legend.title = "Kendall tau", tl.cex = 18,pch.cex = 1, insig = "blank")
# ggsave(filename = "./Images/AIM_ggcorrplot_Exp_post2nd.pdf")
ggcorrplot::ggcorrplot(corr = cor.AIM, p.mat = cor.AIM.pmat, title = "AIM post 2nd dose", legend.title = "Kendall tau", tl.cex = 18,pch.cex = 1) #, insig = "blank"
# ggsave(filename = "./Images/AIM_ggcorrplot_Exp_post2nd_full.pdf")



#' ------------------ B tetramer analyses --------------------------
#' 

subsetData <- subset(mergedData,  timeCategory == "Baseline")
twoSampleBar(data = subsetData, xData = "Prior.COVID.infection.", yData = "Btet_RBD_FreqParent", fillParam = "Prior.COVID.infection.", title = "Baseline", 
             yLabel = "RBD+ (% CD19)", nonparam = T)
# ggsave(filename = "./Images/Btet_frequency_baseline.pdf", width=5)

twoSampleBar(data = subsetData, xData = "Prior.COVID.infection.", yData = "Btet_RBDhiIgDloIgGhi", fillParam = "Prior.COVID.infection.", title = "Baseline", 
             yLabel = "IgG+ IgD- (% RBD)", nonparam = T)
# ggsave(filename = "./Images/Btet_IgGhi_baseline.pdf", width=5)


subsetData <- subset(mergedData,  timeCategory != "two Weeks"& timeCategory != "2 wks post 2nd dose")

prePostTime(subsetData, xData = "timeCategory", yData="Btet_RBD_FreqParent", fillParam = "Prior.COVID.infection.", groupby="Record.ID", title = "RBD-reactive B cells", 
            xLabel = " ", yLabel = "RBD+ (% CD19)", repMeasures = F, exponential=F, newform = T) + scale_y_continuous(trans = "log10", limits = c(0.02, 5))
# ggsave(filename = "./Images/Btet_frequency_overTime.pdf")

bartlett.test(Btet_RBD_FreqParent ~ timeCategory, data=subset(subsetData, Prior.COVID.infection. == 'No'))
kruskal_test(formula = Btet_RBD_FreqParent ~ timeCategory, data=subset(subsetData, Prior.COVID.infection. == 'No')) 
dunn_test(Btet_RBD_FreqParent ~ timeCategory, data=subset(subsetData, Prior.COVID.infection. == 'No'))
bartlett.test(Btet_RBD_FreqParent ~ timeCategory, data=subset(subsetData, Prior.COVID.infection. == 'Yes'))
kruskal_test(formula = Btet_RBD_FreqParent ~ timeCategory, data=subset(subsetData, Prior.COVID.infection. == 'Yes')) 
dunn_test(Btet_RBD_FreqParent ~ timeCategory, data=subset(subsetData, Prior.COVID.infection. == 'Yes'))


subsetData <- subset(mergedData, timeCategory == "Baseline" | timeCategory == "Post 2nd dose")
FC_response2 <- dcast( subsetData, `Record.ID`+`Alias` + `Prior.COVID.infection.` ~`timeCategory`, value.var = c("Btet_RBD_FreqParent")) 
FC_response2$FoldChange <- FC_response2$`Post 2nd dose`/FC_response2$`Baseline`; FC_response2$Cohort <- NULL
FC_response2 <- FC_response2[!is.infinite(FC_response2$FoldChange), ]
FC_response2 %>% group_by(Prior.COVID.infection.) %>% get_summary_stats(type = "common") 
FC_response2 %>% wilcox_test(FoldChange ~ Prior.COVID.infection.)


subsetData <- subset(mergedData,  timeCategory != "two Weeks"& timeCategory != "2 wks post 2nd dose");
prePostTime(subsetData, xData = "timeCategory", yData="Btet_RBD_FreqParent", fillParam = "Prior.COVID.infection.", groupby="Record.ID", title = "RBD-reactive B cells", 
            xLabel = " ", yLabel = "RBD+ (% CD19)", repMeasures = F, exponential=F, newform = T, recentCOVID = T) + scale_y_continuous(trans = "log10", limits = c(0.02, 5))
# ggsave(filename = "./Images/Btet_frequency_overTime_recentCOVID.pdf")



prePostTime(subsetData, xData = "timeCategory", yData="Btet_RBDhiIgDloIgGhi", fillParam = "Prior.COVID.infection.", groupby="Record.ID", title = "Switched RBD+ B cells", 
            xLabel = " ", yLabel = "IgG+ IgD- (% RBD+)", repMeasures = F, exponential=F, newform = T) + scale_y_continuous(breaks=seq(0,100,25), limits = c(0,125))
# ggsave(filename = "./Images/Btet_IgGhi_overTime.pdf")

bartlett.test(Btet_RBDhiIgDloIgGhi ~ timeCategory, data=subset(subsetData, Prior.COVID.infection. == 'No'))
kruskal_test(formula = Btet_RBDhiIgDloIgGhi ~ timeCategory, data=subset(subsetData, Prior.COVID.infection. == 'No')) 
dunn_test(Btet_RBDhiIgDloIgGhi ~ timeCategory, data=subset(subsetData, Prior.COVID.infection. == 'No'))
bartlett.test(Btet_RBDhiIgDloIgGhi ~ timeCategory, data=subset(subsetData, Prior.COVID.infection. == 'Yes'))
tukey_hsd( aov(Btet_RBDhiIgDloIgGhi ~ timeCategory, data = subset(subsetData, Prior.COVID.infection. == 'Yes')))
# kruskal_test(formula = Btet_RBDhiIgDloIgGhi ~ timeCategory, data=subset(subsetData, Prior.COVID.infection. == 'Yes')) 
# dunn_test(Btet_RBDhiIgDloIgGhi ~ timeCategory, data=subset(subsetData, Prior.COVID.infection. == 'Yes'))

prePostTime(subsetData, xData = "timeCategory", yData="Btet_RBDhiIgDloIgGhi", fillParam = "Prior.COVID.infection.", groupby="Record.ID", title = "Switched RBD+ B cells", 
            xLabel = " ", yLabel = "IgG+ IgD- (% RBD+)", repMeasures = F, exponential=F, newform = T, recentCOVID = T) + scale_y_continuous(breaks=seq(0,100,25), limits = c(0,125))
# ggsave(filename = "./Images/Btet_IgGhi_overTime_recentCOVID.pdf")


prePostTime(subsetData, xData = "timeCategory", yData="Btet_RBDhi_DN2_FreqCD19", fillParam = "Prior.COVID.infection.", groupby="Record.ID", title = "DN2 RBD+ B cells", 
            xLabel = " ", yLabel = "DN2 B cells (% CD19+)", repMeasures = F, exponential=F, newform = T) + scale_y_continuous(limits = c(0.001,3),  trans = "log10")
# ggsave(filename = "./Images/Btet_DN2_overTime_fullPlot.pdf")

bartlett.test(Btet_RBDhi_DN2_FreqCD19 ~ timeCategory, data=subset(subsetData, Prior.COVID.infection. == 'No'))
kruskal_test(formula = Btet_RBDhi_DN2_FreqCD19 ~ timeCategory, data=subset(subsetData, Prior.COVID.infection. == 'No')) 
dunn_test(Btet_RBDhi_DN2_FreqCD19 ~ timeCategory, data=subset(subsetData, Prior.COVID.infection. == 'No'))
bartlett.test(Btet_RBDhi_DN2_FreqCD19 ~ timeCategory, data=subset(subsetData, Prior.COVID.infection. == 'Yes'))
kruskal_test(formula = Btet_RBDhi_DN2_FreqCD19 ~ timeCategory, data=subset(subsetData, Prior.COVID.infection. == 'Yes')) 
dunn_test(Btet_RBDhi_DN2_FreqCD19 ~ timeCategory, data=subset(subsetData, Prior.COVID.infection. == 'Yes'))


prePostTime(subsetData, xData = "timeCategory", yData="Btet_RBD_CD71._FreqParent", fillParam = "Prior.COVID.infection.", groupby="Record.ID", title = "CD71+ RBD+ B cells", 
            xLabel = " ", yLabel = "CD71+ IgDlo (% RBD+)", repMeasures = F, exponential=F, newform = T) + scale_y_continuous(breaks=seq(0,100,20),  trans = "identity", limits = c(0, 100))
# ggsave(filename = "./Images/Btet_CD71_overTime.pdf")
bartlett.test(Btet_RBD_CD71._FreqParent ~ timeCategory, data=subset(subsetData, Prior.COVID.infection. == 'No'))
kruskal_test(formula = Btet_RBD_CD71._FreqParent ~ timeCategory, data=subset(subsetData, Prior.COVID.infection. == 'No')) 
dunn_test(Btet_RBD_CD71._FreqParent ~ timeCategory, data=subset(subsetData, Prior.COVID.infection. == 'No'))
bartlett.test(Btet_RBD_CD71._FreqParent ~ timeCategory, data=subset(subsetData, Prior.COVID.infection. == 'Yes'))
kruskal_test(formula = Btet_RBD_CD71._FreqParent ~ timeCategory, data=subset(subsetData, Prior.COVID.infection. == 'Yes')) 
dunn_test(Btet_RBD_CD71._FreqParent ~ timeCategory, data=subset(subsetData, Prior.COVID.infection. == 'Yes'))



prePostTime(subsetData, xData = "timeCategory", yData="Btet_RBD_CD21lo_FreqParent", fillParam = "Prior.COVID.infection.", groupby="Record.ID", title = "CD21lo RBD+ B cells", 
            xLabel = " ", yLabel = "CD21lo (% RBD+)", repMeasures = F, exponential=F, newform = T) + scale_y_continuous(trans="pseudo_log", limits = c(0,100), breaks = c(0,1,10,25,50,100))
  # scale_y_continuous(breaks=seq(0,100,10),  trans = "identity", limits = c(0,60))
# ggsave(filename = "./Images/Btet_CD21lo_overTime.pdf")
bartlett.test(Btet_RBD_CD21lo_FreqParent ~ timeCategory, data=subset(subsetData, Prior.COVID.infection. == 'No'))
kruskal_test(formula = Btet_RBD_CD21lo_FreqParent ~ timeCategory, data=subset(subsetData, Prior.COVID.infection. == 'No')) 
dunn_test(Btet_RBD_CD21lo_FreqParent ~ timeCategory, data=subset(subsetData, Prior.COVID.infection. == 'No'))
bartlett.test(Btet_RBD_CD21lo_FreqParent ~ timeCategory, data=subset(subsetData, Prior.COVID.infection. == 'Yes'))
kruskal_test(formula = Btet_RBD_CD21lo_FreqParent ~ timeCategory, data=subset(subsetData, Prior.COVID.infection. == 'Yes')) 
dunn_test(Btet_RBD_CD21lo_FreqParent ~ timeCategory, data=subset(subsetData, Prior.COVID.infection. == 'Yes'))


prePostTime(subsetData, xData = "timeCategory", yData="Btet_RBD_CD24lo_FreqParent", fillParam = "Prior.COVID.infection.", groupby="Record.ID", title = "CD24lo RBD+ B cells", 
            xLabel = " ", yLabel = "CD24lo (% RBD+)", repMeasures = F, exponential=F, newform = T) + scale_y_continuous(breaks=seq(0,100,20),  trans = "identity", limits = c(0, 120))
# ggsave(filename = "./Images/Btet_CD24lo_overTime.pdf")
bartlett.test(Btet_RBD_CD24lo_FreqParent ~ timeCategory, data=subset(subsetData, Prior.COVID.infection. == 'No'))
kruskal_test(formula = Btet_RBD_CD24lo_FreqParent ~ timeCategory, data=subset(subsetData, Prior.COVID.infection. == 'No')) 
dunn_test(Btet_RBD_CD24lo_FreqParent ~ timeCategory, data=subset(subsetData, Prior.COVID.infection. == 'No'))
bartlett.test(Btet_RBD_CD24lo_FreqParent ~ timeCategory, data=subset(subsetData, Prior.COVID.infection. == 'Yes'))
tukey_hsd( aov(Btet_RBD_CD24lo_FreqParent ~ timeCategory, data=subset(subsetData, Prior.COVID.infection. == 'Yes')) )

subsetData$Btet_RBDhi_mfiRBD <- str_replace(pattern = "n/a", replacement = "", string = subsetData$Btet_RBDhi_mfiRBD)
subsetData$Btet_RBDhi_mfiRBD <- as.numeric(subsetData$Btet_RBDhi_mfiRBD)
prePostTime(subsetData, xData = "timeCategory", yData="Btet_RBDhi_mfiRBD", fillParam = "Prior.COVID.infection.", groupby="Record.ID", title = "RBD+ B cells", 
            xLabel = " ", yLabel = "MFI RBD", repMeasures = F, exponential=F, newform = T) + scale_y_continuous( trans = "log10")
# ggsave(filename = "./Images/Btet_9_overTime.pdf")
bartlett.test(Btet_RBDhi_mfiRBD ~ timeCategory, data=subset(subsetData, Prior.COVID.infection. == 'No'))
kruskal_test(formula = Btet_RBDhi_mfiRBD ~ timeCategory, data=subset(subsetData, Prior.COVID.infection. == 'No')) 
dunn_test(Btet_RBDhi_mfiRBD ~ timeCategory, data=subset(subsetData, Prior.COVID.infection. == 'No'))
bartlett.test(Btet_RBDhi_mfiRBD ~ timeCategory, data=subset(subsetData, Prior.COVID.infection. == 'Yes'))
tukey_hsd( aov(Btet_RBDhi_mfiRBD ~ timeCategory, data=subset(subsetData, Prior.COVID.infection. == 'Yes')) )


subsetData$Btet_RBDhiIgDloIgGhi_mfiRBD <- str_replace(pattern = "n/a", replacement = "", string = subsetData$Btet_RBDhiIgDloIgGhi_mfiRBD)
subsetData$Btet_RBDhiIgDloIgGhi_mfiRBD <- as.numeric(subsetData$Btet_RBDhiIgDloIgGhi_mfiRBD)
prePostTime(subsetData, xData = "timeCategory", yData="Btet_RBDhiIgDloIgGhi_mfiRBD", fillParam = "Prior.COVID.infection.", groupby="Record.ID", title = "RBD+ B cells", 
            xLabel = " ", yLabel = "MFI RBD", repMeasures = F, exponential=F, newform = T) + scale_y_continuous( trans = "log10")
# ggsave(filename = "./Images/Btet_0_overTime.pdf")
# bartlett.test(Btet_RBDhiIgDloIgGhi_mfiRBD ~ timeCategory, data=subset(subsetData, Prior.COVID.infection. == 'No'))
kruskal_test(formula = Btet_RBDhiIgDloIgGhi_mfiRBD ~ timeCategory, data=subset(subsetData, Prior.COVID.infection. == 'No')) 
dunn_test(Btet_RBDhiIgDloIgGhi_mfiRBD ~ timeCategory, data=subset(subsetData, Prior.COVID.infection. == 'No'))
bartlett.test(Btet_RBDhiIgDloIgGhi_mfiRBD ~ timeCategory, data=subset(subsetData, Prior.COVID.infection. == 'Yes'))
tukey_hsd( aov(Btet_RBDhiIgDloIgGhi_mfiRBD ~ timeCategory, data=subset(subsetData, Prior.COVID.infection. == 'Yes')) )






#' ------------------ Tfh analyses --------------------------
#' 
# subsetData <- mergedData[which(mergedData$shortForm == "bL" | mergedData$shortForm == "oW"),]
# prePostTime(data = subsetData, xData = "shortForm", yData="CD4_.Nonnaive.cTfh.ICOS..CD38.._FreqParent", fillParam = "Prior.COVID.infection.", groupby="Record.ID", 
# title = "cTfh", xLabel = "Prior COVID?", yLabel = "ICOS+CD38+ (% cTfh)") 

subsetData <- subset(mergedData,  timeCategory != "two Weeks" & timeCategory != "2 wks post 2nd dose");   
subsetData$timeCategory <- factor(subsetData$timeCategory, levels = c("Baseline", "Post 1st dose", "Pre 2nd dose", "Post 2nd dose","One month post\n2nd dose"))

prePostTime(subsetData, xData = "DPV", yData="CD4_.Nonnaive.cTfh.ICOS..CD38.._FreqParent", fillParam = "Prior.COVID.infection.", groupby="Record.ID", title = "Tfh",
            xLabel = "Prior COVID?", yLabel = "ICOS+CD38+ (% cTfh)", repMeasures = F, exponential=F) +  geom_vline(xintercept = 0,linetype="dashed", alpha=0.5) +  coord_cartesian(xlim = c(-1,12)) 

subsetData$DPV <- subsetData$DPV - difftime(subsetData$Vaccine.2.date, subsetData$Vaccine.1.date, units="days" ) 
prePostTime(subsetData, xData = "DPV", yData="CD4_.Nonnaive.cTfh.ICOS..CD38.._FreqParent", fillParam = "Prior.COVID.infection.", groupby="Record.ID", title = "Tfh",
            xLabel = "Prior COVID?", yLabel = "ICOS+CD38+ (% cTfh)", repMeasures = F, exponential=F) +  geom_vline(xintercept = 0,linetype="dashed", alpha=0.5) +  coord_cartesian(xlim = c(-5,12)) 


prePostTime(subsetData, xData = "timeCategory", yData="CD4_.Nonnaive.cTfh.ICOS..CD38.._FreqParent", fillParam = "Prior.COVID.infection.", groupby="Record.ID", title = "cTfh", 
            xLabel = " ", yLabel = "ICOS+CD38+ (% cTfh)", repMeasures = F, newform=T)  
# ggsave(filename = "./Images/cTfh_responses_bothCohorts_overTime.pdf")
bartlett.test(CD4_.Nonnaive.cTfh.ICOS..CD38.._FreqParent ~ timeCategory, data=subset(subsetData, Prior.COVID.infection. == 'No'))
tukey_hsd(aov(CD4_.Nonnaive.cTfh.ICOS..CD38.._FreqParent ~ timeCategory, data=subset(subsetData, Prior.COVID.infection. == 'No')))
bartlett.test(CD4_.Nonnaive.cTfh.ICOS..CD38.._FreqParent ~ timeCategory, data=subset(subsetData, Prior.COVID.infection. == 'Yes'))
tukey_hsd( aov(CD4_.Nonnaive.cTfh.ICOS..CD38.._FreqParent ~ timeCategory, data=subset(subsetData, Prior.COVID.infection. == 'Yes')) )


prePostTime(subsetData, xData = "timeCategory", yData="CD4_.Nonnaive.cTfh.ICOS..CD38.._FreqParent", fillParam = "Prior.COVID.infection.", groupby="Record.ID", title = "cTfh", 
            xLabel = " ", yLabel = "ICOS+CD38+ (% cTfh)", repMeasures = F, newform=T, recentCOVID = T)  
# ggsave(filename = "./Images/cTfh_responses_bothCohorts_overTime_recentCOVID.pdf")


twoSampleBar(data = subset(mergedData, timeCategory == "Post 1st dose"), xData = "Prior.COVID.infection.", yData = "FCtfh_Vax1", fillParam = "Prior.COVID.infection.", title = "ICOS+CD38+ cTfh", 
             yLabel = "Fold-change at one week")



subsetData = subset(mergedData, mergedData$shortForm == 'oW')
subsetData %>% group_by(Prior.COVID.infection.) %>% shapiro_test(CD4_.Nonnaive.cTfh.ICOS..CD38.._FreqParent)
twoSampleBar(data=subsetData, xData = "Prior.COVID.infection.", yData="CD4_.Nonnaive.cTfh.ICOS..CD38.._FreqParent", fillParam = "Prior.COVID.infection.", title = "cTfh at oneWeek", 
             yLabel = "ICOS+CD38+ (% cTfh)", nonparam = T) + scale_y_continuous(breaks = seq(0,21,3), limits = c(0,20)) 


subsetData <- subset(mergedData,  timeCategory != "two Weeks" & timeCategory != "2 wks post 2nd dose"); 
subsetData$timeCategory <- factor(subsetData$timeCategory, levels = c("Baseline", "Post 1st dose", "Pre 2nd dose", "Post 2nd dose","One month post\n2nd dose"))
prePostTime(data = subsetData, xData = "timeCategory", yData="CD4_.Nonnaive.cTfh.ICOS..CD38...CXCR3._FreqParent", fillParam = "Prior.COVID.infection.", groupby="Record.ID", title = "cTfh", 
            xLabel = " ", yLabel = "CXCR3+ (% ICOS+CD38+ cTfh)", repMeasures = F, newform =  T)  + scale_y_continuous(limits = c(0,80),breaks=seq(0,100,10)) 
# ggsave(filename = "./Images/cTfh_responses_CXCR3_bothCohorts_overTime.pdf")
bartlett.test(CD4_.Nonnaive.cTfh.ICOS..CD38...CXCR3._FreqParent ~ timeCategory, data=subset(subsetData, Prior.COVID.infection. == 'No'))
tukey_hsd(aov(CD4_.Nonnaive.cTfh.ICOS..CD38...CXCR3._FreqParent ~ timeCategory, data=subset(subsetData, Prior.COVID.infection. == 'No')))
bartlett.test(CD4_.Nonnaive.cTfh.ICOS..CD38...CXCR3._FreqParent ~ timeCategory, data=subset(subsetData, Prior.COVID.infection. == 'Yes'))
# tukey_hsd( aov(CD4_.Nonnaive.cTfh.ICOS..CD38...CXCR3._FreqParent ~ timeCategory, data=subset(subsetData, Prior.COVID.infection. == 'Yes')) )
kruskal_test(formula = CD4_.Nonnaive.cTfh.ICOS..CD38...CXCR3._FreqParent ~ timeCategory, data=subset(subsetData, Prior.COVID.infection. == 'Yes')) 
# dunn_test(CD4_.Nonnaive.cTfh.ICOS..CD38...CXCR3._FreqParent ~ timeCategory, data=subset(subsetData, Prior.COVID.infection. == 'Yes'))

prePostTime(data = subsetData, xData = "timeCategory", yData="CD4_.Nonnaive.cTfh.ICOS..CD38...CXCR3._FreqParent", fillParam = "Prior.COVID.infection.", groupby="Record.ID", title = "cTfh", 
            xLabel = " ", yLabel = "CXCR3+ (% ICOS+CD38+ cTfh)", repMeasures = F, newform =  T, recentCOVID = T)  + scale_y_continuous(limits = c(0,80),breaks=seq(0,100,10)) 
# ggsave(filename = "./Images/cTfh_responses_CXCR3_bothCohorts_overTime_recentCOVID.pdf")




bivScatter(data1 = subset(mergedData, Prior.COVID.infection. == 'No' & timeCategory == "Post 1st dose" & Tube == "HEP"), name1 = "Naive", 
           data2 = subset(mergedData, Prior.COVID.infection. == 'Yes' & timeCategory == "Post 1st dose" & Tube == "HEP"), name2 = "Experienced", 
           xData = "Age", yData = "CD4_.Nonnaive.cTfh.ICOS..CD38.._FreqParent", fillParam = "Prior.COVID.infection.", title = "Post 1st dose", 
           xLabel = "Age", yLabel = "ICOS+CD38+ (% cTfh)", nonparam = T) +   scale_y_continuous(breaks=seq(0,27,3), limits = c(0,25) )+ scale_x_continuous(limits = c(20,70))
# ggsave(filename = "./Images/Age_correl_HiHiTfh_Vax1.pdf", width=8)
bivScatter(data1 = subset(mergedData, Prior.COVID.infection. == 'No' & timeCategory == "Post 2nd dose" & Tube == "HEP"), name1 = "Naive", 
           data2 = subset(mergedData, Prior.COVID.infection. == 'Yes' & timeCategory == "Post 2nd dose" & Tube == "HEP"), name2 = "Experienced", 
           xData = "Age", yData = "CD4_.Nonnaive.cTfh.ICOS..CD38.._FreqParent", fillParam = "Prior.COVID.infection.", title = "Post 2nd dose", 
           xLabel = "Age", yLabel = "ICOS+CD38+ (% cTfh)", nonparam = T) +   scale_y_continuous(breaks=seq(0,27,3), limits = c(0,25) ) + scale_x_continuous(limits = c(20,70))
# ggsave(filename = "./Images/Age_correl_HiHiTfh_Vax2.pdf", width=8)

bivScatter(data1 = subset(mergedData, Prior.COVID.infection. == 'No' & timeCategory == "Post 1st dose" & Tube == "HEP"), name1 = "Naive", 
           data2 = subset(mergedData, Prior.COVID.infection. == 'Yes' & timeCategory == "Post 1st dose" & Tube == "HEP"), name2 = "Experienced", 
           xData = "Age", yData = "FCtfh_Vax1", fillParam = "Prior.COVID.infection.", title = "Post 1st dose", 
           xLabel = "Age", yLabel = "Fold-change ICOS+CD38+ cTfh", nonparam = T) +scale_y_continuous(breaks=seq(0,30,1), limits = c(0,6) )+ scale_x_continuous(limits = c(20,70))
# ggsave(filename = "./Images/Age_correl_FCtfh_Vax1.pdf", width=8)
bivScatter(data1 = subset(mergedData, Prior.COVID.infection. == 'No' & timeCategory == "Post 2nd dose" & Tube == "HEP"), name1 = "Naive", 
           data2 = subset(mergedData, Prior.COVID.infection. == 'Yes' & timeCategory == "Post 2nd dose" & Tube == "HEP"), name2 = "Experienced", 
           xData = "Age", yData = "FCtfh_Vax2", fillParam = "Prior.COVID.infection.", title = "Post 2nd dose", 
           xLabel = "Age", yLabel = "Fold-change ICOS+CD38+ cTfh", nonparam = T) +   scale_y_continuous(breaks=seq(0,30,1), limits = c(0,6) ) + scale_x_continuous(limits = c(20,70))
# ggsave(filename = "./Images/Age_correl_FCtfh_Vax2.pdf", width=8)


subsetData <- subset(mergedData,  timeCategory != "two Weeks" & timeCategory != "2 wks post 2nd dose");   
subsetData$timeCategory <- factor(subsetData$timeCategory, levels = c("Baseline", "Post 1st dose", "Pre 2nd dose", "Post 2nd dose","One month post\n2nd dose"))

prePostTime(subsetData, xData = "timeCategory", yData="AIM_CD4.CXCR5.PD1..CD38hi.CD69.CD200._FreqParent", fillParam = "Prior.COVID.infection.", groupby="Record.ID", 
            title = "Activated cTfh", xLabel = " ", yLabel = "CD69+CD200+ (% activated cTfh)", repMeasures = F, exponential=F, newform = T) + 
  scale_y_continuous(limits = c(0,10), breaks = seq(0,10,1))
# ggsave(filename = "./Images/AIM_cTfh_69-200_overTime.pdf")
bartlett.test(AIM_CD4.CXCR5.PD1..CD38hi.CD69.CD200._FreqParent ~ timeCategory, data=subset(subsetData, Prior.COVID.infection. == 'No'))
kruskal_test(formula = AIM_CD4.CXCR5.PD1..CD38hi.CD69.CD200._FreqParent ~ timeCategory, data=subset(subsetData, Prior.COVID.infection. == 'No')) 
dunn_test(AIM_CD4.CXCR5.PD1..CD38hi.CD69.CD200._FreqParent ~ timeCategory, data=subset(subsetData, Prior.COVID.infection. == 'No'))
bartlett.test(AIM_CD4.CXCR5.PD1..CD38hi.CD69.CD200._FreqParent ~ timeCategory, data=subset(subsetData, Prior.COVID.infection. == 'Yes'))
kruskal_test(formula = AIM_CD4.CXCR5.PD1..CD38hi.CD69.CD200._FreqParent ~ timeCategory, data=subset(subsetData, Prior.COVID.infection. == 'Yes')) 
tukey_hsd( aov(AIM_CD4.CXCR5.PD1..CD38hi.CD69.CD200._FreqParent ~ timeCategory, data=subset(subsetData, Prior.COVID.infection. == 'Yes')) )

subsetData <- subset(mergedData, timeCategory == "Baseline" | timeCategory == "Post 2nd dose")
FC_response2 <- dcast( subsetData, `Record.ID`+`Alias` + `Prior.COVID.infection.` ~`timeCategory`, value.var = c("AIM_CD4.CXCR5.PD1..CD38hi.CD69.CD200._FreqParent")) 
FC_response2$FoldChange <- FC_response2$`Post 2nd dose`/FC_response2$`Baseline`; FC_response2$Cohort <- NULL
FC_response2 <- FC_response2[!is.infinite(FC_response2$FoldChange), ]
FC_response2 %>% group_by(Prior.COVID.infection.) %>% get_summary_stats(type = "common") 


#' ------------------ B cell analyses --------------------------
#'

subsetData <- subset(mergedData,  timeCategory != "two Weeks" & timeCategory != "2 wks post 2nd dose");   subsetData$timeCategory <- factor(subsetData$timeCategory, levels = c("Baseline", "Post 1st dose", "Pre 2nd dose", "Post 2nd dose","One month post\n2nd dose"))
prePostTime(subsetData, xData = "DPV", yData="CD19_.CD27..CD38._FreqParent", fillParam = "Prior.COVID.infection.", groupby="Record.ID", title = "PB", 
            xLabel = "Prior COVID?", yLabel = "CD27+CD38+ (% CD19)", repMeasures = F) +  geom_vline(xintercept = 0,linetype="dashed", alpha=0.5)+  coord_cartesian(xlim = c(-1,12)) 
subsetData$DPV <- subsetData$DPV - difftime(subsetData$Vaccine.2.date, subsetData$Vaccine.1.date, units="days" ) 
prePostTime(subsetData, xData = "DPV", yData="CD19_.CD27..CD38._FreqParent", fillParam = "Prior.COVID.infection.", groupby="Record.ID", title = "PB", 
            xLabel = "Prior COVID?", yLabel = "CD27+CD38+ (% CD19)", repMeasures = F) +  geom_vline(xintercept = 0,linetype="dashed", alpha=0.5) +  coord_cartesian(xlim = c(-5,12)) 

subsetData <- subset(mergedData,  timeCategory != "two Weeks" & timeCategory != "2 wks post 2nd dose");
subsetData <- subsetData[which(subsetData$Tube == "HEP"),]
subsetData$timeCategory <- factor(subsetData$timeCategory, levels = c("Baseline", "Post 1st dose", "Pre 2nd dose", "Post 2nd dose","One month post\n2nd dose"))
prePostTime(subsetData, xData = "timeCategory", yData="CD19_.CD27..CD38._FreqParent", fillParam = "Prior.COVID.infection.", groupby="Record.ID", title = "Plasmablasts", 
            xLabel = " ", yLabel = "CD27+CD38+ (% CD19)", repMeasures = F, newform = T) #  + ggrepel::geom_text_repel(aes(label=Record.ID),size=2)
# ggsave(filename = "./Images/PB_responses_bothCohorts_overTime.pdf")

bartlett.test(CD19_.CD27..CD38._FreqParent ~ timeCategory, data=subset(subsetData, Prior.COVID.infection. == 'No'))
kruskal_test(formula = CD19_.CD27..CD38._FreqParent ~ timeCategory, data=subset(subsetData, Prior.COVID.infection. == 'No')) 
dunn_test(CD19_.CD27..CD38._FreqParent ~ timeCategory, data=subset(subsetData, Prior.COVID.infection. == 'No'))
bartlett.test(CD19_.CD27..CD38._FreqParent ~ timeCategory, data=subset(subsetData, Prior.COVID.infection. == 'Yes'))
kruskal_test(formula = CD19_.CD27..CD38._FreqParent ~ timeCategory, data=subset(subsetData, Prior.COVID.infection. == 'Yes')) 
dunn_test(CD19_.CD27..CD38._FreqParent ~ timeCategory, data=subset(subsetData, Prior.COVID.infection. == 'Yes'))
# 
# subsetData = subset(mergedData, mergedData$shortForm == 'bL' | mergedData$shortForm=='oW'); subsetData = subset(subsetData, subsetData$Prior.COVID.infection.=="No")
# prePostTime(data=subsetData, xData = "shortForm", yData="CD19_.CD27..CD38._FreqParent", fillParam = "Prior.COVID.infection.", groupby="Record.ID", title = "PB - Vax1", 
#             xLabel = " ", yLabel = "CD27+CD38+ (% CD19)") + scale_y_continuous(breaks = seq(0,21,0.5), limits = c(0,2))  + 
#   scale_x_discrete(labels= c("Baseline","Post\n1st dose"))# + ggrepel::geom_text_repel(aes(label = Record.ID))
# # ggsave(filename = "./Images/PB_vax1_healthy_paired.pdf", width=4)
# 
# 
# subsetData = subset(mergedData, mergedData$shortForm == '3W' | mergedData$shortForm=='4W'); subsetData = subset(subsetData, subsetData$Prior.COVID.infection.=="No")
# prePostTime(data=subsetData, xData = "shortForm", yData="CD19_.CD27..CD38._FreqParent", fillParam = "Prior.COVID.infection.", groupby="Record.ID", 
#             title = "PB - Vax2", xLabel = " ", yLabel = "CD27+CD38+ (% CD19)") + scale_y_continuous(breaks = seq(0,21,0.5), limits = c(0,2))  + 
#   scale_x_discrete(labels= c("Pre\n2nd dose","Post\n2nd dose"))# + ggrepel::geom_text_repel(aes(label = Record.ID))
# # ggsave(filename = "./Images/PB_vax2_healthy_paired.pdf", width=4)
# 
# subsetData = subset(mergedData, mergedData$shortForm == 'bL' | mergedData$shortForm=='oW'); subsetData = subset(subsetData, subsetData$Prior.COVID.infection.=="Yes")
# prePostTime(data=subsetData, xData = "shortForm", yData="CD19_.CD27..CD38._FreqParent", fillParam = "Prior.COVID.infection.", groupby="Record.ID", 
#             title = "PB - Vax1", xLabel = " ", yLabel = "CD27+CD38+ (% CD19)") +    
#   scale_y_continuous(breaks = seq(0,21,0.5), limits = c(0,3))  + scale_x_discrete(labels= c("Baseline","Post\n1st dose"))+ scale_fill_manual(values=c("#B5B2F1")) 
# # ggsave(filename = "./Images/PB_vax1_covid_paired.pdf", width=4)
# 
# subsetData = subset(mergedData, mergedData$shortForm == '3W' | mergedData$shortForm=='4W'); subsetData = subset(subsetData, subsetData$Prior.COVID.infection.=="Yes")
# prePostTime(data=subsetData, xData = "shortForm", yData="CD19_.CD27..CD38._FreqParent", fillParam = "Prior.COVID.infection.", groupby="Record.ID", 
#             title = "PB - Vax2", xLabel = " ", yLabel = "CD27+CD38+ (% CD19)") +    
#   scale_y_continuous(breaks = seq(0,21,0.5), limits = c(0,3))  + scale_x_discrete(labels= c("Pre\n2nd dose","Post\n2nd dose"))+ scale_fill_manual(values=c("#B5B2F1")) 
# # ggsave(filename = "./Images/PB_vax2_covid_paired.pdf", width=4)
# 
# subsetData = subset(mergedData, mergedData$shortForm == 'bL' ); subsetData = subset(subsetData, subsetData$Prior.COVID.infection.=="No")
# univScatter(data = subsetData, xData = "FCtfh_Vax1", yData = "FCPB_Vax1", fillParam = "Prior.COVID.infection.", title = " ", xLabel = "Fold-change Tfh", 
#             yLabel = "Fold-change PB", position = "left")
# subsetData = subset(mergedData, mergedData$shortForm == 'bL' ); subsetData = subset(subsetData, subsetData$Prior.COVID.infection.=="Yes")
# univScatter(data = subsetData, xData = "FCtfh_Vax1", yData = "FCPB_Vax1", fillParam = "Prior.COVID.infection.", title = " ", xLabel = "Fold-change Tfh", 
#             yLabel = "Fold-change PB", position = "left")+ scale_fill_manual(values=c("#B5B2F1")) 
# 
# 
# subsetData = subset(mergedData, mergedData$shortForm == 'oW')
# subsetData %>% group_by(Prior.COVID.infection.) %>% shapiro_test(CD19_.CD27..CD38._FreqParent)
# twoSampleBar(data=subsetData, xData = "Prior.COVID.infection.", yData="CD19_.CD27..CD38._FreqParent", fillParam = "Prior.COVID.infection.", title = "PB post 1st dose", 
#              yLabel = "CD27+CD38+ (% CD19)", nonparam = T) + scale_y_continuous(breaks = seq(0,20,1), limits = c(0,2)) 
# 
# subsetData = subset(mergedData, mergedData$shortForm == '4W')
# subsetData %>% group_by(Prior.COVID.infection.) %>% shapiro_test(CD19_.CD27..CD38._FreqParent)
# twoSampleBar(data=subsetData, xData = "Prior.COVID.infection.", yData="CD19_.CD27..CD38._FreqParent", fillParam = "Prior.COVID.infection.", title = "PB post 2nd dose", 
#              yLabel = "CD27+CD38+ (% CD19)", nonparam = T) + scale_y_continuous(breaks = seq(0,20,1), limits = c(0,2)) 


subsetData <- subset(mergedData, timeCategory != "two Weeks" & timeCategory != "2 wks post 2nd dose");
subsetData <- subsetData[which(subsetData$Tube == "HEP"),] 
subsetData$timeCategory <- factor(subsetData$timeCategory, levels = c("Baseline", "Post 1st dose", "Pre 2nd dose", "Post 2nd dose", "One month post\n2nd dose"))
subsetData$CD19_.CD27..CD38..CD138._FreqParent <- subsetData$CD19_.CD27..CD38..CD138._FreqParent * subsetData$CD19_.CD27..CD38._FreqParent / 100
if( is.na(sum(subsetData$CD19_.CD27..CD38..CD138._FreqParent, na.rm = F)))    # look for NA values, should give NA result if true
{ subsetData <- subsetData[-which(is.na(subsetData$CD19_.CD27..CD38..CD138._FreqParent)), ]  }
prePostTime(data=subsetData, xData = "timeCategory", yData="CD19_.CD27..CD38..CD138._FreqParent", fillParam = "Prior.COVID.infection.", groupby="Record.ID", 
            title = "CD138+ Plasmablasts", xLabel = " ", yLabel = "CD27+CD38+CD138+CD20- (% CD19)", repMeasures = F, newform = T) + 
  scale_y_continuous(breaks = seq(0,2,0.1), limits = c(0,0.6))
# ggsave(filename = "./Images/PB_deepPhenotype_bothCohorts_overTime.pdf")
bartlett.test(CD19_.CD27..CD38..CD138._FreqParent ~ timeCategory, data=subset(subsetData, Prior.COVID.infection. == 'No'))
kruskal_test(formula = CD19_.CD27..CD38..CD138._FreqParent ~ timeCategory, data=subset(subsetData, Prior.COVID.infection. == 'No')) 
dunn_test(CD19_.CD27..CD38..CD138._FreqParent ~ timeCategory, data=subset(subsetData, Prior.COVID.infection. == 'No'))
bartlett.test(CD19_.CD27..CD38..CD138._FreqParent ~ timeCategory, data=subset(subsetData, Prior.COVID.infection. == 'Yes'))
kruskal_test(formula = CD19_.CD27..CD38..CD138._FreqParent ~ timeCategory, data=subset(subsetData, Prior.COVID.infection. == 'Yes')) 
dunn_test(CD19_.CD27..CD38..CD138._FreqParent ~ timeCategory, data=subset(subsetData, Prior.COVID.infection. == 'Yes'))




#' ------------------ Tbet+ CD11c+ B cell analyses --------------------------

subsetData <- subset(mergedData,  timeCategory != "two Weeks" & timeCategory != "2 wks post 2nd dose");  
subsetData$timeCategory <- factor(subsetData$timeCategory, levels = c("Baseline", "Post 1st dose", "Pre 2nd dose", "Post 2nd dose","One month post\n2nd dose"))
prePostTime(subsetData, xData = "timeCategory", yData="CD19_.Nonnaive.B.Tbet..CD11C._FreqParent", fillParam = "Prior.COVID.infection.", groupby="Record.ID", title = "Tbet+CD11c+ B cells", 
            xLabel = " ", yLabel = "Tbet+CD11+ (% nonnavB)", repMeasures = F, newform =T) 
# ggsave(filename = "./Images/CD11cTbet_Bcells_bothCohorts_overTime.pdf")


#' ------------------ GzmB CD 8 analyses --------------------------
subsetData <- subset(mergedData,  timeCategory != "two Weeks"& timeCategory != "2 wks post 2nd dose"); 
subsetData$timeCategory <- factor(subsetData$timeCategory, levels = c("Baseline", "Post 1st dose", "Pre 2nd dose", "Post 2nd dose", "One month post\n2nd dose"))

prePostTime(data = subsetData, xData = "timeCategory", yData="CD8_.CD38.Ki67..GzmB..CD8_FreqParent", fillParam = "Prior.COVID.infection.", groupby="Record.ID", title = "Activated CD8", 
            xLabel = " ", yLabel = "GzmB+ (% CD8+Ki67+CD38+)", repMeasures=F, newform=T) +
  scale_y_continuous(limits = c(0,110), breaks = seq(0,100,10)) # + 
# ggformula::geom_spline(aes_string(x="DPV", y="CD8_.CD38.Ki67..GzmB..CD8_FreqParent"), color="blue",size=1, spar=0.5) 
# ggrepel::geom_text_repel(data = subset(subsetData, Alias == "HV-002"), aes(label=Record.ID),size=2)
# ggsave(filename = "./Images/CD8_CD38hiKi67hi_GzmB.pdf" )

fit <- aov(CD8_.CD38.Ki67..GzmB..CD8_FreqParent ~ timeCategory, data=subset(subsetData, Prior.COVID.infection. == 'No')); tukey_hsd(fit); formatC(tukey_hsd(fit)$p.adj, format="e")
fit <- aov(CD8_.CD38.Ki67..GzmB..CD8_FreqParent ~ timeCategory, data=subset(subsetData, Prior.COVID.infection. == 'Yes')); tukey_hsd(fit); formatC(tukey_hsd(fit)$p.adj, format="e")


prePostTime(data = subsetData, xData = "timeCategory", yData="CD4_.CD38.Ki67..GzmB..CD4_FreqParent", fillParam = "Prior.COVID.infection.", groupby="Record.ID", title = "Activated CD4", 
            xLabel = "Prior COVID?", yLabel = "GzmB+ (% CD4+Ki67+CD38+)", repMeasures=F) #+ 
# ggformula::geom_spline(aes_string(x="DPV", y="CD4_.CD38.Ki67..GzmB..CD4_FreqParent"), color="blue",size=1, spar=0.5) #+ ggrepel::geom_text_repel(aes(label=Record.ID),size=2)

fit <- aov(CD4_.CD38.Ki67..GzmB..CD4_FreqParent ~ timeCategory, data=subset(subsetData, Prior.COVID.infection. == 'No')); tukey_hsd(fit); formatC(tukey_hsd(fit)$p.adj, format="e")
fit <- aov(CD4_.CD38.Ki67..GzmB..CD4_FreqParent ~ timeCategory, data=subset(subsetData, Prior.COVID.infection. == 'Yes')); tukey_hsd(fit); formatC(tukey_hsd(fit)$p.adj, format="e")

# 
# subsetData = subset(mergedData, mergedData$shortForm == 'oW')
# subsetData %>% group_by(Prior.COVID.infection.) %>% shapiro_test(CD8_.CD38.Ki67..GzmB..CD8_FreqParent)
# twoSampleBar(data=subsetData, xData = "Prior.COVID.infection.", yData="CD8_.CD38.Ki67..GzmB..CD8_FreqParent", fillParam = "Prior.COVID.infection.", title = "GzmB+ Activated CD8 post 1st dose", 
#              yLabel = "GzmB+ (% Activated CD8)", nonparam = T) + scale_y_continuous(breaks = seq(0,100,10), limits = c(0,100)) 
# 
# subsetData = subset(mergedData, mergedData$shortForm == '4W')
# subsetData %>% group_by(Prior.COVID.infection.) %>% shapiro_test(CD8_.CD38.Ki67..GzmB..CD8_FreqParent)
# twoSampleBar(data=subsetData, xData = "Prior.COVID.infection.", yData="CD8_.CD38.Ki67..GzmB..CD8_FreqParent", fillParam = "Prior.COVID.infection.", title = "GzmB+ Activated CD8 post 2nd dose", 
#              yLabel = "GzmB+ (% Activated CD8)", nonparam = T) + scale_y_continuous(breaks = seq(0,100,10), limits = c(0,100)) 
# 
# subsetData = subset(mergedData, mergedData$shortForm == 'oM')
# subsetData %>% group_by(Prior.COVID.infection.) %>% shapiro_test(CD8_.CD38.Ki67..GzmB..CD8_FreqParent)
# twoSampleBar(data=subsetData, xData = "Prior.COVID.infection.", yData="CD8_.CD38.Ki67..GzmB..CD8_FreqParent", fillParam = "Prior.COVID.infection.", title = "GzmB+ Activated CD8 at one month", 
#              yLabel = "GzmB+ (% Activated CD8)", nonparam = T) + scale_y_continuous(breaks = seq(0,100,10), limits = c(0,100)) 



#' ------------------ GzmB+Tbet+ CD8 analyses --------------------------
subsetData <- subset(mergedData,  timeCategory != "two Weeks");  
subsetData$timeCategory <- factor(subsetData$timeCategory, levels = c("Baseline", "Post 1st dose", "Pre 2nd dose", "Post 2nd dose", "One month post\n2nd dose"))

prePostTime(data = subsetData, xData = "timeCategory", yData="CD8_.CD38.Ki67..GzmB..CD8_FreqParent", fillParam = "Prior.COVID.infection.", groupby="Record.ID", title = "Activated CD8", 
            xLabel = "Prior COVID?", yLabel = "GzmB+ (% CD8+Ki67+CD38+)", repMeasures=F) #+ 
# ggformula::geom_spline(aes_string(x="DPV", y="CD8_.CD38.Ki67..GzmB..CD8_FreqParent"), color="blue",size=1, spar=0.5) #+ ggrepel::geom_text_repel(aes(label=Record.ID),size=2)

fit <- aov(CD8_.CD38.Ki67..GzmB..CD8_FreqParent ~ timeCategory, data=subset(subsetData, Prior.COVID.infection. == 'No')); tukey_hsd(fit); formatC(tukey_hsd(fit)$p.adj, format="e")
fit <- aov(CD8_.CD38.Ki67..GzmB..CD8_FreqParent ~ timeCategory, data=subset(subsetData, Prior.COVID.infection. == 'Yes')); tukey_hsd(fit); formatC(tukey_hsd(fit)$p.adj, format="e")



#' ------------------ Foxp3+ CD4 analyses --------------------------
subsetData <- subset(mergedData,  timeCategory != "two Weeks"  & timeCategory != "2 wks post 2nd dose");  
subsetData$timeCategory <- factor(subsetData$timeCategory, levels = c("Baseline", "Post 1st dose", "Pre 2nd dose", "Post 2nd dose", "One month post\n2nd dose"))
subsetData <- subsetData[-which(subsetData$Record.ID == "CV-026"),]     # exclude because multiple samples with not-believable Foxp3 stains
subsetData <- subsetData[-which(subsetData$Record.ID == "CV-033"),]     # exclude because missing file for Pre 2nd dose 

prePostTime(data = subsetData, xData = "timeCategory", yData="CD4_.CD38.Ki67..Foxp3._FreqParent", fillParam = "Prior.COVID.infection.", groupby="Record.ID", title = "Activated CD4", 
            xLabel = "Prior COVID?", yLabel = "Foxp3+ (% CD4+Ki67+CD38+)", repMeasures=F) #+ 
# ggformula::geom_spline(aes_string(x="timeCategory", y="CD4_.CD38.Ki67..Foxp3._FreqParent"), color="blue",size=1, spar=0.5) #+ 
# ggrepel::geom_text_repel(aes(label=Record.ID),size=2)

fit <- aov(CD4_.CD38.Ki67..Foxp3._FreqParent ~ timeCategory, data=subset(subsetData, Prior.COVID.infection. == 'No')); tukey_hsd(fit); formatC(tukey_hsd(fit)$p.adj, format="e")
fit <- aov(CD4_.CD38.Ki67..Foxp3._FreqParent ~ timeCategory, data=subset(subsetData, Prior.COVID.infection. == 'Yes')); tukey_hsd(fit); formatC(tukey_hsd(fit)$p.adj, format="e")



#' ------------------ CD71+ IgD- B cell analyses --------------------------
#' 
subsetData <- subset(mergedData,  timeCategory != "two Weeks" & timeCategory != "2 wks post 2nd dose" );  
subsetData$timeCategory <- factor(subsetData$timeCategory, levels = c("Baseline", "Post 1st dose", "Pre 2nd dose", "Post 2nd dose", "One month post\n2nd dose"))
subsetData <- subsetData[which(subsetData$Tube == "HEP"),]
prePostTime(data = subsetData, xData = "timeCategory", yData="CD19_.IgD..CD71._FreqParent", fillParam = "Prior.COVID.infection.", groupby="Record.ID", 
            title = "CD71+ B cells", xLabel = "", yLabel = "IgD-CD71+ (% CD19+)", repMeasures=F, newform = T) #ggrepel::geom_text_repel(aes(label = Alias))
# ggsave(filename = "./Images/CD71hi_Bcells_bothCohorts_overTime.pdf")
bartlett.test(CD19_.IgD..CD71._FreqParent ~ timeCategory, data=subset(subsetData, Prior.COVID.infection. == 'No'))
kruskal_test(formula = CD19_.IgD..CD71._FreqParent ~ timeCategory, data=subset(subsetData, Prior.COVID.infection. == 'No')) 
# dunn_test(CD19_.IgD..CD71._FreqParent ~ timeCategory, data=subset(subsetData, Prior.COVID.infection. == 'No'))
bartlett.test(CD19_.IgD..CD71._FreqParent ~ timeCategory, data=subset(subsetData, Prior.COVID.infection. == 'Yes'))
tukey_hsd(aov(CD19_.IgD..CD71._FreqParent ~ timeCategory, data=subset(subsetData, Prior.COVID.infection. == 'Yes')) ) 
# kruskal_test(formula = CD19_.IgD..CD71._FreqParent ~ timeCategory, data=subset(subsetData, Prior.COVID.infection. == 'Yes')) 
# dunn_test(CD19_.IgD..CD71._FreqParent ~ timeCategory, data=subset(subsetData, Prior.COVID.infection. == 'Yes'))

prePostTime(data = subsetData, xData = "timeCategory", yData="CD19_.IgD..CD71..CD20.CD38._FreqParent", fillParam = "Prior.COVID.infection.", groupby="Record.ID", 
            title = "ASC - Ellebedy subsets", xLabel = "Prior COVID?", yLabel = "CD20loCD38hi (% CD19+IgD-CD71+)", repMeasures=F) #ggrepel::geom_text_repel(aes(label = Alias))

prePostTime(data = subsetData, xData = "timeCategory", yData="CD19_.IgD..CD71..ActivBCells_FreqParent", fillParam = "Prior.COVID.infection.", groupby="Record.ID", 
            title = "ABC - Ellebedy subsets", xLabel = "Prior COVID?", yLabel = "CD20hiCD38lo (% CD19+IgD-CD71+)", repMeasures=F)# +ggrepel::geom_text_repel(aes(label = Alias))



#' ------------------ CD21lo B cell analyses --------------------------
#' 
subsetData <- subset(mergedData,  timeCategory != "two Weeks"  & timeCategory != "2 wks post 2nd dose");  
subsetData <- subsetData[which(subsetData$Tube == "HEP"),]
subsetData <- subset(subsetData, Label != "PHI-398_V2" & Label != "PHI-021_V1" & Label != "HV-079_V1" )   # exclude because failed QC for CD21lo:  CV-028_PHI-398_oW and CV-030_HV-078_4W and CV-022_PHI-021_bL

subsetData$timeCategory <- factor(subsetData$timeCategory, levels = c("Baseline", "Post 1st dose", "Pre 2nd dose", "Post 2nd dose", "One month post\n2nd dose"))
prePostTime(data = subsetData, xData = "timeCategory", yData="CD19_.CD21lo_FreqParent", fillParam = "Prior.COVID.infection.", groupby="Record.ID", 
            title = "CD21lo B cells", xLabel = " ", yLabel = "CD21lo (% CD19)", repMeasures=F, newform = T) + 
  scale_y_continuous(breaks = seq(0,20,1))#ggrepel::geom_text_repel(aes(label = Alias))
# ggsave(filename = "./Images/CD21lo_bothCohorts_overTime.pdf")
# fit <- aov(CD19_.CD21lo_FreqParent ~ timeCategory, data=subset(subsetData, Prior.COVID.infection. == 'No')); tukey_hsd(fit); formatC(tukey_hsd(fit)$p.adj, format="e")
bartlett.test(CD19_.CD21lo_FreqParent ~ timeCategory, data=subset(subsetData, Prior.COVID.infection. == 'No'))
kruskal_test(formula = CD19_.CD21lo_FreqParent ~ timeCategory, data=subset(subsetData, Prior.COVID.infection. == 'No')) 
dunn_test(CD19_.CD21lo_FreqParent ~ timeCategory, data=subset(subsetData, Prior.COVID.infection. == 'No'))
bartlett.test(CD19_.CD21lo_FreqParent ~ timeCategory, data=subset(subsetData, Prior.COVID.infection. == 'Yes'))
tukey_hsd(aov(CD19_.CD21lo_FreqParent ~ timeCategory, data=subset(subsetData, Prior.COVID.infection. == 'Yes')) ) 
# kruskal_test(formula = CD19_.CD21lo_FreqParent ~ timeCategory, data=subset(subsetData, Prior.COVID.infection. == 'Yes')) 
# dunn_test(CD19_.CD21lo_FreqParent ~ timeCategory, data=subset(subsetData, Prior.COVID.infection. == 'Yes'))





#' ------------------ ELISpot analyses --------------------------
#'
subsetData <- subset(mergedData, timeCategory == 'Post 1st dose')
a.1 <- twoSampleBar(data = subsetData, xData = 'Prior.COVID.infection.',yData='Elispot_IgG_S1', fillParam = 'Prior.COVID.infection.',title = "S1", 
             yLabel = "ASC per 1e6 PBMC", nonparam=T) + 
  scale_y_continuous(trans='pseudo_log', limits = c(0,10000), breaks=c(10^(0:4)), labels=trans_format('log10',math_format(10^.x)), minor_breaks =5*10^(0:10))
a.2 <- twoSampleBar(data = subsetData, xData = 'Prior.COVID.infection.',yData='Elispot_IgG_S2', fillParam = 'Prior.COVID.infection.',title = "S2", 
             yLabel = " ", nonparam=T)+ 
  scale_y_continuous(trans='pseudo_log', limits = c(0,10000), breaks=c(10^(0:4)), labels=trans_format('log10',math_format(10^.x)), minor_breaks =5*10^(0:10))
a.3 <- twoSampleBar(data = subsetData, xData = 'Prior.COVID.infection.',yData='Elispot_IgG_RBD', fillParam = 'Prior.COVID.infection.',title = "RBD", 
             yLabel = " ", nonparam=T)+ 
  scale_y_continuous(trans='pseudo_log', limits = c(0,10000), breaks=c(10^(0:4)), labels=trans_format('log10',math_format(10^.x)), minor_breaks =5*10^(0:10))
grid.arrange(a.1,a.2,a.3, nrow=1)
# ggpubr::ggexport(plotlist = list(a.1, a.2, a.3), filename = "./Images/Elispots_IgG_post1stDose_gridarrange.pdf", nrow = 1, width = 12)

a.1 <- twoSampleBar(data = subsetData, xData = 'Prior.COVID.infection.',yData='Elispot_IgA_S1', fillParam = 'Prior.COVID.infection.',title = "S1", 
             yLabel = "ASC per 1e6 PBMC", nonparam=T) + 
  scale_y_continuous(trans='pseudo_log', limits = c(0,10000), breaks=c(10^(0:4)), labels=trans_format('log10',math_format(10^.x)), minor_breaks =5*10^(0:10))
a.2 <- twoSampleBar(data = subsetData, xData = 'Prior.COVID.infection.',yData='Elispot_IgA_S2', fillParam = 'Prior.COVID.infection.',title = "S2", 
             yLabel = " ", nonparam=T)+ 
  scale_y_continuous(trans='pseudo_log', limits = c(0,10000), breaks=c(10^(0:4)), labels=trans_format('log10',math_format(10^.x)), minor_breaks =5*10^(0:10))
a.3 <- twoSampleBar(data = subsetData, xData = 'Prior.COVID.infection.',yData='Elispot_IgA_RBD', fillParam = 'Prior.COVID.infection.',title = "RBD", 
             yLabel = " ", nonparam=T)+ 
  scale_y_continuous(trans='pseudo_log', limits = c(0,10000), breaks=c(10^(0:4)), labels=trans_format('log10',math_format(10^.x)), minor_breaks =5*10^(0:10))
grid.arrange(a.1,a.2,a.3, nrow=1)
# ggpubr::ggexport(plotlist = list(a.1, a.2, a.3), filename = "./Images/Elispots_IgA_post1stDose_gridarrange.pdf", nrow = 1, width = 12 )

a.1 <- twoSampleBar(data = subsetData, xData = 'Prior.COVID.infection.',yData='Elispot_IgM_S1', fillParam = 'Prior.COVID.infection.',title = "S1", 
             yLabel = "ASC per 1e6 PBMC", nonparam=T)+ 
  scale_y_continuous(trans='pseudo_log', limits = c(0,10000), breaks=c(10^(0:4)), labels=trans_format('log10',math_format(10^.x)), minor_breaks =5*10^(0:10))
a.2 <- twoSampleBar(data = subsetData, xData = 'Prior.COVID.infection.',yData='Elispot_IgM_S2', fillParam = 'Prior.COVID.infection.',title = "S2", 
             yLabel = " ", nonparam=T)+ 
  scale_y_continuous(trans='pseudo_log', limits = c(0,10000), breaks=c(10^(0:4)), labels=trans_format('log10',math_format(10^.x)), minor_breaks =5*10^(0:10))
a.3 <- twoSampleBar(data = subsetData, xData = 'Prior.COVID.infection.',yData='Elispot_IgM_RBD', fillParam = 'Prior.COVID.infection.',title = "RBD", 
             yLabel = " ", nonparam=T)+ 
  scale_y_continuous(trans='pseudo_log', limits = c(0,10000), breaks=c(10^(0:4)), labels=trans_format('log10',math_format(10^.x)), minor_breaks =5*10^(0:10))
grid.arrange(a.1,a.2,a.3, nrow=1)
# ggpubr::ggexport(plotlist = list(a.1, a.2, a.3), filename = "./Images/Elispots_IgM_post1stDose_gridarrange.pdf", nrow = 1, width = 12 )

subsetData <- subset(mergedData, timeCategory == 'Post 2nd dose')
a.1 <- twoSampleBar(data = subsetData, xData = 'Prior.COVID.infection.',yData='Elispot_IgG_S1', fillParam = 'Prior.COVID.infection.',title = "S1", 
                    yLabel = "ASC per 1e6 PBMC", nonparam=T) + 
  scale_y_continuous(trans='pseudo_log', limits = c(0,10000), breaks=c(10^(0:4)), labels=trans_format('log10',math_format(10^.x)), minor_breaks =5*10^(0:10))
a.2 <- twoSampleBar(data = subsetData, xData = 'Prior.COVID.infection.',yData='Elispot_IgG_S2', fillParam = 'Prior.COVID.infection.',title = "S2", 
                    yLabel = " ", nonparam=T)+ 
  scale_y_continuous(trans='pseudo_log', limits = c(0,10000), breaks=c(10^(0:4)), labels=trans_format('log10',math_format(10^.x)), minor_breaks =5*10^(0:10))
a.3 <- twoSampleBar(data = subsetData, xData = 'Prior.COVID.infection.',yData='Elispot_IgG_RBD', fillParam = 'Prior.COVID.infection.',title = "RBD", 
                    yLabel = " ", nonparam=T)+ 
  scale_y_continuous(trans='pseudo_log', limits = c(0,10000), breaks=c(10^(0:4)), labels=trans_format('log10',math_format(10^.x)), minor_breaks =5*10^(0:10))
grid.arrange(a.1,a.2,a.3, nrow=1)
# ggpubr::ggexport(plotlist = list(a.1, a.2, a.3), filename = "./Images/Elispots_IgG_post2ndDose_gridarrange.pdf", nrow = 1, width = 12)

a.1 <- twoSampleBar(data = subsetData, xData = 'Prior.COVID.infection.',yData='Elispot_IgA_S1', fillParam = 'Prior.COVID.infection.',title = "S1", 
                    yLabel = "ASC per 1e6 PBMC", nonparam=T) + 
  scale_y_continuous(trans='pseudo_log', limits = c(0,10000), breaks=c(10^(0:4)), labels=trans_format('log10',math_format(10^.x)), minor_breaks =5*10^(0:10))
a.2 <- twoSampleBar(data = subsetData, xData = 'Prior.COVID.infection.',yData='Elispot_IgA_S2', fillParam = 'Prior.COVID.infection.',title = "S2", 
                    yLabel = " ", nonparam=T)+ 
  scale_y_continuous(trans='pseudo_log', limits = c(0,10000), breaks=c(10^(0:4)), labels=trans_format('log10',math_format(10^.x)), minor_breaks =5*10^(0:10))
a.3 <- twoSampleBar(data = subsetData, xData = 'Prior.COVID.infection.',yData='Elispot_IgA_RBD', fillParam = 'Prior.COVID.infection.',title = "RBD", 
                    yLabel = " ", nonparam=T)+ 
  scale_y_continuous(trans='pseudo_log', limits = c(0,10000), breaks=c(10^(0:4)), labels=trans_format('log10',math_format(10^.x)), minor_breaks =5*10^(0:10))
grid.arrange(a.1,a.2,a.3, nrow=1)
# ggpubr::ggexport(plotlist = list(a.1, a.2, a.3), filename = "./Images/Elispots_IgA_post2ndDose_gridarrange.pdf", nrow = 1, width = 12 )

a.1 <- twoSampleBar(data = subsetData, xData = 'Prior.COVID.infection.',yData='Elispot_IgM_S1', fillParam = 'Prior.COVID.infection.',title = "S1", 
                    yLabel = "ASC per 1e6 PBMC", nonparam=T)+ 
  scale_y_continuous(trans='pseudo_log', limits = c(0,10000), breaks=c(10^(0:4)), labels=trans_format('log10',math_format(10^.x)), minor_breaks =5*10^(0:10))
a.2 <- twoSampleBar(data = subsetData, xData = 'Prior.COVID.infection.',yData='Elispot_IgM_S2', fillParam = 'Prior.COVID.infection.',title = "S2", 
                    yLabel = " ", nonparam=T)+ 
  scale_y_continuous(trans='pseudo_log', limits = c(0,10000), breaks=c(10^(0:4)), labels=trans_format('log10',math_format(10^.x)), minor_breaks =5*10^(0:10))
a.3 <- twoSampleBar(data = subsetData, xData = 'Prior.COVID.infection.',yData='Elispot_IgM_RBD', fillParam = 'Prior.COVID.infection.',title = "RBD", 
                    yLabel = " ", nonparam=T)+ 
  scale_y_continuous(trans='pseudo_log', limits = c(0,10000), breaks=c(10^(0:4)), labels=trans_format('log10',math_format(10^.x)), minor_breaks =5*10^(0:10))
grid.arrange(a.1,a.2,a.3, nrow=1)
# ggpubr::ggexport(plotlist = list(a.1, a.2, a.3), filename = "./Images/Elispots_IgM_post2ndDose_gridarrange.pdf", nrow = 1, width = 12 )

subsetData <- subset(mergedData, !is.na(Elispot_IgG_S1) & timeCategory != "Pre 2nd dose" )
# pdf(file = "./Images/Elispots_IgG_S1_prepostTime.pdf")
  prePostTime(data = subsetData, xData = "timeCategory", yData = "Elispot_IgG_S1", fillParam = "Prior.COVID.infection.", 
              groupby = "Alias", title = "S1", xLabel = " ", yLabel = "IgG ASC per 1e6 PBMC", exponential=T)+ 
    scale_y_continuous(trans='pseudo_log', limits = c(0,10000), breaks=c(10^(0:4)), labels=trans_format('log10',math_format(10^.x)), minor_breaks =5*10^(0:10)) 
# dev.off()

subsetData <- subset(mergedData, !is.na(Elispot_IgG_S2) & timeCategory != "Pre 2nd dose" )
# pdf(file = "./Images/Elispots_IgG_S2_prepostTime.pdf")
prePostTime(data = subsetData, xData = "timeCategory", yData = "Elispot_IgG_S2", fillParam = "Prior.COVID.infection.", 
                  groupby = "Alias", title = "S2", xLabel = " ", yLabel = "IgG ASC per 1e6 PBMC", exponential=T)+ 
  scale_y_continuous(trans='pseudo_log', limits = c(0,10000), breaks=c(10^(0:4)), labels=trans_format('log10',math_format(10^.x)), minor_breaks =5*10^(0:10))
# dev.off()

subsetData <- subset(mergedData, !is.na(Elispot_IgG_RBD) & timeCategory != "Pre 2nd dose" )
# pdf(file = "./Images/Elispots_IgG_RBD_prepostTime.pdf")
prePostTime(data = subsetData, xData = "timeCategory", yData = "Elispot_IgG_RBD", fillParam = "Prior.COVID.infection.", 
                  groupby = "Alias", title = "RBD", xLabel = " ", yLabel = "IgG ASC per 1e6 PBMC", exponential=T)+ 
  scale_y_continuous(trans='pseudo_log', limits = c(0,10000), breaks=c(10^(0:4)), labels=trans_format('log10',math_format(10^.x)), minor_breaks =5*10^(0:10))
# dev.off()



bivScatter(data1 = subset(mergedData, Prior.COVID.infection. == 'No' & timeCategory == "Post 2nd dose"), 
           data2 = subset(mergedData, Prior.COVID.infection. == "Yes" & timeCategory == "Post 2nd dose"),
           name1 = "Naive", name2 = "Experienced", xData = "Elispot_IgG_RBD", yData = 'Elispot_IgG_S1', fillParam = 'Prior.COVID.infection.',
           title = "IgG ELISpots post 2nd dose", xLabel = "RBD ASC per 10^6 PBMC", yLabel = "S1 ASC per 10^6 PBMC", statsOff = F) + 
  scale_y_continuous(trans='pseudo_log', limits = c(0,3e4), breaks=c(10^(0:4)), labels=trans_format('log10',math_format(10^.x)), minor_breaks =5*10^(0:10)) + 
  scale_x_continuous(trans='pseudo_log', limits = c(0,1e4), breaks=c(10^(0:4)), labels=trans_format('log10',math_format(10^.x)), minor_breaks =5*10^(0:10)) 
# ggsave(filename = "./Images/Elispots_IgG_S1-vs-RBD_correl_biv.pdf", width=8)

bivScatter(data1 = subset(mergedData, Prior.COVID.infection. == 'No' & timeCategory == "Post 2nd dose"), 
           data2 = subset(mergedData, Prior.COVID.infection. == "Yes" & timeCategory == "Post 2nd dose"),
           name1 = "Naive", name2 = "Experienced", xData = "Elispot_IgG_RBD", yData = 'Elispot_IgG_S2', fillParam = 'Prior.COVID.infection.',
           title = "IgG ELISpots post 2nd dose", xLabel = "RBD ASC per 10^6 PBMC", yLabel = "S2 ASC per 10^6 PBMC", statsOff = F) +
  scale_y_continuous(trans='pseudo_log', limits = c(0,1e3), breaks=c(10^(0:4)), labels=trans_format('log10',math_format(10^.x)), minor_breaks =5*10^(0:10)) + 
  scale_x_continuous(trans='pseudo_log', limits = c(0,1e3), breaks=c(10^(0:4)), labels=trans_format('log10',math_format(10^.x)), minor_breaks =5*10^(0:10)) 
# ggsave(filename = "./Images/Elispots_IgG_S2-vs-RBD_correl_biv.pdf", width = 8 )

#'' ----------------- DPV analysis post 1st vaccination ------------------------
subsetData <- subset(mergedData,  timeCategory != "two Weeks"& timeCategory != "2 wks post 2nd dose");
subsetData <- subsetData[which(subsetData$Tube == "HEP"),] 
# subsetData <- subset(subsetData, Record.ID != "CV-011" & Record.ID != "CV-012" & Record.ID != "CV-005")        # absence of Ki67 stain
subsetData$timeCategory <- factor(subsetData$timeCategory, levels = c("Baseline", "Post 1st dose", "Pre 2nd dose", "Post 2nd dose", "One month post\n2nd dose"))
prePostTime(subsetData, xData = "DPV", yData="Elispot_IgG_S1", fillParam = "Prior.COVID.infection.", groupby="Record.ID", title = "anti-S1 IgG ASC",
            xLabel = "Days", yLabel = "Spots per 10^6 PBMC", repMeasures = F, exponential=F, pathOff = T) + 
  scale_y_continuous(trans='pseudo_log', limits = c(0,1e3), breaks=c(10^(0:4)), labels=trans_format('log10',math_format(10^.x)), minor_breaks =5*10^(0:10)) + 
  coord_cartesian(xlim = c(-1,15)) + geom_vline(xintercept = 0,linetype="dashed", alpha=0.5) 
# ggsave(filename = "./Images/Elispots_IgG-S1_vs_DPV.pdf")
prePostTime(subsetData, xData = "DPV", yData="Elispot_IgG_S2", fillParam = "Prior.COVID.infection.", groupby="Record.ID", title = "anti-S2 IgG ASC",
            xLabel = "Days", yLabel = "Spots per 10^6 PBMC", repMeasures = F, exponential=F, pathOff = T) + 
  scale_y_continuous(trans='pseudo_log', limits = c(0,1e3), breaks=c(10^(0:4)), labels=trans_format('log10',math_format(10^.x)), minor_breaks =5*10^(0:10)) + 
  coord_cartesian(xlim = c(-1,15)) + geom_vline(xintercept = 0,linetype="dashed", alpha=0.5)
# ggsave(filename = "./Images/Elispots_IgG-S2_vs_DPV.pdf")
prePostTime(subsetData, xData = "DPV", yData="Elispot_IgG_RBD", fillParam = "Prior.COVID.infection.", groupby="Record.ID", title = "anti-RBD IgG ASC",
            xLabel = "Days", yLabel = "Spots per 10^6 PBMC", repMeasures = F, exponential=F, pathOff = T) + 
  scale_y_continuous(trans='pseudo_log', limits = c(0,1e3), breaks=c(10^(0:4)), labels=trans_format('log10',math_format(10^.x)), minor_breaks =5*10^(0:10)) + 
  coord_cartesian(xlim = c(-1,15)) + geom_vline(xintercept = 0,linetype="dashed", alpha=0.5)
# ggsave(filename = "./Images/Elispots_IgG-RBD_vs_DPV.pdf")

#'' ----------------- post 2nd vaccination ------------------------
subsetData$DPV <- subsetData$DPV - as.numeric(difftime(subsetData$Vaccine.2.date, subsetData$Vaccine.1.date, units="days" ) )
subsetData <- subset(subsetData, timeCategory == "Post 2nd dose")
prePostTime(subsetData, xData = "DPV", yData="Elispot_IgG_S1", fillParam = "Prior.COVID.infection.", groupby="Record.ID", title = "anti-S1 IgG ASC",
            xLabel = "Days", yLabel = "Spots per 10^6 PBMC", repMeasures = F, exponential=F, pathOff = T) + 
  scale_y_continuous(trans='pseudo_log', limits = c(0,1e3), breaks=c(10^(0:4)), labels=trans_format('log10',math_format(10^.x)), minor_breaks =5*10^(0:10)) + 
  coord_cartesian(xlim = c(-5,15)) + geom_vline(xintercept = 0,linetype="dashed", alpha=0.5) 
# ggsave(filename = "./Images/Elispots_IgG-S1_vs_DPV_vax2.pdf")
prePostTime(subsetData, xData = "DPV", yData="Elispot_IgG_S2", fillParam = "Prior.COVID.infection.", groupby="Record.ID", title = "anti-S2 IgG ASC",
            xLabel = "Days", yLabel = "Spots per 10^6 PBMC", repMeasures = F, exponential=F, pathOff = T) + 
  scale_y_continuous(trans='pseudo_log', limits = c(0,1e3), breaks=c(10^(0:4)), labels=trans_format('log10',math_format(10^.x)), minor_breaks =5*10^(0:10)) + 
  coord_cartesian(xlim = c(-5,15)) + geom_vline(xintercept = 0,linetype="dashed", alpha=0.5)
# ggsave(filename = "./Images/Elispots_IgG-S2_vs_DPV_vax2.pdf")
prePostTime(subsetData, xData = "DPV", yData="Elispot_IgG_RBD", fillParam = "Prior.COVID.infection.", groupby="Record.ID", title = "anti-RBD IgG ASC",
            xLabel = "Days", yLabel = "Spots per 10^6 PBMC", repMeasures = F, exponential=F, pathOff = T) + 
  scale_y_continuous(trans='pseudo_log', limits = c(0,1e3), breaks=c(10^(0:4)), labels=trans_format('log10',math_format(10^.x)), minor_breaks =5*10^(0:10)) + 
  coord_cartesian(xlim = c(-5,15)) + geom_vline(xintercept = 0,linetype="dashed", alpha=0.5)
# ggsave(filename = "./Images/Elispots_IgG-RBD_vs_DPV_vax2.pdf")






subsetData <- subset(mergedData, Prior.COVID.infection. == 'No' & timeCategory == "Post 1st dose")
subsetData <- subsetData[,grep("^Elispot", names(subsetData))]; subsetData <- subsetData[,-grep("Elispot_IgA_S1", names(subsetData))]
temp <- names(subsetData) ; temp <-  do.call(rbind.data.frame, strsplit(temp, split = "_"))
temp <- paste0(temp[,2]," anti-",temp[,3]); names(subsetData) <- temp
cor.elispots <- cor(subsetData, method="kendall" , use="pairwise.complete.obs" )
cor.elispots.pmat <- ggcorrplot::cor_pmat(subsetData, method="kendall", use="pairwise.complete.obs"  )
ggcorrplot::ggcorrplot(corr = cor.elispots, p.mat = cor.elispots.pmat, title = "ELISpots post 1st dose", legend.title = "Kendall tau", tl.cex = 18,pch.cex = 1, insig = "blank")
# ggsave(filename = "./Images/Elispots_ggcorrplot_Naive_post1st.pdf")
ggcorrplot::ggcorrplot(corr = cor.elispots, p.mat = cor.elispots.pmat, title = "ELISpots post 1st dose", legend.title = "Kendall tau", tl.cex = 18,pch.cex = 1) #, insig = "blank"
# ggsave(filename = "./Images/Elispots_ggcorrplot_Naive_post1st_full.pdf")

subsetData <- subset(mergedData, Prior.COVID.infection. == 'Yes' & timeCategory == "Post 1st dose")
subsetData <- subsetData[,grep("^Elispot", names(subsetData))]; 
temp <- names(subsetData) ; temp <-  do.call(rbind.data.frame, strsplit(temp, split = "_"))
temp <- paste0(temp[,2]," anti-",temp[,3]); names(subsetData) <- temp
cor.elispots <- cor(subsetData, method="kendall" , use="pairwise.complete.obs" )
cor.elispots.pmat <- ggcorrplot::cor_pmat(subsetData, method="kendall", use= "pairwise.complete.obs"  )
ggcorrplot::ggcorrplot(corr = cor.elispots, p.mat = cor.elispots.pmat, title = "ELISpots post 1st dose", legend.title = "Kendall tau", tl.cex = 18,pch.cex = 1, insig = "blank") #
# ggsave(filename = "./Images/Elispots_ggcorrplot_Experienced_post1st.pdf")
ggcorrplot::ggcorrplot(corr = cor.elispots, p.mat = cor.elispots.pmat, title = "ELISpots post 1st dose", legend.title = "Kendall tau", tl.cex = 18,pch.cex = 1) #, insig = "blank"
# ggsave(filename = "./Images/Elispots_ggcorrplot_Experienced_post1st_full.pdf")




subsetData <- subset(mergedData, Prior.COVID.infection. == 'No' & timeCategory == "Post 2nd dose")
subsetData <- subsetData[,grep("^Elispot", names(subsetData))]
temp <- names(subsetData) ; temp <-  do.call(rbind.data.frame, strsplit(temp, split = "_"))
temp <- paste0(temp[,2]," anti-",temp[,3]); names(subsetData) <- temp
cor.elispots <- cor(subsetData, method="kendall" , use="pairwise.complete.obs" )
cor.elispots.pmat <- ggcorrplot::cor_pmat(subsetData, method="kendall", use="pairwise.complete.obs"  )
ggcorrplot::ggcorrplot(corr = cor.elispots, p.mat = cor.elispots.pmat, title = "ELISpots post 2nd dose", legend.title = "Kendall tau", tl.cex = 18,pch.cex = 1, insig = "blank")
# ggsave(filename = "./Images/Elispots_ggcorrplot_Naive_post2nd.pdf")
ggcorrplot::ggcorrplot(corr = cor.elispots, p.mat = cor.elispots.pmat, title = "ELISpots post 2nd dose", legend.title = "Kendall tau", tl.cex = 18,pch.cex = 1) #, insig = "blank"
# ggsave(filename = "./Images/Elispots_ggcorrplot_Naive_post2nd_full.pdf")

subsetData <- subset(mergedData, Prior.COVID.infection. == 'Yes' & timeCategory == "Post 2nd dose")
subsetData <- subsetData[,grep("^Elispot", names(subsetData))]
temp <- names(subsetData) ; temp <-  do.call(rbind.data.frame, strsplit(temp, split = "_"))
temp <- paste0(temp[,2]," anti-",temp[,3]); names(subsetData) <- temp
cor.elispots <- cor(subsetData, method="kendall" , use="pairwise.complete.obs" )
cor.elispots.pmat <- ggcorrplot::cor_pmat(subsetData, method="kendall", use= "pairwise.complete.obs"  )
ggcorrplot::ggcorrplot(corr = cor.elispots, p.mat = cor.elispots.pmat, title = "ELISpots post 2nd dose", legend.title = "Kendall tau", tl.cex = 18,pch.cex = 1, insig = "blank") #
# ggsave(filename = "./Images/Elispots_ggcorrplot_Experienced_post2nd.pdf")
ggcorrplot::ggcorrplot(corr = cor.elispots, p.mat = cor.elispots.pmat, title = "ELISpots post 2nd dose", legend.title = "Kendall tau", tl.cex = 18,pch.cex = 1) #, insig = "blank"
# ggsave(filename = "./Images/Elispots_ggcorrplot_Experienced_post2nd_full.pdf")



#' ------------------ Antibody analyses --------------------------
#'
subsetData <- mergedData[-which(mergedData$timeCategory == "two Weeks"),] ; subsetData <- subsetData[-which(subsetData$timeCategory == "2 wks post 2nd dose"),]
subsetData$timeCategory <- factor(subsetData$timeCategory, levels = c("Baseline", "Post 1st dose", "Pre 2nd dose", "Post 2nd dose", "One month post\n2nd dose",
                                                                      "Four months post\n2nd dose"))
a <- linePlot(data = subsetData, xData = 'timeCategory', yData = 'binding_IgG_S1', groupby = 'Alias', xLabel = ' ', yLabel = "anti-S1 IgG titer", 
              title = "anti-S1 IgG titer", colorby = "Prior.COVID.infection.") + 
  geom_hline(yintercept = 25, linetype = "dashed",alpha=0.3) + annotate("text", x=5,y=15,label = "LOD", color="black", alpha=0.2)+
  scale_y_continuous(trans='pseudo_log', limits = c(0,1e8), breaks=c(10^(0:8)), labels=trans_format('log10',math_format(10^.x)), minor_breaks =5*10^(0:10)) # + 
# ggrepel::geom_text_repel(data = subset(mergedData, is.na(timeCategory) ), aes(label = Label))
a
# ggsave(filename = "./Images/BindingAb_S1_IgG_linePlot.pdf")
# plotly::ggplotly(a)

bartlett.test(binding_IgG_S1 ~ timeCategory, data=subset(subsetData, Prior.COVID.infection. == 'No'))
kruskal_test(formula = binding_IgG_S1 ~ timeCategory, data=subset(subsetData, Prior.COVID.infection. == 'No')) 
dunn_test(binding_IgG_S1 ~ timeCategory, data=subset(subsetData, Prior.COVID.infection. == 'No'))
kruskal_test(formula = binding_IgG_S1 ~ timeCategory, data=subset(subsetData, Prior.COVID.infection. == 'Yes')) 
bartlett.test(binding_IgG_S1 ~ timeCategory, data=subset(subsetData, Prior.COVID.infection. == 'Yes'))
dunn_test(binding_IgG_S1 ~ timeCategory, data=subset(subsetData, Prior.COVID.infection. == 'Yes'))

linePlot(data = subsetData, xData = 'timeCategory', yData = 'binding_IgG_S1', groupby = 'Alias', xLabel = ' ', yLabel = "anti-S1 IgG titer", 
         title = "anti-S1 IgG titer", colorby = "Prior.COVID.infection.", recentCOVID = T) + 
  geom_hline(yintercept = 25, linetype = "dashed",alpha=0.3) + annotate("text", x=5,y=15,label = "LOD", color="black", alpha=0.2)+
  scale_y_continuous(trans='pseudo_log', limits = c(0,1e7), breaks=c(10^(0:7)), labels=trans_format('log10',math_format(10^.x)), minor_breaks =5*10^(0:10))
# ggsave(filename = "./Images/BindingAb_S1_IgG_linePlot_recentCOVID.pdf")



twoSampleBar(data = subset(mergedData, timeCategory == "Post 2nd dose"), xData = "Prior.COVID.infection.", yData = "binding_IgG_S1", fillParam = "Prior.COVID.infection.", 
             title = "Post 2nd dose", yLabel = "anti-S1 IgG titer", nonparam = T) +   coord_cartesian(ylim=c(0e1,1e7),) + 
  scale_y_continuous(trans='pseudo_log', breaks=c(10^(0:7)), labels=trans_format('log10',math_format(10^.x)), minor_breaks =5*10^(0:10)) #+ ggrepel::geom_text_repel(aes(label = Label))
# ggsave(filename = "./Images/BindingAb_S1_IgG_post2ndDose.pdf", width=5)

twoSampleBar(data = subset(mergedData, timeCategory == "One month post\n2nd dose"), xData = "Prior.COVID.infection.", yData = "binding_IgG_S1", fillParam = "Prior.COVID.infection.", 
             title = "One month post", yLabel = "anti-S1 IgG titer", nonparam = T) +   coord_cartesian(ylim=c(0e1,1e7),) + 
  scale_y_continuous(trans='pseudo_log', breaks=c(10^(0:7)), labels=trans_format('log10',math_format(10^.x)), minor_breaks =5*10^(0:10)) #+ ggrepel::geom_text_repel(aes(label = Label))
# ggsave(filename = "./Images/BindingAb_S1_IgG_oneMonth.pdf", width=5)


subsetData <- mergedData[-which(mergedData$timeCategory == "two Weeks"),] 
subsetData$timeCategory <- factor(subsetData$timeCategory, levels = c("Baseline", "Post 1st dose", "Pre 2nd dose", "Post 2nd dose", "One month post\n2nd dose")) 
linePlot(data = subsetData, xData = 'timeCategory', yData = 'binding_IgA_S1', groupby = 'Alias', xLabel = ' ', yLabel = "anti-S1 IgA titer", 
         title = "anti-S1 IgA titer", colorby = "Prior.COVID.infection.") + 
  geom_hline(yintercept = 25, linetype = "dashed",alpha=0.3) + annotate("text", x=5,y=15,label = "LOD", color="black", alpha=0.2)+
  # scale_color_manual(name="Prior COVID?",values = c("#FFC26A","#B5B2F1")) + 
  scale_y_continuous(trans='pseudo_log', limits = c(0,1e6), breaks=c(10^(0:7)), labels=trans_format('log10',math_format(10^.x)), minor_breaks =5*10^(0:10)) #+ 
# ggrepel::geom_text_repel(data = subset(mergedData, Alias =="HV-002"), aes(label = Label))
# ggsave(filename = "./Images/BindingAb_S1_IgA_linePlot.pdf")

# bartlett.test(binding_IgA_S1 ~ timeCategory, data=subset(subsetData, Prior.COVID.infection. == 'No'))
kruskal_test(formula = binding_IgA_S1 ~ timeCategory, data=subset(subsetData, Prior.COVID.infection. == 'No')) 
dunn_test(binding_IgA_S1 ~ timeCategory, data=subset(subsetData, Prior.COVID.infection. == 'No'))
bartlett.test(binding_IgA_S1 ~ timeCategory, data=subset(subsetData, Prior.COVID.infection. == 'Yes'))
kruskal_test(formula = binding_IgA_S1 ~ timeCategory, data=subset(subsetData, Prior.COVID.infection. == 'Yes')) 
dunn_test(binding_IgA_S1 ~ timeCategory, data=subset(subsetData, Prior.COVID.infection. == 'Yes'))

subsetData <- mergedData[-which(mergedData$timeCategory == "two Weeks"),] ; 
# subsetData$timeCategory <- factor(subsetData$timeCategory, levels = c("Baseline", "Post 1st dose", "Pre 2nd dose", "Post 2nd dose"))
linePlot(data = subsetData, xData = 'timeCategory', yData = 'IC50_neutAb_log10', groupby = 'Alias', xLabel = ' ', yLabel = "log10 IC50", 
         title = "Neutralizing antibodies", colorby = "Prior.COVID.infection.") + 
  geom_hline(yintercept = 10, linetype = "dashed",alpha=0.3) + annotate("text", x=5,y=15,label = "LOD", color="black", alpha=0.2)+
  scale_y_continuous(trans='pseudo_log', limits = c(0,1e5), breaks=c(10^(0:6)), labels=trans_format('log10',math_format(10^.x)), minor_breaks =5*10^(0:10)) 
# ggsave(filename = "./Images/neutAb_linePlot.pdf")

# fit <- aov(IC50_neutAb_log10 ~ timeCategory, data=subset(subsetData, Prior.COVID.infection. == 'No')); tukey_hsd(fit)
# fit <- aov(IC50_neutAb_log10 ~ timeCategory, data=subset(subsetData, Prior.COVID.infection. == 'Yes')); tukey_hsd(fit)
bartlett.test(IC50_neutAb_log10 ~ timeCategory, data=subset(subsetData, Prior.COVID.infection. == 'No'))
kruskal_test(formula = IC50_neutAb_log10 ~ timeCategory, data=subset(subsetData, Prior.COVID.infection. == 'No')) 
dunn_test(IC50_neutAb_log10 ~ timeCategory, data=subset(subsetData, Prior.COVID.infection. == 'No'))
bartlett.test(IC50_neutAb_log10 ~ timeCategory, data=subset(subsetData, Prior.COVID.infection. == 'Yes'))
kruskal_test(formula = IC50_neutAb_log10 ~ timeCategory, data=subset(subsetData, Prior.COVID.infection. == 'Yes')) 
dunn_test(IC50_neutAb_log10 ~ timeCategory, data=subset(subsetData, Prior.COVID.infection. == 'Yes'))

linePlot(data = subsetData, xData = 'timeCategory', yData = 'IC50_neutAb_log10', groupby = 'Alias', xLabel = ' ', yLabel = "log10 IC50", 
         title = "Neutralizing antibodies", colorby = "Prior.COVID.infection.", recentCOVID = T) + 
  geom_hline(yintercept = 10, linetype = "dashed",alpha=0.3) + annotate("text", x=5,y=15,label = "LOD", color="black", alpha=0.2)+
  scale_y_continuous(trans='pseudo_log', limits = c(0,1e5), breaks=c(10^(0:6)), labels=trans_format('log10',math_format(10^.x)), minor_breaks =5*10^(0:10))
# ggsave(filename = "./Images/neutAb_linePlot_recentCOVID.pdf")



##' Post 1st dose
# twoSampleBar(data = subset(mergedData, timeCategory == "Post 1st dose"), xData = "Prior.COVID.infection.", yData = "IC50_neutAb_log10", fillParam = "Prior.COVID.infection.", 
#              title = "Post 1st dose", yLabel = "Neutralizing antibodies", nonparam = T) +   coord_cartesian(ylim=c(0e1,1e5),) + 
#   scale_y_continuous(trans='pseudo_log', breaks=c(10^(0:7)), labels=trans_format('log10',math_format(10^.x)), minor_breaks =5*10^(0:10)) #+ ggrepel::geom_text_repel(aes(label = Label))

##' Post 2nd dose
twoSampleBar(data = subset(mergedData, timeCategory == "Post 2nd dose"), xData = "Prior.COVID.infection.", yData = "IC50_neutAb_log10", fillParam = "Prior.COVID.infection.", 
             title = "Post 2nd dose", yLabel = "Neutralizing antibodies", nonparam = T) +   coord_cartesian(ylim=c(0e1,1e5),) + 
  scale_y_continuous(trans='pseudo_log', breaks=c(10^(0:7)), labels=trans_format('log10',math_format(10^.x)), minor_breaks =5*10^(0:10)) #+ ggrepel::geom_text_repel(aes(label = Label))
# ggsave(filename = "./Images/neutAb_post2ndDose.pdf", width=5)

twoSampleBar(data = subset(mergedData, timeCategory == "One month post\n2nd dose"), xData = "Prior.COVID.infection.", yData = "IC50_neutAb_log10", fillParam = "Prior.COVID.infection.", 
             title = "One month post", yLabel = "Neutralizing antibodies", nonparam = T) +   coord_cartesian(ylim=c(0e1,1e5),) + 
  scale_y_continuous(trans='pseudo_log', breaks=c(10^(0:7)), labels=trans_format('log10',math_format(10^.x)), minor_breaks =5*10^(0:10)) #+ ggrepel::geom_text_repel(aes(label = Label))
# ggsave(filename = "./Images/neutAb_oneMonthpost2ndDose.pdf", width=5)


twoSampleBar(data = subset(mergedData, timeCategory == "Baseline"), xData = "Prior.COVID.infection.", yData = "binding_IgG_N", fillParam = "Prior.COVID.infection.", 
             title = "Baseline", yLabel = "Anti-N IgG titer", nonparam = T) + coord_cartesian(ylim=c(0e1,1e5),) + 
  scale_y_continuous(trans='pseudo_log', breaks=c(10^(0:7)), labels=trans_format('log10',math_format(10^.x)), minor_breaks =5*10^(0:10)) #+ ggrepel::geom_text_repel(aes(label = Label))
# ggsave(filename = "./Images/bindingAb_N_IgG_baseline.pdf", width=5)

linePlot(data = subsetData, xData = 'timeCategory', yData = 'binding_IgG_N', groupby = 'Alias', xLabel = ' ', yLabel = "anti-N IgG titer", 
         title = "anti-N IgG titer", colorby = "Prior.COVID.infection.") + 
  # scale_color_manual(name="Prior COVID?",values = c("#FFC26A","#B5B2F1")) + #theme(legend.position = 'none') + 
  scale_y_continuous(trans='pseudo_log', limits = c(0,1e5), breaks=c(10^(0:6)), labels=trans_format('log10',math_format(10^.x)), minor_breaks =5*10^(0:10)) #+ ggrepel::geom_text_repel(aes(label = Label))
# ggsave(filename = "./Images/BindingAb_N_IgG_linePlot.pdf", width=7)


linePlot(data = mergedData, xData = 'DPV', yData = 'binding_IgG_S1', groupby = 'Alias', xLabel = 'Days relative to first dose', yLabel = "anti-S1 IgG titer", 
         title = "anti-S1 IgG titer", colorby = "Prior.COVID.infection.") + 
  # scale_color_manual(name="Prior COVID?",values = c("#FFC26A","#B5B2F1")) + 
  scale_y_continuous(trans='pseudo_log', limits = c(0,1e7), breaks=c(10^(0:7)), labels=trans_format('log10',math_format(10^.x)), minor_breaks =5*10^(0:10))# + ggrepel::geom_text_repel(aes(label = Label))
# ggsave(filename = "./Images/BindingAb_S1_IgG_linePlot_contTime.pdf", width=7)


univScatter(data = subset(mergedData, Prior.COVID.infection. == "Yes" & timeCategory == "Baseline"), xData = "binding_IgG_S1", yData = 'FC_IgG_S1_postVax1', 
            fillParam = 'Prior.COVID.infection.',title = "Post 1st dose Anti-S1 IgG", xLabel = "Baseline anti-S1 IgG titer", yLabel = "Fold-change anti-S1 IgG", nonparam = T) + 
  scale_fill_manual(values=c("#B5B2F1")) + 
  scale_x_continuous(trans='pseudo_log', limits = c(1e2,5e5), breaks=c(10^(0:7)), labels=trans_format('log10',math_format(10^.x)), minor_breaks =5*10^(0:10)) +
  scale_y_continuous(trans='pseudo_log', limits = c(0,1e5), breaks=c(10^(0:7)), labels=trans_format('log10',math_format(10^.x)), minor_breaks =5*10^(0:10))
# ggsave(filename = "./Images/BindingAb_S1_IgG_vs_FC_S1binding.pdf")


subsetData <- subset(mergedData, timeCategory == 'Baseline'); subsetData <- subsetData[which(!is.na(subsetData$binding_IgG_S1)),]
subsetData %>% levene_test( binding_IgG_S1 ~ Prior.COVID.infection.)
subsetData %>% shapiro_test( binding_IgG_S1)
subsetData %>%  wilcox_test(binding_IgG_S1 ~ Prior.COVID.infection.)
subsetData %>% group_by(Prior.COVID.infection.) %>% get_summary_stats(binding_IgG_S1)




subsetData <- subset(mergedData, timeCategory == 'Post 2nd dose'); subsetData <- subsetData[which(!is.na(subsetData$binding_IgG_S1)),]
subsetData %>% levene_test( binding_IgG_S1 ~ Prior.COVID.infection.)
subsetData %>% shapiro_test( binding_IgG_S1)
subsetData %>%  wilcox_test(binding_IgG_S1 ~ Prior.COVID.infection.)

subsetData <- subset(mergedData, timeCategory == 'One month post\n2nd dose'); subsetData <- subsetData[which(!is.na(subsetData$binding_IgG_S1)),]
subsetData %>% group_by(Prior.COVID.infection.) %>% get_summary_stats(binding_IgG_S1)
subsetData %>% levene_test( binding_IgG_S1 ~ Prior.COVID.infection.)
subsetData %>% shapiro_test( binding_IgG_S1)
subsetData %>%  wilcox_test(binding_IgG_S1 ~ Prior.COVID.infection.)


subsetData %>% group_by(Prior.COVID.infection.) %>% get_summary_stats(FC_IgG_S1_postVax1)
subsetData %>% shapiro_test(FC_IgG_S1_postVax1)
subsetData %>%  wilcox_test(FC_IgG_S1_postVax1 ~ Prior.COVID.infection.)
subsetData %>% group_by(Prior.COVID.infection.) %>% get_summary_stats(FC_IgG_S1_postVax2)
subsetData %>% shapiro_test(FC_IgG_S1_postVax2)
subsetData %>%  wilcox_test(FC_IgG_S1_postVax2 ~ Prior.COVID.infection.)



subsetData <- subset(mergedData, timeCategory == 'Baseline'); subsetData <- subsetData[which(!is.na(subsetData$binding_IgG_S1)),]
subsetData %>% group_by(Prior.COVID.infection.) %>% get_summary_stats(IC50_neutAb_log10)
subsetData %>% levene_test( IC50_neutAb_log10 ~ Prior.COVID.infection.)
subsetData %>% shapiro_test( IC50_neutAb_log10)
subsetData %>%  wilcox_test(IC50_neutAb_log10 ~ Prior.COVID.infection.)

subsetData <- subset(mergedData, timeCategory == 'Post 1st dose'); subsetData <- subsetData[which(!is.na(subsetData$binding_IgG_S1)),]
subsetData %>% group_by(Prior.COVID.infection.) %>% get_summary_stats(IC50_neutAb_log10)
subsetData %>% levene_test( IC50_neutAb_log10 ~ Prior.COVID.infection.)
subsetData %>% shapiro_test( IC50_neutAb_log10)
subsetData %>%  wilcox_test(IC50_neutAb_log10 ~ Prior.COVID.infection.)

subsetData <- subset(mergedData, timeCategory == 'Post 2nd dose'); subsetData <- subsetData[which(!is.na(subsetData$binding_IgG_S1)),]
subsetData %>% group_by(Prior.COVID.infection.) %>% get_summary_stats(IC50_neutAb_log10)
subsetData %>% levene_test( IC50_neutAb_log10 ~ Prior.COVID.infection.)
subsetData %>% shapiro_test( IC50_neutAb_log10)
subsetData %>%  wilcox_test(IC50_neutAb_log10 ~ Prior.COVID.infection.)



subsetData <- mergedData[which(!is.na(mergedData$FC_Elispot_IgG_S1)),]
subsetData <- subsetData[which(is.finite(subsetData$FC_Elispot_IgG_S1)),]
subsetData %>% group_by(Prior.COVID.infection.) %>% get_summary_stats(FC_Elispot_IgG_S1)



#' ------------------ Correlational analyses --------------------------
#' 

bivScatter(data1 = subset(mergedData, Prior.COVID.infection. == 'No' & timeCategory == "Post 1st dose" & Tube == "HEP"), name1 = "Naive", 
           data2 = subset(mergedData, Prior.COVID.infection. == 'Yes' & timeCategory == "Post 1st dose" & Tube == "HEP"), name2 = "Experienced", 
           xData = "FCActivCD4_Vax1", yData = "FCActivCD8_Vax1", fillParam = "Prior.COVID.infection.", title = "Post 1st dose", 
           xLabel = "Fold-change CD4+Ki67+CD38+", yLabel = "Fold-change CD8+Ki67+CD38+", nonparam = T) + 
  scale_x_continuous(trans='pseudo_log', breaks=c(10^(0:4)), labels=trans_format('log10',math_format(10^.x)), minor_breaks =5*10^(0:10)) + 
  scale_y_continuous(trans='pseudo_log', breaks=c(10^(0:4)), labels=trans_format('log10',math_format(10^.x)), minor_breaks =5*10^(0:10))
# ggsave(filename = "./Images/FCActivCD4_correl_FCActivCD8_Vax1.pdf", width=8)
bivScatter(data1 = subset(mergedData, Prior.COVID.infection. == 'No' & timeCategory == "Post 2nd dose" & Tube == "HEP"), name1 = "Naive", 
           data2 = subset(mergedData, Prior.COVID.infection. == 'Yes' & timeCategory == "Post 2nd dose" & Tube == "HEP"), name2 = "Experienced", 
           xData = "FCActivCD4_Vax2", yData = "FCActivCD8_Vax2", fillParam = "Prior.COVID.infection.", title = "Post 2nd dose", 
           xLabel = "Fold-change CD4+Ki67+CD38+", yLabel = "Fold-change CD8+Ki67+CD38+", nonparam = T) + 
  scale_x_continuous(trans='pseudo_log', breaks=c(10^(0:4)), labels=trans_format('log10',math_format(10^.x)), minor_breaks =5*10^(0:10)) + 
  scale_y_continuous(trans='pseudo_log', limits = c(0,200), breaks=c(10^(0:4)), labels=trans_format('log10',math_format(10^.x)), minor_breaks =5*10^(0:10))
# ggsave(filename = "./Images/FCActivCD4_correl_FCActivCD8_Vax2.pdf", width=8)


subsetData <- subset(mergedData, Prior.COVID.infection. == 'No' & timeCategory == "Post 1st dose")
subsetData <- subsetData[,grep("^FCActiv",names(subsetData))]
cor.elispots <- cor(subsetData, method="kendall" , use="pairwise.complete.obs" )
cor.elispots.pmat <- ggcorrplot::cor_pmat(subsetData, method="kendall", use="pairwise.complete.obs"  )
ggcorrplot::ggcorrplot(corr = cor.elispots, p.mat = cor.elispots.pmat, title = "Fold-change responses - Naive", legend.title = "Kendall tau", insig = "blank", tl.cex = 20)
# ggsave(filename = "./Images/Fold-changes_ggcorrplot_Naive.pdf")

subsetData <- subset(mergedData, Prior.COVID.infection. == 'Yes' & timeCategory == "Post 1st dose")
subsetData <- subsetData[,grep("^FCActiv",names(subsetData))]
cor.elispots <- cor(subsetData, method="kendall" , use="pairwise.complete.obs" )
cor.elispots.pmat <- ggcorrplot::cor_pmat(subsetData, method="kendall", use="pairwise.complete.obs"  )
ggcorrplot::ggcorrplot(corr = cor.elispots, p.mat = cor.elispots.pmat, title = "Fold-change responses - Experienced", legend.title = "Kendall tau", insig = "blank", tl.cex = 20)
# ggsave(filename = "./Images/Fold-changes_ggcorrplot_Experienced.pdf")


bivScatter(data1 = subset(mergedData, Prior.COVID.infection. == 'No' & timeCategory == "Post 1st dose" & Tube == "HEP"), name1 = "Naive", 
           data2 = subset(mergedData, Prior.COVID.infection. == 'Yes' & timeCategory == "Post 1st dose" & Tube == "HEP"), name2 = "Experienced", 
           xData = "Age", yData = "CD4_.CD38.Ki67._FreqParent", fillParam = "Prior.COVID.infection.", title = "Post 1st dose", 
           xLabel = "Age", yLabel = "Ki67+CD38+ (% CD4)", nonparam = T) + 
  scale_x_continuous(limits = c(20,70), breaks=seq(0,100,10)) + scale_y_continuous(limits = c(0,5), breaks=seq(0,10,1))
# ggsave(filename = "./Images/ActivCD4_correl_Age_Vax1.pdf", width=8)

bivScatter(data1 = subset(mergedData, Prior.COVID.infection. == 'No' & timeCategory == "Post 2nd dose" & Tube == "HEP"), name1 = "Naive", 
           data2 = subset(mergedData, Prior.COVID.infection. == 'Yes' & timeCategory == "Post 2nd dose" & Tube == "HEP"), name2 = "Experienced", 
           xData = "Age", yData = "CD4_.CD38.Ki67._FreqParent", fillParam = "Prior.COVID.infection.", title = "Post 2nd dose", 
           xLabel = "Age", yLabel = "Ki67+CD38+ (% CD4)", nonparam = T) + 
  scale_x_continuous(limits = c(20,70), breaks=seq(0,100,10)) + scale_y_continuous(limits = c(0,5), breaks=seq(0,10,1))
# ggsave(filename = "./Images/ActivCD4_correl_Age_Vax2.pdf", width=8)


bivScatter(data1 = subset(mergedData, Prior.COVID.infection. == 'No' & timeCategory == "Post 1st dose" & Tube == "HEP"), name1 = "Naive", 
           data2 = subset(mergedData, Prior.COVID.infection. == 'Yes' & timeCategory == "Post 1st dose" & Tube == "HEP"), name2 = "Experienced", 
           xData = "Age", yData = "CD8_.CD38.Ki67._FreqParent", fillParam = "Prior.COVID.infection.", title = "Post 1st dose", 
           xLabel = "Age", yLabel = "Ki67+CD38+ (% CD8)", nonparam = T) + 
  scale_x_continuous(limits = c(20,70), breaks=seq(0,100,10)) + scale_y_continuous(limits = c(0,5), breaks=seq(0,10,1))
# ggsave(filename = "./Images/ActivCD8_correl_Age_Vax1.pdf", width=8)

bivScatter(data1 = subset(mergedData, Prior.COVID.infection. == 'No' & timeCategory == "Post 2nd dose" & Tube == "HEP"), name1 = "Naive", 
           data2 = subset(mergedData, Prior.COVID.infection. == 'Yes' & timeCategory == "Post 2nd dose" & Tube == "HEP"), name2 = "Experienced", 
           xData = "Age", yData = "CD8_.CD38.Ki67._FreqParent", fillParam = "Prior.COVID.infection.", title = "Post 2nd dose", 
           xLabel = "Age", yLabel = "Ki67+CD38+ (% CD8)", nonparam = T) + 
  scale_x_continuous(limits = c(20,70), breaks=seq(0,100,10)) + scale_y_continuous(limits = c(0,5), breaks=seq(0,10,1))
# ggsave(filename = "./Images/ActivCD8_correl_Age_Vax2.pdf", width=8)



#'  #------------------------------------ Age correlations with activated CD4 responses ------------------------------------------

subsetData <- subset(mergedData, Prior.COVID.infection. == 'No' & timeCategory == "Post 1st dose")
subsetData <- subsetData[,grep( paste(c("CD8_.CD38.Ki67._FreqParent","CD4_.CD38.Ki67._FreqParent","Age"), collapse = "|"), names(subsetData))]  #  "^FCActivCD4",
cor.matrix <- cor(subsetData, method="kendall" , use="pairwise.complete.obs" )
cor.matrix.pmat <- ggcorrplot::cor_pmat(subsetData, method="kendall", use="pairwise.complete.obs"  )
cor.matrix <- as.data.frame(cor.matrix); cor.matrix$Labels <- row.names(cor.matrix); cor.matrix$Prior.COVID <- "No"
cor.matrix <- merge( x = cor.matrix, y = cor.matrix.pmat[,"Age"], by = "row.names"); names(cor.matrix)[grep("y",names(cor.matrix))] <- "Pvalue"

subsetData <- subset(mergedData, Prior.COVID.infection. == 'Yes' & timeCategory == "Post 1st dose")
subsetData <- subsetData[,grep(  paste(c("CD8_.CD38.Ki67._FreqParent","CD4_.CD38.Ki67._FreqParent","Age" ), collapse = "|"), names(subsetData))]  #"^FCActivCD4",
cor.matrix2 <- cor(subsetData, method="kendall" , use="pairwise.complete.obs" )
cor.matrix2.pmat <- ggcorrplot::cor_pmat(subsetData, method="kendall", use="pairwise.complete.obs"  )
cor.matrix2 <- as.data.frame(cor.matrix2); cor.matrix2$Labels <- row.names(cor.matrix2); cor.matrix2$Prior.COVID <- "Yes"
cor.matrix2 <- merge( x = cor.matrix2, y = cor.matrix2.pmat[,"Age"], by = "row.names"); names(cor.matrix2)[grep("y",names(cor.matrix2))] <- "Pvalue"

temp <- as.data.frame(rbind(cor.matrix, cor.matrix2)); temp <- temp[ -grep( paste( c("Age"), collapse = "|"), temp$Row.names),]
# 
# ggplot( data = temp, aes(y = Labels,x = Age, fill = Prior.COVID)) + geom_bar(stat='identity',position = 'dodge',width=0.75) + theme_bw() + 
#   scale_fill_manual(values=c("#FFDFB1", "#B5B2F1")) + xlab("Correlation with Age") + ylab(" ") + theme(axis.text.y = element_text(angle=0, size = 10)) + 
#   ggtitle("Post 1st dose") + geom_vline(xintercept=0, linetype = "dashed") + scale_x_continuous(limits = c(-1,1))
#   

subsetData <- subset(mergedData, Prior.COVID.infection. == 'No' & timeCategory == "Post 2nd dose")
subsetData <- subsetData[,grep( paste(c("CD8_.CD38.Ki67._FreqParent","CD4_.CD38.Ki67._FreqParent","Age" ), collapse = "|"), names(subsetData))]
cor.matrix <- cor(subsetData, method="kendall" , use="pairwise.complete.obs" )
cor.matrix.pmat <- ggcorrplot::cor_pmat(subsetData, method="kendall", use="pairwise.complete.obs"  )
cor.matrix <- as.data.frame(cor.matrix); cor.matrix$Labels <- row.names(cor.matrix); cor.matrix$Prior.COVID <- "No"
cor.matrix <- merge( x = cor.matrix, y = cor.matrix.pmat[,"Age"], by = "row.names"); names(cor.matrix)[grep("y",names(cor.matrix))] <- "Pvalue"

subsetData <- subset(mergedData, Prior.COVID.infection. == 'Yes' & timeCategory == "Post 2nd dose")
subsetData <- subsetData[,grep( paste(c("CD8_.CD38.Ki67._FreqParent","CD4_.CD38.Ki67._FreqParent","Age" ), collapse = "|"), names(subsetData))]
cor.matrix2 <- cor(subsetData, method="kendall" , use="pairwise.complete.obs" )
cor.matrix2.pmat <- ggcorrplot::cor_pmat(subsetData, method="kendall", use="pairwise.complete.obs"  )
cor.matrix2 <- as.data.frame(cor.matrix2); cor.matrix2$Labels <- row.names(cor.matrix2); cor.matrix2$Prior.COVID <- "Yes"
cor.matrix2 <- merge( x = cor.matrix2, y = cor.matrix2.pmat[,"Age"], by = "row.names"); names(cor.matrix2)[grep("y",names(cor.matrix2))] <- "Pvalue"


temp2 <- as.data.frame(rbind(cor.matrix, cor.matrix2)); temp2 <- temp2[ -grep( paste( c("Age"), collapse = "|"), temp2$Row.names),]


# temp[ grep("CD4_.CD38.Ki67._FreqParent", temp$Row.names), "Labels"] <- "Frequency\npost 1st dose"
temp[,"Labels"]  <- "Frequency\npost 1st dose";  temp2[,"Labels"]  <- "Frequency\npost 2nd dose"
# temp2[ grep("CD4_.CD38.Ki67._FreqParent", temp2$Row.names), "Labels"] <- "Frequency\npost 2nd dose"
# temp[ grep("FCActivCD4_Vax1",temp$Row.names), "Labels"]  <- "Fold-change\npost 1st dose"
# temp[ grep("FCActivCD4_Vax2",temp$Row.names), "Labels"]  <- "Fold-change\npost 2nd dose"
temp <- temp[,c("Row.names","Age","Labels","Prior.COVID","Pvalue")];  temp2 <- temp2[,c("Row.names", "Age","Labels","Prior.COVID","Pvalue")]  
temp <- as.data.frame(rbind(temp, temp2))
temp$Labels <- factor(temp$Labels, levels = c("Fold-change\npost 2nd dose", "Fold-change\npost 1st dose","Frequency\npost 2nd dose",   "Frequency\npost 1st dose"))
# temp <- temp[c(1:6,8,7),]   # rearrange last two rows;  needed if FC lines are included in the df
a <- ggplot( data = subset(temp, Row.names == "CD4_.CD38.Ki67._FreqParent"), aes(x = Labels,y = Age, fill = Prior.COVID)) + geom_bar(stat='identity',position = position_dodge(width=0.5), width=0.05, color="black", size=0.1) + 
  geom_point(aes(fill=Prior.COVID, size=Pvalue), pch=21, color="black", stroke=0.2, position = position_dodge(width=0.5)) + 
  theme_bw() + scale_size(range = c(18,1), breaks = c(0,0.05,0.1,0.2,0.7), limits = c(0,0.8), trans = 'pseudo_log') + guides(size = guide_legend(reverse=TRUE)) + 
  scale_fill_manual(values=c("#FFC26A", "#B5B2F1")) + ylab("Kendall's tau vs Age") + xlab(" ") + ggtitle("CD4+Ki67+CD38+") + geom_hline(yintercept=0, linetype = "dashed") +
  theme(axis.text.y = element_text(size = 20, color="black"), plot.title = element_text(size=20), axis.text.x = element_text(size=20, color="black", angle=45, hjust=1,vjust=1), 
        axis.title.x = element_text(size=20, color="black")) + 
  coord_flip() + scale_y_continuous(limits = c(-1,0.5), breaks = seq(-1,1,0.25))
a
# ggsave(filename = "./Images/Age_CD4correlations_lollipop.pdf", width=5, height = 5)
b <- ggplot( data = subset(temp, Row.names == "CD8_.CD38.Ki67._FreqParent"), aes(x = Labels,y = Age, fill = Prior.COVID)) + geom_bar(stat='identity',position = position_dodge(width=0.5), width=0.05, color="black", size=0.1) + 
  geom_point(aes(fill=Prior.COVID, size=Pvalue), pch=21, color="black", stroke=0.2, position = position_dodge(width=0.5)) + 
  theme_bw() + scale_size(range = c(18,1), breaks = c(0,0.05,0.1,0.2,0.7), limits = c(0,0.8), trans = 'pseudo_log') + guides(size = guide_legend(reverse=TRUE)) + 
  scale_fill_manual(values=c("#FFC26A", "#B5B2F1")) + ylab("Kendall's tau vs Age") + xlab(" ") + ggtitle("CD8+Ki67+CD38+") + geom_hline(yintercept=0, linetype = "dashed") +
  theme(axis.text.y = element_text(size = 20, color="black"), plot.title = element_text(size=20), axis.text.x = element_text(size=20, color="black", angle=45, hjust=1,vjust=1), 
        axis.title.x = element_text(size=20, color="black"), legend.text = element_text(size=20), legend.title = element_text(size=20)) + 
  coord_flip() + scale_y_continuous(limits = c(-1,0.5), breaks = seq(-1,1,0.25))
b
# ggsave(filename = "./Images/Age_CD8correlations_lollipop.pdf", width=5, height = 5)
a <- a + theme(legend.position = "none");  b <- b+theme(axis.text.y = element_blank())
grid.arrange(a,b, nrow=1, widths = c(1,1.1))
# ggsave( plot = arrangeGrob(a,b, nrow=1, widths = c(1,1)), filename =  "./Images/Age_Tcellcorrelations_lollipop.pdf", width=9, height=6)



#' # ------------------ CXCL13 analyses --------------------------
#'
subsetData <- subset(mergedData, !is.na(mergedData$CXCL13))
subsetData %>% group_by( timeCategory) %>% get_summary_stats(CXCL13, type = "common") %>% print(n=500)
linePlot(data = mergedData, xData = 'timeCategory', yData = 'CXCL13', groupby = 'Alias', xLabel = ' ', yLabel = "CXCL13 (pg/mL)", 
         title = "Plasma CXCL13", colorby = "Prior.COVID.infection.") + theme(axis.title.x = element_blank()) + 
  scale_color_manual(name="Prior COVID?",values = c("#FFC26A","#B5B2F1")) + 
  scale_y_continuous(limits = c(0,140),breaks=seq(0,140,10))
  # scale_y_continuous(trans='pseudo_log', limits = c(0,1e4), breaks=c(10^(0:7)), labels=trans_format('log10',math_format(10^.x)), minor_breaks =5*10^(0:10))# + ggrepel::geom_text_repel(aes(label = Label))
# ggsave(filename = "./Images/CXCL13plasma_linePlot.pdf")
xData = "timeCategory";  yData = "CXCL13"; fillParam = "Prior.COVID.infection."; groupby = "Prior.COVID.infection."; title  = "CXCL13 as medians"; yLabel = "CXCL13 (pg/mL)"; xLabel = " "
subsetData <- subsetData[order(subsetData$Record.ID, subsetData$shortForm, decreasing = F),]
subsetData.median <- subsetData %>% group_by_('timeCategory', fillParam ) %>% summarize(median = median(.data[[yData]], na.rm=T))
       # colors:   COVID-exp B5B2F1 (purple)             COVID-naive FFC26A (orange)
    ggplot(data=subsetData, aes_string(x=xData, y=yData, fill=fillParam) ) + theme_bw() + 
      geom_path(aes_string(group=groupby, color=fillParam), alpha=0.4) + 
      geom_point(size = 2, pch=21, color="black", alpha=0.2) + facet_wrap(fillParam ) +   # , scales='free'
      scale_color_manual(values=c("#FFC26A", "#B5B2F1")) + 
      scale_fill_manual(values=c("#FFC26A", "#B5B2F1")) +  
      geom_line(data = subsetData.median, aes_string(x = 'timeCategory', y='median', group='1', color=fillParam), alpha=1,size=2) +
      geom_point(data = subsetData.median, aes_string(x = 'timeCategory', y='median'), alpha=1,size=2) +
      ggtitle(title) + ylab(yLabel) + xlab(xLabel)  +
      theme(axis.text = element_text(size=18,hjust = 0.5, color="black"), axis.title = element_text(size=22,hjust = 0.5), 
            plot.title = element_text(size=36,hjust = 0.5), axis.text.x = element_text(angle=45, hjust=1,vjust=1),
            legend.position = "none", strip.text = element_text(size = 24, color="black"), strip.background = element_rect(fill="white")) 



acuteCOVID <- read.csv(file = "D:/COVID/Infection/Analysis/Rscripts/COVIDmergedData.csv")
acuteCOVID$dummy <- "Acute COVID "
# acuteCOVID.cxcl13 <- acuteCOVID[which(!is.na(acuteCOVID$CXCL13..pg.mL.)),]
# t(acuteCOVID.cxcl13$Subject %in% mergedData$Alias)
# acuteCOVID.cxcl13[c(7,20),c("Subject")]               # too few of the acute covid samples' cxcl13 overlap with the vaccinated cxcl13
ggplot(data = subset(acuteCOVID, DPO<30), aes(x = dummy , y = CXCL13..pg.mL.)) + ggbeeswarm::geom_quasirandom( alpha=0.4, color="black", size=3) + 
  ggtitle("COVID-19") + ylab("Plasma CXCL13 (pg/mL)") + theme_bw() +
  theme(axis.text = element_text(color="black",size=18), axis.title = element_text(size=24), axis.text.x = element_text(angle=45, hjust=1,vjust=1),
        plot.title=element_text(size=28), axis.title.x = element_blank()) + 
  scale_y_continuous(breaks = seq(0,150,10), limits = c(0,140))
# ggsave( filename = "./Images/CXCL13plasma_acuteCOVID.pdf", width=3)


twoSampleBar(data = subset(mergedData, timeCategory == "Post 1st dose"), xData = "Prior.COVID.infection.", yData = "CXCL13", fillParam = "Prior.COVID.infection.", 
             title = "CXCL13 after 1st dose", yLabel = "CXCL13 (pg/mL)", batch="none", position = "left", FCplot=F, confInt=F, nonparam=T)
twoSampleBar(data = subset(mergedData, timeCategory == "Post 2nd dose"), xData = "Prior.COVID.infection.", yData = "CXCL13", fillParam = "Prior.COVID.infection.", 
             title = "CXCL13 after 2nd dose", yLabel = "CXCL13 (pg/mL)", batch="none", position = "left", FCplot=F, confInt=F, nonparam=T)


#' #------------------ Avidity analyses --------------------------
#'
subsetData <- subset(mergedData,  timeCategory != "two Weeks" & timeCategory != "2 wks post 2nd dose");   
subsetData$timeCategory <- factor(subsetData$timeCategory, levels = c("Baseline", "Post 1st dose", "Pre 2nd dose", "Post 2nd dose","One month post\n2nd dose"))
a<- linePlot(data = subsetData, xData = 'timeCategory', yData = 'Avidity', groupby = 'Alias', xLabel = " ", yLabel = "IgG avidity (%)", 
         title = "anti-S1 Avidity", colorby = "Prior.COVID.infection.") + theme(axis.title.x = element_blank()) + 
  scale_color_manual(name="Prior COVID?",values = c("#FFC26A","#B5B2F1")) + scale_y_continuous(limits = c(0,110),breaks=seq(0,140,10)) #+ 
#  ggrepel::geom_text_repel(data = subset(subsetData, Alias =="HV-002"), aes(label = Label))
a
# ggsave(filename = "./Images/avidity_lineplot.pdf")
# plotly::ggplotly(a)

bartlett.test(Avidity ~ timeCategory, data=subset(subsetData, Prior.COVID.infection. == 'No'))
tukey_hsd( aov(Avidity ~ timeCategory, data=subset(subsetData, Prior.COVID.infection. == 'No')) )
# dunn_test(Avidity ~ timeCategory, data=subset(subsetData, Prior.COVID.infection. == 'No'))
bartlett.test(Avidity ~ timeCategory, data=subset(subsetData, Prior.COVID.infection. == 'Yes'))
dunn_test(Avidity ~ timeCategory, data=subset(subsetData, Prior.COVID.infection. == 'Yes'))
# tukey_hsd( aov(Avidity ~ timeCategory, data=subset(subsetData, Prior.COVID.infection. == 'Yes')) )

linePlot(data = subsetData, xData = 'timeCategory', yData = 'Avidity', groupby = 'Alias', xLabel = ' ', yLabel = "IgG avidity (%)", 
         title = "anti-S1 Avidity", colorby = "Prior.COVID.infection.", recentCOVID = T) + theme(axis.title.x = element_blank()) + 
  scale_color_manual(name="Prior COVID?",values = c("grey60","grey60")) + scale_y_continuous(limits = c(0,110),breaks=seq(0,140,10)) 
# ggsave(filename = "./Images/avidity_lineplot_recentCOVID.pdf")


#' ------------------ Correlational analyses --------------------------
#' 
# index <- grep(paste(c("^Prior","DPO.covid","DPV","timeCategory", "FreqParent$", "^FC"), collapse = "|"), names(mergedData), value=F)
# mergedData.deident <- mergedData[, index]
# subsetData <- mergedData.deident[which(mergedData.deident$timeCategory=="Post 1st dose"), ]
# GGally::ggduo( data = subsetData, mapping=ggplot2::aes(color=Prior.COVID.infection.))+ theme_bw() + ggplot2::theme(legend.position = "bottom") +    scale_color_manual(values=c("#FFDFB1", "#B5B2F1")) + 
#   scale_fill_manual(values=c("#FFDFB1", "#B5B2F1")) 
# # ggsave("./Images/exploratory_post1st.pdf", width=50, height=50, limitsize = F)
# 
# subsetData <- mergedData.deident[which(mergedData.deident$timeCategory=="Post 2nd dose"), ]
# GGally::ggduo( data = subsetData, mapping=ggplot2::aes(color=Prior.COVID.infection.))+ theme_bw() + ggplot2::theme(legend.position = "bottom") +    scale_color_manual(values=c("#FFDFB1", "#B5B2F1")) + 
#   scale_fill_manual(values=c("#FFDFB1", "#B5B2F1")) 
# # ggsave("./Images/exploratory_post2nd.pdf", width=50, height=50, limitsize = F)

ggplot(data = subset(mergedData, timeCategory == "Post 1st dose"), aes(x = Elispot_IgM_RBD, y=CD19_.CD27..CD38._FreqParent, color=Prior.COVID.infection. )) + geom_point(size=5) + theme_bw()
ggplot(data = subset(mergedData, timeCategory == "Post 1st dose"), aes(x = Elispot_IgG_S1, y=FCActivCD4_Vax1, color=Prior.COVID.infection. )) + geom_point(size=5) + theme_bw()




##'  ------------------------------------ Age correlations with Tfh responses ------------------------------------------
##'  


subsetData <- subset(mergedData, timeCategory == "Post 1st dose")
subsetData <- subsetData[,grep( paste(c("^FC"), collapse = "|"), names(subsetData))]
cor.matrix <- round(cor(subsetData, method="kendall",use="pairwise.complete.obs"), 2)

bivScatter(data1 = subset(mergedData, Prior.COVID.infection. == 'No' & timeCategory == "Post 1st dose" ), name1 = "Naive", 
           data2 = subset(mergedData, Prior.COVID.infection. == 'Yes' & timeCategory == "Post 1st dose"), name2 = "Experienced", 
           xData = "FC_IgG_S1_postVax1", yData = "FCtfh_CXCR3_Vax1", fillParam = "Prior.COVID.infection.", title = "Post 1st dose", 
           xLabel = "Fold-change CD4+Ki67+CD38+", yLabel = "Fold-change ICOS+CD38+ cTfh", nonparam = T) + 
  scale_x_continuous(trans='pseudo_log', breaks=c(10^(0:4)), labels=trans_format('log10',math_format(10^.x)), minor_breaks =5*10^(0:10)) + 
  scale_y_continuous(trans='pseudo_log', breaks=c(10^(0:4)), labels=trans_format('log10',math_format(10^.x)), minor_breaks =5*10^(0:10))



subsetData <- subset(mergedData, Prior.COVID.infection. == 'No' & timeCategory == "Post 1st dose")
subsetData <- subsetData[,grep( paste(c("CD4_.Nonnaive.cTfh.ICOS..CD38.._FreqParent","^FCtfh","Age" ), collapse = "|"), names(subsetData))]
cor.matrix <- cor(subsetData, method="kendall" , use="pairwise.complete.obs" )
cor.matrix.pmat <- ggcorrplot::cor_pmat(subsetData, method="kendall", use="pairwise.complete.obs"  )
cor.matrix <- as.data.frame(cor.matrix); cor.matrix$Labels <- row.names(cor.matrix); cor.matrix$Prior.COVID <- "No"
cor.matrix <- merge( x = cor.matrix, y = cor.matrix.pmat[,"Age"], by = "row.names"); names(cor.matrix)[grep("y",names(cor.matrix))] <- "Pvalue"

subsetData <- subset(mergedData, Prior.COVID.infection. == 'Yes' & timeCategory == "Post 1st dose")
subsetData <- subsetData[,grep( paste(c("CD4_.Nonnaive.cTfh.ICOS..CD38.._FreqParent","^FCtfh","Age" ), collapse = "|"), names(subsetData))]
cor.matrix2 <- cor(subsetData, method="kendall" , use="pairwise.complete.obs" )
cor.matrix2.pmat <- ggcorrplot::cor_pmat(subsetData, method="kendall", use="pairwise.complete.obs"  )
cor.matrix2 <- as.data.frame(cor.matrix2); cor.matrix2$Labels <- row.names(cor.matrix2); cor.matrix2$Prior.COVID <- "Yes"
cor.matrix2 <- merge( x = cor.matrix2, y = cor.matrix2.pmat[,"Age"], by = "row.names"); names(cor.matrix2)[grep("y",names(cor.matrix2))] <- "Pvalue"

temp <- as.data.frame(rbind(cor.matrix, cor.matrix2)); temp <- temp[ -grep( paste( c("Age", "Foxp3", "CXCR3"), collapse = "|"), temp$Row.names),]
# 
# ggplot( data = temp, aes(y = Labels,x = Age, fill = Prior.COVID)) + geom_bar(stat='identity',position = 'dodge',width=0.75) + theme_bw() +
#   scale_fill_manual(values=c("#FFDFB1", "#B5B2F1")) + xlab("Correlation with Age") + ylab(" ") + theme(axis.text.y = element_text(angle=0, size = 10)) +
#   ggtitle("Post 1st dose") + geom_vline(xintercept=0, linetype = "dashed") + scale_x_continuous(limits = c(-1,1))
#   

subsetData <- subset(mergedData, Prior.COVID.infection. == 'No' & timeCategory == "Post 2nd dose")
subsetData <- subsetData[,grep( paste(c("CD4_.Nonnaive.cTfh.ICOS..CD38.._FreqParent","Age" ), collapse = "|"), names(subsetData))]
cor.matrix <- cor(subsetData, method="kendall" , use="pairwise.complete.obs" )
cor.matrix.pmat <- ggcorrplot::cor_pmat(subsetData, method="kendall", use="pairwise.complete.obs"  )
cor.matrix <- as.data.frame(cor.matrix); cor.matrix$Labels <- row.names(cor.matrix); cor.matrix$Prior.COVID <- "No"
cor.matrix <- merge( x = cor.matrix, y = cor.matrix.pmat[,"Age"], by = "row.names"); names(cor.matrix)[grep("y",names(cor.matrix))] <- "Pvalue"

subsetData <- subset(mergedData, Prior.COVID.infection. == 'Yes' & timeCategory == "Post 2nd dose")
subsetData <- subsetData[,grep( paste(c("CD4_.Nonnaive.cTfh.ICOS..CD38.._FreqParent","Age" ), collapse = "|"), names(subsetData))]
cor.matrix2 <- cor(subsetData, method="kendall" , use="pairwise.complete.obs" )
cor.matrix2.pmat <- ggcorrplot::cor_pmat(subsetData, method="kendall", use="pairwise.complete.obs"  )
cor.matrix2 <- as.data.frame(cor.matrix2); cor.matrix2$Labels <- row.names(cor.matrix2); cor.matrix2$Prior.COVID <- "Yes"
cor.matrix2 <- merge( x = cor.matrix2, y = cor.matrix2.pmat[,"Age"], by = "row.names"); names(cor.matrix2)[grep("y",names(cor.matrix2))] <- "Pvalue"


temp2 <- as.data.frame(rbind(cor.matrix, cor.matrix2)); temp2 <- temp2[ -grep( paste( c("Age", "Foxp3", "CXCR3"), collapse = "|"), temp2$Row.names),]
# 
# ggplot( data = temp2, aes(y = Labels,x = Age, fill = Prior.COVID)) + geom_bar(stat='identity',position = 'dodge',width=0.75) + theme_bw() + 
#   scale_fill_manual(values=c("#FFDFB1", "#B5B2F1")) + xlab("Correlation with Age") + ylab(" ") + theme(axis.text.y = element_text(angle=0, size = 10)) + 
#   ggtitle("Post 2nd dose") + geom_vline(xintercept=0, linetype = "dashed") +  scale_x_continuous(limits = c(-1,1))

# temp <- x; temp2 <- y
temp[ grep("CD4_.Nonnaive.cTfh.ICOS..CD38.._FreqParent", temp$Row.names), "Labels"] <- "Frequency\npost 1st dose"
temp2[ grep("CD4_.Nonnaive.cTfh.ICOS..CD38.._FreqParent", temp2$Row.names), "Labels"] <- "Frequency\npost 2nd dose"
temp[ grep("FCtfh_Vax1",temp$Row.names), "Labels"]  <- "Fold-change\npost 1st dose"
temp[ grep("FCtfh_Vax2",temp$Row.names), "Labels"]  <- "Fold-change\npost 2nd dose"
temp <- temp[,c("Row.names","Age","Labels","Prior.COVID","Pvalue")];  temp2 <- temp2[,c("Row.names", "Age","Labels","Prior.COVID","Pvalue")]  
temp <- as.data.frame(rbind(temp, temp2))
temp$Labels <- factor(temp$Labels, levels = c("Fold-change\npost 2nd dose", "Fold-change\npost 1st dose","Frequency\npost 2nd dose",   "Frequency\npost 1st dose"))
# temp <- temp[c(1:6,8,7),]
ggplot( data = temp, aes(x = Labels,y = Age, fill = Prior.COVID)) + geom_bar(stat='identity',position = position_dodge(width=0.5), width=0.05, color="black", size=0.1) + 
  geom_point(aes(fill=Prior.COVID, size=Pvalue), pch=21, color="black", stroke=0.2, position = position_dodge(width=0.5)) + 
  theme_bw() + scale_size(range = c(8,1), breaks = c(0,0.05,0.1,0.2,0.7), limits = c(0,0.8), trans = 'pseudo_log') + guides(size = guide_legend(reverse=TRUE)) + 
  scale_fill_manual(values=c("#FFC26A", "#B5B2F1")) + ylab("Kendall's tau vs Age") + xlab(" ") + ggtitle("ICOS+CD38+ cTfh") + geom_hline(yintercept=0, linetype = "dashed") +
  theme(axis.text.y = element_text(size = 16, color="black"), plot.title = element_text(size=24), axis.text.x = element_text(size=16, color="black", angle=45, hjust=1,vjust=1), 
        axis.title.x = element_text(size=16, color="black")) + 
  coord_flip() + scale_y_continuous(limits = c(-1,0.5), breaks = seq(-1,1,0.25))
# ggsave(filename = "./Images/Age_Tfhcorrelations_lollipop.pdf", width=5, height = 5)

bivScatter(data1 = subset(mergedData, Prior.COVID.infection. == 'No' & timeCategory == "Post 1st dose" ), name1 = "Naive", 
           data2 = subset(mergedData, Prior.COVID.infection. == 'Yes' & timeCategory == "Post 1st dose"), name2 = "Experienced", 
           xData = "FCActivCD4_Vax1", yData = "FCtfh_Vax1", fillParam = "Prior.COVID.infection.", title = "Post 1st dose", 
           xLabel = "Fold-change CD4+Ki67+CD38+", yLabel = "Fold-change ICOS+CD38+ cTfh", nonparam = T) + 
  scale_x_continuous(trans='pseudo_log', breaks=c(10^(0:4)), labels=trans_format('log10',math_format(10^.x)), minor_breaks =5*10^(0:10)) + 
  scale_y_continuous(trans='pseudo_log', breaks=c(10^(0:4)), limits = c(0,20), labels=trans_format('log10',math_format(10^.x)), minor_breaks =5*10^(0:10))
# ggsave(filename = "./Images/FCTfh_correl_FCactivCD4_Vax1.pdf", width=8)

bivScatter(data1 = subset(mergedData, Prior.COVID.infection. == 'No' & timeCategory == "Post 1st dose"), name1 = "Naive", 
           data2 = subset(mergedData, Prior.COVID.infection. == 'Yes' & timeCategory == "Post 1st dose"), name2 = "Experienced", 
           xData = "FCActivCD4_Vax2", yData = "FCtfh_Vax2", fillParam = "Prior.COVID.infection.", title = "Post 2nd dose", 
           xLabel = "Fold-change CD4+Ki67+CD38+", yLabel = "Fold-change ICOS+CD38+ cTfh", nonparam = T) + 
  scale_x_continuous(trans='pseudo_log', breaks=c(10^(0:4)), labels=trans_format('log10',math_format(10^.x)), minor_breaks =5*10^(0:10)) + 
  scale_y_continuous(trans='pseudo_log', breaks=c(10^(0:4)), labels=trans_format('log10',math_format(10^.x)), minor_breaks =5*10^(0:10))
# ggsave(filename = "./Images/FCTfh_correl_FCactivCD4_Vax2.pdf", width=8)




