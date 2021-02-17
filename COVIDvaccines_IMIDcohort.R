library("grid")
library("ggplot2")
library("gplots")
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



#' ------------------ Demographics and SubjectVisit --------------------------
#'

demog.IMID <- read.csv(file = "../../IMIDcohort/Analysis/Data/demog.csv")
demog.IMID$Alias <- NULL 
# demog.IMID <- demog.IMID[-which(demog.IMID$Record.ID == ""),]
demog.IMID.melt <- melt(demog.IMID, id.vars = c("Record.ID"), measure.vars = c("V1","V2","V3","V4","V5","V6","V7"))
demog.IMID.melt$Label <- paste0(demog.IMID.melt$Record.ID, "_", demog.IMID.melt$variable); 
demog.IMID.melt$Visit <- demog.IMID.melt$variable;  demog.IMID.melt$variable <- NULL
names(demog.IMID.melt)[grep("value",names(demog.IMID.melt), value=F)] <- "DPV"
demog.IMID.melt <- demog.IMID.melt[-which(is.na(demog.IMID.melt$DPV)), ]
# make range for classifying continuous days into categorical
x <- subset(demog.IMID.melt, DPV<100 );  hist(x$DPV , breaks=100)

demog.IMID.melt$shortForm <- demog.IMID.melt$timeCategory <- ""
demog.IMID.melt[which(demog.IMID.melt$DPV == 0), "shortForm"] <- "bL"; demog.IMID.melt[which(demog.IMID.melt$DPV == 0),"timeCategory"] <- "Baseline"

# -------- post 1st dose --------------  
demog.IMID.melt[which(demog.IMID.melt$DPV %in% 6:8 ),"shortForm"] <- "oW"; demog.IMID.melt[which(demog.IMID.melt$DPV %in% 6:8 ),"timeCategory"] <- "Post 1st dose"
# demog.IMID.melt[which(demog.IMID.melt$DPV %in% 6:8 & demog.IMID.melt$Visit == 'V2'),"shortForm"] <- "oW"; demog.IMID.melt[which(demog.IMID.melt$DPV %in% 6:8 & demog.IMID.melt$Visit == 'V2'),"timeCategory"] <- "Post 1st dose"

# -------- pre 2nd dose --------------  
demog.IMID.melt[which(demog.IMID.melt$DPV %in% 14 ),"shortForm"] <- "2W"; demog.IMID.melt[which(demog.IMID.melt$DPV %in% 14 ),"timeCategory"] <- "Two Weeks"
demog.IMID.melt[which(demog.IMID.melt$DPV %in% 19:21 ),"shortForm"] <- "3W"; demog.IMID.melt[which(demog.IMID.melt$DPV %in% 19:21 ),"timeCategory"] <- "Pre 2nd dose"

# -------- post 2nd dose --------------  
demog.IMID.melt[which(demog.IMID.melt$DPV %in% 28),"shortForm"] <- "4W"; demog.IMID.melt[which(demog.IMID.melt$DPV %in% 28),"timeCategory"] <- "Post 2nd dose"

demog.IMID.melt[which(demog.IMID.melt$DPV %in% -10),"shortForm"] <- "preCOVID"; demog.IMID.melt[which(demog.IMID.melt$DPV %in% -10),"timeCategory"] <- "preCOVID"
demog.IMID.melt[which(demog.IMID.melt$DPV %in% 50),"shortForm"] <- "postCOVID"; demog.IMID.melt[which(demog.IMID.melt$DPV %in% 50),"timeCategory"] <- "postCOVID"

#' ------------------ total binding ELISA  --------------------------
#'

totBinding.IMID <- read.csv("../../IMIDcohort/Analysis/Data/Scher_Serum _Normalized IgG EPT_02.09.2021.csv")
head(totBinding.IMID, n=10)

#' ------------------ NeutAb   --------------------------
#'

# neutAb <- read.csv("../Data/neutAb.csv")


#' ------------------ NeutAb   --------------------------
#'

# ASCelispot <- read.csv("../Data/ASCelispot.csv")



#' ------------------ Flow Cytometry --------------------------
#'


flowData.freq.IMID <- read.csv(file="../../IMIDcohort/Analysis/Data/ByParent.csv")                                 # flow cytometry ByParent spreadsheet
flowData.freq.IMID <- flowData.freq.IMID[-grep(paste(c("Mean","SD"),collapse = "|"), flowData.freq.IMID$X, value=F), ]         # delete the Mean and SD rows
names(flowData.freq.IMID)[1] <- "fcsFile"                # simple name

flowData.freq.IMID$fcsFile[1:19] <- c("ID06-0.fcs", "ID06-1W.fcs", "ID06-4W.fcs", "ID26-4W.fcs", "ID26-1W.fcs", "ID27-4W.fcs", 
                                      "ID27-1W.fcs", "ID337-0.fcs", "ID337-1W.fcs", "ID793-3W.fcs", "ID793-4W.fcs", "ID1400-4W.fcs",
                                      "ID1400-1W.fcs","SAGA337_preCOVID.fcs", "SAGA337_postCOVID.fcs", "SAGA1081_preCOVID.fcs", "SAGA1081_postCOVID.fcs", "SAGA1414_preCOVID.fcs", "SAGA1414_postCOVID.fcs")

a <- do.call(rbind.data.frame, strsplit(flowData.freq.IMID$fcsFile, paste(c("-", "_", "\\."), collapse="|")))      # split FCS file name
names(a) <- c("Record.ID","shortForm","Tube")

flowData.freq.IMID <- as.data.frame(cbind(a,flowData.freq.IMID))
flowData.freq.IMID$Tube[grep("Unmixed",flowData.freq.IMID$Tube)] <- "fcs"
flowData.freq.IMID[which(is.na(names(flowData.freq.IMID)))] <- NULL
names(flowData.freq.IMID) <- str_replace(names(flowData.freq.IMID),pattern="...Freq..of.Parent....", "_FreqParent")
names(flowData.freq.IMID) <- str_replace(names(flowData.freq.IMID),pattern="Lymphocytes.Single.Cells.Single.Cells.Live.CD3.CD4", "CD4_")
names(flowData.freq.IMID) <- str_replace(names(flowData.freq.IMID),pattern="Lymphocytes.Single.Cells.Single.Cells.Live.CD3.CD8", "CD8_")
names(flowData.freq.IMID) <- str_replace(names(flowData.freq.IMID),pattern="Lymphocytes.Single.Cells.Single.Cells.Live.CD19", "CD19_")
names(flowData.freq.IMID) <- str_replace(names(flowData.freq.IMID),pattern="Lymphocytes.Single.Cells.Single.Cells.Live.CD3", "CD3_")
flowData.freq.IMID$shortForm[which(flowData.freq.IMID$shortForm == '0')] <- "bL"
flowData.freq.IMID$shortForm[which(flowData.freq.IMID$shortForm == '1W')] <- "oW"

# flowData.freq.IMID$categoricalVisit[which(flowData.freq.IMID$categoricalVisit == "oD")] <- "bL"


#' ------------------ merge Data  --------------------------
#'

merge.1 <-  merge(x=demog.IMID, y=demog.IMID.melt, all=T, by=c("Record.ID"))
merge.2 <- merge(x = merge.1, y = flowData.freq.IMID, all=T, by=c("Record.ID","shortForm"))
merge.3 <- merge(x = merge.2, y = totBinding.IMID, all = T, by=c("Label"))


mergedData.IMID <- merge.3

mergedData.IMID$shortForm <- factor(mergedData.IMID$shortForm, levels = c("bL","oW","2W","3W","4W","5W"))
mergedData.IMID$timeCategory <- factor(mergedData.IMID$timeCategory, levels = c("Baseline","Post 1st dose","two Weeks", "Pre 2nd dose","Post 2nd dose","2 wks post 2nd dose", "One month post 2nd dose"))


# subsetData <- subset(mergedData, timeCategory == "Post 1st dose" | timeCategory == "Baseline")
# FC_response <- dcast(subsetData, `Record.ID`+`Alias` ~`timeCategory`, value.var = c("CD4_.Nonnaive.cTfh.ICOS..CD38.._FreqParent"))
# FC_response$FCtfh_Vax1 <- FC_response$`Post 1st dose`/FC_response$Baseline

# subsetData <- subset(mergedData, timeCategory == "Post 1st dose" | timeCategory == "Baseline")
# FC_response2 <- dcast( subsetData, `Record.ID`+`Alias`~`timeCategory`, value.var = c("CD19_..Nonnaive.B.CD27..CD38._FreqParent")) 
# FC_response2$FCPB_Vax1 <- FC_response2$`Post 1st dose`/FC_response2$Baseline; FC_response2$Cohort <- NULL
# FC_response <- merge(x=FC_response, y=FC_response2, by = c("Record.ID","Alias"))
# 
# subsetData <- subset(mergedData, timeCategory == "Post 2nd dose" | timeCategory == "Pre 2nd dose")
# FC_response2 <- dcast( subsetData, `Record.ID`+`Alias`~`timeCategory`, value.var = c("CD4_.Nonnaive.cTfh.ICOS..CD38.._FreqParent")) 
# FC_response2$FCtfh_Vax2 <- FC_response2$`Post 2nd dose`/FC_response2$`Pre 2nd dose`; FC_response2$Cohort <- NULL
# FC_response <- merge(x=FC_response, y=FC_response2, all = T, by = c("Record.ID","Alias"))
# 
# subsetData <- subset(mergedData, timeCategory == "Post 1st dose" | timeCategory == "Baseline" | timeCategory == "oneDay")
# FC_response2 <- dcast( subsetData, `Record.ID`+`Alias`~`timeCategory`, value.var = c("CD4_.CD38.Ki67._FreqParent")) 
# FC_response2$FCActivCD4_Vax1 <- FC_response2$`Post 1st dose`/FC_response2$Baseline; FC_response2$Cohort <- NULL
# FC_response <- merge(x=FC_response, y=FC_response2, all = T, by = c("Record.ID","Alias"))
# FC_response2 <- dcast( subsetData, `Record.ID`+`Alias`~`timeCategory`, value.var = c("CD8_..CD38.Ki67._FreqParent")) 
# FC_response2$FCActivCD8_Vax1 <- FC_response2$`Post 1st dose`/FC_response2$Baseline; FC_response2$Cohort <- NULL
# FC_response <- merge(x=FC_response, y=FC_response2, all = T, by = c("Record.ID","Alias"))
# 
# subsetData <- subset(mergedData, timeCategory == "Post 1st dose" | timeCategory == "Baseline" | timeCategory == "oneDay")
# FC_response2 <- dcast( subsetData, `Record.ID`+`Alias`~`timeCategory`, value.var = c("binding_IgG_S1")) 
# FC_response2$FC_IgG_S1_postVax1 <- FC_response2$`Post 1st dose`/FC_response2$Baseline; FC_response2$Cohort <- NULL
# FC_response <- merge(x=FC_response, y=FC_response2, all = T, by = c("Record.ID","Alias"))
# 
# subsetData <- subset(mergedData, timeCategory == "Post 2nd dose" | timeCategory == "Pre 2nd dose")
# FC_response2 <- dcast( subsetData, `Record.ID`+`Alias`~`timeCategory`, value.var = c("binding_IgG_S1")) 
# FC_response2$FC_IgG_S1_postVax2 <- FC_response2$`Post 2nd dose`/FC_response2$`Pre 2nd dose`; FC_response2$Cohort <- NULL
# FC_response <- merge(x=FC_response, y=FC_response2, all = T, by = c("Record.ID","Alias"))
# 
# subsetData <- subset(mergedData, timeCategory == "Post 2nd dose" | timeCategory == "Post 1st dose")
# FC_response2 <- dcast( subsetData, `Record.ID`+`Alias`~`timeCategory`, value.var = c("Elispot_IgG_S1")) 
# FC_response2$FC_Elispot_IgG_S1 <- FC_response2$`Post 2nd dose`/FC_response2$`Post 1st dose`; FC_response2$Cohort <- NULL
# FC_response <- merge(x=FC_response, y=FC_response2, all = T, by = c("Record.ID","Alias"))
# 
# 
# 
# 
# 
# FC_response <- FC_response[, c("Record.ID","Alias","FCtfh_Vax1","FCtfh_Vax2","FCActivCD4_Vax1","FCActivCD8_Vax1","FCPB_Vax1", "FC_IgG_S1_postVax1", 
#                                "FC_IgG_S1_postVax2", "FC_Elispot_IgG_S1")]
# mergedData <- merge(x = mergedData, y= FC_response, all=T, by = c('Record.ID', 'Alias'))
# 
# mergedData$dummy <- "dummy"

#' ------------------ Final data object --------------------------
#'

saveRDS(mergedData.IMID, file = "mergedData.IMID.Rds")
# write.csv(mergedData.IMID, file = "../Data/mergedData.IMID.csv")




subsetData <- mergedData.IMID[- which(is.na(mergedData.IMID$timeCategory) ),]
subsetData <- subsetData[-which(subsetData$Record.ID == "ID337"),]        # exclude because only day 0 and day 7 timepoints
subsetData <- subsetData[-which(subsetData$Label == "ID06_V2"),]        # exclude because the day 7 timepoint is unique
subsetData <- subset(subsetData, timeCategory == "Baseline" | timeCategory == "Post 1st dose" | timeCategory == 'Post 2nd dose')
subsetData <- subset(subsetData,  timeCategory != "two Weeks");  subsetData$timeCategory <- factor(subsetData$timeCategory, levels = c("Baseline", "Post 1st dose", "Pre 2nd dose", "Post 2nd dose"))
subsetData <- subset(subsetData, Cohort == "HC" | Cohort == "IMID")
levels(subsetData$timeCategory) <- c("Baseline","Baseline","Post 2nd dose", "Post 2nd dose")
prePostTime(subsetData, xData = "timeCategory", yData="CD4_.CD38.Ki67._FreqParent", fillParam = "Cohort", groupby="Record.ID", title = "Activated CD4", 
            xLabel = " ", yLabel = "Ki67+CD38+ (% CD4)", repMeasures = F, exponential=F)  #+ #scale_y_continuous(limits = c(0,5))
#  ggformula::geom_spline(aes_string(x="timeCategory", y="CD4_.CD38.Ki67._FreqParent"), color="blue",size=1, spar=0.5) #+ 
ggrepel::geom_text_repel(aes(label=Record.ID),size=3)
# ggsave(filename = "./Images/IMID.ActivCD4_bothCohorts_overTime.pdf")
# ggsave(filename = "./Images/IMID.ActivCD4_bothCohorts_overTime_noLabel.pdf")


subsetData = subset(mergedData, mergedData$shortForm == 'bL' | mergedData$shortForm=='oW'); subsetData = subset(subsetData, subsetData$Prior.COVID.infection.=="No")
subsetData <- subset(subsetData, Record.ID != "CV-011" & Record.ID != "CV-012" )        # extreme outliers because of absence of Ki67 stain
prePostTime(data=subsetData, xData = "shortForm", yData="CD4_.CD38.Ki67._FreqParent", fillParam = "Prior.COVID.infection.", groupby="Record.ID", title = "Activated CD4", 
            xLabel = " ", yLabel = "Ki67+CD38+ (% CD4)") + #ggrepel::geom_text_repel(aes(label=Alias)) + 
  scale_y_continuous(breaks = seq(0,20,1), limits = c(0,5)) + scale_x_discrete(labels= c("Baseline","Post\n1st dose"))
# ggsave(filename = "./Images/ActivatedCD4_vax1_healthy_paired.pdf", width=4)



prePostTime(subsetData, xData = "timeCategory", yData="CD8_..CD38.Ki67._FreqParent", fillParam = "Cohort", groupby="Record.ID", title = "Activated CD8", 
            xLabel = " ", yLabel = "Ki67+CD38+ (% CD8)", repMeasures = F, exponential=F) + scale_y_continuous(breaks=seq(0,15,2))# + # limits = c(0,6))
#  ggformula::geom_spline(aes_string(x="timeCategory", y="CD4_.CD38.Ki67._FreqParent"), color="blue",size=1, spar=0.5) #+ 
ggrepel::geom_text_repel(aes(label=Record.ID),size=3)
# ggsave(filename = "./Images/IMID.ActivCD8_bothCohorts_overTime.pdf")
# ggsave(filename = "./Images/IMID.ActivCD8_bothCohorts_overTime_noLabel.pdf")


subsetData = subset(mergedData.IMID, shortForm == 'bL' | shortForm=='4W'); subsetData = subset(subsetData, Cohort =="HC")
prePostTime(data=subsetData, xData = "shortForm", yData="CD4_.CD38.Ki67._FreqParent", fillParam = "Cohort", groupby="Record.ID", title = "Activated CD4", 
            xLabel = " ", yLabel = "Ki67+CD38+ (% CD4)") + #ggrepel::geom_text_repel(aes(label=Alias)) + 
  scale_y_continuous(breaks = seq(0,20,1), limits = c(0,5)) + scale_x_discrete(labels= c("Baseline","Post\n1st dose"))
# ggsave(filename = "./Images/IMID.ActivatedCD4_vax1_healthy_paired.pdf", width=4)



subsetData <- mergedData.IMID[- which(is.na(mergedData.IMID$timeCategory) ),]
subsetData <- subsetData[-which(subsetData$Record.ID == "ID337"),]        # exclude because only day 0 and day 7 timepoints
subsetData <- subsetData[-which(subsetData$Label == "ID06_V2"),]        # exclude because the day 7 timepoint is unique
subsetData <- subset(subsetData, timeCategory == "Baseline" | timeCategory == "Post 1st dose" | timeCategory == 'Post 2nd dose')
subsetData <- subset(subsetData,  timeCategory != "two Weeks");  subsetData$timeCategory <- factor(subsetData$timeCategory, levels = c("Baseline", "Post 1st dose", "Pre 2nd dose", "Post 2nd dose"))
subsetData <- subset(subsetData, Cohort == "HC" | Cohort == "IMID")
levels(subsetData$timeCategory) <- c("Baseline","Baseline","Post 2nd dose", "Post 2nd dose")
subsetData <- subsetData[-which(subsetData$Label == "SE01_V2" | subsetData$Record.ID == "ID793"),]
# subsetData = subset(mergedData.IMID, timeCategory == 'Baseline' | timeCategory =='Post 2nd dose'); # subsetData = subset(subsetData, Cohort =="IMID")

# pdf(file = "../../IMIDcohort/Analysis/Images/IMID.ActivCD4_paired.pdf")
prePostTime(data=subsetData, xData = "timeCategory", yData="CD4_.CD38.Ki67._FreqParent", fillParam = "Cohort", groupby="Record.ID", title = "Activated CD4", 
            xLabel = " ", yLabel = "Ki67+CD38+ (% CD4)") + #ggrepel::geom_text_repel(aes(label=Alias)) + 
  scale_y_continuous(breaks = seq(0,20,1), limits = c(0,5)) + scale_x_discrete(labels= c("Baseline","Post\n1st dose"))
# dev.off()


# pdf(file = "../../IMIDcohort/Analysis/Images/IMID.ActivCD8_paired.pdf")
prePostTime(data=subsetData, xData = "timeCategory", yData="CD8_..CD38.Ki67._FreqParent", fillParam = "Cohort", groupby="Record.ID", title = "Activated CD8", 
            xLabel = " ", yLabel = "Ki67+CD38+ (% CD8)") + #ggrepel::geom_text_repel(aes(label=Alias)) + 
  scale_y_continuous(breaks = seq(0,20,1), limits = c(0,15)) + scale_x_discrete(labels= c("Baseline","Post\n1st dose"))
# dev.off()


prePostTime(data = subsetData, xData = "timeCategory", yData="CD4_.Nonnaive.cTfh.ICOS..CD38.._FreqParent", fillParam = "Cohort", groupby="Record.ID", title = "Tfh responses", 
            xLabel = " ", yLabel = "ICOS+CD38+ (% cTfh)", repMeasures = T, exponential=F)  #+ #scale_y_continuous(limits = c(0,5))
#  ggformula::geom_spline(aes_string(x="timeCategory", y="CD4_.CD38.Ki67._FreqParent"), color="blue",size=1, spar=0.5) #+ 
# ggrepel::geom_text_repel(aes(label=Record.ID),size=3)



prePostTime(subsetData, xData = "timeCategory", yData="CD19_..Nonnaive.B.CD27..CD38._FreqParent", fillParam = "Cohort", groupby="Record.ID", title = "PB responses", 
            xLabel = " ", yLabel = "CD27+CD38+ (% Nonnaive CD19)", repMeasures = T, exponential=F)  #+ #scale_y_continuous(limits = c(0,5))
#  ggformula::geom_spline(aes_string(x="timeCategory", y="CD4_.CD38.Ki67._FreqParent"), color="blue",size=1, spar=0.5) #+ 
ggrepel::geom_text_repel(aes(label=Record.ID),size=3)



pdf(file = "../../IMIDcohort/Analysis/Images/IMID.SpikeReactive_paired.pdf")
prePostTime(subsetData, xData = "timeCategory", yData="CD19_..Nonnaive.B.Spike._FreqParent", fillParam = "Cohort", groupby="Record.ID", title = "Spike frequency", 
            xLabel = " ", yLabel = "Spike+ (% Nonnaive CD19)", repMeasures = T, exponential=F)  #+ #scale_y_continuous(limits = c(0,5))
#  ggformula::geom_spline(aes_string(x="timeCategory", y="CD4_.CD38.Ki67._FreqParent"), color="blue",size=1, spar=0.5) #+ 
# ggrepel::geom_text_repel(aes(label=Record.ID),size=3)
dev.off()
















