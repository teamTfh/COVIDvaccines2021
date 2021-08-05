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

# demog <- read.csv(file = "../Data/COVIDvaccinesObserva-CodedReport_DATA_LABELS_2021-03-16_1049.csv", colClasses = 'character')    
# demog <- read.csv(file = "../Data/COVIDvaccinesObserva-CodedReport_DATA_LABELS_2021-05-24_2112.csv", colClasses = 'character')    
demog <- read.csv(file = "../Data/COVIDvaccinesObserva-CodedReport_DATA_LABELS_2021-07-13_1744.csv", colClasses = 'character')    
index <- grep(pattern = "Vaccine.1.date", names(demog), value=F)  # search for vax1date in case its position has shifted due to new columns preceding it
if ( !is.na(strptime(demog[1,index], format = "%m/%d/%Y")) )       # don't know why but format shifts with exports, so will acount for both possibilities here
  demog[,(index-1):ncol(demog)] <- apply(X = demog[,(index-1):ncol(demog)], MARGIN = 2, FUN = function(x) { strptime(x, format = "%m/%d/%Y")} )      # convert date fields to true dates
if ( !is.na(strptime(demog[1,index], format = "%Y-%m-%d")) )
  demog[,(index-1):ncol(demog)] <- apply(X = demog[,(index-1):ncol(demog)], MARGIN = 2, FUN = function(x) { strptime(x, format = "%Y-%m-%d")} )    # convert date fields to true dates

demog <- demog[-which(demog$Notes != ""),]      # exclude anyone with caveats 
demog[which(demog$Alias == "HV-002"), "Prior.COVID.infection."] <- "No"       # this person did not have COVID during the timefrmae of this study

demog$Age <- as.numeric(demog$Age)
demog$Race <- substr(demog$Race, 4, 30)
index <- grep("Blood",names(demog),value=F)[1]  # start at first Blood.draw.date field
for(i in index:ncol(demog))    { demog[,i] <- as.numeric( round(difftime(time1 = demog[,i], time2=demog$Vaccine.1.date, units='days'), digits = 0)) }
names(demog)[grep("Blood.draw",names(demog),value=F)] <- paste0("V",seq(1,length(names(demog[grep("Blood.draw",names(demog),value=F)]))))     # convert to V1-Vx format
demog$DPO.covid <- round(as.numeric(difftime(time1 = demog$Date.of.Onset.of.Symptoms, time2 = demog$Vaccine.1.date, units = 'days')), 0); demog$Date.of.Onset.of.Symptoms <- NULL

demog.melt <- melt(demog, id.vars = c("Record.ID", "Alias"), measure.vars = paste0("V",1:10))
demog.melt$Label <- paste0(demog.melt$Alias, "_", demog.melt$variable); 
demog.melt$Visit <- demog.melt$variable;  demog.melt$variable <- NULL
names(demog.melt)[grep("value",names(demog.melt), value=F)] <- "DPV"

# make range for classifying continuous days into categorical
x <- subset(demog.melt, DPV<100 );  hist(x$DPV , breaks=100)

# Specimen collection cases:
# baseline -  from before vaccine to 1 day after
# first post-vaccine -   from 
# had two subjects who had extra draw at day 14  so there needs to be an extra assignment step to account for visits being 1 off
demog.melt$shortForm <- demog.melt$timeCategory <- ""
demog.melt[grep(pattern="1$", x = demog.melt$Label, value=F),"shortForm"] <- "bL"; demog.melt[grep(pattern="1$",x = demog.melt$Label, value=F),"timeCategory"] <- "Baseline"

# -------- post 1st dose --------------  
#demog.melt$shortForm == "" &                demog.melt$timeCategory == "" & 
demog.melt[grep(pattern="2$", x = demog.melt$Label, value=F),"shortForm"] <- "oW"; demog.melt[grep(pattern="2$", x = demog.melt$Label, value=F),"timeCategory"] <- "Post 1st dose"

demog.melt[which( demog.melt$Visit == 'V3' & demog.melt$DPV == 14),"shortForm"] <- "2W" 
                  # ^ not already assigned      ^ Label ends in '3' meaning Visit 3 only            ^ in the 13-15 day range  
  demog.melt[which(demog.melt$Visit == 'V3'  & demog.melt$DPV == 14),"timeCategory"] <- "two Weeks"

# -------- pre 2nd dose --------------  
demog.melt[which(demog.melt$Visit == 'V3' & demog.melt$DPV %in% 16:28),"shortForm"] <- "3W" 
  demog.melt[which( demog.melt$Visit == 'V3' & demog.melt$DPV %in% 16:28),"timeCategory"] <- "Pre 2nd dose"

demog.melt[which(demog.melt$Visit == 'V4' & demog.melt$DPV %in% 20:22),"shortForm"] <- "3W" 
  demog.melt[which(demog.melt$Visit == 'V4' &  demog.melt$DPV %in% 20:22),"timeCategory"] <- "Pre 2nd dose"
  
  # -------- post 2nd dose --------------  
demog.melt[which( demog.melt$Visit == 'V4' & demog.melt$DPV %in% 26:36),"shortForm"] <- "4W" 
  demog.melt[which( demog.melt$Visit == 'V4' & demog.melt$DPV %in% 26:36),"timeCategory"] <- "Post 2nd dose"
  
demog.melt[which(demog.melt$shortForm == "" & demog.melt$Visit == 'V5' & demog.melt$DPV %in% 26:36),"shortForm"] <- "4W" 
  demog.melt[which(demog.melt$timeCategory == "" & demog.melt$Visit == 'V5' & demog.melt$DPV %in% 26:36),"timeCategory"] <- "Post 2nd dose"
  
# -------- 2 wks post 2nd dose --------------  
demog.melt[which( grep(pattern="6$", x = demog.melt$Label, value=F) & demog.melt$Label == "HV-084_V6"),"shortForm"] <- "5W" 
  demog.melt[which( grep(pattern="6$", x = demog.melt$Label, value=F) & demog.melt$Label == "HV-084_V6"),"timeCategory"] <- "2 wks post 2nd dose"

demog.melt[which(demog.melt$Label== "PHI-053_V3"),"shortForm"] <- '4W'
demog.melt[which(demog.melt$Label== "PHI-053_V3"),"timeCategory"] <- 'Post 2nd dose'

# -------- One month post 2nd dose -------------- 
demog.melt[which( demog.melt$Visit == 'V4' & demog.melt$DPV %in% 40:60 ),"shortForm"] <- "oM" 
demog.melt[which( demog.melt$Visit == 'V4' & demog.melt$DPV %in% 40:60),"timeCategory"] <- "One month post\n2nd dose"

demog.melt[which( demog.melt$Visit == 'V5' & demog.melt$DPV %in% 40:60 ),"shortForm"] <- "oM" 
demog.melt[which( demog.melt$Visit == 'V5' & demog.melt$DPV %in% 40:60),"timeCategory"] <- "One month post\n2nd dose"

demog.melt[which(demog.melt$shortForm == "" & demog.melt$Visit == 'V6' & demog.melt$DPV %in% 40:60 ),"shortForm"] <- "oM" 
demog.melt[which(demog.melt$timeCategory == "" & demog.melt$Visit == 'V6' & demog.melt$DPV %in% 40:60 ),"timeCategory"] <- "One month post\n2nd dose"

demog.melt[which(demog.melt$shortForm == "" & demog.melt$Visit == 'V7' & demog.melt$DPV %in% 40:89 ),"shortForm"] <- "oM" 
demog.melt[which(demog.melt$timeCategory == "" & demog.melt$Visit == 'V7' & demog.melt$DPV %in% 40:89 ),"timeCategory"] <- "One month post\n2nd dose"

# -------- Four month post 2nd dose -------------- 
demog.melt[which( demog.melt$DPV %in% 90:200 ),"shortForm"] <- "4M" 
demog.melt[which( demog.melt$DPV %in% 90:200),"timeCategory"] <- "Four months post\n2nd dose"


#' ------------------ total binding ELISA  --------------------------
#'

totBinding <- read.csv("../Data/totalBindingELISA.csv")

#' ------------------ NeutAb   --------------------------
#'

neutAb <- read.csv("../Data/neutAb.csv")


#' ------------------ NeutAb   --------------------------
#'

ASCelispot <- read.csv("../Data/ASCelispot.csv")

#' ------------------ CXCL13 --------------------------                  
#'

CXCL13 <- read.csv("../Data/CXCL13.csv")


#' ------------------ Avidity --------------------------                  
#'

avidity <- read.csv("../Data/avidity.csv")



#' ------------------ Flow Cytometry --------------------------
#'

flowData.freq <- read.csv(file="../Data/ByParent.csv")                                 # flow cytometry ByParent spreadsheet
flowData.freq <- flowData.freq[-grep(paste(c("Mean","SD"),collapse = "|"), flowData.freq$X, value=F), ]         # delete the Mean and SD rows
names(flowData.freq)[1] <- "fcsFile"                # simple name

flowDataPHIbL.freq <- read.csv(file="../Data/ByParent_PHI Baselines.csv")                                 # flow cytometry ByParent_PHI Baselines spreadsheet
flowDataPHIbL.freq <- flowDataPHIbL.freq[-grep(paste(c("Mean","SD"),collapse = "|"), flowDataPHIbL.freq$X, value=F), ]         # delete the Mean and SD rows
names(flowDataPHIbL.freq)[1] <- "fcsFile"                # simple name

flowDataPHIbL.freq <- flowDataPHIbL.freq[-which(flowDataPHIbL.freq$fcsFile == "CV-015_PHI-048_bL_CPT.fcs"),]           # exclude CV-015_PHI-048_bL -- will use oD as baseline

flowData.freq <- rbind(flowData.freq, flowDataPHIbL.freq) #merging ByParent and ByParent_PHIbaselines 

# flowData.freq <- flowData.freq[-grep("CV-013", flowData.freq$fcsFile, value=F), ] # exclude due to active COVID at time of vaccination so uncertain cohort

a <- do.call(rbind.data.frame, strsplit(flowData.freq$fcsFile, "_"))      # split FCS file name by underscore
names(a) <- c("Record.ID","Alias","categoricalVisit","Tube","overflow")
a$Tube <- substr(a$Tube,1,3)                                  # get rid of any suffix .fcs and have tube type stand alone
flowData.freq <- as.data.frame(cbind(a,flowData.freq))
names(flowData.freq) <- str_replace(names(flowData.freq),pattern="...Freq..of.Parent....", "_FreqParent")
names(flowData.freq) <- str_replace(names(flowData.freq),pattern="Live.CD16..CD14..CD3.CD4", "CD4_")
names(flowData.freq) <- str_replace(names(flowData.freq),pattern="Live.CD16..CD14..CD3.CD8.", "CD8_")
names(flowData.freq) <- str_replace(names(flowData.freq),pattern="Live.CD16..CD14..CD19.", "CD19_")
names(flowData.freq) <- str_replace(names(flowData.freq),pattern="Live.CD16..CD14..CD3", "CD3_")
flowData.freq[grep(pattern="HV-18", x = flowData.freq$Alias, value=F),"Alias"] <- "HV-018"                 # correcting an error in fcs file name format


flowData.freq$categoricalVisit[which(flowData.freq$categoricalVisit == "oD")] <- "bL"
flowData.freq$overflow <- NULL

flowData.freq <- flowData.freq[-which(flowData.freq$fcsFile == "CV-022_PHI-021_bL_HEP_Live.fcs"),]                        # exclude CV-022_PHI-021_bL due to failed QC
# flowData.freq <- flowData.freq[-which(flowData.freq$fcsFile == "CV-015_PHI-048_bL_CPT.fcs"),]                        # exclude because using oneDay draw as baseline
if("PHI-390" %in% flowData.freq$Alias)  {flowData.freq <- flowData.freq[-which(flowData.freq$Alias == "PHI-390"),]   }    # exclude due to missing fcs file for 3w



#' ------------------ Btetramer panel flow Cytometry --------------------------
#'

Btet <- read.csv(file="../Data/ByParent_Btetramer.csv")                                 # flow cytometry Btetramer spreadsheet
Btet <- Btet[-grep(paste(c("Mean","SD"),collapse = "|"), Btet$X, value=F), ]         # delete the Mean and SD rows
Btet$X <- Btet$X.1 <- NULL                # simple name
a <- do.call(rbind.data.frame, strsplit(Btet$SampleID, "_"))      # split FCS file name by underscore
names(a) <- c("Alias","categoricalVisit")
for (i in 1:nrow(a)) 
{
  if(a[i,"categoricalVisit"] == "baseline") { a[i, "categoricalVisit"] <- "Baseline" }
  if(a[i,"categoricalVisit"] == "oneWeek") { a[i, "categoricalVisit"] <- "Post 1st dose" }
  if(a[i,"categoricalVisit"] == "Pre2nd") { a[i, "categoricalVisit"] <- "Pre 2nd dose" }
  if(a[i,"categoricalVisit"] == "Post2nd") { a[i, "categoricalVisit"] <- "Post 2nd dose" }
  if(a[i,"categoricalVisit"] == "oneMonth") { a[i, "categoricalVisit"] <- "One month post\n2nd dose" }
}

Btet <- as.data.frame(cbind(a,Btet))
names(Btet) <- str_replace(names(Btet),pattern="...Freq..of.Parent....", "_FreqParent")
names(Btet) <- str_replace(names(Btet),pattern="Lymphocytes.Singlets.Live.CD14.CD16..CD19..RBD.decoy..IgD.IgG....Geometric.Mean..RBD..1", "RBDhiIgDloIgGhi_mfiRBD")
names(Btet) <- str_replace(names(Btet),pattern="Lymphocytes.Singlets.Live.CD14.CD16..CD19..RBD.decoy..IgD.IgG....Geometric.Mean..RBD.", "RBDhiIgDhiIgGlo_mfiRBD")
names(Btet) <- str_replace(names(Btet),pattern="Lymphocytes.Singlets.Live.CD14.CD16..CD19..RBD.decoy....Geometric.Mean..RBD.", "RBDhi_mfiRBD")
names(Btet) <- str_replace(names(Btet),pattern="Lymphocytes.Singlets.Live.CD14.CD16..CD19..IgD.IgG._FreqParent.1", "CD19hiIgDloIgGhi")
names(Btet) <- str_replace(names(Btet),pattern="Lymphocytes.Singlets.Live.CD14.CD16..CD19..IgD.IgG._FreqParent", "CD19hiIgDhiIgGlo")
names(Btet) <- str_replace(names(Btet),pattern="Lymphocytes.Singlets.Live.CD14.CD16..CD19..RBD.decoy..CD27lo.MatureB.DN.DN2...Freq..of.CD19.....", "RBDhi_DN2_FreqCD19")
names(Btet) <- str_replace(names(Btet),pattern="Lymphocytes.Singlets.Live.CD14.CD16..CD19..CD27lo.MatureB.DN.DN2...Freq..of.CD19.....", "CD19hi_DN2_FreqCD19")
names(Btet) <- str_replace(names(Btet),pattern="Lymphocytes.Singlets.Live.CD14.CD16..CD19..RBD.decoy..IgD.IgG._FreqParent.1", "RBDhiIgDloIgGhi")
names(Btet) <- str_replace(names(Btet),pattern="Lymphocytes.Singlets.Live.CD14.CD16..CD19..RBD.decoy..IgD.IgG._FreqParent", "RBDhiIgDhiIgGlo")
names(Btet) <- str_replace(names(Btet),pattern="Lymphocytes.Singlets.Live.CD14.CD16..CD19..RBD.decoy..", "RBD_")
names(Btet) <- str_replace(names(Btet),pattern="Lymphocytes.Singlets.Live.CD14.CD16..CD19..", "CD19_")
names(Btet) <- str_replace(names(Btet),pattern="Live.CD16..CD14..CD3", "CD3_")
names(Btet)[4:ncol(Btet)] <- paste0("Btet_",names(Btet)[4:ncol(Btet)])



#' ------------------ AIM panel flow Cytometry --------------------------
#'

AIM <- read.csv(file="../Data/ByParent_AIM.csv")                                 # flow cytometry AIM spreadsheet
AIM <- AIM[-grep(paste(c("Mean","SD"),collapse = "|"), AIM$X, value=F), ]         # delete the Mean and SD rows
AIM <- AIM[-grep(pattern = "Un$", AIM$SampleID),]     # delete the unstim data
AIM$X <- AIM$X.1 <- NULL                # simple name
a <- do.call(rbind.data.frame, strsplit(AIM$SampleID, "_"))      # split FCS file name by underscore
a <- a[,-3]   # delete out the third column
names(a) <- c("Alias","categoricalVisit")
for (i in 1:nrow(a)) 
{
  if(a[i,"categoricalVisit"] == "baseline") { a[i, "categoricalVisit"] <- "Baseline" }
  if(a[i,"categoricalVisit"] == "oneWeek") { a[i, "categoricalVisit"] <- "Post 1st dose" }
  if(a[i,"categoricalVisit"] == "Pre2nd") { a[i, "categoricalVisit"] <- "Pre 2nd dose" }
  if(a[i,"categoricalVisit"] == "Post2nd") { a[i, "categoricalVisit"] <- "Post 2nd dose" }
  if(a[i,"categoricalVisit"] == "oneMonth") { a[i, "categoricalVisit"] <- "One month post\n2nd dose" }
}

AIM <- as.data.frame(cbind(a,AIM))
# AIM <- AIM[, -which(names(AIM) == "Condition" | names(AIM) == "SampleID")]
names(AIM) <- str_replace(names(AIM),pattern="...Freq..of.Parent....", "_FreqParent")
names(AIM) <- str_replace(names(AIM),pattern="Lymphocytes.Single.Cells.Live.CD3.", "")

names(AIM)[3:ncol(AIM)] <- paste0("AIM_",names(AIM)[3:ncol(AIM)])

#' ------------------ merge Data  --------------------------
#'

merge.1 <-  merge(x=demog, y=demog.melt, all=T, by=c("Record.ID","Alias"))
merge.serology <- merge(x = totBinding, y=neutAb, all=T, by=c("Label"))
merge.serology <- merge(x = merge.serology, y=ASCelispot, all=T,  by=c("Label"))
merge.2 <- merge(x = merge.1, y = merge.serology, all = T, by=c("Label"))
merge.3 <- merge(x = merge.2, y = avidity, all=T, by=c("Label"))
merge.4 <- merge(x = merge.3, y = CXCL13, all=T, by=c("Label"))
merge.5 <- merge(x = merge.4, y = flowData.freq, all=T, by.x=c("Alias","shortForm", "Record.ID"), by.y = c("Alias","categoricalVisit","Record.ID"))
merge.6 <- merge(x = merge.5, y = Btet, all=T, by.x = c("Alias", "timeCategory"), by.y = c("Alias", "categoricalVisit"))
merge.7 <- merge(x = merge.6, y=AIM, all=T, by.x = c("Alias", "timeCategory"), by.y = c("Alias", "categoricalVisit"))
mergedData <- merge.7

# mergedData <- merge(x = benchling, y=flowData.freq, by.x = c("Record.ID","Alias","shortForm"), by.y = c("Record.ID","Alias","categoricalVisit") ) 
# mergedData <- merge(x=demog, y=mergedData, by=c("Record.ID","Alias"))           # not elegant but it'll do

flowData.freq[-which(flowData.freq$fcsFile %in% mergedData$fcsFile), 1:5]         # what flow data did not match with final data file? 
# mergedData[-which(mergedData$fcsFile %in% flowData.freq$fcsFile), 1:10]         # what flow data did not match with final data file? 

mergedData$shortForm <- factor(mergedData$shortForm, levels = c("bL","oW","2W","3W","4W","5W","oM" ))# ,"4M"))          
mergedData$timeCategory <- factor(mergedData$timeCategory, levels = c("Baseline","Post 1st dose","two Weeks", "Pre 2nd dose","Post 2nd dose",
                                                                      "2 wks post 2nd dose", "One month post\n2nd dose" ))#, "Four months post\n2nd dose"))

subsetData <- subset(mergedData, timeCategory == "Post 1st dose" | timeCategory == "Baseline")
FC_response <- dcast(subsetData, `Record.ID`+`Alias` ~`timeCategory`, value.var = c("CD4_.Nonnaive.cTfh.ICOS..CD38.._FreqParent"))
FC_response$FCtfh_Vax1 <- FC_response$`Post 1st dose`/FC_response$Baseline
saveNames <- c("FCtfh_Vax1")

subsetData <- subset(mergedData, timeCategory == "Post 1st dose" | timeCategory == "Baseline")
FC_response2 <- dcast( subsetData, `Record.ID`+`Alias`~`timeCategory`, value.var = c("CD19_.CD27..CD38._FreqParent")) 
FC_response2$FCPB_Vax1 <- FC_response2$`Post 1st dose`/FC_response2$Baseline; FC_response2$Cohort <- NULL
FC_response <- merge(x=FC_response, y=FC_response2, by = c("Record.ID","Alias"))
saveNames <- c(saveNames, "FCPB_Vax1")

subsetData <- subset(mergedData, timeCategory == "Post 2nd dose" | timeCategory == "Pre 2nd dose")
FC_response2 <- dcast( subsetData, `Record.ID`+`Alias`~`timeCategory`, value.var = c("CD19_.CD27..CD38._FreqParent")) 
FC_response2$FCPB_Vax2 <- FC_response2$`Post 2nd dose`/FC_response2$`Pre 2nd dose`; FC_response2$Cohort <- NULL
FC_response <- merge(x=FC_response, y=FC_response2, by = c("Record.ID","Alias"))
saveNames <- c(saveNames, "FCPB_Vax2")

subsetData <- subset(mergedData, timeCategory == "Post 2nd dose" | timeCategory == "Pre 2nd dose")
FC_response2 <- dcast( subsetData, `Record.ID`+`Alias`~`timeCategory`, value.var = c("CD4_.Nonnaive.cTfh.ICOS..CD38.._FreqParent")) 
FC_response2$FCtfh_Vax2 <- FC_response2$`Post 2nd dose`/FC_response2$`Pre 2nd dose`; FC_response2$Cohort <- NULL
FC_response <- merge(x=FC_response, y=FC_response2, all = T, by = c("Record.ID","Alias"))
saveNames <- c(saveNames, "FCtfh_Vax2")

subsetData <- subset(mergedData, timeCategory == "Post 1st dose" | timeCategory == "Baseline")
FC_response2 <- dcast( subsetData, `Record.ID`+`Alias`~`timeCategory`, value.var = c("CD4_.Nonnaive.cTfh.ICOS..CD38...CXCR3._FreqParent")) 
FC_response2$FCtfh_CXCR3_Vax1 <- FC_response2$`Post 1st dose`/FC_response2$Baseline; FC_response2$Cohort <- NULL
FC_response <- merge(x=FC_response, y=FC_response2, by = c("Record.ID","Alias"))
saveNames <- c(saveNames, "FCtfh_CXCR3_Vax1")

subsetData <- subset(mergedData, timeCategory == "Post 2nd dose" | timeCategory == "Pre 2nd dose")
FC_response2 <- dcast( subsetData, `Record.ID`+`Alias`~`timeCategory`, value.var = c("CD4_.Nonnaive.cTfh.ICOS..CD38...CXCR3._FreqParent")) 
FC_response2$FCtfh_CXCR3_Vax2 <- FC_response2$`Post 2nd dose`/FC_response2$`Pre 2nd dose`; FC_response2$Cohort <- NULL
FC_response <- merge(x=FC_response, y=FC_response2, all = T, by = c("Record.ID","Alias"))
saveNames <- c(saveNames, "FCtfh_CXCR3_Vax2")

subsetData <- subset(mergedData, timeCategory == "Post 2nd dose" | timeCategory == "Baseline")
FC_response2 <- dcast( subsetData, `Record.ID`+`Alias`~`timeCategory`, value.var = c("CD4_.Nonnaive.cTfh.ICOS..CD38...CXCR3._FreqParent")) 
FC_response2$FCtfh_CXCR3_fullResp <- FC_response2$`Post 2nd dose`/FC_response2$`Baseline`; FC_response2$Cohort <- NULL
FC_response <- merge(x=FC_response, y=FC_response2, all = T, by = c("Record.ID","Alias"))
saveNames <- c(saveNames, "FCtfh_CXCR3_fullResp")



subsetData <- subset(mergedData, timeCategory == "Post 1st dose" | timeCategory == "Baseline" | timeCategory == "oneDay")
FC_response2 <- dcast( subsetData, `Record.ID`+`Alias`~`timeCategory`, value.var = c("CD4_.CD38.Ki67._FreqParent")) 
FC_response2$FCActivCD4_Vax1 <- FC_response2$`Post 1st dose`/FC_response2$Baseline; FC_response2$Cohort <- NULL
FC_response <- merge(x=FC_response, y=FC_response2, all = T, by = c("Record.ID","Alias"))
saveNames <- c(saveNames, "FCActivCD4_Vax1")

FC_response2 <- dcast( subsetData, `Record.ID`+`Alias`~`timeCategory`, value.var = c("CD8_.CD38.Ki67._FreqParent")) 
FC_response2$FCActivCD8_Vax1 <- FC_response2$`Post 1st dose`/FC_response2$Baseline; FC_response2$Cohort <- NULL
FC_response <- merge(x=FC_response, y=FC_response2, all = T, by = c("Record.ID","Alias"))
saveNames <- c(saveNames, "FCActivCD8_Vax1")

subsetData <- subset(mergedData, timeCategory == "Post 2nd dose" | timeCategory == "Pre 2nd dose")
FC_response2 <- dcast( subsetData, `Record.ID`+`Alias`~`timeCategory`, value.var = c("CD4_.CD38.Ki67._FreqParent")) 
FC_response2$FCActivCD4_Vax2 <- FC_response2$`Post 2nd dose`/FC_response2$`Pre 2nd dose`; FC_response2$Cohort <- NULL
FC_response <- merge(x=FC_response, y=FC_response2, all = T, by = c("Record.ID","Alias"))
saveNames <- c(saveNames, "FCActivCD4_Vax2")
FC_response2 <- dcast( subsetData, `Record.ID`+`Alias`~`timeCategory`, value.var = c("CD8_.CD38.Ki67._FreqParent")) 
FC_response2$FCActivCD8_Vax2 <- FC_response2$`Post 2nd dose`/FC_response2$`Pre 2nd dose`; FC_response2$Cohort <- NULL
FC_response <- merge(x=FC_response, y=FC_response2, all = T, by = c("Record.ID","Alias"))
saveNames <- c(saveNames, "FCActivCD8_Vax2")



subsetData <- subset(mergedData, timeCategory == "Post 1st dose" | timeCategory == "Baseline" | timeCategory == "oneDay")
FC_response2 <- dcast( subsetData, `Record.ID`+`Alias`~`timeCategory`, value.var = c("binding_IgG_S1")) 
FC_response2$FC_IgG_S1_postVax1 <- FC_response2$`Post 1st dose`/FC_response2$Baseline; FC_response2$Cohort <- NULL
FC_response <- merge(x=FC_response, y=FC_response2, all = T, by = c("Record.ID","Alias"))
saveNames <- c(saveNames, "FC_IgG_S1_postVax1")

subsetData <- subset(mergedData, timeCategory == "Post 2nd dose" | timeCategory == "Pre 2nd dose")
FC_response2 <- dcast( subsetData, `Record.ID`+`Alias`~`timeCategory`, value.var = c("binding_IgG_S1")) 
FC_response2$FC_IgG_S1_postVax2 <- FC_response2$`Post 2nd dose`/FC_response2$`Pre 2nd dose`; FC_response2$Cohort <- NULL
FC_response <- merge(x=FC_response, y=FC_response2, all = T, by = c("Record.ID","Alias"))
saveNames <- c(saveNames, "FC_IgG_S1_postVax2")

subsetData <- subset(mergedData, timeCategory == "Post 2nd dose" | timeCategory == "Post 1st dose")
FC_response2 <- dcast( subsetData, `Record.ID`+`Alias`~`timeCategory`, value.var = c("Elispot_IgG_S1")) 
FC_response2$FC_Elispot_IgG_S1 <- FC_response2$`Post 2nd dose`/FC_response2$`Post 1st dose`; FC_response2$Cohort <- NULL
FC_response <- merge(x=FC_response, y=FC_response2, all = T, by = c("Record.ID","Alias"))
saveNames <- c(saveNames, "FC_Elispot_IgG_S1")

# subsetData <- subset(mergedData, timeCategory == "Baseline" | timeCategory == "Post 2nd dose")
# FC_response2 <- dcast( subsetData, `Record.ID`+`Alias`~`timeCategory`, value.var = c("Elispot_IgG_S1")) 
# FC_response2$FC_Elispot_IgG_S1 <- FC_response2$`Post 2nd dose`/FC_response2$`Post 1st dose`; FC_response2$Cohort <- NULL
# FC_response <- merge(x=FC_response, y=FC_response2, all = T, by = c("Record.ID","Alias"))
# saveNames <- c(saveNames, "FC_Elispot_IgG_S1")
# 
# subsetData$CD19_.Nonnaive.B.CD27..CD38..CD138._FreqParent <- subsetData$CD19_.Nonnaive.B.CD27..CD38..CD138._FreqParent * subsetData$CD19_.Nonnaive.B.CD27..CD38._FreqParent / 100

subsetData <- subset(mergedData, timeCategory == "Baseline" | timeCategory == "Post 2nd dose")
FC_response2 <- dcast( subsetData, `Record.ID`+`Alias`~`timeCategory`, value.var = c("Btet_RBD_FreqParent")) 
FC_response2$FC_Btet_freq <- FC_response2$`Post 2nd dose`/FC_response2$`Baseline`; FC_response2$Cohort <- NULL
FC_response <- merge(x=FC_response, y=FC_response2, all = T, by = c("Record.ID","Alias"))
saveNames <- c(saveNames, "FC_Btet_freq")

subsetData <- subset(mergedData, timeCategory == "Baseline" | timeCategory == "Post 2nd dose")
FC_response2 <- dcast( subsetData, `Record.ID`+`Alias`~`timeCategory`, value.var = c("Btet_RBDhiIgDloIgGhi")) 
FC_response2$FC_Btet_IgG <- FC_response2$`Post 2nd dose`/FC_response2$`Baseline`; FC_response2$Cohort <- NULL
FC_response <- merge(x=FC_response, y=FC_response2, all = T, by = c("Record.ID","Alias"))
saveNames <- c(saveNames, "FC_Btet_IgG")

subsetData <- subset(mergedData, timeCategory == "Baseline" | timeCategory == "Post 2nd dose")
FC_response2 <- dcast( subsetData, `Record.ID`+`Alias`~`timeCategory`, value.var = c("AIM_CD4.CD69.CD200._FreqParent")) 
FC_response2$FC_AIMCD4_69200 <- FC_response2$`Post 2nd dose`/FC_response2$`Baseline`; FC_response2$Cohort <- NULL
FC_response <- merge(x=FC_response, y=FC_response2, all = T, by = c("Record.ID","Alias"))
saveNames <- c(saveNames, "FC_AIMCD4_69200")

subsetData <- subset(mergedData, timeCategory == "Baseline" | timeCategory == "Post 2nd dose")
FC_response2 <- dcast( subsetData, `Record.ID`+`Alias`~`timeCategory`, value.var = c("AIM_CD8.CD137.IFNg._FreqParent")) 
FC_response2$FC_AIMCD8_137IFNg <- FC_response2$`Post 2nd dose`/FC_response2$`Baseline`; FC_response2$Cohort <- NULL
FC_response <- merge(x=FC_response, y=FC_response2, all = T, by = c("Record.ID","Alias"))
saveNames <- c(saveNames, "FC_AIMCD8_137IFNg")


FC_response <- FC_response[, c("Record.ID","Alias",saveNames)]
mergedData <- merge(x = mergedData, y= FC_response, all=T, by = c('Record.ID', 'Alias'))




#' ------------------ Final data object --------------------------
#'
mergedData <- mergedData[- which(is.na(mergedData$timeCategory) ),]
if ( length(grep(pattern = "^X", names(mergedData))) > 0) {   mergedData <- mergedData[ , -grep(pattern = "^X", names(mergedData)) ]      }  
saveRDS(mergedData, file = "mergedData.Rds")

saveRDS(demog.melt, file = "demog.melt.Rds")


