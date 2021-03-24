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

mergedData <- readRDS(file = "mergedData.Rds")
# colors:   COVID-exp B5B2F1 (purple)             COVID-naive FFDFB1 (orange)
if("PHI-046" %in% mergedData$Alias)  {mergedData <- mergedData[-which(mergedData$Alias == "PHI-046"),] }          # only Moderna recipient, will temporarily exclude
if("HV-083" %in% mergedData$Alias)  {mergedData <- mergedData[-which(mergedData$Alias == "HV-083"),]   }          # exclude due to active COVID at time of vaccination so uncertain cohort    
if("PHI-390" %in% mergedData$Alias)  {mergedData <- mergedData[-which(mergedData$Alias == "PHI-390"),]   }        # exclude due to missing fcs file for 3w


#' ------------------ Cohort description --------------------------
#'

temp <- mergedData %>% group_by(Prior.COVID.infection., timeCategory) %>% get_summary_stats(type = 'common')
# write.csv(temp, file = "summaryStatistics.csv")

##  *****   REMEMBER EXCLUSIONS ABOVE !!!!    *******
##  *****   also deleted CV-034_PHI-058_bL_CPT because it had no CD4 stain    **********

metaData <- mergedData[,c("fcsFile","Alias","Label","timeCategory","Sex","Prior.COVID.infection.","DPO.covid","DPV")]
metaData <- metaData[-which(metaData$Label == "PHI-021_V1"),]
# fileList <- data.frame(original = list.files(path = "../../COVIDvax_export_CD16loCD14lo/"))   # all the fcs files that will go to OMIQ.ai
fileList[! fileList %in% metaData$fcsFile]
metaData$fcsFile[! metaData$fcsFile %in% fileList]

fileList <-  list.files(path = "../../COVIDvax_export_CD16loCD14lo/")
metaData <- metaData[which(metaData$fcsFile %in% fileList),]
metaData <- metaData[,c("fcsFile","Alias","Label","timeCategory","Sex","Prior.COVID.infection.","DPO.covid","DPV")]
omiqID <- read.csv(file = "../../COVIDvax_export_CD16loCD14lo/20210316_metadata.csv")
omiqID$Filename[! omiqID$Filename %in% metaData$fcsFile]

omiqID <- merge(x = omiqID, y = metaData, by.x = "Filename", by.y = "fcsFile")
omiqID <- omiqID[, c(2,1, 3:ncol(omiqID))]
# write.csv(omiqID, file = "../Flow cytometry/COVIDvax_metaData_OMIQ.ai.csv", row.names = F)

#' ------------------ Antibody analyses --------------------------
#'

# subsetData <- subset(mergedData, timeCategory != "two Weeks")
linePlot(data = mergedData, xData = 'timeCategory', yData = 'binding_IgG_S1', groupby = 'Alias', xLabel = ' ', yLabel = "anti-S1 IgG titer", 
         title = "anti-S1 IgG titer", colorby = "Prior.COVID.infection.") + 
  scale_color_manual(name="Prior COVID?",values = c("#FFDFB1","#B5B2F1")) + 
  scale_y_continuous(trans='pseudo_log', limits = c(0,1e7), breaks=c(10^(0:7)), labels=trans_format('log10',math_format(10^.x)), minor_breaks =5*10^(0:10)) #+ 
# ggrepel::geom_text_repel(data = subset(mergedData, Alias =="HV-002"), aes(label = Label))
# ggsave(filename = "./Images/BindingAb_S1_IgG_linePlot.pdf")

twoSampleBar(data = subset(mergedData, timeCategory == "Post 2nd dose"), xData = "Prior.COVID.infection.", yData = "binding_IgG_S1", fillParam = "Prior.COVID.infection.", 
             title = "Post 2nd dose", yLabel = "anti-S1 IgG titer", nonparam = T) +   coord_cartesian(ylim=c(0e1,1e7),) + 
  scale_y_continuous(trans='pseudo_log', breaks=c(10^(0:7)), labels=trans_format('log10',math_format(10^.x)), minor_breaks =5*10^(0:10)) #+ ggrepel::geom_text_repel(aes(label = Label))
# ggsave(filename = "./Images/BindingAb_S1_IgG_post2ndDose.pdf", width=5)

subsetData <- mergedData
subsetData <- subsetData[-which(subsetData$timeCategory == "two Weeks"),] 
subsetData <- subsetData[-which(subsetData$timeCategory == "One month post\n2nd dose"),] 
subsetData$timeCategory <- factor(subsetData$timeCategory, levels = c("Baseline", "Post 1st dose", "Pre 2nd dose", "Post 2nd dose"))
linePlot(data = subsetData, xData = 'timeCategory', yData = 'binding_IgA_S1', groupby = 'Alias', xLabel = ' ', yLabel = "anti-S1 IgA titer", 
         title = "anti-S1 IgA titer", colorby = "Prior.COVID.infection.") + 
  scale_color_manual(name="Prior COVID?",values = c("#FFDFB1","#B5B2F1")) + 
  scale_y_continuous(trans='pseudo_log', limits = c(0,1e5), breaks=c(10^(0:7)), labels=trans_format('log10',math_format(10^.x)), minor_breaks =5*10^(0:10)) #+ 
  # ggrepel::geom_text_repel(data = subset(mergedData, Alias =="HV-002"), aes(label = Label))
# ggsave(filename = "./Images/BindingAb_S1_IgA_linePlot.pdf")

linePlot(data = mergedData, xData = 'timeCategory', yData = 'IC50_neutAb_log10', groupby = 'Alias', xLabel = ' ', yLabel = "log10 IC50", 
         title = "Neutralizing antibodies", colorby = "Prior.COVID.infection.") + 
  scale_color_manual(name="Prior COVID?",values = c("#FFDFB1","#B5B2F1")) + 
  scale_y_continuous(trans='pseudo_log', limits = c(0,1e5), breaks=c(10^(0:6)), labels=trans_format('log10',math_format(10^.x)), minor_breaks =5*10^(0:10))# + ggrepel::geom_text_repel(aes(label = Label))
# ggsave(filename = "./Images/neutAb_linePlot.pdf")

twoSampleBar(data = subset(mergedData, timeCategory == "Post 2nd dose"), xData = "Prior.COVID.infection.", yData = "IC50_neutAb_log10", fillParam = "Prior.COVID.infection.", 
             title = "Post 2nd dose", yLabel = "Neutralizing antibodies", nonparam = T) +   coord_cartesian(ylim=c(0e1,1e5),) + 
  scale_y_continuous(trans='pseudo_log', breaks=c(10^(0:7)), labels=trans_format('log10',math_format(10^.x)), minor_breaks =5*10^(0:10)) #+ ggrepel::geom_text_repel(aes(label = Label))
# ggsave(filename = "./Images/neutAb_post2ndDose.pdf", width=5)

linePlot(data = mergedData, xData = 'timeCategory', yData = 'binding_IgG_N', groupby = 'Alias', xLabel = ' ', yLabel = "anti-N IgG titer", 
         title = "anti-N IgG titer", colorby = "Prior.COVID.infection.") + 
  scale_color_manual(name="Prior COVID?",values = c("#FFDFB1","#B5B2F1")) + #theme(legend.position = 'none') + 
  scale_y_continuous(trans='pseudo_log', limits = c(0,1e5), breaks=c(10^(0:6)), labels=trans_format('log10',math_format(10^.x)), minor_breaks =5*10^(0:10)) #+ ggrepel::geom_text_repel(aes(label = Label))
# ggsave(filename = "./Images/BindingAb_N_IgG_linePlot.pdf", width=7)


linePlot(data = mergedData, xData = 'DPV', yData = 'binding_IgG_S1', groupby = 'Alias', xLabel = 'Days after first vaccine dose', yLabel = "anti-S1 IgG titer", 
         title = "anti-S1 IgG titer", colorby = "Prior.COVID.infection.") + 
  scale_color_manual(name="Prior COVID?",values = c("#FFDFB1","#B5B2F1")) + 
  scale_y_continuous(trans='pseudo_log', limits = c(0,1e7), breaks=c(10^(0:7)), labels=trans_format('log10',math_format(10^.x)), minor_breaks =5*10^(0:10))# + ggrepel::geom_text_repel(aes(label = Label))
# ggsave(filename = "./Images/BindingAb_S1_IgG_linePlot_contTime.pdf", width=7)




subsetData <- subset(mergedData, Prior.COVID.infection. == 'Yes' & timeCategory == "Baseline")
subsetData <- subsetData[,c(24:25, 73:79)]
cor.matrix <- cor(subsetData, method="kendall" , use="pairwise.complete.obs" )
cor.matrix.pmat <- ggcorrplot::cor_pmat(subsetData, method="kendall", use="pairwise.complete.obs"  )
ggcorrplot::ggcorrplot(corr = cor.matrix, p.mat = cor.matrix.pmat, title = "Baseline Ab correlations", legend.title = "Kendall tau", insig = "blank")

univScatter(data = subset(mergedData, Prior.COVID.infection. == "Yes" & timeCategory == "Baseline"), xData = "binding_IgG_S1", yData = 'FC_IgG_S1_postVax1', 
            fillParam = 'Prior.COVID.infection.',title = "Ab response", xLabel = "Baseline anti-S1 IgG titer", yLabel = "Fold-change S1 postVax1", nonparam = T) + 
  scale_fill_manual(values=c("#B5B2F1")) + 
  scale_x_continuous(trans='pseudo_log', limits = c(1e2,1e6), breaks=c(10^(0:7)), labels=trans_format('log10',math_format(10^.x)), minor_breaks =5*10^(0:10)) +
  scale_y_continuous(trans='pseudo_log', limits = c(0,1e5), breaks=c(10^(0:7)), labels=trans_format('log10',math_format(10^.x)), minor_breaks =5*10^(0:10))
# ggsave(filename = "./Images/BindingAb_S1_IgG_vs_FC_S1binding.pdf")



# 
# bivScatter(data1 = subset(mergedData, Prior.COVID.infection. == 'No' & timeCategory == "Post 2nd dose"), data2 = subset(mergedData, Prior.COVID.infection. == "Yes" & timeCategory == "Post 2nd dose"), 
#                        name1 = "Naive", name2 = "Experienced", xData = "binding_IgG_S1", yData = 'IC50_neutAb_log10', fillParam = 'Prior.COVID.infection.', 
#                        title = "Post 2nd dose", xLabel = "anti-S1 IgG titer", yLabel = "NeutAb IC50", statsOff = F)
  

subsetData <- subset(mergedData, timeCategory == 'Baseline'); subsetData <- subsetData[which(!is.na(subsetData$binding_IgG_S1)),]
subsetData %>% levene_test( binding_IgG_S1 ~ Prior.COVID.infection.)
subsetData %>% shapiro_test( binding_IgG_S1)
subsetData %>%  wilcox_test(binding_IgG_S1 ~ Prior.COVID.infection.)
subsetData %>% group_by(Prior.COVID.infection.) %>% get_summary_stats(binding_IgG_S1)
subsetData %>% levene_test( IC50_neutAb_log10 ~ Prior.COVID.infection.)
subsetData %>% shapiro_test( IC50_neutAb_log10)
subsetData %>%  wilcox_test(IC50_neutAb_log10 ~ Prior.COVID.infection.)

subsetData <- subset(mergedData, timeCategory == 'Post 1st dose'); subsetData <- subsetData[which(!is.na(subsetData$binding_IgG_S1)),]
subsetData %>% group_by(Prior.COVID.infection.) %>% get_summary_stats(IC50_neutAb_log10)
subsetData %>% levene_test( IC50_neutAb_log10 ~ Prior.COVID.infection.)
subsetData %>% shapiro_test( IC50_neutAb_log10)
subsetData %>%  wilcox_test(IC50_neutAb_log10 ~ Prior.COVID.infection.)


subsetData %>% group_by(Prior.COVID.infection.) %>% get_summary_stats(FC_IgG_S1_postVax1)
subsetData %>% shapiro_test(FC_IgG_S1_postVax1)
subsetData %>%  wilcox_test(FC_IgG_S1_postVax1 ~ Prior.COVID.infection.)
subsetData %>% group_by(Prior.COVID.infection.) %>% get_summary_stats(FC_IgG_S1_postVax2)
subsetData %>% shapiro_test(FC_IgG_S1_postVax2)
subsetData %>%  wilcox_test(FC_IgG_S1_postVax2 ~ Prior.COVID.infection.)


subsetData <- subset(mergedData, timeCategory == 'Post 2nd dose'); subsetData <- subsetData[which(!is.na(subsetData$binding_IgG_S1)),]
subsetData %>% levene_test( binding_IgG_S1 ~ Prior.COVID.infection.)
subsetData %>% shapiro_test( binding_IgG_S1)
subsetData %>%  t_test(binding_IgG_S1 ~ Prior.COVID.infection.)

subsetData <- mergedData[which(!is.na(mergedData$FC_Elispot_IgG_S1)),]
subsetData <- subsetData[which(is.finite(subsetData$FC_Elispot_IgG_S1)),]
subsetData %>% group_by(Prior.COVID.infection.) %>% get_summary_stats(FC_Elispot_IgG_S1)


#' ------------------ ELISpot analyses --------------------------
#'
subsetData <- subset(mergedData, timeCategory == 'Post 1st dose')
a.1 <- twoSampleBar(data = subsetData, xData = 'Prior.COVID.infection.',yData='Elispot_IgG_S1', fillParam = 'Prior.COVID.infection.',title = "ASC - IgG S1", 
             yLabel = "Spots per 1e6 PBMC", nonparam=T) + 
  scale_y_continuous(trans='pseudo_log', limits = c(0,10000), breaks=c(10^(0:4)), labels=trans_format('log10',math_format(10^.x)), minor_breaks =5*10^(0:10))
a.2 <- twoSampleBar(data = subsetData, xData = 'Prior.COVID.infection.',yData='Elispot_IgG_S2', fillParam = 'Prior.COVID.infection.',title = "ASC - IgG S2", 
             yLabel = " ", nonparam=T)+ 
  scale_y_continuous(trans='pseudo_log', limits = c(0,10000), breaks=c(10^(0:4)), labels=trans_format('log10',math_format(10^.x)), minor_breaks =5*10^(0:10))
a.3 <- twoSampleBar(data = subsetData, xData = 'Prior.COVID.infection.',yData='Elispot_IgG_RBD', fillParam = 'Prior.COVID.infection.',title = "ASC - IgG RBD", 
             yLabel = " ", nonparam=T)+ 
  scale_y_continuous(trans='pseudo_log', limits = c(0,10000), breaks=c(10^(0:4)), labels=trans_format('log10',math_format(10^.x)), minor_breaks =5*10^(0:10))
grid.arrange(a.1,a.2,a.3, nrow=1)
# ggpubr::ggexport(plotlist = list(a.1, a.2, a.3), filename = "./Images/Elispots_IgG_post1stDose_gridarrange.pdf", nrow = 1, width = 12)

a.1 <- twoSampleBar(data = subsetData, xData = 'Prior.COVID.infection.',yData='Elispot_IgA_S1', fillParam = 'Prior.COVID.infection.',title = "ASC - IgA S1", 
             yLabel = "Spots per 1e6 PBMC", nonparam=T) + 
  scale_y_continuous(trans='pseudo_log', limits = c(0,10000), breaks=c(10^(0:4)), labels=trans_format('log10',math_format(10^.x)), minor_breaks =5*10^(0:10))
a.2 <- twoSampleBar(data = subsetData, xData = 'Prior.COVID.infection.',yData='Elispot_IgA_S2', fillParam = 'Prior.COVID.infection.',title = "ASC - IgA S2", 
             yLabel = " ", nonparam=T)+ 
  scale_y_continuous(trans='pseudo_log', limits = c(0,10000), breaks=c(10^(0:4)), labels=trans_format('log10',math_format(10^.x)), minor_breaks =5*10^(0:10))
a.3 <- twoSampleBar(data = subsetData, xData = 'Prior.COVID.infection.',yData='Elispot_IgA_RBD', fillParam = 'Prior.COVID.infection.',title = "ASC - IgA RBD", 
             yLabel = " ", nonparam=T)+ 
  scale_y_continuous(trans='pseudo_log', limits = c(0,10000), breaks=c(10^(0:4)), labels=trans_format('log10',math_format(10^.x)), minor_breaks =5*10^(0:10))
grid.arrange(a.1,a.2,a.3, nrow=1)
# ggpubr::ggexport(plotlist = list(a.1, a.2, a.3), filename = "./Images/Elispots_IgA_post1stDose_gridarrange.pdf", nrow = 1, width = 12 )

a.1 <- twoSampleBar(data = subsetData, xData = 'Prior.COVID.infection.',yData='Elispot_IgM_S1', fillParam = 'Prior.COVID.infection.',title = "ASC - IgM S1", 
             yLabel = "Spots per 1e6 PBMC", nonparam=T)+ 
  scale_y_continuous(trans='pseudo_log', limits = c(0,10000), breaks=c(10^(0:4)), labels=trans_format('log10',math_format(10^.x)), minor_breaks =5*10^(0:10))
a.2 <- twoSampleBar(data = subsetData, xData = 'Prior.COVID.infection.',yData='Elispot_IgM_S2', fillParam = 'Prior.COVID.infection.',title = "ASC - IgM S2", 
             yLabel = " ", nonparam=T)+ 
  scale_y_continuous(trans='pseudo_log', limits = c(0,10000), breaks=c(10^(0:4)), labels=trans_format('log10',math_format(10^.x)), minor_breaks =5*10^(0:10))
a.3 <- twoSampleBar(data = subsetData, xData = 'Prior.COVID.infection.',yData='Elispot_IgM_RBD', fillParam = 'Prior.COVID.infection.',title = "ASC - IgM RBD", 
             yLabel = " ", nonparam=T)+ 
  scale_y_continuous(trans='pseudo_log', limits = c(0,10000), breaks=c(10^(0:4)), labels=trans_format('log10',math_format(10^.x)), minor_breaks =5*10^(0:10))
grid.arrange(a.1,a.2,a.3, nrow=1)
# ggpubr::ggexport(plotlist = list(a.1, a.2, a.3), filename = "./Images/Elispots_IgM_post1stDose_gridarrange.pdf", nrow = 1, width = 12 )

subsetData <- subset(mergedData, timeCategory == 'Post 2nd dose')
a.1 <- twoSampleBar(data = subsetData, xData = 'Prior.COVID.infection.',yData='Elispot_IgG_S1', fillParam = 'Prior.COVID.infection.',title = "ASC - IgG S1", 
                    yLabel = "Spots per 1e6 PBMC", nonparam=T) + 
  scale_y_continuous(trans='pseudo_log', limits = c(0,10000), breaks=c(10^(0:4)), labels=trans_format('log10',math_format(10^.x)), minor_breaks =5*10^(0:10))
a.2 <- twoSampleBar(data = subsetData, xData = 'Prior.COVID.infection.',yData='Elispot_IgG_S2', fillParam = 'Prior.COVID.infection.',title = "ASC - IgG S2", 
                    yLabel = " ", nonparam=T)+ 
  scale_y_continuous(trans='pseudo_log', limits = c(0,10000), breaks=c(10^(0:4)), labels=trans_format('log10',math_format(10^.x)), minor_breaks =5*10^(0:10))
a.3 <- twoSampleBar(data = subsetData, xData = 'Prior.COVID.infection.',yData='Elispot_IgG_RBD', fillParam = 'Prior.COVID.infection.',title = "ASC - IgG RBD", 
                    yLabel = " ", nonparam=T)+ 
  scale_y_continuous(trans='pseudo_log', limits = c(0,10000), breaks=c(10^(0:4)), labels=trans_format('log10',math_format(10^.x)), minor_breaks =5*10^(0:10))
grid.arrange(a.1,a.2,a.3, nrow=1)
# ggpubr::ggexport(plotlist = list(a.1, a.2, a.3), filename = "./Images/Elispots_IgG_post2ndDose_gridarrange.pdf", nrow = 1, width = 12)

a.1 <- twoSampleBar(data = subsetData, xData = 'Prior.COVID.infection.',yData='Elispot_IgA_S1', fillParam = 'Prior.COVID.infection.',title = "ASC - IgA S1", 
                    yLabel = "Spots per 1e6 PBMC", nonparam=T) + 
  scale_y_continuous(trans='pseudo_log', limits = c(0,10000), breaks=c(10^(0:4)), labels=trans_format('log10',math_format(10^.x)), minor_breaks =5*10^(0:10))
a.2 <- twoSampleBar(data = subsetData, xData = 'Prior.COVID.infection.',yData='Elispot_IgA_S2', fillParam = 'Prior.COVID.infection.',title = "ASC - IgA S2", 
                    yLabel = " ", nonparam=T)+ 
  scale_y_continuous(trans='pseudo_log', limits = c(0,10000), breaks=c(10^(0:4)), labels=trans_format('log10',math_format(10^.x)), minor_breaks =5*10^(0:10))
a.3 <- twoSampleBar(data = subsetData, xData = 'Prior.COVID.infection.',yData='Elispot_IgA_RBD', fillParam = 'Prior.COVID.infection.',title = "ASC - IgA RBD", 
                    yLabel = " ", nonparam=T)+ 
  scale_y_continuous(trans='pseudo_log', limits = c(0,10000), breaks=c(10^(0:4)), labels=trans_format('log10',math_format(10^.x)), minor_breaks =5*10^(0:10))
grid.arrange(a.1,a.2,a.3, nrow=1)
# ggpubr::ggexport(plotlist = list(a.1, a.2, a.3), filename = "./Images/Elispots_IgA_post2ndDose_gridarrange.pdf", nrow = 1, width = 12 )

a.1 <- twoSampleBar(data = subsetData, xData = 'Prior.COVID.infection.',yData='Elispot_IgM_S1', fillParam = 'Prior.COVID.infection.',title = "ASC - IgM S1", 
                    yLabel = "Spots per 1e6 PBMC", nonparam=T)+ 
  scale_y_continuous(trans='pseudo_log', limits = c(0,10000), breaks=c(10^(0:4)), labels=trans_format('log10',math_format(10^.x)), minor_breaks =5*10^(0:10))
a.2 <- twoSampleBar(data = subsetData, xData = 'Prior.COVID.infection.',yData='Elispot_IgM_S2', fillParam = 'Prior.COVID.infection.',title = "ASC - IgM S2", 
                    yLabel = " ", nonparam=T)+ 
  scale_y_continuous(trans='pseudo_log', limits = c(0,10000), breaks=c(10^(0:4)), labels=trans_format('log10',math_format(10^.x)), minor_breaks =5*10^(0:10))
a.3 <- twoSampleBar(data = subsetData, xData = 'Prior.COVID.infection.',yData='Elispot_IgM_RBD', fillParam = 'Prior.COVID.infection.',title = "ASC - IgM RBD", 
                    yLabel = " ", nonparam=T)+ 
  scale_y_continuous(trans='pseudo_log', limits = c(0,10000), breaks=c(10^(0:4)), labels=trans_format('log10',math_format(10^.x)), minor_breaks =5*10^(0:10))
grid.arrange(a.1,a.2,a.3, nrow=1)
# ggpubr::ggexport(plotlist = list(a.1, a.2, a.3), filename = "./Images/Elispots_IgM_post2ndDose_gridarrange.pdf", nrow = 1, width = 12 )

subsetData <- subset(mergedData, !is.na(Elispot_IgG_S1) & timeCategory != "Pre 2nd dose" )
# pdf(file = "./Images/Elispots_IgG_S1_prepostTime.pdf")
  prePostTime(data = subsetData, xData = "timeCategory", yData = "Elispot_IgG_S1", fillParam = "Prior.COVID.infection.", 
              groupby = "Alias", title = "ASC IgG S1", xLabel = " ", yLabel = "Spots per 1e6 PBMC", exponential=T)+ 
    scale_y_continuous(trans='pseudo_log', limits = c(0,10000), breaks=c(10^(0:4)), labels=trans_format('log10',math_format(10^.x)), minor_breaks =5*10^(0:10)) 
# dev.off()

subsetData <- subset(mergedData, !is.na(Elispot_IgG_S2) & timeCategory != "Pre 2nd dose" )
# pdf(file = "./Images/Elispots_IgG_S2_prepostTime.pdf")
prePostTime(data = subsetData, xData = "timeCategory", yData = "Elispot_IgG_S2", fillParam = "Prior.COVID.infection.", 
                  groupby = "Alias", title = "ASC IgG S2", xLabel = " ", yLabel = "Spots per 1e6 PBMC", exponential=T)+ 
  scale_y_continuous(trans='pseudo_log', limits = c(0,10000), breaks=c(10^(0:4)), labels=trans_format('log10',math_format(10^.x)), minor_breaks =5*10^(0:10))
# dev.off()

subsetData <- subset(mergedData, !is.na(Elispot_IgG_RBD) & timeCategory != "Pre 2nd dose" )
# pdf(file = "./Images/Elispots_IgG_RBD_prepostTime.pdf")
prePostTime(data = subsetData, xData = "timeCategory", yData = "Elispot_IgG_RBD", fillParam = "Prior.COVID.infection.", 
                  groupby = "Alias", title = "ASC IgG RBD", xLabel = " ", yLabel = "Spots per 1e6 PBMC", exponential=T)+ 
  scale_y_continuous(trans='pseudo_log', limits = c(0,10000), breaks=c(10^(0:4)), labels=trans_format('log10',math_format(10^.x)), minor_breaks =5*10^(0:10))
# dev.off()



bivScatter(data1 = subset(mergedData, Prior.COVID.infection. == 'No' & timeCategory == "Post 2nd dose"), 
           data2 = subset(mergedData, Prior.COVID.infection. == "Yes" & timeCategory == "Post 2nd dose"),
           name1 = "Naive", name2 = "Experienced", xData = "Elispot_IgG_RBD", yData = 'Elispot_IgG_S1', fillParam = 'Prior.COVID.infection.',
           title = "ELISpots post 2nd dose", xLabel = "ASC IgG RBD", yLabel = "ASC IgG S1", statsOff = F) + 
  scale_y_continuous(trans='pseudo_log', limits = c(0,1e4), breaks=c(10^(0:4)), labels=trans_format('log10',math_format(10^.x)), minor_breaks =5*10^(0:10)) + 
  scale_x_continuous(trans='pseudo_log', limits = c(0,1e4), breaks=c(10^(0:4)), labels=trans_format('log10',math_format(10^.x)), minor_breaks =5*10^(0:10)) 
ggsave(filename = "./Images/Elispots_IgG_S1-vs-RBD_correl_biv.pdf", width=8)

bivScatter(data1 = subset(mergedData, Prior.COVID.infection. == 'No' & timeCategory == "Post 2nd dose"), 
           data2 = subset(mergedData, Prior.COVID.infection. == "Yes" & timeCategory == "Post 2nd dose"),
           name1 = "Naive", name2 = "Experienced", xData = "Elispot_IgG_RBD", yData = 'Elispot_IgG_S2', fillParam = 'Prior.COVID.infection.',
           title = "ELISpots post 2nd dose", xLabel = "ASC IgG RBD", yLabel = "ASC IgG S2", statsOff = F) +
  scale_y_continuous(trans='pseudo_log', limits = c(0,1e3), breaks=c(10^(0:4)), labels=trans_format('log10',math_format(10^.x)), minor_breaks =5*10^(0:10)) + 
  scale_x_continuous(trans='pseudo_log', limits = c(0,1e3), breaks=c(10^(0:4)), labels=trans_format('log10',math_format(10^.x)), minor_breaks =5*10^(0:10)) 
ggsave(filename = "./Images/Elispots_IgG_S2-vs-RBD_correl_biv.pdf", width = 8 )

subsetData <- subset(mergedData, Prior.COVID.infection. == 'No' & timeCategory == "Post 1st dose")
subsetData <- subsetData[,grep("^Elispot", names(subsetData))]
subsetData <- subsetData[,-grep("Elispot_IgA_S1", names(subsetData))]
cor.elispots <- cor(subsetData, method="kendall" , use="pairwise.complete.obs" )
cor.elispots.pmat <- ggcorrplot::cor_pmat(subsetData, method="kendall", use="pairwise.complete.obs"  )
ggcorrplot::ggcorrplot(corr = cor.elispots, p.mat = cor.elispots.pmat, title = "ELISpots post 1st dose", legend.title = "Kendall tau", pch.cex = 1, insig = "blank")
# ggsave(filename = "./Images/Elispots_ggcorrplot_Naive_post1st.pdf")
ggcorrplot::ggcorrplot(corr = cor.elispots, p.mat = cor.elispots.pmat, title = "ELISpots post 1st dose", legend.title = "Kendall tau", pch.cex = 1) #, insig = "blank"
ggsave(filename = "./Images/Elispots_ggcorrplot_Naive_post1st_full.pdf")

subsetData <- subset(mergedData, Prior.COVID.infection. == 'Yes' & timeCategory == "Post 1st dose")
subsetData <- subsetData[,grep("^Elispot", names(subsetData))]
cor.elispots <- cor(subsetData, method="kendall" , use="pairwise.complete.obs" )
cor.elispots.pmat <- ggcorrplot::cor_pmat(subsetData, method="kendall", use= "pairwise.complete.obs"  )
ggcorrplot::ggcorrplot(corr = cor.elispots, p.mat = cor.elispots.pmat, title = "ELISpots post 1st dose", legend.title = "Kendall tau", pch.cex = 1, insig = "blank") #
# ggsave(filename = "./Images/Elispots_ggcorrplot_Experienced_post1st.pdf")
ggcorrplot::ggcorrplot(corr = cor.elispots, p.mat = cor.elispots.pmat, title = "ELISpots post 1st dose", legend.title = "Kendall tau", pch.cex = 1) #, insig = "blank"
ggsave(filename = "./Images/Elispots_ggcorrplot_Experienced_post1st_full.pdf")




subsetData <- subset(mergedData, Prior.COVID.infection. == 'No' & timeCategory == "Post 2nd dose")
subsetData <- subsetData[,grep("^Elispot", names(subsetData))]
cor.elispots <- cor(subsetData, method="kendall" , use="pairwise.complete.obs" )
cor.elispots.pmat <- ggcorrplot::cor_pmat(subsetData, method="kendall", use="pairwise.complete.obs"  )
ggcorrplot::ggcorrplot(corr = cor.elispots, p.mat = cor.elispots.pmat, title = "ELISpots post 2nd dose", legend.title = "Kendall tau", pch.cex = 1, insig = "blank")
# ggsave(filename = "./Images/Elispots_ggcorrplot_Naive_post2nd.pdf")
ggcorrplot::ggcorrplot(corr = cor.elispots, p.mat = cor.elispots.pmat, title = "ELISpots post 2nd dose", legend.title = "Kendall tau", pch.cex = 1) #, insig = "blank"
# ggsave(filename = "./Images/Elispots_ggcorrplot_Naive_post2nd_full.pdf")

subsetData <- subset(mergedData, Prior.COVID.infection. == 'Yes' & timeCategory == "Post 2nd dose")
subsetData <- subsetData[,grep("^Elispot", names(subsetData))]
cor.elispots <- cor(subsetData, method="kendall" , use="pairwise.complete.obs" )
cor.elispots.pmat <- ggcorrplot::cor_pmat(subsetData, method="kendall", use= "pairwise.complete.obs"  )
ggcorrplot::ggcorrplot(corr = cor.elispots, p.mat = cor.elispots.pmat, title = "ELISpots post 2nd dose", legend.title = "Kendall tau", pch.cex = 1, insig = "blank") #
# ggsave(filename = "./Images/Elispots_ggcorrplot_Experienced_post2nd.pdf")
ggcorrplot::ggcorrplot(corr = cor.elispots, p.mat = cor.elispots.pmat, title = "ELISpots post 2nd dose", legend.title = "Kendall tau", pch.cex = 1) #, insig = "blank"
# ggsave(filename = "./Images/Elispots_ggcorrplot_Experienced_post2nd_full.pdf")




#' ------------------ Activated T cells analyses --------------------------
#'
subsetData <- subset(mergedData,  timeCategory != "two Weeks"& timeCategory != "2 wks post 2nd dose");
subsetData <- subsetData[which(subsetData$Tube == "HEP"),] 
subsetData <- subset(subsetData, Record.ID != "CV-011" & Record.ID != "CV-012" & Record.ID != "CV-005")        # absence of Ki67 stain
subsetData$timeCategory <- factor(subsetData$timeCategory, levels = c("Baseline", "Post 1st dose", "Pre 2nd dose", "Post 2nd dose", "One month post\n2nd dose"))
prePostTime(subsetData, xData = "DPV", yData="CD4_.CD38.Ki67._FreqParent", fillParam = "Prior.COVID.infection.", groupby="Record.ID", title = "CD4 - post dose 1",
            xLabel = "Days", yLabel = "Ki67+CD38+ (% CD4)", repMeasures = F, exponential=F) +  coord_cartesian(xlim = c(-1,12)) + geom_vline(xintercept = 0,linetype="dashed", alpha=0.5)
# ggsave(filename = "./Images/ActivCD4_bothCohorts_Vax1_continuousTime.pdf")
prePostTime(subsetData, xData = "DPV", yData="CD8_.CD38.Ki67._FreqParent", fillParam = "Prior.COVID.infection.", groupby="Record.ID", title = "CD8 - post dose 1",
            xLabel = "Days", yLabel = "Ki67+CD38+ (% CD8)", repMeasures = F, exponential=F) +  coord_cartesian(xlim = c(-1,12)) + geom_vline(xintercept = 0,linetype="dashed" , alpha=0.5)
# ggsave(filename = "./Images/ActivCD8_bothCohorts_Vax1_continuousTime.pdf")
subsetData$DPV <- subsetData$DPV - as.numeric(difftime(subsetData$Vaccine.2.date, subsetData$Vaccine.1.date, units="days" ) )
prePostTime(subsetData, xData = "DPV", yData="CD4_.CD38.Ki67._FreqParent", fillParam = "Prior.COVID.infection.", groupby="Record.ID", title = "CD4 - post dose 2",
            xLabel = "Days", yLabel = "Ki67+CD38+ (% CD4)", repMeasures = F, exponential=F) +  coord_cartesian(xlim = c(-5,12)) + geom_vline(xintercept = 0,linetype="dashed", alpha=0.5)
# ggsave(filename = "./Images/ActivCD4_bothCohorts_Vax2_continuousTime.pdf")
prePostTime(subsetData, xData = "DPV", yData="CD8_.CD38.Ki67._FreqParent", fillParam = "Prior.COVID.infection.", groupby="Record.ID", title = "CD8 - Post dose 2",
            xLabel = "Days", yLabel = "Ki67+CD38+ (% CD8)", repMeasures = F, exponential=F) +  coord_cartesian(xlim = c(-5,12))  + geom_vline(xintercept = 0,linetype="dashed" , alpha=0.5)
# ggsave(filename = "./Images/ActivCD8_bothCohorts_Vax2_continuousTime.pdf")

data = subsetData; xData = "timeCategory"; yData="CD4_.CD38.Ki67._FreqParent"; fillParam = "Prior.COVID.infection."; groupby="Record.ID"; title = "Activated CD4"; 
xLabel = " "; yLabel = "Ki67+CD38+ (% CD4)"; repMeasures = F; exponential=F; newform = T

prePostTime(data = subsetData, xData = "timeCategory", yData="CD4_.CD38.Ki67._FreqParent", fillParam = "Prior.COVID.infection.", groupby="Record.ID", title = "Activated CD4", 
            xLabel = " ", yLabel = "Ki67+CD38+ (% CD4)", repMeasures = F, exponential=F, newform = T)  #+ 
# ggrepel::geom_text_repel(aes(label=Alias),size=3)
# ggsave(filename = "./Images/ActivCD4_bothCohorts_overTime.pdf")

fit <- aov(CD4_.CD38.Ki67._FreqParent ~ timeCategory, data=subset(subsetData, Prior.COVID.infection. == 'No')); tukey_hsd(fit)
fit <- aov(CD4_.CD38.Ki67._FreqParent ~ timeCategory, data=subset(subsetData, Prior.COVID.infection. == 'Yes')); tukey_hsd(fit)


prePostTime(subsetData, xData = "timeCategory", yData="CD8_.CD38.Ki67._FreqParent", fillParam = "Prior.COVID.infection.", groupby="Record.ID", title = "Activated CD8", 
            xLabel = " ", yLabel = "Ki67+CD38+ (% CD8)", repMeasures = F, exponential=F, newform = T) + scale_y_continuous(breaks=seq(0,10,1), limits = c(0,5))
# ggrepel::geom_text_repel(aes(label=Alias),size=3)
# ggsave(filename = "./Images/ActivCD8_bothCohorts_overTime.pdf")

fit <- aov(CD8_.CD38.Ki67._FreqParent ~ timeCategory, data=subset(subsetData, Prior.COVID.infection. == 'No')); tukey_hsd(fit); formatC(tukey_hsd(fit)$p.adj, format="e")
fit <- aov(CD8_.CD38.Ki67._FreqParent ~ timeCategory, data=subset(subsetData, Prior.COVID.infection. == 'Yes')); tukey_hsd(fit); formatC(tukey_hsd(fit)$p.adj, format="e")

subsetData <- subset(mergedData, Record.ID != "CV-011" & Record.ID != "CV-012" & Record.ID != "CV-005")        # absence of Ki67 stain
subsetData <- subset(subsetData, timeCategory == "Baseline")
twoSampleBar(data = subsetData, xData = "Prior.COVID.infection.", yData = "FCActivCD4_Vax1", fillParam = "Prior.COVID.infection.", title = "Activ CD4", 
             yLabel = "Fold-change at one week", nonparam = T)


# ------------------------------------------------ Healthy side --------------------------------------------------------
subsetData = subset(mergedData, mergedData$shortForm == 'bL' | mergedData$shortForm=='oW'); subsetData = subset(subsetData, subsetData$Prior.COVID.infection.=="No")
subsetData <- subset(subsetData, Record.ID != "CV-011" & Record.ID != "CV-012" )        # extreme outliers because of absence of Ki67 stain
prePostTime(data=subsetData, xData = "shortForm", yData="CD4_.CD38.Ki67._FreqParent", fillParam = "Prior.COVID.infection.", groupby="Record.ID", title = "Activated CD4", 
            xLabel = " ", yLabel = "Ki67+CD38+ (% CD4)") + #ggrepel::geom_text_repel(aes(label=Alias)) + 
  scale_y_continuous(breaks = seq(0,20,1), limits = c(0,5)) + scale_x_discrete(labels= c("Baseline","Post\n1st dose"))
# ggsave(filename = "./Images/ActivatedCD4_vax1_healthy_paired.pdf", width=4)

subsetData = subset(mergedData, mergedData$shortForm == 'bL' | mergedData$shortForm=='oW'); subsetData = subset(subsetData, subsetData$Prior.COVID.infection.=="No")
subsetData <- subset(subsetData, Record.ID != "CV-011" & Record.ID != "CV-012" )        # extreme outliers because of absence of Ki67 stain
prePostTime(data=subsetData, xData = "shortForm", yData="CD8_.CD38.Ki67._FreqParent", fillParam = "Prior.COVID.infection.", groupby="Record.ID", title = "Activated CD8", 
            xLabel = " ", yLabel = "Ki67+CD38+ (% CD8)") + 
  scale_y_continuous(breaks = seq(0,20,1), limits = c(0,6)) + scale_x_discrete(labels= c("Baseline","Post\n1st dose"))
# ggsave(filename = "./Images/ActivatedCD8_vax1_healthy_paired.pdf", width=4)


subsetData = subset(mergedData, mergedData$shortForm == '3W' | mergedData$shortForm=='4W'); subsetData = subset(subsetData, subsetData$Prior.COVID.infection.=="No")
prePostTime(data=subsetData, xData = "shortForm", yData="CD4_.CD38.Ki67._FreqParent", fillParam = "Prior.COVID.infection.", groupby="Record.ID", title = "Activated CD4", 
            xLabel = " ", yLabel = "Ki67+CD38+ (% CD4)") + scale_y_continuous(breaks = seq(0,20,1), limits = c(0,5)) + scale_x_discrete(labels= c("Pre\n2nd dose","Post\n2nd dose"))
# ggsave(filename = "./Images/ActivatedCD4_vax2_healthy_paired.pdf", width=4)

subsetData = subset(mergedData, mergedData$shortForm == '3W' | mergedData$shortForm=='4W'); subsetData = subset(subsetData, subsetData$Prior.COVID.infection.=="No")
prePostTime(data=subsetData, xData = "shortForm", yData="CD8_.CD38.Ki67._FreqParent", fillParam = "Prior.COVID.infection.", groupby="Record.ID", title = "Activated CD8", 
            xLabel = " ", yLabel = "Ki67+CD38+ (% CD8)") + scale_y_continuous(breaks = seq(0,20,1), limits = c(0,6)) + scale_x_discrete(labels= c("Pre\n2nd dose","Post\n2nd dose"))
# ggsave(filename = "./Images/ActivatedCD8_vax2_healthy_paired.pdf", width=4)

# ------------------------------------------------ PHI side --------------------------------------------------------
subsetData = subset(mergedData, mergedData$shortForm == 'bL' | mergedData$shortForm=='oW'); subsetData = subset(subsetData, subsetData$Prior.COVID.infection.=="Yes")
prePostTime(data=subsetData, xData = "shortForm", yData="CD4_.CD38.Ki67._FreqParent", fillParam = "Prior.COVID.infection.", groupby="Record.ID", title = "Activated CD4", 
            xLabel = " ", yLabel = "Ki67+CD38+ (% CD4)") + scale_y_continuous(breaks = seq(0,20,1), limits = c(0,5)) + 
  scale_x_discrete(labels= c("Baseline","Post\n1st dose"))+ scale_fill_manual(values=c("#B5B2F1"))# + ggrepel::geom_text_repel(aes(label=Alias))
# ggsave(filename = "./Images/ActivatedCD4_vax1_covid_paired.pdf", width=4)

subsetData = subset(mergedData, mergedData$shortForm == 'bL' | mergedData$shortForm=='oW'); subsetData = subset(subsetData, subsetData$Prior.COVID.infection.=="Yes")
prePostTime(data=subsetData, xData = "shortForm", yData="CD8_.CD38.Ki67._FreqParent", fillParam = "Prior.COVID.infection.", groupby="Record.ID", title = "Activated CD8", 
            xLabel = " ", yLabel = "Ki67+CD38+ (% CD8)") + 
  scale_y_continuous(breaks = seq(0,20,1), limits = c(0,6)) + scale_x_discrete(labels= c("Baseline","Post\n1st dose"))+ scale_fill_manual(values=c("#B5B2F1")) 
# ggsave(filename = "./Images/ActivatedCD8_vax1_covid_paired.pdf", width=4)


subsetData = subset(mergedData, mergedData$shortForm == '3W' | mergedData$shortForm=='4W'); subsetData = subset(subsetData, subsetData$Prior.COVID.infection.=="Yes")
prePostTime(data=subsetData, xData = "shortForm", yData="CD4_.CD38.Ki67._FreqParent", fillParam = "Prior.COVID.infection.", groupby="Record.ID", title = "Activated CD4", 
            xLabel = " ", yLabel = "Ki67+CD38+ (% CD4)") + scale_y_continuous(breaks = seq(0,20,1), limits = c(0,5)) + 
  scale_x_discrete(labels= c("Pre\n2nd dose","Post\n2nd dose"))+ scale_fill_manual(values=c("#B5B2F1")) 
# ggsave(filename = "./Images/ActivatedCD4_vax2_covid_paired.pdf", width=4)

subsetData = subset(mergedData, mergedData$shortForm == '3W' | mergedData$shortForm=='4W'); subsetData = subset(subsetData, subsetData$Prior.COVID.infection.=="Yes")
prePostTime(data=subsetData, xData = "shortForm", yData="CD8_.CD38.Ki67._FreqParent", fillParam = "Prior.COVID.infection.", groupby="Record.ID", title = "Activated CD8", 
            xLabel = " ", yLabel = "Ki67+CD38+ (% CD8)") + scale_y_continuous(breaks = seq(0,20,1), limits = c(0,6)) + 
  scale_x_discrete(labels= c("Pre\n2nd dose","Post\n2nd dose")) + scale_fill_manual(values=c("#B5B2F1")) 
# ggsave(filename = "./Images/ActivatedCD8_vax2_covid_paired.pdf", width=4)



#' ------------------ CXCL13 analyses --------------------------
#'

linePlot(data = mergedData, xData = 'timeCategory', yData = 'CXCL13', groupby = 'Alias', xLabel = ' ', yLabel = "CXCL13 (pg/mL)", 
         title = "Plasma CXCL13", colorby = "Prior.COVID.infection.") + theme(axis.title.x = element_blank()) + 
  scale_color_manual(name="Prior COVID?",values = c("#FFDFB1","#B5B2F1")) + 
  scale_y_continuous(limits = c(0,140),breaks=seq(0,140,10))
  # scale_y_continuous(trans='pseudo_log', limits = c(0,1e4), breaks=c(10^(0:7)), labels=trans_format('log10',math_format(10^.x)), minor_breaks =5*10^(0:10))# + ggrepel::geom_text_repel(aes(label = Label))
# ggsave(filename = "./Images/CXCL13plasma_linePlot.pdf")

acuteCOVID <- read.csv(file = "D:/COVID/Analysis/Rscripts/COVIDmergedData.csv")
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


#' ------------------ Avidity analyses --------------------------
#'

linePlot(data = mergedData, xData = 'timeCategory', yData = 'Avidity', groupby = 'Alias', xLabel = ' ', yLabel = "IgG avidity (%)", 
         title = "anti-S1 Avidity", colorby = "Prior.COVID.infection.") + theme(axis.title.x = element_blank()) + 
  scale_color_manual(name="Prior COVID?",values = c("#FFDFB1","#B5B2F1")) + scale_y_continuous(limits = c(0,110),breaks=seq(0,140,10))
# ggsave(filename = "./Images/avidity_lineplot.pdf")



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
            xLabel = " ", yLabel = "ICOS+CD38+ (% cTfh)", repMeasures = F, newform=T) #  + ggrepel::geom_text_repel(aes(label=Record.ID),size=2)
# ggsave(filename = "./Images/cTfh_responses_bothCohorts_overTime.pdf")

fit <- aov(CD4_.Nonnaive.cTfh.ICOS..CD38.._FreqParent ~ timeCategory, data=subset(subsetData, Prior.COVID.infection. == 'No')); tukey_hsd(fit); formatC(tukey_hsd(fit)$p.adj, format="e")
fit <- aov(CD4_.Nonnaive.cTfh.ICOS..CD38.._FreqParent ~ timeCategory, data=subset(subsetData, Prior.COVID.infection. == 'Yes')); tukey_hsd(fit); formatC(tukey_hsd(fit)$p.adj, format="e")

  
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
# ggrepel::geom_text_repel(aes(label=Record.ID),size=2)
# ggsave(filename = "./Images/cTfh_responses_CXCR3_bothCohorts_overTime.pdf")
fit <- aov(CD4_.Nonnaive.cTfh.ICOS..CD38...CXCR3._FreqParent ~ timeCategory, data=subset(subsetData, Prior.COVID.infection. == 'No')); tukey_hsd(fit); formatC(tukey_hsd(fit)$p.adj, format="e")
fit <- aov(CD4_.Nonnaive.cTfh.ICOS..CD38...CXCR3._FreqParent ~ timeCategory, data=subset(subsetData, Prior.COVID.infection. == 'Yes')); tukey_hsd(fit); formatC(tukey_hsd(fit)$p.adj, format="e")


bivScatter(data1 = subset(mergedData, Prior.COVID.infection. == 'No' & timeCategory == "Post 1st dose" & Tube == "HEP"), name1 = "Naive", 
           data2 = subset(mergedData, Prior.COVID.infection. == 'Yes' & timeCategory == "Post 1st dose" & Tube == "HEP"), name2 = "Experienced", 
           xData = "Age", yData = "CD4_.Nonnaive.cTfh.ICOS..CD38.._FreqParent", fillParam = "Prior.COVID.infection.", title = "Vaccine dose 1", 
           xLabel = "Age", yLabel = "ICOS+CD38+ (% cTfh)", nonparam = T) +   scale_y_continuous(breaks=seq(0,27,3), limits = c(0,25) )+ scale_x_continuous(limits = c(20,70))
# ggsave(filename = "./Images/Age_correl_HiHiTfh_Vax1.pdf", width=8)
bivScatter(data1 = subset(mergedData, Prior.COVID.infection. == 'No' & timeCategory == "Post 2nd dose" & Tube == "HEP"), name1 = "Naive", 
           data2 = subset(mergedData, Prior.COVID.infection. == 'Yes' & timeCategory == "Post 2nd dose" & Tube == "HEP"), name2 = "Experienced", 
           xData = "Age", yData = "CD4_.Nonnaive.cTfh.ICOS..CD38.._FreqParent", fillParam = "Prior.COVID.infection.", title = "Vaccine dose 2", 
           xLabel = "Age", yLabel = "ICOS+CD38+ (% cTfh)", nonparam = T) +   scale_y_continuous(breaks=seq(0,27,3), limits = c(0,25) ) + scale_x_continuous(limits = c(20,70))
# ggsave(filename = "./Images/Age_correl_HiHiTfh_Vax2.pdf", width=8)

bivScatter(data1 = subset(mergedData, Prior.COVID.infection. == 'No' & timeCategory == "Post 1st dose" & Tube == "HEP"), name1 = "Naive", 
           data2 = subset(mergedData, Prior.COVID.infection. == 'Yes' & timeCategory == "Post 1st dose" & Tube == "HEP"), name2 = "Experienced", 
           xData = "Age", yData = "FCtfh_Vax1", fillParam = "Prior.COVID.infection.", title = "Vaccine dose 1", 
           xLabel = "Age", yLabel = "Fold-change ICOS+CD38+ cTfh", nonparam = T) +scale_y_continuous(breaks=seq(0,30,1), limits = c(0,6) )+ scale_x_continuous(limits = c(20,70))
# ggsave(filename = "./Images/Age_correl_FCtfh_Vax1.pdf", width=8)
bivScatter(data1 = subset(mergedData, Prior.COVID.infection. == 'No' & timeCategory == "Post 2nd dose" & Tube == "HEP"), name1 = "Naive", 
           data2 = subset(mergedData, Prior.COVID.infection. == 'Yes' & timeCategory == "Post 2nd dose" & Tube == "HEP"), name2 = "Experienced", 
           xData = "Age", yData = "FCtfh_Vax2", fillParam = "Prior.COVID.infection.", title = "Vaccine dose 2", 
           xLabel = "Age", yLabel = "Fold-change ICOS+CD38+ cTfh", nonparam = T) +   scale_y_continuous(breaks=seq(0,30,1), limits = c(0,6) ) + scale_x_continuous(limits = c(20,70))
# ggsave(filename = "./Images/Age_correl_FCtfh_Vax2.pdf", width=8)



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

fit <- aov(CD19_.CD27..CD38._FreqParent ~ timeCategory, data=subset(subsetData, Prior.COVID.infection. == 'No')); tukey_hsd(fit); formatC(tukey_hsd(fit)$p.adj, format="e")
fit <- aov(CD19_.CD27..CD38._FreqParent ~ timeCategory, data=subset(subsetData, Prior.COVID.infection. == 'Yes')); tukey_hsd(fit); formatC(tukey_hsd(fit)$p.adj, format="e")


subsetData = subset(mergedData, mergedData$shortForm == 'bL' | mergedData$shortForm=='oW'); subsetData = subset(subsetData, subsetData$Prior.COVID.infection.=="No")
prePostTime(data=subsetData, xData = "shortForm", yData="CD19_.CD27..CD38._FreqParent", fillParam = "Prior.COVID.infection.", groupby="Record.ID", title = "PB - Vax1", 
            xLabel = " ", yLabel = "CD27+CD38+ (% CD19)") + scale_y_continuous(breaks = seq(0,21,0.5), limits = c(0,2))  + 
  scale_x_discrete(labels= c("Baseline","Post\n1st dose"))# + ggrepel::geom_text_repel(aes(label = Record.ID))
# ggsave(filename = "./Images/PB_vax1_healthy_paired.pdf", width=4)


subsetData = subset(mergedData, mergedData$shortForm == '3W' | mergedData$shortForm=='4W'); subsetData = subset(subsetData, subsetData$Prior.COVID.infection.=="No")
prePostTime(data=subsetData, xData = "shortForm", yData="CD19_.CD27..CD38._FreqParent", fillParam = "Prior.COVID.infection.", groupby="Record.ID", 
            title = "PB - Vax2", xLabel = " ", yLabel = "CD27+CD38+ (% CD19)") + scale_y_continuous(breaks = seq(0,21,0.5), limits = c(0,2))  + 
  scale_x_discrete(labels= c("Pre\n2nd dose","Post\n2nd dose"))# + ggrepel::geom_text_repel(aes(label = Record.ID))
# ggsave(filename = "./Images/PB_vax2_healthy_paired.pdf", width=4)



subsetData = subset(mergedData, mergedData$shortForm == 'bL' | mergedData$shortForm=='oW'); subsetData = subset(subsetData, subsetData$Prior.COVID.infection.=="Yes")
prePostTime(data=subsetData, xData = "shortForm", yData="CD19_.CD27..CD38._FreqParent", fillParam = "Prior.COVID.infection.", groupby="Record.ID", 
            title = "PB - Vax1", xLabel = " ", yLabel = "CD27+CD38+ (% CD19)") +    
  scale_y_continuous(breaks = seq(0,21,0.5), limits = c(0,3))  + scale_x_discrete(labels= c("Baseline","Post\n1st dose"))+ scale_fill_manual(values=c("#B5B2F1")) 
# ggsave(filename = "./Images/PB_vax1_covid_paired.pdf", width=4)

subsetData = subset(mergedData, mergedData$shortForm == '3W' | mergedData$shortForm=='4W'); subsetData = subset(subsetData, subsetData$Prior.COVID.infection.=="Yes")
prePostTime(data=subsetData, xData = "shortForm", yData="CD19_.CD27..CD38._FreqParent", fillParam = "Prior.COVID.infection.", groupby="Record.ID", 
            title = "PB - Vax2", xLabel = " ", yLabel = "CD27+CD38+ (% CD19)") +    
  scale_y_continuous(breaks = seq(0,21,0.5), limits = c(0,3))  + scale_x_discrete(labels= c("Pre\n2nd dose","Post\n2nd dose"))+ scale_fill_manual(values=c("#B5B2F1")) 
# ggsave(filename = "./Images/PB_vax2_covid_paired.pdf", width=4)

subsetData = subset(mergedData, mergedData$shortForm == 'bL' ); subsetData = subset(subsetData, subsetData$Prior.COVID.infection.=="No")
univScatter(data = subsetData, xData = "FCtfh_Vax1", yData = "FCPB_Vax1", fillParam = "Prior.COVID.infection.", title = " ", xLabel = "Fold-change Tfh", 
            yLabel = "Fold-change PB", position = "left")
subsetData = subset(mergedData, mergedData$shortForm == 'bL' ); subsetData = subset(subsetData, subsetData$Prior.COVID.infection.=="Yes")
univScatter(data = subsetData, xData = "FCtfh_Vax1", yData = "FCPB_Vax1", fillParam = "Prior.COVID.infection.", title = " ", xLabel = "Fold-change Tfh", 
            yLabel = "Fold-change PB", position = "left")+ scale_fill_manual(values=c("#B5B2F1")) 


subsetData = subset(mergedData, mergedData$shortForm == 'oW')
subsetData %>% group_by(Prior.COVID.infection.) %>% shapiro_test(CD19_.CD27..CD38._FreqParent)
twoSampleBar(data=subsetData, xData = "Prior.COVID.infection.", yData="CD19_.CD27..CD38._FreqParent", fillParam = "Prior.COVID.infection.", title = "PB post 1st dose", 
             yLabel = "CD27+CD38+ (% CD19)", nonparam = T) + scale_y_continuous(breaks = seq(0,20,1), limits = c(0,2)) 

subsetData = subset(mergedData, mergedData$shortForm == '4W')
subsetData %>% group_by(Prior.COVID.infection.) %>% shapiro_test(CD19_.CD27..CD38._FreqParent)
twoSampleBar(data=subsetData, xData = "Prior.COVID.infection.", yData="CD19_.CD27..CD38._FreqParent", fillParam = "Prior.COVID.infection.", title = "PB post 2nd dose", 
             yLabel = "CD27+CD38+ (% CD19)", nonparam = T) + scale_y_continuous(breaks = seq(0,20,1), limits = c(0,2)) 


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
fit <- aov(CD19_.CD27..CD38..CD138._FreqParent ~ timeCategory, data=subset(subsetData, Prior.COVID.infection. == 'No')); tukey_hsd(fit); formatC(tukey_hsd(fit)$p.adj, format="e")
fit <- aov(CD19_.CD27..CD38..CD138._FreqParent ~ timeCategory, data=subset(subsetData, Prior.COVID.infection. == 'Yes')); tukey_hsd(fit); formatC(tukey_hsd(fit)$p.adj, format="e")





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
  scale_y_continuous(limits = c(0,110), breaks = seq(0,100,10))
  # ggformula::geom_spline(aes_string(x="DPV", y="CD8_.CD38.Ki67..GzmB..CD8_FreqParent"), color="blue",size=1, spar=0.5) #+ ggrepel::geom_text_repel(aes(label=Record.ID),size=2)
# ggsave(filename = "./Images/CD8_CD38hiKi67hi_GzmB.pdf" )

fit <- aov(CD8_.CD38.Ki67..GzmB..CD8_FreqParent ~ timeCategory, data=subset(subsetData, Prior.COVID.infection. == 'No')); tukey_hsd(fit); formatC(tukey_hsd(fit)$p.adj, format="e")
fit <- aov(CD8_.CD38.Ki67..GzmB..CD8_FreqParent ~ timeCategory, data=subset(subsetData, Prior.COVID.infection. == 'Yes')); tukey_hsd(fit); formatC(tukey_hsd(fit)$p.adj, format="e")

  
prePostTime(data = subsetData, xData = "timeCategory", yData="CD4_.CD38.Ki67..GzmB..CD4_FreqParent", fillParam = "Prior.COVID.infection.", groupby="Record.ID", title = "Activated CD4", 
            xLabel = "Prior COVID?", yLabel = "GzmB+ (% CD4+Ki67+CD38+)", repMeasures=F) #+ 
# ggformula::geom_spline(aes_string(x="DPV", y="CD4_.CD38.Ki67..GzmB..CD4_FreqParent"), color="blue",size=1, spar=0.5) #+ ggrepel::geom_text_repel(aes(label=Record.ID),size=2)

fit <- aov(CD4_.CD38.Ki67..GzmB..CD4_FreqParent ~ timeCategory, data=subset(subsetData, Prior.COVID.infection. == 'No')); tukey_hsd(fit); formatC(tukey_hsd(fit)$p.adj, format="e")
fit <- aov(CD4_.CD38.Ki67..GzmB..CD4_FreqParent ~ timeCategory, data=subset(subsetData, Prior.COVID.infection. == 'Yes')); tukey_hsd(fit); formatC(tukey_hsd(fit)$p.adj, format="e")


subsetData = subset(mergedData, mergedData$shortForm == 'bL' | mergedData$shortForm=='oW'); subsetData = subset(subsetData, subsetData$Prior.COVID.infection.=="No")
prePostTime(data=subsetData, xData = "shortForm", yData="CD8_.CD38.Ki67..GzmB..CD8_FreqParent", fillParam = "Prior.COVID.infection.", groupby="Record.ID", title = "GzmB+ activated CD8 pre and post 1st dose", 
            xLabel = " ", yLabel = "GzmB+ (% Activated CD8)") + scale_y_continuous(breaks = seq(0,100,10), limits = c(0,100)) + scale_x_discrete(labels= c("Baseline","Post\n1st dose")) 
# ggsave(filename = "./Images/GzmBActivatedCD8_vax1_healthy_paired.pdf", width=4)


subsetData = subset(mergedData, mergedData$shortForm == 'bL' | mergedData$shortForm=='oW'); subsetData = subset(subsetData, subsetData$Prior.COVID.infection.=="Yes")
prePostTime(data=subsetData, xData = "shortForm", yData="CD8_.CD38.Ki67..GzmB..CD8_FreqParent", fillParam = "Prior.COVID.infection.", groupby="Record.ID", title = "GzmB+ activated CD8 pre and post 1st dose", 
            xLabel = " ", yLabel = "GzmB+ (% Activated CD8)") + scale_y_continuous(breaks = seq(0,100,10), limits = c(0,100)) + scale_x_discrete(labels= c("Baseline","Post\n1st dose")) + scale_fill_manual(values=c("#B5B2F1"))
# ggsave(filename = "./Images/GzmBActivatedCD8_vax1_covid_paired.pdf", width=4)


subsetData = subset(mergedData, mergedData$shortForm == '3W' | mergedData$shortForm=='4W'); subsetData = subset(subsetData, subsetData$Prior.COVID.infection.=="No")
prePostTime(data=subsetData, xData = "shortForm", yData="CD8_.CD38.Ki67..GzmB..CD8_FreqParent", fillParam = "Prior.COVID.infection.", groupby="Record.ID", title = "GzmB+ activated CD8 pre and post 2nd dose", 
            xLabel = " ", yLabel = "GzmB+ (% Activated CD8)") + scale_y_continuous(breaks = seq(0,100,10), limits = c(0,100)) + scale_x_discrete(labels= c("Pre\n2nd dose","Post\n2nd dose")) 
# ggsave(filename = "./Images/GzmBActivatedCD8_vax2_healthy_paired.pdf", width=4)


subsetData = subset(mergedData, mergedData$shortForm == '3W' | mergedData$shortForm=='4W'); subsetData = subset(subsetData, subsetData$Prior.COVID.infection.=="Yes")
prePostTime(data=subsetData, xData = "shortForm", yData="CD8_.CD38.Ki67..GzmB..CD8_FreqParent", fillParam = "Prior.COVID.infection.", groupby="Record.ID", title = "GzmB+ activated CD8 pre and post 2nd dose", 
            xLabel = " ", yLabel = "GzmB+ (% Activated CD8)") + scale_y_continuous(breaks = seq(0,100,10), limits = c(0,100)) + scale_x_discrete(labels= c("Pre\n2nd dose","Post\n2nd dose")) + scale_fill_manual(values=c("#B5B2F1"))
# ggsave(filename = "./Images/GzmBActivatedCD8_vax2_covid_paired.pdf", width=4)

subsetData <- subset(mergedData,  timeCategory != "two Weeks");  subsetData$timeCategory <- factor(subsetData$timeCategory, levels = c("Baseline", "Post 1st dose", "Pre 2nd dose", "Post 2nd dose","One month post\n2nd dose"))
prePostTime(subsetData, xData = "timeCategory", yData="CD8_.CD38.Ki67..GzmB..CD8_FreqParent", fillParam = "Prior.COVID.infection.", groupby="Record.ID", title = "GzmB+ Activated CD8", 
            xLabel = "Prior COVID?", yLabel = "GzmB+ (% Activated CD8)", repMeasures = F)  + #ggrepel::geom_text_repel(aes(label=Record.ID),size=2)
# ggsave(filename = "./Images/GzmBActivatedCD8_bothCohorts_overTime.pdf", width=4)

subsetData = subset(mergedData, mergedData$shortForm == 'oW')
subsetData %>% group_by(Prior.COVID.infection.) %>% shapiro_test(CD8_.CD38.Ki67..GzmB..CD8_FreqParent)
twoSampleBar(data=subsetData, xData = "Prior.COVID.infection.", yData="CD8_.CD38.Ki67..GzmB..CD8_FreqParent", fillParam = "Prior.COVID.infection.", title = "GzmB+ Activated CD8 post 1st dose", 
             yLabel = "GzmB+ (% Activated CD8)", nonparam = T) + scale_y_continuous(breaks = seq(0,100,10), limits = c(0,100)) 

subsetData = subset(mergedData, mergedData$shortForm == '4W')
subsetData %>% group_by(Prior.COVID.infection.) %>% shapiro_test(CD8_.CD38.Ki67..GzmB..CD8_FreqParent)
twoSampleBar(data=subsetData, xData = "Prior.COVID.infection.", yData="CD8_.CD38.Ki67..GzmB..CD8_FreqParent", fillParam = "Prior.COVID.infection.", title = "GzmB+ Activated CD8 post 2nd dose", 
             yLabel = "GzmB+ (% Activated CD8)", nonparam = T) + scale_y_continuous(breaks = seq(0,100,10), limits = c(0,100)) 

subsetData = subset(mergedData, mergedData$shortForm == 'oM')
subsetData %>% group_by(Prior.COVID.infection.) %>% shapiro_test(CD8_.CD38.Ki67..GzmB..CD8_FreqParent)
twoSampleBar(data=subsetData, xData = "Prior.COVID.infection.", yData="CD8_.CD38.Ki67..GzmB..CD8_FreqParent", fillParam = "Prior.COVID.infection.", title = "GzmB+ Activated CD8 at one month", 
             yLabel = "GzmB+ (% Activated CD8)", nonparam = T) + scale_y_continuous(breaks = seq(0,100,10), limits = c(0,100)) 



#' ------------------ GzmB+Tbet+ CD8 analyses --------------------------
subsetData <- subset(mergedData,  timeCategory != "two Weeks");  
subsetData$timeCategory <- factor(subsetData$timeCategory, levels = c("Baseline", "Post 1st dose", "Pre 2nd dose", "Post 2nd dose", "One month post\n2nd dose"))

prePostTime(data = subsetData, xData = "timeCategory", yData="CD8_.CD38.Ki67..Tbet.GzmB._FreqParent", fillParam = "Prior.COVID.infection.", groupby="Record.ID", title = "Activated CD8", 
            xLabel = "Prior COVID?", yLabel = "GzmB+Tbet+ (% CD8+Ki67+CD38+)", repMeasures=F) #+ 
# ggformula::geom_spline(aes_string(x="DPV", y="CD8_.CD38.Ki67..GzmB..CD8_FreqParent"), color="blue",size=1, spar=0.5) #+ ggrepel::geom_text_repel(aes(label=Record.ID),size=2)
prePostTime(data = subsetData, xData = "timeCategory", yData="CD8_.CD38.Ki67..GzmB..CD8_FreqParent", fillParam = "Prior.COVID.infection.", groupby="Record.ID", title = "Activated CD8", 
            xLabel = "Prior COVID?", yLabel = "GzmB+ (% CD8+Ki67+CD38+)", repMeasures=F) #+ 
# ggformula::geom_spline(aes_string(x="DPV", y="CD8_.CD38.Ki67..GzmB..CD8_FreqParent"), color="blue",size=1, spar=0.5) #+ ggrepel::geom_text_repel(aes(label=Record.ID),size=2)

fit <- aov(CD8_.CD38.Ki67..GzmB..CD8_FreqParent ~ timeCategory, data=subset(subsetData, Prior.COVID.infection. == 'No')); tukey_hsd(fit); formatC(tukey_hsd(fit)$p.adj, format="e")
fit <- aov(CD8_.CD38.Ki67..GzmB..CD8_FreqParent ~ timeCategory, data=subset(subsetData, Prior.COVID.infection. == 'Yes')); tukey_hsd(fit); formatC(tukey_hsd(fit)$p.adj, format="e")
fit <- aov(CD8_.CD38.Ki67..Tbet.GzmB._FreqParent ~ timeCategory, data=subset(subsetData, Prior.COVID.infection. == 'No')); tukey_hsd(fit); formatC(tukey_hsd(fit)$p.adj, format="e")
fit <- aov(CD8_.CD38.Ki67..Tbet.GzmB._FreqParent ~ timeCategory, data=subset(subsetData, Prior.COVID.infection. == 'Yes')); tukey_hsd(fit); formatC(tukey_hsd(fit)$p.adj, format="e")



#' ------------------ Foxp3+ CD4 analyses --------------------------
subsetData <- subset(mergedData,  timeCategory != "two Weeks" );  
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
# subsetData <- subsetData[-which(subsetData$Record.ID == "CV-033"),]     # exclude because missing file for Pre 2nd dose 

subsetData$timeCategory <- factor(subsetData$timeCategory, levels = c("Baseline", "Post 1st dose", "Pre 2nd dose", "Post 2nd dose", "One month post\n2nd dose"))
subsetData <- subsetData[which(subsetData$Tube == "HEP"),]
prePostTime(data = subsetData, xData = "timeCategory", yData="CD19_.IgD..CD71._FreqParent", fillParam = "Prior.COVID.infection.", groupby="Record.ID", 
            title = "CD71+ B cells", xLabel = "", yLabel = "IgD-CD71+ (% CD19+)", repMeasures=F, newform = T) #ggrepel::geom_text_repel(aes(label = Alias))
# ggsave(filename = "./Images/CD71hi_Bcells_bothCohorts_overTime.pdf")
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
fit <- aov(CD19_.CD21lo_FreqParent ~ timeCategory, data=subset(subsetData, Prior.COVID.infection. == 'No')); tukey_hsd(fit); formatC(tukey_hsd(fit)$p.adj, format="e")
fit <- aov(CD19_.CD21lo_FreqParent ~ timeCategory, data=subset(subsetData, Prior.COVID.infection. == 'Yes')); tukey_hsd(fit); formatC(tukey_hsd(fit)$p.adj, format="e")



#' ------------------ Correlational analyses --------------------------
#' 
index <- grep(paste(c("^Prior","DPO.covid","DPV","timeCategory", "FreqParent$", "^FC"), collapse = "|"), names(mergedData), value=F)
mergedData.deident <- mergedData[, index]
subsetData <- mergedData.deident[which(mergedData.deident$timeCategory=="Post 1st dose"), ]
GGally::ggduo( data = subsetData, mapping=ggplot2::aes(color=Prior.COVID.infection.))+ theme_bw() + ggplot2::theme(legend.position = "bottom") +    scale_color_manual(values=c("#FFDFB1", "#B5B2F1")) + 
  scale_fill_manual(values=c("#FFDFB1", "#B5B2F1")) 
# ggsave("./Images/exploratory_post1st.pdf", width=50, height=50, limitsize = F)

subsetData <- mergedData.deident[which(mergedData.deident$timeCategory=="Post 2nd dose"), ]
GGally::ggduo( data = subsetData, mapping=ggplot2::aes(color=Prior.COVID.infection.))+ theme_bw() + ggplot2::theme(legend.position = "bottom") +    scale_color_manual(values=c("#FFDFB1", "#B5B2F1")) + 
  scale_fill_manual(values=c("#FFDFB1", "#B5B2F1")) 
# ggsave("./Images/exploratory_post2nd.pdf", width=50, height=50, limitsize = F)

ggplot(data = subset(mergedData, timeCategory == "Post 1st dose"), aes(x = Elispot_IgM_RBD, y=CD19_.CD27..CD38._FreqParent, color=Prior.COVID.infection. )) + geom_point(size=5) + theme_bw()
ggplot(data = subset(mergedData, timeCategory == "Post 1st dose"), aes(x = Elispot_IgG_S1, y=FCActivCD4_Vax1, color=Prior.COVID.infection. )) + geom_point(size=5) + theme_bw()






bivScatter(data1 = subset(mergedData, Prior.COVID.infection. == 'No' & timeCategory == "Post 1st dose" & Tube == "HEP"), name1 = "Naive", 
           data2 = subset(mergedData, Prior.COVID.infection. == 'Yes' & timeCategory == "Post 1st dose" & Tube == "HEP"), name2 = "Experienced", 
           xData = "FCActivCD4_Vax1", yData = "FCActivCD8_Vax1", fillParam = "Prior.COVID.infection.", title = "Vaccine dose 1", 
           xLabel = "Fold-change CD4+Ki67+CD38+", yLabel = "Fold-change CD8+Ki67+CD38+", nonparam = T) + 
  scale_x_continuous(trans='pseudo_log', breaks=c(10^(0:4)), labels=trans_format('log10',math_format(10^.x)), minor_breaks =5*10^(0:10)) + 
  scale_y_continuous(trans='pseudo_log', breaks=c(10^(0:4)), labels=trans_format('log10',math_format(10^.x)), minor_breaks =5*10^(0:10))
# ggsave(filename = "./Images/FCActivCD4_correl_FCActivCD8_Vax1.pdf", width=8)
bivScatter(data1 = subset(mergedData, Prior.COVID.infection. == 'No' & timeCategory == "Post 2nd dose" & Tube == "HEP"), name1 = "Naive", 
           data2 = subset(mergedData, Prior.COVID.infection. == 'Yes' & timeCategory == "Post 2nd dose" & Tube == "HEP"), name2 = "Experienced", 
           xData = "FCActivCD4_Vax2", yData = "FCActivCD8_Vax2", fillParam = "Prior.COVID.infection.", title = "Vaccine dose 2", 
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





bivScatter(data1 = subset(mergedData, Prior.COVID.infection. == 'No' & timeCategory == "Post 1st dose" ), name1 = "Naive", 
           data2 = subset(mergedData, Prior.COVID.infection. == 'Yes' & timeCategory == "Post 1st dose"), name2 = "Experienced", 
           xData = "FCActivCD4_Vax1", yData = "FCtfh_Vax1", fillParam = "Prior.COVID.infection.", title = "Vaccine dose 1", 
           xLabel = "Fold-change CD4+Ki67+CD38+", yLabel = "Fold-change ICOS+CD38+ cTfh", nonparam = T) + 
  scale_x_continuous(trans='pseudo_log', breaks=c(10^(0:4)), labels=trans_format('log10',math_format(10^.x)), minor_breaks =5*10^(0:10)) + 
  scale_y_continuous(trans='pseudo_log', breaks=c(10^(0:4)), limits = c(0,20), labels=trans_format('log10',math_format(10^.x)), minor_breaks =5*10^(0:10))
# ggsave(filename = "./Images/FCTfh_correl_FCactivCD4_Vax1.pdf", width=8)


bivScatter(data1 = subset(mergedData, Prior.COVID.infection. == 'No' & timeCategory == "Post 1st dose"), name1 = "Naive", 
           data2 = subset(mergedData, Prior.COVID.infection. == 'Yes' & timeCategory == "Post 1st dose"), name2 = "Experienced", 
           xData = "FCActivCD4_Vax2", yData = "FCtfh_Vax2", fillParam = "Prior.COVID.infection.", title = "Vaccine dose 2", 
           xLabel = "Fold-change CD4+Ki67+CD38+", yLabel = "Fold-change ICOS+CD38+ cTfh", nonparam = T) + 
  scale_x_continuous(trans='pseudo_log', breaks=c(10^(0:4)), labels=trans_format('log10',math_format(10^.x)), minor_breaks =5*10^(0:10)) + 
  scale_y_continuous(trans='pseudo_log', breaks=c(10^(0:4)), labels=trans_format('log10',math_format(10^.x)), minor_breaks =5*10^(0:10))
# ggsave(filename = "./Images/FCTfh_correl_FCactivCD4_Vax2.pdf", width=8)










bivScatter(data1 = subset(mergedData, Prior.COVID.infection. == 'No' & timeCategory == "Post 1st dose" & Tube == "HEP"), name1 = "Naive", 
           data2 = subset(mergedData, Prior.COVID.infection. == 'Yes' & timeCategory == "Post 1st dose" & Tube == "HEP"), name2 = "Experienced", 
           xData = "Age", yData = "CD4_.CD38.Ki67._FreqParent", fillParam = "Prior.COVID.infection.", title = "Vaccine dose 2", 
           xLabel = "Fold-change CD4+Ki67+CD38+", yLabel = "Fold-change CD8+Ki67+CD38+", nonparam = T) + 
  scale_x_continuous(limits = c(20,70), breaks=seq(0,100,10)) + scale_y_continuous(limits = c(0,5), breaks=seq(0,10,1))
# ggsave(filename = "./Images/FCActivCD4_correl_FCActivCD8_Vax2.pdf", width=8)

bivScatter(data1 = subset(mergedData, Prior.COVID.infection. == 'No' & timeCategory == "Post 2nd dose" & Tube == "HEP"), name1 = "Naive", 
           data2 = subset(mergedData, Prior.COVID.infection. == 'Yes' & timeCategory == "Post 2nd dose" & Tube == "HEP"), name2 = "Experienced", 
           xData = "Age", yData = "CD4_.CD38.Ki67._FreqParent", fillParam = "Prior.COVID.infection.", title = "Vaccine dose 2", 
           xLabel = "Fold-change CD4+Ki67+CD38+", yLabel = "Fold-change CD8+Ki67+CD38+", nonparam = T) + 
  scale_x_continuous(limits = c(20,70), breaks=seq(0,100,10)) + scale_y_continuous(limits = c(0,5), breaks=seq(0,10,1))
# ggsave(filename = "./Images/FCActivCD4_correl_FCActivCD8_Vax2.pdf", width=8)


bivScatter(data1 = subset(mergedData, Prior.COVID.infection. == 'No' & timeCategory == "Post 1st dose" & Tube == "HEP"), name1 = "Naive", 
           data2 = subset(mergedData, Prior.COVID.infection. == 'Yes' & timeCategory == "Post 1st dose" & Tube == "HEP"), name2 = "Experienced", 
           xData = "Age", yData = "CD8_.CD38.Ki67._FreqParent", fillParam = "Prior.COVID.infection.", title = "Vaccine dose 2", 
           xLabel = "Fold-change CD4+Ki67+CD38+", yLabel = "Fold-change CD8+Ki67+CD38+", nonparam = T) + 
  scale_x_continuous(limits = c(20,70), breaks=seq(0,100,10)) + scale_y_continuous(limits = c(0,5), breaks=seq(0,10,1))
# ggsave(filename = "./Images/FCActivCD4_correl_FCActivCD8_Vax2.pdf", width=8)

bivScatter(data1 = subset(mergedData, Prior.COVID.infection. == 'No' & timeCategory == "Post 2nd dose" & Tube == "HEP"), name1 = "Naive", 
           data2 = subset(mergedData, Prior.COVID.infection. == 'Yes' & timeCategory == "Post 2nd dose" & Tube == "HEP"), name2 = "Experienced", 
           xData = "Age", yData = "CD8_.CD38.Ki67._FreqParent", fillParam = "Prior.COVID.infection.", title = "Vaccine dose 2", 
           xLabel = "Fold-change CD4+Ki67+CD38+", yLabel = "Fold-change CD8+Ki67+CD38+", nonparam = T) + 
  scale_x_continuous(limits = c(20,70), breaks=seq(0,100,10)) + scale_y_continuous(limits = c(0,5), breaks=seq(0,10,1))
# ggsave(filename = "./Images/FCActivCD4_correl_FCActivCD8_Vax2.pdf", width=8)



subsetData <- subset(mergedData, Prior.COVID.infection. == 'No' & timeCategory == "Post 1st dose")
subsetData <- subsetData[,grep( paste(c("CD38.Ki67", "^FCA", "Age" ), collapse = "|"), names(subsetData))]
cor.matrix <- cor(subsetData, method="kendall" , use="pairwise.complete.obs" )
cor.matrix.pmat <- ggcorrplot::cor_pmat(subsetData, method="kendall", use="pairwise.complete.obs"  )
cor.matrix <- as.data.frame(cor.matrix); cor.matrix$Labels <- row.names(cor.matrix); cor.matrix$Prior.COVID <- "No"

subsetData <- subset(mergedData, Prior.COVID.infection. == 'Yes' & timeCategory == "Post 1st dose")
subsetData <- subsetData[,grep( paste(c("CD38.Ki67", "^FCA" , "Age"), collapse = "|"), names(subsetData))]
cor.matrix2 <- cor(subsetData, method="kendall" , use="pairwise.complete.obs" )
cor.matrix2 <- as.data.frame(cor.matrix2); cor.matrix2$Labels <- row.names(cor.matrix2); cor.matrix2$Prior.COVID <- "Yes"

cor.matrix <- as.data.frame(rbind(cor.matrix, cor.matrix2)); cor.matrix <- cor.matrix[ -grep( paste( c("Age", "Foxp3"), collapse = "|"), rownames(cor.matrix)),]

ggplot( data = cor.matrix, aes(y = Labels,x = Age, fill = Prior.COVID)) + geom_bar(stat='identity',position = 'dodge',width=0.75) + theme_bw() + 
  scale_fill_manual(values=c("#FFDFB1", "#B5B2F1")) + xlab("Correlation with Age") + ylab(" ") + theme(axis.text.y = element_text(angle=0, size = 10)) + 
  ggtitle("Post 1st dose") + geom_vline(xintercept=0, linetype = "dashed") + scale_x_continuous(limits = c(-1,1))
  
subsetData <- subset(mergedData, Prior.COVID.infection. == 'No' & timeCategory == "Post 2nd dose")
subsetData <- subsetData[,grep( paste(c("CD38.Ki67", "^FCA", "Age" ), collapse = "|"), names(subsetData))]
cor.matrix <- cor(subsetData, method="kendall" , use="pairwise.complete.obs" )
cor.matrix.pmat <- ggcorrplot::cor_pmat(subsetData, method="kendall", use="pairwise.complete.obs"  )
cor.matrix <- as.data.frame(cor.matrix); cor.matrix$Labels <- row.names(cor.matrix); cor.matrix$Prior.COVID <- "No"

subsetData <- subset(mergedData, Prior.COVID.infection. == 'Yes' & timeCategory == "Post 2nd dose")
subsetData <- subsetData[,grep( paste(c("CD38.Ki67", "^FCA" , "Age"), collapse = "|"), names(subsetData))]
cor.matrix2 <- cor(subsetData, method="kendall" , use="pairwise.complete.obs" )
cor.matrix2 <- as.data.frame(cor.matrix2); cor.matrix2$Labels <- row.names(cor.matrix2); cor.matrix2$Prior.COVID <- "Yes"

cor.matrix <- as.data.frame(rbind(cor.matrix, cor.matrix2)); cor.matrix <- cor.matrix[ -grep( paste( c("Age", "Foxp3"), collapse = "|"), rownames(cor.matrix)),]

ggplot( data = cor.matrix, aes(y = Labels,x = Age, fill = Prior.COVID)) + geom_bar(stat='identity',position = 'dodge',width=0.75) + theme_bw() + 
  scale_fill_manual(values=c("#FFDFB1", "#B5B2F1")) + xlab("Correlation with Age") + ylab(" ") + theme(axis.text.y = element_text(angle=0, size = 10)) + 
  ggtitle("Post 2nd dose") + geom_vline(xintercept=0, linetype = "dashed") +  scale_x_continuous(limits = c(-1,1))





subsetData <- subset(mergedData, Prior.COVID.infection. == 'Yes' & timeCategory == "Post 1st dose")
subsetData <- subsetData[,grep( paste(c("CD38.Ki67", "^FCA" , "Age"), collapse = "|"), names(subsetData))]
cor.matrix <- cor(subsetData, method="kendall" , use="pairwise.complete.obs" )
cor.matrix.pmat <- ggcorrplot::cor_pmat(subsetData, method="kendall", use="pairwise.complete.obs"  )










