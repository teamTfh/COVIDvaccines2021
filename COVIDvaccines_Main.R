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

mergedData <- readRDS(file = "mergedData.Rds")
# colors:   COVID-exp B5B2F1 (purple)             COVID-naive FFDFB1 (orange)
mergedData <- mergedData[-which(mergedData$Alias == "PHI-046"),]                    # only Moderna recipient, will temporarily exclude
#' ------------------ Cohort description --------------------------
#'

temp <- mergedData %>% group_by(Prior.COVID.infection., timeCategory) %>% get_summary_stats(type = 'common')
# write.csv(temp, file = "summaryStatistics.csv")

#' ------------------ Antibody analyses --------------------------
#'

# subsetData <- subset(mergedData, timeCategory != "two Weeks")
linePlot(data = mergedData, xData = 'timeCategory', yData = 'binding_IgG_S1', groupby = 'Alias', xLabel = 'Time after first vaccine dose', yLabel = "anti-S1 IgG titer", 
         title = "anti-S1 IgG titer", colorby = "Prior.COVID.infection.") + 
  scale_color_manual(name="Prior COVID?",values = c("#FFDFB1","#B5B2F1")) + 
  scale_y_continuous(trans='pseudo_log', limits = c(0,1e7), breaks=c(10^(0:7)), labels=trans_format('log10',math_format(10^.x)), minor_breaks =5*10^(0:10))# + ggrepel::geom_text_repel(aes(label = Label))
# ggsave(filename = "./Images/BindingAb_S1_IgG_linePlot.pdf")

linePlot(data = mergedData, xData = 'timeCategory', yData = 'IC50_neutAb_log10', groupby = 'Alias', xLabel = 'Time after first vaccine dose', yLabel = "log10 IC50", 
         title = "Neutralizing antibodies", colorby = "Prior.COVID.infection.") + 
  scale_color_manual(name="Prior COVID?",values = c("#FFDFB1","#B5B2F1")) + 
  scale_y_continuous(trans='pseudo_log', limits = c(0,1e5), breaks=c(10^(0:6)), labels=trans_format('log10',math_format(10^.x)), minor_breaks =5*10^(0:10))# + ggrepel::geom_text_repel(aes(label = Label))
# ggsave(filename = "./Images/neutAb_linePlot.pdf")

linePlot(data = mergedData, xData = 'timeCategory', yData = 'binding_IgG_N', groupby = 'Alias', xLabel = 'Time after first vaccine dose', yLabel = "anti-N IgG titer", 
         title = "anti-N IgG titer", colorby = "Prior.COVID.infection.") + 
  scale_color_manual(name="Prior COVID?",values = c("#FFDFB1","#B5B2F1")) + #theme(legend.position = 'none') + 
  scale_y_continuous(trans='pseudo_log', limits = c(0,1e5), breaks=c(10^(0:6)), labels=trans_format('log10',math_format(10^.x)), minor_breaks =5*10^(0:10)) #+ ggrepel::geom_text_repel(aes(label = Label))
# ggsave(filename = "./Images/BindingAb_N_IgG_linePlot.pdf", width=7)


linePlot(data = mergedData, xData = 'DPV', yData = 'binding_IgG_S1', groupby = 'Alias', xLabel = 'Time after first vaccine dose', yLabel = "anti-S1 IgG titer", 
         title = "anti-S1 IgG titer", colorby = "Prior.COVID.infection.") + 
  scale_color_manual(name="Prior COVID?",values = c("#FFDFB1","#B5B2F1")) + 
  scale_y_continuous(trans='pseudo_log', limits = c(0,1e7), breaks=c(10^(0:7)), labels=trans_format('log10',math_format(10^.x)), minor_breaks =5*10^(0:10))# + ggrepel::geom_text_repel(aes(label = Label))


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





#' ------------------ Activated T cells analyses --------------------------
#'
subsetData <- subset(mergedData, Prior.COVID.infection. == "No" & DPV>15)
a <- subsetData[,c("Alias", "DPV","binding_IgG_S1")]
# temp <- 
b<-dcast( a, `DPV` ~`Alias` , value.var = c("binding_IgG_S1")) 
c<-melt(data = b, id.vars = c("DPV"), measure.vars = c("HV-001","HV-002","HV-003","HV-008", "HV-018", "HV-052", "HV-054", "HV-068", "HV-076", "HV-078", 
                                                       "HV-079", "HV-080", "HV-081", "HV-082", "HV-084", "HV-085", "PHI-010", "PHI-020", "PHI-021"))
c <- c[which(!is.na(c$value)),]
Labels <- as.character(unique(c$variable)); 
Labels <- Labels[-grep(paste(c('HV-008','HV-084'), collapse ="|"), Labels, value=F)]
d <- data.frame(Labels = Labels, DPV = c(21,22,16,20,21,21,21), binding_IgG_S1 = c(67656, 741655, 364629, 346771, 274699, 104608, 230657))
ggplot(data = d, aes(x=DPV, y=binding_IgG_S1, label=Labels, group=Labels)) + geom_point() + theme_bw()


subsetData <- mergedData[- which(is.na(mergedData$timeCategory) ),]
subsetData <- subset(subsetData,  timeCategory != "two Weeks");  subsetData$timeCategory <- factor(subsetData$timeCategory, levels = c("Baseline", "Post 1st dose", "Pre 2nd dose", "Post 2nd dose"))
prePostTime(subsetData, xData = "DPV", yData="CD4_.CD38.Ki67._FreqParent", fillParam = "Prior.COVID.infection.", groupby="Record.ID", title = "Activated CD4",
            xLabel = "Prior COVID?", yLabel = "Ki67+CD38+ (% CD4)", repMeasures = F, exponential=F) +  coord_cartesian(xlim = c(-5,40))

prePostTime(subsetData, xData = "timeCategory", yData="CD4_.CD38.Ki67._FreqParent", fillParam = "Prior.COVID.infection.", groupby="Record.ID", title = "Activated CD4", 
            xLabel = "Prior COVID?", yLabel = "Ki67+CD38+ (% CD4)", repMeasures = F, exponential=F) #+ #+ scale_y_continuous(limits = c(0,5))
#  ggformula::geom_spline(aes_string(x="timeCategory", y="CD4_.CD38.Ki67._FreqParent"), color="blue",size=1, spar=0.5) #+ 
ggrepel::geom_text_repel(aes(label=Alias),size=3)
# ggsave(filename = "./Images/ActivCD4_bothCohorts_overTime.pdf")

prePostTime(subsetData, xData = "timeCategory", yData="CD8_..CD38.Ki67._FreqParent", fillParam = "Prior.COVID.infection.", groupby="Record.ID", title = "Activated CD8", 
            xLabel = "Prior COVID?", yLabel = "Ki67+CD38+ (% CD8)", repMeasures = F, exponential=F) + scale_y_continuous(breaks=seq(0,10,1)) #+ # limits = c(0,6))
#  ggformula::geom_spline(aes_string(x="timeCategory", y="CD4_.CD38.Ki67._FreqParent"), color="blue",size=1, spar=0.5) #+ 
ggrepel::geom_text_repel(aes(label=Alias),size=3)
# ggsave(filename = "./Images/ActivCD8_bothCohorts_overTime.pdf")



subsetData <- subset(mergedData, Record.ID != "CV-011" & Record.ID != "CV-012" )        # extreme outliers because of absence of Ki67 stain
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
prePostTime(data=subsetData, xData = "shortForm", yData="CD8_..CD38.Ki67._FreqParent", fillParam = "Prior.COVID.infection.", groupby="Record.ID", title = "Activated CD8", 
            xLabel = " ", yLabel = "Ki67+CD38+ (% CD8)") + 
  scale_y_continuous(breaks = seq(0,20,1), limits = c(0,6)) + scale_x_discrete(labels= c("Baseline","Post\n1st dose"))
# ggsave(filename = "./Images/ActivatedCD8_vax1_healthy_paired.pdf", width=4)


subsetData = subset(mergedData, mergedData$shortForm == '3W' | mergedData$shortForm=='4W'); subsetData = subset(subsetData, subsetData$Prior.COVID.infection.=="No")
prePostTime(data=subsetData, xData = "shortForm", yData="CD4_.CD38.Ki67._FreqParent", fillParam = "Prior.COVID.infection.", groupby="Record.ID", title = "Activated CD4", 
            xLabel = " ", yLabel = "Ki67+CD38+ (% CD4)") + scale_y_continuous(breaks = seq(0,20,1), limits = c(0,5)) + scale_x_discrete(labels= c("Pre\n2nd dose","Post\n2nd dose"))
# ggsave(filename = "./Images/ActivatedCD4_vax2_healthy_paired.pdf", width=4)

subsetData = subset(mergedData, mergedData$shortForm == '3W' | mergedData$shortForm=='4W'); subsetData = subset(subsetData, subsetData$Prior.COVID.infection.=="No")
prePostTime(data=subsetData, xData = "shortForm", yData="CD8_..CD38.Ki67._FreqParent", fillParam = "Prior.COVID.infection.", groupby="Record.ID", title = "Activated CD8", 
            xLabel = " ", yLabel = "Ki67+CD38+ (% CD8)") + scale_y_continuous(breaks = seq(0,20,1), limits = c(0,6)) + scale_x_discrete(labels= c("Pre\n2nd dose","Post\n2nd dose"))
# ggsave(filename = "./Images/ActivatedCD8_vax2_healthy_paired.pdf", width=4)

# ------------------------------------------------ PHI side --------------------------------------------------------
subsetData = subset(mergedData, mergedData$shortForm == 'bL' | mergedData$shortForm=='oW'); subsetData = subset(subsetData, subsetData$Prior.COVID.infection.=="Yes")
prePostTime(data=subsetData, xData = "shortForm", yData="CD4_.CD38.Ki67._FreqParent", fillParam = "Prior.COVID.infection.", groupby="Record.ID", title = "Activated CD4", 
            xLabel = " ", yLabel = "Ki67+CD38+ (% CD4)") + scale_y_continuous(breaks = seq(0,20,1), limits = c(0,5)) + 
  scale_x_discrete(labels= c("Baseline","Post\n1st dose"))+ scale_fill_manual(values=c("#B5B2F1"))# + ggrepel::geom_text_repel(aes(label=Alias))
# ggsave(filename = "./Images/ActivatedCD4_vax1_covid_paired.pdf", width=4)

subsetData = subset(mergedData, mergedData$shortForm == 'bL' | mergedData$shortForm=='oW'); subsetData = subset(subsetData, subsetData$Prior.COVID.infection.=="Yes")
prePostTime(data=subsetData, xData = "shortForm", yData="CD8_..CD38.Ki67._FreqParent", fillParam = "Prior.COVID.infection.", groupby="Record.ID", title = "Activated CD8", 
            xLabel = " ", yLabel = "Ki67+CD38+ (% CD8)") + 
  scale_y_continuous(breaks = seq(0,20,1), limits = c(0,6)) + scale_x_discrete(labels= c("Baseline","Post\n1st dose"))+ scale_fill_manual(values=c("#B5B2F1")) 
# ggsave(filename = "./Images/ActivatedCD8_vax1_covid_paired.pdf", width=4)


subsetData = subset(mergedData, mergedData$shortForm == '3W' | mergedData$shortForm=='4W'); subsetData = subset(subsetData, subsetData$Prior.COVID.infection.=="Yes")
prePostTime(data=subsetData, xData = "shortForm", yData="CD4_.CD38.Ki67._FreqParent", fillParam = "Prior.COVID.infection.", groupby="Record.ID", title = "Activated CD4", 
            xLabel = " ", yLabel = "Ki67+CD38+ (% CD4)") + scale_y_continuous(breaks = seq(0,20,1), limits = c(0,5)) + 
  scale_x_discrete(labels= c("Pre\n2nd dose","Post\n2nd dose"))+ scale_fill_manual(values=c("#B5B2F1")) 
# ggsave(filename = "./Images/ActivatedCD4_vax2_covid_paired.pdf", width=4)

subsetData = subset(mergedData, mergedData$shortForm == '3W' | mergedData$shortForm=='4W'); subsetData = subset(subsetData, subsetData$Prior.COVID.infection.=="Yes")
prePostTime(data=subsetData, xData = "shortForm", yData="CD8_..CD38.Ki67._FreqParent", fillParam = "Prior.COVID.infection.", groupby="Record.ID", title = "Activated CD8", 
            xLabel = " ", yLabel = "Ki67+CD38+ (% CD8)") + scale_y_continuous(breaks = seq(0,20,1), limits = c(0,6)) + 
  scale_x_discrete(labels= c("Pre\n2nd dose","Post\n2nd dose")) + scale_fill_manual(values=c("#B5B2F1")) 
# ggsave(filename = "./Images/ActivatedCD8_vax2_covid_paired.pdf", width=4)


#' ------------------ Tfh analyses --------------------------
#' 
subsetData <- mergedData[which(mergedData$shortForm == "bL" | mergedData$shortForm == "oW"),]
prePostTime(data = subsetData, xData = "shortForm", yData="CD4_.Nonnaive.cTfh.ICOS..CD38.._FreqParent", fillParam = "Prior.COVID.infection.", groupby="Record.ID", 
            title = "cTfh", xLabel = "Prior COVID?", yLabel = "ICOS+CD38+ (% cTfh)") 
  
subsetData <- mergedData[- which(is.na(mergedData$timeCategory) ),]
subsetData <- subset(subsetData,  timeCategory != "two Weeks");  subsetData$timeCategory <- factor(subsetData$timeCategory, levels = c("Baseline", "Post 1st dose", "Pre 2nd dose", "Post 2nd dose"))
prePostTime(subsetData, xData = "timeCategory", yData="CD4_.Nonnaive.cTfh.ICOS..CD38.._FreqParent", fillParam = "Prior.COVID.infection.", groupby="Record.ID", title = "cTfh", 
            xLabel = "Prior COVID?", yLabel = "ICOS+CD38+ (% cTfh)", repMeasures = F) #  + ggrepel::geom_text_repel(aes(label=Record.ID),size=2)
# ggsave(filename = "./Images/cTfh_responses_bothCohorts_overTime.pdf")
  
  
twoSampleBar(data = subset(mergedData, timeCategory == "Post 1st dose"), xData = "Prior.COVID.infection.", yData = "FCtfh_Vax1", fillParam = "Prior.COVID.infection.", title = "ICOS+CD38+ cTfh", 
             yLabel = "Fold-change at one week")


subsetData = subset(mergedData, mergedData$shortForm == 'bL' | mergedData$shortForm=='oW'); subsetData = subset(subsetData, subsetData$Prior.COVID.infection.=="No")
prePostTime(data=subsetData, xData = "shortForm", yData="CD4_.Nonnaive.cTfh.ICOS..CD38.._FreqParent", fillParam = "Prior.COVID.infection.", groupby="Record.ID", title = "cTfh - Vax1", 
            xLabel = " ", yLabel = "ICOS+CD38+ (% cTfh)") + scale_y_continuous(breaks = seq(0,21,3), limits = c(0,20))  + 
  scale_x_discrete(labels= c("Baseline","Post\n1st dose"))  # +   ggrepel::geom_text_repel(aes(label=Record.ID))
# ggsave(filename = "./Images/Tfh_vax1_healthy_paired.pdf", width=4)
subsetData = subset(mergedData, mergedData$shortForm == '3W' | mergedData$shortForm=='4W'); subsetData = subset(subsetData, subsetData$Prior.COVID.infection.=="No")
prePostTime(data=subsetData, xData = "shortForm", yData="CD4_.Nonnaive.cTfh.ICOS..CD38.._FreqParent", fillParam = "Prior.COVID.infection.", groupby="Record.ID", 
            title = "cTfh - Vax2", xLabel = " ", yLabel = "ICOS+CD38+ (% cTfh)") + 
  scale_y_continuous(breaks = seq(0,21,3), limits = c(0,20))  + scale_x_discrete(labels= c("Baseline","Post\n1st dose"))
# ggsave(filename = "./Images/Tfh_vax2_healthy_paired.pdf", width=4)


subsetData = subset(mergedData, mergedData$shortForm == 'bL' | mergedData$shortForm=='oW'); subsetData = subset(subsetData, subsetData$Prior.COVID.infection.=="Yes")
prePostTime(data=subsetData, xData = "shortForm", yData="CD4_.Nonnaive.cTfh.ICOS..CD38.._FreqParent", fillParam = "Prior.COVID.infection.", groupby="Record.ID", 
            title = "cTfh - Vax1", xLabel = " ", yLabel = "ICOS+CD38+ (% cTfh)") + scale_x_discrete(labels= c("Baseline","Post\n1st dose")) + scale_fill_manual(values=c("#B5B2F1")) +
  scale_y_continuous(breaks = seq(0,20,3), limits = c(0,20)) 
# ggsave(filename = "./Images/Tfh_vax1_covid_paired.pdf", width=4)

subsetData = subset(mergedData, mergedData$shortForm == '3W' | mergedData$shortForm=='4W'); subsetData = subset(subsetData, subsetData$Prior.COVID.infection.=="Yes")
prePostTime(data=subsetData, xData = "shortForm", yData="CD4_.Nonnaive.cTfh.ICOS..CD38.._FreqParent", fillParam = "Prior.COVID.infection.", groupby="Record.ID", 
            title = "cTfh - Vax2", xLabel = " ", yLabel = "ICOS+CD38+ (% cTfh)") + scale_x_discrete(labels= c("Baseline","Post\n1st dose"))+ scale_fill_manual(values=c("#B5B2F1")) +
  scale_y_continuous(breaks = seq(0,20,3), limits = c(0,20)) 
# ggsave(filename = "./Images/Tfh_vax2_covid_paired.pdf", width=4)



subsetData %>% group_by(shortForm) %>% get_summary_stats(type = "mean") %>% print(n=80)


subsetData = subset(mergedData, mergedData$shortForm == 'oW')
subsetData %>% group_by(Prior.COVID.infection.) %>% shapiro_test(CD4_.Nonnaive.cTfh.ICOS..CD38.._FreqParent)
twoSampleBar(data=subsetData, xData = "Prior.COVID.infection.", yData="CD4_.Nonnaive.cTfh.ICOS..CD38.._FreqParent", fillParam = "Prior.COVID.infection.", title = "cTfh at oneWeek", 
            yLabel = "ICOS+CD38+ (% cTfh)", nonparam = T) + scale_y_continuous(breaks = seq(0,20,1), limits = c(0,16)) 

#' ------------------ B cell analyses --------------------------
#'


prePostTime(data=mergedData, xData = "day", yData="CD19_..Nonnaive.B.CD27..CD38._FreqParent", fillParam = "Prior.COVID.infection.", groupby="Record.ID", title = "Plasmablast", 
            xLabel = "Prior COVID?", yLabel = "PB (% CD19)") + 
  ggformula::geom_spline(aes_string(x="day", y="CD19_..Nonnaive.B.CD27..CD38._FreqParent"), color="blue",size=1, spar=0.5) # + ggrepel::geom_text_repel(aes(label=Record.ID),size=2)


subsetData = subset(mergedData, mergedData$shortForm == 'bL' | mergedData$shortForm=='oW'); subsetData = subset(subsetData, subsetData$Prior.COVID.infection.=="No")
prePostTime(data=subsetData, xData = "shortForm", yData="CD19_..Nonnaive.B.CD27..CD38._FreqParent", fillParam = "Prior.COVID.infection.", groupby="Record.ID", title = "PB - Vax1", 
            xLabel = " ", yLabel = "CD27+CD38+ (% nonnavB)") + scale_y_continuous(breaks = seq(0,21,0.5), limits = c(0,2))  + 
  scale_x_discrete(labels= c("Baseline","One Week"))# + ggrepel::geom_text_repel(aes(label = Record.ID))
# ggsave(filename = "./Images/PB_vax1_healthy_paired.pdf", width=4)


subsetData = subset(mergedData, mergedData$shortForm == '3W' | mergedData$shortForm=='4W'); subsetData = subset(subsetData, subsetData$Prior.COVID.infection.=="No")
prePostTime(data=subsetData, xData = "shortForm", yData="CD19_..Nonnaive.B.CD27..CD38._FreqParent", fillParam = "Prior.COVID.infection.", groupby="Record.ID", 
            title = "PB - Vax2", xLabel = " ", yLabel = "CD27+CD38+ (% nonnavB)") + scale_y_continuous(breaks = seq(0,21,0.5), limits = c(0,2))  + 
  scale_x_discrete(labels= c("Baseline","One Week"))# + ggrepel::geom_text_repel(aes(label = Record.ID))
# ggsave(filename = "./Images/PB_vax2_healthy_paired.pdf", width=4)



subsetData = subset(mergedData, mergedData$shortForm == 'bL' | mergedData$shortForm=='oW'); subsetData = subset(subsetData, subsetData$Prior.COVID.infection.=="Yes")
prePostTime(data=subsetData, xData = "shortForm", yData="CD19_..Nonnaive.B.CD27..CD38._FreqParent", fillParam = "Prior.COVID.infection.", groupby="Record.ID", 
            title = "PB - Vax1", xLabel = " ", yLabel = "CD27+CD38+ (% nonnavB)") +    
  scale_y_continuous(breaks = seq(0,21,0.5), limits = c(0,3))  + scale_x_discrete(labels= c("Baseline","One Week"))+ scale_fill_manual(values=c("#B5B2F1")) 
# ggsave(filename = "./Images/PB_vax1_covid_paired.pdf", width=4)

subsetData = subset(mergedData, mergedData$shortForm == '3W' | mergedData$shortForm=='4W'); subsetData = subset(subsetData, subsetData$Prior.COVID.infection.=="Yes")
prePostTime(data=subsetData, xData = "shortForm", yData="CD19_..Nonnaive.B.CD27..CD38._FreqParent", fillParam = "Prior.COVID.infection.", groupby="Record.ID", 
            title = "PB - Vax2", xLabel = " ", yLabel = "CD27+CD38+ (% nonnavB)") +    
  scale_y_continuous(breaks = seq(0,21,0.5), limits = c(0,3))  + scale_x_discrete(labels= c("Baseline","One Week"))+ scale_fill_manual(values=c("#B5B2F1")) 
# ggsave(filename = "./Images/PB_vax2_covid_paired.pdf", width=4)

subsetData = subset(mergedData, mergedData$shortForm == 'bL' ); subsetData = subset(subsetData, subsetData$Prior.COVID.infection.=="No")
univScatter(data = subsetData, xData = "FCtfh_Vax1", yData = "FCPB_Vax1", fillParam = "Prior.COVID.infection.", title = " ", xLabel = "Fold-change Tfh", 
            yLabel = "Fold-change PB", position = "left")
subsetData = subset(mergedData, mergedData$shortForm == 'bL' ); subsetData = subset(subsetData, subsetData$Prior.COVID.infection.=="Yes")
univScatter(data = subsetData, xData = "FCtfh_Vax1", yData = "FCPB_Vax1", fillParam = "Prior.COVID.infection.", title = " ", xLabel = "Fold-change Tfh", 
            yLabel = "Fold-change PB", position = "left")+ scale_fill_manual(values=c("#B5B2F1")) 

subsetData <- subset(mergedData, mergedData$Record.ID != "CV-022")        # baseline for CD11c looks super strange
prePostTime(subsetData, xData = "day", yData="CD19_..Nonnaive.B.Tbet..CD11C._FreqParent", fillParam = "Prior.COVID.infection.", groupby="Record.ID", title = "Tbet+CD11c+", 
            xLabel = "Prior COVID?", yLabel = "Tbet+CD11+ (% nonnavB)") + 
  ggformula::geom_spline(aes_string(x="day", y="CD19_..Nonnaive.B.Tbet..CD11C._FreqParent"), color="blue",size=1, spar=0.5)

  







