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
library("gridExtra")
source('./COVIDvaccines_PlottingFunctions.R')

mergedData <- readRDS(file = "mergedData.Rds")
# colors:   COVID-exp B5B2F1 (purple)             COVID-naive FFDFB1 (orange)

#' ------------------ Antibody analyses --------------------------
#'

# subsetData <- subset(mergedData, timeCategory != "two Weeks")
linePlot(data = mergedData, xData = 'timeCategory', yData = 'binding_IgG_S1', groupby = 'Alias', xLabel = 'Time after first vaccine dose', yLabel = "anti-S1 IgG titer", 
         title = "anti-S1 IgG titer", colorby = "Prior.COVID.infection.") + 
  scale_color_manual(name="Prior COVID?",values = c("#FFDFB1","#B5B2F1")) + scale_y_continuous(trans='log10')  #+ geom_label(aes(label=Label),size=2) #ggrepel::geom_label_repel(aes(label = Label))
ggsave(filename = "./Images/BindingAb_S1_IgG_linePlot.pdf")
linePlot(data = mergedData, xData = 'timeCategory', yData = 'IC50_neutAb_log10', groupby = 'Alias', xLabel = 'Time after first vaccine dose', yLabel = "log10 IC50", 
         title = "Neutralizing antibodies", colorby = "Prior.COVID.infection.") + 
  scale_color_manual(name="Prior COVID?",values = c("#FFDFB1","#B5B2F1")) + scale_y_continuous(trans='log10') # + geom_label(aes(label=Label),size=2) #ggrepel::geom_label_repel(aes(label = Label))
ggsave(filename = "./Images/neutAb_linePlot.pdf")

#' ------------------ Elispots analyses --------------------------
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
ggpubr::ggexport(plotlist = list(a.1, a.2, a.3), filename = "./Images/Elispots_IgG_post1stDose_gridarrange.pdf", nrow = 1, width = 12)

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
ggpubr::ggexport(plotlist = list(a.1, a.2, a.3), filename = "./Images/Elispots_IgA_post1stDose_gridarrange.pdf", nrow = 1, width = 12 )

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
ggpubr::ggexport(plotlist = list(a.1, a.2, a.3), filename = "./Images/Elispots_IgM_post1stDose_gridarrange.pdf", nrow = 1, width = 12 )

subsetData <- subset(mergedData, !is.na(Elispot_IgG_S1) & timeCategory != "Pre 2nd dose" )
pdf(file = "./Images/Elispots_IgG_S1_prepostTime.pdf")
  prePostTime.expon(data = subsetData, xData = "timeCategory", yData = "Elispot_IgG_S1", fillParam = "Prior.COVID.infection.", 
              groupby = "Alias", title = "ASC IgG S1", xLabel = " ", yLabel = "Spots per 1e6 PBMC")+ 
    scale_y_continuous(trans='pseudo_log', limits = c(0,10000), breaks=c(10^(0:4)), labels=trans_format('log10',math_format(10^.x)), minor_breaks =5*10^(0:10))
dev.off()

subsetData <- subset(mergedData, !is.na(Elispot_IgG_S2) & timeCategory != "Pre 2nd dose" )
pdf(file = "./Images/Elispots_IgG_S2_prepostTime.pdf")
prePostTime.expon(data = subsetData, xData = "timeCategory", yData = "Elispot_IgG_S2", fillParam = "Prior.COVID.infection.", 
                  groupby = "Alias", title = "ASC IgG S2", xLabel = " ", yLabel = "Spots per 1e6 PBMC")+ 
  scale_y_continuous(trans='pseudo_log', limits = c(0,10000), breaks=c(10^(0:4)), labels=trans_format('log10',math_format(10^.x)), minor_breaks =5*10^(0:10))
dev.off()

subsetData <- subset(mergedData, !is.na(Elispot_IgG_RBD) & timeCategory != "Pre 2nd dose" )
pdf(file = "./Images/Elispots_IgG_RBD_prepostTime.pdf")
prePostTime.expon(data = subsetData, xData = "timeCategory", yData = "Elispot_IgG_RBD", fillParam = "Prior.COVID.infection.", 
                  groupby = "Alias", title = "ASC IgG RBD", xLabel = " ", yLabel = "Spots per 1e6 PBMC")+ 
  scale_y_continuous(trans='pseudo_log', limits = c(0,10000), breaks=c(10^(0:4)), labels=trans_format('log10',math_format(10^.x)), minor_breaks =5*10^(0:10))
dev.off()


#' ------------------ Activated T cells analyses --------------------------
#'

prePostTime(mergedData, xData = "DPV", yData="CD4_.CD38.Ki67._FreqParent", fillParam = "Prior.COVID.infection.", groupby="Record.ID", title = "Activated CD4", 
            xLabel = "Prior COVID?", yLabel = "Ki67+CD38+ (% CD4)") + 
  ggformula::geom_spline(aes_string(x="DPV", y="CD4_.CD38.Ki67._FreqParent"), color="blue",size=1, spar=0.5) #+ ggrepel::geom_text_repel(aes(label=Record.ID),size=2)


subsetData <- subset(mergedData, Record.ID != "CV-011" & Record.ID != "CV-012" )
twoSampleBar(data = subset(subsetData, timeCategory == "oneWeek"), xData = "Prior.COVID.infection.", yData = "FCActivCD4_Vax1", fillParam = "Prior.COVID.infection.", title = "Activ CD4", 
             yLabel = "Fold-change at one week", nonparam = T)


subsetData = subset(mergedData, mergedData$shortForm == 'bL' | mergedData$shortForm=='oW'); subsetData = subset(subsetData, subsetData$Prior.COVID.infection.=="No")
subsetData <- subset(subsetData, Record.ID != "CV-011" & Record.ID != "CV-012" )
prePostTime(data=subsetData, xData = "shortForm", yData="CD4_.CD38.Ki67._FreqParent", fillParam = "Prior.COVID.infection.", groupby="Record.ID", title = "Activated CD4", 
            xLabel = " ", yLabel = "Ki67+CD38+ (% CD4)") + scale_y_continuous(breaks = seq(0,20,0.5), limits = c(0,3)) + scale_x_discrete(labels= c("Baseline","One Week"))
# ggsave(filename = "./Images/ActivatedCD4_vax1_healthy_paired.pdf", width=4)

subsetData = subset(mergedData, mergedData$shortForm == 'bL' | mergedData$shortForm=='oW'); subsetData = subset(subsetData, subsetData$Prior.COVID.infection.=="No")
subsetData <- subset(subsetData, Record.ID != "CV-011" & Record.ID != "CV-012" )
prePostTime(data=subsetData, xData = "shortForm", yData="CD8_.CD38.Ki67._FreqParent", fillParam = "Prior.COVID.infection.", groupby="Record.ID", title = "Activated CD8", 
            xLabel = " ", yLabel = "Ki67+CD38+ (% CD8)") + 
  scale_y_continuous(breaks = seq(0,20,0.5), limits = c(0,3)) + scale_x_discrete(labels= c("Baseline","One Week"))
# ggsave(filename = "./Images/ActivatedCD8_vax1_healthy_paired.pdf", width=4)


subsetData = subset(mergedData, mergedData$shortForm == '3W' | mergedData$shortForm=='4W'); subsetData = subset(subsetData, subsetData$Prior.COVID.infection.=="No")
subsetData <- subset(subsetData, Record.ID != "CV-011" & Record.ID != "CV-012" )
prePostTime(data=subsetData, xData = "shortForm", yData="CD4_.CD38.Ki67._FreqParent", fillParam = "Prior.COVID.infection.", groupby="Record.ID", title = "Activated CD4", 
            xLabel = " ", yLabel = "Ki67+CD38+ (% CD4)") + scale_y_continuous(breaks = seq(0,20,0.5), limits = c(0,3)) + scale_x_discrete(labels= c("Baseline","One Week"))
# ggsave(filename = "./Images/ActivatedCD4_vax2_healthy_paired.pdf", width=4)

subsetData = subset(mergedData, mergedData$shortForm == '3W' | mergedData$shortForm=='4W'); subsetData = subset(subsetData, subsetData$Prior.COVID.infection.=="No")
subsetData <- subset(subsetData, Record.ID != "CV-011" & Record.ID != "CV-012" )
prePostTime(data=subsetData, xData = "shortForm", yData="CD8_.CD38.Ki67._FreqParent", fillParam = "Prior.COVID.infection.", groupby="Record.ID", title = "Activated CD8", 
            xLabel = " ", yLabel = "Ki67+CD38+ (% CD8)") + 
  scale_y_continuous(breaks = seq(0,20,0.5), limits = c(0,3)) + scale_x_discrete(labels= c("Baseline","One Week"))
# ggsave(filename = "./Images/ActivatedCD8_vax2_healthy_paired.pdf", width=4)


subsetData = subset(mergedData, mergedData$shortForm == 'bL' | mergedData$shortForm=='oW'); subsetData = subset(subsetData, subsetData$Prior.COVID.infection.=="Yes")
subsetData <- subset(subsetData, Record.ID != "CV-011" & Record.ID != "CV-012" )
prePostTime(data=subsetData, xData = "shortForm", yData="CD4_.CD38.Ki67._FreqParent", fillParam = "Prior.COVID.infection.", groupby="Record.ID", title = "Activated CD4", 
            xLabel = " ", yLabel = "Ki67+CD38+ (% CD4)") + # ggrepel::geom_text_repel(aes(label=Record.ID),size=2) + 
  scale_y_continuous(breaks = seq(0,20,0.5), limits = c(0,3)) + scale_x_discrete(labels= c("Baseline","One Week"))
# ggsave(filename = "./Images/ActivatedCD4_vax1_covid_paired.pdf", width=4)                                                                             

subsetData = subset(mergedData, mergedData$shortForm == 'bL' | mergedData$shortForm=='oW'); subsetData = subset(subsetData, subsetData$Prior.COVID.infection.=="Yes")
subsetData <- subset(subsetData, Record.ID != "CV-011" & Record.ID != "CV-012" )
prePostTime(data=subsetData, xData = "shortForm", yData="CD8_.CD38.Ki67._FreqParent", fillParam = "Prior.COVID.infection.", groupby="Record.ID", title = "Activated CD8", 
            xLabel = " ", yLabel = "Ki67+CD38+ (% CD8)") + # ggrepel::geom_text_repel(aes(label=Record.ID),size=2) 
  scale_y_continuous(breaks = seq(0,20,0.5), limits = c(0,3)) + scale_x_discrete(labels= c("Baseline","One Week"))
# ggsave(filename = "./Images/ActivatedCD8_vax1_covid_paired.pdf", width=4)


subsetData = subset(mergedData, mergedData$shortForm == '3W' | mergedData$shortForm=='4W'); subsetData = subset(subsetData, subsetData$Prior.COVID.infection.=="Yes")
subsetData <- subset(subsetData, Record.ID != "CV-011" & Record.ID != "CV-012" )
prePostTime(data=subsetData, xData = "shortForm", yData="CD4_.CD38.Ki67._FreqParent", fillParam = "Prior.COVID.infection.", groupby="Record.ID", title = "Activated CD4", 
            xLabel = " ", yLabel = "Ki67+CD38+ (% CD4)")  + scale_y_continuous(breaks = seq(0,20,0.5), limits = c(0,3)) + scale_x_discrete(labels= c("Baseline","One Week"))
# ggsave(filename = "./Images/ActivatedCD4_vax2_covid_paired.pdf", width=4)

subsetData = subset(mergedData, mergedData$shortForm == '3W' | mergedData$shortForm=='4W'); subsetData = subset(subsetData, subsetData$Prior.COVID.infection.=="Yes")
subsetData <- subset(subsetData, Record.ID != "CV-011" & Record.ID != "CV-012" )
prePostTime(data=subsetData, xData = "shortForm", yData="CD8_.CD38.Ki67._FreqParent", fillParam = "Prior.COVID.infection.", groupby="Record.ID", title = "Activated CD8", 
            xLabel = " ", yLabel = "Ki67+CD38+ (% CD8)") +  
  scale_y_continuous(breaks = seq(0,20,0.5), limits = c(0,6)) + scale_x_discrete(labels= c("Baseline","One Week"))
                                                      # ^ scale previously cut off several data points 
# ggsave(filename = "./Images/ActivatedCD8_vax2_covid_paired.pdf", width=4)


#' ------------------ Tfh analyses --------------------------
#' 
subsetData <- mergedData[which(mergedData$shortForm == "bL" | mergedData$shortForm == "oW"),]
prePostTime(data = subsetData, xData = "shortForm", yData="CD4_.Nonnaive.cTfh.ICOS..CD38.._FreqParent", fillParam = "Prior.COVID.infection.", groupby="Record.ID", title = "cTfh", 
            xLabel = "Prior COVID?", yLabel = "ICOS+CD38+ (% cTfh)") 
  
prePostTime(mergedData, xData = "day", yData="CD4_.Nonnaive.cTfh.ICOS..CD38.._FreqParent", fillParam = "Prior.COVID.infection.", groupby="Record.ID", title = "cTfh", 
            xLabel = "Prior COVID?", yLabel = "ICOS+CD38+ (% cTfh)") + 
  ggformula::geom_spline(aes_string(x="day", y="CD4_.Nonnaive.cTfh.ICOS..CD38.._FreqParent"), color="blue",size=1, spar=0.5) #  + ggrepel::geom_text_repel(aes(label=Record.ID),size=2)

twoSampleBar(data = subset(mergedData, timeCategory == "oneWeek"), xData = "Prior.COVID.infection.", yData = "FCtfh_Vax1", fillParam = "Prior.COVID.infection.", title = "ICOS+CD38+ cTfh", 
             yLabel = "Fold-change at one week")


subsetData = subset(mergedData, mergedData$shortForm == 'bL' | mergedData$shortForm=='oW'); subsetData = subset(subsetData, subsetData$Prior.COVID.infection.=="No")
prePostTime(data=subsetData, xData = "shortForm", yData="CD4_.Nonnaive.cTfh.ICOS..CD38.._FreqParent", fillParam = "Prior.COVID.infection.", groupby="Record.ID", title = "cTfh - Vax1", 
            xLabel = " ", yLabel = "ICOS+CD38+ (% cTfh)") + scale_y_continuous(breaks = seq(0,21,3), limits = c(0,20))  + 
  scale_x_discrete(labels= c("Baseline","One Week"))  # +   ggrepel::geom_text_repel(aes(label=Record.ID))
# ggsave(filename = "./Images/Tfh_vax1_healthy_paired.pdf", width=4)
subsetData = subset(mergedData, mergedData$shortForm == '3W' | mergedData$shortForm=='4W'); subsetData = subset(subsetData, subsetData$Prior.COVID.infection.=="No")
prePostTime(data=subsetData, xData = "shortForm", yData="CD4_.Nonnaive.cTfh.ICOS..CD38.._FreqParent", fillParam = "Prior.COVID.infection.", groupby="Record.ID", 
            title = "cTfh - Vax2", xLabel = " ", yLabel = "ICOS+CD38+ (% cTfh)") + 
  scale_y_continuous(breaks = seq(0,21,3), limits = c(0,20))  + scale_x_discrete(labels= c("Baseline","One Week"))
# ggsave(filename = "./Images/Tfh_vax2_healthy_paired.pdf", width=4)


subsetData = subset(mergedData, mergedData$shortForm == 'bL' | mergedData$shortForm=='oW'); subsetData = subset(subsetData, subsetData$Prior.COVID.infection.=="Yes")
prePostTime(data=subsetData, xData = "shortForm", yData="CD4_.Nonnaive.cTfh.ICOS..CD38.._FreqParent", fillParam = "Prior.COVID.infection.", groupby="Record.ID", 
            title = "cTfh - Vax1", xLabel = " ", yLabel = "ICOS+CD38+ (% cTfh)") + scale_x_discrete(labels= c("Baseline","One Week"))
#+ scale_y_continuous(breaks = seq(0,20,1), limits = c(0,20)) 
# ggsave(filename = "./Images/Tfh_vax1_covid_paired.pdf", width=4)

subsetData = subset(mergedData, mergedData$shortForm == '3W' | mergedData$shortForm=='4W'); subsetData = subset(subsetData, subsetData$Prior.COVID.infection.=="Yes")
prePostTime(data=subsetData, xData = "shortForm", yData="CD4_.Nonnaive.cTfh.ICOS..CD38.._FreqParent", fillParam = "Prior.COVID.infection.", groupby="Record.ID", 
            title = "cTfh - Vax2", xLabel = " ", yLabel = "ICOS+CD38+ (% cTfh)") + scale_x_discrete(labels= c("Baseline","One Week"))
#+ scale_y_continuous(breaks = seq(0,20,1), limits = c(0,20)) 
# ggsave(filename = "./Images/Tfh_vax2_covid_paired.pdf", width=4)



subsetData %>% group_by(shortForm) %>% get_summary_stats(type = "mean") %>% print(n=80)


subsetData = subset(mergedData, mergedData$shortForm == 'oW')
subsetData %>% group_by(Prior.COVID.infection.) %>% shapiro_test(CD4_.Nonnaive.cTfh.ICOS..CD38.._FreqParent)
twoSampleBar(data=subsetData, xData = "Prior.COVID.infection.", yData="CD4_.Nonnaive.cTfh.ICOS..CD38.._FreqParent", fillParam = "Prior.COVID.infection.", title = "cTfh at oneWeek", 
            yLabel = "ICOS+CD38+ (% cTfh)", nonparam = T) + scale_y_continuous(breaks = seq(0,20,1), limits = c(0,16)) 

#' ------------------ B cell analyses --------------------------
#'
prePostTime(data=mergedData, xData = "day", yData="CD19_.Nonnaive.B.CD27..CD38._FreqParent", fillParam = "Prior.COVID.infection.", groupby="Record.ID", title = "Plasmablast", 
            xLabel = "Prior COVID?", yLabel = "PB (% CD19)") + 
  ggformula::geom_spline(aes_string(x="day", y="CD19_.Nonnaive.B.CD27..CD38._FreqParent"), color="blue",size=1, spar=0.5) # + ggrepel::geom_text_repel(aes(label=Record.ID),size=2)


subsetData = subset(mergedData, mergedData$shortForm == 'bL' | mergedData$shortForm=='oW'); subsetData = subset(subsetData, subsetData$Prior.COVID.infection.=="No")
prePostTime(data=subsetData, xData = "shortForm", yData="CD19_.Nonnaive.B.CD27..CD38._FreqParent", fillParam = "Prior.COVID.infection.", groupby="Record.ID", title = "PB - Vax1", 
            xLabel = " ", yLabel = "CD27+CD38+ (% nonnavB)") + scale_y_continuous(breaks = seq(0,21,0.5), limits = c(0,2))  + 
  scale_x_discrete(labels= c("Baseline","One Week"))# + ggrepel::geom_text_repel(aes(label = Record.ID))
# ggsave(filename = "./Images/PB_vax1_healthy_paired.pdf", width=4)


subsetData = subset(mergedData, mergedData$shortForm == '3W' | mergedData$shortForm=='4W'); subsetData = subset(subsetData, subsetData$Prior.COVID.infection.=="No")
prePostTime(data=subsetData, xData = "shortForm", yData="CD19_.Nonnaive.B.CD27..CD38._FreqParent", fillParam = "Prior.COVID.infection.", groupby="Record.ID", 
            title = "PB - Vax2", xLabel = " ", yLabel = "CD27+CD38+ (% nonnavB)") + scale_y_continuous(breaks = seq(0,21,0.5), limits = c(0,2))  + 
  scale_x_discrete(labels= c("Baseline","One Week"))# + ggrepel::geom_text_repel(aes(label = Record.ID))
# ggsave(filename = "./Images/PB_vax2_healthy_paired.pdf", width=4)



subsetData = subset(mergedData, mergedData$shortForm == 'bL' | mergedData$shortForm=='oW'); subsetData = subset(subsetData, subsetData$Prior.COVID.infection.=="Yes")
prePostTime(data=subsetData, xData = "shortForm", yData="CD19_.Nonnaive.B.CD27..CD38._FreqParent", fillParam = "Prior.COVID.infection.", groupby="Record.ID", 
            title = "PB - Vax1", xLabel = " ", yLabel = "CD27+CD38+ (% nonnavB)") +    
  scale_y_continuous(breaks = seq(0,21,0.5), limits = c(0,3))  + scale_x_discrete(labels= c("Baseline","One Week"))
# ggsave(filename = "./Images/PB_vax1_covid_paired.pdf", width=4)

subsetData = subset(mergedData, mergedData$shortForm == '3W' | mergedData$shortForm=='4W'); subsetData = subset(subsetData, subsetData$Prior.COVID.infection.=="Yes")
prePostTime(data=subsetData, xData = "shortForm", yData="CD19_.Nonnaive.B.CD27..CD38._FreqParent", fillParam = "Prior.COVID.infection.", groupby="Record.ID", 
            title = "PB - Vax2", xLabel = " ", yLabel = "CD27+CD38+ (% nonnavB)") +    
  scale_y_continuous(breaks = seq(0,21,0.5), limits = c(0,3))  + scale_x_discrete(labels= c("Baseline","One Week"))
# ggsave(filename = "./Images/PB_vax2_covid_paired.pdf", width=4)

subsetData = subset(mergedData, mergedData$shortForm == 'bL' ); subsetData = subset(subsetData, subsetData$Prior.COVID.infection.=="No")
univScatter(data = subsetData, xData = "FCtfh_Vax1", yData = "FCPB_Vax1", fillParam = "Prior.COVID.infection.", title = " ", xLabel = "Fold-change Tfh", 
            yLabel = "Fold-change PB", position = "left")
subsetData = subset(mergedData, mergedData$shortForm == 'bL' ); subsetData = subset(subsetData, subsetData$Prior.COVID.infection.=="Yes")
univScatter(data = subsetData, xData = "FCtfh_Vax1", yData = "FCPB_Vax1", fillParam = "Prior.COVID.infection.", title = " ", xLabel = "Fold-change Tfh", 
            yLabel = "Fold-change PB", position = "left")

#Tbet+ CD11c+ B cells

subsetData <- subset(mergedData, mergedData$Record.ID != "CV-022")        # baseline for CD11c looks super strange
prePostTime(subsetData, xData = "day", yData="CD19_.Nonnaive.B.Tbet..CD11C._FreqParent", fillParam = "Prior.COVID.infection.", groupby="Record.ID", title = "Tbet+CD11c+", 
            xLabel = "Prior COVID?", yLabel = "Tbet+CD11+ (% nonnavB)") + 
  ggformula::geom_spline(aes_string(x="day", y="CD19_.Nonnaive.B.Tbet..CD11C._FreqParent"), color="blue",size=1, spar=0.5)

subsetData = subset(mergedData, mergedData$shortForm == 'bL' | mergedData$shortForm=='oW'); subsetData = subset(subsetData, subsetData$Prior.COVID.infection.=="No")
prePostTime(data=subsetData, xData = "shortForm", yData="CD19_.Nonnaive.B.Tbet..CD11C._FreqParent", fillParam = "Prior.COVID.infection.", groupby="Record.ID", title = "Tbet+ CD11c+ B cells - inex Vax1", 
            xLabel = " ", yLabel = "Tbet+CD11c+ (% nonnavB)") + scale_y_continuous(breaks = seq(0,21,0.5), limits = c(0,6))  + 
  scale_x_discrete(labels= c("Baseline","One Week")) + #ggrepel::geom_text_repel(aes(label = Record.ID))
# ggsave(filename = "./Images/Tbet+ CD11c+ B cells_vax1_healthy_paired.pdf", width=4)

subsetData = subset(mergedData, mergedData$shortForm == 'bL' | mergedData$shortForm=='oW'); subsetData = subset(subsetData, subsetData$Prior.COVID.infection.=="Yes")
prePostTime(data=subsetData, xData = "shortForm", yData="CD19_.Nonnaive.B.Tbet..CD11C._FreqParent", fillParam = "Prior.COVID.infection.", groupby="Record.ID", title = "Tbet+ CD11c+ B cells - exp Vax1", 
            xLabel = " ", yLabel = "Tbet+CD11c+ (% nonnavB)") + scale_y_continuous(breaks = seq(0,21,0.5), limits = c(0,6))  + 
  scale_x_discrete(labels= c("Baseline","One Week")) + #ggrepel::geom_text_repel(aes(label = Record.ID))
# ggsave(filename = "./Images/Tbet+ CD11c+ B cells_vax1_covid_paired.pdf", width=4)
  
  
subsetData = subset(mergedData, mergedData$shortForm == '3W' | mergedData$shortForm=='4W'); subsetData = subset(subsetData, subsetData$Prior.COVID.infection.=="No")
prePostTime(data=subsetData, xData = "shortForm", yData="CD19_.Nonnaive.B.Tbet..CD11C._FreqParent", fillParam = "Prior.COVID.infection.", groupby="Record.ID", 
            title = "Tbet+ CD11c+ B cells - Vax2 \ inexp", xLabel = " ", yLabel = "Tbet+CD11c+ (% nonnavB)") + scale_y_continuous(breaks = seq(0,21,0.5), limits = c(0,6))  + 
  scale_x_discrete(labels= c("Baseline","One Week"))# + ggrepel::geom_text_repel(aes(label = Record.ID))
# ggsave(filename = "./Images/Tbet+ CD11c+ B cells__vax2_healthy_paired.pdf", width=4)

subsetData = subset(mergedData, mergedData$shortForm == '3W' | mergedData$shortForm=='4W'); subsetData = subset(subsetData, subsetData$Prior.COVID.infection.=="Yes")
prePostTime(data=subsetData, xData = "shortForm", yData="CD19_.Nonnaive.B.Tbet..CD11C._FreqParent", fillParam = "Prior.COVID.infection.", groupby="Record.ID", 
            title = "Tbet+ CD11c+ B cells - Vax2 \ exp", xLabel = " ", yLabel = "Tbet+CD11c+ (% nonnavB)") + scale_y_continuous(breaks = seq(0,21,0.5), limits = c(0,6))  + 
  scale_x_discrete(labels= c("Baseline","One Week"))# + ggrepel::geom_text_repel(aes(label = Record.ID))
# ggsave(filename = "./Images/Tbet+ CD11c+ B cells__vax2_covid_paired.pdf", width=4)


#' ------------------ GzmB CD 8 analyses --------------------------

prePostTime(mergedData, xData = "DPV", yData="CD8_.CD38.Ki67..GzmB..CD8_FreqParent", fillParam = "Prior.COVID.infection.", groupby="Record.ID", title = "GzmB+ Activated CD8", 
            xLabel = "Prior COVID?", yLabel = "GzmB+ (% Activated CD8)") + 
  ggformula::geom_spline(aes_string(x="DPV", y="CD8_.CD38.Ki67..GzmB..CD8_FreqParent"), color="blue",size=1, spar=0.5) #+ ggrepel::geom_text_repel(aes(label=Record.ID),size=2)


subsetData = subset(mergedData, mergedData$shortForm == 'bL' | mergedData$shortForm=='oW'); subsetData = subset(subsetData, subsetData$Prior.COVID.infection.=="No")
prePostTime(data=subsetData, xData = "shortForm", yData="CD8_.CD38.Ki67..GzmB..CD8_FreqParent", fillParam = "Prior.COVID.infection.", groupby="Record.ID", title = "GzmB+ Activated CD8 -- Vax 1Inexperienced", 
            xLabel = " ", yLabel = "GzmB+ (% Activated CD8)") + scale_y_continuous(breaks = seq(0,100,10), limits = c(0,100)) + scale_x_discrete(labels= c("Baseline","One Week")) 
# ggsave(filename = "./Images/GzmBActivatedCD8_vax1_healthy_paired.pdf", width=4)


subsetData = subset(mergedData, mergedData$shortForm == 'bL' | mergedData$shortForm=='oW'); subsetData = subset(subsetData, subsetData$Prior.COVID.infection.=="Yes")
prePostTime(data=subsetData, xData = "shortForm", yData="CD8_.CD38.Ki67..GzmB..CD8_FreqParent", fillParam = "Prior.COVID.infection.", groupby="Record.ID", title = "GzmB+ Activated CD8 -- Vax 1Experienced", 
            xLabel = " ", yLabel = "GzmB+ (% Activated CD8)") + scale_y_continuous(breaks = seq(0,100,10), limits = c(0,100)) + scale_x_discrete(labels= c("Baseline","One Week"))
# ggsave(filename = "./Images/GzmBActivatedCD8_vax1_covid_paired.pdf", width=4)


subsetData = subset(mergedData, mergedData$shortForm == '3W' | mergedData$shortForm=='4W'); subsetData = subset(subsetData, subsetData$Prior.COVID.infection.=="No")
prePostTime(data=subsetData, xData = "shortForm", yData="CD8_.CD38.Ki67..GzmB..CD8_FreqParent", fillParam = "Prior.COVID.infection.", groupby="Record.ID", title = "GzmB+ Activated CD8 -- Vax 2 Inexperienced", 
            xLabel = " ", yLabel = "GzmB+ (% Activated CD8)") + scale_y_continuous(breaks = seq(0,100,10), limits = c(0,100)) + scale_x_discrete(labels= c("Baseline","One Week")) 
# ggsave(filename = "./Images/GzmBActivatedCD8_vax2_healthy_paired.pdf", width=4)


subsetData = subset(mergedData, mergedData$shortForm == '3W' | mergedData$shortForm=='4W'); subsetData = subset(subsetData, subsetData$Prior.COVID.infection.=="Yes")
prePostTime(data=subsetData, xData = "shortForm", yData="CD8_.CD38.Ki67..GzmB..CD8_FreqParent", fillParam = "Prior.COVID.infection.", groupby="Record.ID", title = "GzmB+ Activated CD8 --Vax 2 Experienced", 
            xLabel = " ", yLabel = "GzmB+ (% Activated CD8)") + scale_y_continuous(breaks = seq(0,100,10), limits = c(0,100)) + scale_x_discrete(labels= c("Baseline","One Week"))
# ggsave(filename = "./Images/GzmBActivatedCD8_vax2_covid_paired.pdf", width=4)

linePlot(data = mergedData, xData = 'timeCategory', yData = 'CD8_.CD38.Ki67..GzmB..CD8_FreqParent', groupby = 'Alias', xLabel = 'Time after first vaccine dose', yLabel = "GzmB+ (% Activated CD8)", 
         title = "GzmB+ Activated CD8", colorby = "Prior.COVID.infection.") + 
  scale_color_manual(name="Prior COVID?",values = c("#FFDFB1","#B5B2F1")) + scale_y_continuous(breaks = seq(0,100,10), limits = c(0,100))  #+ geom_label(aes(label=Label),size=2) #ggrepel::geom_label_repel(aes(label = Label))


