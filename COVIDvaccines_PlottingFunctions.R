library("grid")
library("ggplot2")
library("gplots")
library("gridExtra")
library("viridis")
library("reshape2")
library("gridExtra")
library("RColorBrewer")
library("scales")
library("MASS")
library("ggrepel")
library("pheatmap")
library("e1071")

## ------------------------------------------- PLOTTING FUNCTIONS ------------------------------------------------------------------------------


twoSampleBarMelted <- function (data, xData, yData, fillParam, title, yLabel)
{
  set.seed(101)
  fit <- lm(data[,yData] ~ data[,xData])
  pValue <- summary(fit)$coefficients[2,"Pr(>|t|)"];  CI <- confint(fit)[2,]; CI <- round(CI,2)
  if (pValue < 0.01)
  {
    annotationInfo <- paste0("P = ", formatC(pValue, format="e",digits=1), "\n", "Mean 95%CI: [",CI[1], ", ",CI[2], "]")
  }
  if (pValue >= 0.01)
  {
    annotationInfo <- paste0("P = ", round(pValue, 2),"\n", "Mean 95%CI: \n[",CI[1], ", ",CI[2], "]")
  }
  my_grob = grobTree(textGrob(annotationInfo, x=0.03,  y=0.85, hjust=0, gp=gpar(col="black", fontsize=28)))
  
  if( ! is.factor(data[,xData])) {  data[,xData] <- factor(data[,xData])    }
  for (i in 1:length(levels(data[,xData])))
  {
    data[, paste(levels(data[,xData])[i])] <- mean(data[which(data[,xData] == levels(data[,xData])[i]), yData], na.rm=T)   # now mean calculated for each level of xData
  }
  
  overTime <- aggregate( x = data[,yData], by= list(variable = data$variable, TimeCategory = data$TimeCategory, Cohort = data$Cohort), FUN=mean, na.rm = T)
  names(overTime)[which(names(overTime) == 'x')] <- yData

  
  return (
    ggplot(data=data, aes_string(x=xData, y=yData, fill=fillParam, width=0.6)) + scale_fill_manual(values = c("#7FAEDB", "#FFB18C")) +  
      geom_bar(data=overTime, aes_string(x=xData, y=yData), position = position_dodge(), stat = "identity") + 
      geom_point(size=7, pch=21, fill="black", color="white", alpha=0.5, position = position_jitter(width=0.1)) + 
      ggtitle(title) + ylab(yLabel) +  theme_bw() +
      theme(axis.text = element_text(size=28,hjust = 0.5), axis.title = element_text(size=28,hjust = 0.5), axis.title.x = element_blank(), plot.title = element_text(size=32,hjust = 0.5)) + 
      annotation_custom(my_grob) + theme(legend.position = "none")
  )
}

twoSampleBar <- function (data, xData, yData, fillParam, title, yLabel, batch="none", position = "left", FCplot=F, confInt=F, ttest=T, nonparam=F)
{
  set.seed(102)
  if (ttest == T && batch == "none" && nonparam == F)
  {
    justforttest <- data[, c(xData,yData)]
    fit <- rstatix::t_test(justforttest, formula = as.formula(paste(colnames(justforttest)[2], "~", paste(colnames(justforttest)[1]), sep = "") ))
    pValue <- fit$p
  }
  if (ttest == T && batch == "none" && nonparam == T)
  {
    justforttest <- data[, c(xData,yData)]
    fit <- rstatix::wilcox_test(justforttest, formula = as.formula(paste(colnames(justforttest)[2], "~", paste(colnames(justforttest)[1]), sep = "") ))
    pValue <- fit$p
  }
  if (ttest == F)
  {
    if (batch == "none")  {
      fit <- wilcox_test(data = data, as.formula(paste(yData, xData, sep="~") ));  pValue <- fit$p  }
    if (batch != "none")  {  fit <- lm(data[,yData] ~ data[,xData] + data[,batch]); pValue <- summary(fit)$coefficients[2,"Pr(>|t|)"];  CI <- confint(fit)[2,]; CI <- round(CI,2)    }
  }
  if (! is.nan(pValue) && confInt == T)
  {
    if (pValue < 0.01)
    {    annotationInfo <- paste0("P = ", formatC(pValue, format="e",digits=1), "\n", "Mean 95%CI: [",CI[1], ", ",CI[2], "]")     }
    if (pValue >= 0.01)
    {    annotationInfo <- paste0("P = ", round(pValue, 2),"\n", "Mean 95%CI: \n[",CI[1], ", ",CI[2], "]")      }
    if(position == "left")  { my_grob = grobTree(textGrob(annotationInfo, x=0.03,  y=0.9, hjust=0, gp=gpar(col="black", fontsize=28)))   }
    if(position == "right")  { my_grob = grobTree(textGrob(annotationInfo, x=0.65,  y=0.9, hjust=0, gp=gpar(col="black", fontsize=28)))   }
    if(position == "none") { my_grob = grobTree(textGrob(annotationInfo, x=10,  y=10, hjust=0, gp=gpar(col="black", fontsize=1)))   }
  }
  if (! is.nan(pValue) && confInt == F)
  {
    if (pValue < 0.01)
    {    annotationInfo <- paste0("P = ", formatC(pValue, format="e",digits=1))     }
    if (pValue >= 0.01)
    {    annotationInfo <- paste0("P = ", round(pValue, 3))      }
    if(position == "left")  { my_grob = grobTree(textGrob(annotationInfo, x=0.03,  y=0.93, hjust=0, gp=gpar(col="black", fontsize=28)))   }
    if(position == "right")  { my_grob = grobTree(textGrob(annotationInfo, x=0.65,  y=0.93, hjust=0, gp=gpar(col="black", fontsize=28)))   }
    if(position == "none") { my_grob = grobTree(textGrob(annotationInfo, x=10,  y=10, hjust=0, gp=gpar(col="black", fontsize=1)))   }
    
  }
  
  if( ! is.factor(data[,xData])) {  data[,xData] <- factor(data[,xData])    }

  overTime <- aggregate(x = data[,yData], by= list( data[,xData]), FUN=mean, na.rm = T)
  names(overTime)[which(names(overTime) == 'Group.1')] <- xData
  names(overTime)[which(names(overTime) == 'x')] <- yData
  # colors:   COVID-exp B5B2F1 (purple)             COVID-naive FFDFB1 (orange)
  if (FCplot == T)

    return (
      ggplot(data=data, aes_string(x=xData, y=yData, fill=fillParam, width=0.6)) + scale_fill_manual(values = c("#FFDFB1", "#B5B2F1")) + 
        geom_hline(yintercept=1, linetype="dashed", color = "black", size=0.5) + 
        geom_bar(data=overTime, aes_string(x=xData, y=yData), position = position_dodge(), stat = "identity", color='black',size=0.1) + 
        geom_point(size=7, pch=21, fill="black", color="white", alpha=0.5, position = position_jitter(width=0.15)) + 
        ggtitle(title) + ylab(yLabel) +  theme_bw() +
        theme(axis.text = element_text(size=28,color="black"), axis.title = element_text(size=28,hjust = 0.5), axis.title.x = element_blank(), 
              plot.title = element_text(size=32,hjust = 0.5), axis.text.x = element_text(hjust = 0.5, vjust=0.5, angle=0)) + 
        annotation_custom(my_grob) + theme(legend.position = "none")
    )
  if (FCplot == F)
    return (
      ggplot(data=data, aes_string(x=xData, y=yData, fill=fillParam, width=0.6)) + scale_fill_manual(values = c("#FFDFB1", "#B5B2F1")) + 
        geom_bar(data=overTime, aes_string(x=xData, y=yData), position = position_dodge(), stat = "identity", color='black',size=0.1) + 
        # geom_point(size=7, pch=21, fill="black", color="white", alpha=0.5, position = position_jitter(width=0.15)) + 
        ggbeeswarm::geom_quasirandom(size=7, pch=21, fill="black",color="white",alpha=0.5, width = 0.2) + 
        ggtitle(title) + ylab(yLabel) +  theme_bw() +
        theme(axis.text = element_text(size=28,color="black"), axis.title = element_text(size=28,hjust = 0.5), axis.title.x = element_blank(), 
              plot.title = element_text(size=32,hjust = 0.5), axis.text.x = element_text(hjust = 0.5, vjust=0.5, angle=0)) + 
        annotation_custom(my_grob) + theme(legend.position = "none")
    )
}


univScatter <- function(data, xData, yData, fillParam, title, xLabel, yLabel, position = "left")
{
  pearson <- round(cor(data[,xData], data[,yData], method = "pearson", use = "complete.obs"), 2)
  pValue <- cor.test(data[,xData], data[,yData], method="pearson")
  if (pValue$p.value < 0.01)   {    annotationInfo <- paste0("Pearson r = ", pearson,"\n","P = ", formatC(pValue$p.value, format="e", digits=1))    }
  if (pValue$p.value >= 0.01) {    annotationInfo <- paste0("Pearson r = ", pearson,"\n","P = ", round(pValue$p.value,2))   }
  if(position == "left")  { my_grob = grobTree(textGrob(annotationInfo, x=0.05,  y=0.9, hjust=0, gp=gpar(col="black", fontsize=28)))   }
  if(position == "right")  { my_grob = grobTree(textGrob(annotationInfo, x=0.45,  y=0.9, hjust=0, gp=gpar(col="black", fontsize=28)))   }
  if(position == "none") { my_grob = grobTree(textGrob(annotationInfo, x=10,  y=10, hjust=0, gp=gpar(col="black", fontsize=1)))   }
  # my_grob = grobTree(textGrob(annotationInfo, x=0.05,  y=0.9, hjust=0, gp=gpar(col="black", fontsize=28)))
  return (
    ggplot(data ) + 
      geom_smooth(aes_string(x=xData, y=yData, color=fillParam, fill=fillParam), method='lm') + 
      geom_point(aes_string(x=xData, y=yData, fill=fillParam), size=8, color="black", pch=21) + theme_bw() + 
      scale_color_manual(values=c("black")) + 
      scale_fill_manual(values=c("grey90")) + 
      ggtitle(title) + ylab(yLabel) + xlab(xLabel)  +
      theme(axis.text = element_text(size=28,hjust = 0.5, color="black"), axis.title = element_text(size=28,hjust = 0.5), plot.title = element_text(size=32,hjust = 0.5)) + 
      annotation_custom(my_grob) + theme(legend.position = "none")     )
}

bivScatter <- function(data1, data2, name1, name2, xData, yData, fillParam, title, xLabel, yLabel, statsOff = F)
{
  if (statsOff == F)
  {  
    pearson1 <- round(cor(data1[,xData], data1[,yData], method = "pearson", use = "complete.obs"), 2)
    pValue1 <- cor.test(data1[,xData], data1[,yData], method="pearson")
    pearson2 <- round(cor(data2[,xData], data2[,yData], method = "pearson", use = "complete.obs"), 2)
    pValue2 <- cor.test(data2[,xData], data2[,yData], method="pearson")
    if (pValue1$p.value < 0.01 | pValue2$p.value < 0.01)   {    
      annotationInfo1 <- paste0(name1, " Pearson r = ", pearson1,"\n","P = ", formatC(pValue1$p.value, format="e", digits=1) )
      annotationInfo2 <- paste0("\n", name2, " Pearson r = ", pearson2,"\n","P = ", formatC(pValue2$p.value, format="e", digits=1) )  }
    if (pValue1$p.value >= 0.01 & pValue2$p.value >= 0.01) {    
      annotationInfo1 <- paste0(name1, " Pearson r = ", pearson1,"\n","P = ", round(pValue1$p.value,2) )
      annotationInfo2 <- paste0("\n",name2," Pearson r = ", pearson2, "\n","P = ", round(pValue2$p.value,2))   }
    my_grob1 = grobTree(textGrob(annotationInfo1, x=0.05,  y=0.88, hjust=0, gp=gpar(col="#7FAEDB", fontsize=28)))
    my_grob2 = grobTree(textGrob(annotationInfo2, x=0.05,  y=0.75, hjust=0, gp=gpar(col="#FFB18C", fontsize=28)))
    return (
      ggplot() + 
        geom_point(data = data1, aes_string(x=xData, y=yData, fill= fillParam), size=8, color="black", pch=21) + theme_bw() + 
        geom_point(data = data2, aes_string(x=xData, y=yData, fill=fillParam), size=8, color="black", pch=21) + theme_bw() + 
        geom_smooth(data = data1, aes_string(x=xData, y=yData, color=fillParam, fill=NA), method='lm', fullrange=T) +
        geom_smooth(data = data2, aes_string(x=xData, y=yData, color=fillParam, fill=NA), method='lm', fullrange=T) +
        scale_color_manual(values=c("#FFB18C", "#7FAEDB")) + 
        scale_fill_manual(values=c("#FFB18C", "#7FAEDB")) + 
        ggtitle(title) + ylab(yLabel) + xlab(xLabel)  +
        theme(axis.text = element_text(size=28,hjust = 0.5,color="black"), axis.title = element_text(size=28,hjust = 0.5), plot.title = element_text(size=32,hjust = 0.5)) + 
        annotation_custom(my_grob1) + annotation_custom(my_grob2) + 
        theme(legend.position = "none")     )
  }
  if (statsOff == T)
  {  
    return (
      ggplot() + 
        geom_point(data = data1, aes_string(x=xData, y=yData, fill= fillParam), size=8, color="black", pch=21) + theme_bw() + 
        geom_point(data = data2, aes_string(x=xData, y=yData, fill=fillParam), size=8, color="black", pch=21) + theme_bw() + 
        geom_smooth(data = data1, aes_string(x=xData, y=yData, color=fillParam, fill=NA), method='lm', fullrange=T) +
        geom_smooth(data = data2, aes_string(x=xData, y=yData, color=fillParam, fill=NA), method='lm', fullrange=T) +
        scale_color_manual(values=c("#FFB18C", "#7FAEDB")) + 
        scale_fill_manual(values=c("#FFB18C", "#7FAEDB")) + 
        ggtitle(title) + ylab(yLabel) + xlab(xLabel)  +
        theme(axis.text = element_text(size=28,hjust = 0.5,color="black"), axis.title = element_text(size=28,hjust = 0.5), plot.title = element_text(size=32,hjust = 0.5)) + 
        #annotation_custom(my_grob1) + annotation_custom(my_grob2) + 
        theme(legend.position = "none")     )
  }
}

prePostTime <- function(data, xData, yData, fillParam, groupby, title, xLabel, yLabel, repMeasures=T, exponential=F)
{
  data[,fillParam] <- factor(data[,fillParam])
  
  if(repMeasures==F  && exponential==F)              # original, not repeated measures, permits several cohorts, intended for continuous variable
  {
    targets <- which(table(data$Record.ID) > 0)         
    subsetData <- data[ which( data$Record.ID %in% names(targets)   ), ]
    subsetData <- subsetData[order(subsetData$Record.ID, subsetData$shortForm, decreasing = F),]
    if(length(levels(as.factor(data[,fillParam])))>1)
    {  
      return(          # colors:   COVID-exp B5B2F1 (purple)             COVID-naive FFDFB1 (orange)
        ggplot(data=subsetData, aes_string(x=xData, y=yData, fill=fillParam) ) + theme_bw() + 
          geom_path(aes_string(group=groupby), color="grey70", alpha=0.95) + 
          geom_point(size = 5, pch=21, color="black", alpha=0.5) + facet_wrap(fillParam ) +   # , scales='free'
          scale_color_manual(values=c("#FFDFB1", "#B5B2F1")) + 
          scale_fill_manual(values=c("#FFDFB1", "#B5B2F1")) + 
          ggtitle(title) + ylab(yLabel) + xlab(xLabel)  +
          theme(axis.text = element_text(size=18,hjust = 0.5, color="black"), axis.title = element_text(size=22,hjust = 0.5), 
                plot.title = element_text(size=36,hjust = 0.5), axis.text.x = element_text(angle=45, hjust=1,vjust=1),
                legend.position = "none", strip.text = element_text(size = 24, color="black"), strip.background = element_rect(fill="white"))  
      )
    }
  }
  
  ##  ----------------------------------------------------------------------------------------------------------------
  if(repMeasures==F  && length(levels(data[,fillParam]))==2  && exponential==T)        # not repeated measures, two cohorts, exponential scale
  {
    targets <- which(table(data$Record.ID) > 0)        
    subsetData <- data[ which( data$Record.ID %in% names(targets)   ), ]
    subsetData <- subsetData[order(subsetData$Record.ID, subsetData$shortForm, decreasing = F),]
    if(length(levels(as.factor(data[,fillParam])))>1)     
    {  
      return(          # colors:   COVID-exp B5B2F1 (purple)             COVID-naive FFDFB1 (orange)
        ggplot(data=subsetData, aes_string(x=xData, y=yData, fill=fillParam) ) + theme_bw() + 
          geom_path(aes_string(group=groupby), color="grey70", alpha=0.95) + 
          geom_point(size = 8, pch=21, color="black", alpha=0.5) + facet_wrap(fillParam ) +   # , scales='free'
          scale_color_manual(values=c("#FFDFB1", "#B5B2F1")) + 
          scale_fill_manual(values=c("#FFDFB1", "#B5B2F1")) + 
          ggtitle(title) + ylab(yLabel) + xlab(xLabel)  +
          theme(axis.text = element_text(size=18,hjust = 0.5, color="black"), axis.title = element_text(size=22,hjust = 0.5), 
                plot.title = element_text(size=36,hjust = 0.5), axis.text.x = element_text(angle=45, hjust=1,vjust=1),
                legend.position = "none", strip.text = element_text(size = 24, color="black"), strip.background = element_rect(fill="white")) + 
          scale_y_continuous(trans='pseudo_log', limits = c(0,10000), breaks=c(10^(0:4)), labels=trans_format('log10',math_format(10^.x)), minor_breaks =5*10^(0:10))
        )
    }
  }
  
  if(repMeasures==F && length(levels(data[,fillParam]))==2  && exponential==F)         # not repeated measures, two cohorts, normal scale
  {
    targets <- which(table(data$Record.ID) > 0)       
    subsetData <- data[ which( data$Record.ID %in% names(targets)   ), ]
    subsetData <- subsetData[order(subsetData$Record.ID, subsetData$shortForm, decreasing = F),]
    if(length(levels(as.factor(data[,fillParam])))>1)
    {  
      return(          # colors:   COVID-exp B5B2F1 (purple)             COVID-naive FFDFB1 (orange)
        ggplot(data=subsetData, aes_string(x=xData, y=yData, fill=fillParam) ) + theme_bw() + 
          geom_path(aes_string(group=groupby), color="grey70", alpha=0.95) + 
          geom_point(size = 8, pch=21, color="black", alpha=0.5) + facet_wrap(fillParam ) +   # , scales='free'
          scale_color_manual(values=c("#FFDFB1", "#B5B2F1")) + 
          scale_fill_manual(values=c("#FFDFB1", "#B5B2F1")) + 
          ggtitle(title) + ylab(yLabel) + xlab(xLabel)  +
          theme(axis.text = element_text(size=18,hjust = 0.5, color="black"), axis.title = element_text(size=22,hjust = 0.5), 
                plot.title = element_text(size=36,hjust = 0.5), axis.text.x = element_text(angle=45, hjust=1,vjust=1),
                legend.position = "none", strip.text = element_text(size = 24, color="black"), strip.background = element_rect(fill="white")) 
      )
    }
  }
  
  ##  ----------------------------------------------------------------------------------------------------------------
  if(repMeasures==T && length(levels(data[,fillParam]))==1 && exponential==T)          #  repeated measures, one cohort, exponential scale
  {
    targets <- which(table(data$Record.ID) > 1)         # show only the repeated measures values
    subsetData <- data[ which( data$Record.ID %in% names(targets)   ), ]; subsetData[,xData] <- factor(subsetData[,xData])
    justforttest <- subsetData[, c(xData,yData)]
    justforttest[,xData] <- factor(justforttest[,xData])
    fit <- rstatix::t_test(justforttest, formula = as.formula(paste(colnames(justforttest)[2], "~", paste(colnames(justforttest)[1]), sep = "") ), paired = T)
    pValue <- fit$p
    annotationInfo <- paste0("P = ", formatC(pValue, format="e",digits=1))     
    my_grob = grobTree(textGrob(annotationInfo, x=0.03,  y=0.93, hjust=0, gp=gpar(col="black", fontsize=24))) 
    return( 
      ggplot(data=subsetData, aes_string(x=xData, y=yData,fill=fillParam) ) + theme_bw() + 
        geom_path(aes_string(group=groupby), color="grey70", alpha=0.5) + 
        geom_point(size = 8, color="black", shape=21, alpha=0.5)  +  
        scale_color_manual(values=c("#FFDFB1")) + 
        scale_fill_manual(values=c("#FFDFB1")) + 
        ggtitle(title) + ylab(yLabel) + xlab(xLabel)  +
        theme(axis.text = element_text(size=20,hjust = 0.5, color="black"), axis.title = element_text(size=28,hjust = 0.5), plot.title = element_text(size=28,hjust = 0.5), 
              legend.position = "none", strip.text = element_text(size = 24, color="black"), strip.background = element_rect(fill="white")) + annotation_custom(my_grob) + 
        scale_y_continuous(trans='pseudo_log', limits = c(0,10000), breaks=c(10^(0:4)), labels=trans_format('log10',math_format(10^.x)), minor_breaks =5*10^(0:10))
    )
  }
  
  if(repMeasures==T && length(levels(data[,fillParam]))==1 && exponential==F)          #  repeated measures, one cohort, normal scale
  {
    targets <- which(table(data$Record.ID) > 1)         # show only the repeated measures values
    subsetData <- data[ which( data$Record.ID %in% names(targets)   ), ]; subsetData[,xData] <- factor(subsetData[,xData])
    justforttest <- subsetData[, c(xData,yData)]
    justforttest[,xData] <- factor(justforttest[,xData])
    fit <- rstatix::t_test(justforttest, formula = as.formula(paste(colnames(justforttest)[2], "~", paste(colnames(justforttest)[1]), sep = "") ), paired = T)
    pValue <- fit$p
    annotationInfo <- paste0("P = ", formatC(pValue, format="e",digits=1))     
    my_grob = grobTree(textGrob(annotationInfo, x=0.03,  y=0.93, hjust=0, gp=gpar(col="black", fontsize=24))) 
    return( 
      ggplot(data=subsetData, aes_string(x=xData, y=yData,fill=fillParam) ) + theme_bw() + 
        geom_path(aes_string(group=groupby), color="grey70", alpha=0.5) + 
        geom_point(size = 8, color="black", shape=21, alpha=0.5)  + # facet_wrap(fillParam ) + 
        scale_color_manual(values=c("#FFDFB1")) + 
        scale_fill_manual(values=c("#FFDFB1")) + 
        ggtitle(title) + ylab(yLabel) + xlab(xLabel)  +
        theme(axis.text = element_text(size=20,hjust = 0.5, color="black"), axis.title = element_text(size=28,hjust = 0.5), plot.title = element_text(size=28,hjust = 0.5), 
              legend.position = "none", strip.text = element_text(size = 24, color="black"), strip.background = element_rect(fill="white")) + annotation_custom(my_grob)
    )
  }
  
  ##  ----------------------------------------------------------------------------------------------------------------
  if(repMeasures==T & length(levels(data[,fillParam]))==2 && exponential==T)           #  repeated measures, two cohorts, exponential scale
  {
    targets <- which(table(data$Record.ID) > 1)         # show only the repeated measures values
    subsetData <- data[ which( data$Record.ID %in% names(targets)   ), ]; subsetData[,xData] <- factor(subsetData[,xData])
    justforttest <- subsetData[, c(xData,yData,fillParam)]
    justforttest[,xData] <- factor(justforttest[,xData])
    temp <- justforttest %>% group_by(eval(as.name(fillParam))) %>% levene_test(eval(as.name(yData)) ~ eval(as.name(xData)))
    if( any(temp$p < 0.05)) { justforttest[,yData] <- log(justforttest[,yData]+1)}
    
    fit <- justforttest %>% group_by(eval(as.name(fillParam))) %>% 
      rstatix::t_test(formula = as.formula(paste(colnames(justforttest)[2], "~", paste(colnames(justforttest)[1]), sep = "") ), paired = T)
    
    temp <- factor(data[,fillParam])
    data.1 <- subsetData[which(subsetData[,fillParam] == levels(temp)[1]),]
    pValue <- fit$p[1]
    annotationInfo <- paste0("P = ", formatC(pValue, format="e",digits=1))     
    my_grob = grobTree(textGrob(annotationInfo, x=0.03,  y=0.93, hjust=0, gp=gpar(col="black", fontsize=24))) 
    a.1 <- ggplot(data=data.1, aes_string(x=xData, y=yData,fill=fillParam) ) + theme_bw() + 
      geom_path(aes_string(group=groupby), color="grey70", alpha=0.5) + 
      geom_point(size = 8, color="black", shape=21, alpha=0.5)  +  
      scale_color_manual(values=c("#FFDFB1")) + 
      scale_fill_manual(values=c("#FFDFB1")) + 
      ggtitle(title) + ylab(yLabel) + xlab(xLabel)  +
      theme(axis.text = element_text(size=20,hjust = 0.5, color="black"), axis.title = element_text(size=28,hjust = 0.5), plot.title = element_text(size=28,hjust = 0.5), 
            axis.text.x = element_text(angle=45,hjust=1,vjust=1),
            legend.position = "none", strip.text = element_text(size = 24, color="black"), strip.background = element_rect(fill="white")) + annotation_custom(my_grob) + 
      scale_y_continuous(trans='pseudo_log', limits = c(0,10000), breaks=c(10^(0:4)), labels=trans_format('log10',math_format(10^.x)), minor_breaks =5*10^(0:10))
    
    data.2 <-subsetData[which(subsetData[,fillParam] == levels(temp)[2]),]
    pValue <- fit$p[2]
    annotationInfo <- paste0("P = ", formatC(pValue, format="e",digits=1))     
    my_grob = grobTree(textGrob(annotationInfo, x=0.03,  y=0.93, hjust=0, gp=gpar(col="black", fontsize=24))) 
    a.2 <- ggplot(data=data.2, aes_string(x=xData, y=yData,fill=fillParam) ) + theme_bw() + 
      geom_path(aes_string(group=groupby), color="grey70", alpha=0.5) + 
      geom_point(size = 8, color="black", shape=21, alpha=0.5)  +   
      scale_color_manual(values=c("#B5B2F1")) + 
      scale_fill_manual(values=c("#B5B2F1")) + 
      ggtitle(title) +  xlab(xLabel)  + # ylab(yLabel) +
      theme(axis.text = element_text(size=20,hjust = 0.5, color="black"), axis.title = element_text(size=28,hjust = 0.5), plot.title = element_text(size=28,hjust = 0.5), 
            axis.text.x = element_text(angle=45,hjust=1,vjust=1), axis.title.y = element_blank(), 
            legend.position = "none", strip.text = element_text(size = 24, color="black"), strip.background = element_rect(fill="white")) + annotation_custom(my_grob) + 
      scale_y_continuous(trans='pseudo_log', limits = c(0,10000), breaks=c(10^(0:4)), labels=trans_format('log10',math_format(10^.x)), minor_breaks =5*10^(0:10))
    return(grid.arrange(a.1,a.2, nrow=1))
  }
  
  if(repMeasures==T & length(levels(data[,fillParam]))==2 && exponential==F)           #  repeated measures, two cohorts, normal scale
  {
    targets <- which(table(data$Record.ID) > 1)         # show only the repeated measures values
    subsetData <- data[ which( data$Record.ID %in% names(targets)   ), ]; subsetData[,xData] <- factor(subsetData[,xData])
    justforttest <- subsetData[, c(xData,yData,fillParam)]
    justforttest[,xData] <- factor(justforttest[,xData])
    temp <- justforttest %>% group_by(eval(as.name(fillParam))) %>% levene_test(eval(as.name(yData)) ~ eval(as.name(xData)))
    if( any(temp$p < 0.05)) { justforttest[,yData] <- log(justforttest[,yData]+1)}
    range.data <- range(data[,yData])
    range.data <- round(c(range.data[1]*0.25, range.data[2]*1.5),digits = 0)
    breaks <- seq(range.data[1], range.data[2], round((range.data[2]-range.data[1])/5, digits=0))
    if( length(levels(data[,xData])) == 2 )             # two time point Pre-Post, so do paired statistics
    {
      fit <- justforttest %>% group_by(eval(as.name(fillParam))) %>% 
        rstatix::t_test(formula = as.formula(paste(colnames(justforttest)[2], "~", paste(colnames(justforttest)[1]), sep = "") ), paired = T)
      
      temp <- factor(data[,fillParam])
      data.1 <- subsetData[which(subsetData[,fillParam] == levels(temp)[1]),]
      pValue <- fit$p[1]
      annotationInfo <- paste0("P = ", formatC(pValue, format="e",digits=1))     
      my_grob = grobTree(textGrob(annotationInfo, x=0.03,  y=0.93, hjust=0, gp=gpar(col="black", fontsize=24))) 
      a.1 <- ggplot(data=data.1, aes_string(x=xData, y=yData,fill=fillParam) ) + theme_bw() + 
        geom_path(aes_string(group=groupby), color="grey70", alpha=0.5) + 
        geom_point(size = 8, color="black", shape=21, alpha=0.5)  + # facet_wrap(fillParam ) + 
        scale_color_manual(c("#FFDFB1")) + 
        scale_fill_manual(values="#2861BC") + #c("#FFDFB1")) + 
        ggtitle(title) + ylab(yLabel) + xlab(xLabel)  +
        theme(axis.text = element_text(size=20,hjust = 0.5, color="black"), axis.title = element_text(size=28,hjust = 0.5), plot.title = element_text(size=28,hjust = 0.5), 
              axis.text.x = element_text(angle=45,hjust=1,vjust=1), legend.position = "none", strip.text = element_text(size = 24, color="black"), 
              strip.background = element_rect(fill="white")) + annotation_custom(my_grob) + scale_y_continuous(limits=range.data, breaks=breaks)
      
      data.2 <-subsetData[which(subsetData[,fillParam] == levels(temp)[2]),]
      pValue <- fit$p[2]
      annotationInfo <- paste0("P = ", formatC(pValue, format="e",digits=1))     
      my_grob = grobTree(textGrob(annotationInfo, x=0.03,  y=0.93, hjust=0, gp=gpar(col="black", fontsize=24))) 
      a.2 <- ggplot(data=data.2, aes_string(x=xData, y=yData,fill=fillParam) ) + theme_bw() + 
        geom_path(aes_string(group=groupby), color="grey70", alpha=0.5) + 
        geom_point(size = 8, color="black", shape=21, alpha=0.5)  + # facet_wrap(fillParam ) + 
        scale_color_manual(values=c("#B5B2F1")) + 
        scale_fill_manual(values="orange") + #c("#B5B2F1")) + 
        ggtitle(title) +  xlab(xLabel)  + # ylab(yLabel) +
        theme(axis.text = element_text(size=20,hjust = 0.5, color="black"), axis.title = element_text(size=28,hjust = 0.5), plot.title = element_text(size=28,hjust = 0.5), 
              axis.text.x = element_text(angle=45,hjust=1,vjust=1), axis.title.y = element_blank(), legend.position = "none", strip.text = element_text(size = 24, color="black"), 
              strip.background = element_rect(fill="white")) + annotation_custom(my_grob) + scale_y_continuous(limits=range.data, breaks=breaks)  
      return(grid.arrange(a.1,a.2, nrow=1))
    }  
    if( length(levels(data[,xData])) > 2 )              # >2 time points, so hold off on statistics here
    {
      temp <- factor(data[,fillParam])
      data.1 <- subsetData[which(subsetData[,fillParam] == levels(temp)[1]),]
      a.1 <- ggplot(data=data.1, aes_string(x=xData, y=yData,fill=fillParam) ) + theme_bw() + 
        geom_path(aes_string(group=groupby), color="grey70", alpha=0.5) + 
        geom_point(size = 8, color="black", shape=21, alpha=0.5)  + 
        scale_color_manual(values=c("#FFDFB1")) + 
        scale_fill_manual(values=c("#FFDFB1")) + 
        ggtitle(title) + ylab(yLabel) + xlab(xLabel)  +
        theme(axis.text = element_text(size=20,hjust = 0.5, color="black"), axis.title = element_text(size=28,hjust = 0.5), plot.title = element_text(size=28,hjust = 0.5), 
              axis.text.x = element_text(angle=45,hjust=1,vjust=1),
              legend.position = "none", strip.text = element_text(size = 24, color="black"), strip.background = element_rect(fill="white"))  
      
      data.2 <-subsetData[which(subsetData[,fillParam] == levels(temp)[2]),]
      a.2 <- ggplot(data=data.2, aes_string(x=xData, y=yData,fill=fillParam) ) + theme_bw() + 
        geom_path(aes_string(group=groupby), color="grey70", alpha=0.5) + 
        geom_point(size = 8, color="black", shape=21, alpha=0.5)  +  
        scale_color_manual(values=c("#B5B2F1")) + 
        scale_fill_manual(values=c("#B5B2F1")) + 
        ggtitle(title) +  xlab(xLabel)  + # ylab(yLabel) +
        theme(axis.text = element_text(size=20,hjust = 0.5, color="black"), axis.title = element_text(size=28,hjust = 0.5), plot.title = element_text(size=28,hjust = 0.5), 
              axis.text.x = element_text(angle=45,hjust=1,vjust=1), axis.title.y = element_blank(), 
              legend.position = "none", strip.text = element_text(size = 24, color="black"), strip.background = element_rect(fill="white"))   
      return(grid.arrange(a.1,a.2, nrow=1))
    }  
  }

}

prePostTimeAveraged <- function(data, title, xLabel, yLabel)
{
  subsetData <- data
  subsetData <- subsetData[ which( !is.na(subsetData$value) ),]
  overTime <- aggregate( subsetData$value, by= list(subsetData$variable, subsetData$TimeCategory, subsetData$Cohort), FUN=mean, na.rm = T); names(overTime)[4] <- "Mean"
  overTimeSD <- aggregate( subsetData$value, by= list(subsetData$variable, subsetData$TimeCategory, subsetData$Cohort), FUN=sd, na.rm = T); names(overTimeSD)[4] <- "SD"
  temp <- merge(overTime, overTimeSD, by = c("Group.1","Group.2", "Group.3"))
  #overTimeN <- t(as.matrix(table(subsetData$Cohort, subsetData$TimeCategory)))
  #overTimeSD$SE <- overTimeSD$x / sqrt(as.vector(overTimeN[1:nrow(overTime)]))              
  #overTime$SE <- overTimeSD$SE
  temp2 <- as.data.frame(unlist(t(as.matrix(table(subsetData$Cohort, subsetData$TimeCategory)))));  names(temp2)[3] <- "N"
  temp3 <- merge(temp, temp2, by.x = c("Group.2","Group.3"), by.y = c("Var1","Var2"))
  overTime <- temp3
  overTime$SE <- overTime$SD / sqrt(overTime$N)
  # temp <- colnames(overTime); temp[which(temp == "x")] <- "Ave"; colnames(overTime) <- temp
  annotationInfo1 <- paste0(unique(overTime$Group.3)[1])  # aPD1 group should get the peach FFB18C color
  my_grob1 = grobTree(textGrob(annotationInfo1, x=0.05,  y=0.93, hjust=0, gp=gpar(col="#FFB18C", fontsize=28)))
  annotationInfo2 <- paste0(unique(overTime$Group.3)[2])  # healthy group should get the blue 7FAEDB color
  my_grob2 = grobTree(textGrob(annotationInfo2, x=0.05,  y=0.83, hjust=0, gp=gpar(col="#7FAEDB", fontsize=28)))
  
  return(
    ggplot(data=overTime, aes(x=Group.2, y=Mean, group = Group.3, color=Group.3)) + theme_bw() +   # group.1 = variable, Group.2 = TimeCategory, Group.3 = Cohort
      geom_ribbon(aes(x=Group.2, ymin=Mean-SE, ymax=Mean+SE, fill=Group.3), alpha=0.1) + 
      geom_line(aes(size=3)) +  
      geom_point(aes(size=4, fill=Group.3), pch=21) + 
      scale_color_manual(values=c("#7FAEDB", "#FFB18C")) +     
      scale_fill_manual(values=c("#7FAEDB", "#FFB18C")) + 
      ggtitle(title) + ylab(yLabel) + xlab(xLabel)  +
      theme(axis.text = element_text(size=24,hjust = 0.5, color="black"), axis.title = element_text(size=28,hjust = 0.5), plot.title = element_text(size=36,hjust = 0.5), 
            axis.text.x = element_text(angle = 45, hjust=1), legend.position = "none")  + annotation_custom(my_grob1) + annotation_custom(my_grob2)
  )
}


linePlot <- function(data, xData, yData, groupby, colorby, title, xLabel, yLabel)
{
  subsetData <- subset(data, !is.na(mergedData[,yData]));  #subsetData$subsetData[order(subsetData$Label)]
  return(
    ggplot(subsetData, aes_string(x=xData, y=yData, group=groupby)) + 
    geom_point(size=2,alpha=0.5) + geom_line(aes_string(group=groupby, color=colorby), size=1, alpha=0.75) +  
    xlab(xLabel) + ylab(yLabel) + theme_bw() + 
    ggtitle(title) + #scale_x_discrete(labels= c("Baseline","1 week","2 weeks","3 weeks","4 weeks")) + 
    theme(axis.text = element_text(color="black",size=18), axis.title = element_text(size=24), axis.text.x = element_text(angle=45, hjust=1,vjust=1),
          plot.title=element_text(size=28) )
  )
  
}


