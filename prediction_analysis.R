library(deSolve)
library(data.table)
library(lhs)
library(foreach)
library(doParallel)
library(ggplot2)
library(randomForest)

# Use models to predict clustering
# Use LHS results
# m = dependent variable (mean arc distance)
# LHS parameters = predictors
# -> identify the most important predictors


#### Data_all - 1 species ####

data_all <- out_dt_last_steps[time==simLength]
colnames(data_all)
setnames(data_all, old=c("rMob","rDup"), new=c("rMob_1","rDup_1"))
## Get data ##
data_importance_all <- NULL
lm1_coef <- NULL
## Arc distance prediction with random forest
for(trt in 1:4){ # treatment loop
  print(paste0("Treatment ",trt))
  data_predict <- data_all[time==simLength][x11_1>1E-6][treatment==trt]
  ## Data transformation ##
  # Log #
  #data_predict[, m_1 := log10(rMob_2)]
  data_predict[, rMob_1 := log10(rMob_1)]
  data_predict[, rDup_1 := log10(rDup_1)]
  data_predict[, rUpt_1 := log10(rUpt_1)]
  # data_predict[, rClear := log10(rClear)]
  # data_predict[, fcost := log10(fcost)]
  # data_predict[, deathBasal := log10(deathBasal)]
  # data_predict[, deathIncrease := log10(deathIncrease)]
  # data_predict[, deathPeriod := log10(deathPeriod)]
  # Scale #
  data_predict[, rMob_1 := scale(rMob_1,center=TRUE,scale=TRUE)]
  data_predict[, rDup_1 := scale(rDup_1,center=TRUE,scale=TRUE)]
  data_predict[, rUpt_1 := scale(rUpt_1,center=TRUE,scale=TRUE)]
  data_predict[, rClear := scale(rClear,center=TRUE,scale=TRUE)]
  data_predict[, fcost := scale(fcost,center=TRUE,scale=TRUE)]
  data_predict[, deathBasal := scale(deathBasal,center=TRUE,scale=TRUE)]
  data_predict[, deathIncrease := scale(deathIncrease,center=TRUE,scale=TRUE)]
  data_predict[, deathPeriod := scale(deathPeriod,center=TRUE,scale=TRUE)]
  
  ## Predictors ##
  # if(trt==1){
  #   predictors <- c("rMob_1","rDup_1","rUpt_1","rClear","fcost","deathBasal") # no atb
  # }
  # if(trt==2){
  #   predictors <- c("rMob_1","rDup_1","rUpt_1","rClear","fcost","deathBasal","deathIncrease") # cst atb
  # }
  # if(trt==3 | trt==4){
  #   predictors <- c("rMob_1","rDup_1","rUpt_1","rClear","fcost","deathBasal","deathIncrease","deathPeriod") # synchro or altern atb
  # }
  # if(trt==5){
  #   predictors <- c("rMob_1","rDup_1","rUpt_1","rClear","fcost","deathBasal","deathIncrease","deathPeriod","lifePeriod") # uneven synchro atb
  # }
  #predictors <- c("rMob_1","rDup_1","rUpt_1","rClear","fcost","deathBasal","deathIncrease","deathPeriod")
  #predictors <- c("rDup_1","rUpt_1","rClear","fcost","deathBasal","deathIncrease","deathPeriod") # if rMob=0
  predictors <- c("rMob_1","rDup_1","rClear","fcost","deathBasal","deathIncrease","deathPeriod")
  
  ## Random forest ##
  formula_predictors <- paste(predictors, collapse="+")
  #predict_formula <- as.formula(paste0("m_1 ~", formula_predictors)) # arc distance
  predict_formula <- as.formula(paste0("pab_cond_1 ~", formula_predictors)) # conditional probability co-mobilization
  rf <- randomForest(predict_formula, data=data_predict, mtry=3, importance=TRUE, na.action=na.omit)
  ## Predictor importance (rf) ##
  data_importance <- as.data.frame(importance(rf))
  colnames(data_importance) <- c("IncMSE","IncNodePurity")
  data_importance$predictor <- rownames(data_importance)
  data_importance$treatment <- trt
  ## Store importance ##
  if(is.null(data_importance)){
    data_importance_all <- data_importance
  }else{
    data_importance_all <- rbind(data_importance_all, data_importance)
  }
  
  ## Linear model 1 ##
  lm1 <- lm(predict_formula, data=data_predict)
  
  ## Linear model 2 (interactions) ##
  # predictors_interactions <- c()
  # for(x in 1:(length(predictors)-1)){
  #   for(y in (x+1):length(predictors)){
  #     predictors_interactions <- c(predictors_interactions, paste0(predictors[x],"*",predictors[y]))
  #   }
  # }
  # predictors_interactions <- paste(predictors_interactions, collapse="+")
  # predict_formula_lm2 <- as.formula(paste0("m_1 ~", formula_predictors,"+",predictors_interactions))
  # lm2 <- lm(predict_formula_lm2, data=data_predict)
  
  ## Store lm coefficients ##
  if(is.null(lm1_coef)){
    rows <- c("(intercept)", predictors)
    cols <- c("predictor","trt_1","trt_2","trt_3","trt_4")
    lm1_coef <- as.data.frame(matrix(nrow=length(rows), ncol=length(cols)))
    colnames(lm1_coef) <- cols
    lm1_coef$predictor <- rows
  }
  lm1_coef[, paste0("trt_",trt)] <- lm1$coefficients
} # end treatment loop
rm(trt,rf,data_predict,data_importance,predictors,formula_predictors,predict_formula)
rm(cols,i,nbPop,rows) #predictors_interactions
setDT(data_importance_all)

## Prediction ##
data_predict[, lm1_pred := predict(lm1, data_predict)]
#data_predict[, lm1_resi := m_1 - lm1_pred]
data_predict[, lm1_resi := pab_cond_1 - lm1_pred]
P <- ggplot(data_predict)
P <- P + geom_histogram(aes(x=lm1_resi))
plot(P)

## lm1_coef to long format ##
setDT(lm1_coef)
lm1_coef_long <- melt(
  lm1_coef,
  measure.vars=paste0("trt_",1:4),
  variable.name="treatment",
  value.name="value"
)


# Plots #
figure_path <- "./Figures"
fig_specs <- list(
  "width" = 10,
  "heigh" = 8,
  "dpi" = 200
)

# Importance IncMSE (from random forest) #
P <- ggplot(data_importance_all)
P <- P + ggtitle(paste0(""))
P <- P + geom_bar(aes(x=predictor, y=IncMSE, fill=factor(treatment)), stat="identity", position="dodge",alpha=0.8)
P <- P + scale_x_discrete(name="")
P <- P + scale_y_continuous()
P <- P + scale_fill_hue(labels=c("relaxed","constant","periodic","alternating"))
P <- P + labs(fill="pressure")
P <- P + theme_light()
P <- P + theme(axis.text.x = element_text(angle = 45, vjust = 0.5, hjust=1, size = 10))
plot(P)
ggsave("importance_predictors_IncMSE_.tiff", width=fig_specs$width, height=fig_specs$heigh, units="cm", path=figure_path, dpi=fig_specs$dpi)

# lm coefficients #
P <- ggplot(lm1_coef_long[predictor!="(intercept)"][treatment=="trt_3"])
P <- ggplot(lm1_coef_long[predictor!="(intercept)"])
P <- P + ggtitle(paste0(""))
P <- P + geom_bar(aes(x=predictor, y=value, fill=factor(treatment)), stat="identity", position="dodge")
P <- P + scale_x_discrete(name="")
P <- P + scale_y_continuous(name="coefficient", limits=c(-0.1,0.05)) # m
P <- P + scale_y_continuous(name="coefficient", limits=c(-5E-6,1E-5), breaks=seq(-1E-5,2E-5,5E-6)) # pab_cond
#P <- P + scale_y_continuous(name="coefficient", limits=c(-0.00001,0.00001))
#P <- P + scale_fill_hue(labels=c("relaxed","constant","periodic","alternating"))
P <- P + scale_fill_viridis_d(labels=c("relaxed","constant","periodic","alternating"))
P <- P + labs(fill="pressure")
P <- P + theme_light()
P <- P + theme(axis.text.x = element_text(angle = 45, vjust = 0.5, hjust=1, size = 10))
plot(P)
ggsave("lm1_coef.tiff", width=fig_specs$width, height=fig_specs$heigh, units="cm", path=figure_path, dpi=fig_specs$dpi)

# m=f(rMob) - data_all #
P <- ggplot(data_all[time==simLength][x11_1>1E-6])
P <- P + ggtitle("")
P <- P + geom_point(aes(x=log10(rMob_1), y=m_1, color=factor(treatmentType)), alpha=0.5)
P <- P + scale_x_continuous(name="excision rate (log)")
P <- P + scale_y_continuous(name="arc distance")
P <- P + scale_color_hue(labels=c("relaxed","constant","synchro","altern"))
P <- P + scale_color_viridis_d(labels=c("relaxed","constant","periodic","alternating"))
P <- P + labs(color="pressure")
P <- P + theme_light()
plot(P)
ggsave("m=f(rMob_1)_points.tiff", width=fig_specs$width, height=fig_specs$heigh, units="cm", path=figure_path, dpi=fig_specs$dpi)

# m=f(rMob) - data_all - smooth #
P <- ggplot(data_all[time==simLength][x11_1>1E-6])
P <- P + ggtitle("")
P <- P + geom_smooth(aes(x=log10(rMob_1), y=m_1, color=factor(treatmentType), fill=factor(treatmentType)), size=0.5)
P <- P + scale_x_continuous(name="excision rate (log)")
P <- P + scale_y_continuous(name="arc distance", breaks=seq(0,0.25,0.05))
#P <- P + scale_color_hue(labels=c("relaxed","constant","synchro","altern"))
#P <- P + scale_fill_hue(labels=c("relaxed","constant","synchro","altern"))
P <- P + scale_fill_viridis_d(labels=c("relaxed","constant","periodic","alternating"))
P <- P + scale_color_viridis_d(labels=c("relaxed","constant","periodic","alternating"))
P <- P + labs(fill="pressure",color="pressure")
P <- P + theme_light()
plot(P)
ggsave("m=f(rMob_1)_smooth.tiff", width=fig_specs$width, height=fig_specs$heigh, units="cm", path=figure_path, dpi=fig_specs$dpi)

# m boxplot #
P <- ggplot(data_all[time==simLength][x11_1>1E-6])#[m_1<0.2])
P <- P + ggtitle("arc distance distribution")
P <- P + geom_boxplot(aes(x=factor(treatmentType), y=m_1, color=factor(treatmentType)))
P <- P + scale_x_discrete(name="", labels=c("relaxed","constant","synchro","altern"))
P <- P + scale_y_continuous(name="arc distance")
#P <- P + scale_color_hue(labels=c("relaxed","constant","synchro","altern"))
P <- P + scale_color_viridis_d(labels=c("relaxed","constant","periodic","alternating"))
P <- P + labs(color="pressure")
P <- P + theme_light()
plot(P)
ggsave("m_boxplot.tiff", width=fig_specs$width, height=fig_specs$heigh, units="cm", path=figure_path, dpi=fig_specs$dpi)

# m violin #
P <- ggplot(data_all[time==simLength][x11_1>1E-6])#[m_1<0.2])
P <- P + ggtitle("arc distance distribution")
P <- P + geom_violin(aes(x=factor(treatmentType), y=m_1, fill=factor(treatmentType)), scale="width")
P <- P + scale_x_discrete(name="", labels=c("relaxed","constant","synchro","altern"))
P <- P + scale_y_continuous(name="arc distance")
#P <- P + scale_fill_hue(labels=c("relaxed","constant","synchro","altern"))
P <- P + scale_fill_viridis_d(labels=c("relaxed","constant","periodic","alternating"))
P <- P + labs(fill="pressure")
P <- P + theme_light()
plot(P)
ggsave("m_violin.tiff", width=fig_specs$width, height=fig_specs$heigh, units="cm", path=figure_path, dpi=fig_specs$dpi)

# m density #
P <- ggplot(data_all[time==simLength][x11_1>1E-6])#[m_1<0.2])
P <- P + ggtitle("arc distance distribution")
P <- P + geom_density(aes(x=m_1, fill=factor(treatmentType)), alpha=0.5)
P <- P + scale_y_continuous(name="arc distance")
P <- P + scale_fill_hue(labels=c("relaxed","constant","synchro","altern"))
P <- P + scale_fill_viridis_d(labels=c("relaxed","constant","periodic","alternating"))
P <- P + labs(fill="pressure")
P <- P + theme_light()
plot(P)
ggsave("m_density.tiff", width=fig_specs$width, height=fig_specs$heigh, units="cm", path=figure_path, dpi=fig_specs$dpi)


# mean_upt_x00_dna11_1 #
P <- ggplot(data_all[time==simLength][x11_1>1E-6])#[m_1<0.2])
P <- P + ggtitle("mean_upt_x00_dna11 distribution")
P <- P + geom_violin(aes(x=factor(treatmentType), y=log10(mean_upt_x00_dna11_1), fill=factor(treatmentType)))
P <- P + scale_x_discrete(name="", labels=c("relaxed","constant","synchro","altern"))
P <- P + scale_y_continuous(name="mean_upt_x00_dna11")
#P <- P + scale_color_hue(labels=c("relaxed","constant","synchro","altern"))
P <- P + scale_fill_viridis_d(labels=c("relaxed","constant","periodic","alternating"))
P <- P + labs(fill="treatment")
P <- P + theme_light()
plot(P)
ggsave("mean_upt_x00_dna11_1_distrib.tiff", width=fig_specs$width, height=fig_specs$heigh, units="cm", path=figure_path, dpi=fig_specs$dpi)

# Dna11_1 #
P <- ggplot(data_all[time==simLength][x11_1>1E-6])#[m_1<0.2])
P <- P + ggtitle("mean_dna11 distribution")
P <- P + geom_violin(aes(x=factor(treatmentType), y=log10(Dna11_1), fill=factor(treatmentType)))
P <- P + scale_x_discrete(name="", labels=c("relaxed","constant","synchro","altern"))
P <- P + scale_y_continuous(name="Dna11")#, breaks=seq(-15,15,5), limits=c(-18,10))
P <- P + scale_fill_viridis_d(labels=c("relaxed","constant","periodic","alternating"))
P <- P + labs(color="treatment")
P <- P + theme_light()
plot(P)
ggsave("dna11_distrib.tiff", width=fig_specs$width, height=fig_specs$heigh, units="cm", path=figure_path, dpi=fig_specs$dpi)

# x11_1 #
P <- ggplot(data_all[time==simLength])#[m_1<0.2])
P <- P + ggtitle("mean_x11 distribution")
P <- P + geom_violin(aes(x=factor(treatmentType), y=x11_1, fill=factor(treatmentType)), scale="width")
P <- P + scale_x_discrete(name="", labels=c("relaxed","constant","synchro","altern"))
P <- P + scale_y_continuous(name="x11")
P <- P + scale_fill_viridis_d(labels=c("relaxed","constant","periodic","alternating"))
P <- P + labs(color="treatment")
P <- P + theme_light()
plot(P)
ggsave("x11_distrib.tiff", width=fig_specs$width, height=fig_specs$heigh, units="cm", path=figure_path, dpi=fig_specs$dpi)

# m=f(cost*db*dup) - data_all #
colnames(data_all)
P <- ggplot(data_all[time==simLength][x11_1>1E-6])
P <- P + ggtitle("")
P <- P + geom_point(aes(x=log10(deathBasal*fcost*rDup_1/deathIncrease), y=m_1, color=factor(treatmentType)), alpha=0.7)
P <- P + scale_x_continuous(name="cost*db*dup/di")
P <- P + scale_y_continuous(name="arc distance")
P <- P + scale_color_hue(labels=c("relaxed","constant","synchro","altern"))
P <- P + scale_color_viridis_d(labels=c("relaxed","constant","periodic","alternating"))
P <- P + labs(color="pressure")
P <- P + theme_light()
plot(P)
ggsave("m=f(c,db,dup,di).tiff", width=fig_specs$width, height=fig_specs$heigh, units="cm", path=figure_path, dpi=fig_specs$dpi)

# m=f(cost) - data_all #
colnames(data_all)
P <- ggplot(data_all[time==simLength][x11_1>1E-6])
P <- P + ggtitle("")
P <- P + geom_point(aes(x=log10(fcost), y=m_1, color=factor(treatmentType)), alpha=0.7)
P <- P + scale_x_continuous(name="cost")
P <- P + scale_y_continuous(name="arc distance")
P <- P + scale_color_hue(labels=c("relaxed","constant","synchro","altern"))
P <- P + scale_color_viridis_d(labels=c("relaxed","constant","periodic","alternating"))
P <- P + labs(color="pressure")
P <- P + theme_light()
plot(P)



#### Data_all - 2 species - scenarios ####

## Get data ##
data_all <- out_dt_last_steps[time==simLength]
data_importance_all <- NULL
lm1_coef_all <- NULL
for(scen in 1:3){
  for(trt in 1:4){ # treatment loop
    data_predict <- data_all[time==simLength][x11_2>1E-6][scenario==scen][treatmentType==trt]
    data_predict[, mean_N2N1_ratio := mean_N_2 / (mean_N_1 + mean_N_2)]
    # Log #
    data_predict[, rMob_2 := log10(rMob_2)]
    data_predict[, rDup_2 := log10(rDup_2)]
    data_predict[, rUpt_2 := log10(rUpt_2)]
    data_predict[, rClear := log10(rClear)]
    data_predict[, fcost := log10(fcost)]
    data_predict[, deathBasal := log10(deathBasal)]
    data_predict[, deathIncrease := log10(deathIncrease)]
    data_predict[, deathPeriod := log10(deathPeriod)]
    # Scale #
    data_predict[, rMob_2 := scale(rMob_2,center=TRUE,scale=TRUE)]
    data_predict[, rDup_2 := scale(rDup_2,center=TRUE,scale=TRUE)]
    data_predict[, rUpt_2 := scale(rUpt_2,center=TRUE,scale=TRUE)]
    data_predict[, rClear := scale(rClear,center=TRUE,scale=TRUE)]
    data_predict[, fcost := scale(fcost,center=TRUE,scale=TRUE)]
    data_predict[, deathBasal := scale(deathBasal,center=TRUE,scale=TRUE)]
    data_predict[, deathIncrease := scale(deathIncrease,center=TRUE,scale=TRUE)]
    data_predict[, deathPeriod := scale(deathPeriod,center=TRUE,scale=TRUE)]
    
    ## Predictors ##
    # if(trt==1){
    #   predictors <- c("rMob_2","rDup_2","rUpt_2","rClear","fcost","deathBasal") # no atb
    # }
    # if(trt==2){
    #   predictors <- c("rMob_2","rDup_2","rUpt_2","rClear","fcost","deathBasal","deathIncrease") # cst atb
    # }
    # if(trt==3 | trt==4){
    #   predictors <- c("rMob_2","rDup_2","rUpt_2","rClear","fcost","deathBasal","deathIncrease","deathPeriod") # synchro or altern atb
    # }
    # if(trt==5){
    #   predictors <- c("rMob_2","rDup_2","rUpt_2","rClear","fcost","deathBasal","deathIncrease","deathPeriod","lifePeriod") # uneven synchro atb
    # }
    predictors <- c("rMob_2","rDup_2","rUpt_2","rClear","fcost","deathBasal","deathIncrease","deathPeriod")
    ## Random forest ##
    formula_predictors <- paste(predictors, collapse="+")
    predict_formula <- as.formula(paste0("mean_N2N1_ratio ~", formula_predictors))
    rf <- randomForest(predict_formula, data=data_predict, mtry=3, importance=TRUE, na.action=na.omit)
    ## Predictor importance (rf) ##
    data_importance <- as.data.frame(importance(rf))
    colnames(data_importance) <- c("IncMSE","IncNodePurity")
    data_importance$predictor <- rownames(data_importance)
    data_importance$treatment <- trt
    data_importance$scenario <- scen
    ## Store importance ##
    if(is.null(data_importance)){
      data_importance_all <- data_importance
    }else{
      data_importance_all <- rbind(data_importance_all, data_importance)
    }
    
    ## Linear model 1 ##
    lm1 <- lm(predict_formula, data=data_predict)
    ## Linear model 2 (interactions) ##
    # predictors_interactions <- c()
    # for(x in 1:(length(predictors)-1)){
    #   for(y in (x+1):length(predictors)){
    #     predictors_interactions <- c(predictors_interactions, paste0(predictors[x],"*",predictors[y]))
    #   }
    # }
    # predictors_interactions <- paste(predictors_interactions, collapse="+")
    # predict_formula_lm2 <- as.formula(paste0("m_1 ~", formula_predictors,"+",predictors_interactions))
    # lm2 <- lm(predict_formula_lm2, data=data_predict)
    ## Store lm coefficients ##
    
    rows <- c("(intercept)", predictors)
    cols <- c("scenario","treatment","predictor","coef")
    lm1_coef <- as.data.frame(matrix(nrow=length(rows), ncol=length(cols)))
    colnames(lm1_coef) <- cols
    lm1_coef$scenario <- scen
    lm1_coef$treatment <- trt
    lm1_coef$predictor <- rows
    lm1_coef$coef <- lm1$coefficients
    if(is.null(lm1_coef_all)){
      lm1_coef_all <- lm1_coef
    }else{
      lm1_coef_all <- rbind(lm1_coef_all, lm1_coef)
    }
  } # end treatment loop
} # end scenario loop
rm(scen,trt,rf,data_predict,data_importance,predictors,formula_predictors,predict_formula,lm1_coef,lm1)
rm(cols,rows)


# Figures #
figure_path <- "C:\\Users\\Gabriel Carvalho\\Documents\\Multi-resistance\\Figures"
fig_specs <- list(
  "width" = 11,
  "heigh" = 9,
  "dpi" = 200
)

# Importance IncMSE #
setDT(data_importance_all)
scen <- 1 # scenario
P <- ggplot(data_importance_all[scenario==scen])
P <- P + ggtitle(paste0("scenario ", scen))
P <- P + geom_bar(aes(x=predictor, y=IncMSE, fill=factor(treatment)), stat="identity", position="dodge",alpha=0.8)
P <- P + scale_x_discrete(name="")
P <- P + scale_y_continuous()
#P <- P + scale_fill_hue(labels=c("relaxed","constant","periodic","alterning"))
P <- P + scale_fill_viridis_d(labels=c("relaxed","constant","periodic","alternating"))
P <- P + labs(fill="treatment")
P <- P + theme_light()
P <- P + theme(axis.text.x = element_text(angle = 45, vjust = 0.5, hjust=1, size = 10))
plot(P)
ggsave(paste0("importance_predictors_IncMSE_scen_",scen,".tiff"), width=fig_specs$width, height=fig_specs$heigh, units="cm", path=figure_path, dpi=fig_specs$dpi)

# lm coefficients #
setDT(lm1_coef_all)
scen <- 1 # scenario
P <- ggplot(lm1_coef_all[predictor!="(intercept)"][scenario==scen])
P <- P + ggtitle(paste0("scenario ", scen))
P <- P + geom_bar(aes(x=predictor, y=coef, fill=factor(treatment)), stat="identity", position="dodge")
P <- P + scale_x_discrete(name="")
P <- P + scale_y_continuous(name="coefficient")#, limits=c(-0.1,0.05))
P <- P + scale_fill_viridis_d(labels=c("relaxed","constant","periodic","alternating"))
P <- P + labs(fill="pressure")
P <- P + theme_light()
P <- P + theme(axis.text.x = element_text(angle = 45, vjust = 0.5, hjust=1, size = 10))
plot(P)
ggsave(paste0("lm1_coef_N2pred_scen",scen,".tiff"), width=fig_specs$width, height=fig_specs$heigh, units="cm", path=figure_path, dpi=fig_specs$dpi)

# m boxplot #
P <- ggplot(data_all[time==simLength][x11_2>1E-6])
P <- ggplot(data_all[time==simLength][x11_2>1E-6][scenario==scen])
P <- P + ggtitle("mean_N_2/(mean_N_1+mean_N_2)")
P <- P + geom_boxplot(aes(x=factor(scenario), y=mean_N_2/(mean_N_1+mean_N_2), color=factor(treatmentType)))
P <- P + scale_x_discrete(name="", labels=c("scen.1","scen.2","scen.3"))
P <- P + scale_y_continuous(name="mean_N_2/(mean_N_1+mean_N_2)")
P <- P + scale_color_viridis_d(labels=c("relaxed","constant","periodic","alternating"))
P <- P + labs(color="pressure")
P <- P + theme_light()
plot(P)
ggsave("m_boxplot_scen.tiff", width=fig_specs$width, height=fig_specs$heigh, units="cm", path=figure_path, dpi=fig_specs$dpi)
ggsave(paste0("m_boxplot_scen",scen,".tiff"), width=fig_specs$width, height=fig_specs$heigh, units="cm", path=figure_path, dpi=fig_specs$dpi)

# m violin #
P <- ggplot(data_all[time==simLength][x11_2>1E-6])
P <- ggplot(data_all[time==simLength][x11_2>1E-6][scenario==scen])
P <- P + ggtitle("mean_N_2/(mean_N_1+mean_N_2)")
P <- P + geom_violin(aes(x=factor(scenario), y=mean_N_2/(mean_N_1+mean_N_2), fill=factor(treatmentType)), scale="width")
P <- P + scale_x_discrete(name="", labels=c("scen.1","scen.2","scen.3"))
P <- P + scale_y_continuous(name="mean_N_2/(mean_N_1+mean_N_2)")
P <- P + scale_fill_viridis_d(labels=c("relaxed","constant","periodic","alternating"))
P <- P + labs(fill="pressure")
P <- P + theme_light()
plot(P)
ggsave("m_violin_scen.tiff", width=fig_specs$width, height=fig_specs$heigh, units="cm", path=figure_path, dpi=fig_specs$dpi)
ggsave(paste0("m_violin_scen",scen,".tiff"), width=fig_specs$width, height=fig_specs$heigh, units="cm", path=figure_path, dpi=fig_specs$dpi)

# m density #
P <- ggplot(data_all[time==simLength][x11_2>1E-6][scenario==scen])
P <- P + ggtitle("mean_N_2/(mean_N_1+mean_N_2)")
P <- P + geom_density(aes(x=mean_N_2/(mean_N_1+mean_N_2), fill=factor(treatmentType)), alpha=0.5)
P <- P + scale_y_continuous(name="arc distance")
P <- P + scale_fill_viridis_d(labels=c("relaxed","constant","periodic","alternating"))
P <- P + labs(fill="pressure")
P <- P + theme_light()
plot(P)
ggsave(paste0("m_density_scen",scen,".tiff"), width=fig_specs$width, height=fig_specs$heigh, units="cm", path=figure_path, dpi=fig_specs$dpi)



# N2=f(rMob_2) - data_all #
P <- ggplot(data_all[time==simLength][scenario==scen])
P <- P + ggtitle(paste0("scenario ", scen))
P <- P + geom_point(aes(x=log10(rMob_2), y=mean_N_2/(mean_N_1+mean_N_2), color=factor(treatmentType)), alpha=0.6)
P <- P + scale_x_continuous(name="excision rate (log)")
P <- P + scale_y_continuous(name="N2/(N1+N2)")
P <- P + scale_color_hue(labels=c("relaxed","constant","periodic","alterning"))
P <- P + scale_color_viridis_d(labels=c("relaxed","constant","periodic","alternating"))
P <- P + labs(color="treatment")
P <- P + theme_light()
plot(P)
ggsave(paste0("N2=f(rMob_2)_scen_",scen,".tiff"), width=fig_specs$width, height=fig_specs$heigh, units="cm", path=figure_path, dpi=fig_specs$dpi)

# N2=f(rMob_2) bis #
P <- ggplot(data_all[time==simLength][scenario==scen])
P <- P + ggtitle(paste0("scenario ", scen))
P <- P + geom_smooth(aes(x=log10(rMob_2), y=mean_N_2/(mean_N_1+mean_N_2), fill=factor(treatmentType), color=factor(treatmentType)), alpha=0.4)
P <- P + scale_x_continuous(name="excision rate (log)")
P <- P + scale_y_continuous(name="N2/(N1+N2)")
P <- P + geom_hline(yintercept=0.5, linetype="dashed", color="grey20")
#P <- P + scale_fill_hue(labels=c("relaxed","constant","periodic","alterning"))
P <- P + scale_fill_viridis_d(labels=c("relaxed","constant","periodic","alternating"))
P <- P + scale_color_viridis_d(labels=c("relaxed","constant","periodic","alternating"))
P <- P + labs(fill="pressure",color="pressure")
P <- P + theme_light()
plot(P)
ggsave(paste0("N2=f(rMob_2)_smooth_scen_",scen,".tiff"), width=fig_specs$width, height=fig_specs$heigh, units="cm", path=figure_path, dpi=fig_specs$dpi)

# N2=f(fcost) - data_all #
P <- ggplot(data_all[time==simLength][scenario==scen])
P <- P + ggtitle(paste0("scenario ", scen))
P <- P + geom_point(aes(x=log10(fcost), y=mean_N_2/(mean_N_1+mean_N_2), color=as.character(treatmentType)))
P <- P + scale_x_continuous(name="ARG cost (log)")
P <- P + scale_y_continuous(name="N2/(N1+N2)")
#P <- P + scale_fill_hue(labels=c("constant","synchro","altern"))
P <- P + scale_fill_hue(labels=c("relaxed","constant","periodic","alterning"))
P <- P + labs(color="treatment")
P <- P + theme_light()
plot(P)
ggsave(paste0("N2=f(cost)_scen_",scen,".tiff"), width=fig_specs$width, height=fig_specs$heigh, units="cm", path=figure_path, dpi=fig_specs$dpi)




#### Data_all - 2 patches ####

data_all <- out_dt_last_steps[time==simLength]
colnames(data_all)
setnames(data_all, old=c("rMob","rDup","rUpt"), new=c("rMob_1","rDup_1","rUpt_1"))
## Get data ##
data_importance_all <- NULL
lm1_coef <- NULL
## Arc distance prediction with random forest
data_predict <- data_all[time==simLength][x11_env1_sp1>1E-6]
## Data transformation ##
# Log #
#data_predict[, m_1 := log10(rMob_2)]
data_predict[, rMob_1 := log10(rMob_1)]
data_predict[, rDup_1 := log10(rDup_1)]
data_predict[, rUpt_1 := log10(rUpt_1)]
# data_predict[, rClear := log10(rClear)]
# data_predict[, fcost := log10(fcost)]
# data_predict[, deathBasal := log10(deathBasal)]
# data_predict[, deathIncrease := log10(deathIncrease)]
# data_predict[, G := log10(G)]
# data_predict[, K2 := log10(K2)]
# Scale #
data_predict[, rMob_1 := scale(rMob_1,center=TRUE,scale=TRUE)]
data_predict[, rDup_1 := scale(rDup_1,center=TRUE,scale=TRUE)]
data_predict[, rUpt_1 := scale(rUpt_1,center=TRUE,scale=TRUE)]
data_predict[, rClear := scale(rClear,center=TRUE,scale=TRUE)]
data_predict[, fcost := scale(fcost,center=TRUE,scale=TRUE)]
data_predict[, deathBasal := scale(deathBasal,center=TRUE,scale=TRUE)]
data_predict[, deathIncrease := scale(deathIncrease,center=TRUE,scale=TRUE)]
data_predict[, G := scale(G,center=TRUE,scale=TRUE)]
#data_predict[, K2 := scale(K2,center=TRUE,scale=TRUE)]

## Predictors ##
predictors <- c("rMob_1","rDup_1","rUpt_1","rClear","fcost","deathBasal","deathIncrease","G","K2")
predictors <- c("rMob_1","rDup_1","rClear","fcost","deathBasal","deathIncrease","G")
## Random forest ##
formula_predictors <- paste(predictors, collapse="+")
predict_formula <- as.formula(paste0("m_env1_sp1 ~", formula_predictors))
predict_formula <- as.formula(paste0("pab_cond_env1_sp1 ~", formula_predictors))
rf <- randomForest(predict_formula, data=data_predict, mtry=3, importance=TRUE, na.action=na.omit)
## Predictor importance (rf) ##
data_importance <- as.data.frame(importance(rf))
colnames(data_importance) <- c("IncMSE","IncNodePurity")
data_importance$predictor <- rownames(data_importance)
## Store importance ##
if(is.null(data_importance_all)){
  data_importance_all <- data_importance
}else{
  data_importance_all <- rbind(data_importance_all, data_importance)
}

## Linear model 1 ##
lm1 <- lm(predict_formula, data=data_predict)

## Linear model 2 (interactions) ##
# predictors_interactions <- c()
# for(x in 1:(length(predictors)-1)){
#   for(y in (x+1):length(predictors)){
#     predictors_interactions <- c(predictors_interactions, paste0(predictors[x],"*",predictors[y]))
#   }
# }
# predictors_interactions <- paste(predictors_interactions, collapse="+")
# predict_formula_lm2 <- as.formula(paste0("m_1 ~", formula_predictors,"+",predictors_interactions))
# lm2 <- lm(predict_formula_lm2, data=data_predict)

## Store lm coefficients ##
if(is.null(lm1_coef)){
  rows <- c("(intercept)", predictors)
  cols <- c("predictor","coef")
  lm1_coef <- as.data.frame(matrix(nrow=length(rows), ncol=length(cols)))
  colnames(lm1_coef) <- cols
  lm1_coef$predictor <- rows
  lm1_coef[, "coef"] <- lm1$coefficients
}

rm(rf,data_predict,data_importance,predictors,formula_predictors,predict_formula)
rm(cols,i,nbPop,rows)
setDT(data_importance_all)
setDT(lm1_coef)

# Plots #
figure_path <- "./Figures"
fig_specs <- list(
  "width" = 11,
  "heigh" = 9,
  "dpi" = 200
)

# Importance IncMSE #
P <- ggplot(data_importance_all)
P <- P + ggtitle(paste0(""))
P <- P + geom_bar(aes(x=predictor, y=IncMSE), stat="identity", position="dodge",alpha=0.8)
P <- P + scale_x_discrete(name="")
P <- P + scale_y_continuous()
P <- P + theme_light()
P <- P + theme(axis.text.x = element_text(angle = 45, vjust = 0.5, hjust=1, size = 10))
plot(P)
ggsave("importance_predictors_IncMSE.tiff", width=fig_specs$width, height=fig_specs$heigh, units="cm", path=figure_path, dpi=fig_specs$dpi)

# lm coefficients #
P <- ggplot(lm1_coef[predictor!="(intercept)"])
P <- P + ggtitle(paste0(""))
P <- P + geom_bar(aes(x=predictor, y=coef), stat="identity", position="dodge")
P <- P + scale_x_discrete(name="")
P <- P + scale_y_continuous(name="coefficient")#, limits=c(-0.1,0.05))
P <- P + theme_light()
P <- P + theme(axis.text.x = element_text(angle = 45, vjust = 0.5, hjust=1, size = 10))
plot(P)
ggsave("lm1_coef.tiff", width=fig_specs$width, height=fig_specs$heigh, units="cm", path=figure_path, dpi=fig_specs$dpi)

# m boxplot #
P <- ggplot(data_all[time==simLength][x11_env1_sp1>1E-6])#[m_1<0.2])
P <- P + ggtitle("arc distance distribution")
P <- P + geom_boxplot(aes(x=1, y=m_env1_sp1))
P <- P + scale_y_continuous(name="arc distance")
P <- P + theme_light()
plot(P)
ggsave("m_boxplot.tiff", width=fig_specs$width, height=fig_specs$heigh, units="cm", path=figure_path, dpi=fig_specs$dpi)

# m violin #
P <- ggplot(data_all[time==simLength][x11_env1_sp1>1E-6])#[m_1<0.2])
P <- P + ggtitle("arc distance distribution")
P <- P + geom_violin(aes(x=1, y=m_env1_sp1), scale="width")
P <- P + scale_y_continuous(name="arc distance")
P <- P + theme_light()
plot(P)
ggsave("m_violin.tiff", width=fig_specs$width, height=fig_specs$heigh, units="cm", path=figure_path, dpi=fig_specs$dpi)

# Dna11_1 #
P <- ggplot(data_all[time==simLength][x11_env1_sp1>1E-6])#[m_1<0.2])
P <- P + ggtitle("mean_dna11 distribution")
P <- P + geom_violin(aes(x=1, y=log10(Dna11_env1_sp1)))
P <- P + scale_y_continuous(name="Dna11")#, breaks=seq(-15,15,5), limits=c(-18,10))
P <- P + theme_light()
plot(P)
ggsave("dna11_distrib.tiff", width=fig_specs$width, height=fig_specs$heigh, units="cm", path=figure_path, dpi=fig_specs$dpi)

# x11_1 #
P <- ggplot(data_all[time==simLength])#[m_1<0.2])
P <- P + ggtitle("mean_x11 distribution")
P <- P + geom_violin(aes(x=1, y=x11_env1_sp1))
P <- P + scale_y_continuous(name="x11")
P <- P + theme_light()
plot(P)
ggsave("x11_distrib.tiff", width=fig_specs$width, height=fig_specs$heigh, units="cm", path=figure_path, dpi=fig_specs$dpi)

# m=f(G) #
P <- ggplot(data_all[time==simLength][x11_env1_sp1>1E-6])#[m_1<0.2])
P <- P + ggtitle("m=f(G)")
P <- P + geom_point(aes(x=log10(G), y=m_env1_sp1))
P <- P + scale_y_continuous(name="m")
P <- P + theme_light()
plot(P)
ggsave("m=f(G).tiff", width=fig_specs$width, height=fig_specs$heigh, units="cm", path=figure_path, dpi=fig_specs$dpi)

# m=f(mob) #
P <- ggplot(data_all[time==simLength][x11_env1_sp1>1E-6])#[m_1<0.2])
P <- P + ggtitle("m=f(mob)")
P <- P + geom_point(aes(x=log10(rMob_1), y=m_env1_sp1))
P <- P + scale_y_continuous(name="m")
P <- P + theme_light()
plot(P)
ggsave("m=f(mob).tiff", width=fig_specs$width, height=fig_specs$heigh, units="cm", path=figure_path, dpi=fig_specs$dpi)

# m=f(dup) #
P <- ggplot(data_all[time==simLength][x11_env1_sp1>1E-6])#[m_1<0.2])
P <- P + ggtitle("m=f(dup)")
P <- P + geom_point(aes(x=log10(rDup_1), y=m_env1_sp1))
P <- P + scale_y_continuous(name="m")
P <- P + theme_light()
plot(P)
ggsave("m=f(dup).tiff", width=fig_specs$width, height=fig_specs$heigh, units="cm", path=figure_path, dpi=fig_specs$dpi)

# m=f(cost) #
P <- ggplot(data_all[time==simLength][x11_env1_sp1>1E-6])#[m_1<0.2])
P <- P + ggtitle("m=f(cost)")
P <- P + geom_point(aes(x=log10(fcost), y=m_env1_sp1))
P <- P + scale_y_continuous(name="m")
P <- P + theme_light()
plot(P)
ggsave("m=f(cost).tiff", width=fig_specs$width, height=fig_specs$heigh, units="cm", path=figure_path, dpi=fig_specs$dpi)



#### Plots ####
figure_path <- "~/Documents/Article-fabric/Article - Multi-resistance/R_model_JP/Figures"
colnames(data_predict)

# Importance IncMSE #
P <- ggplot(data_importance)
P <- P + geom_bar(aes(x=predictor, y=IncMSE), stat="identity")
P <- P + scale_x_discrete(name="")
P <- P + scale_y_continuous()
P <- P + theme_light()
P <- P + theme(axis.text.x = element_text(angle = 45, vjust = 0.5, hjust=1, size = 10))
plot(P)
ggsave("importance_predictors_IncMSE.tiff", width=9, height=8, units="cm", path=figure_path, dpi=150)

# Importance IncNodePurity #
P <- ggplot(data_importance)
P <- P + geom_bar(aes(x=predictor, y=IncNodePurity), stat="identity")
P <- P + scale_x_discrete(name="")
P <- P + scale_y_continuous()
P <- P + theme_light()
P <- P + theme(axis.text.x = element_text(angle = 45, vjust = 0.5, hjust=1))
plot(P)
ggsave("importance_predictors_IncNodePurity.tiff", width=9, height=8, units="cm", path=figure_path, dpi=150)

# m=f(rMob) - data_predict #
P <- ggplot(data_predict)
P <- P + geom_point(aes(x=(10^rMob)^2, y=10^m_1))
P <- P + scale_x_continuous(name="(DNA excision rate)^2")
P <- P + scale_y_continuous(name="mean arc distance")
P <- P + theme_light()
plot(P)
ggsave("m=f(rMob).tiff", width=10, height=8, units="cm", path=figure_path, dpi=150)

# m=f(rMob) - out_dt #
colnames(out_dt)
P <- ggplot(out_dt[time==simLength][x11_1>1E-6])
P <- P + geom_point(aes(x=rMob^2, y=m_1, color=WTin))
P <- P + scale_color_viridis_c()
P <- P + scale_x_continuous(name="(DNA excision rate)^2")
P <- P + scale_y_continuous(name="mean arc distance")
P <- P + theme_light()
plot(P)
ggsave("m=f(rMob).tiff", width=10, height=8, units="cm", path=figure_path, dpi=150)

P <- ggplot(out_dt[time==simLength][x11_1>1E-6])
P <- P + geom_point(aes(x=rMob, y=m_1, color=as.character(WTin)))
P <- P + geom_line(aes(x=rMob, y=m_1, color=as.character(WTin)))
P <- P + scale_x_continuous(name="DNA excision rate")
P <- P + scale_y_continuous(name="mean arc distance")
P <- P + theme_light()
plot(P)


## 2 species ##
# m=f(rMob_2) - data_predict #
P <- ggplot(data_predict)
P <- P + geom_point(aes(x=log10(rMob_2), y=mean_N2N1_ratio))
P <- P + theme_light()
plot(P)
ggsave("m=f(rMob_2).tiff", width=10, height=8, units="cm", path=figure_path, dpi=150)

# m=f(deathBasal) - data_predict #
P <- ggplot(data_predict)
P <- P + geom_point(aes(x=deathBasal, y=mean_N2N1_ratio))
P <- P + theme_light()
plot(P)
ggsave("m=f(deathBasal).tiff", width=10, height=8, units="cm", path=figure_path, dpi=150)

# m=f(fcost) - data_predict #
P <- ggplot(data_predict)
P <- P + geom_point(aes(x=log10(fcost), y=mean_N2N1_ratio))
P <- P + theme_light()
plot(P)
ggsave("m=f(fcost).tiff", width=10, height=8, units="cm", path=figure_path, dpi=150)

# m=f(lifePeriod) - data_predict #
P <- ggplot(data_predict)
P <- P + geom_point(aes(x=lifePeriod, y=mean_N2N1_ratio))
P <- P + theme_light()
plot(P)
ggsave("m=f(lifePeriod).tiff", width=10, height=8, units="cm", path=figure_path, dpi=150)


