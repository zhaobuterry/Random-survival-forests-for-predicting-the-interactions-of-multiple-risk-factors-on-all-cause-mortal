library("randomForestSRC")
library("ggRandomForests")
library("tidyverse")
library(survival)
library(ggplot2)
library(cowplot)
library(dplyr)
library(MASS)
library(LTRCtrees)
library(rpart.plot)
library(scales)

##############################
#### Load the NHANES DATA ####
##############################

Data <- nhanes_merged_dataset %>%
  # Define new columns with more legible codenames
  mutate(mortality_status = MORTSTAT
         , time_to_death = PERMTH_INT
         , weights = WTINT2YR
         , cluster = paste(nhanes_merged_dataset$SDMVPSU, nhanes_merged_dataset$SDMVSTRA)
         , cycles = SDDSRVYR
         , age = RIDAGEYR
         , gender = relevel(factor(if_else(RIAGENDR == 1, "_male", "_female"))
                            , ref = "_male")
         , race = relevel(factor(case_when(RIDRETH1 == 1 ~ "_mexican_american"
                                           , RIDRETH1 == 2 ~ "_other_hispanic"
                                           , RIDRETH1 == 3 ~ "_whites"
                                           , RIDRETH1 == 4 ~ "_blacks"
                                           , RIDRETH1 == 5 ~ "_other" ))
                          , ref = "_whites")) %>%
  # Select the pertinent variables after removing highly correlated indicators
  dplyr::select("mortality_status", "time_to_death",  "age", "gender", "race","BMXBMI","BPXPLS","VNRFPI","VNEGFR","VNHOMAIR","LBXCRP","LBXTC","LBXWBCSI","LBXSAL","LBXSTR","LBXGLU","LBXSAPSI","LBDHDD","LBXSBU","VNAVEBPXSY","VNAVEBPXDI") %>%
  # Exclude participants who have missing data
  na.omit(.)%>%
  # Exclude participants with no follow-up data
  filter(time_to_death != 0)

# Create the variable using age as time scale
Data$End <- Data$age + Data$time_to_death/12
Data$Start <- 0

##############################
##### Variable Selection #####
##############################

# Costruct random survival forest model using age as time sclae
set.seed(1)
RSF <- rfsrc(Surv(End, mortality_status) ~ ., data = train,
                       ntree = 500, mtry=5, nodesize = 2000, importance = T)

# Select needed variables based on VIMP
plot(RSF)

############################################################
##### Visualization of Risk Groups Using Survival Tree #####
############################################################

# Fit survival tree model
set.seed(2)
LTRCART.obj_VNEGFR <- LTRCART(Surv(age,End, mortality_status) ~ gender+ race +VNEGFR, 
                              data = Data, control = rpart::rpart.control(cp = 0.00001, minbucket = 1000))
#Visulaize Survival Tree
rpart.plot(LTRCART.obj_VNEGFR, digits=-2, tweak = 1.2, roundint=FALSE)

############################################################
##### Determination of HR Using Random Survival Forest #####
############################################################

# Construct Random Survival Forest model and calculate the Ensemble Mortality
set.seed(1)
RSF.obj_5 <- rfsrc(Surv(time_to_death, mortality_status) ~ age+ gender+ VNEGFR+ VNRFPI+ LBXGLU+ LBXWBCSI, 
                           data = Data, ntree = 500, mtry=3, nodesize = 200, importance = T)
RSF.obj_5.pred <- predict.rfsrc(RSF.obj_5, Data)
pred.scores <- RSF.obj_5.pred$predicted
pred.scores.norm <- (pred.scores-min(pred.scores))/(max(pred.scores)-min(pred.scores))
cv <- coxph(Surv(time_to_death, mortality_status) ~ pred.scores.norm, 
            data = cbind(Data[,c(1,2)], pred.scores.norm))

# Calculate the hazard ratio for the reference person
median_person <- Data[1,]
median_person$VNEGFR <- quantile(Data$VNEGFR, .5)
median_person$LBXWBCSI <- quantile(Data$LBXWBCSI, .5)
median_person$LBXGLU <- quantile(Data$LBXGLU, .5)
median_person$VNRFPI <- quantile(Data$VNRFPI, .5)
median_person$age <- quantile(Data$age, .5)

median_person.pred <- predict.rfsrc(RSF.obj_5, median_person)
pred.scores.median_person <- median_person.pred$predicted
median_person.norm <- (pred.scores.median_person - min(pred.scores))/(max(pred.scores)-min(pred.scores))
HR.median_person <- cv$coefficients*median_person.norm


# Design the needed value grids
VNEGFR <- seq(from = quantile(Data$VNEGFR, .05), 
              to = quantile(Data$VNEGFR, .95), 
              length.out = 2500)
LBXWBCSI_reference <- rep(quantile(Data$LBXWBCSI, .5), 2500)
LBXGLU_reference <- rep(quantile(Data$LBXGLU, .5), 2500)
VNRFPI_reference <- rep(quantile(Data$VNRFPI, .5), 2500)
gender.list <- unique(Data$gender)
age.list <- c(19, 46, 60, 70, 80)
```

# Calculate the HR on the designed grids
Hazard <- vector()
for (i in 1:5){
  for (j in 1:2){
    Test <- as.data.frame(cbind(VNEGFR,LBXWBCSI_reference,
                                LBXGLU_reference,VNRFPI_reference))
    Test$age <- rep(age.list[i],2500)
    Test$gender <- rep(gender.list[j],2500)
    colnames(Test) <- c("VNEGFR", "LBXWBCSI", "LBXGLU", "VNRFPI", "age","gender")
    levels(Test$gender) <- levels(Data$gender)
    RSF.obj_5.test.pred <- predict.rfsrc(RSF.obj_5, Test)
    pred.scores.test <- RSF.obj_5.test.pred$predicted
    pred.scores.norm.test <- (pred.scores.test-min(pred.scores))/(max(pred.scores)-min(pred.scores))
    HR.test <- cv$coefficients*pred.scores.norm.test
    Hazard_norm <- HR.test-HR.median_person
    Hazard <- cbind(Hazard,Hazard_norm)
  }
}

# Plot the hazard ratio on the test set
Hazard_results <- as.data.frame(cbind(VNEGFR,Hazard))
colnames(Hazard_results) <- c("VNEGFR", "19_male", "19_female", "46_male", 
                              "46_female", "60_male", "60_female", "70_male", 
                              "70_female", "80_male", "80_female")

df <- Hazard_results %>%
  gather(key = "variable", value = "value", -VNEGFR)
ggplot(df, aes(x = VNEGFR, y = value)) + 
  geom_line(aes(color = variable, linetype = variable)) + 
  scale_color_manual(values = c("darkgray","darkgray", "darkred","darkred",
                                "blue","blue","darkgreen","darkgreen",
                                "black","black"))+
  scale_linetype_manual(values=c("dashed","solid","dashed","solid", "dashed","solid",
                                 "dashed","solid", "dashed","solid"))+
  theme_bw() +
  theme(text = element_text(size = 15))  +  
  theme(axis.text = element_text(size = 15)) + 
  geom_hline(yintercept = 0, linetype="twodash", color = "red") +
  geom_vline(xintercept = quantile(Data$VNEGFR, .5), linetype="twodash", color = "blue") +
  xlab("Estimated Glomerular Filtration Rate (mL/min/1.73 m2)") + ylab("log Hazard ratios")


##############################################
##### Visualization of HR Using Heatmaps #####
##############################################

# Set the value grids
VNEGFR <- seq(from = quantile(Data$VNEGFR, .05), 
              to = quantile(Data$VNEGFR, .95), length.out = 50)
LBXGLU <- seq(from = quantile(Data$LBXGLU, .05), 
              to = quantile(Data$LBXGLU, .95), length.out = 50)
Test.base <- expand.grid(VNEGFR = VNEGFR, LBXGLU = LBXGLU)
LBXWBCSI.list <- c(quantile(Data$LBXWBCSI, .05), 
                   quantile(Data$LBXWBCSI, .5), 
                   quantile(Data$LBXWBCSI, .95))
VNRFPI.list <- c(quantile(Data$VNRFPI, .05), 
                 quantile(Data$VNRFPI, .5), 
                 quantile(Data$VNRFPI, .95))

# Plot the hazard ratio using heatmaps
plot_list <- list()
count <- 0
for (i in 1:3){
  for (j in 1:3){
    count <- count+1
    age <- rep(age.list[2], 2500)
    VNRFPI <- rep(VNRFPI.list[i], 2500)
    LBXWBCSI <- rep(LBXWBCSI.list[j], 2500)
    Test <- cbind(VNRFPI, Test.base, LBXWBCSI, age)
    colnames(Test) <- c("VNRFPI", "VNEGFR", "LBXGLU", "LBXWBCSI", "age")
    Test$gender=rep(gender.list[1],2500)
    levels(Test$gender) <- levels(Data$gender)
    
    RSF.obj_5.combine.pred.test <- predict.rfsrc(RSF.obj_5, Test)
    all.pred.scores.test <- RSF.obj_5.combine.pred.test$predicted
    all.pred.scores.norm.test <- (all.pred.scores.test - min(pred.scores))/(max(pred.scores)-min(pred.scores))
    HR.all.test <- cv$coefficients*all.pred.scores.norm.test
    Hazard_all_norm <- HR.all.test-HR.median_person
    
    p <- ggplot(Test, aes(VNEGFR, LBXGLU, fill=Hazard_all_norm )) + 
      geom_tile() +
      scale_fill_gradient2(low="green", mid="yellow", high="red", 
                           midpoint=0.4, 
                           breaks=seq(-0.4,1.6,0.4),
                           limits=c(-0.4,1.6))+
      theme(legend.position="none")
    plot_list[[count]] = ggplotGrob(p)
  }
}

cowplot::plot_grid(plotlist = plot_list)
