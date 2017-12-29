library(car)
library(glmnet)
library(pROC)

full_data <- read.csv("diabetic.data.csv", header=T)
cleaned_data <- read.csv("readmission.csv", header=T)
cleaned_data$readmitted <- as.factor(ifelse(cleaned_data$readmitted == '<30', 1, 0))

cleaned_data[which(cleaned_data$race == "?"), 3] <- NA
cleaned_data[which(cleaned_data$gender == "Unknown/Invalid"), 4] <- NA
cleaned_data <- na.omit(cleaned_data)
cleaned_data$race <- factor(cleaned_data$race)
cleaned_data$gender <- factor(cleaned_data$gender)
cleaned_data$not_unique <- as.factor(as.integer(duplicated(cleaned_data$patient_nbr)))

X <- model.matrix(readmitted ~ . -encounter_id -patient_nbr, cleaned_data)[,-1]
Y <- cleaned_data$readmitted

set.seed(10)
fit1.cv <- cv.glmnet(X, Y, alpha=1, family="binomial", nfolds = 10, type.measure = "deviance")
coef1.1se <- coef(fit1.cv, s = "lambda.1se")  
coef1.1se <- coef1.1se[which(coef1.1se != 0), ] 
fit1.1se.glm.insig <- glm(readmitted ~ time_in_hospital + num_medications + number_emergency + 
                     number_inpatient + number_diagnoses + metformin + insulin + 
                     diabetesMed + disch_disp_modified + age_mod + diag1_mod + 
                     diag2_mod + diag3_mod + not_unique, family="binomial", data=cleaned_data)
fit1.1se.glm.Anova.insig <- Anova(fit1.1se.glm.insig)
fit1.1se.glm <- glm(readmitted ~ num_medications + number_emergency + 
                      number_inpatient + number_diagnoses + metformin + insulin + 
                      diabetesMed + disch_disp_modified + age_mod + diag1_mod + 
                      diag2_mod + diag3_mod + not_unique, family="binomial", data=cleaned_data)
fit1.1se.glm.Anova <- Anova(fit1.1se.glm)
fit1.1se.glm.roc <- roc(cleaned_data$readmitted, fit1.1se.glm$fitted.values, plot=T, col="blue")
saveRDS(fit1.cv, 'fit1_cv.rds')
saveRDS(fit1.1se.glm.insig, 'fit1_1se_glm_insig.rds')
saveRDS(fit1.1se.glm.Anova.insig, 'fit1_1se_glm_Anova_insig.rds')
saveRDS(fit1.1se.glm, 'fit1_1se_glm.rds')
saveRDS(fit1.1se.glm.Anova, 'fit1_1se_glm_Anova.rds')
saveRDS(fit1.1se.glm.roc, 'fit1_1se_glm_roc.rds')

set.seed(11)
fit2.cv <- cv.glmnet(X, Y, alpha=1, family="binomial", nfolds = 10, type.measure = "class")
coef2.min <- coef(fit2.cv, s = "lambda.min")  
coef2.min <- coef2.min[which(coef2.min != 0), ] 
fit2.min.glm <- glm(readmitted ~ time_in_hospital + number_inpatient + number_diagnoses + 
                      disch_disp_modified + not_unique, family="binomial", data=cleaned_data)
fit2.min.glm.Anova <- Anova(fit2.min.glm)
fit2.min.glm.roc <- roc(cleaned_data$readmitted, fit2.min.glm$fitted.values, plot=T, col="blue")
saveRDS(fit2.cv, 'fit2_cv.rds')
saveRDS(fit2.min.glm, 'fit2_min_glm.rds')
saveRDS(fit2.min.glm.Anova, 'fit2_min_glm_Anova.rds')
saveRDS(fit2.min.glm.roc, 'fit2_min_glm_roc.rds')

set.seed(12)
fit3.cv <- cv.glmnet(X, Y, alpha=1, family="binomial", nfolds = 10, type.measure = "auc")
coef3.1se <- coef(fit3.cv, s = "lambda.1se")  
coef3.1se <- coef3.1se[which(coef3.1se != 0), ] 
fit3.1se.glm.insig <- glm(readmitted ~ time_in_hospital + num_medications + number_emergency + 
                      number_inpatient + number_diagnoses + A1Cresult + metformin + insulin + 
                      diabetesMed + disch_disp_modified + age_mod + diag1_mod + diag2_mod + 
                      diag3_mod + not_unique, family="binomial", data=cleaned_data)
fit3.1se.glm.Anova.insig <- Anova(fit3.1se.glm.insig)
fit3.1se.glm <- fit1.1se.glm
fit3.1se.glm.Anova <- fit1.1se.glm.Anova
fit3.1se.glm.roc <- fit1.1se.glm.roc
saveRDS(fit3.cv, 'fit3_cv.rds')
saveRDS(fit3.1se.glm.insig, 'fit3_1se_glm_insig.rds')
saveRDS(fit3.1se.glm.Anova.insig, 'fit3_1se_glm_Anova_insig.rds')
saveRDS(fit3.1se.glm, 'fit3_1se_glm.rds')
saveRDS(fit3.1se.glm.Anova, 'fit3_1se_glm_Anova.rds')
saveRDS(fit3.1se.glm.roc, 'fit3_1se_glm_roc.rds')

set.seed(13)
fit4.cv <- cv.glmnet(X, Y, alpha=0.99, family="binomial", nfolds = 10, type.measure = "deviance")
coef4.1se <- coef(fit4.cv, s = "lambda.1se")  
coef4.1se <- coef4.1se[which(coef4.1se != 0), ]
fit4.1se.glm.insig <- glm(readmitted ~ time_in_hospital + num_medications + number_emergency + 
                      number_inpatient + number_diagnoses + insulin + diabetesMed + disch_disp_modified + 
                      age_mod + diag1_mod + diag3_mod + not_unique, family="binomial", data=cleaned_data)
fit4.1se.glm.Anova.insig <- Anova(fit4.1se.glm.insig)
fit4.1se.glm <- glm(readmitted ~ number_emergency + number_inpatient + number_diagnoses + 
                            insulin + diabetesMed + disch_disp_modified + age_mod + diag1_mod + 
                            diag3_mod + not_unique, family="binomial", data=cleaned_data)
fit4.1se.glm.Anova <- Anova(fit4.1se.glm)
fit4.1se.glm.roc <- roc(cleaned_data$readmitted, fit4.1se.glm$fitted.values, plot=T, col="blue")
saveRDS(fit4.cv, 'fit4_cv.rds')
saveRDS(fit4.1se.glm.insig, 'fit4_1se_glm_insig.rds')
saveRDS(fit4.1se.glm.Anova.insig, 'fit4_1se_glm_Anova_insig.rds')
saveRDS(fit4.1se.glm, 'fit4_1se_glm.rds')
saveRDS(fit4.1se.glm.Anova, 'fit4_1se_glm_Anova.rds')
saveRDS(fit4.1se.glm.roc, 'fit4_1se_glm_roc.rds')
