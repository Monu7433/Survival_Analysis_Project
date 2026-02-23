# ============================================================
# SURVIVAL ANALYSIS OF BREAST CANCER PATIENTS (METABRIC)
# FINAL THESIS VERSION
# ============================================================



# ============================================================
# 1. LOAD LIBRARIES
# ============================================================

library(janitor)

library(dplyr)

library(survival)

library(survminer)

library(ggplot2)

library(tableone)

library(corrplot)



# ============================================================
# 2. IMPORT DATASET
# ============================================================

df <- read.csv(file.choose())



# ============================================================
# 3. DATA CLEANING
# ============================================================

df <- clean_names(df)

df <- df %>%
  mutate(across(where(is.character), ~ na_if(.x, "")))


df_main <- df %>%
  select(
    
    patient_id,
    
    age_at_diagnosis,
    
    tumor_size,
    
    tumor_stage,
    
    lymph_nodes_examined_positive,
    
    neoplasm_histologic_grade,
    
    er_status,
    
    pr_status,
    
    her2_status,
    
    chemotherapy,
    
    hormone_therapy,
    
    radio_therapy,
    
    overall_survival_months,
    
    overall_survival_status
    
  )



df_main <- df_main %>%
  filter(
    
    !is.na(overall_survival_months),
    
    !is.na(overall_survival_status)
    
  )



df_main <- df_main %>%
  mutate(
    
    overall_survival_status =
      ifelse(overall_survival_status=="Deceased",1,0)
    
  )



# ============================================================
# 4. CONVERT TO FACTOR
# ============================================================

df_main$er_status <- factor(df_main$er_status)

df_main$pr_status <- factor(df_main$pr_status)

df_main$her2_status <- factor(df_main$her2_status)

df_main$chemotherapy <- factor(df_main$chemotherapy)

df_main$hormone_therapy <- factor(df_main$hormone_therapy)

df_main$radio_therapy <- factor(df_main$radio_therapy)

df_main$tumor_stage <- factor(df_main$tumor_stage)

df_main$neoplasm_histologic_grade <- factor(df_main$neoplasm_histologic_grade)



# ============================================================
# 5. CREATE GROUP VARIABLES
# ============================================================

df_main$nodes_group <- cut(
  
  df_main$lymph_nodes_examined_positive,
  
  breaks=c(-1,0,3,10,100),
  
  labels=c("0","1-3","4-10",">10")
  
)



df_main$tsize_group <- cut(
  
  df_main$tumor_size,
  
  breaks=c(0,20,50,200),
  
  labels=c("<20","20-50",">50")
  
)



# ============================================================
# 6. EDA SUMMARY
# ============================================================

summary(df_main)

colSums(is.na(df_main))


CreateTableOne(
  
  vars=c(
    
    "age_at_diagnosis",
    
    "tumor_size",
    
    "lymph_nodes_examined_positive",
    
    "er_status",
    
    "pr_status",
    
    "her2_status",
    
    "tumor_stage",
    
    "neoplasm_histologic_grade",
    
    "chemotherapy",
    
    "hormone_therapy",
    
    "radio_therapy"
    
  ),
  
  data=df_main
  
)



# ============================================================
# 7. CORRELATION
# ============================================================

cor_matrix <- cor(
  
  df_main[,c(
    
    "age_at_diagnosis",
    
    "tumor_size",
    
    "lymph_nodes_examined_positive",
    
    "overall_survival_months")],
  
  use="complete.obs")



corrplot(cor_matrix)



# ============================================================
# 8. CREATE SURVIVAL OBJECT
# ============================================================

survival_object <- Surv(
  
  df_main$overall_survival_months,
  
  df_main$overall_survival_status)



# ============================================================
# 9. KAPLAN MEIER
# ============================================================

km <- survfit(survival_object~1,data=df_main)
km
km_er <- survfit(survival_object~er_status,data=df_main)

km_pr <- survfit(survival_object~pr_status,data=df_main)

km_her2 <- survfit(survival_object~her2_status,data=df_main)

km_stage <- survfit(survival_object~tumor_stage,data=df_main)

km_grade <- survfit(survival_object~neoplasm_histologic_grade,data=df_main)

km_nodes <- survfit(survival_object~nodes_group,data=df_main)

km_size <- survfit(survival_object~tsize_group,data=df_main)

km_chemo <- survfit(survival_object~chemotherapy,data=df_main)

km_hormone <- survfit(survival_object~hormone_therapy,data=df_main)

km_radio <- survfit(survival_object~radio_therapy,data=df_main)



# ============================================================
# 10. LOG RANK TEST
# ============================================================

survdiff(survival_object~er_status,data=df_main)

survdiff(survival_object~pr_status,data=df_main)

survdiff(survival_object~her2_status,data=df_main)

survdiff(survival_object~tumor_stage,data=df_main)

survdiff(survival_object~nodes_group,data=df_main)



# ============================================================
# 11. UNIVARIABLE COX
# ============================================================

summary(coxph(survival_object~age_at_diagnosis,data=df_main))

summary(coxph(survival_object~tumor_size,data=df_main))

summary(coxph(survival_object~lymph_nodes_examined_positive,data=df_main))

summary(coxph(survival_object~tumor_stage,data=df_main))

summary(coxph(survival_object~neoplasm_histologic_grade,data=df_main))

summary(coxph(survival_object~er_status,data=df_main))

summary(coxph(survival_object~pr_status,data=df_main))

summary(coxph(survival_object~her2_status,data=df_main))

summary(coxph(survival_object~chemotherapy,data=df_main))

summary(coxph(survival_object~hormone_therapy,data=df_main))

summary(coxph(survival_object~radio_therapy,data=df_main))



# ============================================================
# 12. MULTIVARIABLE COX
# ============================================================



#Model

cox_multi <- coxph(
  
  survival_object~
    
    age_at_diagnosis+
    
    tumor_size+
    
    lymph_nodes_examined_positive+
    
    neoplasm_histologic_grade+
    
    tumor_stage+
    
    er_status+
    
    pr_status+
    
    her2_status+
    
    chemotherapy+
    
    hormone_therapy+
    
    radio_therapy,
  
  data=df_main)



summary(cox_multi)


# =====================================
# PH ASSUMPTION TEST
# =====================================

ph_test <- cox.zph(cox_multi)

ph_test

# Plot for visual check
plot(ph_test)



# ============================================================
# 13. EXPORT ALL PLOTS INTO PDF
# ============================================================

pdf("FINAL_METABRIC_ANALYSIS.pdf", width=8, height=6)



# EDA

print(ggplot(df_main, aes(age_at_diagnosis))+geom_histogram())

print(ggplot(df_main, aes(tumor_size))+geom_histogram())

print(ggplot(df_main, aes(lymph_nodes_examined_positive))+geom_histogram())

print(ggplot(df_main, aes(er_status))+geom_bar())



# Correlation

corrplot(cor_matrix)



# Survival plots

print(ggsurvplot(km)$plot)

print(ggsurvplot(km_er)$plot)

print(ggsurvplot(km_pr)$plot)

print(ggsurvplot(km_her2)$plot)

print(ggsurvplot(km_stage)$plot)

print(ggsurvplot(km_nodes)$plot)

print(ggsurvplot(km_size)$plot)

print(ggsurvplot(km_chemo)$plot)

print(ggsurvplot(km_hormone)$plot)

print(ggsurvplot(km_radio)$plot)



dev.off()



# ============================================================
# END OF PROJECT
# ============================================================