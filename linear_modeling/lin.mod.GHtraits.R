library(lme4)
library(ggplot2)
library(lmerTest)

#data <- read.table("GH_IL_data_summary.txt", header=TRUE)
data <- read.delim("~SuppDataset4_SummaryTraits.txt")
names(data)
summary(data)

# -----

# First trait: pct_iso.alkanes
# Take a look at it...
qplot(pct_iso.alkanes,data=data,geom="histogram")

# The dataset is skewed; look at a qq plot
qqnorm(data$pct_iso.alkanes)
qqline(data$pct_iso.alkanes)

# There is some deviation; let's try some transformations
sq_trait <- sqrt(data$pct_iso.alkanes)
log_trait <- log(data$pct_iso.alkanes+1)
rec_trait <- 1/(data$pct_iso.alkanes)
asin_trait <- asin(data$pct_iso.alkanes/100)

# Look at histograms and qq plots for these transformed data sets
p <- ggplot(data, aes(x=sq_trait))
p + geom_histogram()

p <- ggplot(data, aes(x=log_trait))
p + geom_histogram()

# Look at qq plots of the transformed data sets
qqnorm(sq_trait)
qqline(sq_trait)

qqnorm(log_trait)
qqline(log_trait)

qqnorm(rec_trait)
qqline(rec_trait)

qqnorm(asin_trait)
qqline(asin_trait)

# Let's formally test deviation of transformations from a true normal distribution, using Shapiro test
shapiro.test(data$pct_iso.alkanes)
shapiro.test(sq_trait)
shapiro.test(log_trait)
shapiro.test(rec_trait)
shapiro.test(asin_trait)

# The log transformation shows the least deviation from a normal distribution, so we'll use it

# Model our trait using a mixed effect linear model -- no need for model selection
pct_iso.alkanes_model <- lmer(log_trait ~ il + (1|Plot), data=data)

# Heritability values
h2.trans <- lmer(log_trait ~ (1|il) + (1|Plot), data=data) #Transformed
h2.trans
h2.model <- lmer(data$pct_iso.alkanes ~ (1|il) + (1|Plot), data=data) #Not transformed
h2.model

# Create function "modelcheck"
modelcheck <- function(model,h=8,w=10.5) {
  rs <- residuals(model)
  fv <- fitted(model)
  quartz(h=h,w=w)
  plot(rs~fv)
  quartz(h=h,w=w)
  plot(sqrt(abs(rs))~fv)
  quartz(h=h,w=w)
  qqnorm(rs)
  qqline(rs)
}

modelcheck(pct_iso.alkanes_model)

# Everything looks good! Now, find the p values for the ILs deviating from the cvM82 parent
summary(pct_iso.alkanes_model)

# There appear to be many significant ILs. Write out data.
summary <- summary(pct_iso.alkanes_model)$coef
write.table(summary, file="pct_iso.alkanes_model_summary.txt")
# -------------------------------------------------------
# Next trait: n-alkane CPI
qplot(n.alkane_CPI,data=data,geom="histogram")

# QQ plot
qqnorm(data$n.alkane_CPI)
qqline(data$n.alkane_CPI)

# Data transformations
sq_nalkCPI <- sqrt(data$n.alkane_CPI)
log_nalkCPI <- log(data$n.alkane_CPI+1)
rec_nalkCPI <- 1/(data$n.alkane_CPI)
asin_nalkCPI <- asin(data$n.alkane_CPI/100)

# QQ plots
qqnorm(sq_nalkCPI)
qqline(sq_nalkCPI)

qqnorm(log_nalkCPI)
qqline(log_nalkCPI)

qqnorm(rec_nalkCPI)
qqline(rec_nalkCPI)

qqnorm(asin_nalkCPI)
qqline(asin_nalkCPI)

# Formally test deviation from normal distribution
shapiro.test(data$n.alkane_CPI)
shapiro.test(sq_nalkCPI)
shapiro.test(log_nalkCPI)
shapiro.test(rec_nalkCPI)
shapiro.test(asin_nalkCPI)

# Use log transformation

nalkCPI_model <- lmer(log_nalkCPI ~ il + (1|Plot), data=data)

modelcheck(nalkCPI_model)

# Everything looks good! Now, find the p values for the ILs deviating from the cvM82 parent
summary(nalkCPI_model)

# There appear to be many significant ILs. Write out data.
summary <- summary(nalkCPI_model)$coef
write.table(summary, file="nalkaneCPI_model_summary.txt")

# Heritability values
h2.trans <- lmer(log_nalkCPI ~ (1|il) + (1|Plot), data=data) #Transformed
h2.trans
h2.model <- lmer(data$n.alkane_CPI ~ (1|il) + (1|Plot), data=data) #Not transformed
h2.model
# -------------------------------------------------------
# Next trait: n-alkane ACL
qplot(n.alkane_ACL,data=data,geom="histogram")

# QQ plot
qqnorm(data$n.alkane_ACL)
qqline(data$n.alkane_ACL)

# Data transformations
sq_nalkACL <- sqrt(data$n.alkane_ACL)
log_nalkACL <- log(data$n.alkane_ACL+1)
rec_nalkACL <- 1/(data$n.alkane_ACL)
asin_nalkACL <- asin(data$n.alkane_ACL/100)

# QQ plots
qqnorm(sq_nalkACL)
qqline(sq_nalkACL)

qqnorm(log_nalkACL)
qqline(log_nalkACL)

qqnorm(rec_nalkACL)
qqline(rec_nalkACL)

qqnorm(asin_nalkACL)
qqline(asin_nalkACL)

# Formally test deviation from normal distribution
shapiro.test(data$n.alkane_ACL)
shapiro.test(sq_nalkACL)
shapiro.test(log_nalkACL)
shapiro.test(rec_nalkACL)
shapiro.test(asin_nalkACL)

# Mixed effect linear model
nalkACL_model <- lmer(asin_nalkACL ~ il + (1|Plot), data=data)

# Function model check
modelcheck(nalkACL_model)

# Everything looks good! Now, find the p values for the ILs deviating from the cvM82 parent
summary(nalkACL_model)

# There appear to be many significant ILs. Write out data.
summary <- summary(nalkACL_model)$coef
write.table(summary, file="nalkACL_model_summary.txt")

# Heritability values
h2.trans <- lmer(asin_nalkACL ~ (1|il) + (1|Plot), data=data) #Transformed
h2.trans
h2.model <- lmer(data$n.alkane_ACL ~ (1|il) + (1|Plot), data=data) #Not transformed
h2.model
# -------------------------------------------------------

# Next trait: anteiso-alkane CPI
qplot(anteiso_CPI,data=data,geom="histogram")

# QQ plot
qqnorm(data$anteiso_CPI)
qqline(data$anteiso_CPI)

# Data transformations
sq_anteCPI <- sqrt(data$anteiso_CPI)
log_anteCPI <- log(data$anteiso_CPI+1)
rec_anteCPI <- 1/(data$anteiso_CPI)
asin_anteCPI <- asin(data$anteiso_CPI/100)

# QQ plots
qqnorm(sq_anteCPI)
qqline(sq_anteCPI)

qqnorm(log_anteCPI)
qqline(log_anteCPI)

qqnorm(rec_anteCPI)
qqline(rec_anteCPI)

qqnorm(asin_anteCPI)
qqline(asin_anteCPI)

# Formally test deviation from normal distribution
shapiro.test(data$anteiso_CPI)
shapiro.test(sq_anteCPI)
shapiro.test(log_anteCPI)
shapiro.test(rec_anteCPI)
shapiro.test(asin_anteCPI)

#Use no transformation
# Mixed effect linear model
anteCPI_model <- lmer(data$anteiso_CPI ~ il + (1|Plot), data=data)

# Function model check
modelcheck(anteCPI_model)

# Everything looks good! Now, find the p values for the ILs deviating from the cvM82 parent
summary(anteCPI_model)

# There appear to be many significant ILs. Write out data.
summary <- summary(anteCPI_model)$coef
write.table(summary, file="anteisoalkaneCPI_model_summary.txt")

# Heritability values
h2.model <- lmer(data$anteiso_CPI ~ (1|il) + (1|Plot), data=data) #Not transformed
h2.model
# -------------------------------------------------------
# Next trait: iso-alkane CPI
qplot(iso_CPI,data=data,geom="histogram")

# QQ plot
qqnorm(data$iso_CPI)
qqline(data$iso_CPI)

# Data transformations
sq_isoCPI <- sqrt(data$iso_CPI)
log_isoCPI <- log(data$iso_CPI+1)
rec_isoCPI <- 1/(data$iso_CPI)
asin_isoCPI <- asin(data$iso_CPI/100)

# QQ plots
qqnorm(sq_isoCPI)
qqline(sq_isoCPI)

qqnorm(log_isoCPI)
qqline(log_isoCPI)

qqnorm(rec_isoCPI)
qqline(rec_isoCPI)

qqnorm(asin_isoCPI)
qqline(asin_isoCPI)

# Formally test deviation from normal distribution
shapiro.test(data$iso_CPI)
shapiro.test(sq_isoCPI)
shapiro.test(log_isoCPI)
shapiro.test(rec_isoCPI)
shapiro.test(asin_isoCPI)

# Sqrt transformation
# Mixed effect linear model
isoCPI_model <- lmer(sq_isoCPI ~ il + (1|Plot), data=data)

# Function model check
modelcheck(isoCPI_model)

# Everything looks good! Now, find the p values for the ILs deviating from the cvM82 parent
summary(isoCPI_model)

# There appear to be many significant ILs. Write out data.
summary <- summary(isoCPI_model)$coef
write.table(summary, file="isoalkaneCPI_model_summary.txt")

# Heritability values
h2.trans <- lmer(sq_isoCPI ~ (1|il) + (1|Plot), data=data) #Transformed
h2.trans
h2.model <- lmer(data$iso_CPI ~ (1|il) + (1|Plot), data=data) #Not transformed
h2.model
# -----------------------------------------------------------
# Next trait: Methyl index
qplot(Methyl_Index,data=data,geom="histogram")

# QQ plot
qqnorm(data$Methyl_Index)
qqline(data$Methyl_Index)

# Data transformations
sq_mindex <- sqrt(data$Methyl_Index)
log_mindex <- log(data$Methyl_Index+1)
rec_mindex <- 1/(data$Methyl_Index)
asin_mindex <- asin(data$Methyl_Index/100)

# QQ plots
qqnorm(sq_mindex)
qqline(sq_mindex)

qqnorm(log_mindex)
qqline(log_mindex)

qqnorm(rec_mindex)
qqline(rec_mindex)

qqnorm(asin_mindex)
qqline(asin_mindex)

# Formall test deviation from normal distribution
shapiro.test(data$Methyl_Index)
shapiro.test(sq_mindex)
shapiro.test(log_mindex)
shapiro.test(rec_mindex)
shapiro.test(asin_mindex)

# sqrt transformation
# Mixed effect linear model
mindex_model <- lmer(sq_mindex ~ il + (1|Plot), data=data)

# Function model check
modelcheck(mindex_model)

# Looks good; find p values for ILs deviating from cvM82 parent
summary(mindex_model)

# Many significant ILs -- write out data
summary <- summary(mindex_model)$coef
write.table(summary, file="MethylIndex_model_summary.txt")

# Heritability values
h2.trans <- lmer(sq_mindex ~ (1|il) + (1|Plot), data=data) #Transformed
h2.trans
h2.model <- lmer(data$Methyl_Index ~ (1|il) + (1|Plot), data=data) #Not transformed
h2.model
# -----------------------------------------------------------
# Next trait: percent anteiso-alkanes
qplot(pct_anteiso.alkanes,data=data,geom="histogram")

# QQ plot
qqnorm(data$pct_anteiso.alkanes)
qqline(data$pct_anteiso.alkanes)

# Data transformations
sq_pctanteisoalk <- sqrt(data$pct_anteiso.alkanes)
log_pctanteisoalk <- log(data$pct_anteiso.alkanes+1)
rec_pctanteisoalk <- 1/(data$pct_anteiso.alkanes)
asin_pctanteisoalk <- asin(data$pct_anteiso.alkanes/100)

# QQ plots
qqnorm(sq_pctanteisoalk)
qqline(sq_pctanteisoalk)

qqnorm(log_pctanteisoalk)
qqline(log_pctanteisoalk)

qqnorm(rec_pctanteisoalk)
qqline(rec_pctanteisoalk)

qqnorm(asin_pctanteisoalk)
qqline(asin_pctanteisoalk)

# Formally test deviation from normal distribution
shapiro.test(data$pct_anteiso.alkanes)
shapiro.test(sq_pctanteisoalk)
shapiro.test(log_pctanteisoalk)
shapiro.test(rec_pctanteisoalk)
shapiro.test(asin_pctanteisoalk)

# Log transformation
# Mixed effect linear model
pctanteisoalk_model <- lmer(log_pctanteisoalk ~ il + (1|Plot), data=data)

# Function model check
modelcheck(pctanteisoalk_model)

# Looks good; find p values for ILs deviating from cvM82 parent
summary(pctanteisoalk_model)

# Many significant ILs -- write out data
summary <- summary(pctanteisoalk_model)$coef
write.table(summary, file="pctanteisoalkane_model_summary.txt")

# Heritability values
h2.trans <- lmer(log_pctanteisoalk ~ (1|il) + (1|Plot), data=data) #Transformed
h2.trans
h2.model <- lmer(data$pct_anteiso.alkanes ~ (1|il) + (1|Plot), data=data) #Not transformed
h2.model

#----------------
# Next trait: anteiso_ACL
# Take a look at it...
qplot(anteiso_ACL,data=data,geom="histogram")

# The dataset is skewed; look at a qq plot
qqnorm(data$anteiso_ACL)
qqline(data$anteiso_ACL)

# There is some deviation; let's try some transformations
sq_trait <- sqrt(data$anteiso_ACL)
log_trait <- log(data$anteiso_ACL+1)
rec_trait <- 1/(data$anteiso_ACL)
asin_trait <- asin(data$anteiso_ACL/100)

# QQ plots
qqnorm(sq_trait)
qqline(sq_trait)

qqnorm(log_trait)
qqline(log_trait)

qqnorm(rec_trait)
qqline(rec_trait)

qqnorm(asin_trait)
qqline(asin_trait)

# Formally test deviation from normal distribution
shapiro.test(data$anteiso_ACL)
shapiro.test(sq_trait)
shapiro.test(log_trait)
shapiro.test(rec_trait)
shapiro.test(asin_trait)

# No transformation
# Model our trait using a mixed effect linear model
anteiso_ACL_model <- lmer(data$anteiso_ACL ~ il + (1|Plot), data=data)

modelcheck(anteiso_ACL_model)

# Everything looks good! Now, find the p values for the ILs deviating from the cvM82 parent
summary(anteiso_ACL_model)

# There appear to be many significant ILs. Write out data.
summary <- summary(anteiso_ACL_model)$coef
write.table(summary, file="anteiso_ACL_model_summary.txt")

# Heritability values
h2.model <- lmer(data$anteiso_ACL ~ (1|il) + (1|Plot), data=data) #Not transformed
h2.model
# --------
# Next trait: iso_ACL
qplot(iso_ACL,data=data,geom="histogram")

# The dataset is skewed; look at a qq plot
qqnorm(data$iso_ACL)
qqline(data$iso_ACL)

# There is some deviation; let's try some transformations
sq_trait <- sqrt(data$iso_ACL)
log_trait <- log(data$iso_ACL+1)
rec_trait <- 1/(data$iso_ACL)
asin_trait <- asin(data$iso_ACL/100)

# QQ plots
qqnorm(sq_trait)
qqline(sq_trait)

qqnorm(log_trait)
qqline(log_trait)

qqnorm(rec_trait)
qqline(rec_trait)

qqnorm(asin_trait)
qqline(asin_trait)

# Formally test deviation from normal distribution
shapiro.test(data$iso_ACL)
shapiro.test(sq_trait)
shapiro.test(log_trait)
shapiro.test(rec_trait)
shapiro.test(asin_trait)

# Arcsine transformation
# Model our trait using a mixed effect linear model
iso_ACL_model <- lmer(asin_trait ~ il + (1|Plot), data=data)

modelcheck(iso_ACL_model)

# Everything looks good! Now, find the p values for the ILs deviating from the cvM82 parent
summary(iso_ACL_model)

# There appear to be many significant ILs. Write out data.
summary <- summary(iso_ACL_model)$coef
write.table(summary, file="iso_ACL_model_summary.txt")

# Heritability values
h2.trans <- lmer(asin_trait ~ (1|il) + (1|Plot), data=data) #Transformed
h2.trans
h2.model <- lmer(data$iso_ACL ~ (1|il) + (1|Plot), data=data) #Not transformed
h2.model


# --------

# Next trait: eps_iC31
qplot(eps13C_i.C31,data=data,geom="histogram")

# The dataset is skewed; look at a qq plot
qqnorm(data$eps13C_i.C31)
qqline(data$eps13C_i.C31)

# There is some deviation; let's try some transformations
sq_trait <- sqrt(abs(data$eps13C_i.C31))
log_trait <- log(abs(data$eps13C_i.C31)+1)
rec_trait <- 1/(abs(data$eps13C_i.C31))
asin_trait <- asin(abs(data$eps13C_i.C31)/100)

# QQ plots
qqnorm(sq_trait)
qqline(sq_trait)

qqnorm(log_trait)
qqline(log_trait)

qqnorm(rec_trait)
qqline(rec_trait)

qqnorm(asin_trait)
qqline(asin_trait)

# Formally test deviation from normal distribution
shapiro.test(data$eps13C_i.C31)
shapiro.test(sq_trait)
shapiro.test(log_trait)
shapiro.test(rec_trait)
shapiro.test(asin_trait)

# Arcsine transformation
# Model our trait using a mixed effect linear model
eps13C_i.C31_model <- lmer(asin_trait ~ il + (1|Plot), data=data)

modelcheck(eps13C_i.C31_model)

# Everything looks good! Now, find the p values for the ILs deviating from the cvM82 parent
summary(eps13C_i.C31_model)

# There appear to be many significant ILs. Write out data.
summary <- summary(eps13C_i.C31_model)$coef
write.table(summary, file="eps13C_i.C31_model_summary.txt")

# Heritability values
h2.trans <- lmer(asin_trait ~ (1|il) + (1|Plot), data=data) #Transformed
h2.trans
h2.model <- lmer(data$eps13C_i.C31 ~ (1|il) + (1|Plot), data=data) #Not transformed
h2.model
# --------

# Next trait: eps13C_n.C31
qplot(eps13C_n.C31,data=data,geom="histogram")

# The dataset is skewed; look at a qq plot
qqnorm(data$eps13C_n.C31)
qqline(data$eps13C_n.C31)

# There is some deviation; let's try some transformations
sq_trait <- sqrt(abs(data$eps13C_n.C31))
log_trait <- log(abs(data$eps13C_n.C31)+1)
rec_trait <- 1/(abs(data$eps13C_n.C31))
asin_trait <- asin(abs(data$eps13C_n.C31)/100)

# QQ plots
qqnorm(sq_trait)
qqline(sq_trait)

qqnorm(log_trait)
qqline(log_trait)

qqnorm(rec_trait)
qqline(rec_trait)

qqnorm(asin_trait)
qqline(asin_trait)

# Formally test deviation from normal distribution
shapiro.test(data$eps13C_n.C31)
shapiro.test(sq_trait)
shapiro.test(log_trait)
shapiro.test(rec_trait)
shapiro.test(asin_trait)

# Arcsine transformation
# Model our trait using a mixed effect linear model
eps13C_n.C31_model <- lmer(asin_trait ~ il + (1|Plot), data=data)

modelcheck(eps13C_n.C31_model)

# Everything looks good! Now, find the p values for the ILs deviating from the cvM82 parent
summary(eps13C_n.C31_model)

# There appear to be many significant ILs. Write out data.
summary <- summary(eps13C_n.C31_model)$coef
write.table(summary, file="eps13C_n.C31_model_summary.txt")

# Heritability values
h2.trans <- lmer(asin_trait ~ (1|il) + (1|Plot), data=data) #Transformed
h2.trans
h2.model <- lmer(data$eps13C_n.C31 ~ (1|il) + (1|Plot), data=data) #Not transformed
h2.model
# -------

# Next trait: eps13C_i.C33
qplot(D_iC33,data=data,geom="histogram")

# The dataset is skewed; look at a qq plot
qqnorm(data$eps13C_i.C33)
qqline(data$eps13C_i.C33)

# There is some deviation; let's try some transformations
sq_trait <- sqrt(abs(data$eps13C_i.C33))
log_trait <- log(abs(data$eps13C_i.C33)+1)
rec_trait <- 1/(abs(data$eps13C_i.C33))
asin_trait <- asin(abs(data$eps13C_i.C33)/100)

# QQ plots
qqnorm(sq_trait)
qqline(sq_trait)

qqnorm(log_trait)
qqline(log_trait)

qqnorm(rec_trait)
qqline(rec_trait)

qqnorm(asin_trait)
qqline(asin_trait)

# Formally test deviation from normal distribution
shapiro.test(data$eps13C_i.C33)
shapiro.test(sq_trait)
shapiro.test(log_trait)
shapiro.test(rec_trait)
shapiro.test(asin_trait)

# sqrt transformation
# Model our trait using a mixed effect linear model
eps13C_i.C33_model <- lmer(sq_trait ~ il + (1|Plot), data=data)

modelcheck(eps13C_i.C33_model)

# Everything looks good! Now, find the p values for the ILs deviating from the cvM82 parent
summary(eps13C_i.C33_model)

# There appear to be many significant ILs. Write out data.
summary <- summary(eps13C_i.C33_model)$coef
write.table(summary, file="D_iC33_model_summary.txt")

h2.trans <- lmer(sq_trait ~ (1|il) + (1|Plot), data=data) #Transformed
h2.trans
h2.model <- lmer(data$eps13C_i.C33 ~ (1|il) + (1|Plot), data=data) #Not transformed
h2.model
# -------

# Next trait: eps13C_n.C33
qplot(D_nC33,data=data,geom="histogram")

# The dataset is skewed; look at a qq plot
qqnorm(data$eps13C_n.C33)
qqline(data$eps13C_n.C33)

# There is some deviation; let's try some transformations
sq_trait <- sqrt(abs(data$eps13C_n.C33))
log_trait <- log((data$eps13C_n.C33)+1)
rec_trait <- 1/(abs(data$eps13C_n.C33))
asin_trait <- asin((data$eps13C_n.C33)/100)

# QQ plots
qqnorm(sq_trait)
qqline(sq_trait)

qqnorm(log_trait)
qqline(log_trait)

qqnorm(rec_trait)
qqline(rec_trait)

qqnorm(asin_trait)
qqline(asin_trait)

# Formally test deviation from normal distribution
shapiro.test(data$eps13C_n.C33)
shapiro.test(sq_trait)
shapiro.test(log_trait)
shapiro.test(rec_trait)
shapiro.test(asin_trait)

# No transformation
# Model our trait using a mixed effect linear model
D_nC33_model <- lmer(data$eps13C_n.C33 ~ il + (1|Plot), data=data)

modelcheck(eps13C_n.C33_model)

# Everything looks good! Now, find the p values for the ILs deviating from the cvM82 parent
summary(eps13C_n.C33_model)

# There appear to be many significant ILs. Write out data.
summary <- summary(eps13C_n.C33_model)$coef
write.table(summary, file="eps13C_n.C33_model_summary.txt")

# Heritability values
h2.model <- lmer(data$eps13C_n.C33 ~ (1|il) + (1|Plot), data=data) #Not transformed
h2.model
# -------

# Next trait: n31.i31
qplot(n31.i31,data=data,geom="histogram")

# The dataset is skewed; look at a qq plot
qqnorm(data$n31.i31)
qqline(data$n31.i31)

# There is some deviation; let's try some transformations
sq_trait <- sqrt(abs(data$n31.i31))
log_trait <- log(abs(data$n31.i31)+1)
rec_trait <- 1/(abs(data$n31.i31))
asin_trait <- asin(abs(data$n31.i31)/100)

# QQ plots
qqnorm(sq_trait)
qqline(sq_trait)

qqnorm(log_trait)
qqline(log_trait)

qqnorm(rec_trait)
qqline(rec_trait)

qqnorm(asin_trait)
qqline(asin_trait)

# Formally test deviation from normal distribution
shapiro.test(data$n31.i31)
shapiro.test(sq_trait)
shapiro.test(log_trait)
shapiro.test(rec_trait)
shapiro.test(asin_trait)

# No transformation
# Model our trait using a mixed effect linear model
n31.i31_model <- lmer(data$n31.i31 ~ il + (1|Plot), data=data)

modelcheck(n31.i31_model)

# Everything looks good! Now, find the p values for the ILs deviating from the cvM82 parent
summary(n31.i31_model)

# There appear to be many significant ILs. Write out data.
summary <- summary(n31.i31_model)$coef
write.table(summary, file="n31.i31_model_summary.txt")

# Heritability values
h2.model <- lmer(data$n31.i31 ~ (1|il) + (1|Plot), data=data) #Not transformed
h2.model
# -------

# Next trait: n33.i33
qplot(n33.i33,data=data,geom="histogram")

# The dataset is skewed; look at a qq plot
qqnorm(data$n33.i33)
qqline(data$n33.i33)

# There is some deviation; let's try some transformations
sq_trait <- sqrt(abs(data$n33.i33))
log_trait <- log(abs(data$n33.i33)+1)
rec_trait <- 1/(abs(data$n33.i33))
asin_trait <- asin(abs(data$n33.i33)/100)

# QQ plots
qqnorm(sq_trait)
qqline(sq_trait)

qqnorm(log_trait)
qqline(log_trait)

qqnorm(rec_trait)
qqline(rec_trait)

qqnorm(asin_trait)
qqline(asin_trait)

# Formally test deviation from normal distribution
shapiro.test(data$n33.i33)
shapiro.test(sq_trait)
shapiro.test(log_trait)
shapiro.test(rec_trait)
shapiro.test(asin_trait)

# reciprocal transformation
# Model our trait using a mixed effect linear model
n33.i33_model <- lmer(sq_trait ~ il + (1|Plot), data=data)

modelcheck(n33.i33_model)

# Everything looks good! Now, find the p values for the ILs deviating from the cvM82 parent
summary(n33.i33_model)

# There appear to be many significant ILs. Write out data.
summary <- summary(n33.i33_model)$coef
write.table(summary, file="n33.i33_model_summary.txt")

# Heritability values
h2.trans <- lmer(sq_trait ~ (1|il) + (1|Plot), data=data) #Transformed
h2.trans
h2.model <- lmer(data$n33.i33 ~ (1|il) + (1|Plot), data=data) #Not transformed
h2.model
# -------

# Next trait: Sum_Alk.area
qplot(Sum_Alk.area,data=data,geom="histogram")

# The dataset is skewed; look at a qq plot
qqnorm(data$Sum_Alk.area)
qqline(data$Sum_Alk.area)

# There is some deviation; let's try some transformations
sq_trait <- sqrt(data$Sum_Alk.area)
log_trait <- log(data$Sum_Alk.area+1)
rec_trait <- 1/(data$Sum_Alk.area)
asin_trait <- asin(data$Sum_Alk.area/100)

# QQ plots
qqnorm(sq_trait)
qqline(sq_trait)

qqnorm(log_trait)
qqline(log_trait)

qqnorm(rec_trait)
qqline(rec_trait)

qqnorm(asin_trait)
qqline(asin_trait)

# Formally test deviation from normal distribution
shapiro.test(data$Sum_Alk.area)
shapiro.test(sq_trait)
shapiro.test(log_trait)
shapiro.test(rec_trait)
shapiro.test(asin_trait)

# sqrt transformation
# Model our trait using a mixed effect linear model
Sum_Alk.area_model <- lmer(sq_trait ~ il + (1|Plot), data=data)

modelcheck(Sum_Alk.area_model)

# Everything looks good! Now, find the p values for the ILs deviating from the cvM82 parent
summary(Sum_Alk.area_model)

# There appear to be many significant ILs. Write out data.
summary <- summary(Sum_Alk.area_model)$coef
write.table(summary, file="Sum_Alk.area_model_summary.txt")

# Heritability values
h2.trans <- lmer(sq_trait ~ (1|il) + (1|Plot), data=data) #Transformed
h2.trans
h2.model <- lmer(data$Sum_Alk.area ~ (1|il) + (1|Plot), data=data) #Not transformed
h2.model
