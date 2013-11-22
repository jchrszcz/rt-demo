# uncomment and install these if first time using
# install.packages("Rcpp")
# install.packages("devtools")
# install.packages("assist")
# install.packages("signals")
# library(devtools); install_github("gemmR", "jchrszcz", subdir = "gemmR")

# load these every time, but have to install first
library(Rcpp)
library(signals)
library(assist)
library(gemmR)

# load data into workspace, have to download from github

load("demo.RData")
source("vzfilter.R")

# run filtering

subjects <- unique(df$Subj)
df$new.vz.filter <- df$log.rt

for (i in 1:length(subjects)) {
  temp <- df$new.vz.filter[df$Subj == subjects[i] & !is.na(df$new.vz.filter)]
  temp <- est.norm.smoothing.splines(temp)
  df$new.vz.filter[df$Subj == subjects[i] & !is.na(df$new.vz.filter)] <- temp
}

# repeated-measures ANOVA on raw RTs

new.df <- with(df, aggregate(RT ~ Center + Difference + Gain + Subj, FUN = median))
raw.mod <- lm(RT ~ Center * Difference * Gain + Subj, data = new.df)
anova(raw.mod)

# repeated-measures ANOVA on filtered RTs

vz.df <- with(df, aggregate(new.vz.filter ~ Center + Difference + Gain + Subj, FUN = median))
vz.mod <- lm(new.vz.filter ~ Center * Difference * Gain + Subj, data = vz.df)
anova(vz.mod)

# GeMM on RTs

# scaling the RTs just makes them more amenable to search, maintains orders
new.df$RT <- as.vector(scale(new.df$RT))

# run gemm
a <- gemm(RT ~ Center * Difference * Gain + Subj, data = new.df)
print(a)