# RT demo file

# 1. vanZandt's smoothing spline filter

load("demo.RData")
source("vzfilter.R")

subjects <- unique(df$Subj)
df$new.vz.filter <- df$log.rt

for (i in 1:length(subjects)) {
  temp <- df$new.vz.filter[df$Subj == subjects[i] & !is.na(df$new.vz.filter)]
  temp <- est.norm.smoothing.splines(temp)
  df$new.vz.filter[df$Subj == subjects[i] & !is.na(df$new.vz.filter)] <- temp
}

