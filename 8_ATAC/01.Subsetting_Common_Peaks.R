rm(list = ls())

INPUT_DIR="/home/cascione/cascione@Bertoni/Cascione/ATACSeq_Verrilio/results/"
setwd(dir = INPUT_DIR)

#ORIGINAL SERVER LOCATION df <- read.table("./bwa/merged_library/macs2/broad_peak/consensus/consensus_peaks.mLb.clN.boolean.txt", header = TRUE)
#READ THE TXT FILE df <- read.table("~/ATACSEQ_Verrillo/consensus_peaks.mLb.clN.boolean.txt", header = TRUE)
#SAVE AS RDS FILE saveRDS(object = df, file = "~/ATACSEQ_Verrillo/consensus_peaks.mLb.clN.boolean.rds")

df <- readRDS(file = "~/ATACSEQ_Verrillo/consensus_peaks.mLb.clN.boolean.rds")
dim(df)
df$PREEXPOSURE <- ifelse(df$MMCB_PREEXPOSURE_REP1.mLb.clN.bool | df$MMCB_PREEXPOSURE_REP2.mLb.clN.bool | df$MMCB_PREEXPOSURE_REP3.mLb.clN.bool, yes = TRUE, no = FALSE)
df$TWOWEEKS <- ifelse(df$MMCB_2WEEKS_REP1.mLb.clN.bool | df$MMCB_2WEEKS_REP2.mLb.clN.bool | df$MMCB_2WEEKS_REP3.mLb.clN.bool, yes = TRUE, no = FALSE)
df$ONEMONTH <- ifelse(df$MELANOMA_1MONTH_REP1.mLb.clN.bool | df$MELANOMA_1MONTH_REP2.mLb.clN.bool | df$MELANOMA_1MONTH_REP3.mLb.clN.bool, yes = TRUE, no = FALSE)

#Regions detected in PREEXPOSURE but not present in both TWOWEEKS and ONEMONTH are excluded
df$PREEXPOSURE[df$PREEXPOSURE & !(df$TWOWEEKS & df$ONEMONTH)] <- FALSE
table(df$PREEXPOSURE | df$TWOWEEKS | df$ONEMONTH)

dim(df)
df <- df[(df$PREEXPOSURE | df$TWOWEEKS | df$ONEMONTH),]; dim(df)
saveRDS(object = df, file = "~/ATACSEQ_Verrillo/consensus_peaks.forDifferentialAccessibilityAnalysis.rds")
