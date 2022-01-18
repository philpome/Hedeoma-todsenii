hed.all <- read.csv("hedeoma-gccs-all-samples.csv",header=TRUE)
library(gtsummary)
all.sum <- hed.all %>% select(Status, CoD, DeathAge)
tbl1 <- all.sum %>%
     tbl_summary(
         statistic = list(all_continuous() ~ "{mean} ({sd})",
         all_categorical() ~ "{n} / {N} ({p}%)"),
         digits = all_continuous() ~ 2
     )
library(flextable)
library(webshot)
tbl1 %>%
     as_gt() %>%
     gt::gtsave(filename = "hed-all-gccs-sum.png")
hed.all %>% 
  group_by(CoD) %>%
  summarize(mean=mean(DeathAge, na.rm=TRUE))
