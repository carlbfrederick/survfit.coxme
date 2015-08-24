#Some Examples to test with

library(haven)
library(survival)
library(coxme)
library(dplyr)

#From Allison: Cox Regression, (Modified) Table 6.2, columns 8:10  ----
arrests <- read_dta("http://statisticalhorizons.com/wp-content/uploads/arrests.dta")

(m0 <- coxph(Surv(length, arrind)~
               fin+age+male+married+paro+numprop+crimprop+numarst+edcomb+strata(race), 
             data=arrests, y = TRUE, x=TRUE))

(m1 <- coxme(Surv(length, arrind)~
               fin+age+male+married+paro+numprop+crimprop+numarst+edcomb+strata(race)+(1|id), 
             data=arrests, y = TRUE, x=TRUE))

survfit(m0)
