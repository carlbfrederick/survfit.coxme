#Some Examples to test with

library(haven)
library(survival)
library(coxme)
library(dplyr)

#. Allison: Cox Regression, Table 6.2, columns 8:10----
arrests <- read_dta("http://statisticalhorizons.com/wp-content/uploads/arrests.dta")

(m1 <- coxme(Surv(length, arrind)~
               fin+age+white+male+married+paro+numprop+crimprop+numarst+edcomb+(1|id), 
             data=arrests, y = TRUE, x=TRUE))

