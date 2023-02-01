
source ("fcns/fcns.R")
source ("fcns/screening_fcns.R")
source ("fcns/CFR.R")
source("fcns/cost_fcn.R")


library ("rootSolve")
library("ggpubr")
library("gridExtra")
library("wesanderson")
library("gtable")
library("quantmod")
library("reshape2")
library("FinCal")
library(cowplot) 
library("readxl")
library("dplyr")
library(forcats)
library("scales")
## Model parameter definition

BR <- 11.1/52000    #

npop <- 2
B <- c(BR,BR)
miu <- c(BR,BR)
rec <- c(0.5, 0.5)
ro <- c(15,15)
ef <- 0.95
cov <- c(0.9,0.9)
N_t <- 100000

# Economic parameters

p_healthcare <- 275.96       # 2012 GBP
p_productiviy <- 794.43      # 2012 GBP

path_model <- "fcns/measles_model.R"

QALY_death <- data.frame(matrix(0, ncol=2, nrow=nrow(CFR)))
QALY_death[,1] <- CFR$Age_group
QALY_death[,2]  <- c(26.98730376,	26.83265656,	26.4035706,	25.90614236,	25.32948671,	24.66098475,	23.88600777,	22.98759705,	15.42456622)

case_qaly <- 0.019# data.frame(Value = 0.019, Low = 0.016, High = 0.022)


#############################################################################################
#                         Gettting parameters up to date 
############################################################################################
#

# CPI conversion factors
getSymbols("GBRCPIALLMINMEI", src="FRED")
avg.cpi <- apply.yearly(GBRCPIALLMINMEI,mean)

ratio2012 <- as.numeric(avg.cpi["2022"])/ as.numeric(avg.cpi["2012"])

# Conversion rates: GBP/USD 
#getFX(c("GBP/USD"), from = Sys.Date() - 30, to = Sys.Date(),env = parent.frame(), verbose = FALSE, warning = TRUE, auto.assign = TRUE)

#exchange_rate <- as.numeric(GBPUSD[Sys.Date()-1])


p_healthcare <- p_healthcare * ratio2012 
p_productiviy <- p_productiviy * ratio2012


##################################################################################################
#                   Screening different pop size and mixing patterns                            #
#################################################################################################

###  Running conditions parameters ####

cov <-c(0.9,0)

#Time 
years <- 250
tt <- seq(0, years*52, length.out =years*52)

antivax_perc <- seq(from=1e-5, to=1.000001, by=0.01)
coef  <- seq(from=0, to=0.5, by=0.005)


#data_out <- screening(years = years,tt = tt,antivax_perc = antivax_perc,coef = coef,path_model = path_model,npop = npop,B = B,miu = miu,cov = cov,rec = rec,ro = ro,ef = ef,N_t = N_t,I0_t = I0_t,R0_t = R0_t)
#saveRDS(data_out, "screening0-02.rds")

data_out <- readRDS("screening0-05.rds")

inflation <- read_excel("./inflation.xls")

#############################################################################################
#                    Function for costs data frame 
############################################################################################
#

coef_toextract <- 0.292 #0.0543
pos_coef <- which.min(abs(coef-coef_toextract ))

# Function that gets cumulative values at year 20 
# @param av: Anti-vaccination percentage 
# @param value_qaly: GBP value of QUALY (13 000 or 30 000)
# @param inf: Inflation rate: 1% or 3% 


cum_total_costs <- function (av, value_qaly, inf){#,CFR,age_distr,QALY_death,case_qualy,p_productiviy,p_healthcare,inflation){
  
  av_toextract <- av
 
  ## -------------------------------------------------------- ##
  ##                  Extracting data frame 
  ## ------------------------------------------------------ ##
  
  pos_av <- which.min(abs(antivax_perc-av_toextract))
 
  
  df <- data.frame(data_out[[pos_coef]][[pos_av]])
  
  ## Baseline scenatio Antivax == 0 ##
  baseline <- data.frame(data_out[[1]][[pos_av]])
  baseline$inc.2 <- df$inc.2[1]
  
  
  costs_total <- costs_output(df = df,CFR = CFR,QALY_death = QALY_death,
                              years_total = 20,age_distr = age_distr,case_qualy =case_qualy,value_qaly = value_qaly,
                              p_productiviy = p_productiviy,p_healthcare = p_healthcare )
  
  costs_baseline <- costs_output(df=baseline,CFR = CFR,QALY_death = QALY_death,
                                 years_total = 20,age_distr = age_distr,case_qualy =case_qualy,value_qaly = value_qaly,
                                 p_productiviy = p_productiviy,p_healthcare = p_healthcare)
  
  costs <- costs_total
  costs[,1:13] <- costs_total[,1:13] - costs_baseline[,1:13]
  
  # Looop to include inflation in the costs if required: 
  
  # for (i in 1:nrow(costs)){
  #   
  #   rate <- inflation  %>% filter(Year ==  costs[i,]$year)
  #   
  #   if(inf ==1){
  #     
  #     costs[i,5:13] <- costs[i,5:13] *rate$Total
  #   }
  #   
  #   else if (inf == 3 ){
  #     costs[i,5:13] <- costs[i,5:13] *rate$Total_upper
  #   }
  # }
  # 
  #Cumulative values
  costs <- costs %>% group_by(population) %>% mutate(across(c(1:13), .fns = list(cum = ~cumsum(.) )))
  
  last <-costs %>% filter(year ==20)
  last$av <- av_toextract
  
  return (last)
  
}


## Estimating the cumulative costs for all anti-vaccination population sizes:     #### 

qaly <- 30000
inf <- 1

big_df <- data.frame(matrix (0, nrow=0,ncol=29))

for( i in 1:length(antivax_perc)){
  
  temp <-  cum_total_costs(av=antivax_perc[i], value_qaly = qaly,inf = inf)#,CFR = CFR,
  #age_distr = age_distr,QALY_death = QALY_death,case_qualy = case_qualy,
  # p_productiviy = p_productiviy,p_healthcare =  p_healthcare)
  
  big_df <- rbind (big_df, temp)
  
  
}

per_capita  <- big_df[,14:29]  %>% filter(population =="Total") %>%
  mutate(across(c(3:15), .fns = list(pc = ~.x/N_t ))) %>%
  mutate(across(c(3:15), .fns = list(pa = ~.x/(av*N_t) )))

npv <-  per_capita 
npv[,c(7:15,20:28,34:41)] <- npv[,c(7:15,20:28,34:41)]/(1+0.03)^20

x <- npv %>% filter (av == 0.25001)

# To text
x$cost_cum_pa
x$cost_high_cum_pa
x$cost_low_cum_pa

x$cost_cum_pc
x$cost_high_cum_pc
x$cost_low_cum_pc


x$cost_cum /N_t * 60e6
x$cost_low_cum /N_t * 60e6
x$cost_high_cum /N_t * 60e6

write.csv(npv, "./npv_inlf1_qualy30.csv")

npv<- read.csv(file="./npv_inlf1_qualy30.csv")

ggplot(data =  npv, aes(x=av*100,y=cost_cum_pa))+ 
  geom_ribbon(aes(ymin= cost_low_cum_pa, ymax=cost_high_cum_pa), fill ="#f4f5ef") +
  geom_line(size=1.2, color = "#4b5335") +
   xlab("Anti-vaccination population percentage") + ylab("Welfare loss per-antivaccinator (GBP)") +
  xlim(0,100) + theme_bw(base_size=15)+
  theme(strip.background = element_rect(fill = NA),
        panel.border = element_blank(),
        axis.line = element_line(),
        legend.text.align = 0)
 


