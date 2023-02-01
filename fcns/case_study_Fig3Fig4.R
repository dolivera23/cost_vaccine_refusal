
## ----------------------------------------------------------------------------------------- ##
##                 File to estimate the costs for a European population.
##                 Specific mixing coefficient and antivaccination population are needed. 
## ---------------------------------------------------------------------------------------- ##

# Needed functions
source ("fcns/fcns.R")               # Equations and solve functions 
source ("fcns/screening_fcns.R")     # output analysis function 
source ("fcns/CFR.R")                # CFR estimation per age group
source("fcns/cost_fcn.R")            # Functions to estimate costs based on simulation outputs
source("fcns/fcn_plots.R")           # Plotting functions 

# Needed packages 
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
value_qaly <-30000

#############################################################################################
#                         Gettting parameters up to date 
############################################################################################
#

# CPI conversion factors
getSymbols("GBRCPIALLMINMEI", src="FRED")
avg.cpi <- apply.yearly(GBRCPIALLMINMEI,mean)

ratio2012 <- as.numeric(avg.cpi["2022"])/ as.numeric(avg.cpi["2012"])

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


#############################################################################################
#                     Extracting incidence and costs 
############################################################################################
#
av_toextract <- 0.03
coef_toextract <- 0.292 #0.0543

av_pop <- N_t*av_toextract
pv_pop <- N_t * (1-av_toextract)

## -------------------------------------------------------- ##
##                  Extracting data frame 
## ------------------------------------------------------ ##

pos_av <- which.min(abs(antivax_perc-av_toextract))
pos_coef <- which.min(abs(coef-coef_toextract ))

df<- data.frame(data_out[[pos_coef]][[pos_av]])

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

# 
# inflation <- read_excel("./fcns/inflation.xls")
# 
# for (i in 1:nrow(costs)){
#   
#   rate <- inflation  %>% filter(Year ==  costs[i,]$year)     
#   
#   costs[i,5:13] <- costs[i,5:13] *rate$Total         #rate$Total = 1%  rate$Total_upper = 3%
# }


#Cumulative values
costs <- costs %>% group_by(population) %>% mutate(cum_inc = cumsum(cases),
                                                   cum_death = cumsum(deaths),
                                                   cum_death_low = cumsum(deaths_low),
                                                   cum_death_high = cumsum(deaths_high),
                                                   cum_costs = cumsum(cost),
                                                   cum_costs_low = cumsum(cost_low),
                                                   cum_costs_high = cumsum(cost_high))

#############################################################################################
#                                        Plots  
############################################################################################
#

# Figure 3:  Incidence and deaths plots 


plot_inc <- plot_incidence (df=costs, N_t, yscale =60e6)
plot_d <- plot_deaths (df=costs,N_t, yscale=60e6)

plot_deaths_total <- plot_cumu_deaths(df=costs, N_t, yscale=60e6,last_year=20) + theme(legend.position = "bottom")
plot_inc_total <-   plot_cumu_inc (df=costs, N_t, yscale=60e6,last_year=20)

legend_cum <- get_legend(plot_deaths_total)


top <- plot_grid(plot_inc,
                 plot_inc_total + theme(legend.position = "none"),
                 plot_d,
                 plot_deaths_total+theme(legend.position = "none"),
                 labels = c("A","B","C","D"),
                 ncol=2)

plot_grid (top, legend_cum, ncol=1, rel_heights = c(1,0.1))


# Figure 4: Costs 

plot_costs30 <- plot_bar_costs (costs, N_t, 60e6)
plot_cum30 <- plot_cum_costs(costs,N_t,yscale = 60e6)

lab1 <- get_legend(plot_costs30)
lab2 <- get_legend (plot_cum30)

right <- plot_grid(plot_costs30 + theme(legend.position = "none"),
                   #plot_costs13 + theme (legend.position = "none"),# + scale_y_continuous(labels = comma,limits = c(0,2e8)),
                   plot_cum30+ theme (legend.position = "none"),
                  #plot_cum13 + theme (legend.position = "none"),# + scale_y_continuous(labels = comma,limits = c(0,1000000000)),
                   labels = c("A","B"),
                   ncol=1)
left <- plot_grid(lab1,lab2,ncol=1)


plot_grid(right,left, nrow=1, rel_widths = c(1,0.4))


### Supplementary: Per capita plots ####

costs <- costs  %>% mutate (pc_death = 0,
                            pc_death_low =0,
                            pc_death_high = 0,
                            pc_cost = 0,
                            pc_cost_low = 0,
                            pc_cost_high = 0) 

pro <- costs %>% filter(population == "Pro vaccination") %>%
  mutate(pc_inc = cases/pv_pop,
         pc_death = deaths/pv_pop,
         pc_death_low = deaths_low/pv_pop,
         pc_death_high = deaths_high/pv_pop,
         pc_cost = cost/pv_pop,
         pc_cost_low = cost_low/pv_pop,
         pc_cost_high = cost_high/pv_pop)

anti <- costs %>% filter (population == "Anti vaccination")%>%
  mutate(pc_inc = cases/av_pop,
         pc_death = deaths/av_pop,
         pc_death_low = deaths_low/av_pop,
         pc_death_high = deaths_high/av_pop,
         pc_cost = cost/av_pop,
         pc_cost_low = cost_low/av_pop,
         pc_cost_high = cost_high/av_pop)

pc_costs <- rbind (pro,anti)

# Incidence
plot_inc_pc <- ggplot(pc_costs, aes(x=year, y=pc_inc*1e6)) +
  geom_line(size=1.2,color="#00758F") +
   xlab("Time (years)") + ylab("Yearly Incidence per \n million individuals") +
  theme_bw(base_size = 13) +
  theme(strip.background = element_rect(fill = NA),
        panel.border = element_blank(),
        axis.line = element_line(),
        legend.text.align = 0)+
  facet_wrap(~population,scales = "free_y")

#Deaths 
plot_deaths_pc <- ggplot(pc_costs, aes(x=year, y=pc_death*1e6)) +
  geom_ribbon(aes(ymin = pc_death_low*1e6, ymax = pc_death_high*1e6), fill = "#efeddc") +
  #geom_errorbar(aes(ymin= (deaths_low/N_t)*yscale, ymax=(deaths_high/N_t)*yscale), width=.7,color= "#efeddc")+  
  geom_point(size=3, shape=21, fill="white",color="#efeddc")+
  geom_line(color ="#8c8340",size=1.2) +
xlab("Time (years)") + ylab("Yearly Deaths per \n million individuals") +
  theme_bw(base_size = 13) +
  theme(strip.background = element_rect(fill = NA),
        panel.border = element_blank(),
        axis.line = element_line(),
        legend.text.align = 0)+
  facet_wrap(~population,scales = "free_y")

# Costs
plot_costs_pc <- ggplot(pc_costs, aes(x=year, y=pc_cost*1e6)) +
  geom_ribbon(aes(ymin = pc_cost_low*1e6, ymax = pc_cost_high*1e6), fill = "#faedea") +
  #geom_errorbar(aes(ymin= (deaths_low/N_t)*yscale, ymax=(deaths_high/N_t)*yscale), width=.7,color= "#efeddc")+  
  geom_point(size=3, shape=21, fill="white",color="#f1c9c1")+
  geom_line(color ="#bb452a",size=1.2) +
  xlab("Time (years)") + ylab("Yearly Costs (GBP) per \n million individuals") +
  theme_bw(base_size = 13) +
  theme(panel.border = element_blank(),
        axis.line = element_line(),
        legend.text.align = 0)+
  facet_wrap(~population,scales = "free_y") +
  scale_y_continuous(labels = comma)

plot_grid (plot_inc_pc,
           plot_deaths_pc,
           plot_costs_pc,
           ncol =1,
           labels = c("A","B", "C"),
           align = 'v', axis = 'tb')


## =========================================================================== ##
##                       To document                                           ##
## ==========================================================================  ##

# Total Incidence 

last_cumu <- costs %>% filter (year ==20) 
total_incs <- last_cumu %>% filter (population =="Total") %>% select(cum_inc)
total_deah <- last_cumu %>% filter (population =="Total") %>% select(cum_death)

# Max costs

tots <- costs %>% filter(population =="Total") %>% slice(which.max(cost))
(tots$cost/N_t)*60e6

tots <- costs %>% filter(population =="Total")

x<- tots$qaly_death/tots$cost
x



NPV <-  tots$cum_costs_low[20]/(1+0.03)^20
(NPV/N_t)*60
