library('readxl')
library ('ggplot2')
library('reshape2')
library('Hmisc')
library('ggthemes')

cases <- read_excel("cases_age.xlsx",sheet = "Cases")         # Reported cases 
deaths <- read_excel("cases_age.xlsx",sheet = "Deaths")       # Reported deaths by age group 

total_cases <- rowSums(cases[, 2:19])
total_deaths <- rowSums(deaths[, 2:19])

age_distr <- binconf(rowSums(cases[, 2:20]), sum(rowSums(cases[, 2:20])), alpha=0.05,return.df=TRUE) 
age_distr$Age_group <- cases$Age_group


cases_long <- melt(cases, id.vars = "Age_group")
cases_long$Age_group <- as.factor(cases_long$Age_group)#



cases_long$Age_group <- factor(cases_long$Age_group , levels = c("< 1", "1-4", "5-9", "10-14",  "15-19", "20-24", "25-29","30-34", "> 35"))


plot_age <- ggplot(cases_long, aes(x=Age_group, y= value/sum(value), fill=variable)) +   geom_bar(stat="identity")  +  
  geom_errorbar(data=age_distr, aes(x= Age_group,ymin= Lower, ymax=Upper,group=1), width=.1 , inherit.aes = FALSE)+
  geom_point(data=age_distr, aes(x=Age_group, y = PointEst,group=1),inherit.aes =  FALSE) +
  scale_fill_hue(l=40, c=35, name= "Year") + xlab("Age group") + ylab("Frequency ") +theme_bw(base_family = 15) +
  theme(strip.background = element_rect(fill = NA),
        panel.border = element_blank(),
        axis.line = element_line(),
        legend.text.align = 0)



CFR <- binconf(total_deaths, total_cases, alpha=0.05,return.df=TRUE) 
CFR$Age_group <- cases$Age_group

CFR$Age_group <- factor(CFR$Age_group , levels = c("< 1", "1-4", "5-9", "10-14",  "15-19", "20-24", "25-29","30-34", "> 35"))

plot_CFR<- ggplot(CFR, aes(x=Age_group, y=PointEst)) + geom_errorbar(aes(ymin= Lower, ymax=Upper), width=.5) +  
  geom_point(size=2, shape=21, fill="white") + xlab("Age group") + ylab("CFR") + theme_bw(base_family = 15) + 
  theme(strip.background = element_rect(fill = NA),
        panel.border = element_blank(),
        axis.line = element_line(),
        legend.text.align = 0)


#plot_grid(plot_age,plot_CFR,labels = c("A","B"), nrow=1)
