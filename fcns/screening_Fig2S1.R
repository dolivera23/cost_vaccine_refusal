
source ("fcns/fcns.R")
source ("fcns/screening_fcns.R")
library ("rootSolve")
library("ggpubr")
library("gridExtra")
library("wesanderson")
library("gtable")
library("ggthemes")
library("cowplot")
library("dplyr")


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


path_model <- "fcns/measles_model.R"

##################################################################################################
#                   Finding steady state for initial conditions                                  #
#################################################################################################

years_ss <- 200000

ini <- initial_cond(years=years_ss, path_model = path_model,npop=npop,B=B,mi=miu,cov=cov,rec=rec,ro=ro,ef=ef, N_t=N_t)

I0_t <- sum(as.numeric(ini[1,4:5]))
R0_t <- sum(as.numeric(ini[1,6:7]))

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
#saveRDS(data_out, "screening0-05.rds")

data_out <- readRDS("screening0-05.rds")

##################################################################################################
#                  Extracting first outbreak ending time for all                                #
#################################################################################################

# First outbreak position for each coef and antivax population
pos_fo <- data.frame(matrix(0, nrow = length(coef),ncol=length(antivax_perc)))   #  [coef,antivx]

for (i in 1:length(coef)){
  
  temp_av <- vector(mode = "list", length = length(antivax_perc)) 
  
  for (j in 1: length(antivax_perc)){
    
    pos<- first_outbreak_time(data_out[[i]][[j]])
    
    pos_fo[i,j] <- pos +1 
    
  }
  
}

##################################################################################################
#                           Figure 2: Dilution effect of vaccination:                           #
#################################################################################################

## Base line scenario for comparaison: 
x<- data_out[[1]][[1]]
baseline <- data.frame(x)


# Extracting data for especific Antivaxxer perecentages

## ax_$$ List with two objects: One data frame for pro vaccination pop[1] and one data frame for anti vaccination pop[2]
## Dara frame have 3 columns: Time, extra cases during first outbreak, populations

ax_1 <- perc_extraction(toextract = 0.05,antivax_perc = antivax_perc,coef = coef,N_t = N_t,tt = tt,data_out = data_out,baseline = baseline,pos_fo = pos_fo)
ax_2 <- perc_extraction (toextract = 0.02,antivax_perc = antivax_perc,coef = coef,N_t = N_t,tt = tt,data_out = data_out,baseline = baseline,pos_fo = pos_fo)
ax_3 <- perc_extraction(toextract =  0.1,antivax_perc = antivax_perc,coef = coef,N_t = N_t,tt = tt,data_out = data_out,baseline = baseline,pos_fo = pos_fo)
ax_4<- perc_extraction(toextract =  0.15,antivax_perc = antivax_perc,coef = coef,N_t = N_t,tt = tt,data_out = data_out,baseline = baseline,pos_fo = pos_fo)


####        Plots          ####
plot_ax_1 <- plot_totalCases(ax_1[[1]],ax_1[[2]], "5%",N_t,y_up = 30,show_coef =  0.292) +  
  theme(legend.position = "none", plot.title = element_text(hjust=0.5)) 

plot_ax_2 <- plot_totalCases(ax_2[[1]],ax_2[[2]], "2%",N_t,y_up = 30,show_coef =  0.292) +
  theme(legend.direction = "horizontal",plot.title = element_text(hjust=0.5),legend.justification="center") 
  
plot_ax_3 <- plot_totalCases(ax_3[[1]],ax_3[[2]], "10%",N_t,y_up = 30,show_coef =  0.292) +
  theme(legend.position = "none",plot.title = element_text(hjust=0.5))

plot_ax_4 <- plot_totalCases(ax_4[[1]],ax_4[[2]], "15%",N_t,y_up = 30,show_coef =  0.292) +  
  theme(legend.position = "none",plot.title = element_text(hjust=0.5))


legend1 <- get_legend(plot_ax_2)
legend2 <- get_legend(plot_anti)


# Grouped plot:
x_tittle <- ggdraw() + draw_label("Mixing coefficient",size = 15) #ggplot() +labs(title = "Anti vaccination population") + theme(plot.title = element_text(hjust=0.5))
y_tittle <- ggdraw() + draw_label("Additional yearly cases \n per 10 000 people", angle= 90, size=15)
top_plots <- plot_grid(plot_ax_2 + theme(legend.position = "none",axis.title.y = element_blank(),axis.title.x = element_blank(),axis.text.x = element_blank()),
                       plot_ax_1  + theme (axis.text.y = element_blank(),axis.title.y = element_blank(),axis.text.x = element_blank(),axis.title.x  = element_blank()),
                       plot_ax_3  + theme (axis.title.y = element_blank(),axis.title.x = element_blank()),
                       plot_ax_4 + theme (axis.text.y = element_blank(),axis.title.y = element_blank(),axis.title.x = element_blank()),
                       rel_widths = c(1, 1,1 ,1),align = "h", nrow=2)
right  <- plot_grid(top_plots, x_tittle,legend1, ncol = 1, rel_heights = c(1,0.1,0.05))
total_plot <-plot_grid(y_tittle, right, nrow =1, rel_widths = c(0.07,1))


## Pro- anti plots  

# Colors 
cols <- c("0.02" = "#00AFBB", "0.05" ="#E7B800", "0.1" = "#FC4E07", "0.15" = "#293352")
cols_pro <- c("0.98" = "#00AFBB", "0.95" ="#E7B800", "0.9" = "#FC4E07", "0.85" = "#293352")


list_anti <- list(ax_1[[2]], 0.05, ax_2[[2]], 0.02, ax_3[[2]],0.1, ax_4[[2]], 0.15)

plot_anti <- plot_cases_many(list_anti,N_t) +  geom_vline(xintercept = 0.292,aes(colour="grey85"), linetype="dashed") +
    ggtitle("")+  theme_cowplot(font_size = 12) +
  scale_color_manual(values=cols, labels= c("2%","5%","10%","15%"))+
  theme (legend.direction = "horizontal",legend.position = "bottom",legend.justification="center") 

list_pro <- list(ax_1[[1]], 0.95, ax_2[[1]], 0.98, ax_3[[1]],0.9, ax_4[[1]], 0.85) 

plot_pro <- plot_cases_many(list_pro,N_t)  +  geom_vline(xintercept = 0.292,aes(colour="grey85"), linetype="dashed")+ 
  ggtitle ("") + theme_cowplot(font_size = 12) +
  scale_color_manual(values=cols_pro, labels= c("15%","10%","5%","2%")) +
  theme (legend.positio = "none") 


legend_pro <- get_legend (plot_anti)

#Horizontal 
two_plots <- plot_grid(plot_anti + theme(legend.position = "none", axis.title.y = element_blank(),axis.title.x = element_blank()),
                       plot_pro + theme (axis.title.y = element_blank(),axis.title.x = element_blank() ), ncol=2 )
total_two_right <- plot_grid(legend_pro ,two_plots, x_tittle, ncol = 1, rel_heights = c(0.05,1,0.1))

total_two <- plot_grid(y_tittle, total_two_right, nrow =1, rel_widths = c(0.07,1))

#Vertical

two_vertical <-  plot_grid(plot_pro, 
                           plot_anti+ theme(legend.position = "none"),
                           legend_pro,
                           rel_heights = c(1,1,0.1),
                           nrow=3)

two_vertical


### Results to document ####

#5%
total_1 <- ax_1[[1]]
total_1$X2 <- ax_1[[1]]$X2 + ax_1[[2]]$X2

val_1 <- total_3$X2/total_1$X2


#10%
total_3 <- ax_3[[1]]
total_3$X2 <- ax_3[[1]]$X2 + ax_3[[2]]$X2
val_3 <- total_4$X2/total_1$X2

#15%
total_4 <- ax_4[[1]]
total_4$X2 <- ax_4[[1]]$X2 + ax_4[[2]]$X2
val_4 <- (total_4$X2[59]/N_t) *10000

##################################################################################################
#                           Figure S1 : First outbreak dynamics                                 #
#################################################################################################

coef_list <- c (0.0,0.02,0.05,0.1,0.4)
av_list <- c(0.01,0.05,0.1,0.3,0.5)

sup_plot <- plot_dynamics_fo (coef_list,av_list, N_t,antivax_perc,coef,pos_fo,data_out)+xlim(0,17)

