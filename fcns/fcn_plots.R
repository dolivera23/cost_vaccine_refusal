
##################################################################################################
#                                             Plots                                              #
#################################################################################################

# 1. Plot for Adittional annual cases per 10.000 people (PER CAPITA IN EACH POP)
# @param  pList List comes in format of [database1, anti-vaccination perc1, database2, anti-vaccination perc2 ...]
# @patam N_t total population 
plot_cases_many <- function (pList,N_t){
  
  total<- length(pList)/2
  
  db <- data.frame (matrix(ncol = 3, nrow = 0))     # Three columns: Time, extra cases, population
  
  colnames(db) <- colnames(pList[[1]])
  c= 1
  
  for (i in 1:total){
    
    perc <- pList[[c+1]]
    
    to_bind <- pList[[c]]
    to_bind[,"X2"] <- (to_bind[,"X2"] /(N_t*(perc))) *10000
    to_bind[,"V3"] <- perc
    
    db <- rbind(db, to_bind)
    
    c <- c + 2
    
  }
  
  db[,3] <- as.factor(db[,3])
  
  p_out <- ggplot(data=db, aes(x = X1, y = X2, colour=V3)) + geom_line( size=1) +
    xlab("Mixing coefficient") +  ylab("Additional annual cases \n  per 10 000 people") + 
    labs(colour= "Anti vaccination population") + 
    return(p_out)
  
}

# 2. Plot additional yearly cases for overall population and show shaded are by pop
# @param dbp: Pro vax data fram. Three columns: coef,extra casesfirst outbreak, pop
# @param dba: Anti vax data fram. Three columns: coef,extra casesfirst outbreak, pop
# @param y_up: Upper limit for y xis
# @param show_coef: Coefficient vertical line value 
plot_totalCases <- function (dbp,dba,title,N_t, y_up, show_coef) {
  
  binded <- rbind(dbp, dba)
  colnames(binded)<- c("Coef", "Cases", "Population")
  binded[,"Cases"] <- (binded[,"Cases"]/N_t)*10000
  
  r<- ggplot(binded, aes(x=Coef, y=Cases, fill=Population)) + geom_area() + 
    ggtitle (title) + ylab("Additional yearly cases \n per 10 000 people") + xlab("Mixing coeficient") +
    scale_fill_manual(name="", values = wes_palette("Moonrise2", n = 2),labels = c("Pro Vaccination", "Anti Vaccination")) +
    ylim(0,y_up) + 
    theme_cowplot(font_size = 15) + 
    geom_vline(xintercept = show_coef,aes(colour="grey85"), linetype="dashed")
  
  return(r)
} 


# 3. Plot dynamics of first outbreak for a combination of coef and antivax populations
# @param coef_list List of coeficients to plot
# @param av_list List of antivaccination population to plot
#
# @out Plot of first 10 year dynamics. With vetical line showing first outbreak time. 

plot_dynamics_fo <- function (coef_list,av_list, N_t,antivax_perc,coef,pos_fo,data_out){
  
  
  toplot <- data.frame(matrix(0,nrow=0,ncol=5))
  fo_plot <- data.frame(matrix(0,nrow=0,ncol=3))
  
  for (i in 1:length(coef_list)) {
    
    p_coef <- coef_list[i]
    
    for (j in 1:length(av_list)){
      
      p_av <- av_list[j]
      
      av <- which.min(abs(antivax_perc-p_av))
      c <-  which.min(abs(antivax_perc-p_coef))
      
      df <- data_out[[c]][[av]] 
      df <- df[,c("t","inc.1","inc.2")]
      
      df_long = melt(df, id.vars = "t")
      colnames(df_long) <- c("t", "Population","value")
      df_long$Population <- factor(df_long$Population, levels = c("inc.2", "inc.1"))
      levels(df_long$Population) <- c("Anti Vaccination", "Pro vaccination")
      
      df_long$antiVaccination <- p_av*100#paste(p_av*100)
      df_long$mixing <-p_coef
      
      toplot <- rbind(toplot,df_long)
      
      # Firs outbreak information 
      temp <- data.frame(matrix(0,nrow=1,ncol=0))
      temp$fo <- pos_fo[c,av]/52
      temp$antiVaccination <- p_av*100#paste(p_av*100, "%")
      temp$mixing <-p_coef
      
      fo_plot <- rbind(fo_plot,temp)
    }
  }
  
  fo_plot <- fo_plot %>% group_by(antiVaccination,mixing) %>% summarise(mean_x = mean(fo))

  
  plot <- ggplot(toplot, aes(x = t/52, y = (value/N_t)*10000, fill = Population)) + geom_area()+#geom_line(size=1) + 
    xlab("Time (years)") + ylab("Weekly incidence per \n 10000 individuals")+
    scale_fill_manual(values = rev(wes_palette("Moonrise2", n = 2))) + xlim(0,10) +
    theme_bw(base_size = 14) +
    geom_vline(data  = fo_plot, aes(xintercept = mean_x), linetype="dashed",colour="#34568B") +
    facet_grid(antiVaccination~mixing,scales = "free")+
    theme(legend.position = "bottom",
          strip.background = element_rect(fill = NA),
          panel.border = element_blank(),
          axis.line = element_line(),
          legend.text.align = 0)
  
  return(plot)
  
}

# 4. Plot yearly incidence for the overal population.
plot_incidence <- function (df, N_t, yscale){
  
  df <- df %>% filter(population == "Total") %>% select (year, cases)
  
  plot_out <- ggplot(df, aes(x=year, y=(cases/N_t)*yscale,color=population) )+ 
    geom_line(size=1.2, color="#00758F") + 
    labs (y= "Additional yearly cases", x="Years")+
    theme_bw(base_family = 15)+
    theme(strip.background = element_rect(fill = NA),
          panel.border = element_blank(),
          axis.line = element_line(),
          legend.text.align = 0)
  
  color ="#CCC591"
  
  return (plot_out)
  
}

# 5. Plot yearly deaths for the overal population. With upper and lower bounds shadowed 
plot_deaths <- function (df,N_t,yscale){
  
  df <- df %>% filter(population == "Total") %>% select (year,deaths, deaths_high,deaths_low)
  
  plot_out <- ggplot(df, aes(x=year, y=(deaths/N_t)*yscale)) + 
    geom_ribbon(aes(ymin = (deaths_low/N_t)*yscale, ymax = (deaths_high/N_t)*yscale), fill = "#efeddc") +
    #geom_errorbar(aes(ymin= (deaths_low/N_t)*yscale, ymax=(deaths_high/N_t)*yscale), width=.7,color= "#efeddc")+  
    geom_point(size=3, shape=21, fill="white",color="#efeddc")+
    geom_line(color ="#8c8340",size=1.2) + 
    labs (y= "Additional yearly deaths", x="Years")+
    theme_bw(base_family = 13)+
    theme(strip.background = element_rect(fill = NA),
          panel.border = element_blank(),
          axis.line = element_line(),
          legend.text.align = 0)
  
  return (plot_out)
  
}

# 6. Bar plot of cumulative incidence  per population: Pro vaccination and anti-vaccination 
plot_cumu_inc <- function (df, N_t, yscale,last_year){
  
  
  df <- df %>% filter(population != "Total" &  year== last_year) %>%
    select (year,cum_inc,population)#,cost, population)
  
  df$population <- factor(df$population)
  df$population <- fct_rev(df$population)
  
  plot <- ggplot(data = df, aes(x = population, y = (cum_inc/N_t)*yscale), fill = population)+
    geom_col(aes(fill = population))+
    labs(y = "Total cases", x = "", fill = "Population") +
    theme_bw(base_size = 13) +
    theme(strip.background = element_rect(fill = NA),
          panel.border = element_blank(),
          axis.line = element_line(),
          legend.text.align = 0) +
    scale_y_continuous(labels = comma)+
    scale_fill_manual(values = wes_palette("Moonrise2", n = 2))
  
  return(plot)
  
}

# 7. Bar plot with confidence interval for cumulative deaths per population: Pro-Vaccination and anti-vaccination 
plot_cumu_deaths <- function (df, N_t, yscale,last_year){
  
  
  df <- df %>% filter(population != "Total" &  year== last_year) %>%
    select (year,cum_death, cum_death_low, cum_death_high,population)#,cost, population)
  
  df$population <- factor(df$population)
  df$population <- fct_rev(df$population)
  
  plot <- ggplot(data = df, aes(x = population, y = (cum_death/N_t)*yscale), fill = population)+
    geom_col(aes(fill = population))+
    geom_errorbar(data=df, aes(x= population ,ymin= (cum_death_low/N_t)*yscale, ymax=(cum_death_high/N_t)*yscale,group=1), width=.1 , inherit.aes = FALSE)+
    geom_point(data=df, aes(x=population, y = (cum_death/N_t)*yscale,group=1),inherit.aes =  FALSE)+
    labs(y = "Total deaths", x = "", fill = "Population") +
    theme_bw(base_size = 13) +
    theme(strip.background = element_rect(fill = NA),
          panel.border = element_blank(),
          axis.line = element_line(),
          legend.text.align = 0) +
    scale_y_continuous(labels = comma)+
    scale_fill_manual(values = wes_palette("Moonrise2", n = 2))
  
  return(plot)
  
}


#8.  Bar plot of average costs differentiated by cost type: QUALY deaths, QALY morbidity.
#    healthcarecosts, productivity costs. 

plot_bar_costs <- function (df, N_t, yscale){
  
  my_col <- c( "#b2b2b2",
               "#f4e1d2",
               "#f18973",
               "#c97b69")
  
  library(scales)
  
  bound <- df  %>% filter(population == "Total") %>% select (year,cost, cost_low,cost_high) 
  
  df <- df %>% filter(population == "Total")  %>% 
    ungroup() %>% 
    select (year,qaly_death,qaly_dis,cost_healthcare,cost_prod)#,cost, population)
  
  long_df <- melt(df,id.vars = "year")
  
  plot <- ggplot(data = long_df, aes(x = year, y = (value/N_t)*yscale, 
                                     fill = variable))+ 
    geom_col()+
    geom_errorbar(data=bound, aes(x= year ,ymin= (cost_low/N_t)*yscale, ymax=(cost_high/N_t)*yscale,group=1 ),color="grey53", width=.1 , inherit.aes = FALSE)+
    geom_point(data=bound, aes(x=year, y = (cost/N_t)*yscale,group=1),inherit.aes =  FALSE)+
    labs(y = "Total costs (GBP)", x = "Year", fill = "Costs") +
    theme_bw(base_size = 13) +
    theme(strip.background = element_rect(fill = NA),
          panel.border = element_blank(),
          axis.line = element_line(),
          legend.text.align = 0) +
    scale_y_continuous(labels = comma)+
    scale_fill_manual(values= my_col, labels= c("Deaths", "Disease burden","Healthcare","Productivity loss"))
  
  return(plot)
  
}


# 9. Plot of cumulative costs over time showns for the average costs and discretised by population 
plot_cum_costs <- function (df,N_t,yscale){
  
  df <- df %>% filter(population != "Total") 
  
  df$population <- fct_rev(df$population)
  
  plot <- ggplot(data = df, aes(x = year, y = (cum_costs/N_t)*yscale, 
                                fill = population, label = round((cum_costs/N_t)*yscale,2)))+
    geom_area() +
    labs(y = "Cumulative costs (GBP)", x = "Time", fill = "Population") +
    scale_fill_manual(values = wes_palette("Moonrise2", n = 2)) +
    theme_bw(base_size = 13) +
    theme(strip.background = element_rect(fill = NA),
          panel.border = element_blank(),
          axis.line = element_line(),
          legend.text.align = 0) +
    scale_y_continuous(labels = comma)
  
  return(plot)
  
}

