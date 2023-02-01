## Screening functions


##################################################################################################
#                   Finding steady state for initial conditions                                  #
#################################################################################################

initial_cond <- function(years,path_model,npop,B,miu,cov,rec,ro,ef, N_t){
  
  tt <- seq(0, years*52, length.out =years*52)
  
  N<- c(0.5*N_t, 0.5*N_t)
  alpha <- trans_matrix(N[1],N[2],0.5)
  
  ss<- runModel( path_model, tt, np = npop, N = N, B=B,miu = miu, cov= cov,rec = rec,ro = ro,ef=ef,alpha = alpha, I0 = c(1,0),R0 = c(0,0))
  ss <- data.frame(ss)
  
  ini <- ss[52*years,]
  
  return(ini) 
  
}

##################################################################################################
#                   Screening different pop size and mixing patterns                            #
#################################################################################################

## @out   Data_out[[coef]][[Antivax_perc]]

screening <- function (years, tt,antivax_perc, coef,path_model,npop,B,miu,cov,rec,ro,ef, N_t, I0_t,R0_t){
  
  data_out <- vector(mode = "list", length = length(coef))  ## Data_out[[coef]][[Antivax_perc]]
  
  for (j in 1:length(coef)) {
    
    one_coef <- vector(mode = "list", length = length(antivax_perc))
    
    for (i in 1:length(antivax_perc)){
      
      av <- antivax_perc[i]         # Anti vaxxer percentage
      N <- c(N_t*(1-av), N_t*av)    # Total population in each patch
      alpha <- trans_matrix(N[1],N[2],coef[j])
      I0 <- c(I0_t*(1-av), I0_t*av)
      R0 <- c(R0_t*(1-av), R0_t*av)
      
      
      out <- runModel( path_model, tt, np = npop, N = N, B=B,miu = miu, cov= cov,rec = rec,ro = ro,ef=ef,alpha = alpha, I0 = I0,R0 = R0)
      out <- data.frame(out)
      
      one_coef[[i]] <- out  
    }
    data_out[[j]] <- one_coef
  }
  
  return (data_out)
}


##################################################################################################
#                           Extracting first outbreak ending time                                #
#################################################################################################

# Function that return the position of the first minimum after the first peak: 
# The position where the first outbreak ends

first_outbreak_time <- function (db){
  
  total <- db[,"inc.1"] + db[,"inc.2"]
  
  d1<- diff(total)
  d2<- diff(sign(d1))
  out <- which(d2!= 0)[2] 
  
  return (out)
}


##################################################################################################
#                           Extracting data by antivax percentage                                #
#################################################################################################

# Function to extract the extra cases at the first oubreak for each popultion. 
# @param toextract: Anti vaccinatio poercentaje to extract 
# @param antivax_perc: Vector wuth all anti vacccination percentages
# @param coef: Vector with all mixing coefficents.
# @param data_out: Object will ALL simulations 
#
#@return out: List with two objects: One data frame for pro vaccination pop[1] and one data frame for anti vaccination pop[2]

perc_extraction <- function (toextract, antivax_perc, coef, N_t, tt, data_out, baseline, pos_fo){
  
  pos <- which.min(abs(antivax_perc-toextract))
  
  #First outbreak 
  fo_pro_df <-  data.frame(matrix(0,nrow= length(coef),ncol=2))                # Three columns: Coef, Time First outbreak and Incidence/person compared to baseline scenario
  fo_pro_df[,1] <- coef
  fo_anti_df <-  data.frame(matrix(0,nrow= length(coef),ncol=2))
  fo_anti_df[,1] <- coef
  
  
  # Filling the maxtrices
  for (i in 1:length(coef)) {
    
    pos_out <- pos_fo[i,pos]                    # Position of first outbreak 
    
    BL <- baseline[pos_out,"inc.1"]             # Incidence at baseline scenario at time of first outbreak
    BL_pro <- BL*(1-antivax_perc[pos])
    BL_anti <- BL* antivax_perc[pos]
    
    out <- data_out[[i]][[pos]]                 # Simulation results 
    
    fo_pro_df [i,2] <-((out[pos_out,"It.1"]/out[pos_out,"t"]) - BL_pro)*52     #((out[pos_out,8])- (out[pos_out,1]*BL_pro))/(out[pos_out,1]*N_t) #
    fo_anti_df [i,2] <-((out[pos_out,"It.2"]/out[pos_out,"t"]) - BL_anti)*52    #(out[pos_out,9])- (out[pos_out,1]*BL_anti)/(out[pos_out,1]*N_t)#
    
  }
  
  fo_pro_df [, 3] <- "Pro Vaccine"
  fo_anti_df [,3] <- "Anti Vaccine"
  
  out <- vector(mode = "list", length = 2)
  
  out[[1]] <- fo_pro_df
  out[[2]] <- fo_anti_df
  
  return (out)
  
}


##################################################################################################
#                           Extracting data by coefficient                                      #
#################################################################################################
# Function to extract the extra cases at the first oubreak for each popultion. 
# @param coef_selected: Mixing coeficient to extract 
# @param antivax_perc: Vector wuth all anti vacccination percentages
# @param coef: Vector with all mixing coefficents.
# @param data_out: Object will ALL simulations 
#
#@return out: List with two objects: One data frame for pro vaccination pop[1] and one data frame for anti vaccination pop[2]


data_by_coef <- function(coef_selected, data_out,antivax_perc,coef, pos_fo,baseline,tt, N_t){
  
  db_fo <- data.frame(matrix(0,nrow= length(antivax_perc),ncol=3))
  
  db_fo[,1] <- antivax_perc
  
  colnames(db_fo) <- c("Antivacc_perc", "Pro Vaccine", "Anti Vaccine")
  
  pos <- which.min(abs(coef - coef_selected))
  
  for (i in 1:length(antivax_perc)){
    
    pos_out <- pos_fo[pos,i]
    
    BL <- baseline[pos_out,"inc.1"]
    BL_pro <- BL*(1-antivax_perc[i])
    BL_anti <- BL* antivax_perc[i]
    
    out<- data_out[[pos]][[i]]
    
    pro_fo <- ((out[pos_out,"It.1"]/out[pos_out,"t"]) - BL_pro)*52  
    anti_fo<- ((out[pos_out,"It.2"]/out[pos_out,"t"]) - BL_anti)*52
    
    db_fo[i,2:3] <- c(pro_fo, anti_fo)
  }
  
  return(db_fo)
  
}

##################################################################################################
#                          Deaths by age                                                        #
#################################################################################################

deaths_by_age <- function (df_year,age_distr,CFR,year){
  
  out <-   vector(mode = "list", length = 4)
  
  names <- c("Age_group","cases", "value", "value_low", "value_high")
  
  out_pro <- data.frame(matrix(0, nrow=nrow(age_distr), ncol=5))
  colnames(out_pro) <- names
  
  out_pro$Age_group <- age_distr$Age_group
  out_pro$cases <- age_distr$PointEst * df_year$inc.1
  out_pro$value <- out_pro$cases * CFR$PointEst
  out_pro$value_low <- out_pro$cases * CFR$Lower
  out_pro$value_high <- out_pro$cases * CFR$Upper
  
  out[[1]] <- out_pro
  
  out_anti <- data.frame(matrix(0, nrow=nrow(age_distr), ncol=5))
  colnames(out_anti) <- names
  
  out_anti$Age_group <- age_distr$Age_group
  out_anti$cases <- age_distr$PointEst * df_year$inc.2
  out_anti$value <- out_anti$cases * CFR$PointEst
  out_anti$value_low <- out_anti$cases * CFR$Lower
  out_anti$value_high <- out_anti$cases * CFR$Upper
  
  out[[2]] <- out_anti
  
  
  out_total <- data.frame(matrix(0, nrow=nrow(age_distr), ncol=5))
  colnames(out_total) <- names
  out_total$Age_group <- age_distr$Age_group
  out_total[,2:5] <- out_pro[,2:5] + out_anti[,2:5]   
  
  out[[3]] <- out_total 
  
  
  return (out) 
  
}

## Function that converst incidence data from weeks to years

to_years <- function (df){
  
  df_yearly <- df
  df_yearly[, "year"] <- ceiling(df[,"t"] /52)
  df_yearly <- aggregate(df_yearly[,c("inc.1","inc.2","total")], by=list(Category=df_yearly$year), FUN=sum)
  
}


##################################################################################################
#                          Plotting functions                                                   #
#################################################################################################


plot_yearly <- function (df, N_t, ylabel,av_pop,pv_pop, col_min, col_max){
  
  df[which(df$Population == "Total"), col_min:col_max] <- (df[which(df$Population == "Total"), col_min:col_max]/N_t) * 1e6
  df[which(df$Population =="Anti vaccination"), col_min:col_max] <- (df[which(df$Population == "Anti vaccination"), col_min:col_max]/av_pop) * 1e6
  df[which(df$Population == "Pro vaccination"), col_min:col_max] <- (df[which(df$Population == "Pro vaccination"), col_min:col_max]/pv_pop) * 1e6
  
  plot_out <- ggplot(df, aes(x=year, y=value, col= Population)) + geom_errorbar(aes(ymin= value_low, ymax=value_high), width=.7) +  
    geom_point(size=3, shape=21, fill="white") + xlab("Time (years)") + ylab(paste(ylabel," per million\n individuals")) +
    theme_tufte(base_family = 20) + scale_color_manual(values = wes_palette("Moonrise2", n = 3))+
    geom_line(size=1.2) 
  
  return (plot_out)
  
}

plot_casesEach <- function(dbp,dba, title, N_a, N_p){
  
  r <- ggplot() + geom_line( data=dbp, aes(x = X1, y = dbp[,2]/N_p, colour="Pro Vaccination"), size=1) +  geom_line( data=dba, aes(x = X1, y = dba[,2]/N_a, colour="Anti Vaccination"), size=1) +
    xlab("Mix Coef") + ylab("Additional annual \n cases per person")+ ggtitle (title)  +
    scale_colour_manual("Population", breaks = c("Pro Vaccination", "Anti Vaccination"),values = wes_palette("Moonrise2", n = 2))
  return(r)
}

plot_total_AV<- function (db, title, N_t){
  
  probs1 <- melt(db, id.vars = "Antivacc_perc")
  probs1[,"value"] <- probs1[,"value"] / N_t
  
  probs1$variable <- factor(probs1$variable, levels = rev(levels(probs1$variable)))
  r <-  ggplot(probs1, aes(x = Antivacc_perc, y = value, fill = variable)) + geom_area() + 
    ggtitle (title) +xlab("Antivaxer proportion") + ylab("Additional annual \n cases per person") +
    scale_fill_manual(values = wes_palette("Moonrise2", n = 2))
  
  return(r)
} 

plot_dynamics <- function (df, N_t, av){
  
  out <-  vector(mode = "list", length = 3)
  
  #av_pop <- N_t*av
  #pv_pop <- N_t*(1-av)
  
  # df <- df[,c("t","inc.1","inc.2")]
  # df[, "total"] <- ((df[,"inc.1"] + df[,"inc.2"]) /N_t)*1e6
  # df[,"inc.1"] <- (df[,"inc.1"] /(pv_pop))*1e6
  # df[,"inc.2"] <- (df[,"inc.2"]/(av_pop))*1e6
  
  
  df <- df[,c("t","inc.1","inc.2")]
  df[, "total"] <- (df[,"inc.1"] + df[,"inc.2"])
  
  
  df_long = melt(df, id.vars = "t")
  colnames(df_long) <- c("t", "Population","value")
  df_long$Population <- factor(df_long$Population, levels = c("inc.2", "inc.1","total"))
  levels(df_long$Population) <- c("Anti Vaccination", "Pro vaccination", "Total ")
  
  plot <- ggplot(df_long, aes(x = t/52, y = value, col = Population)) + geom_line(size=1) + xlab("Time (years)") + ylab("Weekly incidence per \n million individuals")+
    scale_color_manual(values = wes_palette("Moonrise2", n = 3)) + xlim(0,20)
  
  
  
  return(out)
  
}

plot_dynamics_onevar <-function (df,av_pop, pv_pop, N_t, ylabel, value){
  
  
  df[which(df$population == "Total"), value] <- (df[which(df$population == "Total"), value]/N_t) * 1e6
  df[which(df$population =="Anti vaccination"), value] <- (df[which(df$population == "Anti vaccination"), value]/av_pop) * 1e6
  df[which(df$population == "Pro vaccination"), value] <- (df[which(df$population == "Pro vaccination"), value]/pv_pop) * 1e6
  
  y_col <- which(colnames(df) == value)
  
  out <- ggplot(df, aes(x = year, y = df[,y_col], col = population)) + geom_line(size=1) + xlab("Time (years)") + ylab(paste( ylabel, " per \n million individuals"))+
    scale_color_manual(values = wes_palette("Moonrise2", n = 3)) +  theme_tufte(base_family = 20) 
  return(out)
  
}


plot_dynamics_onevar <-function (df,av_pop, pv_pop, N_t, ylabel, value){
  
  
  df[which(df$population == "Total"), value] <- (df[which(df$population == "Total"), value]/N_t) * 1e6
  df[which(df$population =="Anti vaccination"), value] <- (df[which(df$population == "Anti vaccination"), value]/av_pop) * 1e6
  df[which(df$population == "Pro vaccination"), value] <- (df[which(df$population == "Pro vaccination"), value]/pv_pop) * 1e6
  
  y_col <- which(colnames(df) == value)
  
  out <- ggplot(df, aes(x = year, y = df[,y_col], col = population)) + geom_line(size=1) + xlab("Time (years)") + ylab(paste( ylabel, " per \n million individuals"))+
    scale_color_manual(values = wes_palette("Moonrise2", n = 3)) +  theme_tufte(base_family = 20) 
  return(out)
  
}




