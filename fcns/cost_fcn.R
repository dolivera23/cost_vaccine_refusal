
# Function that returns costs for a specific population in an specific year 
# @param age_distr: Proportion of cases per age group. 
# @param names: Column names of data-frame
# @param incidence: Yearly incidence
# @param population: "Pro Vaccination", "Anti Vaccination" or "Total"
# @param year: Year to estimate the costs 
# @param CFR: Caste fatality ratio per age group with lower and upper bounds. 
# @param QALY_death: Moratiliyy QALY per age group.
# @param value_qaly: Cost in GBP per QALY
# @param case_qaly: QALY lost due to morbidity during disease (short term)
# @param p_productiviy: Productivity loss Cost in GBP per case 
# @param p_healthcare: Heatlh care costs in GBP per case 


cost_year <- function (age_distr,names,incidence,population,year,
                       CFR,QALY_death,value_qaly,case_qaly,p_healthcare,p_productiviy){
  
  out <- data.frame(matrix(0, nrow=nrow(age_distr), ncol=14))
  colnames(out) <- c("Age_group", names[3:15])

  
  out$Age_group <- age_distr$Age_group
  out$cases <- age_distr$PointEst *incidence
  
  out$deaths_low <- out$cases * CFR$Lower
  out$deaths <- out$cases * CFR$PointEst
  out$deaths_high <- out$cases * CFR$Upper
  
  out$qaly_death_low <- out$deaths_low * QALY_death[,2] * value_qaly
  out$qaly_death <- out$deaths * QALY_death[,2] * value_qaly
  out$qaly_death_high <- out$deaths_high * QALY_death[,2] * value_qaly
  
  
  out$qaly_dis <- out$cases * case_qaly * value_qaly
  
  out$cost_healthcare <- out$cases * p_healthcare
  out$cost_prod <- out$cases * p_productiviy
  
  sub_total <- out$qaly_dis  +  out$cost_healthcare + out$cost_prod
  
  out$cost_low <- out$qaly_death_low + sub_total
  out$cost <-  out$qaly_death  + sub_total
  out$cost_high <- out$qaly_death_high  + sub_total
  
  
  total <- data.frame(t(colSums(out[,2:14])))
  total$population <- population
  total$year <- year
  
  return  <- total
  
  
}


# Function that returns yearly costs per population (anti and pro vaccination).
# @param df: Data Frame with weekly dynamics (output model)
# @param CFR: Caste fatality ratio per age group with lower and upper bounds. 
# @param QALY_death: Moratiliyy QALY per age group.
# @param years_total: Number of years analysed (20 years default)
# @param age_distr: Proportion of cases per age group. 
# @param case_qaly: QALY lost due to morbidity during disease (short term)
# @param value_qaly: Cost in GBP per QALY
# @param p_productiviy: Productivity loss Cost in GBP per case 
# @param p_healthcare: Heatlh care costs in GBP per case 
costs_output <- function (df,CFR,QALY_death, years_total=20, age_distr, case_qualy, value_qaly,
                          p_productiviy,p_healthcare){
  
  ## -------------------------------------------------------- ##
  ##                   Output
  ## ------------------------------------------------------ ##
  
  # Data frames with complete information 
  df_out <- data.frame(matrix(0,nrow=0, ncol =15))
  names<- c("year","population", "cases", "deaths_low", "deaths", "deaths_high", "qaly_death_low", "qaly_death",
            "qaly_death_high", "qaly_dis","cost_healthcare", "cost_prod", "cost_low", "cost","cost_high")
  colnames(df_out) <- names
  
  ## -------------------------------------------------------- ##
  ##                   Yearly Incidence 
  ## ------------------------------------------------------ ##
  
  # Estimating yearly incidence and compared to baseline scenario. 
  df_yearly <- df
  df_yearly$total <- df_yearly$inc.1 + df_yearly$inc.2 
  
  df_yearly[, "year"] <- ceiling(df[,"t"] /52)
  
  df_yearly <- aggregate(df_yearly[,c("inc.1","inc.2","total")], by=list(Category=df_yearly$year), FUN=sum)
  
  df_yearly <- df_yearly[2:(years_total+1),]   # Selectin only years we are interested in 
  
  ## -------------------------------------------------------- ##
  ##                   Values per age group and year
  ## ------------------------------------------------------ ##
  
  
  # Estimating values for each group
  
  for (i in 1:years_total){
    
    df_year <- df_yearly[i,]
    
    # Pro vaccination 
    total_pro <- cost_year(age_distr = age_distr,names = names,incidence =df_year$inc.1,
                           population = "Pro vaccination",year = i,CFR = CFR,QALY_death = QALY_death,
                           value_qaly = value_qaly,case_qaly = case_qaly,p_healthcare = p_healthcare,
                           p_productiviy =p_productiviy)
    
    # Anti Vaccination
    
    total_anti <- cost_year(age_distr = age_distr,names = names,incidence =df_year$inc.2,
                            population = "Anti vaccination",year = i,CFR = CFR,QALY_death = QALY_death,
                            value_qaly = value_qaly,case_qaly = case_qaly,p_healthcare = p_healthcare,
                            p_productiviy =p_productiviy)
    
    #Total
    total_total  <- cost_year(age_distr = age_distr,names = names,incidence =df_year$total,
                              population = "Total",year = i,CFR = CFR,QALY_death = QALY_death,
                              value_qaly = value_qaly,case_qaly = case_qaly,p_healthcare = p_healthcare,
                              p_productiviy =p_productiviy)
    
    df_out<- rbind(df_out,total_pro)
    df_out<- rbind(df_out,total_anti)
    df_out<- rbind(df_out,total_total)
    
    
  }
  return (df_out) 
}


