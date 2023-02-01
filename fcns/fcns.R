##################################################################################################
#                                    Dynamic model functions                                    #
#################################################################################################

##  This function creates the model and returns a data frame with the model outputs.
##   @param path  Path to the model to be run
##   @param tt    Vector of the time points where the model is going to be run

runModel <- function (path, tt, np, N, rec, ro, I0, alpha, ...) {
  ##-------- Packages required -------- ##
  
  library("odin")
  library("here")
  library("ggplot2")
  library("reshape2")
  
  # Path to  model
  path_mult <- here::here(path)
  
  
  # Model generator
  model_generator <- odin::odin(path_mult)
  mod <- model_generator(np = np,N = N, rec = rec, ro = ro, alpha = alpha,I0 = I0,...)
  
  # Running the model
  out_model <- mod$run(tt, maxsteps=10000000, method="ode45")
  

  out <- mod$transform_variables(out_model)

  
  #out_model <- data.frame(out_model)
  
  return (out)
  
}


# Function that calculare the transmission matrix for a 2 population model
#   @param N1   Population size patch 1
#   @param N2   Population Size patch 2
#   @param G    Moving pattern parameter. It is assumed both populations have the same moving patterns. Value between 0 and 0.5

trans_matrix <- function (N1,N2,G) {
  
  a <- matrix(0, nrow = 2,  ncol = 2,byrow = TRUE) 
  
  a[1,1] <- (1-G)^2/((1-G)*N1 + G*N2) + G^2/(G*N1 + (1-G)*N2)
  a[1,2] <- ((1-G)*G)/((1-G)*N1 + G*N2) + ((1-G)*G)/(G*N1 + (1-G)*N2)
  a[2,1] <- ((1-G)*G)/((1-G)*N1 + G*N2) + ((1-G)*G)/(G*N1 + (1-G)*N2)
  a[2,2] <- (1-G)^2/((1-G)*N2 + G*N1) + G^2/(G*N2 + (1-G)*N1)

  return(a)
  
}


