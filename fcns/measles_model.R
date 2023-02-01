#Math model for two population. Population can vary in size. 

np<- user(2)                      #  Number of populations to be modelled 

# -------- Differential equations for the SIR model -------- #
dim(S) <- np
dim(I) <- np
dim(R) <- np
deriv(S[1:np]) <- B[i]*(1-(ef*cov[i]))*N[i] - lambda[i]*S[i]- miu[i] *S[i]         
deriv(I[1:np]) <- lambda[i]*S[i] - rec[i]*I[i] - miu[i]*I[i]
deriv(R[1:np]) <- B[i]*ef*cov[i]*N[i] + rec[i]*I[i]  - miu[i]* R[i]             

# -------- Initial conditions -------- #

initial(S[1:np]) <- N[i] - I0[i] - R0[i]
initial(I[1:np]) <- I0[i]
initial(R[1:np]) <- R0[i]

# -------- Outbreak size --------#

dim(It) <- np
deriv(It[1:np]) <- lambda[i]*S[i]
initial (It[1:np]) <- 0


inc[] <- lambda[i]*S[i]

# -------- Parameters -------- #

# User parameters
N[] <- user()                #  Population in each patch 
B[] <- user()                #  Birth rate in each patch
miu[] <- user()              #  Death rate in each patch
rec[] <-user()               #  Recovery rate in each patch 
ro[] <- user()               #  Reproductive number in each patch  
cov[] <- user()              #  Vaccination coverage in each patch 
ef <- user()                 #  Vaccine effectiveness


# Force of infection
dim(beta) <- np
beta[] <- ro[i] * (rec[i] + miu[i])

#  Individual rate of infection on each patch   
alpha[,] <- user()           #  Strength of transmission between populations (patches)


dim(rows) <- c(np,np)
rows[,] <- alpha[i,j]*I[j]
lambda[] <- (beta[i]) * sum(rows[i,])

# Initial conditions
I0[] <- user()
R0[] <- user()

# -------- Dimensions -------- #
dim(rec) <- np
dim(cov) <- np
dim(N) <- np
dim(B) <- np
dim(miu) <- np
dim (alpha) <- c(np,np)
dim(ro) <-np
dim(lambda) <- np
dim (I0) <- np
dim(R0) <- np
dim(inc) <- np



# -------- Outputs -------- #

output(inc[]) <- TRUE
#output(inc[]) <- inc[i]

config(base) <- "measles_model"