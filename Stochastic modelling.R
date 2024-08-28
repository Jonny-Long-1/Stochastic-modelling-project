library(tidyverse)
library(Hmisc)

species <- c("TnB","TcmB","TemB","TeffB","CD19+B","TnP","TcmP","TemP","TeffP","CD19+P")

#Defining the stoicheometric matrix

#The reactions representing the natural elimination rate of CAR-T cells in the blood 
react_elim <- c(-1,rep(0,9),
                0,-1,rep(0,8),
                0,0,-1,rep(0,7),
                rep(0,3),-1,rep(0,6)) %>% 
  matrix(nrow=10) %>% t()

colnames(react_elim) <- species

#The reactions representing the differentiation of CAR-T cells in the blood
react_diff_B <- c(-1,1,rep(0,8),
                  0,-1,1,rep(0,7),
                  0,0,-1,1,rep(0,6)) %>% 
  matrix(nrow = 10) %>% t()

colnames(react_diff_B) <- species

#The reactions representing the natural proliferation rate of species in the blood
react_nat_prolif_B <- c(1,rep(0,9),
                    0,1,rep(0,8),
                    0,0,1,rep(0,7),
                    rep(0,4),1,rep(0,5)) %>% 
  matrix(nrow = 10) %>% t()

colnames(react_nat_prolif_B) <- species

#The reactions representing tumor cells in the blood being killed by CAR-T cells
react_kill_B <- c(rep(0,4),-1,rep(0,5),
                  rep(0,4),-1,rep(0,5),
                  rep(0,4),-1,rep(0,5),
                  rep(0,4),-1,rep(0,5)) %>% 
  matrix(nrow=10) %>% t()

colnames(react_kill_B)<- species

#The reactions representing the proliferation of CAR-T cells in the blood induced by the presence of tumor cells
react_induce_prolif_B <- c(1,rep(0,9),
                           0,1,rep(0,8),
                           0,0,1,rep(0,7)) %>% 
  matrix(nrow=10) %>% t()

#The reactions representing the migration of species between compartments
react_migration <- c(-1,rep(0,4),1,rep(0,4),
                      1,rep(0,4),-1,rep(0,4),
                     0,-1,rep(0,4),1,rep(0,3),
                     0,1,rep(0,4),-1,rep(0,3),
                     0,0,-1,rep(0,4),1,0,0,
                     0,0,1,rep(0,4),-1,0,0,
                     rep(0,3),-1,rep(0,4),1,0,
                     rep(0,3),1,rep(0,4),-1,0,
                     rep(0,4),1,rep(0,4),-1) %>%
                    #Tumor cells are assumed to be unable to migrate back to peripheral tissue
  matrix(nrow = 10) %>% t()

colnames(react_migration) <- species

#The reactions representing CAR-T cell differentiation in the peripheral tissue
react_diff_P <- c(rep(0,5),-1,1,rep(0,3),
                  rep(0,6),-1,1,0,0,
                  rep(0,7),-1,1,0) %>% 
  matrix(nrow=10) %>% t()

colnames(react_diff_P) <- species

#The reactions representing natural species proliferation in peripheral tissue
react_nat_prolif_P <- c(rep(0,5),1,rep(0,4),
                        rep(0,6),1,rep(0,3),
                        rep(0,7),1,0,0,
                        rep(0,9),1) %>% 
  matrix(nrow = 10) %>% t()

colnames(react_nat_prolif_P) <- species

#The reactions representing tumor cells in peripheral tissue being killed by CAR-T cells
react_kill_P <- c(rep(0,9),-1,
                  rep(0,9),-1,
                  rep(0,9),-1,
                  rep(0,9),-1) %>% 
  matrix(nrow=10) %>% t()

#The reactions representing the proliferation of CAR-T cells in peripheral tissue induced by the presence of tumor cells
react_induce_prolif_P <- c(rep(0,5),1,rep(0,4),
                           rep(0,6),1,rep(0,3),
                           rep(0,7),1,0,0) %>% 
  matrix(nrow = 10) %>% t()

#All the reactions taking place in the blood
blood_reactions <- rbind(react_elim,react_diff_B,react_nat_prolif_B,react_kill_B,react_induce_prolif_B)

#All the reactions taking place in the peripheral tissue
p_tissue_reactions <- rbind(react_diff_P,react_nat_prolif_P,react_kill_P,react_induce_prolif_P)

#The entire stoicheometric matrix
s_mat <- rbind(blood_reactions,react_migration,p_tissue_reactions)

scale <- 5*10^4

s_mat <- scale*s_mat

#Assumptions
blood_vol <- 5*10^6
p_tissue_vol <- 1.75*10^6
tumor_density <- (2/3)*10^6
patient_weight <- 70 #kgs


#Defining the parameters that will be used in the rate equations.
#For these values T cells are defined in cells/ul, the tumor is defined as a volume in ul, and time is defined in days

#Blood
  #T cell elimination rates
ke1 <- 0.0104
ke2 <- 0.0104
ke3 <- 0.0104
ke4 <- 0.518

  #Differentiation
k12_B <- 0.14 #Naive to cm
k23_B <- 0.191  #cm to em
k34_B <- 0.355  #em to effector

  #Homeostatic proliferation
kp1_B <- 0.0005
kp2_B <- 0.007
kp3_B <- 0.007

  #Tumor cell proliferation
k5_B <- 0
#-------------------------------------------------------------------------
                            #Units
K0_B <- 3.85*10^6           #ul
K0_B <- (K0_B * tumor_density)

  #Tumor elimination
Vmax51_B <- 2570 #Naive     #uL(tumor_cell)·day^-1·(T-cells·μL(compartment)^-1)^-1
Vmax51_B <- Vmax51_B * (tumor_density/blood_vol)

Vmax52_B <- 4040 #cm        #uL(tumor_cell)·day^-1·(T-cells·μL(compartment)^-1)^-1
Vmax52_B <- Vmax52_B * (tumor_density/blood_vol)

Vmax53_B <- 3780 #em        #uL(tumor_cell)·day^-1·(T-cells·μL(compartment)^-1)^-1
Vmax53_B <- Vmax53_B * (tumor_density/blood_vol)

Vmax54_B <- 4240 #effector  #uL(tumor_cell)·day^-1·(T-cells·μL(compartment)^-1)^-1
Vmax54_B <- Vmax54_B * (tumor_density/blood_vol)


#Tumor cell volume at half maximum killing rate
KM5_B <- 276*10^3           #μL
KM5_B <- KM5_B *tumor_density

  #Tumor-induced T cell proliferation
Vmax1_B <- 8.46*10^-6       #(cells·μL(compartment)^-1) ·day^-1 ·uL(tumor)^-1  
Vmax1_B <- Vmax1_B * (blood_vol/tumor_density)

KM1_B <- 1.13               #T-cells * ul^-1
KM1_B <- KM1_B * blood_vol
#^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

#Migration
kT_BP <- 0.11    #T-cell blood to peripheral tissue
kT_PB <- 0.176   #T-cell peripheral tissue to blood
kCD_BP <- 0     #Tumor blood to peripheral tissue (0)
kCD_PB <- 0.001 #Tumor peripheral tissue to blood

#Peripheral tissue
  #Differentiation
k12_P <- 0.14 #Naive to cm
k23_P <- 0.191  #cm to em
k34_P <- 0.355  #em to effector

  #Homeostatic proliferation
kp1_P <- 0.0005
kp2_P <- 0.007
kp3_P <- 0.007
    #Tumor cell proliferation
k5_P <- 0.0023
#--------------------------------------------------------------------------------
                            #units:

K0_P <- 1.15*10^6           #ul
K0_P <- K0_P * tumor_density

#Tumor elimination
Vmax51_P <- 2570 #Naive     #uL(tumor_cell)·day^-1·(T-cells·μL(compartment)^-1)^-1
Vmax51_P <- Vmax51_P *(tumor_density/p_tissue_vol)

Vmax52_P <- 4040 #cm        #uL(tumor_cell)·day^-1·(T-cells·μL(compartment)^-1)^-1
Vmax52_P <- Vmax52_P *(tumor_density/p_tissue_vol)

Vmax53_P <- 3780 #em        #uL(tumor_cell)·day^-1·(T-cells·μL(compartment)^-1)^-1
Vmax53_P <- Vmax53_P *(tumor_density/p_tissue_vol)

Vmax54_P <- 4240 #effector  #uL(tumor_cell)·day^-1·(T-cells·μL(compartment)^-1)^-1
Vmax54_P <- Vmax54_P * (tumor_density/p_tissue_vol)

#Tumor cell volume at half maximum killing rate
KM5_P <- 276*10^3           #μL
KM5_P <- KM5_P * tumor_density

#Tumor-induced T cell proliferation
Vmax1_P <- 8.46*10^-6       #(cells·μL(compartment)^-1) ·day^-1 ·uL(tumor)^-1  
Vmax1_P <- Vmax1_P *(p_tissue_vol/tumor_density)

KM1_P <- 1.13               #T-cells * ul^-1
KM1_P <- KM1_P *p_tissue_vol

#^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

#Initial values
#species <- c("TnB","TcmB","TemB","TeffB","CD19+B","TnP","TcmP","TemP","TeffP","CD19+P")
state_0  <- c((2*10^6)*patient_weight,0,0,0,0,0,0,0,0,(50*10^3)*tumor_density)
N <- 10^6 #The number of reactions that will occur (jumps)
t <- rep(0,N) #Time at each jump
dt <- rep(0,N) #Size of each time step

states <- matrix(0,nrow=N,ncol=10) #The state matrix
colnames(states) <- species

states[1,] <- state_0

for(i in 2:N){                 #Loop is as long as the number of jumps given
  #Defining species numbers at the beginning of each time step
  TnB <- states[i-1,1]
  TcmB <- states[i-1,2]
  TemB <- states[i-1,3]
  TeffB <- states[i-1,4]
  CD19B <- states[i-1,5]
  TnP <- states[i-1,6]
  TcmP <- states[i-1,7]
  TemP <- states[i-1,8]
  TeffP <- states[i-1,9]
  CD19P <- states[i-1,10]
  
  #Defining a rate for each reaction
  
  #Elimination
  re_A <- ke1*TnB
  re_B <- ke2*TcmB
  re_C <- ke3*TemB
  re_D <- ke4*TeffB
  
  elim_rates <- c(re_A,re_B,re_C,re_D)
  
  #Differentiation in blood
  re_E <- k12_B*TnB
  re_F <- k23_B*TcmB
  re_G <- k34_B*TemB
  
  diff_B_rates <- c(re_E,re_F,re_G)
  
  #Proliferation in blood
  re_H <- kp1_B*TnB
  re_I <- kp2_B*TcmB
  re_J <- kp3_B*TemB
    #tumor
  re_K <- k5_B*(1-CD19B/K0_B)*CD19B
  
  prolif_B_rates <- c(re_H,re_I,re_J,re_K)
  
  #Tumor cell death in blood
  re_L <- (Vmax51_B*TnB*CD19B)/(KM5_B+CD19B)
  re_M <- (Vmax52_B*TcmB*CD19B)/(KM5_B+CD19B)
  re_N <- (Vmax53_B*TemB*CD19B)/(KM5_B+CD19B)
  re_O <- (Vmax54_B*TeffB*CD19B)/(KM5_B+CD19B)
  
  tumor_death_B_rates <- c(re_L,re_M,re_N,re_O)
  
  #Tumor cell induced proliferation in blood
  re_P <- (Vmax1_B*CD19B*TnB)/(KM1_B+TnB)
  re_Q <- (Vmax1_B*CD19B*TcmB)/(KM1_B+TnB)
  re_R <- (Vmax1_B*CD19B*TemB)/(KM1_B+TnB)
  
  induce_prolif_B_rates <- c(re_P,re_Q,re_R)
  
  #Migration
  re_S <- kT_BP*TnB
  re_T <- kT_PB*TnP
  re_U <- kT_BP*TcmB
  re_V <- kT_PB*TcmP
  re_W <- kT_BP*TemB
  re_X <- kT_PB*TemP
  re_Y <- kT_BP*TeffB
  re_Z <- kT_PB*TeffP
  re_A2 <- kCD_PB*CD19P
  
  migration_rates <- c(re_S,re_T,re_U,re_V,re_W,re_X,re_Y,re_Z,re_A2)
  
  #Differentiation in peripheral tissue
  re_B2 <- k12_P*TnP
  re_C2 <- k23_P*TcmP
  re_D2 <- k34_P*TemP
  
  diff_P_rates <- c(re_B2,re_C2,re_D2)
  
  #Natural proliferation in peripheral tissue
  re_E2 <- kp1_P*TnP
  re_F2 <- kp2_P*TcmP
  re_G2 <- kp3_P*TemP
    #tumor
  re_H2 <- k5_P*(1-CD19P/K0_P)*CD19P
  
  prolif_P_rates <- c(re_E2,re_F2,re_G2,re_H2)
  
  #Tumor cell death in peripheral tissue
  re_I2 <- (Vmax51_P*TnP*CD19P)/(KM5_P+CD19P)
  re_J2 <- (Vmax52_P*TcmP*CD19P)/(KM5_P+CD19P)
  re_K2 <- (Vmax53_P*TemP*CD19P)/(KM5_P+CD19P)
  re_L2 <- (Vmax54_P*TeffP*CD19P)/(KM5_P+CD19P)
  
  tumor_death_P_rates <- c(re_I2,re_J2,re_K2,re_L2)
  
  #Tumor cell induced proliferation in peripheral tissue
  re_M2 <- (Vmax1_P*CD19P*TnP)/(KM1_P+TnP)
  re_N2 <- (Vmax1_P*CD19P*TcmP)/(KM1_P+TnP)
  re_O2 <- (Vmax1_P*CD19P*TemP)/(KM1_P+TnP)
  
  induce_prolif_P_rates <- c(re_M2,re_N2,re_O2)
  
  #Combining all of these rates
  blood_rates <- c(elim_rates,diff_B_rates,prolif_B_rates,tumor_death_B_rates,induce_prolif_B_rates)
  periph_rates <- c(diff_P_rates,prolif_P_rates,tumor_death_P_rates,induce_prolif_P_rates)
  
  all_rates <- c(blood_rates,migration_rates,periph_rates)
  tot_rate <- sum(all_rates)
  
  dt[i-1]<- (-log(runif(1)))/tot_rate
  t[i]<-t[i-1]+dt[i-1] #Update time
  
  choice<-sample.int(length(all_rates),1,prob=all_rates) #Choose one of the reactions at random to occur, proportional to the rate of each reaction
  states[i,]<-states[i-1,]+s_mat[choice,] #Update state matrix
} #End of loop for each jump


#Generating graphs

plotting <- cbind(states,t) %>% data.frame() %>% mutate(t=t)

sub <- seq(1,N,length.out=N/100)
sub_plotting<- plotting[sub,]

print(ggplot(data=sub_plotting, aes(x=t,y=TnB)) + geom_line())
print(ggplot(data=sub_plotting, aes(x=t,y=TcmB)) + geom_line())
print(ggplot(data=sub_plotting, aes(x=t,y=TemB)) + geom_line())
print(ggplot(data=sub_plotting, aes(x=t,y=TeffB)) + geom_line())
print(ggplot(data=sub_plotting, aes(x=t,y=CD19.B)) + geom_line())
print(ggplot(data=sub_plotting, aes(x=t,y=TnP)) + geom_line())
print(ggplot(data=sub_plotting, aes(x=t,y=TcmP)) + geom_line())
print(ggplot(data=sub_plotting, aes(x=t,y=TemP)) + geom_line())
print(ggplot(data=sub_plotting, aes(x=t,y=TeffP)) + geom_line())
print(ggplot(data=sub_plotting, aes(x=t,y=CD19.P)) + geom_line())
