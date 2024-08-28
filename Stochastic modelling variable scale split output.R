library(tidyverse)
library(Hmisc)

species <- c("TnB","TcmB","TemB","TeffB","CD19+B","TnP","TcmP","TemP","TeffP","CD19+P")

#defining the scales altering the rates and number of cells involved in each reaction

scales <- c(1,#A
            1,#B
            1,#C
            1,#D
            1,#E
            1,#F
            1,#G
            1,#H
            1,#I
            1,#J
            100,#K
            100,#L
            100,#M
            100,#N
            100,#O
            100,#P
            100,#Q
            100,#R
            1,#S
            1,#T
            1,#U
            1,#V
            1,#Q
            1,#R
            1,#S
            1,#T
            1,#U
            1,#V
            1,#W
            1,#X
            1,#Y
            1,#Z
            100,#A2
            1,#B2
            1,#C2
            1,#D2
            1,#E2
            1,#F2
            1,#G2
            100,#H2
            100,#I2
            100,#J2
            100,#K2
            100,#L2
            100,#M2
            100,#N2
            100 #O2
)

scales <- 10*scales

#Defining the stoicheometric matrix

#The reactions representing the natural elimination rate of CAR-T cells in the blood 
react_elim <- c(-scales[1],rep(0,9), #A
                0,-scales[2],rep(0,8),#B
                0,0,-scales[3],rep(0,7),#C
                rep(0,3),-scales[4],rep(0,6)) %>% #D
  matrix(nrow=10) %>% t()

colnames(react_elim) <- species

#The reactions representing the differentiation of CAR-T cells in the blood
react_diff_B <- c(-scales[5],scales[5],rep(0,8),#E
                  0,-scales[6],scales[6],rep(0,7),#F
                  0,0,-scales[7],scales[7],rep(0,6)) %>% #G
  matrix(nrow = 10) %>% t()

colnames(react_diff_B) <- species

#The reactions representing the natural proliferation rate of species in the blood
react_nat_prolif_B <- c(scales[8],rep(0,9),#H
                    0,scales[9],rep(0,8),#I
                    0,0,scales[10],rep(0,7),#J
                    rep(0,4),scales[11],rep(0,5)) %>% #K
  matrix(nrow = 10) %>% t()

colnames(react_nat_prolif_B) <- species

#The reactions representing tumor cells in the blood being killed by CAR-T cells
react_kill_B <- c(rep(0,4),-scales[12],rep(0,5), #L
                  rep(0,4),-scales[13],rep(0,5), #M
                  rep(0,4),-scales[14],rep(0,5),#N
                  rep(0,4),-scales[15],rep(0,5)) %>% #O
  matrix(nrow=10) %>% t()

colnames(react_kill_B)<- species

#The reactions representing the proliferation of CAR-T cells in the blood induced by the presence of tumor cells
react_induce_prolif_B <- c(scales[16],rep(0,9), #P
                           0,scales[17],rep(0,8), #Q
                           0,0,scales[18],rep(0,7)) %>% #R
  matrix(nrow=10) %>% t()

#The reactions representing the migration of species between compartments
react_migration <- c(-scales[19],rep(0,4),scales[19],rep(0,4), #S
                      scales[20],rep(0,4),-scales[20],rep(0,4), #T
                     0,-scales[21],rep(0,4),scales[21],rep(0,3), #U
                     0,scales[22],rep(0,4),-scales[22],rep(0,3), #V
                     0,0,-scales[23],rep(0,4),scales[23],0,0,#W
                     0,0,scales[24],rep(0,4),-scales[24],0,0,#X
                     rep(0,3),-scales[25],rep(0,4),scales[25],0,#Y
                     rep(0,3),scales[26],rep(0,4),-scales[26],0, #Z
                     rep(0,4),scales[27],rep(0,4),-scales[27]) %>% #A2
                    #Tumor cells are assumed to be unable to migrate back to peripheral tissue
  matrix(nrow = 10) %>% t()

colnames(react_migration) <- species

#The reactions representing CAR-T cell differentiation in the peripheral tissue
react_diff_P <- c(rep(0,5),-scales[28],scales[28],rep(0,3), #B2
                  rep(0,6),-scales[29],scales[29],0,0, #C2
                  rep(0,7),-scales[30],scales[30],0) %>% #D2
  matrix(nrow=10) %>% t()

colnames(react_diff_P) <- species

#The reactions representing natural species proliferation in peripheral tissue
react_nat_prolif_P <- c(rep(0,5),scales[31],rep(0,4), #E2
                        rep(0,6),scales[32],rep(0,3),#F2
                        rep(0,7),scales[33],0,0,#G2
                        rep(0,9),scales[34]) %>% #H2
  matrix(nrow = 10) %>% t()

colnames(react_nat_prolif_P) <- species

#The reactions representing tumor cells in peripheral tissue being killed by CAR-T cells
react_kill_P <- c(rep(0,9),-scales[35], #I2
                  rep(0,9),-scales[36], #J2
                  rep(0,9),-scales[37], #K2
                  rep(0,9),-scales[38]) %>% #L2 
  matrix(nrow=10) %>% t()

#The reactions representing the proliferation of CAR-T cells in peripheral tissue induced by the presence of tumor cells
react_induce_prolif_P <- c(rep(0,5),scales[39],rep(0,4), #M2
                           rep(0,6),scales[40],rep(0,3), #N2
                           rep(0,7),scales[41],0,0) %>%  #O2
  matrix(nrow = 10) %>% t()

#All the reactions taking place in the blood
blood_reactions <- rbind(react_elim,react_diff_B,react_nat_prolif_B,react_kill_B,react_induce_prolif_B)

#All the reactions taking place in the peripheral tissue
p_tissue_reactions <- rbind(react_diff_P,react_nat_prolif_P,react_kill_P,react_induce_prolif_P)

#The entire stoicheometric matrix
s_mat <- rbind(blood_reactions,react_migration,p_tissue_reactions)

#Assumptions
blood_vol <- 5*10^6
p_tissue_vol <- 1.75*10^6
tumor_density <- (1/2)*10^7 #Assuming a tumor cell volume of 200 um^3 and a maximally dense tumor
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
#state0_base <- 10^-2*c(5.46,5.24,3.25,40.6,0,6.66,42.6,43.2,132,(2.57*10^6))
#state0_adjust <- c(rep(blood_vol,4),tumor_density,rep(p_tissue_vol,4),tumor_density)
#state_0 <- state0_base*state0_adjust

state_0  <- c((2*10^6)*patient_weight,0,0,0,0,0,0,0,0,(1*10^3)*tumor_density)


N <- 10^5 #The number of reactions that will occur (jumps) between each data output
N2 <- 10^1 #The number of times that the data is output throughout a run

t <- rep(0,(N*N2+1)) #Time at each jump
dt <- rep(0,(N*N2+1)) #Size of each time step

states <- matrix(0,nrow=(N*N2+1),ncol=10) #The state matrix
colnames(states) <- species

states[1,] <- state_0

#for running for a particularly long length of time, outputting data as you go

save_curve <- function(path,iteration,title){
  fn=paste0(path,"/",iteration,"/",title,".png")
  ggsave(filename=fn,create.dir = TRUE)
}

path <- "Reactions6Naive2-6Tumor1-3_split" #Change based on order of reactions, dose of CAR-T cells/kg and tumor size in ul. Create a directory prior to running

for(j in 1:N2){  
  for(i in 1:N){                 #Loop is as long as the number of jumps given
    #Defining species numbers at the beginning of each time step
    i<- i+1+((j-1)*N)
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
    re_A <- ke1*TnB/scales[1]
    re_B <- ke2*TcmB/scales[2]
    re_C <- ke3*TemB/scales[3]
    re_D <- ke4*TeffB/scales[4]
    
    elim_rates <- c(re_A,re_B,re_C,re_D)
    
    #Differentiation in blood
    re_E <- k12_B*TnB/scales[5]
    re_F <- k23_B*TcmB/scales[6]
    re_G <- k34_B*TemB/scales[7]
    
    diff_B_rates <- c(re_E,re_F,re_G)
    
    #Proliferation in blood
    re_H <- kp1_B*TnB/scales[8]
    re_I <- kp2_B*TcmB/scales[9]
    re_J <- kp3_B*TemB/scales[10]
      #tumor
    re_K <- (k5_B*(1-CD19B/K0_B)*CD19B)/scales[11]
    
    prolif_B_rates <- c(re_H,re_I,re_J,re_K)
    
    #Tumor cell death in blood
    re_L <- ((Vmax51_B*TnB*CD19B)/(KM5_B+CD19B))/scales[12]
    re_M <- ((Vmax52_B*TcmB*CD19B)/(KM5_B+CD19B))/scales[13]
    re_N <- ((Vmax53_B*TemB*CD19B)/(KM5_B+CD19B))/scales[14]
    re_O <- ((Vmax54_B*TeffB*CD19B)/(KM5_B+CD19B))/scales[15]
    
    tumor_death_B_rates <- c(re_L,re_M,re_N,re_O)
    
    #Tumor cell induced proliferation in blood
    re_P <- ((Vmax1_B*CD19B*TnB)/(KM1_B+TnB))/scales[16]
    re_Q <- ((Vmax1_B*CD19B*TcmB)/(KM1_B+TnB))/scales[17]
    re_R <- ((Vmax1_B*CD19B*TemB)/(KM1_B+TnB))/scales[18]
    
    induce_prolif_B_rates <- c(re_P,re_Q,re_R)
    
    #Migration
    re_S <- kT_BP*TnB/scales[19]
    re_T <- kT_PB*TnP/scales[20]
    re_U <- kT_BP*TcmB/scales[21]
    re_V <- kT_PB*TcmP/scales[22]
    re_W <- kT_BP*TemB/scales[23]
    re_X <- kT_PB*TemP/scales[24]
    re_Y <- kT_BP*TeffB/scales[25]
    re_Z <- kT_PB*TeffP/scales[26]
    re_A2 <- (kCD_PB*CD19P)/scales[27]
    
    migration_rates <- c(re_S,re_T,re_U,re_V,re_W,re_X,re_Y,re_Z,re_A2)
    
    #Differentiation in peripheral tissue
    re_B2 <- k12_P*TnP/scales[28]
    re_C2 <- k23_P*TcmP/scales[29]
    re_D2 <- k34_P*TemP/scales[30]
    
    diff_P_rates <- c(re_B2,re_C2,re_D2)
    
    #Natural proliferation in peripheral tissue
    re_E2 <- kp1_P*TnP/scales[31]
    re_F2 <- kp2_P*TcmP/scales[32]
    re_G2 <- kp3_P*TemP/scales[33]
      #tumor
    re_H2 <- (k5_P*(1-CD19P/K0_P)*CD19P)/scales[34]
    
    prolif_P_rates <- c(re_E2,re_F2,re_G2,re_H2)
    
    #Tumor cell death in peripheral tissue
    re_I2 <- ((Vmax51_P*TnP*CD19P)/(KM5_P+CD19P))/scales[35]
    re_J2 <- ((Vmax52_P*TcmP*CD19P)/(KM5_P+CD19P))/scales[36]
    re_K2 <- ((Vmax53_P*TemP*CD19P)/(KM5_P+CD19P))/scales[37]
    re_L2 <- ((Vmax54_P*TeffP*CD19P)/(KM5_P+CD19P))/scales[38]
    
    tumor_death_P_rates <- c(re_I2,re_J2,re_K2,re_L2)
    
    #Tumor cell induced proliferation in peripheral tissue
    re_M2 <- ((Vmax1_P*CD19P*TnP)/(KM1_P+TnP))/scales[39]
    re_N2 <- ((Vmax1_P*CD19P*TcmP)/(KM1_P+TnP))/scales[40]
    re_O2 <- ((Vmax1_P*CD19P*TemP)/(KM1_P+TnP))/scales[41]
    
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
  
  plotting <- cbind(states,t) %>% data.frame() #%>% mutate(t=t)
  
  sub <- seq(1,N*j,length.out=(N*j)/100)
  sub_plotting<- plotting[sub,]
  
  ggplot(data=sub_plotting, aes(x=t,y=TnB)) + geom_line()
  save_curve(path, j,"TnB")
  
  tcmb_plot<-ggplot(data=sub_plotting, aes(x=t,y=TcmB)) + geom_line()
  save_curve(path,j,"TcmB")
  
  temb_plot<-ggplot(data=sub_plotting, aes(x=t,y=TemB)) + geom_line()
  save_curve(path,j,"TemB")
  
  teffb_plot<-ggplot(data=sub_plotting, aes(x=t,y=TeffB)) + geom_line()
  save_curve(path,j,"TeffB")
  
  tumorb_plot<-ggplot(data=sub_plotting, aes(x=t,y=CD19.B)) + geom_line()
  save_curve(path,j,"CD19B")
  
  tnp_plot<-ggplot(data=sub_plotting, aes(x=t,y=TnP)) + geom_line()
  save_curve(path,j,"Tnp")
  
  tcmp_plot<-ggplot(data=sub_plotting, aes(x=t,y=TcmP)) + geom_line()
  save_curve(path,j,"Tcmp")
  
  temp_plot<-ggplot(data=sub_plotting, aes(x=t,y=TemP)) + geom_line()
  save_curve(path,j,"TemP")
  
  teffp_plot<- ggplot(data=sub_plotting, aes(x=t,y=TeffP)) + geom_line()
  save_curve(path,j,"TeffP")
  
  tumorp_plot<-ggplot(data=sub_plotting, aes(x=t,y=CD19.P)) + geom_line()
  save_curve(path,j,"CD19P")
  
  csvfn <- paste0(path,"/",j,"/data.csv")
  write.csv(sub_plotting,csvfn)
}

#print(ggplot(data=sub_plotting, aes(x=t,y=TnB)) + geom_line())
#print(ggplot(data=sub_plotting, aes(x=t,y=TcmB)) + geom_line())
#print(ggplot(data=sub_plotting, aes(x=t,y=TemB)) + geom_line())
#print(ggplot(data=sub_plotting, aes(x=t,y=TeffB)) + geom_line())
#print(ggplot(data=sub_plotting, aes(x=t,y=CD19.B)) + geom_line())
#print(ggplot(data=sub_plotting, aes(x=t,y=TnP)) + geom_line())
#print(ggplot(data=sub_plotting, aes(x=t,y=TcmP)) + geom_line())
#print(ggplot(data=sub_plotting, aes(x=t,y=TemP)) + geom_line())
#print(ggplot(data=sub_plotting, aes(x=t,y=TeffP)) + geom_line())
#print(ggplot(data=sub_plotting, aes(x=t,y=CD19.P)) + geom_line())
