
# Author: Sheri Tubach
# Date: 2/28/2026
# Manuscript: Informing Hepatitis C Intervention Implementation: Cost-Effectiveness Modeling for Rural Populations of Persons Who Inject Drugs
# Description: HCV compartment model and cost-effectiveness analysis using ode in DeSolve. 


library("deSolve")
library("dplyr")
library("tidyverse")

#----MODEL with MAT at 25% and 50%, SSP at 25% and 50% coverage with Increased Testing and Treatment----

SIS_HCV_model <- function(time, state, parameters, ...) 
{with(as.list(c(state, parameters)), {
 
  N  <- S + I + R + DC + HCC
  
  #FOI Calculations to account for the spontaneous recovery (rec.rate) and for the effect of MAT and/or SSP on reducing the FOI
  
  MAT.coverage <- MAT.rate / (MAT.rate + MATleave.rate)
  SSP.coverage <- SSP.rate / (SSP.rate + SSPleave.rate)
  
  FOI <- beta * shared * ((I + R + DC + HCC)/N) #CHANGE TO INCLUDE R
  effective.FOI <- FOI * (1 - MAT.coverage * MAT.eff) * (1 - SSP.coverage * SSP.eff)
  
  #RNA POC test modeled. This account for RNA tests in S population that are negative for cost calculations
  RNA_tests_S <- test.rate * S
  RNA_tests_I <- test.rate * I
  RNA_tests_DC <- test.rate * DC
  RNA_tests_HCC <- test.rate * HCC
  
  Total_RNA_test <- RNA_tests_S + RNA_tests_I + RNA_tests_DC + RNA_tests_HCC
  treated <- treat.rate * (RNA_tests_I + RNA_tests_DC + RNA_tests_HCC + R)      
  
  
  #DIFFERENTIAL EQUATIONS
  
  # Susceptible Compartment (S)
  dS <- ent.rate - effective.FOI * S * (1 - rec.rate) + (test.rate * treat.rate * SVR.rate * I) +
    (R * treat.rate * SVR.rate)- 
    (Mu1.rate + Mu2.rate + Mu3.rate) * S              
  
  # Infectious Compartment (I)
  dI <- effective.FOI * S * (1 - rec.rate) - (test.rate * treat.rate * SVR.rate * I) - 
    (test.rate * treat.rate * (1 - SVR.rate) * I) -
    (dc.rate * I) - (hcc.rate * I) - (Mu1.rate + Mu2.rate + Mu3.rate) * I 
  
  dR <- (test.rate * treat.rate * (1 - SVR.rate) * I) - (R * treat.rate * SVR.rate) - 
    (dc.rate * R) - (hcc.rate * R) - (Mu1.rate + Mu2.rate + Mu3.rate) * R 
  
  # Complications Compartments
  dDC <- (dc.rate * I) + (dc.rate * R) - (dcm.rate * DC) -(dclt.rate * DC) 
  
  dHCC <- (hcc.rate * I) + (hcc.rate * R)- (hccm.rate * HCC) - (hcclt.rate * HCC) 
  
  
    # Health Outcomes, DALYs, and Cost Calculations
  mortality_C <- dcm.rate * DC + hccm.rate * HCC 
  prevalence <- (I + R + DC + HCC)
  prevalence_rate <- prevalence / N                        
  incidence <- effective.FOI * S   
  
  Discount <- 1/(1+ 0.03)^time 
  
  YLD_HCV <- (I + R) * HDV_HCV
  YLD_DC <- (DC * HDV_DC)
  YLD_HCC <- (HCC * HDV_HCC)
  YLL <- mortality_C * life_loss
  DALY <- (YLD_HCV + YLD_DC + YLD_HCC + YLL)
  
  Test_Cost <- Discount * (Total_RNA_test * Test_cost_RNA + RNA_Maint)          
  Treat_Cost <- Discount * ((treat.rate * (RNA_tests_I + RNA_tests_DC + RNA_tests_HCC) * treat_Cost_Param) +
    ((R * treat.rate) * retreat_Cost_Param))
                              
  
  SSP_Cost <- Discount * (SSP.coverage * N * SSP_Cost_Param) 
  SSP_Imp <- ifelse(time < 1 & SSP.rate > 0, SSP_Implementation_Cost, 0) #SSP implementation cost at time = 0
  MAT_Cost <- Discount * (MAT.coverage * N * MAT_Cost_Param)                        
  DC_Cost <- Discount * (DC_Cost_p * DC)
  HCC_Cost <- Discount * (HCC_Cost_p * HCC)
  
  LT_Cost <- Discount * (dclt.rate * DC + hcclt.rate * HCC) * LT_Cost_p
  
  # Return the derivatives
  list(c(dS, dI, dR, dDC, dHCC), prevalence = prevalence, prevalence_rate = prevalence_rate, treated = treated, 
       incidence = incidence, mortality_C = mortality_C, Test_Cost = Test_Cost, Treat_Cost = Treat_Cost,  
       MAT_Cost = MAT_Cost, SSP_Cost = SSP_Cost, SSP_Imp = SSP_Imp,
      DC_Cost = DC_Cost, HCC_Cost = HCC_Cost, LT_Cost = LT_Cost, DALY = DALY)
  
  
})
}

initial_state <- c(S = 11763, I = 14368, R = 0, DC = 6, HCC = 3)

time <- seq(0, 20, by = 1) 

fixed_parameters <- list(
  beta = 0.01,
  shared = 60,    
  MAT.eff = 0.50,
  SSP.eff = 0.56,
  rec.rate = 0.25,       
  ent.rate = 1560,       
  RNA_Prev = 0.437,
  SSPleave.rate = 0.43,
  MATleave.rate = 0.39,
  SVR.rate = 0.93,       
  Mu1.rate = 0.07,       
  Mu2.rate = 0.003,     
  Mu3.rate = 0.0264,      
  dc.rate = 0.00336,       
  hcc.rate = 0.00174,      
  dcm.rate = 0.14,       
  hccm.rate = 0.55,
  dclt.rate = 0.023,
  hcclt.rate = 0.040,
  Test_cost_RNA = 92,    
  RNA_Maint = 2778,      
  treat_Cost_Param = 26083,
  retreat_Cost_Param = 72579,
  MAT_Cost_Param = 3431,   
  SSP_Cost_Param = 983,     
  SSP_Implementation_Cost = 12077, 
  DC_Cost_p = 24790,
  HCC_Cost_p = 52183,
  LT_Cost_p = 1017800,
  HDV_PWID = 0.088,        
  HDV_HCV = 0.05,         
  HDV_DC = 0.236,     
  HDV_HCC = 0.451,
  life_loss = 28)

# Labels and values for intervention levels

test.rate.labels <- c("44%", "90%")
treat.rate.labels <- c("27%", "80%")
SSP.rate.labels <- c("0","25%", "50%")
MAT.rate.labels <- c("0","25%", "50%")


test.rate.values <- c(0.44, 0.90)
treat.rate.values <- c(0.27, 0.80)
SSP.rate.values <- c(0, 0.25, 0.50) 
MAT.rate.values <- c(0, 0.25, 0.50)

#Combining all combination of intervention levels

parameter_combinations <- expand.grid(test.rate = test.rate.values, treat.rate = treat.rate.values,
                                      MAT.rate = MAT.rate.values, SSP.rate = SSP.rate.values
                                      )
MSresults_20 <- list()
for (i in 1:nrow(parameter_combinations)) {
 parameters <- c(fixed_parameters,
                  test.rate = parameter_combinations$test.rate[i],
                  treat.rate = parameter_combinations$treat.rate[i],
                  MAT.rate = parameter_combinations$MAT.rate[i],
                  SSP.rate = parameter_combinations$SSP.rate[i])
  
  
  # Model uses ode solver from deSolve
  MSoutput_20 <- ode(y = initial_state, times = time, func = SIS_HCV_model, parms = parameters)
  
  MSoutput_20 <- as.data.frame(MSoutput_20)
  
  
  MSoutput_20$`Test_Rate` <- test.rate.labels[which(test.rate.values == parameter_combinations$test.rate[i])]
  MSoutput_20$`Treat_Rate` <- treat.rate.labels[which(treat.rate.values == parameter_combinations$treat.rate[i])]
  MSoutput_20$`MAT_Rate` <- MAT.rate.labels[which(MAT.rate.values == parameter_combinations$MAT.rate[i])]
  MSoutput_20$`SSP_Rate` <- SSP.rate.labels[which(SSP.rate.values == parameter_combinations$SSP.rate[i])]
  

  MSoutput_20$prevalence <- round(MSoutput_20$prevalence, 0)
  MSoutput_20$incidence <- round(MSoutput_20$incidence, 0)
  MSoutput_20$mortality_C <- round(MSoutput_20$mortality_C, 0)
  
  MSresults_20[[i]] <- MSoutput_20

}

MS_results_20 <- bind_rows(MSresults_20)

# Formatting results for export

MSHCV_20 <- MS_results_20 %>% 
  pivot_wider(
    id_cols = time,           # Use 'time' as the identifier variable
    names_from = c(MAT_Rate, SSP_Rate, Test_Rate, Treat_Rate),
    names_glue = "{.value}_MAT_{MAT_Rate}_SSP_{SSP_Rate}_Test_{Test_Rate}_Treat_{Treat_Rate}",
    values_from = c(S, I, R, DC, HCC, prevalence, prevalence_rate, incidence)
  )

All_HCV_Sum <- "C:\\Users\\stubach\\OneDrive - State of Kansas, OITS\\Documents\\Dissertation\\DATA\\R program\\HCV CEA\\Final Commented Scripts\\Revision\\HCV_Interventions_v7.csv"
write.csv(MSHCV_20, file = All_HCV_Sum, row.names = FALSE)

#Calculating cost summary measures

MS_summary_20 <- MS_results_20 %>%
  group_by(Test_Rate, Treat_Rate, MAT_Rate, SSP_Rate) %>%
  summarize(
    Total_Incidence = sum(incidence),
    Total_Death = sum(mortality_C),
    Total_Treated = sum(treated),
    Total_DALY = sum(DALY),
    Total_MAT_Cost = sum(MAT_Cost),
    Total_Test_Cost = sum(Test_Cost),
    Total_Treat_Cost = sum(Treat_Cost),
    Total_Comp_Cost = sum(DC_Cost + HCC_Cost),
    Total_LT_Cost = sum(LT_Cost),
    Total_SSP_Cost = sum(SSP_Cost),
    Total_SSPI_Cost = sum(SSP_Imp),
    Total_SSPIMP_Cost = Total_SSP_Cost + Total_SSPI_Cost,
    Total_Cost = (Total_Test_Cost + Total_Treat_Cost + Total_Comp_Cost + Total_LT_Cost 
                  + Total_SSPIMP_Cost + Total_MAT_Cost), # sum of all equations
      ) %>%
  ungroup()

#Extracting DALYs and Cost from basecase to calculate incremental cost and DALYs averted

reference_summary_20 <- MS_summary_20 %>%
  filter(Test_Rate == "44%", Treat_Rate == "27%", MAT_Rate == "0", SSP_Rate == "0") %>%
  select(reference_DALY = Total_DALY,
         reference_cost = Total_Cost,
         reference_inc = Total_Incidence)

reference_DALY <- reference_summary_20$reference_DALY[1]
reference_cost <- reference_summary_20$reference_cost[1]
reference_inc <- reference_summary_20$reference_inc[1]

HCV_Final_20 <- MS_summary_20 %>%
  mutate(
    Incremental_Cost = Total_Cost - reference_cost,
    IC_yr = (Incremental_Cost / 20),
    DALY_Averted = reference_DALY - Total_DALY,
    DA_yr = DALY_Averted / 20,
    ICER_DALY = Incremental_Cost / DALY_Averted,
    HCV_Averted = reference_inc - Total_Incidence
    )

HCV_cost_sum <- "C:\\Users\\stubach\\OneDrive - State of Kansas, OITS\\Documents\\Dissertation\\DATA\\R program\\HCV CEA\\Final Commented Scripts\\Revision\\HCV_cost_summary_v7.csv"

write.csv(HCV_Final_20, file = HCV_cost_sum, row.names = FALSE)









