#Author: Sheri Tubach
#Date: 9/22/2025
#Manuscript: Informing Hepatitis C Intervention Implementation: Cost-Effectiveness Modeling for Rural Populations of Persons Who Inject Drugs
#Description: HCV compartment model and cost-effectiveness analysis using ode in DeSolve.
#Probabilistic sensitivity analysis with replacement, 95% confidence Ellipse calculating using ConfidenceEllipse package, and 
#cost-effectiveness acceptability curve.

library("deSolve")
library("tidyverse")
library("ConfidenceEllipse")

#----MODEL with MAT at 25% and 50%, SSP at 25% and 50% coverage with Increased Testing and Treatment

SIS_HCV_model <- function(time, state, parameters, ...) 
{with(as.list(c(state, parameters)), {
 
  N  <- S + I + C
  
  #FOI Calculations to account for the spontaneous recovery (rec.rate) and for the effect of MAT and/or SSP on reducing the FOI
    
  FOI <- FOI.rate * (1-rec.rate) * I/(S+I)                
  MAT.coverage <- MAT.rate / (MAT.rate + MATleave.rate)   
  SSP.coverage <- SSP.rate / (SSP.rate + SSPleave.rate)   
  effective.FOI <- FOI * (1 - MAT.coverage * MAT.eff) * (1 - SSP.coverage * SSP.eff)
  
  #RNA POC test modeled. This account for RNA tests in S population that are negative for cost calculations
  RNA_tests_S <- test.rate * S
  RNA_tests_I <- test.rate * I
  RNA_tests_C <- test.rate * C
  
  Total_RNA_test <- RNA_tests_S + RNA_tests_I + RNA_tests_C
  
  treated <- (RNA_tests_I * treat.rate) + (RNA_tests_C * treat.rate) 
  
  #DIFFERENTIAL EQUATIONS
  
  # Susceptible Compartment (S) 
  dS <- ent.rate - 
    effective.FOI * S +                               #FOI from S to I  
    (RNA_tests_I * treat.rate * SVR.rate) -            #Cessation for MAT users that leave the model
    (Mu1.rate + Mu2.rate + Mu3.rate) * S              #Exit from S compartment for cessation, mortality
  
  # Infectious Compartment (I)
  dI <- effective.FOI * S -                           #FOI from S to I
    (RNA_tests_I * treat.rate * SVR.rate) -         #Treatment and disease clearance still PWID move from I to S
    (dc.rate * I) -                                   #Develop complications move to C compartment
    (hcc.rate * I)-
    (Mu1.rate * I) -                                  #Infected PWID cessation move to FP - former PWID
    ((Mu2.rate + Mu3.rate) * I)                       # Exit from I compartment, mortality
  
  # Complications Compartment
  dC <- (dc.rate * I) + (hcc.rate * I) - 
    (dcm.rate * C) - (hccm.rate * C)
  
  #Health Outcomes, DALYs, and Cost Calculations
  mortality_C <- dcm.rate * C + hccm.rate * C        
  prevalence <- (I + C)                              
  prevalence_rate <- prevalence / N                  
  incidence <- effective.FOI * S                     
  DC <- dc.rate * (I)                                 
  HCC <- hcc.rate * (I)                              
  
  #DALY calculation taking into account different disability weights from HCV, DC, and HCC
  YLD_HCV <- prevalence * HDV_HCV                    
  YLD_DC <- (dc.rate * I) * HDV_DC                   
  YLD_HCC <- (hcc.rate * I) * HDV_HCC                
  YLL <- mortality_C * life_loss                              
  DALY <- YLD_HCV + YLD_DC + YLD_HCC + YLL           
  
  Discount <- 1/(1+ 0.03)^time                                                 
  Test_Cost <- Discount * (Total_RNA_test * Test_cost_RNA + RNA_Maint)          
  Treat_Cost <- Discount * (treated * treat_Cost_Param)    
  SSP_Cost <- Discount * ((S + I) * SSP.coverage * SSP_Cost_Param) 
  SSP_Imp <- ifelse(time == 1 & SSP.rate > 0, SSP_Implementation_Cost, 0)   #SSP implementation cost at time = 0
  MAT_Cost <- Discount * (MAT.coverage * (S + I) * MAT_Cost_Param)                        
  DC_Cost <- Discount * (DC_Cost_p * dc.rate * I)
  HCC_Cost <- Discount * (HCC_Cost_p * hcc.rate * I)
  
  
  # Return the derivatives
  list(c(dS, dI, dC), prevalence = prevalence, prevalence_rate = prevalence_rate, treated = treated, 
       incidence = incidence, mortality_C = mortality_C, Test_Cost = Test_Cost, Treat_Cost = Treat_Cost, DC = DC, HCC = HCC, 
       MAT_Cost = MAT_Cost, SSP_Cost = SSP_Cost, SSP_Imp = SSP_Imp,
       DC_Cost = DC_Cost, HCC_Cost = HCC_Cost, DALY = DALY)
  
  
})
}

#Start of the probabilistic sensitivity analysis

set.seed(123)
n_sim <- 5000

   #This is a progress bar to monitor progress
pb <- txtProgressBar(min = 0, max = nrow(parameter_combinations) * n_sim, style = 3)
counter <- 0

PSA_results <- list()

# Labels and values for intervention levels

test.rate.labels <- c("44%", "90%")
treat.rate.labels <- c("27%", "80%")
SSP.rate.labels <- c("0","25%", "50%")
MAT.rate.labels <- c("0","25%", "50%")

test.rate.values <- c(0.44, 0.90)
treat.rate.values <- c(0.27, 0.80)
SSP.rate.values <- c(0, 0.25, 0.50) 
MAT.rate.values <- c(0, 0.25, 0.50)

# Generate all combinations of interventions

parameter_combinations <- expand.grid(test.rate = test.rate.values, treat.rate = treat.rate.values,
                                      MAT.rate = MAT.rate.values, SSP.rate = SSP.rate.values
)

initial_state <- c(S = 11763, I = 14368, C = 9)

time <- seq(0, 20, by = 1) 

param_log_list <- list()

for (i in 1:nrow(parameter_combinations)) {
  
  # Extract parameters for this iteration
  
  test_rate <- parameter_combinations$test.rate[i]
  treat_rate <- parameter_combinations$treat.rate[i]
  MAT_rate <- parameter_combinations$MAT.rate[i]
  SSP_rate <- parameter_combinations$SSP.rate[i]
  
  for (j in 1:n_sim) {
   

# Define the model parameters
    
parameters <- list(
  FOI.rate = runif(1, 0.32, 0.73),                                         
  MAT.eff = rlnorm(1, meanlog = -0.7525, sdlog = 0.123),            
  SSP.eff = rlnorm(1, meanlog = -0.9419, sdlog = 0.3406),           
  SSPleave.rate = runif(1, 0.23, 0.63),                             
  MATleave.rate = runif(1, 0.06, 0.80),                             
  MAT.cess = runif(1, 0.20, 0.60),                                  
  rec.rate = runif(1, 0.22, 0.30),                                  
  ent.rate = runif (1, 1040, 1820), 
  RNA_Prev = runif (1,0.407, 0.467),
  SVR.rate = runif(1, 0.91, 0.95),                                  
  Mu1.rate = runif(1, 0.04, 0.14),                                   
  Mu2.rate = runif(1, 0.0125, 0.04025),                             
  Mu3.rate = runif(1, 0.0079, 0.0081),                                
  dc.rate = runif(1, 0.00105,	0.00513),                             
  hcc.rate = runif(1, 0.0001,	0.003),                               
  dcm.rate = runif(1, 0.11, 0.17),                                  
  hccm.rate = runif(1, 0.31, 0.79),                                 
  Test_cost_RNA = 92,                                               
  RNA_Maint = 2778,                                                 
  treat_Cost_Param = 26083,                                         
  MAT_Cost_Param = 3431,                                            
  SSP_Cost_Param = rgamma(1, shape = 52.5, scale = 18.5),           
  SSP_Implementation_Cost = rgamma(1, shape = 40.1, scale = 271.4), 
  DC_Cost_p = rgamma(1, shape = 11.52, scale = 2150),               
  HCC_Cost_p = rgamma(1, shape = 10.43, scale = 5000),              
  HDV_HCV = rbeta(1, 6.32, 118.2),                                  
  HDV_DC = rbeta(1, 27.94, 88.99),                                  
  HDV_HCC = rbeta(1, 19.29, 23.41),                                 
  life_loss = 28,                                                   
  test.rate = test_rate,                                            
  treat.rate = treat_rate,                                          
  MAT.rate = MAT_rate,                                              
  SSP.rate = SSP_rate)                                              

# Loop for the simulations
param_row <- parameters
param_row$sim_id <- j
param_log_list[[j]] <- as.data.frame(param_row)
  
  
  # Model uses ode solver from deSolve
  PSAoutput <- ode(y = initial_state, times = time, func = SIS_HCV_model, parms = parameters)
  
  # Add Identifying Info about simulation number and test, treat, MAT, and SSP coverage
  df_output <- as.data.frame(PSAoutput)
  df_output$sim_id <- j
  df_output$test_rate <- test_rate
  df_output$treat_rate <- treat_rate
  df_output$MAT_rate <- MAT_rate
  df_output$SSP_rate <- SSP_rate
  
  # Map values to their descriptive labels
  df_output$`Test_Rate` <- test.rate.labels[which(test.rate.values == parameter_combinations$test.rate[i])]
  df_output$`Treat_Rate` <- treat.rate.labels[which(treat.rate.values == parameter_combinations$treat.rate[i])]
  df_output$`MAT_Rate` <- MAT.rate.labels[which(MAT.rate.values == parameter_combinations$MAT.rate[i])]
  df_output$`SSP_Rate` <- SSP.rate.labels[which(SSP.rate.values == parameter_combinations$SSP.rate[i])]
  
  
  PSA_results[[length(PSA_results) + 1]] <- df_output
  
 
  param_log_list[[j]] <- as.data.frame(parameters)
  
  # Sets the counter and progress bar to monitor program
  counter <- counter + 1
  setTxtProgressBar(pb, counter)
  
  }
}

close(pb)


# Combine all simulation outputs into one data frame

PSA_results_df <- do.call(rbind, PSA_results)


# Calculating cost summary measures

PSA_summary <- PSA_results_df %>%
  group_by(Test_Rate, Treat_Rate, MAT_Rate, SSP_Rate, sim_id) %>%
  summarize(
    Total_Incidence = sum(incidence),
    Total_Death = sum(mortality_C),
    Total_Treated = sum(treated),
    Total_Comp = sum(DC + HCC),
    Total_DALY = sum(DALY),
    Total_Test_Cost = sum(Test_Cost),
    Total_Treat_Cost = sum(Treat_Cost),
    Total_Comp_Cost = sum(DC_Cost + HCC_Cost),
    Total_MAT_Cost = sum(MAT_Cost),
    Total_SSP_Cost = sum(SSP_Cost),
    Total_SSPI_Cost = sum(SSP_Imp),
    Total_SSPIMP_Cost = Total_SSP_Cost + Total_SSPI_Cost,
    Total_Cost = (Total_Test_Cost + Total_Treat_Cost + Total_Comp_Cost 
                     + Total_SSPIMP_Cost + Total_MAT_Cost), 
    .groups = "drop")

# Extracting DALYs and Cost from basecase to calculate incremental cost and DALYs averted

reference_summary <- PSA_summary %>%
  filter(Test_Rate == "44%", Treat_Rate == "27%", MAT_Rate == "0", SSP_Rate == "0") %>%
  select(sim_id, reference_DALY = Total_DALY,
         reference_cost = Total_Cost)

reference_DALY <- reference_summary$reference_DALY
reference_cost <- reference_summary$reference_cost

HCV_PSA <- PSA_summary %>%
  mutate(
    Incremental_Cost = Total_Cost - reference_cost,
    DALY_Averted = reference_DALY - Total_DALY,
    ICER_DALY = Incremental_Cost / DALY_Averted
    
    )


#Calculate average yearly incremental cost and DALYs averted

HCV_PSA$Incr_Cost_Yr <- HCV_PSA$Incremental_Cost / 20
HCV_PSA$DALY_Yr <- HCV_PSA$DALY_Averted / 20


# Define the nine scenarios for further analysis

scenarios <- data.frame(
  Test_Rate  = rep("90%", 9 ),
  Treat_Rate = rep("80%", 9),
  SSP_Rate   = c("0", "25%", "50%","0", "0", "25%", "50%", "25%", "50%"), 
  MAT_Rate   = c("0", "0", "0", "25%", "50%", "25%", "25%", "50%", "50%") 
)

combined_results <- data.frame()

for (i in seq_len(nrow(scenarios))) {
  s <- scenarios[i, ]
  
  # Force each to be a character string, not a data.frame column
  test_val  <- as.character(s$Test_Rate)
  treat_val <- as.character(s$Treat_Rate)
  ssp_val   <- as.character(s$SSP_Rate)
  mat_val   <- as.character(s$MAT_Rate)
  
  # Scenario label
  label <- paste0("T", s$Test_Rate, 
                  "R", s$Treat_Rate, 
                  "S", gsub("%", "", s$SSP_Rate), 
                  "M", s$MAT_Rate)
  
  
  filtered_data <- HCV_PSA %>%
    filter(Test_Rate  == test_val,
           Treat_Rate == treat_val,
           SSP_Rate   == ssp_val,
           MAT_Rate   == mat_val)
  
  
  cat("Scenario:", label, "\n")
  cat("Rows:", nrow(filtered_data), "\n")
  
  # Calculates the 95% Confidence ellipse
  ellipse <- confidence_ellipse(filtered_data, x = DALY_Yr, y = Incr_Cost_Yr, conf_level = 0.95) 
  
  # Calculates the CEAC
  wtp <- seq(0, 100000, by = 1000)
  ceac_prob <- sapply(wtp, function(w) {
    nmb <- filtered_data$DALY_Yr * w - filtered_data$Incr_Cost_Yr
    mean(nmb > 0)
  })
  
  # This section prepares labeled data frames for cost-effectiveness analysis (CEA) outputs:
  # ceac: stores willingness-to-pay (WTP) thresholds and probability of cost-effectiveness for CEAC plotting.
  # filtered_data_labeled: adds scenario labels and formatting for incremental costs and DALYs averted.
  # ellipse_labeled: formats 95% confidence ellipse data with scenario and type labels.
  # ceac_labeled: creates a labeled version of CEAC 
  # Type and scenario allows for filtering for plotting in Excel.
  
   ceac <- data.frame(
    DALY_Yr = NA,
    Incr_Cost_Yr = NA,
    Type = "CEAC",
    Scenario = label,
    WTP = wtp,
    Probability_CE = ceac_prob
  )
  
  filtered_data_labeled <- filtered_data %>%
    mutate(
      WTP = NA,
      Probability_CE = NA,
      Type = "DALY_Incremental_Cost",
      Scenario = label
    )
  
  ellipse_labeled <- ellipse %>%
    mutate(
      Type = "Ellipse",
      Scenario = label,
      WTP = NA,
      Probability_CE = NA
    )
  
  ceac_labeled <- ceac %>%
    mutate(
      DALY_Yr = NA,
      Incr_Cost_Yr = NA,
      Type = "CEAC",
      Scenario = label
    )
  
  scenario_combined <- bind_rows(filtered_data_labeled, ellipse_labeled, ceac_labeled)
  
  combined_results <- bind_rows(combined_results, scenario_combined)
  
}

write.csv(combined_results,
          file = "C:/Users/stubach/OneDrive - State of Kansas, OITS/Documents/Dissertation/DATA/New Analysis/Final Analysis/PSA_Ellipse_CEAC_NEW.csv",
          row.names = FALSE)

