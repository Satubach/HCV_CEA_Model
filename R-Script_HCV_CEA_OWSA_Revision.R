
# Author: Sheri Tubach
# Date: 9/22/2025
# Manuscript: Informing Hepatitis C Intervention Implementation: Cost-Effectiveness Modeling for Rural Populations of Persons Who Inject Drugs
# Description: HCV compartment model and cost-effectiveness analysis using ode in DeSolve. One-way sensitivity analysis

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

# Sets the initial parameters
base_initial <- c(S = 11763, I = 14368, R = 0, DC = 6, HCC = 3)
base_time <- seq(0, 20, 1)


base_parameters <- list(
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
  life_loss = 28,
  test.rate = 0.44,
  treat.rate = 0.27,
  MAT.rate = 0,
  SSP.rate = 0
  )


intervention_list <- list(
  "Test_Treat" = list(test.rate = 0.9, treat.rate = 0.8, SSP.rate = 0, MAT.rate = 0),
  "SSP"        = list(test.rate = 0.9, treat.rate = 0.8, SSP.rate = 0.25, MAT.rate = 0),
  "MAT"        = list(test.rate = 0.9, treat.rate = 0.8, SSP.rate = 0, MAT.rate = 0.25),
  "MAT_SSP"    = list(test.rate = 0.9, treat.rate = 0.8, SSP.rate = 0.25, MAT.rate = 0.25)
)


# Define variables and ranges for the sensitivity analysis
owsa_list <- list(
  
    "HCV/PWID prevalence (initial_state)" = list (
      values = list (
        low = c(S = 15500, I = 10632, R = 0, DC = 4, HCC = 2),
        high = c(S = 8364, I = 17764, R = 0, DC = 7, HCC = 3)
      ),
    type = "initial"
  ),
  
  "shared" = list(values = c(low = 30, high = 90), type = "parameter"),
  "Mu1.rate" = list(values = c(low = 0.04, high = 0.14), type = "parameter"),
  "treat_Cost_Param" = list(values = c(low = 16954, high = 45712), type = "parameter"),
  "retreat_Cost_Param" = list(values = c(low = 46885, high = 97376), type = "parameter"),
  "SSPleave.rate" = list(values = c(low = 0.23, high = 0.63), type = "parameter"),
  "MATleave.rate" = list(values = c(low = 0.06, high = 0.8), type = "parameter"),
  "ent.rate" = list(values = c(low = 1040, high = 1820), type = "parameter"),
  "Time horizon" = list(values = list(low = seq(0, 10, 1), high = seq(0, 50, 1)), type = "time")
)


summarize_run <- function(df) {
  df%>%
  summarize(
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
                  + Total_SSPIMP_Cost + Total_MAT_Cost) 
  ) }

# Runs the SIS HCV model for Baseline using ode solver from deSolve. 

baseline_out <- ode(y = base_initial, times = base_time, func = SIS_HCV_model, parms = base_parameters)
baseline_df  <- as.data.frame(baseline_out)
baseline_summary <- summarize_run(baseline_df) %>%
  mutate(Scenario = "Baseline", Param = "Base", Level = "Base")

all_results <- list()

for (scenario_name in names(intervention_list)) {
  
  # Apply intervention parameters
  scenario_parameters <- modifyList(base_parameters, intervention_list[[scenario_name]])
  
  # Runs the SIS HCV model for interventions using ode solver from deSolve and calculates incremental costs and DALYs averted from baselin.
  
  scenario_out <- ode(y = base_initial, times = base_time, func = SIS_HCV_model, parms = scenario_parameters)
  scenario_df  <- as.data.frame(scenario_out)
  scenario_summary <- summarize_run(scenario_df) %>%
    mutate(Scenario = scenario_name, Param = "Base", Level = "Base")
  
  inc_cost <- scenario_summary$Total_Cost - baseline_summary$Total_Cost
  daly_averted <- baseline_summary$Total_DALY - scenario_summary$Total_DALY
  icer <- inc_cost / daly_averted
  
  all_results[[length(all_results) + 1]] <- scenario_summary %>%
    transmute(
      Scenario,
      Comparison = "ScenarioBase_vs_TrueBaseline",
      Param, Level,
      Total_Cost, Total_DALY,
      Incremental_Cost = inc_cost,
      DALYs_Averted = daly_averted,
      ICER = icer
    )
  
  # One-Way Sensitivity Analysis (OWSA):
  # Iterates over each variable in 'owsa_list' and tests different values (initial conditions, parameters, or time).
  # Runs the SIS HCV model with ode() and summarizes results and compares them to the baseline.
  
  for (var_name in names(owsa_list)) {
    var <- owsa_list[[var_name]]
    
    for (level in names(var$values)) {
      
      current_initial    <- base_initial
      current_time       <- base_time
      current_parameters <- scenario_parameters
      
      if (var$type == "initial") current_initial <- var$values[[level]] 
      if (var$type == "parameter") current_parameters[[var_name]] <- var$values[[level]] 
      if (var$type == "time") current_time <- var$values[[level]]
      
      run_out <- ode(y = current_initial, times = current_time, func = SIS_HCV_model, parms = current_parameters)
      run_df  <- as.data.frame(run_out)
      run_sum <- summarize_run(run_df) %>%
        mutate(Scenario = scenario_name, Param = var_name, Level = level)
      
      inc_cost <- run_sum$Total_Cost - baseline_summary$Total_Cost
      daly_averted <- baseline_summary$Total_DALY - run_sum$Total_DALY
      icer <-  inc_cost / daly_averted
      
      all_results[[length(all_results) + 1]] <- run_sum %>%
        transmute(
          Scenario,
          Comparison = paste0("OWSA_vs_TrueBaseline_", var_name),
          Param, Level,
          Total_Cost, Total_DALY,
          Incremental_Cost = inc_cost,
          DALYs_Averted = daly_averted,
          ICER = icer
        )
    }
  }
}

owsa_icers <- bind_rows(
  tibble(
    Scenario = "Baseline",
    Comparison = "TrueBaseline",
    Param = "Base",
    Level = "Base",
    Total_Cost = baseline_summary$Total_Cost,
    Total_DALY  = baseline_summary$Total_DALY,
    Incremental_Cost = 0,
    DALYs_Averted = 0,
    ICER = NA_real_
  ),
  bind_rows(all_results)
)

ref_ICER <- owsa_icers %>%
  filter(Level == "Base" & Param == "Base") %>%
  select(Scenario, ICER_Ref = ICER)

# Summarize OWSA ICERs 
owsa_ICERs <- owsa_icers %>%
  filter(grepl("OWSA_vs_TrueBaseline", Comparison)) %>%
  select(Scenario, Parameter = Param, Level, ICER) %>%
  left_join(ref_ICER, by = "Scenario") %>%
  group_by(Scenario, Parameter) %>%
  summarize(
    ICER_Low  = if(any(Level == "low"))  ICER[Level == "low"] else NA_real_,
    ICER_High = if(any(Level == "high")) ICER[Level == "high"] else NA_real_,
    ICER_Ref  = first(ICER_Ref),
    Min_ICER  = if(all(is.na(c(ICER_Low, ICER_High)))) NA_real_ else min(ICER_Low, ICER_High, na.rm = TRUE),
    Max_ICER  = if(all(is.na(c(ICER_Low, ICER_High)))) NA_real_ else max(ICER_Low, ICER_High, na.rm = TRUE),
    Range     = Max_ICER - Min_ICER,
    Abs_Diff  = if(all(is.na(c(ICER_Low, ICER_High)))) NA_real_ else max(abs(ICER_Low - ICER_Ref), abs(ICER_High - ICER_Ref), na.rm = TRUE)
  ) %>%
  ungroup()
   

  write.csv(owsa_icers, "C:\\Users\\stubach\\OneDrive - State of Kansas, OITS\\Documents\\Dissertation\\DATA\\R program\\HCV CEA\\Final Commented Scripts\\Revision\\HCV_OWSA_results_REV2.csv", 
          row.names = FALSE) 
 
  write.csv(owsa_ICERs, "C:\\Users\\stubach\\OneDrive - State of Kansas, OITS\\Documents\\Dissertation\\DATA\\R program\\HCV CEA\\Final Commented Scripts\\Revision\\HCV_OWSA_plot_REV2.csv", 
            row.names = FALSE)
  