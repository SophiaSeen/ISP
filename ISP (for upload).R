library("writexl")
library("mvtBinaryEP")
library("FinancialMath")

num_simulations <-100000
  
i_success <- 0.535
ii_success <- 0.351
iii_success <- 0.673
nda_success <- 0.957

correlation <-0.4

patient_number <- 2000
indiv_cost <- .2
dose_cost <- indiv_cost*patient_number*0.25
half_dose_cost <- dose_cost/2

i_cost <- -8.5
ii_cost <- -7.3
iii_cost <- -15.7
nda_cost <- -0.3

i_months <- 32
ii_months <- 44
iii_months <-37
nda_months <- 13
months_for_approval <- i_months + ii_months + iii_months + nda_months

societal_gain <- 0.029810 + 0.022470

npv_trial = c()
QALY_trial = c()
Societal_trial = c()
Approval_number_trial= c()
invested_trial= c()
cost_for_successful_trial = c()
ROI_trial= c()
for (j in 1:num_simulations){
  z <- ranMvnXch(correlation, 100, nRep = 1, seed = NULL)
  i_to_ii <- sum(z >= qnorm(1-i_success, 0, 1))
  
  z <- ranMvnXch(correlation, i_to_ii, nRep = 1, seed = NULL)
  ii_to_iii <- sum(z >= qnorm(1-ii_success, 0, 1))
  
  z <- ranMvnXch(correlation, ii_to_iii, nRep = 1, seed = NULL)
  iii_to_nda <- sum(z >= qnorm(1-iii_success, 0, 1))
  
  z <- ranMvnXch(correlation, iii_to_nda, nRep = 1, seed = NULL)
  nda_to_approval <- sum(z >= qnorm(1-nda_success, 0, 1))
  
  Approval_number_trial <- c(nda_to_approval,Approval_number_trial)
  
  npv <- NPV(cf0=-i_cost*100,cf=c(-ii_cost*i_to_ii,-iii_cost*ii_to_iii,-nda_cost*iii_to_nda,dose_cost*nda_to_approval,dose_cost*nda_to_approval,dose_cost*nda_to_approval,dose_cost*nda_to_approval,dose_cost*nda_to_approval,dose_cost*nda_to_approval,dose_cost*nda_to_approval,dose_cost*nda_to_approval,dose_cost*nda_to_approval,dose_cost*nda_to_approval,half_dose_cost*nda_to_approval),times=c(i_months,i_months+ii_months,i_months+ii_months+iii_months,months_for_approval,months_for_approval+12*1,months_for_approval+12*2,months_for_approval+12*3,months_for_approval+12*4,months_for_approval+12*5,months_for_approval+12*6,months_for_approval+12*7,months_for_approval+12*8,months_for_approval+12*9,months_for_approval+12*9.5),i=1.1^(1/12)-1)
  
  npv_trial <- c(npv, npv_trial)
  
  invested <- -NPV(cf0=i_cost*100,cf=c(ii_cost*i_to_ii,iii_cost*ii_to_iii,nda_cost*iii_to_nda),times=c(i_months,i_months+ii_months,i_months+ii_months+iii_months),i=1.1^(1/12)-1)
  
  invested_trial <- c(invested,invested_trial)
  
  cost_for_successful <- invested/nda_to_approval
  
  cost_for_successful_trial <- c(cost_for_successful,cost_for_successful_trial)
  
  ROI <- npv/invested*100
  
  ROI_trial <- c(ROI, ROI_trial)
  
  q <- rnorm(nda_to_approval, mean=0.5, sd=0.2832)
  
  QALY <- patient_number*0.25*10.5*q
  
  QALY_gained <- sum(QALY)
  
  QALY_trial <-c(QALY_gained,QALY_trial)
  
  SQ <- societal_gain*QALY_gained
  
  Societal <- NPV(cf0=0,cf=c(SQ,SQ,SQ,SQ,SQ,SQ,SQ,SQ,SQ,SQ,SQ/2),times=c(months_for_approval, months_for_approval+12*1, months_for_approval+12*2, months_for_approval+12*3, months_for_approval+12*4, months_for_approval+12*5, months_for_approval+12*6, months_for_approval+12*7, months_for_approval+12*8, months_for_approval+12*9, months_for_approval+12*9.5),i=1.1^(1/12)-1)
  
  Societal_trial <-c(Societal,Societal_trial)
}

dataframe_npv_trial <- data.frame(npv_trial)
dataframe_invested_trial <- data.frame(invested_trial)
dataframe_cost_for_succuessful_trial <- data.frame(cost_for_successful_trial)
dataframe_ROI_trial <- data.frame(ROI_trial)
dataframe_QALY_trial <- data.frame(QALY_trial)
dataframe_Societal_trial <- data.frame(Societal_trial)
dataframe_Approval_number_trial <- data.frame(Approval_number_trial )
df_parameters <-data.frame (paramenter=c('patient number','phase 1 success','phase 2 success','phase 3 success','phase 4 success','correlation','drug price'), Chosen=c(patient_number,i_success,ii_success,iii_success,nda_success,correlation,indiv_cost))
df_results <- data.frame (measurement=c('mean approved','sd approved', 'mean npv','sd npv', 'mean invested', 'sd invested','mean cost of successful','sd cost of successful','mean ROI', 'sd ROI', 'mean QALY', 'sd QALY','mean societal','sd societal'), Results=c(mean(Approval_number_trial),sd(Approval_number_trial), mean(npv_trial),sd(npv_trial),mean(cost_for_successful_trial),sd(cost_for_successful_trial),mean(invested_trial),sd(invested_trial),mean(ROI_trial),sd(ROI_trial),mean(QALY_trial),sd(QALY_trial),mean(Societal_trial),sd(Societal_trial)))
names <- list('Parameters'=df_parameters,'Results' = df_results, 'Number of drugs approved' = dataframe_Approval_number_trial, 'NPV' = dataframe_npv_trial,'Invested' = dataframe_invested_trial,'Cost for Successful' = dataframe_cost_for_succuessful_trial, 'ROI' = dataframe_ROI_trial, 'QALY'= dataframe_QALY_trial,'Societal'= dataframe_Societal_trial)
write_xlsx(names, "path for file")
