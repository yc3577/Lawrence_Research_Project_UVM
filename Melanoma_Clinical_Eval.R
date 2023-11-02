# R 3.6.2
library(openxlsx)
library('survival');
library('survminer')
library(ggfortify)
library('stringr')
library(ConsensusClusterPlus);
library(ade4);

rm(list=ls())
######################### load data ###########################

########################## consensus clustering ###############
consensus_result_root = '../Data/ConsensusClustering_SelfNorm_500_0.8/'
clusterAlg_params = c("pam");
distance_params = c("spearman");
MAXK = 6;

newpath_base = '../Results/'
dir.create(newpath_base);

newpath_base = '../Results/Clinical_Evaluation/'
dir.create(newpath_base);

########################################## KM plots ###########
skcm_clinical_fn = paste('../Data/Meta/','skcm_tcga_pan_can_atlas_2018/data_clinical_patient.csv',sep = '')
skcm_clinical_data = read.csv(skcm_clinical_fn,header=TRUE,row.names = 1)
skcm_clinical_data$OS_STATUS = as.character(skcm_clinical_data$OS_STATUS)
skcm_clinical_data$PFS_STATUS = as.character(skcm_clinical_data$PFS_STATUS)
skcm_clinical_data$OS_STATUS[skcm_clinical_data$OS_STATUS=='0:LIVING']=0
skcm_clinical_data$OS_STATUS[skcm_clinical_data$OS_STATUS=='1:DECEASED']=1
skcm_clinical_data$OS_STATUS = as.numeric(skcm_clinical_data$OS_STATUS)
skcm_clinical_data$PFS_STATUS[skcm_clinical_data$PFS_STATUS=='0:CENSORED']=0
skcm_clinical_data$PFS_STATUS[skcm_clinical_data$PFS_STATUS=='1:PROGRESSION']=1
skcm_clinical_data$PFS_STATUS = as.numeric(skcm_clinical_data$PFS_STATUS)

uvm_clinical_fn = paste('../Data/Meta/','uvm_tcga_pan_can_atlas_2018/data_clinical_patient.csv',sep = '')
uvm_clinical_data = read.csv(uvm_clinical_fn,header=TRUE,row.names = 1)
uvm_clinical_data$OS_STATUS = as.character(uvm_clinical_data$OS_STATUS)
uvm_clinical_data$PFS_STATUS = as.character(uvm_clinical_data$PFS_STATUS)
uvm_clinical_data$OS_STATUS[uvm_clinical_data$OS_STATUS=='0:LIVING']=0
uvm_clinical_data$OS_STATUS[uvm_clinical_data$OS_STATUS=='1:DECEASED']=1
uvm_clinical_data$OS_STATUS = as.numeric(uvm_clinical_data$OS_STATUS)
uvm_clinical_data$PFS_STATUS[uvm_clinical_data$PFS_STATUS=='0:CENSORED']=0
uvm_clinical_data$PFS_STATUS[uvm_clinical_data$PFS_STATUS=='1:PROGRESSION']=1
uvm_clinical_data$PFS_STATUS = as.numeric(uvm_clinical_data$PFS_STATUS)

for(cls_ii in 1:length(clusterAlg_params)){
  
  clusterAlg_param = clusterAlg_params[cls_ii];
  
  for(dist_jj in 1:length(distance_params)){
    
    distance_param = distance_params[dist_jj];
    
    
    cat('\n')
    print("########################################################################################################")
    print(newpath_base)
    print("########################################################################################################")
    cat('\n')
    
    for (kk in c(256)){ #### input
      
      result_summary_fn = paste(consensus_result_root,'/', clusterAlg_param, '-', distance_param, '-', kk, '_subtyping_summary.csv', sep = '')
      result_summary = read.csv(result_summary_fn,header = T, row.names = 1)
      skcm_result_summary = merge(skcm_clinical_data,result_summary,by="row.names")
      uvm_result_summary = merge(uvm_clinical_data,result_summary,by="row.names")
      
      for(nn in c(2:MAXK)){
        skcm_formular = paste('Surv(OS_MONTHS, OS_STATUS) ~ ',paste("K",nn,"Label",sep = '.'),sep = '')
        skcm_fit <- survfit(as.formula(skcm_formular), data = skcm_result_summary)
        plot_skcm = ggsurvplot(
          skcm_fit,                     # survfit object with calculated statistics.
          data = skcm_result_summary,  # data used to fit survival curves.
          risk.table = TRUE,       # show risk table.
          pval = TRUE,             # show p-value of log-rank test.
          # conf.int = TRUE,         # show confidence intervals for
          # point estimaes of survival curves.
          # xlim = c(0,2000),        # present narrower X axis, but not affect
          # survival estimates.
          # break.time.by = 500,     # break X axis in time intervals by 500.
          # ggtheme = theme_minimal(), # customize plot and risk table with a theme.
          risk.table.y.text.col = T, # colour risk table text annotations.
          risk.table.y.text = FALSE, # show bars instead of names in text annotations,
          xlab = "Month",
          ylab = "Probability"
          # in legend of risk table
        )

        pdf(paste(newpath_base,'/',"TCGA_SKCM_ConsensusCluster_KK=",as.character(nn), "_OS.pdf",sep = ''),onefile=F)
        print(plot_skcm)
        dev.off()
        
        skcm_formular = paste('Surv(PFS_MONTHS, PFS_STATUS) ~ ',paste("K",nn,"Label",sep = '.'),sep = '')
        skcm_fit <- survfit(as.formula(skcm_formular), data = skcm_result_summary)
        plot_skcm = ggsurvplot(
          skcm_fit,                     # survfit object with calculated statistics.
          data = skcm_result_summary,  # data used to fit survival curves.
          risk.table = TRUE,       # show risk table.
          pval = TRUE,             # show p-value of log-rank test.
          # conf.int = TRUE,         # show confidence intervals for
          # point estimaes of survival curves.
          # xlim = c(0,2000),        # present narrower X axis, but not affect
          # survival estimates.
          # break.time.by = 500,     # break X axis in time intervals by 500.
          # ggtheme = theme_minimal(), # customize plot and risk table with a theme.
          risk.table.y.text.col = T, # colour risk table text annotations.
          risk.table.y.text = FALSE, # show bars instead of names in text annotations,
          xlab = "Month",
          ylab = "Probability"
          # in legend of risk table
        )
        
        pdf(paste(newpath_base,'/',"TCGA_SKCM_ConsensusCluster_KK=",as.character(nn), "_PFS.pdf",sep = ''),onefile=F)
        print(plot_skcm)
        dev.off()
        
        
        uvm_formular = paste('Surv(OS_MONTHS, OS_STATUS) ~ ',paste("K",nn,"Label",sep = '.'),sep = '')
        uvm_fit <- survfit(as.formula(uvm_formular), data = uvm_result_summary)
        plot_uvm = ggsurvplot(
          uvm_fit,                     # survfit object with calculated statistics.
          data = uvm_result_summary,  # data used to fit survival curves.
          risk.table = TRUE,       # show risk table.
          pval = TRUE,             # show p-value of log-rank test.
          # conf.int = TRUE,         # show confidence intervals for
          # point estimaes of survival curves.
          # xlim = c(0,2000),        # present narrower X axis, but not affect
          # survival estimates.
          # break.time.by = 500,     # break X axis in time intervals by 500.
          # ggtheme = theme_minimal(), # customize plot and risk table with a theme.
          risk.table.y.text.col = T, # colour risk table text annotations.
          risk.table.y.text = FALSE, # show bars instead of names in text annotations,
          xlab = "Month",
          ylab = "Probability"
          # in legend of risk table
        )
        
        pdf(paste(newpath_base,'/',"TCGA_UVM_ConsensusCluster_KK=",as.character(nn), "_OS.pdf",sep = ''),onefile=F)
        print(plot_uvm)
        dev.off()
        
        uvm_formular = paste('Surv(PFS_MONTHS, PFS_STATUS) ~ ',paste("K",nn,"Label",sep = '.'),sep = '')
        uvm_fit <- survfit(as.formula(uvm_formular), data = uvm_result_summary)
        plot_uvm = ggsurvplot(
          uvm_fit,                     # survfit object with calculated statistics.
          data = uvm_result_summary,  # data used to fit survival curves.
          risk.table = TRUE,       # show risk table.
          pval = TRUE,             # show p-value of log-rank test.
          # conf.int = TRUE,         # show confidence intervals for
          # point estimaes of survival curves.
          # xlim = c(0,2000),        # present narrower X axis, but not affect
          # survival estimates.
          # break.time.by = 500,     # break X axis in time intervals by 500.
          # ggtheme = theme_minimal(), # customize plot and risk table with a theme.
          risk.table.y.text.col = T, # colour risk table text annotations.
          risk.table.y.text = FALSE, # show bars instead of names in text annotations,
          xlab = "Month",
          ylab = "Probability"
          # in legend of risk table
        )
        
        pdf(paste(newpath_base,'/',"TCGA_UVM_ConsensusCluster_KK=",as.character(nn), "_PFS.pdf",sep = ''),onefile=F)
        print(plot_uvm)
        dev.off()
      }
    }
  }
}
