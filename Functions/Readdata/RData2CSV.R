#define directory and files
directory_DR = "/net/lunaria/work/rdamiano/legacy/Righetti_et_al_2019_sciadv/Data/DERIVED_Model_Frames/Output_1/"
file_GAM = "dfs_monthly_trsp_gam_pa_gridded_gr_bg(sit_ove).RData"
file_GLM = "dfs_monthly_trsp_glm_pa_gridded_gr_bg(sit_ove).RData"
file_RF = "dfs_monthly_trsp_rf_pa_gridded_gr_bg(sit_ove).RData"

#load data
GAMdata = get(load(paste(directory_DR,file_GAM,sep="")))
GLMdata = get(load(paste(directory_DR,file_GLM,sep="")))
RFdata = get(load(paste(directory_DR,file_RF,sep="")))

#store each month of each model
month = 1:12
folder = c('GAM', 'GLM', 'RF')

dir_out = "/net/kryo/work/ursho/PhD/Projects/Biomes/Scripts/Biomes_PROOCE/Data/00Probabilities/GAM/"

for (m in month){
  
  #file_out = paste(dir_out,"GAM/dfs_month_",m,"_trsp_gam_pa_gridded_gr_bg(sit_ove).csv",sep="")
  #write.csv(GAMdata[[m]], file = file_out)
  
  file_out = paste(dir_out,"GLM/dfs_month_",m,"_trsp_glm_pa_gridded_gr_bg(sit_ove).csv",sep="")
  write.csv(GLMdata[[m]], file = file_out) 
  
  file_out = paste(dir_out,"RF/dfs_month_",m,"_trsp_rf_pa_gridded_gr_bg(sit_ove).csv",sep="")
  write.csv(RFdata[[m]], file = file_out)
  print(paste("Done with month",m))
  
}


