library(flowCore)
library(flowWorkspace)
library(cipheCytoExploreR)
library(dplyr) 
options(java.parameters = c("-XX:+UseConcMarkSweepGC", "-Xmx8192m"))
library(xlsx)
library(rJava)
file_to_pick_from_1 <- sample(list.files("S:\\Mcyto\\Experiments\\R&D\\2022\\George-Alehandro_equalsampling\\shiny_infinity_flow\\FCS_background_corrected\\split\\", pattern = "Plate1",full.names = T),5)
file_to_pick_from_2 <- sample(list.files("S:\\Mcyto\\Experiments\\R&D\\2022\\George-Alehandro_equalsampling\\shiny_infinity_flow\\FCS_background_corrected\\split\\", pattern = "Plate2",full.names = T),5)
file_to_pick_from_3 <- sample(list.files("S:\\Mcyto\\Experiments\\R&D\\2022\\George-Alehandro_equalsampling\\shiny_infinity_flow\\FCS_background_corrected\\split\\", pattern = "Plate3",full.names = T),5)

all_files_to_pick <- c(file_to_pick_from_1, file_to_pick_from_2,file_to_pick_from_3)
test <- read.FCS("S:\\Mcyto\\Experiments\\R&D\\2022\\George-Alehandro_equalsampling\\shiny_infinity_flow\\FCS_background_corrected\\split\\Plate1_equally_sampled_Specimen_001_A2_A02_selection_new_gate_0.fcs_target_Armenian Hamster IgG.fcs")
gs_pop_get_stats(test, type = pop.quantiles)
file <- "S:\\Mcyto\\Experiments\\R&D\\2022\\George-Alehandro_equalsampling\\shiny_infinity_flow\\FCS_background_corrected\\split\\"
gs<- cyto_setup("S:\\Mcyto\\Experiments\\R&D\\2022\\George-Alehandro_equalsampling\\shiny_infinity_flow\\FCS_background_corrected\\split\\")
first.quartile <- function(fr){
  chnls <- colnames(fr)
  res <- matrixStats::colQuantiles(exprs(fr), probs = 0.25)
  names(res) <- chnls
  res
}
third.quartile <- function(fr){
  chnls <- colnames(fr)
  res <- matrixStats::colQuantiles(exprs(fr), probs = 0.75)
  names(res) <- chnls
  res
}
mediane <-function(fr){
  chnls <- colnames(fr)
  res <- apply(exprs(fr),2,median)
  names(res) <- chnls
  res
}
standard_deviation <-function(fr){
  chnls <- colnames(fr)
  res <- apply(exprs(fr),2,sd)
  names(res) <- chnls
  res
}
first.decile <- function(fr){
  chnls <- colnames(fr)
  
  for (i in 1:length(chnls)){
    res <- c(res,quantile(fr@exprs[,i], prob = 0.1, type = 5))
  }
  names(res) <- chnls
  res
}
ninth.decile <- function(fr){
  chnls <- colnames(fr)
  
  for (i in 1:length(chnls)){
    res <- c(res,quantile(fr@exprs[,i], prob = 0.9, type = 5))
  }
  names(res) <- chnls
  res
}
modefunc <- function(x){
  tabresult <- tabulate(x)
  themode <- which(tabresult == max(tabresult))
  if(sum(tabresult == max(tabresult))>1) themode <- NA
  return(themode)
}
mode.df <- function(fr){
  chnls <- colnames(fr)
  res<-apply(exprs(fr),2,modefunc)
  names(res) <- chnls
  res
}
mean_pops <-function(fr){
  chnls <- colnames(fr)
  res <- apply(exprs(fr),2,mean)
  names(res) <- chnls
  res
}


mode_stats = gs_pop_get_stats(gs, type = mode.df)
mode_stats$pop <- NULL
standard_deviation_stats = gs_pop_get_stats(gs, type = standard_deviation)
standard_deviation_stats$pop <- NULL
mean_stats = gs_pop_get_stats(gs, type = mean_pops)
mfi_stats = gs_pop_get_stats(gs, type = pop.MFI)
mfi_stats$pop <- NULL
mediane_stats = gs_pop_get_stats(gs, type = mediane)
mediane_stats$pop <- NULL
first_decile_stats = gs_pop_get_stats(gs, type = first.decile)
first_decile_stats$pop <- NULL
ninth_decile_stats = gs_pop_get_stats(gs, type = ninth.decile)
ninth_decile_stats$pop <- NULL
first_quartile_stats = gs_pop_get_stats(gs, type = first.quartile)
first_quartile_stats$pop <- NULL
third_quartile_stats = gs_pop_get_stats(gs, type = third.quartile)
third_quartile_stats$pop <- NULL
options(java.parameters = "-Xmx8000m")
jgc <- function()
{
  .jcall("java/lang/System", method = "gc")
}  
write.xlsx(mode_stats, file = paste(
  getwd(),
  paste0('here_test_local', '_log_file.xlsx'),
  sep = '/'), sheetName = "mode", append = FALSE)
write.xlsx(mode_stats, file = paste(
  getwd(),
  paste0('here_test_local', '_log_file.xlsx'),
  sep = '/'), sheetName = "mode", append = FALSE)
jgc()
write.xlsx(standard_deviation_stats, file = paste(
  getwd(),
  paste0('here_test_local', '_log_file.xlsx'),
  sep = '/'), 
  sheetName="SD", append=TRUE)

write.xlsx(mfi_stats, file = paste(
  getwd(),
  paste0('here_test_local', '_log_file.xlsx'),
  sep = '/'), 
  sheetName="MFI", append=TRUE)
jgc()
write.xlsx(mediane_stats, file = paste(
  getwd(),
  paste0('here_test_local', '_log_file.xlsx'),
  sep = '/'), 
  sheetName="Mediane", append=TRUE)
jgc()
write.xlsx(first_decile_stats, file = paste(
  getwd(),
  paste0('here_test_local_2', '_log_file.xlsx'),
  sep = '/'), 
  sheetName="First decile", append=F)
jgc()
write.xlsx(ninth_decile_stats, file = paste(
  getwd(),
  paste0('here_test_local_2', '_log_file.xlsx'),
  sep = '/'), 
  sheetName="Ninth decile", append=TRUE)
jgc()

write.xlsx(first_quartile_stats, file = paste(
  getwd(),
  paste0('here_test_local_2', '_log_file.xlsx'),
  sep = '/'), 
  sheetName="First quartile", append=TRUE)
jgc()
write.xlsx(third_quartile_stats, file = paste(
  getwd(),
  paste0('here_test_local_2', '_log_file.xlsx'),
  sep = '/'), 
  sheetName="Third quartile", append=TRUE)
jgc()