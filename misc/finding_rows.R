original_a3 <- read.FCS("S:\\Mcyto\\Experiments\\R&D\\2022\\George-Alehandro_InfinityFlow\\DATA\\18-Oct-12_CONTEXT_KD\\FCS Plate 1\\18-Oct-12_CONTEXT_KD\\Plaque 1\\Specimen_001_A3_A03.fcs")
annotated_a3 <- read.FCS("S:\\Mcyto\\Experiments\\R&D\\2022\\George-Alehandro_gating\\Plate1_cleaning_test\\Specimen_001_A3_A03.fcs")
infinity_a3 <- read.FCS("S:\\Mcyto\\Experiments\\R&D\\2022\\George-Alehandro_InfinityFlow\\output_annotated_by_gating_10_files\\FCS\\split\\Specimen_001_A3_A03_selection_testing_target_CD3.fcs")
View(original_a3@exprs)
View(annotated_a3@exprs)
View(infinity_a3@exprs)

all_files <- list.files("C:\\Users\\gsaad\\Desktop\\split\\", full.names = T)
for (file in all_files){
  t<-file
  file <- read.FCS(file)@exprs
  plot(file[,'Event ID'], file[,'G-PE-A'], xlab = t)
  
}
