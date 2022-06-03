cd_3 <- read.FCS("S:\\Mcyto\\Experiments\\R&D\\2022\\George-Alehandro_InfinityFlow\\output_annotated_by_gating_10_files\\FCS\\split\\Specimen_001_A3_A03_selection_testing_target_CD3.fcs")
plot(cd_3@exprs[,'UMAP1'],cd_3@exprs[,'UMAP2'], col = as.numeric(cut(cd_3_gpe,breaks = 10)))
cd_3_gpe <- cd_3@exprs[,'G-PE-A']
original <- read.FCS("S:\\Mcyto\\Experiments\\R&D\\2022\\George-Alehandro_gating\\selection_real_data\\Specimen_001_A3_A03_selection_testing.fcs")@exprs
plot(original[,'Event ID'], original[,'G-PE-A'])
View(original)
View(cd_3@exprs)
table(original[,'Event ID'])


new_subset <- read.FCS("S:\\Mcyto\\Experiments\\R&D\\2022\\George-Alehandro_equalsampling\\shiny_infinity_flow\\FCS_background_corrected\\split\\Plate3_equally_sampled_Specimen_001_B9_B09_selection_plate_3_until_c12__0_target_Rat IgG2b, k.fcs")@exprs
plot(new_subset[,'Event ID'], new_subset[,"I-A/I-E.XGBoost"], pch =1)
plot(new_subset[,'Event ID'], new_subset[,"G-PE-A"], pch =1)
df = subset(new_subset,  new_subset[,'Event ID']== 1)
plot(hist(df[,'G-PE-A']))



not_background_corrected_trans <- read.FCS("S:\\Mcyto\\Experiments\\R&D\\2022\\George-Alehandro_InfinityFlow\\output_annotated_by_gating_10_files\\FCS\\split\\Specimen_001_A3_A03_selection_testing_target_CD3.fcs")@exprs
plot(not_background_corrected_trans[,'Event ID'], not_background_corrected_trans[,'G-PE-A'])
df = subset(not_background_corrected_trans,  not_background_corrected_trans[,'Event ID']== 1)
plot(hist(df[,'G-PE-A']))

not_background_corrected_not_trans <- read.FCS("C:\\Users\\gsaad\\Desktop\\Specimen_001_A3_A03_selection_testing_target_CD3epsilon.fcs")@exprs
plot(not_background_corrected_not_trans[,'Event ID'], not_background_corrected_not_trans[,'G-PE-A'])
df = subset(not_background_corrected_not_trans,  not_background_corrected_not_trans[,'Event ID']== 4)
plot(hist(df[,'G-PE-A']))
View(not_background_corrected_not_trans)
View(not_background_corrected_trans)

done_here <- read.FCS("S:\\Mcyto\\Experiments\\R&D\\2022\\George-Alehandro_InfinityFlow\\shiny_test_infinity_flow\\FCS\\split\\Specimen_001_A3_A03_selection_testing_target_rIgM.fcs")@exprs
plot(done_here[,'Event ID'], done_here[,'G-PE-A'])
df = subset(done_here,  done_here[,'Event ID']== 6)
plot(hist(df[,'G-PE-A']))
View(done_here)
View(not_background_corrected_trans)

one_file <- read.FCS("C:\\Users\\gsaad\\Desktop\\gating\\selection_one_gate_for_a3_0\\Specimen_001_A3_A03_selection_one_gate_for_a3_0.fcs")
#plot(cd_3@exprs[,'Event ID'],cd_3@exprs[,'G-PE-A'])
df = subset(one_file@exprs,  one_file@exprs[,'Event ID']== 1)
plot(hist(df[,'G-PE-A']))
table(one_file@exprs[,'Event ID'])
in_many_files <- read.FCS("C:\\Users\\gsaad\\Desktop\\gating\\selection_with_many_files_0\\Specimen_001_A3_A03_selection_with_many_files_0.fcs")
#plot(cd_3@exprs[,'Event ID'],cd_3@exprs[,'G-PE-A'])
df = subset(in_many_files@exprs,  in_many_files@exprs[,'Event ID']== 1)
plot(hist(df[,'G-PE-A']))
table(in_many_files@exprs[,'Event ID'])
test <- list.files("C:\\Users\\gsaad\\Desktop\\gating\\selection_with_many_files_0\\", full.name = T)
for (a in test){
  print(a)
  b <- read.FCS(a)
  print(table(b@exprs[,'Event ID']))
  df = subset(b@exprs,  b@exprs[,'Event ID']== 1)
  plot(hist(df[,'G-PE-A']))
}
##FOR TWO files
in_two_files <- read.FCS("C:\\Users\\gsaad\\Desktop\\gating\\selection_2_files_0\\Specimen_001_A3_A03_selection_2_files_0.fcs")
#plot(cd_3@exprs[,'Event ID'],cd_3@exprs[,'G-PE-A'])
df = subset(in_two_files@exprs,  in_two_files@exprs[,'Event ID']== 1)
plot(hist(df[,'G-PE-A']))
table(in_two_files@exprs[,'Event ID'])

##FoR THREE FILES
in_3_files <- read.FCS("C:\\Users\\gsaad\\Desktop\\gating\\selection_in_3_files_0\\Specimen_001_A3_A03_selection_in_3_files_0.fcs")
#plot(cd_3@exprs[,'Event ID'],cd_3@exprs[,'G-PE-A'])
df = subset(in_3_files@exprs,  in_3_files@exprs[,'Event ID']== 1)
plot(hist(df[,'G-PE-A']))
table(in_3_files@exprs[,'Event ID'])

#For three files but in the middle
in_mid_files <- read.FCS("C:\\Users\\gsaad\\Desktop\\gating\\selection_selection_in_mid_0\\Specimen_001_A3_A03_selection_selection_in_mid_0.fcs")
#plot(cd_3@exprs[,'Event ID'],cd_3@exprs[,'G-PE-A'])
df = subset(in_mid_files@exprs,  in_mid_files@exprs[,'Event ID']== 1)
plot(hist(df[,'G-PE-A']))
table(in_mid_files@exprs[,'Event ID'])

not_corr <- read.FCS("S:\\Mcyto\\Experiments\\R&D\\2022\\George-Alehandro_gating\\selection_new_gate_0\\Specimen_001_A3_A03_selection_new_gate_0.fcs")
plot(not_corr@exprs[,'Event ID'],not_corr@exprs[,'G-PE-A'])
df = subset(not_corr@exprs,  not_corr@exprs[,'Event ID']== 1)
plot(hist(df[,'G-PE-A']))
table(not_corr@exprs[,'Event ID'])

#corr <- read.FCS("S:\\Mcyto\\Experiments\\R&D\\2022\\George-Alehandro_gating\\selection_PLATE1_0\\Specimen_001_A3_A03_selection_correction_plate_1_0.fcs")@exprs
corr <- read.FCS("C:\\Users\\gsaad\\Desktop\\gating\\selection_in_3_files_0\\Specimen_001_A3_A03_selection_in_3_files_0.fcs")
plot(corr@exprs[,'Event ID'],corr@exprs[,'G-PE-A'],pch = 16)
df = subset(corr@exprs,  corr@exprs[,'Event ID']== 1)
#Create a function to generate a continuous color palette
rbPal <- colorRampPalette(c('red','blue'))

#This adds a column of color values
# based on the y values
corr@exprs$Col <- rbPal(10)[as.numeric(cut(corr@exprs[,'G-PE-A'],breaks = 10))]
plot(corr@epxrs[,'UMAP1'],corr@epxrs[,'UMAP2'],pch = 20,col = dat$Col)
plot(hist(df[,'G-PE-A']))
table(corr@exprs[,'Event ID'])

