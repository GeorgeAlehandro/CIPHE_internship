


before_infinity_flow <- read.FCS("S:\\Mcyto\\Experiments\\R&D\\2022\\George-Alehandro_gating\\selection_PLATE1_0\\Specimen_001_A3_A03_selection_correction_plate_1_0.fcs")@exprs
before_infinity_flow <- as.data.frame(before_infinity_flow)
plot(before_infinity_flow[,'Event ID'], before_infinity_flow[,'G-PE-A'])

after_inf_not_bg <- read.FCS("S:\\Mcyto\\Experiments\\R&D\\2022\\George-Alehandro_InfinityFlow\\second_test_10_plate1\\FCS\\split\\Specimen_001_A3_A03_selection_correction_plate_1_0_target_CD3epsilon.fcs")@exprs
after_inf_not_bg <- as.data.frame(after_inf_not_bg)
plot(after_inf_not_bg[,'Event ID'], after_inf_not_bg[,'G-PE-A'])



after_inf <- read.FCS("S:\\Mcyto\\Experiments\\R&D\\2022\\George-Alehandro_InfinityFlow\\second_test_10_plate1\\FCS_background_corrected\\split\\Specimen_001_A3_A03_selection_correction_plate_1_0_target_CD3epsilon.fcs")@exprs
after_inf <- as.data.frame(after_inf)
rbPal <- colorRampPalette(c('red','blue'))

#This adds a column of color values
# based on the y values
after_inf$Col <- rbPal(10)[as.numeric(cut(after_inf[,'G-PE-A'],breaks = 10))]

plot(after_inf[,'UMAP1'],after_inf[,'UMAP2'],pch = 20,col = after_inf[,'Event ID'])
plot(after_inf[,'Event ID'], after_inf[,'G-PE-A'])

for (file in list.files("S:\\Mcyto\\Experiments\\R&D\\2022\\George-Alehandro_InfinityFlow\\second_test_10_plate1\\FCS_background_corrected\\split\\", full.names = T)){
  View(read.FCS("S:\\Mcyto\\Experiments\\R&D\\2022\\George-Alehandro_InfinityFlow\\second_test_10_plate1\\FCS_background_corrected\\split\\Specimen_001_A1_A01_selection_new_gate_0_target_Blank.fcs")@exprs)
  print(file)
}
View(read.FCS("S:\\Mcyto\\Experiments\\R&D\\2022\\George-Alehandro_InfinityFlow\\second_test_10_plate1\\FCS_background_corrected\\split\\Specimen_001_A1_A01_selection_new_gate_0_target_Blank.fcs")@exprs)

interesting <- read.FCS("S:\\Mcyto\\Experiments\\R&D\\2022\\George-Alehandro_InfinityFlow\\second_test_10_plate1\\FCS_background_corrected\\split\\Specimen_001_A4_A04_selection_new_gate_0_target_CD80.fcs")@exprs
interesting <- as.data.frame(interesting)
plot(interesting[,'Event ID'], interesting[,'G-PE-A'])


#5allina njarreb CD80
interesting_original <- read.FCS("S:\\Mcyto\\Experiments\\R&D\\2022\\George-Alehandro_gating\\selection_PLATE1_0\\Specimen_001_A4_A04_selection_new_gate_0.fcs")@exprs
interesting_original <- as.data.frame(interesting_original)
plot(interesting_original[,'Event ID'], interesting_original[,'G-PE-A'])
