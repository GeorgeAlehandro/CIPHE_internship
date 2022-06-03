#Decompensate

compensated <- read.FCS("C:\\Users\\gsaad\\Desktop\\gating\\selection_test_13_5_1\\Specimen_001_A3_A03_selection_test_13_5_1.fcs")
decompensate(compensated)

decompensated_file <- decompensate(compensated)
write.FCS(decompensated_file, 'decompensated_test')


#TEST PLOT ON OUTPUT OF INFINITY ON DECOMPENSATED DATA
after_infinity_decompensated <- read.FCS("S:\\Mcyto\\Experiments\\R&D\\2022\\George-Alehandro_InfinityFlow\\shiny_test_infinity_flow\\FCS\\split\\decompensated_test_target_rIgM.fcs")@exprs

#after_infinity_decompensated <- read.FCS("C:\\Users\\gsaad\\Desktop\\Specimen_001_A3_A03_selection_testing_target_CD3epsilon.fcs")@exprs
plot(after_infinity_decompensated[,'Event ID'], after_infinity_decompensated[,'G-PE-A'])
df = subset(after_infinity_decompensated,  after_infinity_decompensated[,'Event ID']== 1)
plot(hist(df[,'G-PE-A']))


decompensate(decompensated_file)

#TEST PLOT ONT OUTPUT OF INFINITY ON COMPENSATED DATA
after_infinity_compensated <- read.FCS("S:\\Mcyto\\Experiments\\R&D\\2022\\George-Alehandro_InfinityFlow\\shiny_test_infinity_flow\\FCS\\split\\Specimen_001_A3_A03_selection_test_13_5_1_target_rIgMt.fcs")@exprs

#after_infinity_decompensated <- read.FCS("C:\\Users\\gsaad\\Desktop\\Specimen_001_A3_A03_selection_testing_target_CD3epsilon.fcs")@exprs
plot(after_infinity_compensated[,'Event ID'], after_infinity_compensated[,'G-PE-A'])
df = subset(after_infinity_compensated,  after_infinity_compensated[,'Event ID']== 1)

