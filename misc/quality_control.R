library(flowCore)
library(ggplot2)
original <- read.FCS("C:\\Users\\gsaad\\Desktop\\Plate2_equally_sampled_Specimen_001_E10_E10_selection_plate_2_until_f12_0_target_F4-80.fcs")@exprs
sumbsampled <- read.FCS("S:\\Mcyto\\Experiments\\R&D\\2022\\George-Alehandro_equalsampling\\shiny_infinity_flow\\FCS_background_corrected\\split\\Plate2_equally_sampled_Specimen_001_EE10_E10_selection_plate2_until_F12_36__rea_0_target_F4-80.fcs")@exprs
gated <- read.FCS("S:\\Mcyto\\Experiments\\R&D\\2022\\George-Alehandro_gating\\selection_plate_2\\Specimen_001_E10_E10_selection_plate2_until_F12_36__rea_0.fcs")@exprs
plot(original[,'Event ID'], original[,'F4/80.XGBoost'], pch =1)
df = subset(original,  original[,'Event ID']== 3)
plot(hist(df[,'G-PE-A']))
#after infinity flow
after_inf <- read.FCS("S:\\Mcyto\\Experiments\\R&D\\2022\\George-Alehandro_equalsampling\\shiny_infinity_flow\\FCS_background_corrected\\split\\Plate1_equally_sampled_Specimen_001_A3_A03_selection_correction_plate_1_0.fcs_target_CD3epsilon.fcs")@exprs
plot(after_inf[,'Event ID'], after_inf[,'G-PE-A'], pch =1)
plot(after_inf[,'Event ID'], after_inf[,"CD3epsilon.XGBoost_bgc"], pch =1)
df = subset(after_inf,  after_inf[,'Event ID']== 9)
plot(hist(df[,'G-PE-A']))
after_inf = as.data.frame(after_inf)
ordering = after_inf[order(after_inf[,"CD3epsilon.XGBoost_bgc"], decreasing = TRUE),]
#ordering[,'Event ID']
ggplot(after_inf, aes(x=UMAP1, y=UMAP2)) +
  geom_point(aes(col = `G-PE-A`))+
  scale_colour_gradient2(  low = ("blue"),
                           mid = "yellow",
                           high = ("red"),
                           midpoint = quantile(after_inf$`G-PE-A`, probs = 0.95)/2)+
  geom_point(data = after_inf[after_inf$`G-PE-A`>quantile(after_inf$`G-PE-A`, probs = 0.95),], col = "red")

ggplot(after_inf, aes(x=UMAP1, y=UMAP2)) +
  geom_point(aes(col = CD3epsilon.XGBoost_bgc))+
  scale_colour_gradient2(  low = ("blue"),
                           mid = "yellow",
                           high = ("red"),
                           midpoint = quantile(after_inf$CD3epsilon.XGBoost_bgc, probs = 0.85)/2)+
  geom_point(data = after_inf[after_inf$CD3epsilon.XGBoost_bgc>quantile(after_inf$CD3epsilon.XGBoost_bgc, probs = 0.85),], col = "red")
after_inf[,'Event ID'] = as.factor(after_inf[,'Event ID'])
ggplot(after_inf, aes(x=UMAP1, y=UMAP2)) +
  geom_point(aes(col = `Event ID`))
#Equally sampled A3
View(read.FCS("S:\\Mcyto\\Experiments\\R&D\\2022\\George-Alehandro_equalsampling\\output\\plate1\\equally_sampled_Specimen_001_A3_A03_selection_correction_plate_1_0.fcs.fcs")@exprs)
after_gating <- (read.FCS("S:\\Mcyto\\Experiments\\R&D\\2022\\George-Alehandro_gating\\new_gating\\plate1\\Specimen_001_H12_H12_selection_23_5_plate_1_0.fcs")@exprs)
org <-(read.FCS("S:\\Mcyto\\Experiments\\R&D\\2022\\George-Alehandro_InfinityFlow\\DATA\\18-Oct-12_CONTEXT_KD\\FCS Plate 3\\18-Oct-12_CONTEXT_KD\\Plaque 3\\Specimen_001_F12_F12.fcs")@exprs)
new <-(read.FCS("S:\\Mcyto\\Experiments\\R&D\\2022\\George-Alehandro_equalsampling\\output\\equally_sampled_Specimen_001_F12_F12_selection_after_d12_plate3__0.fcs")@exprs)
