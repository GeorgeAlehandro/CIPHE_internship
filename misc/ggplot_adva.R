library(ggplot2)

iso <- read.FCS("U:\\shiny_infinity_flow\\FCS\\split\\Plate1_Specimen_001_A3_A03_selection_correction_plate_1_0_target_CD3epsilon.fcs")@exprs
plot(iso[,"Event ID"], iso[,"G-PE-A"])
test_high_CD3 = iso[which (iso[,"G-PE-A"]> 2),]
plot(test_high_CD3[,"Event ID"], test_high_CD3[,"G-PE-A"])
plot(iso[,'Event ID'], iso[,"CD3epsilon.XGBoost"])
df = subset(iso,  iso[,'Event ID']== 1)
plot(hist(df[,"CD3epsilon.XGBoost_bgc"]))
plot(iso[,'UMAP1'],iso[,'UMAP2'], pch = 16)
iso = as.data.frame(iso)
iso[,'Event ID'] = as.factor(iso[,'Event ID'])
ggplot(iso, aes(x=UMAP1, y=UMAP2)) +
  geom_point(aes(col = `G-PE-A`))

iso <- read.FCS("U:\\shiny_infinity_flow\\FCS_background_corrected\\split\\Plate1_Specimen_001_A3_A03_selection_correction_plate_1_0_target_CD3epsilon.fcs")@exprs
plot(iso[,"Event ID"], iso[,"G-PE-A"])
test_high_CD3 = iso[which (iso[,"G-PE-A"]> 2),]
plot(test_high_CD3[,"Event ID"], test_high_CD3[,"G-PE-A"])
plot(iso[,'Event ID'], iso[,"CD3epsilon.XGBoost"])
df = subset(iso,  iso[,'Event ID']== 1)
plot(hist(df[,"CD3epsilon.XGBoost_bgc"]))
plot(iso[,'UMAP1'],iso[,'UMAP2'], pch = 16)
iso = as.data.frame(iso)
con[,'Event ID'] = as.factor(con[,'Event ID'])
ggplot(con, aes(x=UMAP1, y=UMAP2)) +
  geom_point(aes(col = `Event ID`))
con = as.data.frame(con)
ggplot(con, aes(x=UMAP1, y=UMAP2)) +
  geom_point(aes(col = `CD3epsilon.XGBoost_bgc`))+
  scale_colour_gradient2(  low = ("blue"),
                           mid = "yellow",
                           high = ("red"),
                           midpoint = 300)+
  geom_point(data = con[con$`CD3epsilon.XGBoost_bgc`>406.0723,], col = "red")


ggplot(, aes(x=UMAP1, y=UMAP2)) +
  geom_point(aes(col = `G-PE-A`))+
  scale_colour_gradient2(  low = ("blue"),
                           mid = "yellow",
                           high = ("red"),
                           midpoint = quantile(iso$`G-PE-A`, probs = 0.85)/2)+
  geom_point(data = iso[iso$`G-PE-A`>quantile(iso$`G-PE-A`, probs = 0.85),], col = "red")
quantile(iso$`CD30.XGBoost_bgc`, probs = 0.85)
ggplot(iso, aes(x=UMAP1, y=UMAP2)) +
  geom_point(aes(col = `CD30.XGBoost_bgc`))+
  geom_point(col = "blue")
kruskal.test(`CD30.XGBoost_bgc` ~ `Event ID`, data = iso)
