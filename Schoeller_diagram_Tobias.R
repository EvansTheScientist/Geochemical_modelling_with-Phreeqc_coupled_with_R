#+eval=FALSE

### Skript 7 # Projektarbeit T.Schnepper # 03.2021

# Line-Plot / Custom Schoeller Diagramm


library(tidyverse)


## Test-Daten simulieren

# Zwei Proben, Ionenkonzentrationen liegen mit Einheit mg/l vor

df_mgl <- Tidy_data_Greece

colorBlindBlack8  <- c("#E69F00", "#56B4E9", "#009E73", 
                       "#F0E442", "#0072B2", "#D55E00", "#CC79A7")

#df_mgl <- data.frame("Sample" = c("1", "2"),
#                      "Mg" = c(42.4, 27.0),
#                      "Ca" = c(92.0, 71.0),
#                      "Na" = c(45.5, 32.1),
#                      "K" = c(3.7, 6.4),
#                      "HCO3" = c(330.7, 330.7),
#                      "Cl" = c(138.65, 75.1),
#                      "SO4" = c(51.23, 64.76),
#                      "pH" = c(7.9, 8.2),
#                      "TDS" = c(688.8, 502.4))


#  meq/l aus mg/l berechnen: Massenkonzentration / Molare Masse * Valenz

#df_meql <- df_mgl %>% summarize("Sample" = Sample,
#                      "Mg" = round(Mg / 24.31 * 2, digits = 2),
#                      "Ca" = round(Ca / 40.08 * 2, digits = 2),
#                      "Na" = round(Na / 22.99 * 1, digits = 2),
#                      "K" = round(K / 39.10 * 1, digits = 2),
#                      "HCO3" = round(HCO3 / 61.02 * 1, digits = 2),
#                      "Cl" = round(Cl / 35.45 * 1, digits = 2),
#                      "SO4" = round(SO4 / 96.06 * 2, digits = 2),
#                      "pH" = pH,
#                      "TDS" = TDS)
                      
df_meql <- df_mgl %>% summarize("Date" = Date,
                                "ID" = ID,          
                      "Mg" = round(Mg * 2, digits = 2),
                      "Ca" = round(Ca * 2, digits = 2),
                      "Na" = round(Na  * 1, digits = 2),
                      "K" = round(K  * 1, digits = 2),
                      "HCO3" = round(HCO3  * 1, digits = 2),
                      "Cl" = round(Cl  * 1, digits = 2),
                      "SO4" = round(SO4  * 2, digits = 2),
                      "pH" = pH) %>%
  filter(complete.cases(.) & ID == "ΥΝP_222" | 
           complete.cases(.) & ID == "MWS")

df_meql$ID_ext <- paste(df_meql$ID, df_meql$Date)


# Datensatz für das Schoeller-Diagramm vorbereiten. Die ausgewählten Ionen können individuell festgelegt werden.
# Die Notation ist wichtig, um die Ionen später im richtigen Format auf der x-Achse darzustellen.
df_schoeller <- df_meql %>% summarise("ID" = ID,
                                      "ID_ext" = ID_ext,
                                      "Mg^2^+phantom()" = Mg,
                                      "Ca^2^+phantom()" = Ca,
                                      "Ca^2^+phantom()~+~Mg^2^+phantom()" = Ca + Mg,
                                      "Na^+phantom()~+~K^+phantom()" = Na + K,
                                      "Cl^-phantom()" = Cl,
                                      "HCO[3]^-phantom()" = HCO3,
                                      "SO[4]^2^-phantom()" = SO4)

# Tabelle aus dem wide in das long format umwandeln
df_schoeller_gather <- df_schoeller %>% gather(Ion, value, -ID, -ID_ext) 

df_schoeller_gather$Ion <- factor(df_schoeller_gather$Ion,
                                 levels = c("Mg^2^+phantom()",
                                            "Ca^2^+phantom()",
                                            "Ca^2^+phantom()~+~Mg^2^+phantom()",
                                            "Na^+phantom()~+~K^+phantom()",
                                            "Cl^-phantom()",
                                            "HCO[3]^-phantom()",
                                            "SO[4]^2^-phantom()"),
                                 labels=c("Mg^2^+phantom()",
                                          "Ca^2^+phantom()",
                                          "Ca^2^+phantom()~+~Mg^2^+phantom()",
                                          "Na^+phantom()~+~K^+phantom()",
                                          "Cl^-phantom()",
                                          "HCO[3]^-phantom()",
                                          "SO[4]^2^-phantom()"))

## Schoeller diagramm als Linienplot
df_schoeller_gather %>% ggplot(aes(Ion, value, colour = ID, group = ID_ext, fill=ID)) +
  geom_point(size=5, pch=21, colour="black") + # Punkt-Plot
  theme_bw()+ # Layout anpassen
  scale_y_log10(limits= c(0.1, 50)) + # Y-Achse formatieren
  annotation_logticks(sides = "lr", size = 1) + # Logs anpassen
  labs(y=expression(paste("Concentrations", " (meq kgw" ^-{1}, ")"))) +
  xlab("") +
  scale_x_discrete(labels = rlang::parse_exprs) +  #"Übersetzt" die Schreibweise in die korrekte Ionenbeschriftung
  guides(x.sec = "axis", y.sec = "axis", labels = NULL)+
  theme(
    axis.title = element_text(size = 15, face = "bold", family="verdana"),
    axis.text = element_text(size = 15, face = "bold", family="verdana"),
    axis.line = element_line(size = 1),
    legend.text = element_text(size = 15, face = "bold", family="verdana"),
    legend.position = c(0.18, 0.2),
    legend.title = element_text(size = 15, face = "bold", family="verdana"),
    legend.background = element_blank()) + # Achsenbeschriftungen
  scale_color_manual(values=c("#e8b56f", "#78aaee"),
                     name = "Samples",
                     labels = c("Mine water", "Groundwater"))+
scale_fill_manual(name="Samples", values=c("#e8b56f", "#78aaee"), labels=c("Mine water", "Groundwater"))+
  geom_line(size=0.5)  # Linien-Plot 
