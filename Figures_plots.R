setwd("C:/Users/aivanova/Documents/Orbitrap Data/settings comparison")
library(ggplot2)
library(dplyr)
library(tidyr)
library(tidyverse)
library(patchwork)
library(forcats)
library(grid)
library(mzR)
library(msdata)

#####################################################################################################################################
# FIGURE 3
# INTENSITY MODULATION PATTERNS
color_180 <- c("#882255", "#CC6677")
color_342AB <- c("#88CCEE", "#44AA99")
color_342BF <- c("#44AA99", "#332288")

real_precursors_180 <- c(181.064141, 181.085921)
real_precursors_342BF <- c(343.081231, 343.12885)
real_precursors_342AB <- c(343.074706, 343.081231)


formulas_180 <- c("180E (181.0641)", "180G (181.0859)")
formulas_342AB <- c("342A (343.0747)", '342B (343.0812)')
formulas_342BF <- c('342B (343.0812)', "342F (343.1288)")

massDeviation <- function(mz1, mz2) { abs(mz1-mz2)/((mz2+mz1)/2)*1000000 }
findPeak <- function(mz, pl, dev, type) {
  match <- -1
  all <- sapply(pl[, 1], function(x) {massDeviation(mz, x)})
  matches <- which(all<=dev)
  if (length(matches)==1 || (type=="sum" & length(matches)>1)) {
    match <- matches
  } else if (length(matches)>1) {
    if (type=="maxInt") {
      match <- matches[which.max(pl[matches, 2])]
    } else if (type=="nearest") {
      match <- matches[which.min(matches)]
    }
  }
  match
}

isoWinPrecursorPlot_data <- function(precursors, labels, ms2, scan_vector, dev, matchType="maxInt") {
  o <- c()
  l <- c()
  for (p in 1:length(precursors)) {
    tmplist <- c()
    for (i in scan_vector) {
      pl <- peaks(ms2, i)
      index <- findPeak(precursors[p], pl, dev, matchType)
      if (length(index)>1) {
        tmplist <- rbind(tmplist, append(append(i, mean(pl[index, ][, 1])), sum(pl[index, ][, 2])))
      } else if (length(index)==1 & index!=-1) {
        tmplist <- rbind(tmplist, append(i, pl[index, ]))
      }
    }
    if (!is.null(tmplist)) {
      tmplist[, 2] <- rep(mean(tmplist[, 2]), nrow(tmplist))
      tmplist <- cbind(tmplist, tmplist[, 3]/max(tmplist[, 3]))
      if (!is.null(labels)) {
        l <- append(l, rep(labels[p], nrow(tmplist)))
      } else {
        l <- append(l, rep(as.character(precursors[p]), nrow(tmplist)))
      }
      colnames(tmplist) <- c("Scan", "mz", "Intensity", "Normalized_Intensity")
      o <- rbind(o, tmplist)
    } else {
      if (!is.null(formulas)) {
        cat("WARNING: Precursor ", formulas[p], " is not found in any scans!\n")
      } else {
        cat("WARNING: Precursor ", precursors[p], " is not found in any scans!\n")
      }
    }
  }
  data <- data.frame(
    scan <- o[, 1],
    mz <- as.character(o[, 2]),
    intensity <- o[, 3],
    normalizedIntensity <- o[, 4],
    labels <- l
  )
  colnames(data) <- c("Scan", "mz", "Intensity", "Normalized_Intensity", "Peak")
  data 
}

scientific_10 <- function(x) {
  sci_x <- scales::scientific_format()(x)
  parse(text=gsub("e\\+", "%*%10^", replace(sci_x, sci_x=="0e+00", "0")))
}


getRowNumber <- function(l) {
  if (!is.null(nrow(l))) {
    nrow(l)
  } else {
    if (length(l)>=1) {
      1
    } else {
      0
    }
  }
}

def_180_el <- openMSfile("./180_settings_comparison/180E+G/1-180E+G_default.mzML")
ms2_180_el <- isoWinPrecursorPlot_data(real_precursors_180, formulas_180, def_180_el, c(1:81), 6)

def_el_180EG <- ggplot(data=ms2_180_el, aes(x=Scan, y=Normalized_Intensity, color=Peak)) +
  geom_point() + geom_line() +
  theme(text = element_text(size = 15),
        plot.title = element_text(size = 17),
        axis.title = element_text(size = 13),
        panel.background = element_blank(), 
        axis.line = element_line(colour = "black"),
        legend.position = 'right') +
  scale_color_manual(values=color_180 ) +
  #scale_y_continuous(labels = scientific_10)+
  ggtitle('Orbitrap Elite')+
  ylab("Normalized intensity")+
  xlab("Scan number")+
  labs(color = "Isobaric compound")


## 342 A+B

def_342AB_el <- openMSfile("./342_settings_comparison/342A+B/342A+B_3microsc.mzML")
ms2_342AB_el <- isoWinPrecursorPlot_data(real_precursors_342AB, formulas_342AB, def_342AB_el, c(1:81), 8)

def_el_342AB <- ggplot(data=ms2_342AB_el, aes(x=Scan, y=Normalized_Intensity, color=Peak)) +
  geom_point() + geom_line() +
  theme(text = element_text(size = 15),
        plot.title = element_text(size = 17),
        axis.title = element_text(size = 13),
        panel.background = element_blank(), 
        axis.line = element_line(colour = "black"),
        legend.position = 'right') +
  scale_color_manual(values=color_342AB) +
  #scale_y_continuous(labels = scientific_10)+
  ggtitle('Orbitrap Elite: 342A+B')+
  ylab("Normalized intensity")+
  xlab("Scan number")+
  labs(color = "Isobaric compound")

## 342B+F
def_342BF_el <- openMSfile("./342_settings_comparison/342B+F/1-342B+F_default.mzML")
ms2_342BF_el <- isoWinPrecursorPlot_data(real_precursors_342BF, formulas_342BF, def_342BF_el, c(1:81), 8)

def_el_342BF <- ggplot(data=ms2_342BF_el, aes(x=Scan, y=Normalized_Intensity, color=Peak)) +
  geom_point() + geom_line() +
  theme(text = element_text(size = 15),
        plot.title = element_text(size = 17),
        axis.title = element_text(size = 13),
        panel.background = element_blank(), 
        axis.line = element_line(colour = "black"),
        legend.position = 'right') +
  scale_color_manual(values=color_342BF) +
  #scale_y_continuous(labels = scientific_10)+
  ggtitle('Orbitrap Elite: 342B+342F')+
  ylab("Normalized intensity")+
  xlab("Scan number")+
  labs(color = "Isobaric compound")


## 2 Orbitrap Exploris
def_180_ex <- openMSfile("./Orbi_Exploris/180_settings_comparison/180E+G/1-180E+G_default_29.mzML")

ms2_180_ex <- isoWinPrecursorPlot_data(real_precursors_180, formulas_180, def_180_ex, c(82:163), 5)

def_ex_180EG <- ggplot(data=ms2_180_ex, aes(x=Scan, y=Normalized_Intensity, color=Peak)) +
  geom_point() + geom_line() +
  theme(text = element_text(size = 15),
        plot.title = element_text(size = 17),
        axis.title = element_text(size = 13),
        panel.background = element_blank(), 
        axis.line = element_line(colour = "black"),
        legend.position = 'right') +
  scale_color_manual(values=color_180 ) +
  #scale_y_continuous(labels = scientific_10)+
  ggtitle('Orbitrap Exploris')+
  ylab("Normalized intensity")+
  xlab("Scan number")+
  labs(color = "Isobaric compound")

## 342 A+B

def_342AB_ex <- openMSfile("./Orbi_Exploris/342_settings_comparison/342A+B/1-342A+B_default_03.mzML")

ms2_342AB_ex <- isoWinPrecursorPlot_data(real_precursors_342AB, formulas_342AB, def_342AB_ex, c(81:163), 4)

def_ex_342AB <- ggplot(data=ms2_342AB_ex, aes(x=Scan, y=Normalized_Intensity, color=Peak)) +
  geom_point() + geom_line() +
  theme(text = element_text(size = 15),
        plot.title = element_text(size = 17),
        axis.title = element_text(size = 13),
        panel.background = element_blank(), 
        axis.line = element_line(colour = "black"),
        legend.position = 'right') +
  scale_color_manual(values=color_342AB) +
  #scale_y_continuous(labels = scientific_10)+
  ggtitle('Orbitrap Exploris: 342A+342B')+
  ylab("Normalized intensity")+
  xlab("Scan number")+
  labs(color = "Isobaric compound")

## 342B+F
def_342BF_ex <- openMSfile("./Orbi_Exploris/342_settings_comparison/342B+F/1-342B+F_default_05.mzML")
ms2_342BF_ex <- isoWinPrecursorPlot_data(real_precursors_342BF, formulas_342BF, def_342BF_ex, c(82:162), 5)

def_ex_342BF <- ggplot(data=ms2_342BF_ex, aes(x=Scan, y=Normalized_Intensity, color=Peak)) +
  geom_point() + geom_line() +
  theme(text = element_text(size = 15),
        plot.title = element_text(size = 17),
        axis.title = element_text(size = 13),
        panel.background = element_blank(), 
        axis.line = element_line(colour = "black"),
        legend.position = 'right') +
  scale_color_manual(values=color_342BF) + 
  #scale_y_continuous(labels = scientific_10)+
  ggtitle('Orbitrap Exploris')+
  ylab("Normalized intensity")+
  xlab("Scan number")+
  labs(color = "Isobaric compound")

theme_combplot_top <- theme(axis.title = element_blank(), ##for combined plot
                            axis.text.y = element_blank(),##for combined plot
                            axis.ticks.y = element_blank(),
                            axis.line.y = element_blank())##for combined plot)

theme_combplot_bottom <- theme(axis.title.y = element_blank(), ##for combined plot
                               axis.text.y = element_blank(),##for combined plot
                               axis.ticks.y = element_blank(),
                               axis.line.y = element_blank())


Figure3 <- ((def_el_342AB  + theme(axis.title = element_blank())) + (def_el_180EG + theme_combplot_top) + (def_el_342BF + theme_combplot_top)) /
  ((def_ex_342AB + theme(axis.title = element_blank())) + (def_ex_180EG + theme_combplot_top)+ (def_ex_342BF+ theme_combplot_top))+
  plot_layout(guides = "collect") & 
  theme(legend.position = 'right',
        legend.title = element_blank(),
        plot.title = element_blank())

##################################################################################################################################################
# FIGURE 4
#SIMILARITY SCORE BAR CHARTS HORIZONTAL

color_180EG <- c("180G" = "#CC6677", "180E" = "#882255")
color_342AB <- c("342A" = "#88CCEE", "342B"="#44AA99")
color_342BF <- c("342B"="#44AA99", "342F"= "#332288")

levels = c("Default", "0.01 m/z step", "0.04 m/z step", "0.1 m/z step", "1 microscan", "3 microscans", "R = 15 000", "R = 30 000",
           "R = 60 000", "0.4 m/z isolation width", "0.7 m/z isolation width", "2 m/z isolation width", 
           "50% AGC", "200% AGC", "NCE = 35", "NCE = 50", "NCE = 65")

El_180EG <- read.csv("./180_settings_comparison/180E+G/180EG_comparison_table_inprog.csv")
settings <- read.csv("./settings_name_180EG.csv")
levels <- settings$Settings

El_180EG <- El_180EG %>%
  select(Instrument, settings, time, int_thres_mix, starts_with("Isobar_")) %>%
  pivot_longer(cols = starts_with("Isobar_"),
               names_to = c("Isobar", ".value"), 
               names_pattern = "Isobar_(\\d+)_(.*)")%>%
  mutate(Isobar = recode(Isobar, `1` = "180E", `2` = "180G"))%>%
  left_join(settings, by = "settings") %>%
  mutate(unassign_int = 1 - assign_int)

El_180EG$Settings <- factor(El_180EG$Settings, levels = (levels))

score_180EG_El <-
  ggplot(El_180EG) + 
  geom_bar(aes(fill=(Isobar), y=score, x=Settings), position="dodge", stat="identity", width = 0.6) +

  scale_x_discrete(guide = guide_axis(angle = 45))+
  theme(text = element_text(size = 11),
        axis.title = element_text(size = 14),
        axis.text.x = element_text(angle = 0),
        panel.background = element_blank())+ 
  scale_fill_manual(values=color_180EG, guide = guide_legend(title = "Isobar"), breaks=c("180E","180G"))+
  theme(axis.text.x = element_text(angle = 0, vjust = 1, hjust=0.5))+
  labs( y = "Similarity score")


El_342AB <- read.csv("./342_settings_comparison/342A+B/342AB_comparison_table_inprog.csv")
settings <- read.csv("./settings_name.csv")
levels <- settings$Settingss

El_342AB <- El_342AB %>%
  select(Instrument, settings, time, int_thres_mix, starts_with("Isobar_")) %>%
  pivot_longer(cols = starts_with("Isobar_"),
               names_to = c("Isobar", ".value"), 
               names_pattern = "Isobar_(\\d+)_(.*)")%>%
  mutate(Isobar = recode(Isobar, `1` = "342A", `2` = "342B"))%>%
  left_join(settings, by = "settings") %>%
  mutate(unassign_int = 1 - assign_int)%>%
  mutate(incorr_int = 0)

El_342AB$Settings <- factor(El_342AB$Settings, levels = levels)

score_342AB_El <-
  ggplot(El_342AB) + 
  geom_bar(aes(fill=(Isobar), y=score, x=Settings), position="dodge", stat="identity", width = 0.6) +
  scale_x_discrete(guide = guide_axis(angle = 45))+
  theme(text = element_text(size = 15),
        axis.title = element_text(size = 15),
        axis.text.x = element_text(angle = 0, size = 13),
        panel.background = element_blank())+
  scale_fill_manual(values=color_342AB, guide = guide_legend(title = "Isobar"), breaks = c("342A", "342B"))+
  theme(axis.text.x = element_text(angle = 0, vjust = 1, hjust=0.5))+
  labs( y = "Similarity score")

El_342BF <- read.csv("./342_settings_comparison/342B+F/342BF_comparison_table_inprog.csv")
settings <- read.csv("./settings_name_342BF.csv")
levels <- settings$Settings

El_342BF <- El_342BF %>%
  select(Instrument, settings, time, int_thres_mix, starts_with("Isobar_")) %>%
  pivot_longer(cols = starts_with("Isobar_"),
               names_to = c("Isobar", ".value"), 
               names_pattern = "Isobar_(\\d+)_(.*)")%>%
  mutate(Isobar = recode(Isobar, `1` = "342B", `2` = "342F"))%>%
  left_join(settings, by = "settings") %>%
  mutate(unassign_int = 1 - assign_int)

El_342BF$Settings <- factor(El_342BF$Settings, levels = levels)

score_342BF_El <-
  ggplot(El_342BF) + 
  geom_bar(aes(fill=(Isobar), y=score, x=Settings), position="dodge", stat="identity", width = 0.6) +
  #coord_flip()+
  scale_x_discrete(guide = guide_axis(angle = 45))+
  theme(text = element_text(size = 15),
        axis.title = element_text(size = 15),
        #axis.ticks.y=element_blank(),
        axis.text.x = element_text(angle = 0, size = 13),
        panel.background = element_blank())+
  scale_fill_manual(values=color_342BF, guide = guide_legend(title = "Isobar"), breaks = c("342B", "342F"))+ #no guide_legend(reverse = TRUE)
  theme(axis.text.x = element_text(angle = 0, vjust = 1, hjust=0.5))+
  labs( y = "Similarity score")

Ex_180EG <- read.csv("C:/Users/aivanova/Documents/Orbitrap Data/settings comparison/Orbi_Exploris/180_settings_comparison/180EG_comparison_table_inprog.csv")
settings <- read.csv("C:/Users/aivanova/Documents/Orbitrap Data/settings comparison/settings_name_180EG_exp.csv")

levels <- settings$Settings

Ex_180EG <- Ex_180EG %>%
  select(Instrument, settings, time, int_thres_mix, starts_with("Isobar_")) %>%
  pivot_longer(cols = starts_with("Isobar_"),
               names_to = c("Isobar", ".value"), 
               names_pattern = "Isobar_(\\d+)_(.*)")%>%
  mutate(Isobar = recode(Isobar, `1` = "180E", `2` = "180G"))%>%
  left_join(settings, by = "settings") %>%
  mutate(unassign_int = 1 - assign_int)

Ex_180EG$Settings <- factor(Ex_180EG$Settings, levels = levels)

score_180EG_Ex <-
  ggplot(Ex_180EG) + 
  geom_bar(aes(fill=(Isobar), y=score, x=Settings), position="dodge", stat="identity", width = 0.6) +
  #coord_flip()+
  scale_x_discrete(guide = guide_axis(angle = 45))+
  theme(text = element_text(size = 15),
        axis.title = element_text(size = 15),
        #axis.ticks.y=element_blank(),
        axis.text.x = element_text(angle = 0, size = 13),
        panel.background = element_blank())+ 
  scale_fill_manual(values=color_180EG, guide = guide_legend(title = "Isobar"), breaks=c("180E","180G"))+
  theme(axis.text.x = element_text(angle = 0, vjust = 1, hjust=0.5))+
  labs( y = "Similarity score")

Ex_342AB <- read.csv("./Orbi_Exploris/342_settings_comparison/342A+B/342AB_comparison_table_inprog.csv")
settings <- read.csv("./settings_name_342AB_ex.csv")
levels <- settings$Settings

Ex_342AB <- Ex_342AB %>%
  select(Instrument, settings, time, int_thres_mix, starts_with("Isobar_")) %>%
  pivot_longer(cols = starts_with("Isobar_"),
               names_to = c("Isobar", ".value"), 
               names_pattern = "Isobar_(\\d+)_(.*)")%>%
  mutate(Isobar = recode(Isobar, `1` = "342A", `2` = "342B"))%>%
  left_join(settings, by = "settings") %>%
  mutate(unassign_int = 1 - assign_int)


Ex_342AB$Settings <- factor(Ex_342AB$Settings, levels = levels)

score_342AB_Ex <-
  ggplot(Ex_342AB) + 
  geom_bar(aes(fill=(Isobar), y=score, x=Settings), position="dodge", stat="identity", width = 0.6) +
  #coord_flip()+
  scale_x_discrete(guide = guide_axis(angle = 45))+
  theme(text = element_text(size = 15),
        axis.title = element_text(size = 15),
        #axis.ticks.y=element_blank(),
        axis.text.x = element_text(angle = 0, size = 13),
        panel.background = element_blank())+
  scale_fill_manual(values=color_342AB, guide = guide_legend(title = "Isobar"), breaks = c("342A", "342B"))+ #no guide_legend(reverse = TRUE)
  theme(axis.text.x = element_text(angle = 0, vjust = 1, hjust=0.5))+
  labs( y = "Similarity score")

Ex_342BF <- read.csv("./Orbi_Exploris/342_settings_comparison/342B+F/342BF_comparison_table_inprog.csv")
settings <- read.csv("./settings_name_342AB_ex.csv")
levels <- settings$Settings

Ex_342BF <- Ex_342BF %>%
  select(Instrument, settings, time, int_thres_mix, starts_with("Isobar_")) %>%
  pivot_longer(cols = starts_with("Isobar_"),
               names_to = c("Isobar", ".value"), 
               names_pattern = "Isobar_(\\d+)_(.*)")%>%
  mutate(Isobar = recode(Isobar, `1` = "342B", `2` = "342F"))%>%
  left_join(settings, by = "settings") %>%
  mutate(unassign_int = 1 - assign_int)

Ex_342BF$Settings <- factor(Ex_342BF$Settings, levels = levels)

score_342BF_Ex <-
  ggplot(Ex_342BF) + 
  geom_bar(aes(fill=(Isobar), y=score, x=Settings), position="dodge", stat="identity", width = 0.6) +
  #coord_flip()+
  scale_x_discrete(guide = guide_axis(angle = 45))+
  theme(text = element_text(size = 15),
        axis.title = element_text(size = 15),
        #axis.ticks.y=element_blank(),
        axis.text.x = element_text(angle = 0, size = 13),
        panel.background = element_blank())+ 
  scale_fill_manual(values=color_342BF, guide = guide_legend(title = "Isobar"), breaks = c("342B", "342F"))+ #no guide_legend(reverse = TRUE)
  theme(axis.text.x = element_text(angle = 0, vjust = 1, hjust=0.5))+
  labs( y = "Similarity score")


score_180EG_Ex_mod <- score_180EG_Ex +  theme(
  text = element_text(size = 13),
  axis.text.y=element_blank(),
  axis.title.y = element_blank(),
  axis.text.x=element_blank(),
  axis.ticks.y = element_blank(),
  axis.title.x = element_blank())

score_180EG_El_mod <- score_180EG_El +
  theme(axis.text.y=element_blank(),
        axis.ticks.y = element_blank(),
        axis.title.y = element_blank(),
        text = element_text(size = 13),
        axis.text.x = element_text(size = 13),
        axis.title.x = element_blank(),
        legend.title = element_blank())

score_180EG <- (score_180EG_Ex_mod/score_180EG_El_mod) + 
  plot_layout(guides = "collect") & theme(legend.position = "bottom")


score_342AB_El_mod <- score_342AB_El+
  theme(#axis.text.x=element_blank(),
    axis.title.y = element_text(size = 13),
    text = element_text(size = 13),
    axis.text.x = element_text(size = 13),
    axis.title.x = element_blank(),
    legend.title = element_blank())


score_342AB_Ex_mod <- score_342AB_Ex+
  theme(axis.text.x=element_blank(),
        axis.title.y = element_text(size = 13),
        text = element_text(size = 13),
        axis.title.x = element_blank())

score_342AB <- (score_342AB_Ex_mod/score_342AB_El_mod) + 
  plot_layout(guides = "collect")& theme(legend.position = "bottom")


score_342BF_El_mod <-score_342BF_El+ 
theme(axis.text.y=element_blank(),
            axis.ticks.y = element_blank(),
            axis.title.y = element_blank(),
            text = element_text(size = 13),
            axis.text.x = element_text(size = 13),
            axis.title.x = element_blank(),
            legend.title = element_blank())

score_342BF_Ex_mod <- score_342BF_Ex +
  theme(axis.text.y=element_blank(),
        text = element_text(size = 13),
        axis.title.y = element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.y = element_blank(),
        axis.title.x = element_blank())
score_342BF <- (score_342BF_Ex_mod/score_342BF_El_mod) +
  plot_layout(guides = "collect")& theme(legend.position = "bottom")

Figure_4 <- (score_342AB | score_180EG |score_342BF) & theme(legend.title  = element_blank())

################################################################################################################################################
#FIGURE 5
#INTENSITY DISCTRIBUTION BAR CHARTS HORIZONTAL
color_180EG_int <- c("180E correctly assigned intensity" = "#882255", "180G correctly assigned intensity" = "#CC6677",
                     "180E unassigned intensity" = "darkgrey",  "180G unassigned intensity" = "darkgrey",
                     "180E incorrectly assigned intensity" = "darkred", "180G incorrectly assigned intensity" = "darkred")

color_342AB_int <- c("342A correctly assigned intensity" = "#88CCEE", "342B correctly assigned intensity" = "#44AA99",
                     "342A unassigned intensity" = "darkgrey",  "342B unassigned intensity" = "darkgrey",
                     "342A incorrectly assigned intensity" = "darkred", "342B incorrectly assigned intensity" = "darkred")

color_342BF_int <- c("342B correctly assigned intensity" = "#44AA99", "342F correctly assigned intensity" = "#332288",
                     "342B unassigned intensity" = "darkgrey",  "342F unassigned intensity" = "darkgrey",
                     "342B incorrectly assigned intensity" = "darkred", "342F incorrectly assigned intensity" = "darkred")
##180EG Elite

El_180EG_resh_int <- El_180EG %>%
  select(settings, Isobar, score, Settings, ends_with("_int")) %>%
  pivot_longer(cols = ends_with("_int"), names_to = "type", names_pattern = "^(.*)_int$", values_to = "value")%>%
  mutate(int_type = interaction(Isobar, type))


int_type_lvl <- c("180E.incorr", "180G.incorr", "180E.unassign", "180G.unassign", "180E.assign", "180G.assign")
El_180EG_resh_int$int_type <- factor(El_180EG_resh_int$int_type, levels = int_type_lvl, 
                                     labels = c("180E incorrectly assigned intensity", "180G incorrectly assigned intensity",
                                                "180E unassigned intensity", "180G unassigned intensity",
                                                "180E correctly assigned intensity", "180G correctly assigned intensity"))

plots_set_180EG <- as.list(c())

for (i in 1:(length(unique(El_180EG_resh_int$Settings)))){
  setting = unique(El_180EG_resh_int$Settings)[i]
  plot = 
    ggplot(as.data.frame(El_180EG_resh_int[El_180EG_resh_int$Settings == setting, ]), aes(fill = int_type, y=value, x=(Isobar)))+
    geom_bar(stat = "identity", position = "stack")+
    scale_fill_manual(values= color_180EG_int)+
    ylim(0, 1.7)+
    labs(x = setting)+
    scale_x_discrete(guide = guide_axis(angle = 45))+
    theme(text = element_text(size = 12), 
          #axis.text.x=element_text(angle = 45, vjust = 1, hjust=0.5),
          axis.text = element_blank(),
          axis.ticks.x=element_blank(),
          axis.ticks.y=element_blank(),
          panel.background = element_blank(),
          axis.title = element_blank(),
          legend.position = "none"
    )
  #+
    #coord_flip()+
  #  theme(#axis.title.y = element_text(angle = 0, hjust = 0.5, vjust = 0.5), 
  #    axis.title.x =  element_text(angle = 45, hjust = 1, vjust = 1))
  plots_set_180EG[[i]] = plot
}

plots_comb_180EG_int_El_noset <- (plots_set_180EG[[1]]|plots_set_180EG[[2]]|plots_set_180EG[[3]]|plots_set_180EG[[4]]|
                              plots_set_180EG[[9]]|plots_set_180EG[[13]]|plots_set_180EG[[6]]|plots_set_180EG[[11]]|
                              plots_set_180EG[[15]]|plots_set_180EG[[17]]|plots_set_180EG[[5]]|plots_set_180EG[[10]]|plots_set_180EG[[7]]|
                              plots_set_180EG[[8]]|plots_set_180EG[[12]]|
                              plots_set_180EG[[14]]|plots_set_180EG[[16]]) + 
  plot_layout(guides = "collect")& theme(axis.title.x = element_blank(),
        plot.margin = 
          margin(t=0.1, l=0, r=0, b=0, unit = "cm"))


##180 E+G EXPLORIS

Ex_180EG_resh_int <- Ex_180EG %>%
  select(settings, Isobar, score, Settings, ends_with("_int")) %>%
  pivot_longer(cols = ends_with("_int"), names_to = "type", names_pattern = "^(.*)_int$", values_to = "value")%>%
  mutate(int_type = interaction(Isobar, type))


int_type_lvl <- c("180E.incorr", "180G.incorr", "180E.unassign", "180G.unassign", "180E.assign", "180G.assign")
Ex_180EG_resh_int$int_type <- factor(Ex_180EG_resh_int$int_type, levels = int_type_lvl, 
                                     labels = c("180E incorrectly assigned intensity", "180G incorrectly assigned intensity",
                                                "180E unassigned intensity", "180G unassigned intensity",
                                                "180E correctly assigned intensity", "180G correctly assigned intensity"))
plots_set_180EG_ex <- as.list(c())

for (i in 1:(length(unique(Ex_180EG_resh_int$Settings)))){
  setting = unique(Ex_180EG_resh_int$Settings)[i]
  plot = 
    ggplot(as.data.frame(Ex_180EG_resh_int[Ex_180EG_resh_int$Settings == setting, ]), aes(fill = int_type, y=value, x=(Isobar)))+
    geom_bar(stat = "identity", position = "stack")+
    scale_fill_manual(values= color_180EG_int)+
    ylim(0, 1.7)+
    labs(x = setting)+
    scale_x_discrete(guide = guide_axis(angle = 45))+
    theme(text = element_text(size = 12), 
          #axis.text.x=element_text(angle = 45, vjust = 1, hjust=0.5),
          axis.text = element_blank(),
          axis.ticks.x=element_blank(),
          axis.ticks.y=element_blank(),
          panel.background = element_blank(),
          axis.title = element_blank(),
          legend.position = "none"
    )
  #+
  #coord_flip()+
  #  theme(#axis.title.y = element_text(angle = 0, hjust = 0.5, vjust = 0.5), 
  #    axis.title.x =  element_text(angle = 45, hjust = 1, vjust = 1))
  plots_set_180EG_ex[[i]] = plot
}

plots_comb_180EG_int_Ex_noset <- (plots_set_180EG_ex[[1]]|plots_set_180EG_ex[[2]]|plots_set_180EG_ex[[3]]|plots_set_180EG_ex[[4]]|
                                    plots_set_180EG_ex[[7]]|plots_set_180EG_ex[[11]]|plots_set_180EG_ex[[6]]|plots_set_180EG_ex[[9]]|
                                    plots_set_180EG_ex[[13]]|plots_set_180EG_ex[[17]]|plots_set_180EG_ex[[5]]|plots_set_180EG_ex[[8]]|
                                    plots_set_180EG_ex[[16]]|plots_set_180EG_ex[[15]]|plots_set_180EG_ex[[10]]|
                                    plots_set_180EG_ex[[12]]|plots_set_180EG_ex[[14]]) + 
  plot_layout(guides = "collect") & theme(axis.title.x = element_blank(),
                                          plot.margin = margin(t=0.1, l=0, r=0, b=0, unit = "cm"))

#### 342A+B Elite

El_342AB_resh_int <- El_342AB %>%
  select(settings, Isobar, score, Settings, ends_with("_int")) %>%
  pivot_longer(cols = ends_with("_int"), names_to = "type", names_pattern = "^(.*)_int$", values_to = "value")%>%
  mutate(int_type = interaction(Isobar, type))

int_type_lvl <- c("342A.incorr", "342B.incorr", "342A.unassign", "342B.unassign", "342A.assign", "342B.assign")
El_342AB_resh_int$int_type <- factor(El_342AB_resh_int$int_type, levels = int_type_lvl, 
                                     labels = c("342A incorrectly assigned intensity", "342B incorrectly assigned intensity",
                                                "342A unassigned intensity", "342B unassigned intensity",
                                                "342A correctly assigned intensity", "342B correctly assigned intensity"))

plots_set <- list(c())
for (i in 1:(length(unique(El_342AB_resh_int$Settings)))){
  setting = unique(El_342AB_resh_int$Settings)[i]
  plot = 
    ggplot(as.data.frame(El_342AB_resh_int[El_342AB_resh_int$Settings == setting, ]), aes(fill = int_type, y=value, x=(Isobar)))+
    geom_bar(stat = "identity", position = "stack")+
    scale_fill_manual(values= color_342AB_int)+
    ylim(0, 1.7)+
    labs(x = setting)+
    scale_x_discrete(guide = guide_axis(angle = 45))+
    theme(text = element_text(size = 12), 
          #axis.text.x=element_text(angle = 45, vjust = 1, hjust=0.5),
          axis.text = element_blank(),
          axis.ticks.x=element_blank(),
          axis.ticks.y=element_blank(),
          panel.background = element_blank(),
          axis.title = element_blank(),
          legend.position = "none"
    )
  #+
  #coord_flip()+
  #  theme(#axis.title.y = element_text(angle = 0, hjust = 0.5, vjust = 0.5), 
  #    axis.title.x =  element_text(angle = 45, hjust = 1, vjust = 1))
  plots_set[[i]] = plot
}
i = 1
setting = unique(El_342AB_resh_int$Settings)[i]

plots_set[[i]] <- 
  ggplot(as.data.frame(El_342AB_resh_int[El_342AB_resh_int$Settings == setting, ]), aes(fill = int_type, y=value, x=(Isobar)))+
  geom_bar(stat = "identity", position = "stack")+
  scale_fill_manual(values= color_342AB_int)+
  ylim(0, 1.7)+
  labs(x = setting)+
  scale_x_discrete(guide = guide_axis(angle = 45))+
  theme(text = element_text(size = 13), 
        axis.text.x=element_blank(),
        #axis.text = element_blank(),
        axis.ticks.x=element_blank(),
        #axis.ticks.y=element_blank(),
        panel.background = element_blank(),
        axis.title = element_blank(),
        legend.position = "none"
  )#+
  #coord_flip()+
  #theme(#axis.title.y = element_text(angle = 0, hjust = 0.5, vjust = 0.5), 
  #  axis.title.x = element_text(angle = 45, hjust = 1, vjust = 1))



plots_comb_342AB_int_El <- (plots_set[[1]]|plots_set[[2]]|plots_set[[3]]|plots_set[[4]]|plots_set[[7]]|
                              plots_set[[11]]|plots_set[[6]]|plots_set[[9]]|plots_set[[13]]|plots_set[[17]]|
                              plots_set[[5]]|plots_set[[8]]|plots_set[[15]]|plots_set[[16]]|plots_set[[10]]|
                              plots_set[[12]]|plots_set[[14]] + theme(legend.position = "none")) + 
  plot_layout(guides = "collect") & theme(axis.title.x = element_blank(), legend.position = "none",
                                                                            plot.margin = 
                                                                              margin(t=0.1, l=0, r=0, b=0, unit = "cm"))


# 342 A+B Exploris

Ex_342AB_resh_int <- Ex_342AB %>%
  select(settings, Isobar, score, Settings, ends_with("_int")) %>%
  pivot_longer(cols = ends_with("_int"), names_to = "type", names_pattern = "^(.*)_int$", values_to = "value")%>%
  mutate(int_type = interaction(Isobar, type))

int_type_lvl <- c("342A.incorr", "342B.incorr", "342A.unassign", "342B.unassign", "342A.assign", "342B.assign")
Ex_342AB_resh_int$int_type <- factor(Ex_342AB_resh_int$int_type, levels = int_type_lvl, 
                                     labels = c("342A incorrectly assigned intensity", "342B incorrectly assigned intensity",
                                                "342A unassigned intensity", "342B unassigned intensity",
                                                "342A correctly assigned intensity", "342B correctly assigned intensity"))

plots_set_342AB_ex <- list(c())

for (i in 1:(length(unique(Ex_342AB_resh_int$Settings)))){
  setting = unique(Ex_342AB_resh_int$Settings)[i]
  plot = 
    ggplot(as.data.frame(Ex_342AB_resh_int[Ex_342AB_resh_int$Settings == setting, ]), aes(fill = int_type, y=value, x=(Isobar)))+
    geom_bar(stat = "identity", position = "stack")+
    scale_fill_manual(values= color_342AB_int)+
    ylim(0, 1.74)+
    labs(x = setting)+
    scale_x_discrete(guide = guide_axis(angle = 45))+
    theme(text = element_text(size = 13), 
          #axis.text.x=element_text(angle = 45, vjust = 1, hjust=0.5),
          axis.text = element_blank(),
          axis.ticks.x=element_blank(),
          axis.ticks.y=element_blank(),
          panel.background = element_blank(),
          axis.title = element_blank(),
          legend.position = "none"
    )#+
    #coord_flip()+
    #theme(#axis.title.y = element_text(angle = 0, hjust = 0.5, vjust = 0.5), 
      #axis.title.x = element_text(angle = 45, hjust = 1, vjust = 1))
  plots_set_342AB_ex[[i]] = plot
}
i = 1
setting = unique(Ex_342AB_resh_int$Settings)[i]

plots_set_342AB_ex[[i]] <- 
  ggplot(as.data.frame(Ex_342AB_resh_int[Ex_342AB_resh_int$Settings == setting, ]), aes(fill = int_type, y=value, x=(Isobar)))+
  geom_bar(stat = "identity", position = "stack")+
  scale_fill_manual(values= color_342AB_int)+
  ylim(0, 1.74)+
  labs(x = setting)+
  scale_x_discrete(guide = guide_axis(angle = 45))+
  theme(text = element_text(size = 13), 
        axis.text.x=element_blank(),
        #axis.text = element_blank(),
        axis.ticks.x=element_blank(),
        #axis.ticks.y=element_blank(),
        panel.background = element_blank(),
        axis.title = element_blank(),
        legend.position = "none"
  )#+
  # #coord_flip()+
  # theme(#axis.title.y = element_text(angle = 0, hjust = 0.5, vjust = 0.5), 
  #   axis.title.x = element_text(angle = 45, hjust = 1, vjust = 1))


plots_comb_342AB_int_Ex_noset <- (plots_set_342AB_ex[[1]]|plots_set_342AB_ex[[2]]|plots_set_342AB_ex[[3]]|plots_set_342AB_ex[[4]]|plots_set_342AB_ex[[7]]|
                                    plots_set_342AB_ex[[11]]|plots_set_342AB_ex[[6]]|plots_set_342AB_ex[[9]]|
                                    plots_set_342AB_ex[[13]]|plots_set_342AB_ex[[17]]|plots_set_342AB_ex[[5]]|plots_set_342AB_ex[[8]]|
                                    plots_set_342AB_ex[[16]]|plots_set_342AB_ex[[15]]|plots_set_342AB_ex[[10]]|
                                    plots_set_342AB_ex[[12]]|plots_set_342AB_ex[[14]] + theme(legend.position = "none")) + 
  plot_layout(guides = "collect") & theme(axis.title.x = element_blank(),
                                          plot.margin = 
                                            margin(t=0.1, l=0, r=0, b=0, unit = "cm"), 
                                          legend.title = element_blank(),
                                          legend.position = "none")


#342B+F Elite

El_342BF_resh_int <- El_342BF %>%
  select(settings, Isobar, score, Settings, ends_with("_int")) %>%
  pivot_longer(cols = ends_with("_int"), names_to = "type", names_pattern = "^(.*)_int$", values_to = "value")%>%
  mutate(int_type = interaction(Isobar, type))


int_type_lvl <- c("342B.incorr", "342F.incorr", "342B.unassign", "342F.unassign", "342B.assign", "342F.assign")
El_342BF_resh_int$int_type <- factor(El_342BF_resh_int$int_type, levels = int_type_lvl, 
                                     labels = c("342B incorrectly assigned intensity", "342F incorrectly assigned intensity",
                                                "342B unassigned intensity", "342F unassigned intensity",
                                                "342B correctly assigned intensity", "342F correctly assigned intensity"))

plots_set_342BF <- as.list(c())

for (i in 1:(length(unique(El_342BF_resh_int$Settings)))){
  setting = unique(El_342BF_resh_int$Settings)[i]
  plot = 
    ggplot(as.data.frame(El_342BF_resh_int[El_342BF_resh_int$Settings == setting, ]), aes(fill = int_type, y=value, x=(Isobar)))+
    geom_bar(stat = "identity", position = "stack")+
    scale_fill_manual(values= color_342BF_int)+
    ylim(0, 1.7)+
    labs(x = setting)+
    scale_x_discrete(guide = guide_axis(angle = 45))+
    theme(text = element_text(size = 12), 
          #axis.text.x=element_text(angle = 45, vjust = 1, hjust=0.5),
          axis.text = element_blank(),
          axis.ticks.x=element_blank(),
          axis.ticks.y=element_blank(),
          panel.background = element_blank(),
          axis.title = element_blank(),
          legend.position = "none"
    )
  #+
  #coord_flip()+
  #  theme(#axis.title.y = element_text(angle = 0, hjust = 0.5, vjust = 0.5), 
  #    axis.title.x =  element_text(angle = 45, hjust = 1, vjust = 1))
  plots_set_342BF[[i]] = plot
}

plots_comb_342BF_int_El_noleg <- (plots_set_342BF[[1]]|plots_set_342BF[[2]]|plots_set_342BF[[3]]|plots_set_342BF[[4]]|plots_set_342BF[[7]]|
                                    plots_set_342BF[[11]]|plots_set_342BF[[6]]|plots_set_342BF[[9]]|
                                    plots_set_342BF[[13]]|plots_set_342BF[[5]]|plots_set_342BF[[8]]|
                                    plots_set_342BF[[15]]|plots_set_342BF[[16]]|plots_set_342BF[[10]]|
                                    plots_set_342BF[[12]]|plots_set_342BF[[14]] + 
                                    theme(legend.position = "none")) + plot_layout(guides = "collect") & 
  theme(axis.title.x = element_blank(),
        plot.margin = 
          margin(t=0.1, l=0, r=0, b=0, unit = "cm"), 
        legend.position = "none")

# 342B+F Exploris

Ex_342BF_resh_int <- Ex_342BF %>%
  select(settings, Isobar, score, Settings, ends_with("_int")) %>%
  pivot_longer(cols = ends_with("_int"), names_to = "type", names_pattern = "^(.*)_int$", values_to = "value")%>%
  mutate(int_type = interaction(Isobar, type))

int_type_lvl <- c("342B.incorr", "342F.incorr", "342B.unassign", "342F.unassign", "342B.assign", "342F.assign")
Ex_342BF_resh_int$int_type <- factor(Ex_342BF_resh_int$int_type, levels = int_type_lvl, 
                                     labels = c("342B incorrectly assigned intensity", "342F incorrectly assigned intensity",
                                                "342B unassigned intensity", "342F unassigned intensity",
                                                "342B correctly assigned intensity", "342F correctly assigned intensity"))

plots_set_342BF_ex <- as.list(c())

for (i in 1:(length(unique(Ex_342BF_resh_int$Settings)))){
  setting = unique(Ex_342BF_resh_int$Settings)[i]
  plot = 
    ggplot(as.data.frame(Ex_342BF_resh_int[Ex_342BF_resh_int$Settings == setting, ]), aes(fill = int_type, y=value, x=(Isobar)))+
    geom_bar(stat = "identity", position = "stack")+
    scale_fill_manual(values= color_342BF_int)+
    ylim(0, 1.7)+
    labs(x = setting)+
    scale_x_discrete(guide = guide_axis(angle = 45))+
    theme(text = element_text(size = 12), 
          #axis.text.x=element_text(angle = 45, vjust = 1, hjust=0.5),
          axis.text = element_blank(),
          axis.ticks.x=element_blank(),
          axis.ticks.y=element_blank(),
          panel.background = element_blank(),
          axis.title = element_blank(),
          legend.position = "none"
    )
  #+
  #coord_flip()+
  #  theme(#axis.title.y = element_text(angle = 0, hjust = 0.5, vjust = 0.5), 
  #    axis.title.x =  element_text(angle = 45, hjust = 1, vjust = 1))
  plots_set_342BF_ex[[i]] = plot
}

plots_comb_342BF_int_Ex_noset <- (plots_set_342BF_ex[[1]]|plots_set_342BF_ex[[2]]|plots_set_342BF_ex[[3]]|plots_set_342BF_ex[[4]]|plots_set_342BF_ex[[7]]|
                                    plots_set_342BF_ex[[11]]|plots_set_342BF_ex[[6]]|
                                    plots_set_342BF_ex[[9]]|plots_set_342BF_ex[[13]]|plots_set_342BF_ex[[17]]|plots_set_342BF_ex[[5]]|
                                    plots_set_342BF_ex[[8]]|plots_set_342BF_ex[[16]]|plots_set_342BF_ex[[15]]|
                                    plots_set_342BF_ex[[10]]|plots_set_342BF_ex[[12]]|plots_set_342BF_ex[[14]] + theme(legend.position = "none")) + 
  plot_layout(guides = "collect") & theme(axis.title.x = element_blank(),
                                          plot.margin = 
                                            margin(t=0.1, l=0, r=0, b=0, unit = "cm"))

#### COMBINING INTENSITY PLOTS ###################

Figure5 <- (plots_comb_342AB_int_Ex_noset/plots_comb_342AB_int_El)|(plots_comb_180EG_int_Ex_noset/plots_comb_180EG_int_El_noset)|(plots_comb_342BF_int_Ex_noset / plots_comb_342BF_int_El_noleg)

###############################################################################################################################
# FIGURE 6

data_combined <- rbind(El_180EG, El_342AB, El_342BF, Ex_180EG, Ex_342AB, Ex_342BF) %>%
  group_by(Instrument, Settings)%>%
  mutate(average_sim_score = mean(score),
         average_time = mean(time),
         Settings = factor(Settings, levels = levels),
         Instrument = if_else(Instrument == "Elite", "LIT-Orbitrap", "Q-Orbitrap"))

Figure_6 <- ggplot(data_combined, aes(x = Settings, y = average_time, group = Instrument))+
  geom_point(aes(color = average_sim_score, shape = Instrument), size = 6)+
  geom_line(linetype = 2)+
  theme(text = element_text(size = 14), 
        panel.background = element_blank())+
  scale_x_discrete(guide = guide_axis(angle = 45))+
  scale_colour_gradientn(colours = c('#FD7446', '#332288'),
                         rescaler = ~ scales::rescale_mid(.x, mid = 0.7),
                         name = "Average similarity score")+
  ylab("Time, s")+
  ylim(0, 800)