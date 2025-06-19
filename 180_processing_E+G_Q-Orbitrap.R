setwd("C:/Users/aivanova/Documents/Orbitrap Data/settings comparison/Orbi_Exploris")
library(OrgMassSpecR) ## for spectra similarity calculations
library(dplyr) ##for fuzzy join via mz as numeric
library(fuzzyjoin) ## for fuzzy join
library(stringr) ##to extract settings name
source("processing.R")

#Identify comparison table to sum up results - done once for a pair of isobars
file_names <-list.files("./180_settings_comparison/180E+G/", pattern = ".mzML", full.names = TRUE)

##load the table in progress
comparison_table_180EG <- read.csv("./180_settings_comparison/180EG_comparison_table_inprog.csv")

##Which file are you working on?
##1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17
iter <- 8 #<- TO CHANGE!!! 
comparison_table_180EG$settings[iter] <- str_match(file_names[iter], "180E\\+G_([^_]+)_")[,2]
str_match(file_names[iter], "180E\\+G_([^_]+)_")[,2]

#0. Load reference spectra and separate by settings - 
#done once for a pair of isobars

int_thres <- 0.003 #relative intensity threshold for imported reference spectra

ref_180E <- read.csv("./180_settings_comparison/180E_ref.csv")
ref_180E_defNCE <- ref_180E %>%
  filter(NCE == 45) %>%
  filter(norm_int > int_thres)

ref_180E_35NCE <- ref_180E %>%
  filter(NCE == 35) %>%
  filter(norm_int  > int_thres)

ref_180E_50NCE <- ref_180E %>%
  filter(NCE == 50) %>%
  filter(norm_int  > int_thres)

ref_180E_65NCE <- ref_180E %>%
  filter(NCE == 65) %>%
  filter(norm_int  > int_thres)

ref_180G <- read.csv("./180_settings_comparison/180G_ref.csv")
ref_180G_defNCE <- ref_180G %>%
  filter(NCE == 45) %>%
  filter(norm_int > int_thres)

ref_180G_35NCE <- ref_180G %>%
  filter(NCE == 35) %>%
  filter(norm_int > int_thres)

ref_180G_50NCE <- ref_180G %>%
  filter(NCE == 50) %>%
  filter(norm_int > int_thres)

ref_180G_65NCE <- ref_180G %>%
  filter(NCE == 65) %>%
  filter(norm_int > int_thres)

# Get reference spectra that contains only theoretical mz and intensities

# !!! change NCE if necessary - for 3 data files with NCE  35, 50 and 65
ref_180E_defNCE_int <- ref_180E_defNCE %>%
  select(mz = mz_theoretical, int_ref = norm_int)
ref_180G_defNCE_int <- ref_180G_defNCE %>%
  select(mz = mz_theoretical, int_ref = norm_int)


# 1. Load and process shifting window data
dev <- 2 ##starting mass deviation in ppm to combine between scans'

msExp <- openMSfile(file_names[iter])

# Find all centroided MS2 scan number / set it yourself
ms2Scan <- c(1:137) # 82-162 for default, 42-82 for 0.04 Da step, 35-51 for 0.1 Da step nd 1-137 for 2 Da wide window
#ms2Scan <- getMS2Scan(msExp)

# Find the scan with maximum TIC and use it as the merging center.
TIC_max <- which.max(sapply(ms2Scan, function(x) {header(msExp, x)$totIonCurrent}))

# Only take peaks with at least 0.003 relative intensity (relative to TIC)
int_thres_mix <- 0.003
comparison_table_180EG$int_thres_mix[iter] <- int_thres_mix
filtered_spectra <- getSpectrumLists(msExp, ms2Scan, int_thres_mix)

# The template chimeric MS2 is an averaged spectrum from scan [TIC_max-5, TIC_max+5]
# Mass deviation for merging: 6ppm
template_chimera <- averageSpectra(filtered_spectra, (TIC_max - 5):(TIC_max + 5), dev*2)


# 2. Load and process blank
blankExp <- openMSfile("./blnk_180_separNCE_04.mzML") #for separate NCE
blankExp <- openMSfile("./180_settings_comparison/blank_180_stepNCE_03.mzML") #for stepped NCE
# Find all centroided MS2 scan number
ms2Scan_blank <- getMS2Scan(blankExp)

#for different NCEs than default:
NCE = 65 #35 50 OR 65
ms2Scan_blank <- which(header(blankExp)$collisionEnergy == NCE)

# Get an averaged blank spectrum, only merge peaks with at least 0.005 relative intensity
# Mass deviation for merging: 6ppm
avg_blank <- averageSpectra(getSpectrumLists(blankExp, ms2Scan_blank, 0.01), ms2Scan_blank, dev*2)
# Remove blank peaks from template chimeric MS2
# Mass deviation for peak matching: 3ppm
template_chimera <- removeBlank_2(template_chimera, avg_blank, dev)
#remove rows and columns for [M+1+H]+
template_chimera <- template_chimera%>%
  as.data.frame()%>%
  filter(mz < 181.8)
# Get the normalized intensities for peaks (middle of run)
int_norm <- as.data.frame(template_chimera) %>%
  mutate(norm_int = intensity/max(intensity)) 


# 3. Get intensity matrix and correlation matrix
int_matrix <- getIntensityMatrix(template_chimera[, 1], getSpectrumLists(msExp, ms2Scan, thres=NULL), dev)
corr_matrix <- getCorrMatrix(int_matrix)

# 4. Take from correlation matrix the correlation between precursors and fragments only
corr_matrix_prec <- as.data.frame(corr_matrix[ ,c(ncol(corr_matrix) - 1, ncol(corr_matrix))]) %>%
  filter(if_any(everything(), ~ . > 0))

##For NCE 50, use the highest intensity fragment as a proxy for the precursors, fragment: 180E fragment:102.0549 (column 9), 180G fragment: 135.0441 (col 12)
corr_matrix_prec <- corr_matrix[ ,c(9, 12)]%>%
  as.data.frame()%>%
  filter(if_any(everything(), ~ . > 0))

## for NCE 65,  use the highest intensity fragment as a proxy for the precursors. 180E fragment: 56.049 (column 1), 180G fragment: 135.0441 (col 12)
corr_matrix_prec <- corr_matrix[ ,c(1, 12)]%>%
  as.data.frame()%>%
  filter(if_any(everything(), ~ . > 0))

# 5. Filter to which precursor the fragment belongs better - shared fragments are assigned to one precursor only
reconstr_spec_1 <- data.frame(matrix(ncol = 2, nrow = 0))
reconstr_spec_2 <- data.frame(matrix(ncol = 2, nrow = 0))

for (i in 1:nrow(corr_matrix_prec)){
  if(corr_matrix_prec[i, 1] > corr_matrix_prec[i, 2]){
    reconstr_spec_1[i, 1] = as.numeric(rownames(corr_matrix_prec)[i])
    reconstr_spec_1[i, 2] = as.numeric((corr_matrix_prec[,1])[i])
  } else {
    #reconstr_spec_2 = rbind(reconstr_spec_2, as.numeric(rownames(corr_matrix_prec)[i]))
    reconstr_spec_2[i, 1] = as.numeric(rownames(corr_matrix_prec)[i])
    reconstr_spec_2[i, 2] = as.numeric((corr_matrix_prec[,2])[i])
  }
}

# Get and normalize intensities for picked peaks (retrieved from template with averaged 10 scans around the highest TICs)
reconstr_spec_1 <- reconstr_spec_1 %>%
  na.omit()%>%
  setNames(c("mz", "cor_val")) %>%
  difference_left_join(int_norm, by = "mz", max_dist = 1e-4, distance_col = NULL) %>%
  select(mz = mz.y, norm_int, cor_val)%>%
  mutate(norm_int = norm_int/max(norm_int))%>%
  filter(norm_int > 0.003)


reconstr_spec_2 <- reconstr_spec_2%>%
  na.omit()%>%
  setNames(c("mz", "cor_val")) %>%
  difference_left_join(int_norm, by = "mz", max_dist = 2e-4, distance_col = NULL) %>%
  select(mz = mz.x, norm_int, cor_val)%>%
  mutate(norm_int = norm_int/max(norm_int))%>%
  filter(norm_int > 0.003)

##compare reference and reconstructed spectra - join with a large mass error to calculate mz errors between measured and theoretical mz
##for isobar 1
comp_180E <- ref_180E_defNCE_int %>%
  difference_full_join(reconstr_spec_1, by = "mz", max_dist = 0.003, distance_col = NULL)%>%
  select(mz_theor = mz.x, mz_measur = mz.y, int_reconstr = norm_int, int_ref, cor_val)%>%
  mutate(mz_err_Da = mz_measur-mz_theor)
#optionally - plot absolute error in Da to check if it has clear log pattern
plot(x = comp_180E$mz_measur, y = comp_180E$mz_err_Da)

model_1 <- lm(comp_180E$mz_err_Da ~ log(comp_180E$mz_measur))
summary(model_1)
comp_180E <- comp_180E %>%
  mutate(mz_corr = mz_measur - predict(model_1),
         mz_err_corr = (mz_corr-mz_theor)/mz_corr*1e6,
         mz_cor_err_Da = mz_corr-mz_theor)

#optionally - check destriburtion of mz errors after recalibration
plot(x = comp_180E$mz_measur, y = comp_180E$mz_err_corr)
## save the recalibrated reconstructed spectra
reconstr_spec_1_t <- as.data.frame(na.omit(select(comp_180E, mz_corr, int_reconstr)))

#for isobar 2
comp_180G <-ref_180G_defNCE_int %>%
  difference_full_join(reconstr_spec_2, by = "mz", max_dist = 0.003, distance_col = NULL)%>%  
  select(mz_theor = mz.x, mz_measur = mz.y, int_reconstr = norm_int, int_ref, cor_val)%>%
  mutate(mz_err_Da = mz_measur-mz_theor)

plot(x = comp_180G$mz_measur, y = comp_180G$mz_err_Da)

##change the mz value fix model type: linear, log, or quadratic
model_2 <- lm(comp_180G$mz_err_Da ~ comp_180G$mz_measur)

summary(model_2)
comp_180G <- comp_180G %>%
  mutate(mz_corr = mz_measur - predict(model_2),
         mz_err_corr = (mz_corr-mz_theor)/mz_corr*1e6,
         mz_cor_err_Da = mz_corr-mz_theor)
#optionally - check distriburtion of mz errors after recalibration
plot(x = comp_180G$mz_measur, y = comp_180G$mz_err_corr)
reconstr_spec_2_t <- as.data.frame(na.omit(select(comp_180G, mz_corr, int_reconstr)))

##fill out comparison table
###get the run time
comparison_table_180EG$time[iter] <- as.numeric(header(msExp)$retentionTime[length(ms2Scan)] - header(msExp)$retentionTime[1])


##fill out comparison table for isobar 1
#number of cor assigned peaks
comparison_table_180EG$Isobar_1_corr_peaks[iter] <- sum(!is.na(comp_180G$int_reconstr) & !is.na(comp_180G$int_ref))
#number of peaks missing in reconstr spectra
comparison_table_180EG$Isobar_1_unassign_peaks[iter] <- 
  length(comp_180E$int_ref) - sum(!is.na(comp_180E$int_reconstr)& !is.na(comp_180E$mz_theor))
#sum int of correctly assigned peaks compare to ref intensity = recall
comparison_table_180EG$Isobar_1_assign_int[iter] <-
  sum(comp_180E$int_ref[!is.na(comp_180E$int_reconstr) & !is.na (comp_180E$int_ref)])/sum(comp_180E$int_ref, na.rm = TRUE)
#number of incorrectly assigned peaks
comparison_table_180EG$Isobar_1_incorr_peaks[iter] <- sum(is.na(comp_180E$mz_theor), na.rm = TRUE)
#spectrum similarity 
comparison_table_180EG$Isobar_1_score[iter] <- 
  SpectrumSimilarity(reconstr_spec_1_t, ref_180E_defNCE_int, t = 0.00018, b = 0, xlim = c(50, 190))
#average errors before recalibration
comparison_table_180EG$Isobar_1_aver_err_mDa[iter] <- mean(comp_180E$mz_err_Da, na.rm = TRUE)*1000
comparison_table_180EG$Isobar_1_aver_err_mDa_sd[iter] <- sd(comp_180E$mz_err_Da, na.rm = TRUE)*1000
comparison_table_180EG$Isobar_1_aver_err_ppm[iter] <- mean((comp_180E$mz_err_Da/comp_180E$mz_theor*10^6), na.rm = TRUE)
comparison_table_180EG$Isobar_1_aver_err_ppm_sd[iter] <- sd((comp_180E$mz_err_Da/comp_180E$mz_theor*10^6), na.rm = TRUE)
#average correl value between prec and fragments
comparison_table_180EG$Isobar_1_cor_aver[iter] <- mean(comp_180E$cor_val, na.rm = TRUE)
#sum intensity of incorrectly assigned peaks compare to ref spectra
comparison_table_180EG$Isobar_1_incorr_int[iter] <- sum(comp_180E$int_reconstr[is.na(comp_180E$mz_theor)])/sum(ref_180E_defNCE_int$int_ref)
# precision
comparison_table_180EG$Isobar_1_precision[iter] <-  
  sum(comp_180E$int_reconstr[!is.na(comp_180E$int_reconstr) & !is.na (comp_180E$int_ref)])/sum(comp_180E$int_reconstr[!is.na(comp_180E$int_reconstr)])

##fill put comparison table for isobar 2
comparison_table_180EG$Isobar_2_corr_peaks[iter] <- sum(!is.na(comp_180G$int_reconstr) & !is.na(comp_180G$int_ref))
comparison_table_180EG$Isobar_2_unassign_peaks[iter] <- length(comp_180G$int_ref[!is.na(comp_180G$int_ref)]) - sum(!is.na(comp_180G$int_reconstr[!is.na(comp_180G$mz_theor)]))
comparison_table_180EG$Isobar_2_assign_int[iter] <- sum(comp_180G$int_ref[!is.na(comp_180G$int_reconstr) & !is.na (comp_180G$int_ref)])/sum(comp_180G$int_ref, na.rm = TRUE)
comparison_table_180EG$Isobar_2_incorr_peaks[iter] <- sum(is.na(comp_180G$mz_theor), na.rm = TRUE)
comparison_table_180EG$Isobar_2_score[iter] <- SpectrumSimilarity(reconstr_spec_2_t, ref_180G_defNCE_int, t = 0.00018, b = 0, xlim = c(50, 190))
comparison_table_180EG$Isobar_2_aver_err_mDa[iter] <- mean(comp_180G$mz_err_Da, na.rm = TRUE)*1000
comparison_table_180EG$Isobar_2_aver_err_mDa_sd[iter] <- sd(comp_180G$mz_err_Da, na.rm = TRUE)*1000
comparison_table_180EG$Isobar_2_aver_err_ppm[iter] <- mean((comp_180G$mz_err_Da/comp_180G$mz_theor*10^6), na.rm = TRUE)
comparison_table_180EG$Isobar_2_aver_err_ppm_sd[iter] <- sd((comp_180G$mz_err_Da/comp_180G$mz_theor*10^6), na.rm = TRUE)
comparison_table_180EG$Isobar_2_cor_aver[iter] <- mean(comp_180G$cor_val, na.rm = TRUE)
comparison_table_180EG$Isobar_2_incorr_int[iter] <- sum(comp_180G$int_reconstr[is.na(comp_180G$mz_theor)])/sum(int_norm$norm_int)
comparison_table_180EG$Isobar_2_precision[iter] <-  
  sum(comp_180G$int_reconstr[!is.na(comp_180G$int_reconstr) & !is.na (comp_180G$int_ref)])/sum(comp_180G$int_reconstr[!is.na(comp_180G$int_reconstr)])
write.csv(comparison_table_180EG, "C:/Users/aivanova/Documents/Orbitrap Data/settings comparison/Orbi_Exploris/180_settings_comparison/180EG_comparison_table_inprog.csv")
##Iteration is finished!

mirror_plot_data_180E <- bind_rows(
  comp_180E %>%
    transmute(mz = mz_theor, int_norm = -int_ref, type = "Reference"),
  comp_180E %>%
    transmute(mz = mz_corr, int_norm = int_reconstr, type = "Reconstructed"))

write.csv(mirror_plot_data_180E, file = paste0("./180_settings_comparison/180E+G_", 
                                               str_match(file_names[iter], "180E\\+G_(.*?)\\.mzML")[,2], "_180Epeakstable.csv" ))

mirror_plot_data_180G <- bind_rows(
  comp_180G %>%
    transmute(mz = mz_theor, int_norm = -int_ref, type = "Reference"),
  comp_180G %>%
    transmute(mz = mz_corr, int_norm = int_reconstr, type = "Reconstructed"))

write.csv(mirror_plot_data_180G, file = paste0("./342_settings_comparison/180E+G_", 
                                               str_match(file_names[iter], "180E\\+G_(.*?)\\.mzML")[,2], "_180Gpeakstable.csv" ))

color_180 <- c("#882255", "#CC6677")


m_p_E <- ggplot(mirror_plot_data_180E, aes(x=mz, y=int_norm)) +
  geom_segment(aes(x=mz, xend=mz, y=0, yend=int_norm, color = type), linewidth = 1) +
  xlab("m/z") +
  ylab("Relative intensity") +
  annotate("text", y = c(1.1, -1.1), x = 120, label = c("180E: Reconstructed", "180E: Reference"), color = c("#882255", "#bca5ad"), size = 5)+
  geom_hline(yintercept=0, linewidth = 0.5)+
  geom_vline(xintercept = 50, linewidth = 0.5)+
  theme(axis.title = element_text(size=12),
        axis.text = element_text(size=12),
        panel.grid=element_blank(),
        panel.background = element_blank(),
        legend.position = 'none')+
  scale_color_manual(values=c("#882255", "#bca5ad"))+
  xlim(50, 200)

m_p_G <- ggplot(mirror_plot_data_180G, aes(x=mz, y=int_norm)) +
  geom_segment(aes(x=mz, xend=mz, y=0, yend=int_norm, color = type), linewidth = 1) +
  xlab("m/z") +
  ylab("Relative intensity") +
  annotate("text", y = c(1.1, -1.1), x = 120, label = c("180G: Reconstructed", "180G: Reference"), color = c("#CC6677", "#bfa5a7"), size = 5)+
  geom_hline(yintercept=0, linewidth = 0.5)+
  theme(axis.title.x = element_text(size=12),
        axis.title.y = element_blank(),
        axis.text.x = element_text(size=12),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        panel.grid=element_blank(),
        panel.background = element_blank(),
        legend.position = 'none')+
  scale_color_manual(values= c("#CC6677", "#bfa5a7"))+
  xlim(50, 200)

library(patchwork)
m_p_180EG_Ex <- m_p_E+m_p_G
m_p_180EG_Ex 
save(m_p_180EG_Ex, file = paste0("./180_settings_comparison/180E+G_", 
                                 str_match(file_names[iter], "180E\\+G_(.*?)\\.mzML")[,2], "_180EG_mirror.RData" ))
