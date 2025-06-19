setwd("C:/Users/aivanova/Documents/Orbitrap Data/settings comparison")
library(OrgMassSpecR) ## for spectra similarity calculations
library(dplyr) ##for fuzzy join via mz as numeric
library(fuzzyjoin) ## for fuzzy join
library(stringr) ##to extract settings name
source("processing.R") #<- for functions

#Identify comparison table to sum up results - done once for a pair of isobars

file_names <-list.files("./342_settings_comparison/342B+F/", pattern = ".mzML", full.names = TRUE)
file_names <- file_names[-c(grep(file_names, pattern = "MS1"))]

##load the table in progress
comparison_table_342BF <- read.csv("./342_settings_comparison/342B+F/342BF_comparison_table_inprog.csv")

##Which file are you working on?
##1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17
iter <- 15 #<- TO CHANGE!!! 
comparison_table_342BF$settings[iter] <- str_match(file_names[iter], "342B\\+F_(.*?)\\.mzML")[,2]
str_match(file_names[iter], "342B\\+F_(.*?)\\.mzML")[,2]

#0. Load reference spectra and separate by settings - 
#done once for a pair of isobars

int_thres <- 0.003 #relative intensity threshold for imported reference spectra

ref_342B <- read.csv("./342_settings_comparison/342B_ref.csv")
ref_342B_defNCE <- ref_342B %>%
  filter(NCE == 45) %>%
  filter(relative.intensity > int_thres)

ref_342B_35NCE <- ref_342B %>%
  filter(NCE == 35) %>%
  filter(relative.intensity > int_thres)

ref_342B_50NCE <- ref_342B %>%
  filter(NCE == 50) %>%
  filter(relative.intensity > int_thres)

ref_342B_65NCE <- ref_342B %>%
  filter(NCE == 65) %>%
  filter(relative.intensity > int_thres)

ref_342F <- read.csv("./342_settings_comparison/342F_ref_cor.csv")
ref_342F_defNCE <- ref_342F %>%
  filter(NCE == 45) %>%
  filter(relative_intensity > int_thres)

ref_342F_35NCE <- ref_342F %>%
  filter(NCE == 35) %>%
  filter(relative_intensity > int_thres)

ref_342F_50NCE <- ref_342F %>%
  filter(NCE == 50) %>%
  filter(relative_intensity > int_thres)

ref_342F_65NCE <- ref_342F %>%
  filter(NCE == 65) %>%
  filter(relative_intensity > int_thres)

# Get reference spectra that contains only theoretical mz and intensities

# !!! change NCE if necessary - for 3 data files with NCE  35, 50 and 65
ref_342B_defNCE_int <- ref_342B_defNCE %>%
  select(mz = mz_theoretical, int_ref = relative.intensity)
ref_342F_defNCE_int <- ref_342F_defNCE %>%
  select(mz = mz_theoretical, int_ref = relative_intensity)


# 1. Load and process shifting window data
dev <- 3 ##starting mass deviation in ppm to combine between scans'

msExp <- openMSfile(file_names[iter])

# Find all centroided MS2 scan number
ms2Scan <- c(1:81) #1:41 for 0.04 Da step, 1:17 for 0.1 Da step and 1:137 for 2 Da wide window


# Find the scan with maximum TIC and use it as the merging center.
TIC_max <- which.max(sapply(ms2Scan, function(x) {header(msExp, x)$totIonCurrent}))

# Only take peaks with at least 0.003 relative intensity (relative to TIC)
int_thres_mix <- 0.003
comparison_table_342BF$int_thres_mix[iter] <- int_thres_mix
filtered_spectra <- getSpectrumLists(msExp, ms2Scan, int_thres_mix)

# The template chimeric MS2 is an averaged spectrum from scan [TIC_max-5, TIC_max+5]
# Mass deviation for merging: 6ppm
template_chimera <- averageSpectra(filtered_spectra, (TIC_max-5):(TIC_max+5), dev*2)


# 2. Load and process blank
blankExp <- openMSfile("./342_settings_comparison/blank342_stepped_NCE_3_1.mzML") # <= for stepped NCE
blankExp <- openMSfile("./342_settings_comparison/blank342_separate_NCE_3_2.mzML") # <= for separate NCEs

# Find all centroided MS2 scan number 
ms2Scan_blank <- getMS2Scan(blankExp)

#for different NCEs than default stepped NCE:
NCE = 65 #35 50 OR 65
ms2Scan_blank <- which(header(blankExp)$collisionEnergy == NCE)


# Get an averaged blank spectrum, only merge peaks with at least 0.01 relative intensity
# Mass deviation for merging: 6ppm
avg_blank <- averageSpectra(getSpectrumLists(blankExp, ms2Scan_blank, 0.01), ms2Scan_blank, dev*2)
# Remove blank peaks from template chimeric MS2
# Mass deviation for peak matching: 3ppm
template_chimera <- removeBlank_2(template_chimera, avg_blank, dev)
#remove rows and columns for [M+1+H]+
template_chimera <- template_chimera%>%
  as.data.frame()%>%
  filter(mz < 344)
# Get the normalized intensities for peaks (middle of run)
int_norm <- as.data.frame(template_chimera) %>%
  mutate(norm_int = intensity/max(intensity)) 


# 3. Get intensity matrix and correlation matrix
int_matrix <- getIntensityMatrix(template_chimera[, 1], getSpectrumLists(msExp, ms2Scan, thres=NULL), dev)
corr_matrix <- getCorrMatrix(int_matrix)

# 4. Take from correlation matrix the correlation between precursors and fragments only, retain only positive correlations
corr_matrix_prec <- as.data.frame(corr_matrix[ ,c(ncol(corr_matrix) - 1, ncol(corr_matrix))]) %>%
  filter(if_any(everything(), ~ . > 0))

## For NCE 50, use the highest intensity fragment as a proxy for the precursor. 342B fragment: 108.0570 (column 8), 342F fragment: 207.07642 (column 43)
#corr_matrix_prec <-  as.data.frame(corr_matrix[ ,c(ncol(corr_matrix), 44)]) %>%
#filter(if_any(everything(), ~ . > 0))

## for NCE 65,  use the highest intensity fragment as a proxy for the precursor. 342B fragment: 108.0570 (column 8), 342F fragment: 179.05769 (col 58)
#corr_matrix_prec <- as.data.frame(corr_matrix[ ,c(16, 58)])%>%
#  filter(if_any(everything(), ~ . > 0))

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
  difference_left_join(int_norm, by = "mz", max_dist = 1e-4, distance_col = NULL) %>%
  select(mz = mz.x, norm_int, cor_val)%>%
  mutate(norm_int = norm_int/max(norm_int))%>%
  filter(norm_int > 0.003)

##compare reference and reconstructed spectra - join with a large mass error to calculate mz errors between measured and theoretical mz
##for isobar 1
comp_342B <- ref_342B_defNCE_int %>%
  difference_full_join(reconstr_spec_1, by = "mz", max_dist = 0.003, distance_col = NULL)%>%
  select(mz_theor = mz.x, mz_measur = mz.y, int_reconstr = norm_int, int_ref, cor_val)%>%
  mutate(mz_err_Da = mz_measur-mz_theor)

#optionally - plot absolute error in Da to check if it has clear quadratic pattern
plot(x = comp_342B$mz_measur, y = comp_342B$mz_err_Da)

#recalibration using quadratic model
model_1 <- lm(comp_342B$mz_err_Da ~(comp_342B$mz_measur+ I(comp_342B$mz_measur^2)))

comp_342B <- comp_342B %>%
  mutate(mz_corr = mz_measur - predict(model_1),
         mz_err_corr = (mz_corr-mz_theor)/mz_corr*1e6,
         mz_cor_err_Da = mz_corr-mz_theor)

#optionally - check destriburtion of mz errors after recalibration
plot(x = comp_342B$mz_measur, y = comp_342B$mz_err_corr)
## save the recalibrated reconstructed spectra
reconstr_spec_1_t <- as.data.frame(na.omit(select(comp_342B, mz_corr, int_reconstr)))

#for isobar 2
comp_342F <-ref_342F_defNCE_int %>%
  difference_full_join(reconstr_spec_2, by = "mz", max_dist = 0.003, distance_col = NULL)%>%  
  select(mz_theor = mz.x, mz_measur = mz.y, int_reconstr = norm_int, int_ref, cor_val)%>%
  mutate(mz_err_Da = mz_measur-mz_theor)
plot(x = comp_342F$mz_measur, y = comp_342F$mz_err_Da)

##change the mz value fix model type: linear, log, or quadratic
model_2 <- lm(comp_342F$mz_err_Da ~ (comp_342F$mz_measur+ I(comp_342F$mz_measur^2)))

comp_342F <- comp_342F %>%
  mutate(mz_corr = mz_measur - predict(model_2),
         mz_err_corr = (mz_corr-mz_theor)/mz_corr*1e6,
         mz_cor_err_Da = mz_corr-mz_theor)
#optionally - check distriburtion of mz errors after recalibration
plot(x = comp_342F$mz_measur, y = comp_342F$mz_err_corr)
reconstr_spec_2_t <- as.data.frame(na.omit(select(comp_342F, mz_corr, int_reconstr)))

##fill out comparison table
###get the run time
comparison_table_342BF$time[iter] <- as.numeric(header(msExp)$retentionTime[length(ms2Scan)] - header(msExp)$retentionTime[1])

##SAVING RESULTS
##fill out comparison table  for isobar 1

#number of cor assigned peaks
comparison_table_342BF$Isobar_1_corr_peaks[iter] <-  sum(!is.na(comp_342B$int_reconstr[!is.na(comp_342B$mz_theor)]))
#number of peaks missing in reconstr spectra
comparison_table_342BF$Isobar_1_unassign_peaks[iter] <- 
  length(comp_342B$int_ref) - sum(!is.na(comp_342B$int_reconstr) & !is.na(comp_342B$mz_theor))
#sum int of correctly assigned peaks compare to ref intensity = recall
comparison_table_342BF$Isobar_1_assign_int[iter] <-
  sum(comp_342B$int_ref[!is.na(comp_342B$int_reconstr) & !is.na (comp_342B$int_ref)])/sum(comp_342B$int_ref, na.rm = TRUE)
#number of incorrectly assigned peaks
comparison_table_342BF$Isobar_1_incorr_peaks[iter] <- sum(is.na(comp_342B$mz_theor), na.rm = TRUE)
#sum intensity of incorrectly assigned peaks in ref spectra
comparison_table_342BF$Isobar_1_incorr_int[iter] <-  
  sum(comp_342B$int_reconstr[is.na(comp_342B$mz_theor)])/sum(ref_342B_defNCE_int$int_ref)
#spectrum similarity 
comparison_table_342BF$Isobar_1_score[iter] <- 
  SpectrumSimilarity(reconstr_spec_1_t, ref_342B_defNCE_int, t = 0.0003, b = 0, xlim = c(50, 350))
#average errors before recalibration
comparison_table_342BF$Isobar_1_aver_err_mDa[iter] <- mean(comp_342B$mz_err_Da, na.rm = TRUE)*1000
comparison_table_342BF$Isobar_1_aver_err_mDa_sd[iter] <- sd(comp_342B$mz_err_Da, na.rm = TRUE)*1000
comparison_table_342BF$Isobar_1_aver_err_ppm[iter] <- mean((comp_342B$mz_err_Da/comp_342B$mz_theor*10^6), na.rm = TRUE)
comparison_table_342BF$Isobar_1_aver_err_ppm_sd[iter] <- sd((comp_342B$mz_err_Da/comp_342B$mz_theor*10^6), na.rm = TRUE)
#average correl value between prec and fragments
comparison_table_342BF$Isobar_1_cor_aver[iter] <- mean(comp_342B$cor_val, na.rm = TRUE)
#precision
comparison_table_342AB$Isobar_1_precision[iter] <-  
  sum(comp_342B$int_reconstr[!is.na(comp_342B$int_reconstr) & !is.na (comp_342B$int_ref)])/sum(comp_342B$int_reconstr[!is.na(comp_342B$int_reconstr)])



##fill put comparison table for isobar 2 - same parameters
comparison_table_342BF$Isobar_2_corr_peaks[iter] <-  sum(!is.na(comp_342F$int_reconstr[!is.na(comp_342F$mz_theor)]))
comparison_table_342BF$Isobar_2_unassign_peaks[iter] <- length(comp_342F$int_ref) - sum(!is.na(comp_342F$int_reconstr))
comparison_table_342BF$Isobar_2_assign_int[iter] <- sum(comp_342F$int_ref[!is.na(comp_342F$int_reconstr) & !is.na (comp_342F$int_ref)])/sum(comp_342F$int_ref, na.rm = TRUE)
comparison_table_342BF$Isobar_2_incorr_peaks[iter] <- sum(is.na(comp_342F$mz_theor), na.rm = TRUE)
comparison_table_342BF$Isobar_2_incorr_int[iter] <- sum(comp_342F$int_reconstr[is.na(comp_342F$mz_theor)])/sum(ref_342F_defNCE_int$int_ref)
comparison_table_342BF$Isobar_2_score[iter] <- SpectrumSimilarity(reconstr_spec_2_t, ref_342F_defNCE_int, t = 0.0003, b = 0, xlim = c(50, 350))
comparison_table_342BF$Isobar_2_aver_err_mDa[iter] <- mean(comp_342F$mz_err_Da, na.rm = TRUE)*1000
comparison_table_342BF$Isobar_2_aver_err_mDa_sd[iter] <- sd(comp_342F$mz_err_Da, na.rm = TRUE)*1000
comparison_table_342BF$Isobar_2_aver_err_ppm[iter] <- mean((comp_342F$mz_err_Da/comp_342F$mz_theor*10^6), na.rm = TRUE)
comparison_table_342BF$Isobar_2_aver_err_ppm_sd[iter] <- sd((comp_342F$mz_err_Da/comp_342F$mz_theor*10^6), na.rm = TRUE)
comparison_table_342BF$Isobar_2_cor_aver[iter] <- mean(comp_342F$cor_val, na.rm = TRUE)
comparison_table_342AB$Isobar_2_precision[iter] <-  
  sum(comp_342F$int_reconstr[!is.na(comp_342F$int_reconstr) & !is.na (comp_342F$int_ref)])/sum(comp_342F$int_reconstr[!is.na(comp_342F$int_reconstr)])


#save results
write.csv(comparison_table_342BF, "./342_settings_comparison/342B+F/342BF_comparison_table_inprog.csv")

##create and save data files for mirror plots 

mirror_plot_data_342B<- bind_rows(
  comp_342B %>%
    transmute(mz = mz_theor, int_norm = -int_ref, type = "Reference"),
  comp_342B %>%
    transmute(mz = mz_corr, int_norm = int_reconstr, type = "Reconstructed"))

write.csv(mirror_plot_data_342B, file = 
            paste0("./342_settings_comparison/342B+F_", str_match(file_names[iter], "342B\\+F_(.*?)\\.mzML")[,2], "_342Bpeakstable.csv" ))

mirror_plot_data_342F <- bind_rows(
  comp_342F %>%
    transmute(mz = mz_theor, int_norm = -int_ref, type = "Reference"),
  comp_342F %>%
    transmute(mz = mz_corr, int_norm = int_reconstr, type = "Reconstructed"))

write.csv(mirror_plot_data_342F, file = 
            paste0("./342_settings_comparison/342B+F_", str_match(file_names[iter], "342B\\+F_(.*?)\\.mzML")[,2], "_342Fpeakstable.csv" ))

m_p_B <- ggplot(mirror_plot_data_342B, aes(x=mz, y=int_norm)) +
  geom_segment(aes(x=mz, xend=mz, y=0, yend=int_norm, color = type), linewidth = 1) +
  xlab("m/z") +
  ylab("Relative intensity") +
  annotate("text", y = c(1.1, -1.1), x = 200, label = c("342B: Reconstructed", "342B: Reference"), color = c("#44AA99", "#97b1ab"), size = 5)+
  geom_hline(yintercept=0, linewidth = 0.5)+
  geom_hline(yintercept=0, linewidth = 0.5)+
  geom_vline(xintercept = 50, linewidth = 0.5)+
  theme(axis.title = element_text(size=12),
        axis.text = element_text(size=12),
        panel.grid=element_blank(),
        panel.background = element_blank(),
        legend.position = 'none')+
  scale_color_manual(values= c("#44AA99", "#97b1ab"))+
  xlim(50, 350)
color_342BF <- c("#44AA99", "#332288")

m_p_F <- ggplot(mirror_plot_data_342F, aes(x=mz, y=int_norm)) +
  geom_segment(aes(x=mz, xend=mz, y=0, yend=int_norm, color = type), linewidth = 1) +
  xlab("m/z") +
  ylab("Relative intensity") +
  annotate("text", y = c(1.1, -1.1), x = 200, label = c("342F: Reconstructed", "342F: Reference"), color = c("#332288", "#afa8ba"), size = 5)+
  geom_hline(yintercept=0, linewidth = 0.5)+
  theme(axis.title.x = element_text(size=12),
        axis.title.y = element_blank(),
        axis.text.x = element_text(size=12),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        panel.grid=element_blank(),
        panel.background = element_blank(),
        legend.position = 'none')+
  scale_color_manual(values= c("#332288", "#afa8ba"))+
  xlim(50, 350)
m_p_342BF_El <- m_p_B+m_p_F
m_p_342BF_El 
save(m_p_342BF_El, file = paste0("./342_settings_comparison/342B+F_", str_match(file_names[iter], "342B\\+F_(.*?)\\.mzML")[,2], "_342BFmirror.RData" ))
##Iteration is finished!
