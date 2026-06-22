data <- read.csv("/Users/marikaclark/dark_matter_atx/data/Lifespan_Study_selected_data (1).csv")
library(tidyverse)
head(data)
colnames(data)

sample_ids <- data %>%
  select(Sample, Cell_line, DNAmethylation_sampleID, RNAseq_sampleID, mtDNA_sequencing_runBarcode, WGS_sampleID)

sample_ids_clean <- sample_ids %>%
  mutate(across(where(is.character), ~ na_if(., "")))

head(sample_ids_clean)
head(data)
subject_pairing <- sample_ids_clean %>%
  left_join(data %>% select(Sample), by = "Sample") %>%
  group_by(Cell_line) %>%
  summarise(
    has_RNA = any(!is.na(RNAseq_sampleID)),
    has_WGS = any(!is.na(WGS_sampleID)),
    has_mtDNA = any(!is.na(mtDNA_sequencing_runBarcode)),
    has_DNAmethyl = any(!is.na(DNAmethylation_sampleID)),
    n_RNA_timepoints = sum(!is.na(RNAseq_sampleID)),
    n_WGS_timepoints = sum(!is.na(WGS_sampleID)),
    .groups = "drop"
  ) %>%
  mutate(
    sex = if_else(Cell_line %in% c("HC2", "HC6", "SURF1_3"), "Female", "Male")
  )

subject_pairing <- subject_pairing[-nrow(subject_pairing),]
write.csv(subject_pairing, "/Users/marikaclark/dark_matter_atx/data_availability/PICARD_data_availability.csv")

age_range <- data %>%
  group_by(Cell_line)%>%
  distinct(data$Donor_age)%>%
  ungroup()
age_range
View(data)
