```r
#############################################################
# Required Packages
#############################################################
library(dplyr)

#setwd("/home/adminsb/pc100_80_42_95/SF_analysis/")
setwd("/home/adminsb/hpc/space2/metalysisAML/FI_SF_analysis_starfusion/analysis_rfus/novel_fusion/")

# Data Import and Case Normalisation (case vs control fusion)
# Remove fusions detected in control samples from the case dataset

cnt <- read.delim("../control77_combined_FI-star.tsv")
cnt_fuz <- data.frame(unique(cnt$X.FusionName))

case <- read.delim("../case801_combined_FI-star.tsv")

ft0 <- case %>% 
  dplyr::filter(!X.FusionName %in% cnt_fuz[[1]]) %>% 
  dplyr::distinct(sample_name, X.FusionName, LeftBreakpoint, RightBreakpoint, .keep_all = TRUE)

write.table(ft0, "case801_combined_FI-star_normalised_FT0.tsv", row.names = F, sep = "\t")


# ============================================================
# DATA FILTERING FOR ANALYSIS OF NOVEL AND KNOWN FUSION SPECTRUM
# ============================================================

# FT1: Initial filtering based on fusion read support, annotation,
# large anchor support and exclusion of common normal/readthrough
# and immunoglobulin, T-cell receptor, HLA, ribosomal and mitochondrial genes

ft1 <- ft0 %>%
  dplyr::filter(
    JunctionReadCount >= 5, SpanningFragCount >= 2, FFPM >= 0.1,
    LargeAnchorSupport == "YES",
    !grepl("readthrough|GTEx|BodyMap|Normal", annots, ignore.case = TRUE),
    !grepl("^IG[HKL]|^TRA|^TRB|^TRD|^TRG|^HLA|^RPL|^RPS|^MT-|^RNU|^RNVU|^RN7SL", LeftGene, ignore.case = TRUE),
    !grepl("^IG[HKL]|^TRA|^TRB|^TRD|^TRG|^HLA|^RPL|^RPS|^MT-|^RNU|^RNVU|^RN7SL", RightGene, ignore.case = TRUE)
  ) %>%
  dplyr::distinct(sample_name, X.FusionName, LeftBreakpoint, RightBreakpoint, .keep_all = TRUE)

write.table(ft1, "novel_fusions_FT1.tsv", sep = "\t", quote = FALSE, row.names = FALSE, na = "")


# FT2: Additional filtering based on fusion abundance,
# maximum fusion allele ratio, microhomology near breakpoint
# and in-frame protein fusion status

ft2 <- ft1 %>%
  dplyr::filter(
    FFPM >= 0.1, pmax(FAR_left, FAR_right, na.rm = TRUE) >= 0.2,
    num_microh_near_brkpt <= 2 | microh_brkpt_dist > 10000,
    PROT_FUSION_TYPE == "INFRAME"
  ) %>%
  dplyr::distinct(sample_name, X.FusionName, LeftBreakpoint, RightBreakpoint, .keep_all = TRUE)

write.table(ft2, "novel_fusions_FT2.tsv", sep = "\t", quote = FALSE, row.names = FALSE, na = "")


# FT3: Remove fusions already reported in FusionGDB 3.0

known_fuzGDB <- read.table("../fuslist_fusionGDB3.0.txt", header = F, stringsAsFactors = F)

ft3 <- ft2 %>%  
  dplyr::filter(!X.FusionName %in% known_fuzGDB[[1]]) %>% 
  dplyr::distinct(sample_name, X.FusionName, LeftBreakpoint, RightBreakpoint, .keep_all = TRUE)

write.table(ft3, "novel_fusions_FT3.tsv", row.names = F, sep = "\t")


# FT4: Further remove fusions present in the curated FusionGDB 3.0 list

lil_db_llm <- read.table("../fuslist_cur_fuzGDB3.0.txt", header = T, stringsAsFactors = F)

ft4 <- ft3 %>%  
  dplyr::filter(!X.FusionName %in% lil_db_llm[[1]]) %>% 
  dplyr::distinct(sample_name, X.FusionName, LeftBreakpoint, RightBreakpoint, .keep_all = TRUE)

write.table(ft4, "novel_fusions_FT4.tsv", row.names = F, sep = "\t")


# KNOWN FUSIONS: FILTERING
# Extract known fusions from the initial filtered dataset (FT1)

known_fusions <- read.table("../fuslist_cur_fuzGDB3.0.txt", header = F, stringsAsFactors = F)

df_known_fusion <- ft1 %>%  
  dplyr::filter(X.FusionName %in% known_fusions[[1]])

write.table(df_known_fusion, "known_fusion_data.txt", row.names = F, sep = "\t")


# ============================================================
# SAMPLE SUMMARY: COUNT SAMPLES AND FUSIONS AT EACH FILTERING STEP
# cnt = control samples; csf = case samples
# ============================================================

cnt <- read.delim("../control77_combined_FI-star.tsv")
csf_raw <- read.delim("../case801_combined_FI-star.tsv")

csf_nor <- read.delim("./case801_combined_FI-star_normalised_FT0.tsv")
csf_ft1 <- read.delim("./novel_fusions_FT1.tsv")
csf_ft2 <- read.delim("./novel_fusions_FT2.tsv")
csf_ft3 <- read.delim("./novel_fusions_FT3.tsv")
csf_ft4 <- read.delim("./novel_fusions_FT4.tsv")

kfuz <- read.table("known_fusion_data.txt", header = T, stringsAsFactors = F)


# Print Overview of Data
# Summarise the number of unique fusions and samples
# across raw, normalised and filtered datasets

cat("# Fusions:==>", "\n", 
    "Raw Control Sample's Unique Fusion: ", length(unique(cnt$X.FusionName)), "\n",
    "Raw AML Unique Fusion: ", length(unique(csf_raw$X.FusionName)), "\n",
    "Normalised AML Unique Fusion: ", length(unique(csf_nor$X.FusionName)), "\n", 
    "Filter 1 Applied; Fusions: ", length(unique(csf_ft1$X.FusionName)), "\n",
    "Filter 2 Applied; Fusions: ", length(unique(csf_ft2$X.FusionName)), "\n",
    "Filter 3 Applied; Fusions: ", length(unique(csf_ft3$X.FusionName)), "\n",
    "Filter 4 Applied; Fusions: ", length(unique(csf_ft4$X.FusionName)), "\n",
    "Filter Novel Fusios; Fusions: ", length(unique(kfuz$X.FusionName)), "\n",
    "# Samples:==>", "\n", 
    "Raw Control Unique Samples (nControl=77) : ", length(unique(cnt$sample_name)), "\n",
    "Raw AML Unique Samples (nCase=801): ", length(unique(csf_raw$sample_name)), "\n",
    "Normalised AML Unique Samples: ", length(unique(csf_nor$sample_name)), "\n",
    "Filter 1 Applied; Samples: ", length(unique(csf_ft1$sample_name)), "\n",
    "Filter 2 Applied; Samples: ", length(unique(csf_ft2$sample_name)), "\n",
    "Filter 3 Applied; Samples: ", length(unique(csf_ft3$sample_name)), "\n",
    "Filter 4 Applied; Samples: ", length(unique(csf_ft4$sample_name)), "\n",
    "Filter Novel Fusions; Samples: ", length(unique(kfuz$sample_name)), "\n",
    sep = "")


# ============================================================
# KNOWN FUSION SUMMARY
# Count known fusions according to protein fusion type
# (in-frame and frameshift) across fusions and samples
# ============================================================

t1 <- kfuz %>%  
  dplyr::distinct(X.FusionName, LeftGene, RightGene, LeftBreakpoint, RightBreakpoint, 
                  sample_name, .keep_all = T)

# Count known fusions according to protein fusion type
t1 %>% 
  dplyr::count(X.FusionName, PROT_FUSION_TYPE, sort = T) %>% 
  dplyr::count(PROT_FUSION_TYPE)

# Count samples according to protein fusion type
t1 %>% 
  dplyr::count(sample_name, PROT_FUSION_TYPE, sort = T) %>% 
  dplyr::count(PROT_FUSION_TYPE)
```
