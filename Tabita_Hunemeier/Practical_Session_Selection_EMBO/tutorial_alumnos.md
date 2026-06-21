# EMBO Practical Course: Genomic Diversity & Natural Selection Scan (Worksheet)

This worksheet guides you through the analysis of genomic data in humans and canines to identify genomic regions under natural selection. The tutorial is divided into two main parts:

1. **Part 1: Human Genomic Diversity and Natural Selection**: Investigating selection signatures in the candidate gene ***EDAR*** (associated with ectodermal traits in East Asians and Native Americans) using population differentiation ($F_{ST}$, PBS) and haplotype-based metrics (EHH, iHS, XP-EHH).
2. **Part 2: Genomic Selection Scan in Canines**: Identifying the selective sweep at the ***IGF1*** body-size locus by comparing small vs. large dog breeds using PCA, PCAdapt (outlier scan), and haplotype homozygosity methods (XP-nSL and Rsb).

---

# Part 1: Human Genomic Diversity and Natural Selection

## 1. Background and Dataset

### Goal
Our goal is to explore approaches and methods which seek to identify regions of the genome with signatures of natural selection. We will use real genomic data and two classes of tests: one based on population differentiation ($F_{ST}$ / PBS) and another based on extended haplotype homozygosity (EHH / iHS / XP-EHH).

### Dataset
Whole-genome sequencing data from the 1000 Genomes Project Phase III. The full database can be accessed via:
<ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/>

### Data Pre-processing
We will analyze a pre-processed dataset for chromosome 2 corresponding to individuals sampled from the African (AFR: 504 individuals), European (EUR: 503 individuals), and East Asian (EAS: 504 individuals) populations. In this dataset, INDELs, singletons, and SNPs with MAF < 0.05 have been removed. The pairwise $F_{ST}$ was then estimated using `vcftools`.

All data files are located in the `input/` directory:
- `input/Part_1_HumanDiversity/AFR_EAS.weir.fst` (Fst between Africans and East Asians)
- `input/Part_1_HumanDiversity/AFR_EUR.weir.fst` (Fst between Africans and Europeans)
- `input/Part_1_HumanDiversity/EAS_EUR.weir.fst` (Fst between East Asians and Europeans)
- `input/Part_1_HumanDiversity/Chr2_EDAR_LWK_500K.recode.vcf` (Phased African haplotypes around *EDAR*)
- `input/Part_1_HumanDiversity/Chr2_EDAR_CHS_500K.recode.vcf` (Phased East Asian haplotypes around *EDAR*)
- `input/Part_1_HumanDiversity/Chr2_NAM_EAS.weir.fst` (Fst between Native Americans and East Asians in candidate region)
- `input/Part_1_HumanDiversity/Chr2_NAM_EUR.weir.fst` (Fst between Native Americans and Europeans in candidate region)
- `input/Part_1_HumanDiversity/Chr2_EUR_EAS.weir.fst` (Fst between Europeans and East Asians in candidate region)

---

## 2. Genetic Differentiation ($F_{ST}$ and PBS)

### Investigating the Candidate Gene *EDAR*
The human Ectodysplasin A receptor gene, or ***EDAR***, is part of the EDA signaling pathway which specifies prenatally the location, size, and shape of ectodermal appendages (such as hair follicles, teeth, and glands). *EDAR* is a textbook example of positive selection in East Asians. A specific non-synonymous variant, **rs3827760** (chr2:109,513,601 A>G), results in a Val370Ala substitution and is strongly associated with thicker hair shafts and shovel-shaped incisors. Another hypothesis states that *EDAR* acted along with *FADS* and *VDR* in the Beringia Standstill, allowing Native American ancestors to survive in extreme arctic environments.

### Questions for Students

1. **The estimate of $F_{ST}$ by the Weir and Cockerham metric can sometimes generate negative values and "NA". What does that mean? How can this interfere with the results?**
   * *Answer*:
2. **The $F_{ST}$ values observed between pairs of populations for the SNP rs3827760 (position 109,513,601) fall within which distribution quantiles of $F_{ST}$ values for the studied chromosome? Can they be considered outliers?**
   * *Answer*:
3. **From the observed $F_{ST}$ values between population pairs and the significance estimates, what can we say about the rs3827760 SNP differentiation between populations?**
   * *Answer*:
4. **Discuss how these results justify performing another type of analysis based on PBS (Population Branch Statistics).**
   * *Answer*:
5. **What does the PBS analysis reveal? What is the difference between PBS and $F_{ST}$ analysis?**
   * *Answer*:

---

### R Code Exercise: Pairwise $F_{ST}$ Calculation
Write the R code necessary to perform the following:
1. Read the pairwise $F_{ST}$ files from `input/`.
2. Filter duplicate SNP positions and exclude NA values.
3. Align the datasets by overlapping positions.
4. Set negative $F_{ST}$ values to zero.
5. Check Fst values at position `109513601`.
6. Calculate distribution quantiles to determine if rs3827760 is an outlier.
7. Plot pairwise Fst around `109513601` in a 10kb window, highlighting the candidate SNP.

**Write your R code here:**
```R
# ============================================================
# Part 1: Human Genomic Diversity and Natural Selection
# Pairwise FST analysis around EDAR
# Working directory: .../Practical_Session_Selection_EMBO/input
# ============================================================

# Install packages once if needed:
# install.packages(c("dplyr", "readr", "ggplot2", "tidyr"))

library(dplyr)
library(readr)
library(ggplot2)
library(tidyr)

# Candidate EDAR SNP: rs3827760
edar_pos <- 109513601

# ------------------------------------------------------------
# 1. Read the pairwise FST files
# ------------------------------------------------------------

afr_eas <- read_tsv(
  "Part_1_HumanDiversity/AFR_EAS.weir.fst",
  show_col_types = FALSE
)

afr_eur <- read_tsv(
  "Part_1_HumanDiversity/AFR_EUR.weir.fst",
  show_col_types = FALSE
)

eas_eur <- read_tsv(
  "Part_1_HumanDiversity/EAS_EUR.weir.fst",
  show_col_types = FALSE
)

# Check columns. The expected FST column is WEIR_AND_COCKERHAM_FST
colnames(afr_eas)
head(afr_eas)
head(afr_eur)
head(eas_eur)
# ------------------------------------------------------------
# 2. Remove duplicated positions and missing FST values
#    Set negative FST values to zero
# ------------------------------------------------------------

afr_eas_clean <- afr_eas %>%
  distinct(CHROM, POS, .keep_all = TRUE) %>%
  filter(!is.na(WEIR_AND_COCKERHAM_FST)) %>%
  mutate(FST_AFR_EAS = pmax(WEIR_AND_COCKERHAM_FST, 0)) %>%
  select(CHROM, POS, FST_AFR_EAS)

afr_eur_clean <- afr_eur %>%
  distinct(CHROM, POS, .keep_all = TRUE) %>%
  filter(!is.na(WEIR_AND_COCKERHAM_FST)) %>%
  mutate(FST_AFR_EUR = pmax(WEIR_AND_COCKERHAM_FST, 0)) %>%
  select(CHROM, POS, FST_AFR_EUR)

eas_eur_clean <- eas_eur %>%
  distinct(CHROM, POS, .keep_all = TRUE) %>%
  filter(!is.na(WEIR_AND_COCKERHAM_FST)) %>%
  mutate(FST_EAS_EUR = pmax(WEIR_AND_COCKERHAM_FST, 0)) %>%
  select(CHROM, POS, FST_EAS_EUR)

# ------------------------------------------------------------
# 3. Keep only SNP positions present in all three comparisons
# ------------------------------------------------------------

fst_merged <- afr_eas_clean %>%
  inner_join(afr_eur_clean, by = c("CHROM", "POS")) %>%
  inner_join(eas_eur_clean, by = c("CHROM", "POS"))

# Number of shared SNPs
nrow(fst_merged)

# ------------------------------------------------------------
# 4. Extract FST values for EDAR SNP rs3827760
# ------------------------------------------------------------

edar_fst <- fst_merged %>%
  filter(POS == edar_pos)

edar_fst

# If this returns zero rows, check nearby positions:
fst_merged %>%
  filter(POS >= edar_pos - 5, POS <= edar_pos + 5)

# ------------------------------------------------------------
# 5. Calculate FST distribution quantiles
# ------------------------------------------------------------

fst_quantiles <- fst_merged %>%
  summarise(
    AFR_EAS_Q50 = quantile(FST_AFR_EAS, 0.50),
    AFR_EAS_Q90 = quantile(FST_AFR_EAS, 0.90),
    AFR_EAS_Q95 = quantile(FST_AFR_EAS, 0.95),
    AFR_EAS_Q99 = quantile(FST_AFR_EAS, 0.99),

    AFR_EUR_Q50 = quantile(FST_AFR_EUR, 0.50),
    AFR_EUR_Q90 = quantile(FST_AFR_EUR, 0.90),
    AFR_EUR_Q95 = quantile(FST_AFR_EUR, 0.95),
    AFR_EUR_Q99 = quantile(FST_AFR_EUR, 0.99),

    EAS_EUR_Q50 = quantile(FST_EAS_EUR, 0.50),
    EAS_EUR_Q90 = quantile(FST_EAS_EUR, 0.90),
    EAS_EUR_Q95 = quantile(FST_EAS_EUR, 0.95),
    EAS_EUR_Q99 = quantile(FST_EAS_EUR, 0.99)
  )

fst_quantiles

# A clearer quantile table: 50%, 90%, 95%, 99%, 99.9%
quantile_table <- data.frame(
  Quantile = c("50%", "90%", "95%", "99%", "99.9%"),
  AFR_EAS = quantile(
    fst_merged$FST_AFR_EAS,
    probs = c(0.50, 0.90, 0.95, 0.99, 0.999)
  ),
  AFR_EUR = quantile(
    fst_merged$FST_AFR_EUR,
    probs = c(0.50, 0.90, 0.95, 0.99, 0.999)
  ),
  EAS_EUR = quantile(
    fst_merged$FST_EAS_EUR,
    probs = c(0.50, 0.90, 0.95, 0.99, 0.999)
  )
)

quantile_table

# ------------------------------------------------------------
# 6. Calculate percentile rank of EDAR SNP
# ------------------------------------------------------------

edar_percentiles <- data.frame(
  Comparison = c("AFR vs EAS", "AFR vs EUR", "EAS vs EUR"),
  EDAR_FST = c(
    edar_fst$FST_AFR_EAS,
    edar_fst$FST_AFR_EUR,
    edar_fst$FST_EAS_EUR
  ),
  Percentile = c(
    mean(fst_merged$FST_AFR_EAS <= edar_fst$FST_AFR_EAS) * 100,
    mean(fst_merged$FST_AFR_EUR <= edar_fst$FST_AFR_EUR) * 100,
    mean(fst_merged$FST_EAS_EUR <= edar_fst$FST_EAS_EUR) * 100
  )
)

edar_percentiles

# ------------------------------------------------------------
# 7. Plot pairwise FST in a 10 kb window around EDAR
# ------------------------------------------------------------

window_size <- 10000

fst_window <- fst_merged %>%
  filter(
    POS >= edar_pos - window_size,
    POS <= edar_pos + window_size
  )

fst_long <- fst_window %>%
  pivot_longer(
    cols = c(FST_AFR_EAS, FST_AFR_EUR, FST_EAS_EUR),
    names_to = "Comparison",
    values_to = "FST"
  ) %>%
  mutate(
    Comparison = recode(
      Comparison,
      FST_AFR_EAS = "AFR vs EAS",
      FST_AFR_EUR = "AFR vs EUR",
      FST_EAS_EUR = "EAS vs EUR"
    )
  )

ggplot(fst_long, aes(x = POS, y = FST)) +
  geom_point(alpha = 0.7, size = 1.5) +
  geom_vline(
    xintercept = edar_pos,
    linetype = "dashed",
    linewidth = 0.8
  ) +
  facet_wrap(~ Comparison, scales = "free_y") +
  labs(
    title = "Pairwise FST around the EDAR locus",
    subtitle = "10 kb window around rs3827760 (chr2:109,513,601)",
    x = "Genomic position on chromosome 2",
    y = "Weir and Cockerham FST"
  ) +
  theme_classic() +
  theme(
    axis.line.y = element_line(linewidth = 0.8),
    axis.title = element_text(face = "bold")
  )


```

---

### R Code Exercise: Population Branch Statistics (PBS)
Write the R code necessary to:
1. Estimate the Population Branch Statistic for East Asians ($PBS_{EAS}$) using the AFR, EAS, and EUR populations.
2. Convert negative branch lengths to zero.
3. Check the PBS value for the candidate SNP rs3827760.
4. Calculate distribution quantiles to determine if it is an outlier.
5. Plot PBS values around the candidate SNP in a 10kb window.

**Write your R code here:**
```R
# 
# ============================================================
# Population Branch Statistic (PBS) for East Asians
# ============================================================

# PBS transforms FST into branch lengths:
# T = -log(1 - FST)
#
# PBS_EAS = (T_AFR_EAS + T_EAS_EUR - T_AFR_EUR) / 2

# Avoid log(0) if any FST value is exactly 1
# Values extremely close to 1 are capped slightly below 1
fst_merged_pbs <- fst_merged %>%
  mutate(
    FST_AFR_EAS_safe = pmin(FST_AFR_EAS, 0.999999),
    FST_AFR_EUR_safe = pmin(FST_AFR_EUR, 0.999999),
    FST_EAS_EUR_safe = pmin(FST_EAS_EUR, 0.999999),

    T_AFR_EAS = -log(1 - FST_AFR_EAS_safe),
    T_AFR_EUR = -log(1 - FST_AFR_EUR_safe),
    T_EAS_EUR = -log(1 - FST_EAS_EUR_safe),

    PBS_EAS_raw = (T_AFR_EAS + T_EAS_EUR - T_AFR_EUR) / 2,

    # Set negative estimated branch lengths to zero
    PBS_EAS = pmax(PBS_EAS_raw, 0)
  )

# View the new PBS columns
head(fst_merged_pbs)

# ------------------------------------------------------------
# 1. PBS value at the EDAR candidate SNP rs3827760
# ------------------------------------------------------------

edar_pbs <- fst_merged_pbs %>%
  filter(POS == edar_pos) %>%
  select(
    CHROM, POS,
    FST_AFR_EAS, FST_AFR_EUR, FST_EAS_EUR,
    T_AFR_EAS, T_AFR_EUR, T_EAS_EUR,
    PBS_EAS_raw, PBS_EAS
  )

edar_pbs

print(edar_pbs, width = Inf)
# ------------------------------------------------------------
# 2. PBS distribution quantiles
# ------------------------------------------------------------

pbs_quantiles <- quantile(
  fst_merged_pbs$PBS_EAS,
  probs = c(0.50, 0.90, 0.95, 0.99, 0.999),
  na.rm = TRUE
)

pbs_quantiles

# Make the quantiles easier to read
pbs_quantile_table <- data.frame(
  Quantile = c("50%", "90%", "95%", "99%", "99.9%"),
  PBS_EAS = as.numeric(pbs_quantiles)
)

pbs_quantile_table

# ------------------------------------------------------------
# 3. Percentile rank of EDAR PBS value
# ------------------------------------------------------------

edar_pbs_percentile <- mean(
  fst_merged_pbs$PBS_EAS <= edar_pbs$PBS_EAS
) * 100

edar_pbs_percentile

# ------------------------------------------------------------
# 4. Plot PBS_EAS in a 10 kb window around rs3827760
# ------------------------------------------------------------

window_size <- 10000

pbs_window <- fst_merged_pbs %>%
  filter(
    POS >= edar_pos - window_size,
    POS <= edar_pos + window_size
  )

ggplot(pbs_window, aes(x = POS, y = PBS_EAS)) +
  geom_point(alpha = 0.7, size = 1.7) +
  geom_vline(
    xintercept = edar_pos,
    linetype = "dashed",
    linewidth = 0.8
  ) +
  geom_point(
    data = pbs_window %>% filter(POS == edar_pos),
    size = 3
  ) +
  labs(
    title = "PBS for East Asians around the EDAR locus",
    subtitle = "10 kb window around rs3827760 (chr2:109,513,601)",
    x = "Genomic position on chromosome 2",
    y = "PBS_EAS"
  ) +
  theme_classic() +
  theme(
    axis.line.y = element_line(linewidth = 0.8),
    axis.title = element_text(face = "bold")
  )


```

---

## 3. Extended Haplotype Homozygosity (EHH)

### Extended Haplotype Homozygosity (EHH) and Haplotype Sweeps
Different approaches can detect genomic signatures of selection at different timescales. More recent selection signals can be detected from haplotype-based tests. Positive selection causes a rapid rise in the frequency of the selected allele, such that recombination does not have enough time to break down the haplotype on which the mutation arose. This creates a signature of **Extended Haplotype Homozygosity (EHH)** extending over a long physical distance.

### Questions for Students

1. **How is the haplotype profile of genetic variants under recent positive selection?**
   * *Answer*:
2. **What is the profile of ancestral and derived haplotypes of the rs3827760 SNP in AFR and EAS?**
   * *Answer*:
3. **The iHS score observed for the SNP rs3827760 falls within which distribution quantiles of iHS values for the studied chromosome? Can it be considered an outlier? How can we make this analysis more robust?**
   * *Answer*:
4. **What information does the XP-EHH analysis add about natural selection in the candidate SNP?**
   * *Answer*:

---

### R Code Exercise: EHH & Furcation Trees
Write the R code necessary to:
1. Convert the VCF databases to `rehh` format using `data2haplohh()`.
2. Estimate the EHH decay for rs3827760 in both populations.
3. Plot the EHH decay and furcation trees for both AFR and EAS.

**Write your R code here:**
```R
# 
# ============================================================
# Part 1: Human Diversity
# EHH decay and furcation trees around EDAR rs3827760
# ============================================================

# Install once if needed:
# install.packages("rehh")

library(rehh)

# Candidate EDAR SNP
edar_pos <- 109513601

# ------------------------------------------------------------
# 1. Convert the phased VCF files into rehh haplohh objects
# ------------------------------------------------------------

# LWK = African population
haplohh_afr <- data2haplohh(
  hap_file = "Part_1_HumanDiversity/Chr2_EDAR_LWK_500K.recode.vcf",
  haplotype.in.columns = FALSE,
  polarize_vcf = FALSE,
  verbose = TRUE
)

# CHS = East Asian population
haplohh_eas <- data2haplohh(
  hap_file = "Part_1_HumanDiversity/Chr2_EDAR_CHS_500K.recode.vcf",
  haplotype.in.columns = FALSE,
  polarize_vcf = FALSE,
  verbose = TRUE
)


# 2. Find the marker name at the EDAR position
# Marker name
edar_marker <- "rs3827760"

# 3. Calculate EHH
ehh_afr <- calc_ehh(haplohh_afr, mrk = edar_marker)
ehh_eas <- calc_ehh(haplohh_eas, mrk = edar_marker)
# 3. Calculate EHH
ehh_afr <- calc_ehh(haplohh_afr, mrk = edar_marker)
ehh_eas <- calc_ehh(haplohh_eas, mrk = edar_marker)

plot(ehh_afr, main = "AFR")
plot(ehh_eas, main = "EAS")

```

---

### R CodeExercise: iHS & XP-EHH (Window-based)
Write the R code necessary to:
1. Perform a genome-wide scan of homozygosity using `scan_hh()` for AFR and EAS.
2. Calculate iHS scores for both populations using `ihh2ihs()`.
3. Check the iHS score at rs3827760 and generate a single-site iHS plot in EAS.
4. Create a function to estimate the average absolute iHS in sliding windows (50 SNPs/40 step) and plot the results.
5. Estimate cross-population XP-EHH between EAS and AFR using `ies2xpehh()`, calculate window-based averages, and plot them.

**Write your R code here:**
```r
# iHS and XP-EHH analysis

library(rehh)
library(dplyr)
library(ggplot2)

# focal SNP
edar_marker <- "rs3827760"

# ------------------------------------------------------------
# 1. Genome-wide EHH scan
# ------------------------------------------------------------

scan_afr <- scan_hh(haplohh_afr)
scan_eas <- scan_hh(haplohh_eas)

# ------------------------------------------------------------
# 2. Calculate iHS
# ------------------------------------------------------------

ihs_afr <- ihh2ihs(scan_afr)
ihs_eas <- ihh2ihs(scan_eas)

# ------------------------------------------------------------
# 3. iHS at rs3827760 and single-site plot in EAS
# ------------------------------------------------------------

ihs_eas[ihs_eas$MARKER == edar_marker, ]

plot(
  ihs_eas$POSITION,
  ihs_eas$IHS,
  pch = 16,
  cex = 0.5,
  xlab = "Position on chromosome 2",
  ylab = "iHS",
  main = "iHS scan in EAS"
)

abline(v = 109513601, lty = 2)

####################################
# ------------------------------------------------------------
# iHS and XP-EHH analysis around EDAR
# ------------------------------------------------------------

library(rehh)
library(dplyr)
library(ggplot2)

edar_marker <- "rs3827760"
edar_pos <- 109513601

# ------------------------------------------------------------
# 1. Genome-wide EHH scan
# ------------------------------------------------------------

scan_afr <- scan_hh(haplohh_afr)
scan_eas <- scan_hh(haplohh_eas)

# ------------------------------------------------------------
# 2. Calculate iHS
# ------------------------------------------------------------

ihs_afr <- ihh2ihs(scan_afr)
ihs_eas <- ihh2ihs(scan_eas)

# Extract the actual iHS tables
ihs_afr_table <- ihs_afr$ihs
ihs_eas_table <- ihs_eas$ihs

# Check column names
colnames(ihs_afr_table)
colnames(ihs_eas_table)

# ------------------------------------------------------------
# 3. iHS at rs3827760 and iHS plot in EAS
# ------------------------------------------------------------

# iHS value at EDAR
ihs_eas_table[ihs_eas_table$MARKER == edar_marker, ]

# Plot iHS across the EAS region
plot(
  ihs_eas_table$POSITION,
  ihs_eas_table$IHS,
  pch = 16,
  cex = 0.5,
  xlab = "Position on chromosome 2",
  ylab = "iHS",
  main = "iHS scan in EAS"
)

abline(v = edar_pos, lty = 2)

# ------------------------------------------------------------
# 4. Average absolute iHS in sliding windows
#    50 SNPs per window, step = 40 SNPs
# ------------------------------------------------------------

window_ihs <- function(ihs_df, window_size = 50, step = 40) {
  
  ihs_df <- ihs_df %>%
    filter(!is.na(IHS)) %>%
    arrange(POSITION)
  
  starts <- seq(
    1,
    nrow(ihs_df) - window_size + 1,
    by = step
  )
  
  bind_rows(lapply(starts, function(i) {
    
    x <- ihs_df[i:(i + window_size - 1), ]
    
    data.frame(
      MID_POS = mean(x$POSITION),
      MEAN_ABS_IHS = mean(abs(x$IHS))
    )
  }))
}

window_ihs_afr <- window_ihs(ihs_afr_table)
window_ihs_eas <- window_ihs(ihs_eas_table)

# Plot windowed mean absolute iHS in EAS
ggplot(window_ihs_eas, aes(x = MID_POS / 1e6, y = MEAN_ABS_IHS)) +
  geom_line() +
  geom_vline(xintercept = edar_pos / 1e6, linetype = "dashed") +
  labs(
    title = "Windowed mean absolute iHS: EAS",
    x = "Position on chromosome 2 (Mb)",
    y = "Mean |iHS|"
  ) +
  theme_classic()

# ------------------------------------------------------------
# 5. XP-EHH: EAS compared with AFR
# ------------------------------------------------------------

xpehh_eas_afr <- ies2xpehh(
  scan_eas,
  scan_afr
)


# ------------------------------------------------------------
# 5. XP-EHH: EAS compared with AFR
# ------------------------------------------------------------

xpehh_eas_afr <- ies2xpehh(scan_eas, scan_afr)

# In this rehh version, the result is already a data frame
xpehh_table <- xpehh_eas_afr

# XP-EHH at EDAR
xpehh_table[xpehh_table$POSITION == edar_pos, ]

# ------------------------------------------------------------
# 6. Sliding-window average XP-EHH
# ------------------------------------------------------------

window_xpehh <- function(xpehh_df, window_size = 50, step = 40) {
  
  xpehh_df <- xpehh_df %>%
    filter(!is.na(XPEHH)) %>%
    arrange(POSITION)
  
  starts <- seq(
    1,
    nrow(xpehh_df) - window_size + 1,
    by = step
  )
  
  bind_rows(lapply(starts, function(i) {
    
    x <- xpehh_df[i:(i + window_size - 1), ]
    
    data.frame(
      MID_POS = mean(x$POSITION),
      MEAN_XPEHH = mean(x$XPEHH)
    )
  }))
}

window_xpehh_eas_afr <- window_xpehh(xpehh_table)

ggplot(window_xpehh_eas_afr, aes(x = MID_POS / 1e6, y = MEAN_XPEHH)) +
  geom_line() +
  geom_vline(xintercept = edar_pos / 1e6, linetype = "dashed") +
  labs(
    title = "Windowed XP-EHH: EAS vs AFR",
    x = "Position on chromosome 2 (Mb)",
    y = "Mean XP-EHH"
  ) +
  theme_classic()
####################################
# ------------------------------------------------------------
# 4. Average absolute iHS in sliding windows
#    50 SNPs per window, step = 40 SNPs
# ------------------------------------------------------------

window_ihs <- function(ihs_df, window_size = 50, step = 40) {
  
  ihs_df <- ihs_df %>%
    filter(!is.na(IHS)) %>%
    arrange(POSITION)
  
  starts <- seq(1, nrow(ihs_df) - window_size + 1, by = step)
  
  bind_rows(lapply(starts, function(i) {
    
    x <- ihs_df[i:(i + window_size - 1), ]
    
    data.frame(
      MID_POS = mean(x$POSITION),
      MEAN_ABS_IHS = mean(abs(x$IHS))
    )
  }))
}

window_ihs_afr <- window_ihs(ihs_afr)
window_ihs_eas <- window_ihs(ihs_eas)

# Plot windowed iHS for EAS
ggplot(window_ihs_eas, aes(x = MID_POS / 1e6, y = MEAN_ABS_IHS)) +
  geom_line() +
  geom_vline(xintercept = 109513601 / 1e6, linetype = "dashed") +
  labs(
    title = "Windowed mean absolute iHS: EAS",
    x = "Position on chromosome 2 (Mb)",
    y = "Mean |iHS|"
  ) +
  theme_classic()

# ------------------------------------------------------------
# 5. XP-EHH: EAS compared with AFR
# ------------------------------------------------------------

xpehh_eas_afr <- ies2xpehh(
  scan_eas,
  scan_afr
)

# Check XP-EHH at EDAR
xpehh_eas_afr[xpehh_eas_afr$MARKER == edar_marker, ]

# Function for sliding-window XP-EHH
window_xpehh <- function(xpehh_df, window_size = 50, step = 40) {
  
  xpehh_df <- xpehh_df %>%
    filter(!is.na(XPEHH)) %>%
    arrange(POSITION)
  
  starts <- seq(1, nrow(xpehh_df) - window_size + 1, by = step)
  
  bind_rows(lapply(starts, function(i) {
    
    x <- xpehh_df[i:(i + window_size - 1), ]
    
    data.frame(
      MID_POS = mean(x$POSITION),
      MEAN_XPEHH = mean(x$XPEHH)
    )
  }))
}

window_xpehh_eas_afr <- window_xpehh(xpehh_eas_afr)

# Plot windowed XP-EHH
ggplot(window_xpehh_eas_afr, aes(x = MID_POS / 1e6, y = MEAN_XPEHH)) +
  geom_line() +
  geom_vline(xintercept = 109513601 / 1e6, linetype = "dashed") +
  labs(
    title = "Windowed XP-EHH: EAS vs AFR",
    x = "Position on chromosome 2 (Mb)",
    y = "Mean XP-EHH"
  ) +
  theme_classic()


```

---

## 4. Native American Selection Analysis

### Background
Hlusko et al. (2018), using morphological data, found a strong selection signal in the *EDAR* gene in Native Americans. Using the additional database from the 1000 Genomes Project (Peruvian samples with over 95% Native American Ancestry, represented as **NAM**), we evaluate genomic signatures of selection at the functional variant rs3827760.

### Questions for Students

1. **Is the functional allele in East Asian at high frequency in other human populations (e.g. Native Americans)?**
   * *Answer*:
2. **Can we identify signatures of natural selection on EDAR in Native Americans using PBS?**
   * *Answer*:
3. **Is selection targeting the same functional variant?**
   * *Answer*:
4. **What is your conclusion based on the data generated?**
   * *Answer*:

---

### R Code Exercise: PBS in Native Americans (NAM)
Write the R code necessary to:
1. Read the pairwise $F_{ST}$ files involving Native Americans (`input/Part_1_HumanDiversity/Chr2_NAM_EAS.weir.fst` and `input/Part_1_HumanDiversity/Chr2_NAM_EUR.weir.fst`) and Europeans-East Asians (`input/Part_1_HumanDiversity/Chr2_EUR_EAS.weir.fst`).
2. Filter duplicates, exclude NA values, and align positions.
3. Convert negative $F_{ST}$ values to zero.
4. Estimate $PBS_{NAM}$ using NAM, EAS, and EUR.
5. Check PBS value at rs3827760, check quantiles, and plot the PBS scan.

**Write your R code here:**
```r

library(dplyr)
library(readr)
library(ggplot2)

edar_pos <- 109513601

# Read FST files
nam_eas <- read_tsv(
  "Part_1_HumanDiversity/Chr2_NAM_EAS.weir.fst",
  show_col_types = FALSE
)

nam_eur <- read_tsv(
  "Part_1_HumanDiversity/Chr2_NAM_EUR.weir.fst",
  show_col_types = FALSE
)

eur_eas <- read_tsv(
  "Part_1_HumanDiversity/Chr2_EUR_EAS.weir.fst",
  show_col_types = FALSE
)

# Keep one row per position, remove NA, rename FST column
clean_fst <- function(df, new_name) {
  df %>%
    distinct(POS, .keep_all = TRUE) %>%
    filter(!is.na(WEIR_AND_COCKERHAM_FST)) %>%
    select(CHROM, POS, FST = WEIR_AND_COCKERHAM_FST) %>%
    rename(!!new_name := FST)
}

nam_eas_clean <- clean_fst(nam_eas, "FST_NAM_EAS")
nam_eur_clean <- clean_fst(nam_eur, "FST_NAM_EUR")
eur_eas_clean <- clean_fst(eur_eas, "FST_EUR_EAS")

# Align SNPs present in all three comparisons
pbs_nam <- nam_eas_clean %>%
  inner_join(nam_eur_clean, by = c("CHROM", "POS")) %>%
  inner_join(eur_eas_clean, by = c("CHROM", "POS")) %>%
  mutate(
    # Negative FST values are set to zero
    FST_NAM_EAS = pmax(FST_NAM_EAS, 0),
    FST_NAM_EUR = pmax(FST_NAM_EUR, 0),
    FST_EUR_EAS = pmax(FST_EUR_EAS, 0),

    # Prevent log(0) if an FST value equals 1
    FST_NAM_EAS_safe = pmin(FST_NAM_EAS, 0.999999),
    FST_NAM_EUR_safe = pmin(FST_NAM_EUR, 0.999999),
    FST_EUR_EAS_safe = pmin(FST_EUR_EAS, 0.999999),

    # Transform FST values into distances
    T_NAM_EAS = -log(1 - FST_NAM_EAS_safe),
    T_NAM_EUR = -log(1 - FST_NAM_EUR_safe),
    T_EUR_EAS = -log(1 - FST_EUR_EAS_safe),

    # PBS for Native Americans
    PBS_NAM_raw = (T_NAM_EAS + T_NAM_EUR - T_EUR_EAS) / 2,

    # Negative branch lengths are not biologically meaningful
    PBS_NAM = pmax(PBS_NAM_raw, 0)
  )

# PBS value at EDAR rs3827760
edar_pbs_nam <- pbs_nam %>%
  filter(POS == edar_pos)

edar_pbs_nam

# Distribution quantiles
pbs_nam_quantiles <- quantile(
  pbs_nam$PBS_NAM,
  probs = c(0.50, 0.90, 0.95, 0.99, 0.999),
  na.rm = TRUE
)

pbs_nam_quantiles

# Percentile of EDAR PBS among all SNPs
edar_pbs_percentile <- mean(
  pbs_nam$PBS_NAM <= edar_pbs_nam$PBS_NAM
) * 100

edar_pbs_percentile

# Plot PBS around EDAR: 10 kb window
pbs_nam_window <- pbs_nam %>%
  filter(
    POS >= edar_pos - 10000,
    POS <= edar_pos + 10000
  )

ggplot(pbs_nam_window, aes(x = POS / 1e6, y = PBS_NAM)) +
  geom_point(size = 1.5, alpha = 0.7) +
  geom_vline(
    xintercept = edar_pos / 1e6,
    linetype = "dashed"
  ) +
  labs(
    title = "PBS scan around EDAR in Native Americans",
    x = "Position on chromosome 2 (Mb)",
    y = "PBS_NAM"
  ) +
  theme_classic()

```

---
---

# Part 2: Genomic Selection Sweep Scan in Canines

## 1. Background and Dataset

The dataset is sourced from the **Dog10K** consortium ([Download Link](https://dog10k.kiz.ac.cn/Home/Download)). The original genomic dataset is a high-coverage phased BCF file containing 1,929 individuals and over 29 million SNPs:
- Original file: `AutoAndXPAR.Dog10K.phased.bcf`
- Metadata table: `dog10K-alignment-sample-table.2022-02-23-v7.txt`

### Sample Selection
For this analysis, we will use a subset of **130 individuals** representing body size extremes:

| Group | Dog Breed (Breed.Type) | Number of Samples |
| :--- | :--- | :---: |
| **Small** (61) | Dachshund | 17 |
| | Toy Fox Terrier | 10 |
| | Pomeranian | 8 |
| | Brussels Griffon | 7 |
| | Yorkshire Terrier | 5 |
| | Shih Tzu | 5 |
| | Maltese | 4 |
| | Pekingese | 3 |
| | Chihuahua | 2 |
| **Large** (69) | Saint Bernard | 13 |
| | Leonberger | 11 |
| | Bernese Mountain Dog | 10 |
| | Greater Swiss Mountain Dog | 10 |
| | Great Pyrenees | 7 |
| | Bullmastiff | 6 |
| | Mastiff | 6 |
| | Newfoundland | 6 |

---

## 2. Preprocessing and Filtering

We filter the massive BCF file to include only our 130 samples and chromosome 15, while removing low-frequency SNPs (MAF < 0.05) that are not informative for breed/size differentiation. The coordinates correspond to the CanFam3 reference genome, where the ***IGF1*** gene is located approximately in the **41.2 Mb to 44.5 Mb** region.

### Bash Code Exercise 1: Extraction & Format Conversion
```bash
# 1. Extract chr15 and filter for samples and MAF >= 0.05
bcftools view \
  -S input/Part_2_CanidDiversity/subset_dogs.txt \
  -r chr15 \
  -q 0.05:minor \
  -O b \
  -o input/Part_2_CanidDiversity/subset_chr15.bcf \
  input/Part_2_CanidDiversity/AutoAndXPAR.Dog10K.phased.bcf

# 2. Index the subset BCF
bcftools index input/Part_2_CanidDiversity/subset_chr15.bcf

# 3. Convert BCF to PLINK binary format (.bed/.bim/.fam)
plink1.9 \
  --bcf input/Part_2_CanidDiversity/subset_chr15.bcf \
  --dog \
  --keep-allele-order \
  --make-bed \
  --out output/subset_chr15
```

> **Premise**: The filtered BCF file contains **177,953 SNPs** on chromosome 15 across the 130 samples. The coordinates correspond to the CanFam3 reference genome, where the ***IGF1*** gene is located approximately in the **41.2 Mb to 44.5 Mb** region.

---

## 3. Population Structure Analysis (PCA)

Before checking for selection outliers, we must examine the genetic structure of our subset. We will run a PCA using **PLINK 1.9** and visualize it in R.

### R Code Exercise 2: PCA Visualization in R
Write the R code necessary to:
1. Load the eigenvectors and eigenvalues generated by PLINK.
2. Merge them with the sample metadata (`input/Part_2_CanidDiversity/sample_info.txt`).
3. Calculate the percentage of variance explained by PC1 and PC2.
4. Generate a scatter plot of PC1 vs. PC2, coloring the points by breed and shaping them by size group (small vs. large).

**Write your R code here:**
```R
# 
# ============================================================
# Part 2: Canid Diversity
# PCA visualization: small vs large dog breeds
# Working directory:
# .../Practical_Session_Selection_EMBO/input
# ============================================================

library(dplyr)
library(readr)
library(ggplot2)

# ------------------------------------------------------------
# 1. Load PLINK PCA eigenvectors and eigenvalues
# ------------------------------------------------------------

# PLINK .eigenvec format:
# Column 1 = FID
# Column 2 = IID
# Remaining columns = PCs

eigenvec <- read.table(
  "Part_2_CanidDiversity/plink_pca.eigenvec",
  header = FALSE,
  stringsAsFactors = FALSE
)

# Give columns meaningful names
colnames(eigenvec) <- c(
  "FID",
  "IID",
  paste0("PC", seq_len(ncol(eigenvec) - 2))
)

# PLINK .eigenval format:
# one eigenvalue per row
eigenval <- scan(
  "Part_2_CanidDiversity/plink_pca.eigenval",
  what = numeric()
)

# Check files
head(eigenvec)
head(eigenval)

# ------------------------------------------------------------
# 2. Calculate variance explained by PC1 and PC2
# ------------------------------------------------------------

variance_explained <- eigenval / sum(eigenval) * 100

pc_variance <- data.frame(
  PC = paste0("PC", seq_along(variance_explained)),
  Variance_explained_percent = variance_explained
)

pc_variance[1:2, ]

pc1_percent <- round(variance_explained[1], 2)
pc2_percent <- round(variance_explained[2], 2)

# ------------------------------------------------------------
# 3–5. Load metadata and merge with PCA coordinates
# ------------------------------------------------------------

sample_info <- read.table(
  "Part_2_CanidDiversity/sample_info.txt",
  header = TRUE,
  sep = "\t",
  stringsAsFactors = FALSE,
  check.names = FALSE
)

# Rename using the actual metadata columns
sample_info_clean <- sample_info %>%
  transmute(
    IID = sampleName,
    Size_group = group,
    Breed = breed
  )

# Check IDs in the PCA file before merging
head(eigenvec[, c("FID", "IID")])
head(sample_info_clean)

# Merge using IID because sample_info has no FID column
pca_data <- eigenvec %>%
  left_join(sample_info_clean, by = "IID")

# Check whether metadata matched
table(is.na(pca_data$Breed))
table(is.na(pca_data$Size_group))

# Keep small and large dogs only for this PCA plot
pca_plot_data <- pca_data %>%
  filter(Size_group %in% c("small", "large")) %>%
  mutate(
    Size_group = factor(
      Size_group,
      levels = c("small", "large")
    )
  )

# ------------------------------------------------------------
# 6. Plot PC1 versus PC2
# ------------------------------------------------------------

ggplot(
  pca_plot_data,
  aes(
    x = PC1,
    y = PC2,
    color = Breed,
    shape = Size_group
  )
) +
  geom_point(size = 3, alpha = 0.85) +
  labs(
    title = "PCA of Canid Genomic Variation",
    subtitle = "Small and large dog breeds",
    x = paste0("PC1 (", pc1_percent, "% variance explained)"),
    y = paste0("PC2 (", pc2_percent, "% variance explained)"),
    color = "Breed",
    shape = "Size group"
  ) +
  theme_classic() +
  theme(
    axis.line = element_line(linewidth = 0.8),
    axis.title = element_text(face = "bold"),
    legend.title = element_text(face = "bold")
  )

```

### Questions for Students
1. **What pattern do you observe along the first principal component (PC1)?**
   * *Answer*:
2. **Why does PC1 capture body size differences in this particular dataset?**
   * *Answer*:

---

## 4. Genomic Outlier Detection using PCAdapt

**PCAdapt** is a method designed to find SNPs that are exceptionally related to population structure (PCs) rather than neutral drift.

### The Role of Linkage Disequilibrium (LD) and Clumping
> [!IMPORTANT]
> **Key Concept**: PCA is highly sensitive to Linkage Disequilibrium (LD). If a region contains many highly correlated markers (due to a selective sweep or low recombination), that single region will dominate the principal components, biasing the PCA and masking other genomic signals (such as the *IGF1* gene sweep).
> 
> To resolve this, we must enable **LD Clumping** in PCAdapt. Thinning out redundant SNPs in strong LD allows the global genomic structure to be correctly computed and helps locate narrow selection sweeps.

### R Code Exercise 3: PCAdapt with LD Clumping
Write the R script necessary to:
1. Load the genotypes into PCAdapt.
2. Execute `pcadapt()` using $K = 2$ components and enable **LD clumping** (for example, with a window size of 500 SNPs and an $r^2$ threshold of 0.1).
3. Merge the resulting p-values with the physical genomic positions from the `.bim` file.
4. Generate a Manhattan Plot of the results, plotting physical position (in Mb) on the X-axis and $-\log_{10}(\text{p-value})$ on the Y-axis.

**Write your R code here:**
```R
# 
# ============================================================
# Part 2: PCAdapt outlier scan with LD clumping
# Working directory:
# .../Practical_Session_Selection_EMBO/input
# ============================================================

# Install once if needed:
#install.packages("pcadapt")

library(pcadapt)
library(dplyr)
library(ggplot2)

# ------------------------------------------------------------
# 1. Load the PLINK genotype data
# ------------------------------------------------------------

plink_prefix <- "Part_2_CanidDiversity/subset_chr15.bed"

# Read the PLINK .bed/.bim/.fam dataset
genotypes <- read.pcadapt(
  input = plink_prefix,
  type = "bed"
)

# ------------------------------------------------------------
# 2. Run PCAdapt using K = 2 and LD clumping
# ------------------------------------------------------------

# LD clumping:
# size = 500 SNPs per window
# threshold = r^2 = 0.1
#
# LD clumping is performed before PCA so highly correlated SNPs
# do not dominate the PCs.

pcadapt_result <- pcadapt(
  input = genotypes,
  K = 2,
  method = "componentwise",
  min.maf = 0.05,
  LD.clumping = list(
    size = 500,
    threshold = 0.1
  )
)

# Inspect result
pcadapt_result

# ------------------------------------------------------------
# 3. Extract p-values and merge with physical positions
# ------------------------------------------------------------

# Read the PLINK BIM file
# Columns: chromosome, SNP ID, genetic distance, physical position, allele1, allele2
bim <- read.table(
  "Part_2_CanidDiversity/subset_chr15.bim",
  header = FALSE,
  stringsAsFactors = FALSE
)

colnames(bim) <- c(
  "CHROM",
  "SNP",
  "GENETIC_DISTANCE",
  "POS",
  "ALLELE1",
  "ALLELE2"
)

# Confirm that number of SNPs matches number of p-values
nrow(bim)
length(pcadapt_result$pvalues)

# Create a data frame for plotting
# Check the dimensions first
dim(pcadapt_result$pvalues)

# Combine the PC1 and PC2 p-values:
# take the smallest p-value for each SNP
pvalue_combined <- apply(
  pcadapt_result$pvalues,
  1,
  min,
  na.rm = TRUE
)

# Create plotting data frame
pcadapt_df <- bim %>%
  mutate(
    pvalue = pvalue_combined,
    POS_MB = POS / 1e6,
    neg_log10_p = -log10(pvalue)
  ) %>%
  filter(!is.na(pvalue), pvalue > 0)

# ------------------------------------------------------------
# 4. Manhattan plot
# ------------------------------------------------------------

ggplot(pcadapt_df, aes(x = POS_MB, y = neg_log10_p)) +
  geom_point(alpha = 0.65, size = 1.2) +
  labs(
    title = "PCAdapt genomic outlier scan: chromosome 15",
    subtitle = "K = 2 PCs with LD clumping (500 SNP window, r² threshold = 0.1)",
    x = "Physical position on chromosome 15 (Mb)",
    y = expression(-log[10](p-value))
  ) +
  theme_classic() +
  theme(
    axis.line = element_line(linewidth = 0.8),
    axis.title = element_text(face = "bold")
  )
  
  
#make the graph pretty

# ------------------------------------------------------------
# Highlight the IGF1 peak around 43.4 Mb
# ------------------------------------------------------------

# Define the IGF1 region to highlight
igf1_center <- 43.4
igf1_window <- 0.5   # highlights SNPs from 42.9 to 43.9 Mb

pcadapt_df <- pcadapt_df %>%
  mutate(
    Region = ifelse(
      POS_MB >= (igf1_center - igf1_window) &
      POS_MB <= (igf1_center + igf1_window),
      "IGF1 peak (43.4 Mb)",
      "Other SNPs"
    ),
    Region = factor(
      Region,
      levels = c("Other SNPs", "IGF1 peak (43.4 Mb)")
    )
  )

# ------------------------------------------------------------
# Manhattan plot with IGF1 peak highlighted
# ------------------------------------------------------------

ggplot(pcadapt_df, aes(x = POS_MB, y = neg_log10_p, color = Region)) +
  
  geom_point(size = 1.2, alpha = 0.75) +
  
  geom_vline(
    xintercept = igf1_center,
    linetype = "dashed",
    linewidth = 0.7
  ) +
  
  scale_color_manual(
    values = c(
      "Other SNPs" = "grey40",
      "IGF1 peak (43.4 Mb)" = "turquoise4"
    )
  ) +
  
  annotate(
    "text",
    x = igf1_center,
    y = max(pcadapt_df$neg_log10_p, na.rm = TRUE) * 0.95,
    label = "IGF1\n43.4 Mb",
    fontface = "bold",
    size = 4,
    vjust = 1
  ) +
  
  labs(
    title = "PCAdapt Genomic Outlier Scan: Chromosome 15",
    subtitle = "K = 2 PCs with LD clumping (500 SNP window, r² threshold = 0.1)",
    x = "Physical position on chromosome 15 (Mb)",
    y = expression(-log[10](p-value)),
    color = NULL
  ) +
  
  theme_classic() +
  theme(
    axis.line = element_line(linewidth = 0.8),
    axis.title = element_text(face = "bold"),
    legend.position = "top",
    legend.direction = "horizontal",
    legend.text = element_text(size = 11),
    legend.background = element_rect(
      fill = "white",
      color = "black",
      linewidth = 0.4
    )
  )
```

### Questions for Students
1. **If you ran PCAdapt without enabling LD clumping, a single massive peak at ~61 Mb would dominate the plot, hiding other regions. What is the effect of LD clumping on outlier detection and why is it necessary here?**
   * *Answer*:
2. **Did you detect the outlier peak at the *IGF1* locus? What is the approximate coordinate of the peak and what is its biological significance?**
   * *Answer*:

---

## 5. Cross-Population Selection Scan using XP-nSL

To confirm that the outlier peak on chromosome 15 is indeed driven by a selective sweep, we will perform a haplotype-based selection scan. Specifically, we will run **XP-nSL (Cross-Population Number of Segregating Sites by Length)** to compare the haplotype homozygosity decay between small dogs and large dogs.

### Key Concepts: nSL and XP-nSL
- **nSL**: A within-population selection scan metric similar to iHS. However, instead of measuring haplotype decay in terms of genetic distance (which requires a genetic map), nSL measures distance by counting the number of segregating sites (segregating site count by length). This makes it highly robust to recombination rate variation and suitable for genomes without well-defined genetic maps.
- **XP-nSL**: A cross-population statistic that compares nSL profiles between a target population and a reference population. A high positive score indicates a selective sweep specific to the target population (longer haplotypes around the derived allele).
- **Phased Mode**: Since our input Dog10K BCF file is already phased (containing haplotype data formatted as `0|0`, `1|0`, etc.), we will perform a phased XP-nSL scan. This utilizes the precise haplotype sequences, which provides a significantly stronger selection signal compared to unphased analyses.

### Outgroup Allele Polarization
Haplotype selection scans require knowing which allele is **ancestral** (original) and which is **derived** (new mutant). `selscan` expects a VCF file where `0` is the ancestral allele and `1` is the derived allele.
To polarize our dataset, we use the **gray wolves** in the Dog10K metadata as an outgroup:
- Gray wolves are the evolutionary ancestor of domestic dogs.
- For each SNP, the most common (major) allele in the gray wolf population is designated as the ancestral allele.
- If the ALT allele in the original VCF is the major allele in wolves, we must physically swap the REF/ALT alleles and swap genotypes (`0` becomes `1`, and `1` becomes `0`) for all individuals.

---

### Bash Code Exercise 4: Extracting and Polarizing Alleles
```bash
# 1. Extract wolf samples and combine with dog samples
bcftools query -l input/Part_2_CanidDiversity/AutoAndXPAR.Dog10K.phased.bcf | grep -E '^CLUP' > input/Part_2_CanidDiversity/wolves.txt
cat input/Part_2_CanidDiversity/subset_dogs.txt input/Part_2_CanidDiversity/wolves.txt > input/Part_2_CanidDiversity/dogs_and_wolves.txt

# 2. Extract polymorphic sites for both dogs and wolves
bcftools query -f '%CHROM\t%POS\n' input/Part_2_CanidDiversity/subset_chr15.bcf > input/Part_2_CanidDiversity/subset_chr15_positions.txt
bcftools view \
  -S input/Part_2_CanidDiversity/dogs_and_wolves.txt \
  -T input/Part_2_CanidDiversity/subset_chr15_positions.txt \
  -O z \
  -o input/Part_2_CanidDiversity/subset_chr15_with_wolves.vcf.gz \
  input/Part_2_CanidDiversity/AutoAndXPAR.Dog10K.phased.bcf

# 3. Polarize alleles using Python (major allele in wolves = 0)
python3 scripts/polarize_by_wolves.py

# 4. Re-compress to block gzip format and index polarized VCF
bcftools view input/Part_2_CanidDiversity/subset_chr15_polarized.vcf.gz -O z -o output/subset_chr15_polarized_bgzf.vcf.gz
mv output/subset_chr15_polarized_bgzf.vcf.gz input/Part_2_CanidDiversity/subset_chr15_polarized.vcf.gz
bcftools index input/Part_2_CanidDiversity/subset_chr15_polarized.vcf.gz

# 5. Extract polarized VCFs for small and large dogs separately
bcftools view -S input/Part_2_CanidDiversity/small_dogs.txt -O z -o input/Part_2_CanidDiversity/small_dogs_polarized.vcf.gz input/Part_2_CanidDiversity/subset_chr15_polarized.vcf.gz
bcftools index input/Part_2_CanidDiversity/small_dogs_polarized.vcf.gz

bcftools view -S input/Part_2_CanidDiversity/large_dogs.txt -O z -o input/Part_2_CanidDiversity/large_dogs_polarized.vcf.gz input/Part_2_CanidDiversity/subset_chr15_polarized.vcf.gz
bcftools index input/Part_2_CanidDiversity/large_dogs_polarized.vcf.gz
```


library(readr)
library(dplyr)
library(ggplot2)

# Read SNP-level normalized XP-nSL results
xpnsl <- read_tsv(
  "Part_2_CanidDiversity/xpnsl_phased.xpnsl.out.norm",
  show_col_types = FALSE
)

# Prepare data
xpnsl <- xpnsl %>%
  mutate(
    POS_MB = pos / 1e6,
    Region = ifelse(
      POS_MB >= 43.3 & POS_MB <= 43.5,
      "IGF1 region",
      "Other regions"
    )
  )

# SNP-level XP-nSL plot
ggplot(xpnsl, aes(x = POS_MB, y = norm_xpnsl, color = Region)) +
  geom_point(size = 1, alpha = 0.7) +
  geom_vline(xintercept = 43.4, linetype = "dashed") +
  labs(
    title = "XP-nSL scan: Small dogs vs Large dogs",
    x = "Position on chromosome 15 (Mb)",
    y = "Normalized XP-nSL",
    color = NULL
  ) +
  theme_classic() +
  theme(
    legend.position = "top",
    axis.title = element_text(face = "bold")
  )
---

### Selection Scan Execution and Normalization
We will execute the selection scan using the compiled `selscan` binary, comparing small dogs (target) against large dogs (reference) in phased mode, and normalize the raw scores.

### Bash Code Exercise 5: Running XP-nSL and Normalizing
Write the bash commands necessary to:
1. Execute the `selscan` program in phased XP-nSL mode using the polarized target (small dogs) and reference (large dogs) VCFs. Use 4 threads and output results to `output/xpnsl_phased`.
2. Normalize the raw XP-nSL scores using the `selscan norm` utility.
3. Calculate window-based statistics (fraction of outliers) in 100 Kb non-overlapping windows.

**Write your Bash code here:**
```bash
# 


```

---

### R Code Exercise 6: XP-nSL Haplotype Manhattan Plots
Write the R code necessary to:
1. Load the normalized SNP-level XP-nSL scores from `input/Part_2_CanidDiversity/xpnsl_phased.xpnsl.out.norm`.
2. Load the 100 Kb window-level data from `input/Part_2_CanidDiversity/xpnsl_phased.xpnsl.out.norm.100kb.windows`.
3. Highlight the *IGF1* region (between 41 Mb and 45.5 Mb).
4. Generate two Manhattan plots: one for the raw SNP-level scores (calculating $-\log_{10}(p\text{-value})$ from the normalized Z-scores), and another for the window-based fraction of extreme positive SNPs (`frac_top`).

**Write your R code here:**
```R
# 


```

### Questions for Students
1. **Why do we need to polarize alleles using an outgroup like the gray wolf? What does the `0` vs `1` coding represent in `selscan`?**
   * *Answer*:
2. **Why does the raw SNP-level XP-nSL scan look like a noisy cloud at individual sites? What is the effect of calculating window-based scores (e.g. 100 Kb)?**
   * *Answer*:

---

## 6. Alternative Haplotype Selection Scan: Rsb using `rehh`

As a complementary approach to XP-nSL, we will run **Rsb**, another widely used cross-population EHH-based statistic.

### Key Concepts & Comparison
- **Rsb**: Compares the integrated Extended Haplotype Homozygosity (iHH) between two populations. It is calculated as $\ln(iES_{pop1} / iES_{pop2})$, where $iES$ is the integrated EHH over physical distance (bp). A high positive score indicates selection in population 1 (small dogs), while a negative score indicates selection in population 2 (large dogs).
- **Difference from XP-nSL**: 
  - **XP-nSL** integrates the nSL metric over the number of segregating sites (SNP count). This makes it highly robust to local recombination rate variation.
  - **Rsb** integrates EHH over physical distance (bp). In species with very strong selective sweeps and long-range Linkage Disequilibrium (like domestic dogs), Rsb can produce exceptionally high, clear peaks at sweep loci like *IGF1*.

### Phase Preservation during Outgroup Polarization
A common question in haplotype selection scans is: **Does swapping the REF and ALT alleles (and flipping 0 to 1 and 1 to 0) to align with the outgroup corrupt or destroy the phasing information?** The answer is **No**. In a phased VCF, genotypes are represented as `0|1` or `1|0` to denote which allele lies on which homologous chromosome (haplotype). Swapping REF and ALT alleles and swapping `0` and `1` (converting `0|1` to `1|0` and vice-versa) is a mathematically symmetric operation. It maintains the exact same physical haplotype alignment across all sites on each chromosome, merely updating the label of which allele is ancestral and which is derived. Thus, outgroup polarization is fully compatible with phased haplotype scans and no data is lost.

---

### R Code Exercise 7: Running Rsb in R using `rehh`
Write the R code necessary to:
1. Load the polarized Small Dogs and Large Dogs VCFs into `rehh` using `data2haplohh()`.
2. Compute the haplotype homozygosity scan for both populations using `scan_hh()`.
3. Compute the cross-population Rsb statistic using `ines2rsb()`.
4. Plot the Rsb Manhattan plot using `ggplot2`, highlighting the *IGF1* region between 41 Mb and 45.5 Mb in red.

**Write your R code here:**
```r
# 



library(rehh)
library(dplyr)
library(ggplot2)
library(R.utils)

# ------------------------------------------------------------
# 1. Load polarized VCFs
# ------------------------------------------------------------

hap_small <- data2haplohh(
  hap_file = "Part_2_CanidDiversity/small_dogs_polarized.vcf.gz",
  haplotype.in.columns = FALSE,
  polarize_vcf = FALSE
)

hap_large <- data2haplohh(
  hap_file = "Part_2_CanidDiversity/large_dogs_polarized.vcf.gz",
  haplotype.in.columns = FALSE,
  polarize_vcf = FALSE
)

# ------------------------------------------------------------
# 2. Haplotype homozygosity scans
# ------------------------------------------------------------

scan_small <- scan_hh(hap_small)
scan_large <- scan_hh(hap_large)

# ------------------------------------------------------------
# 3. Rsb: small dogs compared with large dogs
# ------------------------------------------------------------

rsb <- ines2rsb(scan_small, scan_large)

# In this rehh version, rsb is already a data frame
head(rsb)
colnames(rsb)

# ------------------------------------------------------------
# 4. Manhattan plot highlighting IGF1 region
# ------------------------------------------------------------

rsb_plot <- rsb %>%
  filter(!is.na(RSB)) %>%
  mutate(
    POS_MB = POSITION / 1e6,
    Region = ifelse(
      POS_MB >= 41 & POS_MB <= 45.5,
      "IGF1 region",
      "Other regions"
    )
  )

ggplot(rsb_plot, aes(x = POS_MB, y = RSB, color = Region)) +
  geom_point(size = 1, alpha = 0.7) +
  scale_color_manual(
    values = c(
      "Other regions" = "grey50",
      "IGF1 region" = "red"
    )
  ) +
  labs(
    title = "Rsb scan: Small dogs vs Large dogs",
    x = "Position on chromosome 15 (Mb)",
    y = "Rsb",
    color = NULL
  ) +
  theme_classic() +
  theme(
    legend.position = "top",
    axis.title = element_text(face = "bold")
  )

```

---

### Questions for Students
1. **Explain the physical and mathematical difference between Rsb and XP-nSL. Why does Rsb show a much higher, less noisy peak at the *IGF1* locus in dogs compared to XP-nSL?**
   * *Answer*:
2. **Does outgroup polarization of a phased VCF file corrupt or destroy the phasing information? Why or why not?**
   * *Answer*:
