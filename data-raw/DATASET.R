## code to prepare `DATASET` dataset goes here
suppressPackageStartupMessages(library(spatstat.random))
suppressPackageStartupMessages(library(spatstat.geom))
suppressPackageStartupMessages(library(spatstat.explore))
suppressPackageStartupMessages(library(tidyverse))
suppressPackageStartupMessages(library(tictoc))
suppressPackageStartupMessages(library(VectraPolarisData))
suppressPackageStartupMessages(library(tictoc))
suppressPackageStartupMessages(library(purrr))
set.seed(3000)


# Load the data
sample_index = 1

# load data from VectraPolarisData package and process here
spe_ovarian <- HumanOvarianCancerVP()

## Assays slots
assays_slot <- assays(spe_ovarian)
intensities_df <- assays_slot$intensities
nucleus_intensities_df<- assays_slot$nucleus_intensities
rownames(nucleus_intensities_df) <- paste0("nucleus_", rownames(nucleus_intensities_df))
membrane_intensities_df<- assays_slot$membrane_intensities
rownames(membrane_intensities_df) <- paste0("membrane_", rownames(membrane_intensities_df))

# colData and spatialData
colData_df <- colData(spe_ovarian)
spatialCoords_df <- spatialCoords(spe_ovarian)

# clinical data
patient_level_ovarian <- metadata(spe_ovarian)$clinical_data %>%
  # create binary stage variable
  dplyr::mutate(stage_bin = ifelse(stage %in% c("1", "2"), 0, 1))

ovarian <- as.data.frame(cbind(colData_df,
                               spatialCoords_df,
                               t(intensities_df),
                               t(nucleus_intensities_df),
                               t(membrane_intensities_df))
) %>%
  dplyr::rename(cd19 = cd19_opal_480,
                cd68 = cd68_opal_520,
                cd3 = cd3_opal_540,
                cd8 = cd8_opal_650,
                ier3 = ier3_opal_620,
                pstat3 = p_stat3_opal_570,
                ck = ck_opal_780,
                ki67 = ki67_opal_690,
                x = cell_x_position,
                y = cell_y_position) %>%
  dplyr::select(contains("id"), tissue_category, x,y, contains("phenotype"), ck:dapi) %>%
  # remove control subjects
  dplyr::filter(sample_id %in% patient_level_ovarian$sample_id) %>%
  mutate(immune = ifelse(phenotype_cd19 == "CD19+" | phenotype_cd8 == "CD8+" |
                           phenotype_cd3 == "CD3+" | phenotype_cd68 == "CD68+", "immune", "background"),
         immune = factor(immune, levels = c("immune", "background"))) %>%
  select(cell_id, sample_id, x, y, immune, tissue_category, everything())

rm(spe_ovarian, assays_slot, intensities_df, nucleus_intensities_df, membrane_intensities_df, colData_df, spatialCoords_df)


ovarian = ovarian %>% filter(tissue_category == "Tumor")
sample_index = 1
ids = unique(ovarian$sample_id)

ovarian = ovarian %>%
  filter(sample_id == ids[sample_index])

marksvar = "immune"

w = convexhull.xy(ovarian[["x"]], ovarian[["y"]])
pp_obj_ovarian = ppp(ovarian[["x"]], ovarian[["y"]], window = w, marks = ovarian[[marksvar]])

# Save the data


usethis::use_data(pp_obj_ovarian, overwrite = TRUE)
