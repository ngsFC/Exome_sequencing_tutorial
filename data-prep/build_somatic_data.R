# Builds the two files under data/somatic/ (Module 2) from the public GDC Data
# Portal (https://portal.gdc.cancer.gov/projects/TARGET-WT):
#   - every open-access "Masked Somatic Mutation" WXS file of TARGET-WT
#     (one MAF per tumor aliquot), merged into a single MAF;
#   - the clinical fields of the same cases (histology, age, sex, vital status).
#
# This script is NOT run during the tutorial -- it was run once to produce the
# files committed under data/somatic/. Requires: tidyverse, httr, jsonlite.

library(tidyverse)
library(httr)
library(jsonlite)

out_dir <- "data/somatic"
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
gdc <- "https://api.gdc.cancer.gov"

# ---- 1. File manifest ----
filters <- list(op = "and", content = list(
  list(op = "=", content = list(field = "cases.project.project_id", value = "TARGET-WT")),
  list(op = "=", content = list(field = "data_type", value = "Masked Somatic Mutation")),
  list(op = "=", content = list(field = "experimental_strategy", value = "WXS")),
  list(op = "=", content = list(field = "access", value = "open"))
))
resp <- POST(file.path(gdc, "files"), encode = "json",
             body = list(filters = filters, fields = "file_id,cases.submitter_id", size = 200, format = "JSON"))
hits <- content(resp, as = "parsed", simplifyVector = TRUE)$data$hits
cat("MAF files:", nrow(hits), "\n")

# ---- 2. Download and merge the per-aliquot MAFs ----
read_one_maf <- function(file_id) {
  path <- tempfile(fileext = ".maf.gz")
  GET(file.path(gdc, "data", file_id), write_disk(path, overwrite = TRUE))
  read_tsv(path, comment = "#", col_types = cols(.default = col_character()), progress = FALSE)
}
maf <- map_dfr(hits$file_id, read_one_maf)
cat("Mutations:", nrow(maf), "in", n_distinct(maf$Tumor_Sample_Barcode), "tumor aliquots\n")

write_tsv(maf, file.path(out_dir, "TARGET_WT_somatic_mutations.maf.gz"))

# ---- 3. Clinical data for the same cases ----
# Tumor_Sample_Barcode = TARGET-50-<case>-<sample><vial>-..., so the case ID is
# its first three dash-separated fields.
case_ids <- unique(str_extract(maf$Tumor_Sample_Barcode, "^[^-]+-[^-]+-[^-]+"))
resp_cases <- POST(file.path(gdc, "cases"), encode = "json", body = list(
  filters = list(op = "in", content = list(field = "submitter_id", value = case_ids)),
  fields  = "submitter_id,diagnoses.wilms_tumor_histologic_subtype,diagnoses.age_at_diagnosis,demographic.sex_at_birth,demographic.vital_status",
  size    = 200,
  format  = "JSON"
))
cases <- content(resp_cases, as = "parsed", simplifyVector = FALSE)$data$hits

clinical <- map_dfr(cases, function(h) {
  dx  <- h$diagnoses[[1]]
  dem <- h$demographic
  tibble(
    case_submitter_id     = h$submitter_id,
    histologic_subtype    = dx$wilms_tumor_histologic_subtype %||% NA_character_,
    age_at_diagnosis_days = as.numeric(dx$age_at_diagnosis %||% NA),
    sex_at_birth          = dem$sex_at_birth %||% NA_character_,
    vital_status          = dem$vital_status %||% NA_character_
  )
})

# One row per tumor aliquot, the unit maftools works with.
sample_sheet <- maf %>%
  distinct(Tumor_Sample_Barcode) %>%
  mutate(case_submitter_id = str_extract(Tumor_Sample_Barcode, "^[^-]+-[^-]+-[^-]+")) %>%
  left_join(clinical, by = "case_submitter_id") %>%
  arrange(case_submitter_id, Tumor_Sample_Barcode)

write_tsv(sample_sheet, file.path(out_dir, "TARGET_WT_clinical.tsv"))
cat("Wrote", file.path(out_dir, "TARGET_WT_somatic_mutations.maf.gz"), "and", file.path(out_dir, "TARGET_WT_clinical.tsv"), "\n")
print(count(sample_sheet, histologic_subtype))
