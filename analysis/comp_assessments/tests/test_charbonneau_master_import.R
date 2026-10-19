source("analysis/comp_assessments/scripts/database/071_import_charbonneau_master.R")

master <- data.frame(
  Genus = c("Gadus", "Clupea", "Atheresthes", "Merlangius"),
  Species = c("morhua", "harengus", "stomias", "merlangus"),
  Assessor_ID = c("DFO", "DFO-GSL", "AFSC", "ICES-WGSCE"),
  unique_ID = c(
    "DFO_morhua_2J3KL",
    "DFO-GSL_harengus_4T-Spring",
    "AFSC_stomias_BSAI-Females",
    "ICES-WGCSE_merlangus_CS7a"
  ),
  PDF.file = c("cod.pdf", "herring.pdf", "arrowtooth.pdf", "whiting.pdf"),
  stringsAsFactors = FALSE
)

stocks <- data.frame(
  stock_id = c("dfo_cod_2j3kl", "afsc_arrowtooth_flounder_bsai"),
  charbonneau_id = c("DFO_2J3KL_Gadus_morhua", ""),
  authority = c("DFO", "NOAA-AFSC"),
  authority_stock_id = c("2J3KL", "BSAI arrowtooth"),
  scientific_name = c("Gadus morhua", "Atheresthes stomias"),
  common_name = c("Northern cod", "Arrowtooth flounder"),
  area = c("2J3KL", "BSAI"),
  region = c("Newfoundland and Labrador", "Alaska"),
  ocean = c("Northwest Atlantic", "Northeast Pacific"),
  notes = c("Existing stock record.", "Existing aggregate stock."),
  stringsAsFactors = FALSE
)

matches <- data.frame(
  master_id = "DFO_morhua_2J3KL",
  stock_id = "dfo_cod_2j3kl",
  reason = "scientific name and stock area",
  stringsAsFactors = FALSE
)
overlaps <- data.frame(
  master_id = "AFSC_stomias_BSAI-Females",
  stock_id = "afsc_arrowtooth_flounder_bsai",
  reason = "sex-specific record associated with the existing aggregate stock",
  stringsAsFactors = FALSE
)

result <- ck_prepare_master_import(master, stocks, matches = matches, overlaps = overlaps)
stopifnot(
  nrow(result$stocks) == nrow(stocks) + 3L,
  result$matched == 1L,
  result$added == 3L,
  all(vapply(master$unique_ID, function(id) {
    sum(result$stocks$charbonneau_id == id, na.rm = TRUE) == 1L
  }, logical(1))),
  result$stocks$charbonneau_id[result$stocks$stock_id == "dfo_cod_2j3kl"] ==
    "DFO_morhua_2J3KL",
  grepl("Previous Charbonneau identifier retained as alias: DFO_2J3KL_Gadus_morhua",
    result$stocks$notes[result$stocks$stock_id == "dfo_cod_2j3kl"], fixed = TRUE),
  result$stocks$area[result$stocks$charbonneau_id == "DFO-GSL_harengus_4T-Spring"] ==
    "4T-Spring",
  result$stocks$area[result$stocks$charbonneau_id == "AFSC_stomias_BSAI-Females"] ==
    "BSAI-Females",
  result$stocks$authority[result$stocks$charbonneau_id == "AFSC_stomias_BSAI-Females"] ==
    "NOAA-AFSC",
  result$stocks$ocean[result$stocks$charbonneau_id == "AFSC_stomias_BSAI-Females"] ==
    "Northeast Pacific",
  grepl("Potential overlap with existing stock afsc_arrowtooth_flounder_bsai",
    result$stocks$notes[result$stocks$charbonneau_id == "AFSC_stomias_BSAI-Females"], fixed = TRUE),
  grepl("Sex-specific catalogue entry", result$stocks$notes[
    result$stocks$charbonneau_id == "AFSC_stomias_BSAI-Females"], fixed = TRUE),
  grepl("Master unique_ID prefix differs from Assessor_ID ICES-WGSCE",
    result$stocks$notes[result$stocks$charbonneau_id == "ICES-WGCSE_merlangus_CS7a"], fixed = TRUE),
  grepl("Assessor_ID=ICES-WGSCE",
    result$stocks$notes[result$stocks$charbonneau_id == "ICES-WGCSE_merlangus_CS7a"], fixed = TRUE),
  all(result$stocks$ocean %in% c(
    "Northwest Atlantic", "Northeast Atlantic", "Northeast Pacific", "Northwest Pacific"
  )),
  all(result$stocks$common_name[result$stocks$stock_id %in% result$mapping$stock_id[
    result$mapping$mapping == "new metadata row"]] == ""),
  all(result$stocks$authority_stock_id[result$stocks$stock_id %in% result$mapping$stock_id[
    result$mapping$mapping == "new metadata row"]] == "")
)

rerun <- ck_prepare_master_import(
  master, result$stocks, matches = matches, overlaps = overlaps
)
stopifnot(identical(result$stocks, rerun$stocks), rerun$matched == nrow(master), rerun$added == 0L)

bad_matches <- rbind(matches, data.frame(
  master_id = "not_in_master", stock_id = "dfo_cod_2j3kl", reason = "stale",
  stringsAsFactors = FALSE
))
missing_crosswalk_error <- tryCatch({
  ck_prepare_master_import(master, stocks, matches = bad_matches, overlaps = overlaps)
  FALSE
}, error = function(e) grepl("absent from Master.csv", conditionMessage(e), fixed = TRUE))
stopifnot(missing_crosswalk_error)

cat("Charbonneau-Keith catalogue mapping, metadata safeguards and idempotency passed.\n")
