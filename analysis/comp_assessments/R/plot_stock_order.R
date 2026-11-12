# Higher groups follow Betancur-R. et al. (2017), doi:10.1186/s12862-017-0958-3.
# Additional genera/families checked against Eschmeyer's Catalog of Fishes (2026):
# https://researcharchive.calacademy.org/research/ichthyology/catalog/SpeciesByFamily.asp
# Species records: Molva (Lotidae), Argentina (Argentinidae), Anarhichas
# (Anarhichadidae), Anoplopoma (Anoplopomatidae); flatfishes remain Pleuronectidae.
# Scientific names break ties within groups; stocks of a species stay together.
assessment_plot_stocks <- function(assessments, stocks, assessment_ids) {
  groups <- list(
    Clupeiformes = c("Clupea", "Sprattus", "Engraulis"),
    Gadidae = c("Gadus", "Melanogrammus", "Merlangius", "Micromesistius",
                "Pollachius", "Trisopterus"),
    Merlucciidae = "Merluccius",
    Phycidae = "Urophycis",
    Scombriformes = "Scomber",
    Carangidae = "Trachurus",
    Cyclopsettidae = "Paralichthys",
    Scophthalmidae = "Scophthalmus",
    Pleuronectidae = c("Atheresthes", "Eopsetta", "Glyptocephalus", "Hippoglossoides", "Lepidopsetta", "Limanda", "Pleuronectes", "Pseudopleuronectes",
                      "Reinhardtius", "Myzopsetta", "Microstomus"),
    Soleidae = "Solea",
    Sparidae = "Stenotomus",
    Uranoscopiformes = "Ammodytes",
    Moronidae = c("Morone", "Dicentrarchus"),
    Sebastidae = "Sebastes",
    Hexagrammidae = "Pleurogrammus",
    Argentinidae = "Argentina",
    Lotidae = "Molva",
    Anarhichadidae = "Anarhichas",
    Anoplopomatidae = "Anoplopoma"
  )
  group_rank <- setNames(rep(seq_along(groups), lengths(groups)),
                        unlist(groups, use.names = FALSE))
  ids <- unique(as.character(assessment_ids))
  if (anyNA(ids) || any(!ids %in% assessments$assessment_id)) {
    cli::cli_abort("Plot assessments must have assessment metadata.")
  }
  selected <- assessments[match(ids, assessments$assessment_id), c("assessment_id", "stock_id")]
  metadata <- stocks[match(selected$stock_id, stocks$stock_id),
                     c("common_name", "scientific_name", "ocean")]
  info <- cbind(selected, metadata)
  rownames(info) <- NULL
  info$stock <- ifelse(is.na(info$common_name) | !nzchar(info$common_name),
                       info$assessment_id, info$common_name)
  info$stock <- sub(" \\(formerly.*\\)$", "", info$stock)
  info$basin <- factor(info$ocean, levels = c("Northwest Atlantic", "Northeast Atlantic",
                                            "Northeast Pacific", "Northwest Pacific"))
  if (anyNA(info$basin)) {
    cli::cli_abort("Plot stocks must have a classified ocean basin.")
  }
  genus <- sub(" .*", "", info$scientific_name)
  rank <- unname(group_rank[genus])
  if (anyNA(rank)) {
    missing <- unique(info$scientific_name[is.na(rank)])
    cli::cli_abort("Add the missing plot taxonomy for: {paste(missing, collapse = ', ')}.")
  }
  info <- info[order(info$basin, rank, info$scientific_name, info$stock,
                     info$assessment_id, method = "radix"), ]

  # Add local geography where the stock name is generic or uses an area code alone.
  labels <- c(
    ck_dfo_aeglefinus_4x5y = "4X5Y haddock (western Scotian Shelf / Bay of Fundy)",
    ck_dfo_platessoides_4t = "Southern Gulf of St. Lawrence American plaice (4T)",
    ck_nefsc_saw_chrysops_usatl = "Mid-Atlantic / New England scup",
    nafo_cod_3m = "Flemish Cap cod (3M)"
  )
  label <- unname(labels[info$stock_id])
  use_label <- !is.na(label)
  info$stock[use_label] <- label[use_label]
  info
}
