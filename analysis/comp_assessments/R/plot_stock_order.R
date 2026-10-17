# Higher groups follow Betancur-R. et al. (2017), doi:10.1186/s12862-017-0958-3.
# Scientific names break ties within groups; stocks of a species stay together.
assessment_plot_stocks <- function(assessments, stocks, assessment_ids) {
  groups <- list(
    Clupeiformes = c("Clupea", "Sprattus"),
    Gadidae = c("Gadus", "Melanogrammus", "Merlangius", "Micromesistius",
                "Pollachius", "Trisopterus"),
    Merlucciidae = "Merluccius",
    Scombriformes = "Scomber",
    Cyclopsettidae = "Paralichthys",
    Pleuronectidae = c("Atheresthes", "Lepidopsetta", "Limanda", "Pleuronectes"),
    Sebastidae = "Sebastes",
    Hexagrammidae = "Pleurogrammus"
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
  info[order(info$basin, rank, info$scientific_name, info$stock,
             info$assessment_id, method = "radix"), ]
}
