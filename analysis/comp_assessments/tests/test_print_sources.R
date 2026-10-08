root <- file.path("analysis", "comp_assessments")
pkgload::load_all(".", quiet = TRUE)
source(file.path(root, "R", "run_assessment.R"))

all_sources <- data.frame(
  assessment_url = "https://example.test/report%2F2026",
  framework_url = "https://example.test/framework",
  data_url = "https://example.test/data",
  model_url = "https://example.test/model",
  repository_url = "https://example.test/repository",
  stringsAsFactors = FALSE
)
all_lines <- print_sources(all_sources)
stopifnot(
  identical(all_lines[[1L]], "#### Assessment documentation"),
  all(vapply(c(
    "- [Assessment report](https://example.test/report%2F2026)",
    "- [Framework / methodology](https://example.test/framework)",
    "- [Assessment data](https://example.test/data)",
    "- [Native model / fitted run](https://example.test/model)",
    "- [Assessment repository](https://example.test/repository)"
  ), function(line) line %in% all_lines, logical(1)))
)

sparse_sources <- data.frame(
  assessment_url = NA_character_,
  framework_url = "",
  data_url = "   ",
  model_url = "https://example.test/only",
  repository_url = NA_character_,
  stringsAsFactors = FALSE
)
sparse_lines <- print_sources(sparse_sources)
stopifnot(length(sparse_lines) == 2L,
          identical(sparse_lines[[2L]],
                    "- [Native model / fitted run](https://example.test/only)"))

duplicate_sources <- data.frame(
  assessment_url = "https://example.test/shared",
  framework_url = "https://example.test/shared",
  data_url = NA_character_,
  model_url = "",
  repository_url = "https://example.test/repository",
  stringsAsFactors = FALSE
)
duplicate_lines <- print_sources(duplicate_sources)
stopifnot(length(duplicate_lines) == 3L,
          identical(duplicate_lines[[2L]],
                    "- [Assessment report / Framework / methodology](https://example.test/shared)"),
          sum(grepl("https://example.test/shared", duplicate_lines, fixed = TRUE)) == 1L)

no_sources <- data.frame(
  assessment_url = NA_character_,
  framework_url = "",
  data_url = NA_character_,
  model_url = " ",
  repository_url = "",
  stringsAsFactors = FALSE
)
stopifnot(identical(print_sources(no_sources), character()),
          identical(print_sources(data.frame(other = "value")),
                    character()))

cat("Assessment source-link tests passed.\n")
