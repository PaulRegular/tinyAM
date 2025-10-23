
#' TAM color palette helpers
#'
#' @description
#' Provides the TAM package's primary color palette as a vector of
#' hex codes and a convenience function for generating interpolated palettes.
#'
#' @return
#' * `tam_cols()` returns a character vector of hex color codes.
#' * `tam_pal(n)` returns a character vector of length `n` created by linearly
#'   interpolating between the base palette colors.
#'
#' @examples
#' tam_cols()
#' tam_pal(3)
#'
#' @export
#' @md
#' @seealso [grDevices::colorRampPalette()]
tam_cols <- function() {
  c(
    "#003e56",
    "#007187",
    "#028f8d",
    "#2aaa85",
    "#93c25c",
    "#fec124",
    "#ffa31d",
    "#fe8427",
    "#f05c28"
  )
}

#' @rdname tam_cols
#' @param n Number of colors to generate.
#' @export
tam_pal <- function(n) {
  grDevices::colorRampPalette(tam_cols())(n)
}
