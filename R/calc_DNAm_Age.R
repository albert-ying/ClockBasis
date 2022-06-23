#-----------------------------------------------------------------------------
#' Calculation function for mouse DNAm age
#' @importFrom ggthemes geom_rangeframe theme_tufte
#' @importFrom ggplot2 geom_point ggplot_build scale_fill_viridis_c scale_color_viridis_c guides coord_cartesian update_geom_defaults margin
#' @importFrom hrbrthemes theme_ipsum
#' @importFrom grid unit
#' @importFrom ggtext element_markdown
#' @importFrom ggsci scale_color_npg scale_fill_npg
#' @export
#-----------------------------------------------------------------------------

# do_dnam_clock_human = function(data, ...) {
#   return(NULL)
# }

#-----------------------------------------------------------------------------
#' Calculation function for mouse DNAm age
#' @importFrom ggthemes geom_rangeframe theme_tufte
#' @importFrom ggplot2 geom_point ggplot_build scale_fill_viridis_c scale_color_viridis_c guides coord_cartesian update_geom_defaults margin
#' @importFrom hrbrthemes theme_ipsum
#' @importFrom grid unit
#' @importFrom ggtext element_markdown
#' @importFrom ggsci scale_color_npg scale_fill_npg
#' @export
#-----------------------------------------------------------------------------

# do_dnam_clock_mouse = function(data, ...) {

# }

#-----------------------------------------------------------------------------
#' Global calculation function for DNAm age
#' @importFrom ggthemes geom_rangeframe theme_tufte
#' @importFrom ggplot2 geom_point ggplot_build scale_fill_viridis_c scale_color_viridis_c guides coord_cartesian update_geom_defaults margin
#' @importFrom hrbrthemes theme_ipsum
#' @importFrom grid unit
#' @importFrom ggtext element_markdown
#' @importFrom ggsci scale_color_npg scale_fill_npg
#' @export
#-----------------------------------------------------------------------------

# do_dnam_clock = function(data, species = c("human", "mouse"), ...) {
#   if (species == "human") {
#     do_dnam_clock_human(data, ...)
#   } else if (species == "mouse") {
#     do_dnam_clock_mouse(data, ...)
#   } else {
#     stop("Invalid species: only 'human' or 'mouse' are allowed")
#   }
# }

# debug
#-----------------------------------------------------------------------------
if (FALSE) {
  library(ggplot2)
  library(dplyr)
  library(ggRetro)
  library(ohmyggplot)
  library(ggplot2pipes)
  library(paletteer)
  library(hrbrthemes)
  library(ggtext)

  ohmyggplot::oh_my_ggplot()

  {
  oh_my_ggplot()
  tibble(X = c(1:16), y = rep(10,16)) |>
    mutate(X = as.character(X)) |>
    ggplot(aes(X, y)) +
    geom_col(aes(fill = X))
  }

  iris |>
    ggplot(aes(Petal.Length, Petal.Width, fill = Species)) +
    geom_point() +
    geom_smooth(aes(color = Species))

  cars |>
    tibble() |>
    mutate(type = rep(1:10, 5)) |>
    mutate(type = as.character(type)) |>
    ggplot(aes(speed, dist, fill = type)) +
    geom_point(alpha = 0.8)
    # geom_text(data = annot_tb, aes(x, y, label = lab))
  p + better_fill_legend + theme(legend.position = "top")
  p |>
    base_mode()
  p |>
    base_facet("am", guides = "auto", nrow = 2)
  
  base_facet(p2,"am")
}
