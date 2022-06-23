#-----------------------------------------------------------------------------
#' Calculation function for mouse DNAm age
#' @importFrom readr read_tsv
#' @importFrom purrr map
#' @importFrom dplyr arrange
#-----------------------------------------------------------------------------

{
  require(readr)
  require(purrr)
  require(dplyr)
  clock_folder = "data/Clocks/Mouse/"
# [1] "Csaba.tsv"                   "MeerClock.tsv"               "PetkovichClock.tsv"          "ThompsonClock-all.tsv"      
# [5] "ThompsonClock-conserved.tsv" "WangLiver.tsv" 
  meer = read_tsv(paste0(clock_folder, "MeerClock.tsv"))
  petkovich = read_tsv(paste0(clock_folder, "PetkovichClock.tsv"))
  thompson = read_tsv(paste0(clock_folder, "ThompsonClock-all.tsv"))
  wangliver = read_tsv(paste0(clock_folder, "WangLiver.tsv"))
  Meer_intercept = 234.64
  Petkovich_intercept = 0
  Thompson_intercept = 30.3172
  Wang_intercept = 5.827926399

  mouse_clock.ls = list(meer, petkovich, thompson, wangliver)

## reformate the clock
  mouse_clock.ls = map(mouse_clock.ls, ~{
    colnames(.x) = c("CpG", "Weight")
    arrange(.x, CpG)
  })

  clock_intercept.ls = list(
    Meer_intercept,
    Petkovich_intercept,
    Thompson_intercept,
    Wang_intercept
  )

  names(mouse_clock.ls) = c(
    "MeerClock",
    "PetkovichClock",
    "ThompsonClock",
    "WangLiver"
  )

  mouse_mean_me = read_csv("data/thompson_site_mean.csv")

  mouse_imputer.ls = list()
  for (i in 1:length(mouse_clock.ls)) {
    print(i)
    clock = mouse_clock.ls[[i]]
    mouse_imputer = clock |>
      select(CpG) |>
      arrange(CpG) |>
      left_join(mouse_mean_me, by = "CpG") |>
      pull(mean_level)
    mouse_imputer[is.na(mouse_imputer)] = 0.5
    mouse_imputer.ls[[i]] = mouse_imputer
    names(mouse_imputer.ls)[i] = names(mouse_clock.ls)[i]
  }
}

#-----------------------------------------------------------------------------
#' Calculation function for mouse DNAm age

#' @export
#-----------------------------------------------------------------------------

do_dnam_clock_mouse = function(data, ...) {
  return(mouse_imputer.ls)
}

#-----------------------------------------------------------------------------
#' Global calculation function for DNAm age
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
  # remotes::install_github("albert-ying/ClockBasis") 
  # library(ClockBasis)

}
