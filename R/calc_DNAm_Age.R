#-----------------------------------------------------------------------------
#' Load mouse clock internal objects
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

  mouse_clock_intercept.ls = list(
    Meer_intercept,
    Petkovich_intercept,
    Thompson_intercept,
    Wang_intercept
  )

  names(mouse_clock.ls) = names(mouse_clock_intercept.ls) = c(
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
#'
#' @param data A dataframe with the chr_pos id as first column. The rest of the columns should be methylation levels for each sample, with sample names as column names.
#' @param clock A string vector with the names of the clocks to use.
#' @param fuzzy_pos_window A interger indicating a window size that to be used to infer missing methylation value. Default is 100.
#' @importFrom dplyr left_join select pull arrange group_by summarize mutate bind_rows ungroup bind_cols
#' @importFrom tidyr extract
#' @importFrom tibble rownames_to_column
#' @importFrom purrr map_dfr
#' @importFrom fuzzyjoin genome_left_join
#' @return An tibble with sample per row and DNAm age as columns.
#' @examples
#' do_dnam_clock_mouse(df, c("MeerClock", "PetkovichClock", "ThompsonClock", "WangLiver"), 100)
#' @export
#-----------------------------------------------------------------------------

do_dnam_clock_mouse = function(
  data,
  clock = c("MeerClock", "PetkovichClock", "ThompsonClock", "WangLiver"),
  fuzzy_pos_window = 100
) {
  require(dplyr)
  require(fuzzyjoin)
  require(tidyr)
  require(tibble)
  require(purrr)

  colnames(data)[1] = "CpG"
  if (max(pull(data[,2]), na.rm = T) <= 1) {
    message("The data is at 0-1 scale, transforming it to 0-100 scale")
    data = data |>
      mutate_if(is.numeric, \(x){x * 100})
  }
  clock.ls = mouse_clock.ls[clock]
  clock_intercept.ls = mouse_clock_intercept.ls[clock]
  imputer.ls = mouse_imputer.ls[clock]
  map_dfr(2:ncol(data), ~{
    sample_name = colnames(data)[.x]
    message(paste0("Processing ", sample_name))
    me_mat = cbind(data[,1], data[,.x])
    colnames(me_mat) = c("CpG", "level")
    res.ls = list()
    for (i in 1:length(clock.ls)) {
      print(i)
      clock = clock.ls[[i]]
      clock_name = names(clock.ls)[i]
      ## first match
      raw = clock |>
        select(CpG) |>
        left_join(me_mat, by = "CpG") 
      missed = raw[is.na(raw$level),]
      ## 1st fuzzy match + - 3 bp
      missed_pos = missed |>
        extract(CpG, c("chr", "tar_bp"), '(chr.*)_(.*)', remove = F, convert = T) |>
        select(-level) |>
        mutate(start = tar_bp - 3, end = tar_bp + 3)
      me_mat2 = me_mat |>
        extract(CpG, c("chr", "bp"), '(chr.*)_(.*)', remove = T, convert = T) |>
        mutate(start = bp, end = bp)
      fuzzy_joined = genome_left_join(missed_pos, me_mat2, by = c("chr", 'start', 'end'))

      mean_fuzzy = fuzzy_joined |>
        select(CpG, level) |> 
        group_by(CpG) |>
        summarize(level = mean(level, na.rm = T)) |>
        ungroup() 
    
      raw = raw[!is.na(raw$level),] |>
        bind_rows(mean_fuzzy) |>
        arrange(CpG)

      missed = raw[is.na(raw$level),]

      ## 2nd fuzzy match + - 100 bp
      missed_pos = missed |>
        extract(CpG, c("chr", "tar_bp"), '(chr.*)_(.*)', remove = F, convert = T) |>
        select(-level) |>
        mutate(start = tar_bp - fuzzy_pos_window, end = tar_bp + fuzzy_pos_window)
      me_mat2 = me_mat |>
        extract(CpG, c("chr", "bp"), '(chr.*)_(.*)', remove = T, convert = T) |>
        mutate(start = bp, end = bp)
      fuzzy_joined = genome_left_join(missed_pos, me_mat2, by = c("chr", 'start', 'end'))

      mean_fuzzy = fuzzy_joined |>
        select(CpG, level) |> 
        group_by(CpG) |>
        summarize(level = mean(level, na.rm = T)) |>
        ungroup() 
    
      raw_mat = raw[!is.na(raw$level),] |>
        bind_rows(mean_fuzzy) |>
        arrange(CpG) |>
        select(-CpG) |>
        as.matrix()

      missingness = apply(raw_mat, 2, function(x) mean(is.na(x)))

      if (any(missingness > 0)) { ## impute missingness
        # Impute missing values
        imputer = imputer.ls[[i]]
        k = which(is.na(raw_mat), arr.ind = T)
        raw_mat[k] = imputer[k[,1]]
      }
      # Calc meAge
      clock_coef = clock$Weight
      # make table
      if (clock_name == 'MeerClock') {
        bage = (matrix(clock_coef, nrow = 1) %*% raw_mat) + clock_intercept.ls[[i]]
        bage = bage/30.417
      } else if (clock_name == 'PetkovichClock') {
        a = 0.1666
        b = 0.4185 
        c = -1.712
        bage = (matrix(clock_coef, nrow = 1) %*% raw_mat/100) + clock_intercept.ls[[i]]
        bage = ((((bage - c)/a)**(1/b))/30.417)

      } else if (clock_name == 'ThompsonClock') {
        bage = (matrix(clock_coef, nrow = 1) %*% raw_mat/100) + clock_intercept.ls[[i]]
      } else if (clock_name == 'WangLiver') {
        bage = (matrix(clock_coef, nrow = 1) %*% raw_mat/100) + clock_intercept.ls[[i]]
        bage = (2**(bage))/30.417
      }
      res = data.frame(Age = bage[1,], missing = missingness)
      colnames(res) = c(clock_name, paste0(clock_name, "_missing"))
      res.ls[[i]] = res
    }
    res_mat = bind_cols(res.ls) |>
      as.data.frame() |>
      rownames_to_column("sample_id") |>
      mutate(sample_id = sample_name)
    return(res_mat)
  })
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
  remotes::install_github("albert-ying/ClockBasis") 
  library(ClockBasis)
  library(tidyverse)
  do_dnam_clock_mouse(tibble(a = "chr1_1000", lel = 0.5, kkk = 1))
}
