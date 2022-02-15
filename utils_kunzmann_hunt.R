prepare_kunzmann_hunt = function(df){
  # compute new variables
  df$White = !(df$Asian | df$Black |df$HawaiianPacific | df$IndianAlaskan)
  df$smoke_ever = df$smoke_current | df$smoke_former
  df = df[, c("ID", "CaseControl", "smoke_current", "smoke_former", "Gender", 
                         "GerdAtIndex", "White", "smoke_ever", "bmi", "ageatindex",
              "h2r_max", "ppi_max")]
  df$k_age_bin = 
    relevel(cut(df$ageatindex, breaks=c(0, 50, 55, 60, 65, 100)), ref="(50,55]")
  df$k_bmi_bin = cut(df$bmi, breaks=c(0, 25, 30, 35, 100))
  df$h_age_bin = cut(df$ageatindex, breaks=c(0, 50, 60, 70, 100))
  df$h_bmi_bin = cut(df$bmi, breaks=c(0, 30, 100))
  df$h_smoke_any = pmax(
    df$smoke_former,
    df$smoke_current
  )
  df$k_ec = pmax(
    df$h2r_max>0, 
    df$ppi_max>0, 
    df$GerdAtIndex
  ) > 0
  return(df)
}

get_complete_cases = function(df){
  df = prepare_kunzmann_hunt(df)
  complete_cases = !(
    is.na(df$Gender) |
      is.na(df$ageatindex) |
      is.na(df$smoke_former) |
      is.na(df$smoke_current) |
      is.na(df$bmi) |
      is.na(df$White) |
      is.na(df$GerdAtIndex) |
      is.na(df$h2r_max) |
      is.na(df$ppi_max)
  )
  return(df %>% subset(complete_cases) %>% pull(ID) )
}

kunzmann_score = function(df){
  df = prepare_kunzmann_hunt(df)
  score = 0
  score = score + (df$k_age_bin == "(55,60]") * 1.5
  score = score + (df$k_age_bin == "(60,65]") * 2.5
  score = score + (df$k_age_bin == "(65,100]") * 3.5
  score = score + df$Gender * 4
  score = score + (df$k_bmi_bin == "(25,30]") * 1
  score = score + (df$k_bmi_bin == "(30,35]") * 1.5
  score = score + (df$k_bmi_bin == "(35,100]") * 2.5
  score = score + df$smoke_former * 2.5
  score = score + df$smoke_current * 3.5
  score = score + df$k_ec * 1.5
  return(score)
}

hunt_score = function(df){
  df = prepare_kunzmann_hunt(df)
  score = 3.6
  score = score * ifelse(df$Gender, 1.9, 1.)
  score = score * ifelse(df$h_age_bin == "(50,60]", 2.1, 1.)
  score = score * ifelse(df$h_age_bin == "(60,70]", 3.2, 1.)
  score = score * ifelse(df$h_age_bin == "(70,100]", 3.1, 1.)
  score = score * ifelse(df$h_bmi_bin == "(30,100]", 1.8, 1.)
  score = score * ifelse(df$GerdAtIndex, 3.7, 1.)
  score = score * ifelse(df$h_smoke_any, 2.1, 1.)
  return(score/100000)
}