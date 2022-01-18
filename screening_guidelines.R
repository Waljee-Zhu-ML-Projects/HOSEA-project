screening = function(df){
  df$White = !(df$Asian | df$Black | df$HawaiianPacific | df$IndianAlaskan)
  df$smoke_ever = df$smoke_current | df$smoke_former
  screenings = dplyr::mutate(df, 
    ACG2016     = Gender & GerdAtIndex & ( (ageatindex>=50) + (White) + (bmi>30) + (smoke_ever) > 1),
    ACG2022     = GerdAtIndex & ( (Gender) & (ageatindex>=50) + (White) + (bmi>30) + (smoke_ever) > 2),
    ACP2012     = Gender & GerdAtIndex & (ageatindex>=50) & ( (bmi>30) + (smoke_ever) > 0),
    AGA2011     = ( (GerdAtIndex) & (Gender) & (ageatindex>=50) + (White) + (bmi>30) > 1),
    AGA_CPU2022 = ( (GerdAtIndex) & (Gender) & (ageatindex>=50) + (White) + (bmi>30) + (smoke_ever) > 2),
    ASGE2019    = GerdAtIndex & ( (Gender) & (ageatindex>=50) + (bmi>30) + (smoke_ever) > 0),
    BSG2013     = GerdAtIndex & ( (Gender) & (ageatindex>=50) + (White) + (bmi>30) > 2),
    ESGE2017    = GerdAtIndex & ( (Gender) & (ageatindex>=50) + (White) + (bmi>30) > 1)
  )[, -(2:243)]
  return(screenings)
}