# Code to create data for DLM example
# Written by Ander Wilson

rm(list = ls())
library(tidyverse)
library(elevatr)
library(tigris)
library(sf)
library(table1)


#--------------------------------------------------------------
# read in CO birth data

co_birth <- read_csv("~/CAPBOS/CDPHE_data/co_birth/co_birth_0718.csv")

total <- nrow(co_birth)


#--------------------------------------------------------------
# add elevation (new way)
tracts <- tracts(state = "CO", year = 2010) %>%
  select(GEOID = GEOID10, geometry)
cents <- st_centroid(tracts)
elevations <- get_elev_point(cents, units = "feet")
elevations <- elevations %>%
  data.frame() %>%
  mutate(
    fipscoor = substr(GEOID, 3, 5),
    tract = substr(GEOID, 6, 11),
    elev_feet_tract = elevation
  ) %>%
  select(fipscoor, tract, elev_feet_tract)
co_birth <- left_join(x = co_birth, y = elevations, by = c("fipscoor", "tract"))


#--------------------------------------------------------------
# Filter based on:
# estimated gestational age at birth
# complete sex data
# year of conception
# singleton births
# front range counties
# tract elevation <= 6000

co_birth <- co_birth %>%
  filter(
    #EstGest>=37,     # limit to full term
    Sex %in% c("M", "F"), # limit to defined sex
    fipscoor != "999",
    !is.na(tract),
    !is.na(EstGest),
    !is.na(BWGr), # limit to births with a birth weight
    elev_feet_tract <= 6000,
    Plurality == 1
  ) %>% # limit singleton births
  mutate(
    DOB = lubridate::as_date(DOB, format = "%m/%d/%Y"), # convert date to date format
    DOC = DOB - 7 * EstGest,
    YOC = lubridate::year(DOC),
    MOC = lubridate::month(DOC)
  ) %>% #calculate date of conception
  filter(YOC %in% 2007:2017) # limit cohort based on conception year


# front range metro counties
fr <-
  c(
    "069", # larimer
    "123", # weld
    "013", # Boulder
    "031", # city of denver
    "005", # Arapahoe
    "059", # Jefferson County
    "001", # Adams
    "035", # Douglas
    "014", # City and County of Broomfield
    "039", # Elbert County
    "093", # Park County
    "019", # Clear Creek County
    "047", # Gilpin County
    "041", # El Paso County
    "119", # Teller County
    "043", # Fremont County
    "101" # Pueblo County
  )

co_birth <- co_birth %>%
  filter(fipscoor %in% fr)

# remove duplicate vsid
co_birth <- co_birth[!duplicated(co_birth$vsid), ]


#--------------------------------------------------------------
# Select Variables

co_birth <- co_birth %>%
  select(
    vsid,
    DOC,
    BWGr,
    EstGest,
    Sex,
    Plurality,
    MOC,
    YOC,
    MatAge,
    Marital,
    Income,
    MEduc,
    PrenatalCare,
    MotherHeight,
    PriorWeight,
    methnic,
    mracebrg,
    CigsPrePreg,
    CigsTri1,
    CigsTri2,
    CigsTri3,
    fipscoor,
    tract,
    elev_feet_tract
  )

dim(co_birth)


#--------------------------------------------------------------
# Construct analysis variables

co_birth <- co_birth %>%
  mutate(
    hispanic = ifelse(
      methnic >= 200 & methnic < 300,
      "Hispanic",
      "NonHispanic"
    ),
    race = ifelse(
      mracebrg == "01",
      "white",
      ifelse(
        mracebrg == "02",
        "Black",
        ifelse(mracebrg == "03", "AmInd", "AsianPI")
      )
    ),
    PriorWeight = ifelse(PriorWeight %in% c("?", "999", "NA"), NA, PriorWeight),
    PriorWeight = as.numeric(PriorWeight),
    PriorWeight = as.numeric(PriorWeight),
    Income = letters[Income + 1],
    MEduc = letters[MEduc + 1]
  )
co_birth <- co_birth %>%
  mutate(
    Income = recode(
      Income,
      a = "<15k",
      b = "15_24k",
      c = "25_34k",
      d = "35_49k",
      e = "50_74k",
      f = ">75k"
    ),
    MEduc = recode(
      MEduc,
      a = "lsHS",
      b = "lsHS",
      c = "HSdeg",
      d = "HSdeg",
      e = "AssocDeg",
      f = "CollegeDeg",
      g = "AdvDeg",
      h = "AdvDeg",
      i = "AdvDeg"
    ),
    CigsPrePreg = ifelse(CigsPrePreg %in% c("?", NA), NA, CigsPrePreg),
    CigsPrePreg = as.numeric(CigsPrePreg),
    CigsTri1 = ifelse(CigsTri1 %in% c("?", NA), NA, CigsTri1),
    CigsTri1 = as.numeric(CigsTri1),
    CigsTri2 = ifelse(CigsTri2 %in% c("?", NA), NA, CigsTri2),
    CigsTri2 = as.numeric(CigsTri2),
    CigsTri3 = ifelse(CigsTri3 %in% c("?", NA), NA, CigsTri3),
    CigsTri3 = as.numeric(CigsTri3),
    SmkPre = ifelse(
      is.na(CigsPrePreg),
      NA,
      c("N", "Y")[1 + 1 * (0 < CigsPrePreg)]
    ),
    SmkTri1 = ifelse(is.na(CigsTri1), NA, c("N", "Y")[1 + 1 * (0 < CigsTri1)]),
    SmkTri2 = ifelse(is.na(CigsTri2), NA, c("N", "Y")[1 + 1 * (0 < CigsTri2)]),
    SmkTri3 = ifelse(is.na(CigsTri3), NA, c("N", "Y")[1 + 1 * (0 < CigsTri3)]),
    SmkAny = ifelse(
      is.na(SmkTri1) & is.na(SmkTri2) & is.na(SmkTri3),
      NA,
      ifelse(SmkTri1 == "Y" | SmkTri2 == "Y" | SmkTri3 == "Y", "Y", "N")
    ),
    Marital2 = Marital,
    Marital2 = ifelse(Marital == "MR", "CM", Marital2),
    Marital2 = ifelse(Marital == "D", "DMSW", Marital2),
    Marital2 = ifelse(Marital == "MS", "DMSW", Marital2),
    Marital2 = ifelse(Marital == "W", "DMSW", Marital2),
    Marital2 = ifelse(Marital == "U", NA, Marital2),
    Marital2 = ifelse(Marital == "US", NA, Marital2)
  )

# calculate BMI and clean height
co_birth <- co_birth %>%
  mutate(
    MotherHeight = ifelse(
      MotherHeight == "?",
      NA,
      ifelse(
        substring(MotherHeight, 3, 3) == ":",
        MotherHeight,
        paste(
          substring(MotherHeight, 1, 2),
          substring(MotherHeight, 3, 4),
          sep = ":"
        )
      )
    ),
    MotherHeightIn = ifelse(
      is.na(MotherHeight),
      NA,
      as.numeric(substring(MotherHeight, 1, 2)) *
        12 +
        as.numeric(substring(MotherHeight, 4, 5))
    ),
    MotherBMI = (PriorWeight / 2.2) / (MotherHeightIn * 0.0254)^2
  )

# alternative smoking variable
co_birth$Smk <- ifelse(is.na(co_birth$SmkAny), NA, "Never")
co_birth$Smk <- ifelse(co_birth$CigsPrePreg > 0, "Former", co_birth$Smk)
co_birth$Smk <- ifelse(
  co_birth$CigsTri1 > 0 | co_birth$CigsTri2 > 0,
  "CurLow",
  co_birth$Smk
)
co_birth$Smk <- ifelse(
  co_birth$CigsTri1 > 10 | co_birth$CigsTri2 > 10,
  "CurHigh",
  co_birth$Smk
)


# drop variables used to create items above
co_birth <- co_birth %>%
  select(
    vsid,
    DOC,
    BWGr,
    EstGest,
    Sex,
    Plurality,
    MOC,
    YOC,
    MatAge,
    Marital2,
    Income,
    MEduc,
    PrenatalCare,
    MotherHeightIn,
    PriorWeight,
    MotherBMI,
    hispanic,
    race,
    Smk,
    fipscoor,
    tract,
    elev_feet_tract
  )

dim(co_birth)
co_birth %>% drop_na()


#--------------------------------------------------------------
# Add PM2.5 from CMAQ

load(
  "~/CAPBOS/ExposureData/CMAQ_ozone_pm/cmaqdownscaler_colorado_tract_pm25_weeklag.rda"
)

co_birth <- left_join(
  x = co_birth,
  y = cmaqdownscaler_colorado_tract_pm25_weeklag,
  by = join_by(DOC == Date, fipscoor, tract)
)


#--------------------------------------------------------------
# Add temperature from gridmet

load("~/CAPBOS/ExposureData/gridmet_temp/gridmet_tmmx_tract_trimester.rda")


co_birth <- left_join(
  x = co_birth,
  y = gridmet_tmmx_tract_trimester,
  by = join_by(DOC == date, fipscoor, tract)
)


co_birth_dlm_data <- co_birth %>%
  # filter(EstGest <= 42 & EstGest >= 20) %>%
  select(-vsid, -DOC, -Plurality, -tract) %>%
  drop_na()
# hist(co_birth_dlm_overview %>% pull(EstGest), breaks=20:42)

# save file
save(
  co_birth_dlm_data,
  file = "Data/co_dlm_data.rda"
)
