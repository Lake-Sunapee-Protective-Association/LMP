#*****************************************************************
#*      Cary Institute of Ecosystem Studies (Millbrook, NY)      *
#*                                                               *
#* TITLE:   Sunapee_long_term_sampling_01Mar2021.R               *
#* AUTHOR:  Amanda Lindsey, Bethel Steele                        *
#* PURPOSE: Collate long term records of monitoring data         *
#*          collected in lake Sunapee                            *
#* LAST UPDATE: v01Mar2021 - update through 2020                 *
#*          v27Jul2020 - update to use tidyverse, update         *
#*          with data through 2019                               *
#*          no version change 25Sept2020: updated bio and do     *
#*          sections                                             *  
#*                                                               *
#*****************************************************************

# set working directory - change if running from a different location
setwd('C:/Users/steeleb/Dropbox/Lake Sunapee/long term Sunapee data')

library(tidyverse)
library(readxl)
library(ggthemes)


# import sunapee id #
sun_id <- read.csv('raw data/sunapee station ids w SUNSUNGEN.csv', header=T)

#************************************************
#### chem data ####

# import data
raw_chem_2017 <- read_xlsx("raw data/LMP files/LMP 2017/CHEMISTRY.xlsx", 
                      sheet="CHEMISTRY",
                      col_types = 'text') %>% 
  mutate(DATE = as.Date(as.numeric(DATE), origin = '1899-12-30'))
raw_chem_2018 <- read_csv('raw data/LMP files/2018/Sunapee 2018 Chem.csv',
                          col_types = cols(.default = col_character())) %>% 
  mutate(DATE = as.Date(DATE, format = '%d-%b-%y'))
raw_chem_2019 <- read_csv('raw data/LMP files/2019/Sunapee2019Chem.csv',
                          col_types = cols(.default = col_character()))%>% 
  mutate(DATE = as.Date(DATE, format = '%d-%b-%y'))
raw_chem_2020 <- read_csv('raw data/LMP files/2020/Sunapee2020CHEM.csv',
                          col_types = cols(.default = col_character()))%>% 
  mutate(DATE = as.Date(DATE, format = '%d-%b-%y'))

str(raw_chem_2017)
str(raw_chem_2018)
str(raw_chem_2019)
str(raw_chem_2020)

raw_chem <- full_join(raw_chem_2017, raw_chem_2018) %>% 
  full_join(., raw_chem_2019) %>% 
  full_join(., raw_chem_2020)


#format chem columns to numeric
raw_chem <- raw_chem %>% 
  mutate_at(vars(PH, H_ION, ALK, COLOR, TP, COND, TURBIDITY, YEAR, MONTH),
            ~ as.numeric(.))

# remove unneeded variables and rename others
raw_chem <- raw_chem %>% 
  select(-"LAKE", -"TOWN", -"YEAR", -"MONTH", -"H_ION", -"COLOR") %>% 
  rename(date = DATE,
         station =STATION,
         depth_m = Depth,
         layer = LAYER, 
         pH = PH,
         alk_mglCaCO3 = ALK,
         TP_mgl = TP,
         cond_uScm = COND, 
         turb_NTU = TURBIDITY)
str(raw_chem)


# set missing data to NA
range(raw_chem$pH, na.rm = T)
raw_chem$pH [raw_chem$pH==-99] = NA
raw_chem$pH [raw_chem$pH>14] = NA
range(raw_chem$pH, na.rm = T)

range(raw_chem$alk_mglCaCO3, na.rm = T)
raw_chem$alk_mglCaCO3 [raw_chem$alk_mglCaCO3<=-99] = NA
raw_chem$alk_mglCaCO3 [raw_chem$alk_mglCaCO3==99] = NA
range(raw_chem$alk_mglCaCO3, na.rm = T)

range(raw_chem$TP_mgl, na.rm = T)
raw_chem$TP_mgl [raw_chem$TP_mgl<=-99] = NA
#add flag for presumed BDL
raw_chem <- raw_chem %>% 
  mutate(TP_flag = case_when(TP_mgl<0.005 ~ 'BDL',
                             TP_mgl < 0 ~ 'negative value reported, recoded to 0',
                             TRUE ~ '')) %>% 
  mutate(TP_mgl = case_when(TP_mgl< 0 ~ 0,
                            TRUE ~ TP_mgl))
range(raw_chem$TP_mgl, na.rm = T)

range(raw_chem$cond_uScm, na.rm = T)
raw_chem$cond_uScm [raw_chem$cond_uScm<=-99] = NA
range(raw_chem$cond_uScm, na.rm = T)

range(raw_chem$turb_NTU, na.rm = T)
raw_chem$turb_NTU [raw_chem$turb_NTU==-99] = NA
range(raw_chem$turb_NTU, na.rm = T)

#plot to check for funky values
#ph
plot(raw_chem$date, raw_chem$pH)
#flag the data above 10
raw_chem <- raw_chem %>% 
  mutate(ph_flag = case_when(pH>10 ~ 'suspect',
                             TRUE ~ ''))

#alk
plot(raw_chem$date, raw_chem$alk_mglCaCO3)
#remove anomolous point
#no anomolous points

#TP
plot(raw_chem$date, raw_chem$TP_mgl)
#remove anomolous point
#no anomolous points

#cond
plot(raw_chem$date, raw_chem$cond_uScm)
#remove anomolous points above 2500
raw_chem$cond_uScm [raw_chem$cond_uScm>2500] = NA
plot(raw_chem$date, raw_chem$cond_uScm)


#turbidity
plot(raw_chem$date, raw_chem$turb_NTU)
#remove anomolous point
#no anamolous points

#COLOR IS INCOMPLETE
# #COLOR
# plot(raw_chem$date, raw_chem$COLOR)
# #remove anomolous point
# #no anamolous points



# calculate H+ ion from pH (H+ = 10^-pH)
raw_chem$H_M=10^(raw_chem$pH * -1)
raw_chem$pH = NULL

#create vertical dataset
raw_chem_vert <- raw_chem %>% 
  select(-c(Comments:Date_Sta_Lr)) %>% 
  gather(parameter, value, -station, -layer, -depth_m, -date, -ph_flag, -TP_flag) %>% 
  filter(!is.na(value)) %>% 
  mutate(flag = case_when(parameter == 'TP_mgl' & !is.na(value) & TP_flag != '' ~ TP_flag,
                          parameter == 'H_M' & !is.na(value) & ph_flag != '' ~ ph_flag,
                          TRUE ~ '')) %>% 
  select(-ph_flag, -TP_flag) %>% 
  mutate(value = as.numeric(value))


ggplot(raw_chem_vert, aes(x = date, y = value, color = flag)) +
  geom_point() +
  facet_grid(parameter~., scales = 'free_y') +
  theme_bw()

stream_chem_vert <- raw_chem_vert %>% 
  mutate(station = as.numeric(station)) %>% 
  filter(station > 230)
unique(stream_chem_vert$station)
ggplot(stream_chem_vert, aes(x = date, y = value, color = flag)) +
  geom_point() +
  facet_grid(parameter~., scales = 'free_y') +
  theme_bw()

lake_chem_vert <- raw_chem_vert %>% 
  mutate(station = as.numeric(station)) %>% 
  filter(station <= 230)
unique(lake_chem_vert$station)
ggplot(raw_chem_vert, aes(x = date, y = value, color = flag)) +
  geom_point() +
  facet_grid(parameter~., scales = 'free_y') +
  theme_bw()

# export raw_chem

write_csv(stream_chem_vert, 'master files/stream_chem_1986-2020_v01Mar2021.csv')
write_csv(lake_chem_vert, 'master files/lake_chem_1986-2020_v01Mar2021.csv')
write_csv(raw_chem, "master files/raw_chem_all_1986-2020_v01Mar2021.csv")



#************************************************


#### bio data ####

# import data
raw_bio_2017 <- read_xlsx("raw data/LMP files/LMP 2017/BIOLOGY.xlsx", 
                           sheet="BIOLOGY",
                           col_types = 'text',
                          na = '-99') %>% 
  mutate(DATE = as.Date(as.numeric(DATE), origin = '1899-12-30'))
raw_bio_2018 <- read_csv('raw data/LMP files/2018/Sunapee 2018 Bio.csv',
                          col_types = cols(.default = col_character()),
                         na = '-99') %>% 
  mutate(DATE = as.Date(DATE, format = '%m/%d/%Y'))
raw_bio_2019 <- read_csv('raw data/LMP files/2019/Sunapee2019Bio.csv',
                         col_types = cols(.default = col_character()),
                         na = '-99')%>% 
  mutate(DATE = as.Date(DATE, format = '%d-%b-%y'))
raw_bio_2020 <- read_csv('raw data/LMP files/2020/Sunapee2020BIO.csv',
                         col_types = cols(.default = col_character()),
                         na = '-99')%>% 
  mutate(DATE = as.Date(DATE, format = '%d-%b-%y'))
str(raw_bio_2017)
str(raw_bio_2018)
str(raw_bio_2019)
str(raw_bio_2020)

raw_bio <- full_join(raw_bio_2017, raw_bio_2018) %>% 
  full_join(., raw_bio_2019) %>% 
  full_join(., raw_bio_2020)
head(raw_bio)

#format bio columns to numeric
raw_bio <- raw_bio %>% 
  mutate_at(vars(CHL, SD, PCTPHY1, PCTPHY2, PCTPHY3),
            ~ as.numeric(.))

# remove unneeded variables and rename others
raw_bio <- raw_bio %>% 
  select(-"LAKE", -"TOWN", -"YEAR", -"MONTH", -'CompleteDate', -'SampleID', -'ID', - 'Date_Sta') %>% 
  rename(date = DATE,
         station =STATION,
         chla_ugl = CHL,
         secchidepth_m = SD)
str(raw_bio)

# format, QAQC and standardize variables

#chlorophyll-a
range(raw_bio$chla_ugl, na.rm = T)
raw_bio$chla_ugl [raw_bio$chla_ugl==-99.99] = NA
raw_bio$chla_ugl [raw_bio$chla_ugl==-9.99] = NA
range(raw_bio$chla_ugl, na.rm = T)
#chl-a
plot(raw_bio$date, raw_bio$chla_ugl)
#no anomolous points
raw_bio$flag_chla <- NA_character_
raw_bio$flag_chla [raw_bio$chla_ugl==0] = 'suspect'

#secchi depth
range(raw_bio$secchidepth_m, na.rm = T) 
#secchi of 0 seems suspect
plot(raw_bio$date, raw_bio$secchidepth_m)

#remove anomolous point
raw_bio$secchidepth_m [raw_bio$station==30 & raw_bio$date=="1999-07-15" & raw_bio$secchidepth_m>25] = NA #point deeper than maximum depth
plot(raw_bio$date, raw_bio$secchidepth_m)

#flag secchi of 0
raw_bio$flag_secchi <- NA_character_
raw_bio$flag_secchi [raw_bio$secchidepth_m==0] = 'suspect'


# drop phyto data (there are only about 20 observations of phyto data) and create vertical dataset of clha and sd
raw_bio_vert <- raw_bio %>% 
  select(-(NETPHY1:PCTPHY3)) %>% 
  gather(parameter, value, -station, -date, -flag_secchi, -flag_chla) %>% 
  mutate(flag = case_when(parameter == 'chla_ugl' & !is.na(flag_chla) ~ flag_chla,
                          parameter == 'secchidepth_m' & !is.na(flag_secchi) ~ flag_secchi,
                          TRUE ~ NA_character_)) %>% 
  select(-flag_secchi, -flag_chla)
head(raw_bio_vert)

ggplot(raw_bio_vert, aes(x=date, y=value, color = flag)) +
  geom_point() +
  facet_grid(parameter~station) +
  theme_bw()

# filter for deep spots
raw_bio_deep <- raw_bio_vert %>% 
  filter(station == 200 | station == 210 | station == 220 | station == 230)

#plot historical deep spot records
ggplot(raw_bio_deep, aes(x=date, y=value, color = flag)) +
  geom_point() +
  facet_grid(parameter~station) +
  theme_bw()


# export raw_bio
write.csv(raw_bio_vert, file="master files/lake_bio_1986-2020_v01Mar2021.csv", row.names=F, na="")



#### DO data ####
raw_do_2017 <- read_xlsx("raw data/LMP files/LMP 2017/DO_BGSqaqc.xlsx", 
                          sheet="DO",
                          col_types = 'text') %>% 
  mutate(DATE = as.Date(as.numeric(DATE), origin = '1899-12-30'))
raw_do_2018 <- read_csv('raw data/LMP files/2018/Sunapee 2018 DO.csv',
                         col_types = cols(.default = col_character())) %>% 
  mutate(DATE = as.Date(DATE, format = '%d-%b-%y'))
raw_do_2019 <- read_csv('raw data/LMP files/2019/Sunapee2019DO.csv',
                         col_types = cols(.default = col_character()),
                         col_names = c('DATE', 'TIME', 'STATION', 'SENSOR', 'TEMP', 'PCNTSAT','DO',  'DEPTH'),
                        skip = 1)%>% 
  mutate(DATE = as.Date(DATE, format = '%m/%d/%Y'))
raw_do_2020 <- read_csv('raw data/LMP files/2020/Sunapee2020DO.csv',
                        col_types = cols(.default = col_character()),
                        skip = 1,
                        col_names = c('LAKE', 'TOWN','STATION', 'DATE', 'DEPTH', 'TEMP', 'DO', 'PCNTSAT', 'SPC_USCM', 'PH', 'NTU', 'TIME', 'BOTTOMZ', 'WEATHER', 'COMMENTS', 'ID', 'DATESTAID'))%>% 
  mutate(DATE = as.Date(DATE, format = '%d-%b-%y'))
str(raw_do_2017)
str(raw_do_2018)
str(raw_do_2019)
str(raw_do_2020)

raw_do <- full_join(raw_do_2017, raw_do_2018) %>% 
  full_join(., raw_do_2019) %>% 
  full_join(., raw_do_2020)
head(raw_do)
unique(raw_do$Comments)
unique(raw_do$WEATHER)
unique(raw_do$BOTTOMZ)

#select pertinent columns
raw_do <- raw_do %>% 
  select(STATION, DATE, DEPTH, TEMP, DO, PCNTSAT, SPC_USCM, PH, NTU, BOTTOMZ)

#format depth, temp, do columns to numeric
raw_do <- raw_do %>% 
  mutate_at(vars(DEPTH, TEMP, DO, PCNTSAT, SPC_USCM, PH, NTU, BOTTOMZ),
            ~ as.numeric(.))

#quick reality check for anamolous points
range(raw_do$TEMP, na.rm = T)
plot(raw_do$DATE, raw_do$TEMP)

#values where temp less than 4.5 in june and august seem errant, recode
ix = which(raw_do$TEMP<4.5 & as.numeric(format(raw_do$DATE, '%m')) >=4 & as.numeric(format(raw_do$DATE, '%m')) <10)

raw_do$TEMP[ix] = NA_real_
plot(raw_do$DATE, raw_do$TEMP)




range(raw_do$DO, na.rm = T)

#recode do and pctsat where raw do is 0 or less 
raw_do <- raw_do %>% 
  mutate_at(vars('DO', 'PCNTSAT'),
            ~ case_when(DO <=0 ~ NA_real_,
                        TRUE ~ .))
range(raw_do$DO, na.rm = T)

plot(raw_do$DATE, raw_do$DO)

#recode those greater than 80 as NA
raw_do$DO [raw_do$DO>80] = NA_real_
plot(raw_do$DATE, raw_do$DO)


range(raw_do$PCNTSAT, na.rm = T)
plot(raw_do$DATE, raw_do$PCNTSAT)
#recode those greater than 200 as NA
raw_do$PCNTSAT [raw_do$PCNTSAT>200] = NA_real_
plot(raw_do$DATE, raw_do$PCNTSAT)
#flag those GT 120 as suspect
raw_do$do_flag = NA_character_
raw_do$do_flag [raw_do$PCNTSAT>140] = 'DO PCT >140'

ggplot(raw_do, aes(x = DATE, y = DO, color = do_flag)) +
  geom_point()

#plot ph
plot(raw_do$DATE, raw_do$PH)
range(raw_do$PH, na.rm = T)

#plot turb
plot(raw_do$DATE, raw_do$NTU)
range(raw_do$NTU, na.rm = T)
#recode negative to NA
raw_do$NTU [raw_do$NTU < 0] = NA_real_
#LOOKS LIKE THERE ARE SOME PROBLEMATIC DATA HERE, LIKELY SONDE IS IN SEDIMENT where NTU > 150
raw_do <- raw_do %>% 
  mutate_at(vars(TEMP, DO, PCNTSAT, SPC_USCM, NTU, PH),
            ~ case_when(NTU > 150 ~ NA_real_,
                        TRUE ~ .)) %>% 
  mutate(flag = case_when(NTU > 20 ~ 'sensor may be in sediment, turbidity > 20 NTU',
                          TRUE ~ NA_character_))
plot(raw_do$DATE, raw_do$NTU)


# plot cond
plot(raw_do$DATE, raw_do$SPC_USCM)
range(raw_do$SPC_USCM, na.rm = T)


# rename columns
raw_do <- raw_do %>% 
  rename(station = STATION,
         depth_m = DEPTH,
         temp_C = TEMP,
         DO_mgl = DO,
         DO_pctsat = PCNTSAT,
         turb_NTU = NTU,
         pH = PH,
         cond_uscm = SPC_USCM,
         date = DATE) %>% 
  select(-BOTTOMZ)
str(raw_do)

#create vertical dataset
raw_do_vert <- raw_do %>% 
  gather(parameter, value, -station, -depth_m, -date, -do_flag, -flag)

ggplot(raw_do_vert, aes(x = date, y = value, color = flag)) +
  geom_point() +
  facet_grid(parameter ~ station, scales = 'free_y') +
  theme_bw()

ggplot(raw_do_vert, aes(x = date, y = value, color = do_flag)) +
  geom_point() +
  facet_grid(parameter ~ station, scales = 'free_y') +
  theme_bw()

# export raw_do
write.csv(raw_do_vert, file="master files/raw_do_1986-2020_v01Mar2021.csv.csv", row.names=F, na="")

