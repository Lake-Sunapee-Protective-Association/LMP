#*****************************************************************
#*      Cary Institute of Ecosystem Studies (Millbrook, NY)      *
#*                on behalf of                                   *
#*              the Lake Sunapee Protective Assosciation         *
#*                                                               *
#* TITLE:   LSPA_LMP_collation.R                                 *
#* AUTHOR:  Bethel Steele steeleb@caryinstitute.org              *
#* RStudio version: 2022.12                                      *
#* R version: 4.2.2                                              *
#* PURPOSE: Collate long term records of monitoring data         *
#*          collected in the Lake Sunapee watershed by the LSPA  *
#* LAST UPDATE: v25May2023 - fix site type labels and include    *
#*                area lakes in collation                        *
#*              v22Dec2022 - update data through 2022            *
#*              v10May2022 - update through 2021 and add some    *
#*                flags to previous archived data                *
#*              v08Mar2021 - update in Git                       *
#*              v01Mar2021 - update through 2020                 *
#*              v27Jul2020 - update to use tidyverse, update     *
#*          with data through 2019                               *
#*          no version change 25Sept2020: updated bio and do     *
#*          sections                                             *  
#*****************************************************************

# point to data directories
datadir = 'raw data files/'
dumpdir = 'primary files/'
areadump = 'raw data files/area lakes subset/'
figuredump = 'collation code/tempfigs/'

library(tidyverse) #v1.3.2
library(readxl)
library(ggthemes)


# The data from the LSPA come from three different files: Chemistry, Biology, and DO. In the recent past, 
#   there have been fluctuations of where data are stored, so this script collates the data, removes obviously
#   errant values, flags suspect data, and collates all of the data together. At that collation step, further
#   QAQC is completed on the data to recode datat recorded that were obtained with the sensor in the sediment.

# CHEMISTRY DATA ####

# import data
raw_chem_2017 <- read_xlsx(paste0(datadir, "1986-2017/CHEMISTRY.xlsx"), 
                      sheet="CHEMISTRY",
                      col_types = 'text') %>% 
  mutate(DATE = as.Date(as.numeric(DATE), origin = '1899-12-30'))
raw_chem_2018 <- read_csv(paste0(datadir, "2018/Sunapee 2018 Chem.csv"),
                          col_types = cols(.default = col_character())) %>% 
  mutate(DATE = as.Date(DATE, format = '%d-%b-%y'))
raw_chem_2019 <- read_csv(paste0(datadir,'2019/Sunapee2019Chem.csv'),
                          col_types = cols(.default = col_character()))%>% 
  mutate(DATE = as.Date(DATE, format = '%d-%b-%y'))
raw_chem_2020 <- read_csv(paste0(datadir,'2020/Sunapee2020CHEM.csv'),
                          col_types = cols(.default = col_character()))%>% 
  mutate(DATE = as.Date(DATE, format = '%d-%b-%y'))
raw_chem_2021 <- read_xlsx(paste0(datadir,'2021-2022/LMPCHEM_2021_2022.xlsx'),
                          col_types = 'text') %>% 
  mutate(DATE = as.Date(as.numeric(DATE), origin = '1899-12-30'))


### collate chem data ####

comp_chem <- full_join(raw_chem_2017, raw_chem_2018) %>% 
  full_join(., raw_chem_2019) %>% 
  full_join(., raw_chem_2020) %>% 
  full_join(., raw_chem_2021) %>% 
  arrange(DATE)
colnames(comp_chem)

chem_vars = c("Depth","PH","H_ION","COLOR", "ALK","TP","COND","TURBIDITY","Chloride")

#format chem columns to numeric
comp_chem <- comp_chem %>%
  mutate_at(vars(all_of(chem_vars), STATION),
            ~ as.numeric(.))

# pull out sunapee from other data without station ids
unique(comp_chem$LAKE)
unique(comp_chem$STATION)

# grab area lakes
area_chem <- comp_chem %>% 
  filter(!grepl('sunapee', LAKE, ignore.case =T))
write.csv(area_chem, file.path(areadump, 'area_chem_subset.csv'), row.names = F)

# remove unneeded variables and rename others
raw_chem <- comp_chem %>%
  filter(grepl('sunapee', LAKE, ignore.case = T)) %>% #drom non-LSPA LMP sites
  select(-"LAKE", -"TOWN", -"YEAR", -"MONTH",-'H_ION', -'Date_Sta_Lr',-'CompleteDate', -'StationText') %>% #no data in COLOR or H_ION variable
  rename(date = DATE,
         station =STATION,
         depth_m = Depth,
         layer = LAYER, 
         pH = PH,
         alk_mglCaCO3 = ALK,
         TP_mgl = TP,
         cond_uScm = COND, 
         turb_NTU = TURBIDITY,
         cl_mgl = Chloride,
         color_PCU = COLOR,
         org_id = ID,
         org_sampid = 'SampleID') %>% 
  mutate(station = case_when(station < 0 ~ NA_real_,
                             TRUE ~ station)) %>% 
  mutate(site_type = case_when(station <=230 ~ 'lake',
                              station >230 ~ 'tributary',
                              is.na(station) ~ 'new site',
                              TRUE ~ '')) %>% 
  relocate(station, site_type, date)

#idiot check for station labels
stations = raw_chem %>% 
  arrange(station) %>% 
  mutate (station_type = paste(station, site_type)) %>% 
  distinct(station_type)
view(stations)

### baseline chem data clean up ####
# start with a new dataframe
qaqc_chem <- raw_chem

# look at layer info:
unique(qaqc_chem$layer)
qaqc_chem$layer [qaqc_chem$layer=="0"] = "O" #outlet
qaqc_chem$layer [qaqc_chem$layer=="1"] = "I" #presumed transcription error this is a shallow location
qaqc_chem$layer [qaqc_chem$layer=="L"] = "I" #presumed transcription error all others are integrated
qaqc_chem$layer [qaqc_chem$layer=="N"] = "" #this is the layer for many of the early records, recoding to blank
qaqc_chem$layer [qaqc_chem$layer=="5"] = "" #this is in a tributary
unique(qaqc_chem$layer)

# all stations < 110 are I
qaqc_chem$layer [qaqc_chem$station <110 & !is.na(qaqc_chem$station) & qaqc_chem$layer == ''] = 'I'

#recode depth and layer for tributary samples
qaqc_chem$depth_m [qaqc_chem$site_type=="tributary"] = NA_real_
qaqc_chem$layer [qaqc_chem$site_type=="tributary"] = ""

unique(qaqc_chem$layer)
unique(qaqc_chem$depth_m)

### look at ranges and recode data that are obviously NA or errant ####
colnames(qaqc_chem)

# set missing data to NA
range(qaqc_chem$pH, na.rm = T)
qaqc_chem$pH [qaqc_chem$pH==-99] = NA
range(qaqc_chem$pH, na.rm = T)
qaqc_chem$pH [qaqc_chem$pH>14] = NA
range(qaqc_chem$pH, na.rm = T)

range(qaqc_chem$alk_mglCaCO3, na.rm = T)
qaqc_chem$alk_mglCaCO3 [qaqc_chem$alk_mglCaCO3<=-99] = NA
qaqc_chem$alk_mglCaCO3 [qaqc_chem$alk_mglCaCO3==99] = NA
range(qaqc_chem$alk_mglCaCO3, na.rm = T)

range(qaqc_chem$color_PCU, na.rm = T)
qaqc_chem$color_PCU [qaqc_chem$color_PCU<=0] = NA
range(qaqc_chem$color_PCU, na.rm = T)

range(qaqc_chem$TP_mgl, na.rm = T)
qaqc_chem$TP_mgl [qaqc_chem$TP_mgl<=-99] = NA
range(qaqc_chem$TP_mgl, na.rm = T) #still some negative values, will recode to BDL

#add flag for presumed BDL
qaqc_chem <- qaqc_chem %>% 
  mutate(TP_flag = case_when(TP_mgl<0.005 ~ 'BDL',
                             TP_mgl < 0 ~ 'BDL; negative value reported: recoded to 0',
                             TRUE ~ '')) %>% 
  mutate(TP_mgl = case_when(TP_mgl< 0 ~ 0,
                            TRUE ~ TP_mgl))
range(qaqc_chem$TP_mgl, na.rm = T)

range(qaqc_chem$cond_uScm, na.rm = T)
qaqc_chem$cond_uScm [qaqc_chem$cond_uScm<=-99] = NA
range(qaqc_chem$cond_uScm, na.rm = T)
#value for cond >9000 was in middle of lake - recoding, because that is obviously errant. The other values >1000 are in tributarys, which seems reasonable.
qaqc_chem$cond_uScm [qaqc_chem$cond_uScm>9000] = NA
range(qaqc_chem$cond_uScm, na.rm = T)

range(qaqc_chem$turb_NTU, na.rm = T)
qaqc_chem$turb_NTU [qaqc_chem$turb_NTU==-99] = NA
range(qaqc_chem$turb_NTU, na.rm = T)

range(qaqc_chem$cl_mgl, na.rm = T)

### plot to check for funky values ####
#### ph ####
ggplot(qaqc_chem, aes(x = date, y = pH)) +
  facet_grid(site_type~. ) +
  geom_point()

#recode the data in lake above 10
qaqc_chem$pH [qaqc_chem$pH>10] = NA
qaqc_chem$pH_flag = ''
qaqc_chem$pH_flag [qaqc_chem$pH<5.5 & qaqc_chem$site_type == 'lake'] = 'low pH suspect'

#replot
ggplot(qaqc_chem, aes(x = date, y = pH, color = pH_flag)) +
  facet_grid(site_type~. ) +
  geom_point(aes(shape = layer))

ggplot(qaqc_chem, aes(x = as.numeric(format(date, '%j')), y = pH, color = pH_flag)) +
  geom_point() +  
  facet_grid(site_type~. )

qaqc_chem$pH_flag [qaqc_chem$pH>7.5 & qaqc_chem$site_type == 'tributary' & as.numeric(format(qaqc_chem$date, '%j'))>275] = 'high pH suspect'

ggplot(qaqc_chem, aes(x = as.numeric(format(date, '%j')), y = pH, color = pH_flag)) +
  geom_point() +  
  facet_grid(site_type~. )

# calculate H+ ion from pH (H+ = 10^-pH) and drop pH column
qaqc_chem$conc_H_molpl=10^(qaqc_chem$pH * -1)
qaqc_chem$pH = NULL


#### alk ####
ggplot(qaqc_chem, aes(x = date, y = alk_mglCaCO3, color = depth_m)) +
  facet_grid(site_type~. ) +
  geom_point(aes(shape = layer))

#flag 0 values as well as in-lake values >9
qaqc_chem <- qaqc_chem %>% 
  mutate(year = format(date, '%Y')) %>% 
  mutate(alk_flag = case_when(alk_mglCaCO3 == 0 ~ 'low alk suspect',
                              alk_mglCaCO3 > 9 & site_type == 'lake' ~ 'high alk suspect',
                              (year == 1995 | year == 1993 |  year == 1998 | year == 2000) & alk_mglCaCO3 <2.5 & site_type == 'lake' ~ 'low alk suspect',
                              year == 1993 & site_type == 'tributary' & alk_mglCaCO3 < 1.5 ~ 'low alk suspect',
                              TRUE ~ ''))

ggplot(qaqc_chem, aes(x = date, y = alk_mglCaCO3, color = alk_flag)) +
  facet_grid(site_type~. ) +
  geom_point(aes(shape = layer)) +
  scale_x_date(minor_breaks = '1 year')

#### TP ####
ggplot(qaqc_chem, aes(x = date, y = TP_mgl, color = depth_m)) +
  facet_grid(site_type~. , scales = 'free_y') +
  geom_point(aes(shape = layer))

#### cond ####
ggplot(qaqc_chem, aes(x = date, y = cond_uScm, color = depth_m)) +
  facet_grid(site_type~. , scales = 'free_y') +
  geom_point(aes(shape = layer)) +
  scale_x_date(minor_breaks = '1 year')

#remove anomolous point above 2500
qaqc_chem <- qaqc_chem %>% 
  mutate(cond_flag = case_when(cond_uScm>2500 ~ 'high cond suspect',
                                site_type == 'lake' & cond_uScm < 47 ~ 'low cond suspect',
                               site_type == 'lake' & as.numeric(format(date, '%Y')) == 1997 &
                                 cond_uScm <65 ~ 'low cond suspect',
                                TRUE ~ ''))

ggplot(qaqc_chem, aes(x = date, y = cond_uScm, color = cond_flag)) +
  facet_grid(site_type~., scales = 'free_y') +
  geom_point(aes(shape = layer)) +
  scale_x_date(minor_breaks = '1 year')

ggplot(qaqc_chem, aes(x = date, y = cond_uScm, color = cond_flag)) +
  facet_grid(site_type~.) +
  scale_y_continuous(limit = c(0, 50)) +
  geom_point(aes(shape = layer)) +
  scale_x_date(minor_breaks = '1 year')

#flag low tributary values
qaqc_chem <- qaqc_chem %>% 
  mutate(cond_flag = case_when(cond_uScm<9 & site_type == 'tributary' ~ 'low cond suspect',
                               TRUE ~ cond_flag))

ggplot(qaqc_chem, aes(x = date, y = cond_uScm, color = cond_flag)) +
  facet_grid(site_type~. , scales = 'free_y') +
  geom_point(aes(shape = layer))+
  scale_x_date(minor_breaks = '1 year')


#### turbidity ####
ggplot(qaqc_chem, aes(x = date, y = turb_NTU, color = depth_m)) +
  facet_grid(site_type~. , scales = 'free_y') +
  geom_point(aes(shape = layer))

#flag in-lake above 30 and tributary > 500 
qaqc_chem <- qaqc_chem %>% 
  mutate(turb_flag = case_when(turb_NTU >30 & site_type == 'lake' ~ 'high turb suspect unless recent storm', 
                               turb_NTU >500 & site_type == 'tributary' ~ 'suspect unless recent storm', 
                               TRUE ~ ''))

ggplot(qaqc_chem, aes(x = date, y = turb_NTU, color = turb_flag)) +
  facet_grid(site_type~. , scales = 'free_y') +
  geom_point(aes(shape = layer))

#zoom in on values near zero
ggplot(qaqc_chem, aes(x = date, y = turb_NTU, color = turb_flag)) +
  facet_grid(site_type~.) +
  scale_y_continuous(limits = c(0, 5))+
  geom_point(aes(shape = layer))

#### chloride ####
ggplot(qaqc_chem, aes(x = date, y = cl_mgl, color = depth_m)) +
  facet_grid(site_type~. , scales = 'free_y') +
  geom_point(aes(shape = layer))

#### color ####
ggplot(qaqc_chem, aes(x = date, y = color_PCU)) +
  facet_grid(site_type~. , scales = 'free_y') +
  geom_point(aes(shape = layer))

#need to recode values of 99 prior to 2021
qaqc_chem$color_PCU = ifelse(qaqc_chem$color_PCU == 99 & qaqc_chem$date < as.Date('2021-01-01'), NA_real_, qaqc_chem$color_PCU)

ggplot(qaqc_chem, aes(x = date, y = color_PCU)) +
  facet_grid(site_type~. , scales = 'free_y') +
  geom_point(aes(shape = layer))


#### look at comments ####
unique(qaqc_chem$Comments)
# move the storm event into own column
qaqc_chem <- qaqc_chem %>% 
  mutate(gen_flag = case_when(Comments == 'STORM EVENT' ~ 'storm event sampling',
                              TRUE ~ '')) %>% 
  select(-Comments)

### create vertical dataset, eliminate na rows and export master files ####
#create vertical dataset
colnames(qaqc_chem)
qaqc_chem_vert <- qaqc_chem %>% 
  select(-year) %>% 
  pivot_longer(names_to = 'parameter', values_to = 'value',
               cols = c("alk_mglCaCO3", "color_PCU", "TP_mgl", "cond_uScm", "turb_NTU", "cl_mgl", "conc_H_molpl"))

qaqc_chem_vert <- qaqc_chem_vert %>% 
  filter(!is.na(value)) %>% 
  mutate(flag = case_when(parameter == 'TP_mgl' & TP_flag != '' ~ TP_flag,
                          parameter == 'alk_mglCaCO3' & alk_flag != '' ~ alk_flag,
                          parameter == 'cond_uScm' & cond_flag != '' ~ cond_flag,
                          parameter == 'turb_NTU' & turb_flag != '' ~ turb_flag,
                          parameter == 'conc_H_molpl' & pH_flag != '' ~ pH_flag,
                          TRUE ~ '')) %>% 
  mutate(flag = case_when(gen_flag != '' & flag != '' ~ paste(flag, gen_flag, sep = '; '),
                          gen_flag != '' & flag == '' ~ gen_flag,
                          TRUE ~ flag)) %>% 
  select(-TP_flag,-alk_flag, -cond_flag, -turb_flag, -gen_flag, -pH_flag) 

#plot all data with flags
ggplot(qaqc_chem_vert, aes(x = date, y = value, color = flag)) +
  geom_point() +
  facet_grid(parameter~site_type, scales = 'free_y') +
  theme_bw()

#plot only tributary data
tributary_chem_vert <- qaqc_chem_vert %>% 
  filter(site_type == 'tributary') %>% 
  select(-depth_m, -layer)
unique(tributary_chem_vert$station)

ggplot(tributary_chem_vert, aes(x = date, y = value, color = flag)) +
  geom_point() +
  facet_grid(parameter~., scales = 'free_y') +
  theme_bw()

#plot only lake data
lake_chem_vert <- qaqc_chem_vert %>% 
  filter(site_type == 'lake')
unique(lake_chem_vert$station)
ggplot(qaqc_chem_vert, aes(x = date, y = value, color = flag)) +
  geom_point() +
  facet_grid(parameter~., scales = 'free_y') +
  theme_bw()

#apply ODM cv
unique(qaqc_chem_vert$parameter)
qaqc_chem_vert <- qaqc_chem_vert %>% 
  mutate(parameter = case_when(parameter == 'alk_mglCaCO3' ~ 'alkalinity_mglCaCO3',
                               parameter == 'TP_mgl' ~ 'phosphorusTotal_mgl',
                               parameter == 'cond_uScm' ~ 'specificConductance_uScm',
                               parameter == 'turb_NTU' ~ 'turbidity_NTU',
                               parameter == 'cl_mgl' ~ 'chloride_mgl',
                               parameter == 'conc_H_molpl' ~ 'hydrogenDissolved_moll',
                               TRUE ~ parameter))


# BIOLOGICAL DATA ####

#NOTE, IN 2021 DATA STORAGE CHANGES - BOTH BIO AND DO, AS WELL AS SOME OTHER PARAMS ARE STORED IN ONE FILE

# import data
raw_bio_2017 <- read_xlsx(paste0(datadir, "1986-2017/BIOLOGY.xlsx"), 
                           sheet="BIOLOGY",
                           col_types = 'text') %>% 
  mutate(DATE = as.Date(as.numeric(DATE), origin = '1899-12-30'))

raw_bio_2018 <- read_csv(paste0(datadir, "2018/Sunapee 2018 Bio.csv"),
                          col_types = cols(.default = col_character()),
                         na = '-99') %>% 
  mutate(DATE = as.Date(DATE, format = '%m/%d/%Y'))

raw_bio_2019 <- read_csv(paste0(datadir,'2019/Sunapee2019Bio.csv'),
                         col_types = cols(.default = col_character()),
                         na = '-99')%>% 
  mutate(DATE = as.Date(DATE, format = '%d-%b-%y'))

raw_bio_2020 <- read_csv(paste0(datadir,'2020/Sunapee2020BIO.csv'),
                         col_types = cols(.default = col_character()),
                         na = '-99')%>% 
  mutate(DATE = as.Date(DATE, format = '%d-%b-%y'))

raw_bio_2021 <- read_xlsx(file.path(datadir, '2021-2022/LMPBIO_2021_2022.xlsx'),
                          col_types = 'text') %>% 
  mutate(DATE = as.Date(as.numeric(DATE), origin = '1899-12-30'))

#check for character strings
str(raw_bio_2017)
str(raw_bio_2018)
str(raw_bio_2019)
str(raw_bio_2020)
str(raw_bio_2021)

### collate bio data ####
comp_bio <- full_join(raw_bio_2017, raw_bio_2018) %>% 
  full_join(., raw_bio_2019) %>% 
  full_join(., raw_bio_2020) %>% 
  full_join(., raw_bio_2021)
head(comp_bio)
colnames(comp_bio)


# format bio columns to numeric
comp_bio <- comp_bio %>% 
  mutate_at(vars(STATION, CHL, SD, SD_Scope, PCTPHY1, PCTPHY2, PCTPHY3),
            ~ as.numeric(.))

#grab area lakes
area_bio = comp_bio %>% 
  filter(!grepl('sunapee', LAKE, ignore.case = T))
write.csv(area_bio, file.path(areadump, 'area_bio_subset.csv'), row.names = F)

# remove unneeded variables and rename others
raw_bio <- comp_bio %>% 
  filter(grepl('sunapee', LAKE, ignore.case = T)) %>% #drop non LSPA data
  select(-"LAKE", -"TOWN", -"YEAR", -"MONTH", -'CompleteDate',- 'Date_Sta', -'StationText') %>% 
  rename(date = DATE,
         station =STATION,
         chla_ugl = CHL,
         secchidepth_m = SD,
         org_id = ID,
         org_sampid = SampleID,
         phyto_net_a = NETPHY1,
         phyto_pct_a = PCTPHY1,
         phyto_net_b = NETPHY2,
         phyto_pct_b = PCTPHY2,
         phyto_net_c = NETPHY3,
         phyto_pct_c = PCTPHY3) %>% 
  mutate(site_type = case_when(station <=230 ~ 'lake',
                              station >230 ~ 'tributary')) %>% 
  relocate(station, site_type, date)
head(raw_bio)


# check to see if the 2021/2022 SD and SD scope are the same
raw_bio %>% 
  filter(!is.na(SD_Scope)) %>% 
  filter(secchidepth_m != SD_Scope)
#these are all the same
raw_bio %>%
  filter(!is.na(SD_Scope) & is.na(secchidepth_m))
# there are no instances where scope has data that other doesn't, we can drop this column

raw_bio <- raw_bio %>% select(-SD_Scope)

### baseline cleanup and QAQC ####

#save into new dataframe
qaqc_bio <- raw_bio

# filter out oddball tributary samples (1700)
qaqc_bio <- qaqc_bio %>% 
  filter(site_type == 'lake')

##### chlorophyll-a #####
range(qaqc_bio$chla_ugl, na.rm = T)
qaqc_bio$chla_ugl [qaqc_bio$chla_ugl==-99.99] = NA
range(qaqc_bio$chla_ugl, na.rm = T)
qaqc_bio$chla_ugl [qaqc_bio$chla_ugl==-99] = NA
range(qaqc_bio$chla_ugl, na.rm = T)
qaqc_bio$chla_ugl [qaqc_bio$chla_ugl==-9.99] = NA
range(qaqc_bio$chla_ugl, na.rm = T)

# plot chl-a
ggplot(qaqc_bio, aes(x = date, y = chla_ugl)) +
  geom_point() +
  facet_grid(station ~ .)

#no anomolous points
qaqc_bio$flag_chla <- ''
qaqc_bio$flag_chla [qaqc_bio$chla_ugl<1] = 'likely BDL'
qaqc_bio$flag_chla [qaqc_bio$chla_ugl==0] = 'suspect'


# plot chl-a with flags
ggplot(qaqc_bio, aes(x = date, y = chla_ugl, color = flag_chla)) +
  geom_point() +
  facet_grid(station ~ .)


##### secchi depth ####
range(qaqc_bio$secchidepth_m, na.rm = T) 
qaqc_bio$secchidepth_m [qaqc_bio$secchidepth_m==-99] = NA
range(qaqc_bio$secchidepth_m, na.rm = T)
#will need to flag secchi of 0

# plot secchi
ggplot(qaqc_bio, aes(x = date, y = secchidepth_m)) +
  geom_point() 

#remove anomolous point
qaqc_bio$secchidepth_m [qaqc_bio$station==30 & qaqc_bio$date=="1999-07-15" & qaqc_bio$secchidepth_m>25] = NA #point deeper than maximum depth

#plot again
ggplot(qaqc_bio, aes(x = date, y = secchidepth_m)) +
  geom_point() 


#flag secchi of 0 and other repeated values
qaqc_bio$flag_secchi <- ''
qaqc_bio$flag_secchi [qaqc_bio$secchidepth_m==0] = 'suspect'
qaqc_bio$flag_secchi [qaqc_bio$secchidepth_m==2 & qaqc_bio$station == 80] = 'suspect - possible bottom hit'
qaqc_bio$flag_secchi [qaqc_bio$secchidepth_m==4 & (qaqc_bio$station == 20 | qaqc_bio$station == 60) ] = 'suspect - possible bottom hit'
qaqc_bio$flag_secchi [qaqc_bio$secchidepth_m==4.3 & qaqc_bio$station == 20] = 'suspect - possible bottom hit'

#and by station
ggplot(qaqc_bio, aes(x=date, y=secchidepth_m, color = flag_secchi)) +
  geom_point() +
  facet_grid(.~station) +
  theme_bw()

#recode value at site 20 > 10m - that is almost certainly a typo
qaqc_bio$secchidepth_m [qaqc_bio$secchidepth_m>10 & qaqc_bio$station == 20] = NA


#plot again with flags
ggplot(qaqc_bio, aes(x = date, y = secchidepth_m, color = flag_secchi )) +
  geom_point() 

### phyto data into separate file (very few observations) ####
phyto <- qaqc_bio %>% 
  select(station, date, phyto_net_a:phyto_pct_c) %>% 
  filter(!is.na(phyto_net_a))

phyto_net = c('phyto_net_a', 'phyto_net_b', 'phyto_net_c')
phyto_pct = c('phyto_pct_a', 'phyto_pct_b','phyto_pct_c')

phyto <- phyto %>%
  mutate_at(vars(all_of(phyto_pct)),
            ~ case_when(phyto_pct_a < 0 ~ NA_real_,
                        TRUE ~ .)) %>% 
  filter(!is.na(phyto_pct_a))

# export phytos
firstsamp <- format(as.Date(min(phyto$date)), '%Y')
lastsamp <- format(as.Date(max(phyto$date)), '%Y')
#write_csv(phyto, paste0(dumpdir, 'lake/phyto_', firstsamp, '_', lastsamp, '_v', Sys.Date(),'.csv'))

### create vertical dataset for secchi and chla ####
qaqc_bio_vert <- qaqc_bio %>% 
  select(-(phyto_net_a:phyto_pct_c)) %>% 
  gather(parameter, value, -station, -date, -org_id, -org_sampid, -flag_secchi, -flag_chla, -site_type) %>% 
  mutate(flag = case_when(parameter == 'chla_ugl' & !is.na(flag_chla) ~ flag_chla,
                          parameter == 'secchidepth_m' & !is.na(flag_secchi) ~ flag_secchi,
                          TRUE ~ '')) %>% 
  select(-flag_secchi, -flag_chla)
head(qaqc_bio_vert)

ggplot(qaqc_bio_vert, aes(x=date, y=value, color = flag)) +
  geom_point() +
  facet_grid(parameter~station) +
  theme_bw()

# filter for deep spots
qaqc_bio_deep <- qaqc_bio_vert %>% 
  filter(station == 200 | station == 210 | station == 220 | station == 230)

#plot historical deep spot records
ggplot(qaqc_bio_deep, aes(x=date, y=value, color = flag)) +
  geom_point() +
  facet_grid(parameter~station) +
  theme_bw()

#add flag for hight 210 chl-a (May 2022)
qaqc_bio_vert <- qaqc_bio_vert %>% 
  mutate(flag = case_when(parameter == 'chla_ugl' & station == '210' & value > 10 ~ 'suspect',
                          TRUE ~ flag))

# filter for deep spots
qaqc_bio_deep <- qaqc_bio_vert %>% 
  filter(station == 200 | station == 210 | station == 220 | station == 230)

#plot historical deep spot records
ggplot(qaqc_bio_deep, aes(x=date, y=value, color = flag)) +
  geom_point() +
  facet_grid(parameter~station) +
  theme_bw()

#rename with ODM2 CV
unique(qaqc_bio_vert$parameter)
qaqc_bio_vert <- qaqc_bio_vert %>% 
  mutate(parameter = case_when(parameter == 'chla_ugl' ~ 'chlorophyll_a_ugl',
                               parameter == 'secchidepth_m' ~ 'secchiDepth_m',
                               TRUE ~ parameter))

# DO data ####

#NOTE, IN 2021 DATA STORAGE CHANGES - BOTH BIO AND DO, AS WELL AS SOME OTHER PARAMS ARE STORED IN ONE FILE

# load files
raw_do_2017 <- read_xlsx(paste0(datadir, "1986-2017/DO.xlsx"), 
                          sheet="DO",
                          col_types = 'text') %>% 
  mutate(DATE = as.Date(as.numeric(DATE), origin = '1899-12-30'))
head(raw_do_2017)

raw_do_2018 <- read_csv(paste0(datadir, "2018/Sunapee 2018 DO.csv"),
                         col_types = cols(.default = col_character())) %>% 
  mutate(DATE = as.Date(DATE, format = '%d-%b-%y'))
head(raw_do_2018)

raw_do_2019 <- read_csv(paste0(datadir, '2019/Sunapee2019DO.csv'),
                         col_types = cols(.default = col_character()),
                         col_names = c('DATE', 'TIME', 'STATION', 'SENSOR', 'TEMP', 'PCNTSAT','DO',  'DEPTH'),
                        skip = 1)%>% 
  mutate(DATE = as.Date(DATE, format = '%m/%d/%Y'))
head(raw_do_2019)

raw_do_2020 <- read_csv(paste0(datadir, '2020/Sunapee2020DO.csv'),
                        col_types = cols(.default = col_character()),
                        skip = 1,
                        col_names = c('LAKE', 'TOWN','STATION', 'DATE', 'DEPTH', 'TEMP', 'DO', 'PCNTSAT', 'SPC_USCM', 'PH', 'NTU', 'TIME', 'BOTTOMZ', 'WEATHER', 'COMMENTS', 'ID', 'DATESTAID'))%>% 
  mutate(DATE = as.Date(DATE, format = '%d-%b-%y'))
head(raw_do_2020)

raw_do_2021 = read_xlsx(file.path(datadir, '2021-2022/LMPDO_2021_2022.xlsx'),
                        col_types = 'text') %>% 
  mutate(DATE = as.Date(as.numeric(DATE), origin = '1899-12-30')) %>% 
  rename(SPC_USCM = `SPC-uS/cm`,
         PH = pH)
head(raw_do_2021)

#### collate do data ####
comp_do <- full_join(raw_do_2017, raw_do_2018) %>% 
  full_join(., raw_do_2019) %>% 
  full_join(., raw_do_2020) %>% 
  full_join(., raw_do_2021)
head(raw_do)
unique(comp_do$Comments)
unique(comp_do$WEATHER)
unique(comp_do$BOTTOMZ)

#### pull out weather conditions ####
weather_obs <- comp_do %>% 
  select(DATE, STATION, WEATHER) %>% 
  rename(date = DATE,
         station = STATION,
         weather = WEATHER)  %>% 
  mutate(weather = case_when(grepl('not recorded', weather, ignore.case = T) ~ NA_character_,
                             TRUE ~ weather)) %>% 
  filter(!is.na(weather))
weather_obs <- unique(weather_obs) %>% 
  arrange(date) %>% 
  mutate(source = 'DO/Temperature record')

write_csv(weather_obs, paste0(dumpdir, 'weather_observations.csv'))
  

#### pull out bottom z ####
bottom_depth <- comp_do %>% 
  select(DATE, STATION, BOTTOMZ) %>% 
  rename(date = DATE,
         station = STATION,
         bottom_depth_m = BOTTOMZ) %>% 
  mutate(bottom_depth_m = case_when(as.numeric(bottom_depth_m) <= 2.5 ~ NA_real_, #0 and negative are often listed; 2.5 is listed at a deep site, which is incorrect
                                    TRUE ~ as.numeric(bottom_depth_m))) %>% 
  filter(!is.na(bottom_depth_m))

ggplot(bottom_depth, aes(x = date, y = bottom_depth_m)) +
  geom_point() +
  facet_grid(station ~ .)

#some of these look okay, but there are inconsistencies - not pushing these to the primary file at this time.

#select pertinent columns
qaqc_do <- comp_do %>% 
  select(STATION, DATE, DEPTH, TEMP, DO, PCNTSAT, SPC_USCM, PH, NTU, `Chl ug/L`, `Chl RFU`, ID)

#format depth, temp, do columns to numeric
qaqc_do <- qaqc_do %>% 
  mutate_at(vars(DEPTH, TEMP, DO, PCNTSAT, SPC_USCM, PH, NTU, `Chl ug/L`, `Chl RFU`),
            ~ as.numeric(.))

#### look at ranges and recode na and obviously errant data ####
range(qaqc_do$TEMP, na.rm = T)
# negative values are likely incorrect, will come back to this.

range(qaqc_do$DO, na.rm = T)
#negative do needs to be recoded
qaqc_do$DO [qaqc_do$DO < 0] = NA_real_
range(qaqc_do$DO, na.rm = T)
# high values are incorrect to, but need to see them in context.

range(qaqc_do$PCNTSAT, na.rm = T)
#negative do needs to be recoded
qaqc_do$PCNTSAT [qaqc_do$PCNTSAT < 0] = NA_real_
range(qaqc_do$PCNTSAT, na.rm = T)
#some really high values are impossible, need to see in context

range(qaqc_do$NTU, na.rm = T)
# recode negative values
qaqc_do$NTU [qaqc_do$NTU < 0] = NA_real_
range(qaqc_do$NTU, na.rm = T)

range(qaqc_do$PH, na.rm = T)
#0 will have to be recoded, see in context.

range(qaqc_do$SPC_USCM, na.rm = T)
#again, zero

range(qaqc_do$`Chl RFU`, na.rm = T)
#0s

range(qaqc_do$`Chl ug/L`, na.rm = T)
#0s

#### plot to look at data ####

##### do ####
ggplot(qaqc_do, aes(x = DATE, y = DO)) +
  geom_point()

#recode points above 75 - those are wrong, likely transcription errors
qaqc_do$DO [qaqc_do$DO >75] = NA_real_
ggplot(qaqc_do, aes(x = DATE, y = DO)) +
  geom_point()
#couple of oddballs, but look good other than that. 

#plot by year, color by station
years = seq(1986, 2022, by =1)

#Set up pdf device
pdf(file=paste0(figuredump, 'sunapee annual do.pdf'),width=11,height=8.5)
par()
for(i in 1:length(years)){
  DF <- qaqc_do %>% 
    filter(DATE >= paste(years[i], '01', '01', sep='-') &
             DATE < paste(years[i]+1, '01', '01', sep='-'))
  PLOT <- ggplot(DF, aes(x = DATE, y = DO, color = STATION)) +
    geom_point() +
    labs(title = years[i]) +
    coord_cartesian(ylim = c(0,max(qaqc_do$DO, na.rm = T)))
  print(PLOT)
}
dev.off()


# flags added May2022
qaqc_do <- qaqc_do %>% 
  mutate(year = format(DATE, '%Y')) %>% 
  mutate(do_flag = case_when(year == 1991 & DO < 2.5 ~ 'suspect',
                             DATE > as.Date('1992-01-01') & DATE < as.Date('1992-06-01') & DO <7.5 ~ 'suspect',
                             DATE > as.Date('1994-08-01') & DATE < as.Date('1994-09-01') & DO <2.5 ~ 'suspect',
                             DATE > as.Date('1995-01-01') & DATE < as.Date('1995-06-02') & DO <5 ~ 'suspect',
                             DATE > as.Date('1995-07-01') & DATE < as.Date('1995-08-01') & DO <5 ~ 'suspect',
                             DATE > as.Date('1996-06-01') & DATE < as.Date('1996-07-01') & DO <5 ~ 'suspect',
                             DATE > as.Date('1997-07-01') & DATE < as.Date('1997-08-01') & DO <2.5 ~ 'suspect',
                             DATE > as.Date('1998-06-01') & DATE < as.Date('1998-08-01') & DO <5 ~ 'suspect',
                             DATE > as.Date('1999-01-01') & DATE < as.Date('1999-06-01') & DO <5 ~ 'suspect',
                             DATE > as.Date('2001-01-01') & DATE < as.Date('2001-06-01') & DO >12.5 ~ 'suspect',
                             DATE > as.Date('2003-01-01') & DATE < as.Date('2003-06-01') & DO >7.5 ~ 'suspect',
                             TRUE ~ ''))
         
#Set up pdf device and overwrite previous file
pdf(file=paste0(figuredump, 'sunapee annual do.pdf'),width=11,height=8.5)
par()
for(i in 1:length(years)){
  DF <- qaqc_do %>% 
    filter(DATE >= paste(years[i], '01', '01', sep='-') &
             DATE < paste(years[i]+1, '01', '01', sep='-'))
  PLOT <- ggplot(DF, aes(x = DATE, y = DO, color = do_flag, shape = STATION)) +
    geom_point() +
    labs(title = years[i]) +
    coord_cartesian(ylim = c(0,max(qaqc_do$DO, na.rm = T)))
  print(PLOT)
}
dev.off()

##### percent saturation ####
ggplot(qaqc_do, aes(x = DATE, y = PCNTSAT)) +
  geom_point()

#percent saturation gt 250 are errant
qaqc_do$PCNTSAT [qaqc_do$PCNTSAT>250] = NA_real_

ggplot(qaqc_do, aes(x = DATE, y = PCNTSAT)) +
  geom_point()

# add flag to site 210 data on 1994-07-29 seems too high
qaqc_do <- qaqc_do %>% 
  mutate(do_flag = case_when(DATE == as.Date('1994-07-29') & do_flag == '' ~'likely do calibration issue on this date - do sat very high',
                             DATE == as.Date('2008-07-01') & do_flag == '' ~ 'possible do calibration issue on this date - do sat very high at 210',
                             DATE == as.Date('1994-07-29') & do_flag != '' ~  paste(do_flag, 'likely do calibration issue on this date - do sat very high', sep = '; '),
                             DATE == as.Date('2008-07-01') & do_flag != '' ~  paste(do_flag, 'possible do calibration issue on this date - do sat very high at 210', sep = '; '),
                             TRUE ~ do_flag)) 
ggplot(qaqc_do, aes(x = DATE, y = PCNTSAT, color = do_flag)) +
  geom_point()

#print year-by-year pdf
pdf(file=paste0(figuredump, 'sunapee annual dopct.pdf'),width=11,height=8.5)
par()
for(i in 1:length(years)){
  DF <- qaqc_do %>% 
    filter(DATE >= paste(years[i], '01', '01', sep='-') &
             DATE < paste(years[i]+1, '01', '01', sep='-'))
  PLOT <- ggplot(DF, aes(x = DATE, y = PCNTSAT, color = do_flag, shape = STATION)) +
    geom_point() +
    labs(title = years[i]) +
    coord_cartesian(ylim = c(0,max(qaqc_do$PCNTSAT, na.rm = T)))
  print(PLOT)
}
dev.off()

#additional data flags and recoding May 2022
qaqc_do <- qaqc_do %>% 
  mutate(do_flag = case_when(DATE > as.Date('1998-07-01') & DATE < as.Date('1998-07-15') & PCNTSAT <25 ~'suspect',
                             year == 2003 & PCNTSAT > 150 ~ 'suspect',
                             year == 2007 & PCNTSAT > 150 ~ 'suspect',
                             TRUE ~ do_flag))  %>% 
  mutate(PCNTSAT = case_when(DATE >= as.Date('1998-09-01') & DATE < as.Date('1998-10-01') & PCNTSAT < 5 ~ NA_real_, # these should not be zero
                             DATE == as.Date('1999-07-01') & PCNTSAT < 5 ~ NA_real_, # these should not be zero
                             TRUE ~ PCNTSAT))

#overwrite year-by-year pdf
pdf(file=paste0(figuredump, 'sunapee annual dopct.pdf'),width=11,height=8.5)
par()
for(i in 1:length(years)){
  DF <- qaqc_do %>% 
    filter(DATE >= paste(years[i], '01', '01', sep='-') &
             DATE < paste(years[i]+1, '01', '01', sep='-'))
  PLOT <- ggplot(DF, aes(x = DATE, y = PCNTSAT, color = do_flag, shape = STATION)) +
    geom_point() +
    labs(title = years[i]) +
    coord_cartesian(ylim = c(0,max(qaqc_do$PCNTSAT, na.rm = T)))
  print(PLOT)
}
dev.off()

##### temperature ####
ggplot(qaqc_do, aes(x = DATE, y = TEMP)) +
  geom_point()

#drop Feb 2 2022. All temps recorded below zero, there may be something malfunctioning.
remove = qaqc_do %>% 
  mutate(year = format(DATE, '%Y'), month = format(DATE, '%m'), day = format(DATE, '%d')) %>% 
  filter(year == 2022 & month == "02" & day == '02') %>% 
  select(-c(year, month, day))

qaqc_do = anti_join(qaqc_do, remove)

#values where temp less than 2 are incorrect, recoding the do data as well, since those are temp based
ix = which(qaqc_do$TEMP<2 & format(qaqc_do$DATE, '%Y') < 2020)

qaqc_do$TEMP[ix] = NA_real_
qaqc_do$DO[ix] = NA_real_
qaqc_do$PCNTSAT[ix] = NA_real_

ggplot(qaqc_do, aes(x = DATE, y = TEMP)) +
  geom_point()

#recode oddball temp at site 200, 1999-06-28, 10m, <4 degC
qaqc_do$TEMP [qaqc_do$TEMP < 4 & qaqc_do$STATION == 200 & qaqc_do$DEPTH == 10] = NA_real_

ggplot(qaqc_do, aes(x = DATE, y = TEMP)) +
  geom_point()

##### pH ####
ggplot(qaqc_do, aes(x = DATE, y = PH)) +
  geom_point()

#recode 0 to na
qaqc_do$PH = ifelse(qaqc_do$PH == 0, NA_real_, qaqc_do$PH)

# calculate H+ ion from pH (H+ = 10^-pH) and drop pH column
qaqc_do$conc_H_molpl=10^(qaqc_do$PH * -1)
qaqc_do$PH = NULL

##### turbidity ####
ggplot(qaqc_do, aes(x = DATE, y = NTU)) +
  geom_point()

#LOOKS LIKE THERE ARE SOME PROBLEMATIC DATA HERE, LIKELY SONDE IS IN SEDIMENT - want to look at these relative to the DO data to make sure those data are recoded properly.


##### conductivity ####
ggplot(qaqc_do, aes(x = DATE, y = SPC_USCM)) +
  geom_point()

#recode 0 to na 
qaqc_do$SPC_USCM = ifelse(qaqc_do$SPC_USCM == 0, NA_real_, qaqc_do$SPC_USCM)

#AGAIN, PROBABLY SOME PROMLEMATIC DATA HERE

#### chl ugl and rfu ----
ggplot(qaqc_do, aes(x = DATE, y = `Chl ug/L`)) +
  geom_point()

# need to recode zeros on May date
qaqc_do$`Chl ug/L` = ifelse(qaqc_do$`Chl ug/L` == 0, NA_real_, qaqc_do$`Chl ug/L`)
qaqc_do$`Chl RFU` = ifelse(qaqc_do$`Chl RFU` == 0, NA_real_, qaqc_do$`Chl RFU`)

# and flag rfu > 1 as suspect
qaqc_do$chl_flag = ifelse(qaqc_do$`Chl RFU`>1, 'suspect', '')

#### look at 2020 data to see better the relationship between parameters that are suspected in the sediment
do20 <- qaqc_do %>% 
  filter(DATE > '2020-01-01') %>% 
  pivot_longer(cols = c(TEMP:NTU, `Chl ug/L`, `Chl RFU`), names_to = 'variable')

ggplot(do20, aes(x = DATE, y = value, color = DEPTH))+
  geom_point() +
  facet_grid(variable ~ ., scales = 'free_y') +
  theme_bw()

#looks like where turbidity is above 100, the sensor is in the sediment.
qaqc_do <- qaqc_do %>% 
  mutate_at(vars(TEMP:NTU, conc_H_molpl, `Chl ug/L`, `Chl RFU`),
            ~case_when(NTU>100 ~ NA_real_,
                       TRUE ~ .))
#look again
do20 <- qaqc_do %>% 
  filter(DATE > '2020-01-01') %>% 
  pivot_longer(cols = c(TEMP:NTU, `Chl ug/L`, `Chl RFU`), names_to = 'variable')

ggplot(do20, aes(x = DATE, y = value, color = DEPTH))+
  geom_point() +
  facet_grid(variable ~ ., scales = 'free_y') +
  theme_bw()

# flagging anything greater than 10 NTU as suspect data
qaqc_do <- qaqc_do %>% 
  mutate(gen_flag = case_when(NTU>10 ~ 'sonde suspected in sediment',
                              TRUE ~ ''))


#### rename columns ####
qaqc_do <- qaqc_do %>% 
  rename(station = STATION,
         depth_m = DEPTH,
         waterTemperature_degC = TEMP,
         oxygenDissolved_mgl = DO,
         oxygenDissolvedPercentOfSaturation_pct = PCNTSAT,
         turbidity_NTU = NTU,
         specificConductance_uScm = SPC_USCM,
         chlorophyllFluorescence_ugl = `Chl ug/L`,
         chlorophyllFluorescence_RFU = `Chl RFU`,
         date = DATE,
         org_id = ID) %>% 
  mutate(station = as.numeric(station))
str(qaqc_do)

#### create vertical dataset ####
qaqc_do_vert <- qaqc_do %>% 
  select(-year) %>% 
  gather(parameter, value, -station, -depth_m, -date, -do_flag, -chl_flag, -gen_flag, -org_id) %>% 
  filter(!is.na(value)) %>% 
  mutate(site_type = 'lake',
         depth_m = round(depth_m, digits = 1)) %>% 
  mutate(flag = '') %>% 
  mutate(flag = case_when(parameter == 'oxygenDissolved_mgl' ~ do_flag,
                          parameter == 'oxygenDissolvedPercentOfSaturation_pct' ~ do_flag,
                          parameter == 'chlorophyllFluorescence_ugl' | parameter == 'chlorophyllFluorescence_RFU' ~ chl_flag,
                          TRUE ~ '')) %>% 
  select(-do_flag, -chl_flag) 

ggplot(qaqc_do_vert, aes(x = date, y = value, color = flag)) +
  geom_point() +
  facet_grid(parameter ~ station, scales = 'free_y') +
  theme_bw()


# JOIN AND HARMONIZE ALL DATA ####
head(qaqc_chem_vert)
unique(qaqc_chem_vert$parameter)
head(qaqc_bio_vert)
unique(qaqc_bio_vert$parameter)
head(qaqc_do_vert)
unique(qaqc_do_vert$parameter)

primary_file_vert <- qaqc_chem_vert %>% 
  mutate(station = as.numeric(station)) %>% 
  full_join(., qaqc_bio_vert) %>% 
  full_join(., qaqc_do_vert) %>% 
  mutate(station = as.character(station))

head(qaqc_chem)
head(qaqc_bio)
head(qaqc_do)

unique(primary_file_vert$parameter)
unique(primary_file_vert$site_type)

primary_file_tributary <- primary_file_vert %>% 
  filter(site_type == 'tributary') 
primary_file_lake <- primary_file_vert %>% 
  filter(site_type == 'lake')

#### qaqc with whole dataset ####
# plot turbidity and do together
lake_do <- primary_file_lake %>% 
  filter(parameter == 'oxygenDissolved_mgl') %>% 
  rename(DO_mgl = value,
         do_flag = flag) %>% 
  select(-parameter, -layer, - org_id, -org_sampid)
lake_dosat<- primary_file_lake %>% 
  filter(parameter == 'oxygenDissolvedPercentOfSaturation_pct') %>% 
  rename(DO_pctsat = value,
         dosat_flag = flag) %>% 
  select(-parameter, -layer, - org_id, -org_sampid)
lake_turb<- primary_file_lake %>% 
  filter(parameter == 'turbidity_NTU') %>% 
  rename(turb_NTU = value,
         turb_flag = flag) %>% 
  select(-parameter, -layer, - org_id, -org_sampid)
do_turb <- full_join(lake_do, lake_dosat) %>% 
  full_join(., lake_turb)

years_doturb = seq(2001, 2022, by = 1)

pdf(file=paste0(figuredump, 'annual do and turbidity by year.pdf'),width=11,height=8.5)
par()
for(i in 1:length(years_doturb)) {
  PLOT <- do_turb %>%
    filter(as.numeric(station) >= 200) %>%
    filter(date >= as.Date(paste(years_doturb[i], '01', '01', sep='-')) &
             date < as.Date(paste(years_doturb[i]+1, '01', '01', sep='-'))) %>%
    ggplot(., aes(x = date, y = DO_pctsat, size = turb_NTU)) +
    geom_point()+
    facet_grid(station ~ .) +
    theme(legend.position = 'bottom') +
    labs(title = years_doturb[i])
  print(PLOT)
}
dev.off()

#### plot the data together ####
lakevars = unique(primary_file_lake$parameter)
pdf(file=paste0(figuredump, 'sunapee whole record by parameter.pdf'),width=11,height=8.5)
par()
for(i in 1:length(lakevars)) {
  DF <- primary_file_lake %>%
    filter(parameter == lakevars[i])
  PLOT <- ggplot(DF, aes(x = date, y = value, color = station, shape = flag)) +
    geom_point() +
    facet_grid(parameter ~ ., scales = 'free_y') +
    theme(legend.position = 'bottom')
  print(PLOT)
}
dev.off()

# SAVE primary FILE ####
primary_file_vert <- primary_file_vert %>% 
  mutate(gen_flag = case_when(is.na(gen_flag) ~ '',
                              TRUE ~ gen_flag),
         layer = case_when(is.na(layer) ~ '', 
                           TRUE ~ layer)) %>% 
  select(station, date, depth_m, layer, site_type, 
         parameter, value, flag, gen_flag, org_sampid, org_id)

write_csv(primary_file_vert, paste0(dumpdir, 'LSPALMP_1986-2022_v', Sys.Date(), '.csv'))

# COLLATE STATION LOCATIONS AND CREATE PARAMETER SUMMARIES ####

#create primary list of stations in the primary file
station <- unique(primary_file_vert$station) 
station <- tibble(station)

#add tributary/lake identifier as well as site type and sub type for lake sites
station_details <- station %>% 
  mutate(station = as.numeric(station)) %>% 
  mutate(site_type = case_when((station) <= 230 ~ 'lake',
                              TRUE ~ 'tributary')) %>% 
  mutate(sub_site_type = case_when(site_type == 'lake' & (station) >=200 ~ 'deep',
                                   site_type == 'lake' & station <200 ~ 'cove',
                                   TRUE ~ '')) 

# make a summary of first sample, last sample per parameter per site, export as summary file in primary files
parameter_summary <- primary_file_vert %>% 
  group_by(station, parameter) %>% 
  summarize(first_sample = min(date),
            last_sample = max(date),
            n_obs = length(value)) %>% 
  mutate(station = as.numeric(station)) %>% 
  full_join(station_details)

write_csv(parameter_summary, paste0(dumpdir, 'parameter_by_site_sample_summary.csv'))


# summarize year of first and last sample by station to summarize station status
station_summary <- primary_file_vert %>% 
  group_by(station) %>% 
  summarize(first_year = as.numeric(format(min(date), '%Y')),
            last_year = as.numeric(format(max(date), '%Y'))) %>% 
  mutate(station = as.numeric(station))

station_details <- full_join(station_details, station_summary) %>% 
  mutate(status = case_when(last_year >= 2020 ~ 'ongoing',
                            TRUE ~ 'inactive'),
         status = case_when(first_year == 2020 ~ 'temporary',
                            TRUE ~ status))

#add lat long information as available
templake <- read_csv(paste0(datadir, 'station locations/VLAP_LS_AddedDeepStations2020.csv'))
temptributary <- read_csv(paste0(datadir, 'station locations/VLAP_LS_AddedTribStations2020.csv'))
tributary <- read_xls(paste0(datadir, 'station locations/tributary_pts_alt.xls')) %>% 
  select(stream_no, lat_dd, lon_dd) %>% 
  rename(station = stream_no)  %>% 
  mutate(lat_dd = round(lat_dd, digits = 4),
         lon_dd = round(lon_dd, digits = 4))
lake <- read_xls(paste0(datadir, 'station locations/Cove  Deep + Buoy WQ Sample Point GPS Coordinates.xls'),
                 skip = 3) %>% 
  mutate(station = as.numeric(WAYPOINT)) %>% #warnings okay - this is intentional character -> numeric
  filter(!is.na(station)) %>% 
  mutate(lat_dd = round(Y, digits = 4),
         lon_dd = round(X, digits = 4)) %>% 
  select(station, lat_dd, lon_dd)
  
temp <- full_join(templake, temptributary) %>% 
  select(AliasID, LATDECDEG, LONGDECDEG) %>% 
  rename(station = AliasID, 
         lat_dd = LATDECDEG,
         lon_dd = LONGDECDEG) %>% 
  mutate(lat_dd = round(lat_dd, digits = 4),
         lon_dd = round(lon_dd, digits = 4))

station_locs <- full_join(tributary, lake) %>% 
  full_join(temp)

station_details <- full_join(station_details, station_locs)

#now bring in historical site locations from OneStop (NHDES Onestop Mapper)
onestop <- read_delim(paste0(datadir, 'station locations/VLAP_LS_allinWS.txt'), delim = ',') %>%
  filter(startsWith(STATNAME,'SUNAPEE LAKE'))
onestop <- onestop %>% 
  mutate(site_char = substr(STATNAME, 14, 20)) %>% 
  mutate(site_num = as.numeric(site_char)) %>% 
  mutate(site_num = case_when(is.na(site_num) ~ as.numeric(substr(site_char, 1, 4)),
                              TRUE~ site_num)) %>% 
  filter(!is.na(site_num)) %>% 
  select(site_num, LONGITUDE, LATITUDE) %>% 
  rename(station = site_num) %>% 
  mutate(LONGITUDE = round(LONGITUDE, digits = 4),
         LATITUDE = round(LATITUDE, digits = 4))

station_details <- full_join(station_details, onestop) %>% 
  mutate(lat_dd = case_when(is.na(lat_dd) ~ LATITUDE,
                            TRUE ~ lat_dd),
         lon_dd = case_when(is.na(lon_dd) ~ LONGITUDE,
                            TRUE ~ lon_dd)) %>% 
  select(-LATITUDE, -LONGITUDE) %>% 
  filter(!is.na(site_type))

write_csv(station_details, paste0(dumpdir, 'station_location_details.csv'))
