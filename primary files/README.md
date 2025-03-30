# LSPA LMP primary files

This folder contains summary data tables, station location details, and collated, QAQC'd primary files. 

***

# Data files

## [LSPALMP_1986-2023_v2024-01-20.csv](https://github.com/Lake-Sunapee-Protective-Association/LMP/blob/main/primary%20files/LSPALMP_1986-2023_v2024-01-20.csv)

This is the collated file of all LSPA LMP data as collated and QAQC'd in the code found in the [collation code](http://github.com/Lake-Sunapee-Protective-Association/LMP/tree/main/collation%20code) folder from 1986-2023. These data have been collated, QAQC'd to recode obviously errant data, and to flag suspicious data.

| 	Column Name		|	Column Description															|	Units	|	Format 	|
| 	----			|	----																		|	----	|	----	|
|	station	|	station identifier, see [station_location_details.csv](https://github.com/Lake-Sunapee-Protective-Association/LMP/blob/main/primary%20files/station_location_details.csv) for location lat-longs	|	*N/A*	|	numeric	|
|	date	|	date of observation, measurement, or sample	|	*N/A*	|	YYYY-MM-DD	|
|	depth_m	|	depth of observation, measurement, or sample if applicable	|	meters	|	numeric	|
|	layer	|	layer from which an observation, measurement, or sample was taken, if applicable	|	I = integrated, E = epilimnion, M = metalimnion, H = hypolimnion	|	character	|
|	site_type	|	indication of 'lake' site or 'tributary' site	|	*N/A*	| character string	|
|	parameter	| 	parameter for which *value* is reported, see table below for parameter descriptions and units	|	*N/A*	|	character string	|
|	value	|	value of the *parameter* measured at *station* on *date* at *depth_m* or *layer*	|	various, see parameter units below	|	numeric	|
|	flag	|	use flag for *parameter*-*value* reported	|	*N/A*	|	character string, 'BDL' = 'below detection limit', all flags indicate use-with-caution	|	
| gen_flag | any flag for *parameter*-*value* that applies to more than one parameter | *N/A* | character |
|	org_sampid	|	where 'SampleID' was listed in the original raw data file, it was copied to this column	|	*N/A*	|	character	|	
|	org_id	|	where 'ID' was listed in the original raw data file, it was copied to this column	|	*N/A*	|	character	|


Additional information for the parameter column, names use ODM2 variable names and units where applicable:

| 	Parameter Name	|	Parameter Description							|	Parameter Units	|	instrument	|	location of measurement	|
| 	----			|	----											|	----	| 	----	|	----	|
|	alkalinity_milligramCaCO3PerLiter	|	water alkalinity								| 	milligramPerLiter |	VWR SympHony B10P	|	CSC	|
|	chlorophyll_a_microgramPerLiter		|	concentration of chlorophyll-*a* in water		|	microgramPerLiter	|	ThermoScientific Genesys 30	|	CSC	|
|	chloride_milligramPerLiter			|	concentration of chloroide in water				|	milligramPerLiter	|	Orion VersaStar Pro Meter	|	CSC	|
|	hydrogenDissolved_molePerLiter	|	concentration of hydrogen ions in water			|	molePerLiter	|	VWR SympHony B10P (pre-2019), then YSI ProDSS	|	CSC, then *in-situ*	|
|	specificConductance_microSiemenPerCentimeter		|	conductivity of water							|	microSiemenPerCentimeter	|	Orion 3Star Meter (pre-2019), then YSI ProDSS	|	CSC, then *in-situ*	|
|	secchiDepth_meter	|	Secchi depth, note that an indication of 'scope' in *flag* column indicates a scope was used for measurement									| 	meter	|	Secchi disc	|	*in-situ*	|
|	phosphorusTotal_milligramPerLiter			| 	concentration of total phosphorus in water		|	milligramPerLiter	|	ThermoScientific Genesys 30	|	CSC	|
|	turbidity_nephelometricTurbidityUnit		|	turbidity of water								|	nephelometricTurbidityUnit	|	 HF Scientific Micro 100 (pre-2019), then YSI ProDSS	|	CSC, then *in-situ*	|
|	oxygenDissolved_milligramPerLiter			|	dissolved oxygen concentration in water			| 	milligramPerLiter	|	unkown sonde (pre-2019), then YSI ProDSS	|	*in-situ*	|
|	oxygenDissolvedPercentOfSaturation_percent		|	dissolved oxygen saturation in water			| 	percent	|	unkown sonde (pre-2019), then YSI ProDSS	|	*in-situ*	|
|	waterTemperature_degreeCelsius			|	temperature of water							|	degreeCelsius	|	unkown sonde (pre-2019), then YSI ProDSS	|	*in-situ*	|
|	chlorophyllFluorescence_microgramPerLiter			|	chlorophyll estimated from *in situ* fluorescence						|	microGramPerLiter	|	YSI ProDSS	|	*in-situ*	|
|	chlorophyllFluorescence_relativeFluorescenceUnit		|	chlorophyll fluroescence 					|	relativeFluorescenceUnit	|	YSI ProDSS	|	*in-situ*	|
| blue_GreenAlgae_Cyanobacteria_Phycocyanin_microgramPerLiter | BGA phycocyanin estimated from *in situ* fluorescence | microGramPerLiter | YSI ProDSS	|	*in-situ*	|
| blue_GreenAlgae_Cyanobacteria_Phycocyanin_relativeFluorescenceUnit | BGA phycocyanin as relative fluorescence units |relativeFluorescenceUnit  | YSI ProDSS	|	*in-situ*	|

Sample analysis methodology questions should be directed to the LSPA. 'CSC' = LSPA lab at Colby Sawyer College

#### [station_location_details.csv](https://github.com/Lake-Sunapee-Protective-Association/LMP/blob/main/primary%20files/station_location_details.csv)

This file contains all station location information available from the LSPA. Note that not all stations have latitude/longitude information stored with them.

| 	Column Name		|	Column Description															|	Units	|	Format 	|
| 	----			|	----																		|	----	|	----	|
|	station			|	station identifier															| 	*N/A*	|	numeric	|
|	site_type		|	site type indicating whether sample obtained in "lake" or "tributary"			|	*N/A*	|	character	|
|	sub_site_type	| 	"cove" (shallow) or "deep" location; applicable only for "lake" site_type	| 	*N/A*	|	character	|
|	first_year		|	first year station was sampled at											|	*N/A*	|	numeric YYYY	|
|	last_year		|	last year station was sampled at											|	*N/A*	|	numeric YYYY	|
|	status			|	indication of whether station sampling is 'ongoing', 'temporary', or 'inactive' |	*N/A*	|	character	|
|	lat_dd			|	latitude of station															| decimal degrees | numeric, WGS84 |
|	lon_dd			|	longitude of station														| decimal degrees |	numeric, WGS84 |	

Status of 'temporary' is for station(s) at which we do not anticipate continued sampling. These stations will move to 'inactive' once sampling has been completed. 


#### [parameter_by_site_sample_summary.csv](https://github.com/Lake-Sunapee-Protective-Association/LMP/blob/main/primary%20files/parameter_by_site_sample_summary.csv)

This file contains a brief summary of when/how long each parameter was measured at any given site.

| 	Column Name		|	Column Description															|	Units	|	Format 	|
| 	----			|	----																		|	----	|	----	|
|	station			|	station identifier, see [station_location_details.csv](https://github.com/Lake-Sunapee-Protective-Association/LMP/blob/main/primary%20files/station_location_details.csv) for location lat-longs	| 	*N/A*	|	numeric	|
|	parameter		|	parameter reported in primary file (see list below)							| various, see table below	| numeric	|
|	first_sample	| 	date of first sample of the *parameter* at a given *station*				|	*N/A*	|	date, YYYY-MM-DD |
|	last_sample		|	date of first sample of the *parameter* at a given *station*				|	*N/A*	|	date, YYYY-MM-DD |
|	n_obs			| 	total number of observations the given *parameter* at a given *station*		|	count	|	numeric	|
|	site_type		|	site type indicating whether sample obtained in "lake" or "tributary"			|	*N/A*	|	character	|
|	sub_site_type	| 	"cove" (shallow) or "deep" location; applicable only for "lake" site_type	| 	*N/A*	|	character	|


Additional information for the parameter column are provided in the section above


#### [weather_observations.csv](https://github.com/Lake-Sunapee-Protective-Association/LMP/blob/main/primary%20files/weather_observations.csv)

This file contains a summary of unique date/station pairs where weather conditions were reported in the LSPA's LMP files. 

|	Column Name	|	description	|
|	----	|	----	|
|	date	|	date of observation in YYYY-MM-DD format	|
|	station	|	location of observation, often this is the same at all sites, see [station_location_details.csv](https://github.com/Lake-Sunapee-Protective-Association/LMP/blob/main/primary%20files/station_location_details.csv) for location lat-longs	|
|	weather	source	|	raw file stream from which the weather data were collated	|
