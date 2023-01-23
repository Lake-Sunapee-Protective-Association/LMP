# LSPA LMP primary files

This folder contains summary data tables, station location details, and collated, QAQC'd primary files. 

***

# Data files

## [LSPALMP_1986-2022_v2023-01-22.csv](https://github.com/Lake-Sunapee-Protective-Association/LMP/blob/main/primary%20files/LSPALMP_1986-2022_v2023-01-22.csv)

This is the collated file of all LSPA LMP data as collated and QAQC'd in the code found in the [collation code](http://github.com/Lake-Sunapee-Protective-Association/LMP/tree/main/collation%20code) folder from 1986-2022. These data have been collated, QAQC'd to recode obviously errant data, and to flag suspicious data.

| 	Column Name		|	Column Description															|	Units	|	Format 	|
| 	----			|	----																		|	----	|	----	|
|	station	|	station identifier, see [station_location_details.csv](https://github.com/Lake-Sunapee-Protective-Association/LMP/blob/main/primary%20files/station_location_details.csv) for location lat-longs	|	*N/A*	|	numeric	|
|	date	|	date of observation, measurement, or sample	|	*N/A*	|	YYYY-MM-DD	|
|	depth_m	|	depth of observation, measurement, or sample if applicable	|	meters	|	numeric	|
|	layer	|	layer from which an observation, measurement, or sample was taken, if applicable	|	I = integrated, E = epilimnion, M = metalimnion, H = hypolimnion	|	character	|
|	site_type	|	indication of 'lake' site or 'stream' site	|	*N/A*	| character string	|
|	parameter	| 	parameter for which *value* is reported, see table below for parameter descriptions and units	|	*N/A*	|	character string	|
|	value	|	value of the *parameter* measured at *station* on *date* at *depth_m* or *layer*	|	various, see parameter units below	|	numeric	|
|	flag	|	use flag for *parameter*-*value* reported	|	*N/A*	|	character string, 'BDL' = 'below detection limit', all flags indicate use-with-caution	|	
|	org_sampid	|	where 'SampleID' was listed in the original raw data file, it was copied to this column	|	*N/A*	|	character	|	
|	org_id	|	where 'ID' was listed in the original raw data file, it was copied to this column	|	*N/A*	|	character	|


Additional information for the parameter column, names use ODM2 variable names and units where applicable:
| 	Parameter Name	|	Parameter Description							|	Parameter Units	|	instrument	|	location of measurement	|
| 	----			|	----											|	----	| 	----	|	----	|
|	alkalinity_mglCaCO3	|	water alkalinity								| 	milligramsPerLiter |	VWR SympHony B10P	|	CSC	|
|	chlorphyll_a_ugl		|	concentration of chlorophyll-*a* in water		|	microgramsPerLiter	|	ThermoScientific Genesys 30	|	CSC	|
|	chloride_mgl			|	concentration of chloroide in water				|	milligramsPerLiter	|	Orion VersaStar Pro Meter	|	CSC	|
|	dissolvedHydrogen_moll	|	concentration of hydrogen ions in water			|	molesPerLiter	|	VWR SympHony B10P (pre-2019), then YSI ProDSS	|	CSC, then *in-situ*	|
|	specificConductance_uScm		|	conductivity of water							|	microSiemensPerCentimeter	|	Orion 3Star Meter (pre-2019), then YSI ProDSS	|	CSC, then *in-situ*	|
|	secchiDepth_m	|	Secchi depth									| 	meters	|	Secchi disc	|	*in-situ*	|
|	totalPhosphorus_mgl			| 	concentration of total phosphorus in water		|	milligramsPerLiter	|	ThermoScientific Genesys 30	|	CSC	|
|	turbidity_NTU		|	turbidity of water								|	nephelometricTurbidityUnit	|	 HF Scientific Micro 100 (pre-2019), then YSI ProDSS	|	CSC, then *in-situ*	|
|	oxygenDissolved_mgl			|	dissolved oxygen concentration in water			| 	milligramsPerLiter	|	unkown sonde (pre-2019), then YSI ProDSS	|	*in-situ*	|
|	oxygenDissolvedPercentOfSaturation_pct		|	dissolved oxygen saturation in water			| 	percent	|	unkown sonde (pre-2019), then YSI ProDSS	|	*in-situ*	|
|	waterTemperature_degC			|	temperature of water							|	degreesCelsius	|	unkown sonde (pre-2019), then YSI ProDSS	|	*in-situ*	|
|	chlorophyllFluorescence_ugl			|	chlorophyll estimated from *in situ* fluorescence						|	microGramsPerLiter	|	YSI ProDSS	|	*in-situ*	|
|	chlorophyllFluorescence_RFU		|	chlorophyll fluroescence 					|	relativeFluorescenceUnits	|	YSI ProDSS	|	*in-situ*	|

Sample analysis methodology questions should be directed to the LSPA. 'CSC' = LSPA lab at Colby Sawyer College

#### [station_location_details.csv](https://github.com/Lake-Sunapee-Protective-Association/LMP/blob/main/primary%20files/station_location_details.csv)

This file contains all station location information available from the LSPA. Note that not all stations have latitude/longitude information stored with them.

| 	Column Name		|	Column Description															|	Units	|	Format 	|
| 	----			|	----																		|	----	|	----	|
|	station			|	station identifier															| 	*N/A*	|	numeric	|
|	site_type		|	site type indicating whether sample obtained in "lake" or "stream"			|	*N/A*	|	character	|
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
|	site_type		|	site type indicating whether sample obtained in "lake" or "stream"			|	*N/A*	|	character	|
|	sub_site_type	| 	"cove" (shallow) or "deep" location; applicable only for "lake" site_type	| 	*N/A*	|	character	|


Additional information for the parameter column:
| 	Parameter Name	|	Parameter Description							|	Parameter Units	|	instrument	|	location of measurement	|
| 	----			|	----											|	----	| 	----	|	----	|
|	alkalinity_mglCaCO3	|	water alkalinity								| 	milligramsPerLiter |	VWR SympHony B10P	|	CSC	|
|	chlorophyll_a_ugl		|	concentration of chlorophyll-*a* in water		|	microgramsPerLiter	|	ThermoScientific Genesys 30	|	CSC	|
|	chloride_mgl			|	concentration of chloroide in water				|	milligramsPerLiter	|	Orion VersaStar Pro Meter	|	CSC	|
|	dissolvedHydrogen_moll	|	concentration of hydrogen ions in water			|	molesPerLiter	|	VWR SympHony B10P (pre-2019), then YSI ProDSS	|	CSC, then *in-situ*	|
|	specificConductance_uScm		|	specific conductance of water							|	microSiemensPerCentimeter	|	Orion 3Star Meter (pre-2019), then YSI ProDSS	|	CSC, then *in-situ*	|
|	secchiDepth_m	|	Secchi depth									| 	meters	|	Secchi disc	|	*in-situ*	|
|	totalPhosphorus_mgl			| 	concentration of total phosphorus in water		|	milligramsPerLiter	|	ThermoScientific Genesys 30	|	CSC	|
|	turbidity_NTU		|	turbidity of water								|	nephelometricTurbidityUnit	|	 HF Scientific Micro 100 (pre-2019), then YSI ProDSS	|	CSC, then *in-situ*	|
|	oxygenDissolved_mgl			|	dissolved oxygen concentration in water			| 	milligramsPerLiter	|	unkown sonde (pre-2019), then YSI ProDSS	|	*in-situ*	|
|	oxygenDissolvedPercentOfSaturation_pct		|	dissolved oxygen saturation in water			| 	percent	|	unkown sonde (pre-2019), then YSI ProDSS	|	*in-situ*	|
|	waterTemperature_degC			|	temperature of water							|	degreesCelsius	|	unkown sonde (pre-2019), then YSI ProDSS	|	*in-situ*	|
|	chlorophyllFluorescence_ugl			|	chlorophyll estimated from *in situ* fluorescence						|	microGramsPerLiter	|	YSI ProDSS	|	*in-situ*	|
|	chlorophyllFluorescence_RFU		|	chlorophyll fluroescence 					|	relativeFluorescenceUnits	|	YSI ProDSS	|	*in-situ*	|

Sample analysis methodology questions should be directed to the LSPA. 'CSC' = LSPA lab at Colby Sawyer College


#### [weather_observations.csv](https://github.com/Lake-Sunapee-Protective-Association/LMP/blob/main/primary%20files/weather_observations.csv)

This file contains a summary of unique date/station pairs where weather conditions were reported in the LSPA's LMP files. 

|	Column Name	|	description	|
|	----	|	----	|
|	date	|	date of observation in YYYY-MM-DD format	|
|	station	|	location of observation, often this is the same at all sites, see [station_location_details.csv](https://github.com/Lake-Sunapee-Protective-Association/LMP/blob/main/primary%20files/station_location_details.csv) for location lat-longs	|
|	weather	source	|	raw file stream from which the weather data were collated	|
