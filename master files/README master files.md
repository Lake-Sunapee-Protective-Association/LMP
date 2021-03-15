# LSPA LMP master files

This folder contains summary data tables, station location details, and collated, QAQC'd master files. 

***

# Data files

#### [station_location_details.csv](https://github.com/Lake-Sunapee-Protective-Association/LMP/blob/main/master%20files/station_location_details.csv)

This file contains all station location information available from the LSPA and NH VLAP. Note that not all stations have latitude/longitude information stored with them.

| 	Column Name		|	Column Description															|	Units	|	Format 	|
| 	----			|	----																		|	----	|	----	|
|	station			|	station identifier															| 	*N/A*	|	numeric	|
|	site_type		|	site type indicating whether sample obtained in "lake" or "stream"			|	*N/A*	|	character	|
|	sub_site_type	| 	"cove" (shallow) or "deep" location; applicable only for "lake" site_type	| 	*N/A*	|	character	|
|	first_year		|	first year station was sampled												|	*N/A*	|	numeric YYYY	|
|	last_year		|	last year station was sampled												|	*N/A*	|	numeric YYYY	|
|	status			|	indication of whether station sampling is 'ongoing', 'temporary', or 'inactive' |	*N/A*	|	character	|
|	lat_dd			|	latitude of station															| decimal degrees | numeric, WGS84 |
|	lon_dd			|	longitude of station														| decimal degrees |	numeric, WGS84 |	

Status of 'temporary' is for station(s) that we do not anticipate continued sampling at. These stations will move to 'inactive' once sampling has been completed. 


#### [parameter_by_site_sample_summary.csv](https://github.com/Lake-Sunapee-Protective-Association/LMP/blob/main/master%20files/parameter_by_site_sample_summary.csv)

This file contains a brief summary of when/how long each parameter was measured at any given site.

| 	Column Name		|	Column Description															|	Units	|	Format 	|
| 	----			|	----																		|	----	|	----	|
|	station			|	station identifier															| 	*N/A*	|	numeric	|
|	parameter		|	parameter reported in master file (see list below)							| various, see table below	| numeric	|
|	first_sample	| 	date of first sample of the *parameter* at a given *station*				|	*N/A*	|	date, YYYY-MM-DD |
|	last_sample		|	date of first sample of the *parameter* at a given *station*				|	*N/A*	|	date, YYYY-MM-DD |
|	n_obs			| 	total number of observations the given *parameter* at a given *station*		|	count	|	numeric	|
|	site_type		|	site type indicating whether sample obtained in "lake" or "stream"			|	*N/A*	|	character	|
|	sub_site_type	| 	"cove" (shallow) or "deep" location; applicable only for "lake" site_type	| 	*N/A*	|	character	|

Additional information for the parameter column:
| 	Parameter Name	|	Parameter Description							|	Parameter Units	|
| 	----			|	----											|	----	|
|	alk_mglCaCO3	|	water alkalinity								| 	milligramsPerLiter |
|	chla_ugl		|	concentration of chlorophyll-*a* in water		|	microgramsPerLiter	|
|	cl_mgl			|	concentration of chloroide in water				|	milligramsPerLiter	|
|	conc_H_molpL	|	concentration of hydrogen ions in water			|	molesPerLiter	|
|	cond_uScm		|	conductivity of water							|	microSiemensPerCentimeter	|
|	secchidepth_m	|	secchi depth									| 	meters	|
|	TP_mgl			| 	concentration of total phosphorus in water		|	milligramsPerLiter	|
|	turb_NTU		|	turbidity of water								|	nephelometricTurbidityUnit	|
|	DO_mgl			|	dissolved oxygen concentration in water			| 	milligramsPerLiter	|
|	DO_pctsat		|	dissolved oxygen saturation in water			| 	percent	|
|	temp_C			|	temperature of water							|	degreesCelsius	|


