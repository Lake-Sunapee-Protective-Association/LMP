# LSPA LMP master files

This folder contains summary data tables, station location details, and collated, QAQC'd master files. 

***

# Data files

#### [station_location_details.csv](https://github.com/Lake-Sunapee-Protective-Association/LMP/blob/main/master%20files/station_location_details.csv)

This file contains all station location information available from the LSPA and NH VLAP. Note that not all stations have latitude/longitude information stored with them.

| 	Column Name		|	Column Description															|	Units	|	Format 	|
| 	----			|	----																		|	----	|	----	|
|	station			|	station identifier															| 	*N/A*	|	*N/A*	|
|	site_type		|	site type indicating whether sample obtained in "lake" or "stream"			|	*N/A*	|	*N/A*	|
|	sub_site_type	| 	"cove" (shallow) or "deep" location; applicable only for "lake" site_type	| 	*N/A*	|	*N/A*	|
|	first_year		|	first year station was sampled												|	*N/A*	|	'YYYY'	|
|	last_year		|	last year station was sampled												|	*N/A*	|	'YYYY'	|
|	status			|	indication of whether station sampling is 'ongoing', 'temporary', or 'inactive' |	*N/A*	|	*N/A*	|
|	lat_dd			|	latitude of station															| decimal degrees |	WGS84 |
|	lon_dd			|	longitude of station														| decimal degrees |	WGS84 |	

Status of 'temporary' is for station(s) that we do not anticipate continued sampling at. These stations will move to 'inactive' once sampling has been completed. 


#### [parameter_by_site_sample_summary.csv](https://github.com/Lake-Sunapee-Protective-Association/LMP/blob/main/master%20files/parameter_by_site_sample_summary.csv)

This file contains a brief summary of when each parameter was measured at any given site.

| 	Column Name		|	Column Description															|	Units	|	Format 	|
| 	----			|	----																		|	----	|	----	|
station	parameter	first_sample	last_sample	n_obs	site_type	sub_site_type

Additional information for the parameter column:
| 	Column Name		|	Column Description															|	Units	|	Format 	|
| 	----			|	----																		|	----	|	----	|

alk_mglCaCO3
chla_ugl
cl_mgl
conc_H_molpL
cond_uScm
secchidepth_m
TP_mgl
turb_NTU
DO_mgl
DO_pctsat
temp_C


