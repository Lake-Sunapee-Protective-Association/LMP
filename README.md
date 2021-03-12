# Lake Sunapee Protective Association's Longterm Monitoring Program data repository

This repository contains the raw [Lake Sunapee Protective Association (LSPA)](https://www.lakesunapee.org/) Longterm Monitoring Program (LMP) data and collation scripts for all data collected by the LSPA in the watershed 1986-2020.

These data use the ODC Open Database License v1.0, see *ODC license.txt* for details on use and reuse. In addition, we invite you to contact the data owner (the LSPA) with any questions or concerns. Please cite data using the Zenodo DOI associated with this repository 10.5281/zenodo.XXX.

Suggested citation for the MASTER FILES:
**xxx; doi:xxx**

#### Contacts: 

code and repository questions: steeleb@caryinstitute.org, weathersk@caryinstitute.org

data questions: lspa@lakesunapee.org

This repository is maintained by B. Steele of the Cary Institute of Ecosystem Studies (steeleb@caryinstitute.org). 

# Lake Sunapee Watershed and Sampling sites

Below is a map of Lake Sunapee, the lake's watershed, and the sampling sites referenced in these data:

![Lake Sunapee Watershed and Sampling sites](https://github.com/Lake-Sunapee-Protective-Association/LMP/blob/main/master%20files/LMP%20sampling%20map.jpg)

# Data organization

## [raw data files](https://github.com/Lake-Sunapee-Protective-Association/LMP/tree/main/raw%20data%20files)

This folder contains the original collated data files from the LSPA that are organized for integration into the NH Volunteer Lake Assessment Program (VLAP) database. Files have come directly from the LSPA and may contain QAQC errors. These data are stored here for the purposes of transparency and processing, please use the files in the 'MASTER FILES' folder instead of these. For questions about the QAQC process from the raw files, please contact Bethel Steele (steeleb@caryinstitute.org) for information. All QAQC of these files are performed in the R scripts found in the 'collation code' folder.

These data were obtained by request from the LSPA and are archived at the New Hampshire Department of Environmental Services (NHDES) Environmental Monitoring Database (data are also available upon request from the NHDES).

## [collation code](https://github.com/Lake-Sunapee-Protective-Association/LMP/tree/main/collation%20code)

This folder contains the code used to collate and QAQC the data (where needed) to recode obviously errant data. Data from multiple files in the 'raw data files' folder are collated in this script.

## [master files](https://github.com/Lake-Sunapee-Protective-Association/LMP/tree/main/master%20files)

These are the harmonized, collated data that should be used by other researchers. Proper citation is required.


# Data status

The monitoring and data acquisition is ongoing, this repostitory will be updated approximately annually with additional data.
