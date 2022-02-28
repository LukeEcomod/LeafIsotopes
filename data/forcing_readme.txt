Readme for forcing.csv

Contains:

-	Half-hourly environmental data and ET from SMEAR II Hyytiälä forest run
	by Institute of Atmospheric and Earth System Research (INAR) of University of Helsinki.
	For data use permissions see: https://smear.avaa.csc.fi/terms-of-use
	To obtain data from SMEAR -stations visit https://smear.avaa.csc.fi/download

-	d13Ca from Pallas-Sammaltunturi GAW-station (White et al. 2015)
	Contact sylvia.michel@colorado.edu for data use permissions
	Original weekly data interpolated linearly

-	d18Ov from IsoGSM (Yoshimure et al., 2008)
	Converted to VSMOW and interpolated linearly from original 6-hourly data
	For data use permissions see: http://isotope.iis.u-tokyo.ac.jp/~kei/tmp/isogsm2/Readme

-	d18O_prec analyzed from study site monthly precipitation
	Contact katja.rinne-garmston@luke.fi for data use permissions

Columns:

datetime: datetime [UTC + 2.0]
Tair: air temperature [degC]
CO2: mole fraction of atmospheric CO2  [umol mol-1]
H2O: mole fraciton of atmospheric water vapor [mol mol-1]
Prec: precipitation [mm d-1]
P: ambient pressure [Pa]
Ws: volumetric soil moisture in at 5 cm depth in mineral soil (m3 m-3)
Zen: solar zenith angle [rad]
Par: photosynthetically active radiation (umol m-2 s-1)
d13Ca: carbon isotope composition of atmospheric CO2 (permil, VPDB)
d18Ov: oxygen isotope composition of atmospheric CO2 (permil, VSMOW)
d18O_prec: oxygen isotope composition of precipitation (permil, VSMOW)
ET: evapotranspiration [mm d-1]

References:

White, J.W.C., B.H. Vaughn, and S.E. Michel (2015), University of Colorado, Institute of Arctic and Alpine Research (INSTAAR), Stable Isotopic Composition of Atmospheric Carbon Dioxide (13C and 18O) from the NOAA ESRL Carbon Cycle Cooperative Global Air Sampling Network, 1990-2014, Version: 20210204, Path: ftp://aftp.cmdl.noaa.gov/data/trace_gases/co2c13/flask/.

Yoshimura, K., Kanamitsu, M., Noone, D. and Oki, T. (2008) Historical isotope simulation using Reanalysis atmospheric data. Journal of Geophysical Research: Atmospheres, 113.
