# FlexDroughtIndex

FlexDroughtIndex is a package of functions in R that calculates different drought indices by means of a flexible procedure described in: 
- Sergio M. Vicente-Serrano, Fergus Reig, Santiago Begueréa, Ahmed El-Kenawy, Fernando Domínguez-Castro, Magí Franquesa, Luis Gimeno-Sotelo, María Adell, Amar Halifa, Iván Noguera, Miguel Andres-Martin, Cesar Azorín-Molina: An optimal flexible approach to calculate standardized drought indices. Under review

It calculates the following drought indices:
- SPI: Standardized Precipitation Index
- SPEI:  Standardized Precipitation Evapotranspiration Index
- SEDI: Standardized Evapotranspiration Deficit Index
  
Nevertheless, actually serves to standardize any data series that has different characteristics. For example, SPI would normalize any series truncated to 0 at the bottom, SEDI truncated at the top, and SPEI any series ranging from -inf to +inf. The functions have a fundamental advantage over previous approaches as the procedure uses the distribution that gives us a better normalized resulting series for each monthly (and scale) series according to the Shapiro-Wilks test of normality.

## EXAMPLE OF CALCULATION

```r
source("flex_drought_index.r")

precip_aed <- read.table("data/precip_aed.csv",sep=";",header=T) #Monthly precipitation and atmospheric demand in the coordinate 0.25E 42.25N from the Climatic Research Unit dataset from 1901 to 2021
et_aed <- read.table("data/et_aed.csv",sep=";",header=T) #Monthly evaporation and atmospheric demand in the coordinate 32.375E 42.375N from the GLEAM dataset from 1980 to 2020 

#Example of calculation of 3-month SPI
spi3 <- SPI(precip_aed$precip,fr=12,scale=3)#fr refers to the frequiency of the data, in this case 12 per year.

spi3$serie # provides the SPI data
spi3$statistics # provide the AIC values for the different six distributions used for the different monthly series

#Calculation of the 12-month SPEI
spei12 <- SPEI(precip_aed$precip-precip_aed$atmospheric_demand,fr=12,scale=12)

#Calculation of the 1-month SEDI
sedi1 <- SEDI(et_aed$Evaporation-et_aed$Atmospheric_Demand,fr=12,scale=1)

#Calculation of 3-month SPI
#scale is the time scale of the index
#fr is the frequency of the data. E.g., monthly series have a frequency of 12. 
#ref.start and ref.end define the reference period for calculations.
spi3 <- SPI(precip_aed$precip,fr=12,scale=3,ref.start=1,ref.end=48)
```
