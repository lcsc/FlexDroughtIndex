# FlexDroughtIndex

FlexDroughtIndex is a package of functions in R that calculates different drought indices by means of a flexible procedure described in: 
- Sergio M. Vicente-Serrano, Fergus Reig, Santiago Begueréa, Ahmed El-Kenawy, Fernando Domínguez-Castro, Magí Franquesa, Luis Gimeno-Sotelo, María Adell, Amar Halifa, Miguel Andres-Martin, Cesar Azorín-Molina: An optimal flexible approach to calculate standardized drought indices. Under review

It calculates the following drought indices:
- SPI: Standardized Precipitation Index
- SPEI:  Standardized Precipitation Evapotranspiration Index
- SEDI: Standardized Evapotranspiration Deficit Index
- SEDI: Standardized Evapotranspiration Deficit Index
  
Nevertheless, actually serves to standardize any data series that has different characteristics. For example, SPI would normalize any series truncated to 0 at the bottom, SEDI truncated at the top, and SPEI any series ranging from -inf to +inf. The functions have a fundamental advantage over previous approaches as the procedure uses the distribution that gives us a better normalized resulting series for each monthly (and scale) series according to the Shapiro-Wilks test of normality.

## EXAMPLE OF CALCULATION

- SPI(series, scale, fr, ref.start = NULL, ref.end = NULL)
  - scale is the time scale of the index
  - fr is the frequency of the data. E.g., monthly series have a frequency of 12. 
  - ref.start and ref.end define the reference period for calculations.
