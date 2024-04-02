# Author: Sergio M. Vicente-Serrano  <http://lcsc.csic.es>; Laboratorio de Climatología y Servicios Climáticos, IPE, CSIC 
# Fergus Reig Gracia <http://lcsc.csic.es>; Laboratorio de Climatología y Servicios Climáticos, IPE, CSIC 
# Santiago Beguería Portugués <http://lcsc.csic.es>; Laboratorio de Climatología y Servicios Climáticos, EEAD, CSIC 
# Version: 1.0

# This program is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation, either version 3 of the License, or any later version.

# This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more details.

# You should have received a copy of the GNU General Public License along with this program. If not, see <http://www.gnu.org/licenses/> <http://www.gnu.org/licenses/gpl.txt/>.

################################################################

# flex_drought_index.R is a package of functions in R that calculates different drought indices by means of a flexible procedure described in: 
###Sergio M. Vicente-Serrano, Fergus Reig, Santiago Begueréa, Ahmed El-Kenawy, Fernando Domínguez-Castro, Magí Franquesa, Luis Gimeno-Sotelo, María Adell, Amar Halifa, Miguel Andres-Martin, Cesar Azorín-Molina: An optimal flexible approach to calculate standardized drought indices. Under review

# It calculates the following drought indices:
# SPI: Standardized Precipitation Index
# SPEI:  Standardized Precipitation Evapotranspiration Index
# SEDI: Standardized Evapotranspiration Deficit Index
# Nevertheless, actually serves to standardize any data series that has different characteristics. For example, SPI would normalize any series truncated to 0 at the bottom, SEDI truncated at the top, and SPEI any series ranging from -inf to +inf.
# The functions have a fundamental advantage over previous approaches as the procedure uses the distribution that gives us a better normalized resulting series for each monthly (and scale) series according to the Shapiro-Wilks test of normality

#### EXAMPLE OF CALCULATION: SPI(series,scale,fr,ref.start = NULL,ref.end = NULL)
#scale is the time scale of the index
#fr is the frequency of the data. E.g., monthly series have a frequency of 12. 
#ref.start and ref.end define the reference period for calculations.

library(lmom)

#' Calculate shapiro test
#'
#' @param datos_month 
#'
#' @return shapiro test
#' @export
calc_sha = function(datos_month){
    sha <- tryCatch({
        return(shapiro.test(datos_month)$p.value)
    }, error = function(cond) {
        return(0)
    })  
    return(sha)
}

#' Apply function to the given serie (calculate the SPI, SPEI or SEDI)
#'
#' @param function spi_month_data or spei_month_data or sedi_month_data
#' @param serie vector of data
#' @param serie_par vector of reference data
#' @param scale is the time scale of the index
#' @param fr is the frequency of the data. E.g., monthly series have a frequency of 12
#'
#' @return Index serie and parameters
index_data <- function(function_month_data, serie, serie_par, scale, fr){

    name_functions <- c("exp", "gam", "gev", "glo", "gpa", "gno", "ln3", "nor", "pe3", "wei")
    name_fr = paste0("X", seq(1:fr))
    statistics <- array(NA, dim = c(length(name_functions), length(name_fr)), dimnames = list(name_functions, name_fr))

    if (scale>1) {
        serie[scale:length(serie)] <- rowSums(embed(serie, scale), na.rm = FALSE)
        serie[1:(scale-1)] <- NA

        serie_par[scale:length(serie_par)] <- rowSums(embed(serie_par, scale), na.rm = FALSE)
        serie_par[1:(scale-1)] <- NA
    }

    imonth <- 1
    for(imonth in c(1:fr)){
        month <- seq(imonth, length(serie), by = fr)
        month_par <- seq(imonth, length(serie_par), by = fr)
        serie_month <- serie[month]
        serie_month_par <- serie_par[month_par]

        ssi_month <- list()

        ssi_month[["exp"]] <- function_month_data(serie_month, serie_month_par, fun1 = pelexp, fun2 = cdfexp)
        ssi_month[["gam"]] <- function_month_data(serie_month, serie_month_par, fun1 = pelgam, fun2 = cdfgam)
        ssi_month[["gev"]] <- function_month_data(serie_month, serie_month_par, fun1 = pelgev, fun2 = cdfgev)
        ssi_month[["glo"]] <- function_month_data(serie_month, serie_month_par, fun1 = pelglo, fun2 = cdfglo)
        ssi_month[["gpa"]] <- function_month_data(serie_month, serie_month_par, fun1 = pelgpa, fun2 = cdfgpa)
        ssi_month[["gno"]] <- function_month_data(serie_month, serie_month_par, fun1 = pelgno, fun2 = cdfgno)
        ssi_month[["ln3"]] <- function_month_data(serie_month, serie_month_par, fun1 = pelln3, fun2 = cdfln3)
        ssi_month[["nor"]] <- function_month_data(serie_month, serie_month_par, fun1 = pelnor, fun2 = cdfnor)
        ssi_month[["pe3"]] <- function_month_data(serie_month, serie_month_par, fun1 = pelpe3, fun2 = cdfpe3)
        ssi_month[["wei"]] <- function_month_data(serie_month, serie_month_par, fun1 = pelwei, fun2 = cdfwei)

        sha_month <- array(NA, length(ssi_month), dimnames = list(names(ssi_month)))
        sha_month["exp"] <- calc_sha(ssi_month[["exp"]])
        sha_month["gam"] <- calc_sha(ssi_month[["gam"]])
        sha_month["gev"] <- calc_sha(ssi_month[["gev"]])
        sha_month["glo"] <- calc_sha(ssi_month[["glo"]])
        sha_month["gpa"] <- calc_sha(ssi_month[["gpa"]])
        sha_month["gno"] <- calc_sha(ssi_month[["gno"]])
        sha_month["ln3"] <- calc_sha(ssi_month[["ln3"]])
        sha_month["nor"] <- calc_sha(ssi_month[["nor"]])
        sha_month["pe3"] <- calc_sha(ssi_month[["pe3"]])
        sha_month["wei"] <- calc_sha(ssi_month[["wei"]])

        statistics[, imonth] <- sha_month

        max_month <- names(sha_month)[sha_month == max(sha_month, na.rm = FALSE)][1]
        serie[month] <- ssi_month[[max_month]]
    }
    return(list(serie = serie, statistics = statistics))
}

#' SPI: Standardized Precipitation Index
#'
#' @param serie vector of data
#' @param scale is the time scale of the index 
#' @param fr is the frequency of the data. E.g., monthly series have a frequency of 12
#' @param ref.start define the reference period for calculations
#' @param ref.end define the reference period for calculations
#'
#' @return SPI serie and parameters
#' @export
SPI <- function(serie, scale, fr, ref.start = NULL, ref.end = NULL){

    #' SPI for one period
    #'
    #' @param serie_month vector of data
    #' @param serie_month_par vector of reference data
    #' @param fun1 Probability distribution function
    #' @param fun2 Cumulative distribution function
    #'
    #' @return SPI serie
    spi_month_data <- function(serie_month, serie_month_par, fun1, fun2){
        length1 <- length(serie_month)

        no_cero_month <- which(serie_month > 0)
        cero_month <- which(serie_month == 0)
        datos_month <- serie_month[no_cero_month]

        datos_month_par <- serie_month_par[serie_month_par > 0]

        serie_month <- tryCatch({
            l_mom <- samlmu(datos_month_par, nmom = 4, sort.data = TRUE, ratios = TRUE, trim = 0)
            par_data <- fun1(l_mom)
            cdf_data <- (length(cero_month)/(length1+1))+(1-(length(cero_month)/(length1+1)))*fun2(datos_month, para = par_data)
            norm_data <- qnorm(cdf_data)
            serie_month[no_cero_month] <- norm_data
            serie_month[cero_month] <- qnorm(((length(cero_month)+1)/(2*(length1+1))))
            serie_month[serie_month == -Inf] <- -2.88
            serie_month[serie_month == Inf] <- 2.88
            return(serie_month)
        }, error = function(cond) {
            serie_month[1:length(serie_month)] <- NA
            return(serie_month)
        })

        return(serie_month)
    }

    serie_par <- as.numeric(window(serie, ref.start, ref.end))
    spi <- suppressWarnings(index_data(function_month_data = spi_month_data, serie = serie, serie_par = serie_par, scale = scale, fr = fr))

    return(spi)

}


#' SPEI: Standardized Precipitation Evapotranspiration Index
#'
#' @param serie vector of data
#' @param scale is the time scale of the index 
#' @param fr is the frequency of the data. E.g., monthly series have a frequency of 12
#' @param ref.start define the reference period for calculations
#' @param ref.end define the reference period for calculations
#'
#' @return SPEI serie and parameters
#' @export
SPEI <- function(serie, scale, fr, ref.start = NULL, ref.end = NULL){

    #' SPEI for one period
    #'
    #' @param serie_month vector of data
    #' @param serie_month_par vector of reference data
    #' @param fun1 Probability distribution function
    #' @param fun2 Cumulative distribution function
    #'
    #' @return SPEI serie
    spei_month_data <- function(serie_month, serie_month_par, fun1, fun2){
        serie_month <- tryCatch({
            l_mom <- samlmu(serie_month_par, nmom = 4, sort.data = TRUE, ratios = TRUE, trim = 0)
            par_data <- fun1(l_mom)
            cdf_data <- fun2(serie_month, para = par_data)
            norm_data <- qnorm(cdf_data)
            serie_month <- norm_data
            serie_month[serie_month == -Inf] <- -2.88
            serie_month[serie_month == Inf] <- 2.88
            return(serie_month)
        }, error = function(cond) {
            serie_month[1:length(serie_month)] <- NA
            return(serie_month)
        })

        return(serie_month)
    }

    serie_par <- as.numeric(window(serie, ref.start, ref.end))
    spei <- suppressWarnings(index_data(function_month_data = spei_month_data, serie = serie, serie_par = serie_par, scale = scale, fr = fr))

    return(spei)
}

#' SEDI: Standardized Evapotranspiration Deficit Index
#'
#' @param serie vector of data
#' @param scale is the time scale of the index
#' @param fr is the frequency of the data. E.g., monthly series have a frequency of 12
#' @param ref.start define the reference period for calculations
#' @param ref.end define the reference period for calculations
#'
#' @return SEDI serie and parameters
#' @export
SEDI <- function(serie, scale, fr, ref.start = NULL, ref.end = NULL){

    #' SEDI for one period
    #'
    #' @param serie_month vector of data
    #' @param serie_month_par vector of reference data
    #' @param fun1 Probability distribution function
    #' @param fun2 Cumulative distribution function
    #'
    #' @return SEDI serie
    #' @export
    sedi_month_data <- function(serie_month, serie_month_par, fun1, fun2){
        length1 <- length(serie_month)
        no_cero_month <- which(serie_month < 0)
        cero_month <- which(serie_month == 0)
        datos_month <- serie_month[no_cero_month]

        datos_month_par <- serie_month_par[serie_month_par < 0]

        serie_month <- tryCatch({
            l_mom <- samlmu(datos_month_par, nmom = 4, sort.data = TRUE, ratios = TRUE, trim = 0)
            par_data <- fun1(l_mom)
            cdf_data <- (1-(length(cero_month)/(length1+1)))*fun2(datos_month, para = par_data)
            norm_data <- qnorm(cdf_data)
            serie_month[no_cero_month] <- norm_data
            serie_month[cero_month] <- qnorm(1-((length(cero_month)+1)/(2*(length1+1))))
            serie_month[serie_month == -Inf] <- -2.88
            serie_month[serie_month == Inf] <- 2.88
            return(serie_month)
        }, error = function(cond) {
            serie_month[1:length(serie_month)] <- NA
            return(serie_month)
        })

        return(serie_month)
    }

    serie_par <- as.numeric(window(serie, ref.start, ref.end))
    sedi <- suppressWarnings(index_data(function_month_data = sedi_month_data, serie = serie, serie_par = serie_par, scale = scale, fr = fr))
    return(sedi)
}



