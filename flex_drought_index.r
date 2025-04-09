# Author: Sergio M. Vicente-Serrano  <http://lcsc.csic.es>; Laboratorio de Climatología y Servicios Climáticos, IPE, CSIC 
# Fergus Reig Gracia <http://lcsc.csic.es>; Laboratorio de Climatología y Servicios Climáticos, IPE, CSIC 
# Santiago Beguería Portugués <http://lcsc.csic.es>; Laboratorio de Climatología y Servicios Climáticos, EEAD, CSIC 
# Version: 1.0

# This program is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation, either version 3 of the License, or any later version.

# This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more details.

# You should have received a copy of the GNU General Public License along with this program. If not, see <http://www.gnu.org/licenses/> <http://www.gnu.org/licenses/gpl.txt/>.

################################################################

# flex_drought_index.R is a package of functions in R that calculates different drought indices by means of a flexible procedure described in: 
###Sergio M. Vicente-Serrano, Fergus Reig, Santiago Beguería, Ahmed El-Kenawy, Fernando Domínguez-Castro, Magí Franquesa, Luis Gimeno-Sotelo, María Adell, Amar Halifa, Miguel Andres-Martin, Cesar Azorín-Molina: An optimal flexible approach to calculate standardized drought indices. Under review

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

logLik_cdf <- function(cdf_theor, cdf_emp) {
  eps <- 1e-10  # Evitar log(0)
  p_empirical <- pmax(diff(cdf_emp), eps)  # Probabilidad empírica en cada intervalo
  p_theoretical <- pmax(diff(cdf_theor), eps)  # Probabilidad teórica en cada intervalo
  length(cdf_emp)*sum(p_empirical * log(p_theoretical))
}
AIC_cdf <- function(cdf_theor, cdf_emp, k) {
  logL <- logLik_cdf(cdf_theor, cdf_emp)
  2 * k - 2 * logL
}


#' Select the best distribution
#'
#' @param serie vector of data
#'
#' @return Scoring of distributions and selected distribution
#' @export

selec_distrib_AIC <- function(serie){
	lmom_data <- samlmu(serie, nmom = 4, sort.data = TRUE, ratios = TRUE, trim = 0)
	fit_gev <- try(pelgev(lmom_data),silent=TRUE) 	  # Generalized Extreme Value (GEV)
	fit_gev <- if(is(fit_gev, "try-error")){fit_gev <- NA}else{fit_gev <- fit_gev}
	fit_glo <- try(pelglo(lmom_data),silent=TRUE) 	  # Generalized Logistic (GLO)
	fit_glo <- if(is(fit_glo, "try-error")){fit_glo <- NA}else{fit_glo <- fit_glo}
	fit_gpa <- try(pelgpa(lmom_data),silent=TRUE) 	  # Generalized Pareto (GPD)
	fit_gpa <- if(is(fit_gpa, "try-error")){fit_gpa <- NA}else{fit_gpa <- fit_gpa}
	fit_ln3 <- try(pelln3(lmom_data),silent=TRUE) 	  # log-Normal
	fit_ln3 <- if(is(fit_ln3, "try-error")){fit_ln3 <- NA}else{fit_ln3 <- fit_ln3}
	fit_pe3 <- try(pelpe3(lmom_data),silent=TRUE) 	   # Pearson 3
	fit_pe3 <- if(is(fit_pe3, "try-error")){fit_pe3 <- NA}else{fit_pe3 <- fit_pe3}
	fit_wei <- try(pelwei(lmom_data),silent=TRUE) 	  # Weibull
	fit_wei <- if(is(fit_wei, "try-error")){fit_wei <- NA}else{fit_wei <- fit_wei}

	x_vals <- sort(serie)
	cdf_emp <- ecdf(serie)(x_vals)

	cdf_gev <- try(cdfgev(x_vals, fit_gev),silent=TRUE)
	cdf_glo <- try(cdfglo(x_vals, fit_glo),silent=TRUE)
	cdf_gpa <- try(cdfgpa(x_vals, fit_gpa),silent=TRUE)
	cdf_ln3 <- try(cdfln3(x_vals, fit_ln3),silent=TRUE)
	cdf_pe3 <- try(cdfpe3(x_vals, fit_pe3),silent=TRUE)
	cdf_wei <- try(cdfwei(x_vals, fit_wei),silent=TRUE)
	cdf_gev <- if(is(cdf_gev, "try-error")){cdf_gev <- NA}else{cdf_gev <- cdf_gev}
	cdf_glo <- if(is(cdf_glo, "try-error")){cdf_glo <- NA}else{cdf_glo <- cdf_glo}
	cdf_gpa <- if(is(cdf_gpa, "try-error")){cdf_gpa <- NA}else{cdf_gpa <- cdf_gpa}
	cdf_ln3 <- if(is(cdf_ln3, "try-error")){cdf_ln3 <- NA}else{cdf_ln3 <- cdf_ln3}
	cdf_pe3 <- if(is(cdf_pe3, "try-error")){cdf_pe3 <- NA}else{cdf_pe3 <- cdf_pe3}
	cdf_wei <- if(is(cdf_wei, "try-error")){cdf_wei <- NA}else{cdf_wei <- cdf_wei}

	AIC_gev <- try(AIC_cdf(cdf_gev, cdf_emp, 3), silent = TRUE)
	AIC_gev <- if(is(AIC_gev, "try-error")){AIC_gev <- NA} else {AIC_gev <- AIC_gev}
	AIC_glo <- try(AIC_cdf(cdf_glo, cdf_emp, 3), silent = TRUE)
	AIC_glo <- if(is(AIC_glo, "try-error")){AIC_glo <- NA} else {AIC_glo <- AIC_glo}
	AIC_gpa <- try(AIC_cdf(cdf_gpa, cdf_emp, 3), silent = TRUE)
	AIC_gpa <- if(is(AIC_gpa, "try-error")){AIC_gpa <- NA} else {AIC_gpa <- AIC_gpa}
	AIC_ln3 <- try(AIC_cdf(cdf_ln3, cdf_emp, 3), silent = TRUE)
	AIC_ln3 <- if(is(AIC_ln3, "try-error")){AIC_ln3 <- NA} else {AIC_ln3 <- AIC_ln3}
	AIC_pe3 <- try(AIC_cdf(cdf_pe3, cdf_emp, 3), silent = TRUE)
	AIC_pe3 <- if(is(AIC_pe3, "try-error")){AIC_pe3 <- NA} else {AIC_pe3 <- AIC_pe3}
	AIC_wei <- try(AIC_cdf(cdf_wei, cdf_emp, 3), silent = TRUE)
	AIC_wei <- if(is(AIC_wei, "try-error")){AIC_wei <- NA} else {AIC_wei <- AIC_wei}

	AIC_gev <- if(length(cdf_gev) == 1){AIC_gev <- NA}else{AIC_gev <- AIC_gev}
	AIC_glo <- if(length(cdf_glo) == 1){AIC_glo <- NA}else{AIC_glo <- AIC_glo}
	AIC_gpa <- if(length(cdf_gpa) == 1){AIC_gpa <- NA}else{AIC_gpa <- AIC_gpa}
	AIC_ln3 <- if(length(cdf_ln3) == 1){AIC_ln3 <- NA}else{AIC_ln3 <- AIC_ln3}
	AIC_pe3 <- if(length(cdf_pe3) == 1){AIC_pe3 <- NA}else{AIC_pe3 <- AIC_pe3}
	AIC_wei <- if(length(cdf_wei) == 1){AIC_wei <- NA}else{AIC_wei <- AIC_wei}

	AIC <- c(AIC_gev,AIC_glo,AIC_gpa,AIC_ln3,AIC_pe3,AIC_wei)
      min_AIC <- min(AIC, na.rm = TRUE)
      minimum <- min(which(min_AIC == AIC))

	distrib <- c("gev", "glo", "gpa", "ln3", "pe3", "wei")
	selected <- distrib[minimum]
	return(list(AIC = AIC, selected = selected))
}

#' Apply function to the given serie (calculate the SPI, SPEI or SEDI)
#'
#' @param function spi_month_data or spei_month_data or sedi_month_data
#' @param serie vector of data
#' @param serie_par vector of reference data
#' @param scale is the time scale of the index
#' @param fr is the frequency of the data. E.g., monthly series have a frequency of 12
#' @param method 1 = WD (Weighted Distance), 2 = Cramer-von Mises or 3 = Anderson-Darling
#'
#' @return Index serie and parameters
index_data <- function(function_month_data, serie, serie_par, scale, fr){

    name_functions <- c("gev", "glo", "gpa", "ln3", "pe3", "wei")
    name_fr = paste0("X", seq(1:fr))
    statistics <- array(NA, dim = c(length(name_functions), length(name_fr)), dimnames = list(name_functions, name_fr))

    if (scale>1) {
        serie[scale:length(serie)] <- rowSums(embed(serie, scale), na.rm = FALSE)
        serie[1:(scale - 1)] <- NA

        serie_par[scale:length(serie_par)] <- rowSums(embed(serie_par, scale), na.rm = FALSE)
        serie_par[1:(scale - 1)] <- NA
    }

    imonth <- 1
    for(imonth in c(1:fr)){
        month <- seq(imonth, length(serie), by = fr)
        month_par <- seq(imonth, length(serie_par), by = fr)
        serie_month <- serie[month]
        serie_month_par <- serie_par[month_par]

        ssi_month <- list()

        ssi_month[["gev"]] <- function_month_data(serie_month, serie_month_par, fun1 = pelgev, fun2 = cdfgev)
        ssi_month[["glo"]] <- function_month_data(serie_month, serie_month_par, fun1 = pelglo, fun2 = cdfglo)
        ssi_month[["gpa"]] <- function_month_data(serie_month, serie_month_par, fun1 = pelgpa, fun2 = cdfgpa)
        ssi_month[["ln3"]] <- function_month_data(serie_month, serie_month_par, fun1 = pelln3, fun2 = cdfln3)
        ssi_month[["pe3"]] <- function_month_data(serie_month, serie_month_par, fun1 = pelpe3, fun2 = cdfpe3)
        ssi_month[["wei"]] <- function_month_data(serie_month, serie_month_par, fun1 = pelwei, fun2 = cdfwei)

        selec_distrib <- selec_distrib_AIC(serie_month_par[!is.na(serie_month_par)])
        if(!is.na(selec_distrib$selected)){
            max_month <- selec_distrib$selected
            statistics[, imonth] <- selec_distrib$AIC
            serie[month] <- ssi_month[[max_month]]
        }else{
            statistics[, imonth] <- NA
            serie[month] <- NA
        }
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
#' @param method 1 = WD (Weighted Distance), 2 = Cramer-von Mises or 3 = Anderson-Darling
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

        length1_par <- length(serie_month_par)
        no_cero_month_par <- which(serie_month_par > 0)
        cero_month_par <- which(serie_month_par == 0)
        datos_month_par <- serie_month_par[no_cero_month_par]

        serie_month <- tryCatch({
            l_mom <- samlmu(datos_month_par, nmom = 4, sort.data = TRUE, ratios = TRUE, trim = 0)
            par_data <- fun1(l_mom)
            cdf_data <- (length(cero_month_par)/(length1_par+1))+(1-(length(cero_month_par)/(length1_par+1)))*fun2(datos_month, para = par_data)
            norm_data <- qnorm(cdf_data)
            serie_month[no_cero_month] <- norm_data
            serie_month[cero_month] <- qnorm(((length(cero_month_par)+1)/(2*(length1_par+1))))
            serie_month[serie_month == -Inf] <- -2.88
            serie_month[serie_month == Inf] <- 2.88
            return(serie_month)
        }, error = function(cond) {
            serie_month[1:length(serie_month)] <- NA
            return(serie_month)
        })

        return(serie_month)
    }

    if(class(serie) == "numeric"){
        expected_length <- ceiling(length(serie) / fr)
        serie <- c(serie, rep(NA, fr * expected_length - length(serie)))
        serie <- ts(data = serie, start = c(1, 1), end = c(expected_length, fr), frequency = fr)
    }
    serie_par <- suppressWarnings(window(serie, ref.start, ref.end, frequency = fr))
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
#' @param method 1 = WD (Weighted Distance), 2 = Cramer-von Mises or 3 = Anderson-Darling
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
    if(class(serie) == "numeric"){
        expected_length <- ceiling(length(serie) / fr)
        serie <- c(serie, rep(NA, fr * expected_length - length(serie)))
        serie <- ts(data = serie, start = c(1, 1), end = c(expected_length, fr), frequency = fr)
    }
    serie_par <- suppressWarnings(window(serie, ref.start, ref.end, frequency = fr))
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
#' @param method 1 = WD (Weighted Distance), 2 = Cramer-von Mises or 3 = Anderson-Darling
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

        length1_par <- length(serie_month_par)
        no_cero_month_par <- which(serie_month_par < 0)
        cero_month_par <- which(serie_month_par == 0)
        datos_month_par <- serie_month_par[no_cero_month_par]

        serie_month <- tryCatch({
            l_mom <- samlmu(datos_month_par, nmom = 4, sort.data = TRUE, ratios = TRUE, trim = 0)
            par_data <- fun1(l_mom)
            cdf_data <- (1-(length(cero_month_par)/(length1_par+1)))*fun2(datos_month, para = par_data)
            norm_data <- qnorm(cdf_data)
            serie_month[no_cero_month] <- norm_data
            serie_month[cero_month] <- qnorm(1-((length(cero_month_par)+1)/(2*(length1_par+1))))
            serie_month[serie_month == -Inf] <- -2.88
            serie_month[serie_month == Inf] <- 2.88
            return(serie_month)
        }, error = function(cond) {
            serie_month[1:length(serie_month)] <- NA
            return(serie_month)
        })

        return(serie_month)
    }
    if(class(serie) == "numeric"){
        expected_length <- ceiling(length(serie) / fr)
        serie <- c(serie, rep(NA, fr * expected_length - length(serie)))
        serie <- ts(data = serie, start = c(1, 1), end = c(expected_length, fr), frequency = fr)
    }
    serie_par <- suppressWarnings(window(serie, ref.start, ref.end, frequency = fr))
    sedi <- suppressWarnings(index_data(function_month_data = sedi_month_data, serie = serie, serie_par = serie_par, scale = scale, fr = fr))
    return(sedi)
}



