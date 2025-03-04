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

#' Cramer-von Mises (W²) Test Statistic
#'
#' @param probs_teoricas A numeric vector of theoretical probabilities (values between 0 and 1).
#'
#' @return A numeric value representing the Cramer-von Mises test statistic (W²).
#' @export
cvm_con_probs <- function(probs_teoricas) {
  probs_ordenadas=sort(probs_teoricas)
  n <- length(probs_ordenadas)  
  i <- 1:n  
  W2 <- sum((probs_ordenadas - (2 * i - 1) / (2 * n))^2) + 1 / (12 * n)
  return(W2)
}


#' Anderson-Darling (A²) Test Statistic
#'
#' @param probs_teoricas A numeric vector of theoretical probabilities (values between 0 and 1).
#'
#' @return A numeric value representing the Anderson-Darling test statistic (A²).
#' @export
ad_con_probs <- function(probs_teoricas) {
  probs_ordenadas=sort(probs_teoricas)
  n <- length(probs_ordenadas)  
  i <- 1:n  
  A2 <- -n - sum((2 * i - 1) / n * (log(probs_ordenadas) + log(1 - rev(probs_ordenadas))))
  return(A2)
}

#' Calculates the empirical cumulative distribution function (ECDF) values
#' for a given numeric series based on a specified plotting position formula.
#'
#' @param serie vector of data
#' @param a constant used in the plotting position formula
#'
#' @return cumulative probabilities for each element in the original series.
#' @export
plotting_position <- function(serie, a){
	longitud <- length(serie)
	empiricas <- rep(NA, longitud)
	ordenadas <- sort(serie, decreasing=F)
	emp <- empiricas

	for (i in 1:longitud){ 
		emp[i] <- (i - a) / (longitud + 1 - 2 * a)
	}
	for (i in 1:longitud){ 
		orden <- min(which(serie[i] == ordenadas))
		empiricas[i] <- emp[orden]
	} 
	return(empiricas)
}

#' Calculates the weighted distance (WD) between two cumulative distribution functions (CDFs).
#' It is typically used to compare modeled CDF values with empirical CDF values.
#'
#' @param modeladas A numeric vector representing the modeled CDF values 
#' @param empiricas A numeric vector representing the empirical CDF value
#'
#' @return Weighted distance between the CDFs
#' @export
weight_dist <- function(modeladas, empiricas){
		pesos <- 1 / empiricas
		weight_d <- sqrt(sum((pesos * (empiricas - modeladas)) ^ 2))
		return(weight_d)
}

#' Select the best distribution
#'
#' @param serie vector of data
#' @param method 1 = WD (Weighted Distance), 2 = Cramer-von Mises or 3 = Anderson-Darling
#'
#' @return Scoring of distributions and selected distribution
#' @export
selec_distrib_dist <- function(serie, method = 1){
	lmom <- samlmu(serie, nmom = 4, sort.data = TRUE, ratios = TRUE, trim = 0)
	param_gev <- try(pelgev(lmom), silent = TRUE)
	model_gev <- if(is(param_gev, "try-error")){rep(NA, length(serie))} else {cdfgev(serie, para = param_gev)}
	param_glo <- try(pelglo(lmom), silent = TRUE)
	model_glo <- if(is(param_glo, "try-error")){rep(NA, length(serie))} else {cdfglo(serie, para = param_glo)}
	param_gpa <- try(pelgpa(lmom), silent = TRUE)
	model_gpa <- if(is(param_gpa, "try-error")){rep(NA, length(serie))} else {cdfgpa(serie, para = param_gpa)}
	param_ln3 <- try(pelln3(lmom), silent = TRUE)
	model_ln3 <- if(is(param_ln3, "try-error")){rep(NA, length(serie))} else {cdfln3(serie, para = param_ln3)}
	param_pe3 <- try(pelpe3(lmom), silent = TRUE)
	model_pe3 <- if(is(param_pe3, "try-error")){rep(NA, length(serie))} else {cdfpe3(serie, para = param_pe3)}
	param_wei <- try(pelwei(lmom), silent = TRUE)
	model_wei <- if(is(param_wei, "try-error")){rep(NA, length(serie))} else {cdfwei(serie, para = param_wei)}
	empiricas <- plotting_position(serie, 0)

    if(method == 1){
        dist_gev <- weight_dist(model_gev, empiricas)
        dist_glo <- weight_dist(model_glo, empiricas)
        dist_gpa <- weight_dist(model_gpa, empiricas)
        dist_ln3 <- weight_dist(model_ln3, empiricas)
        dist_pe3 <- weight_dist(model_pe3, empiricas)
        dist_wei <- weight_dist(model_wei, empiricas)
        dist <- c(dist_gev, dist_glo, dist_gpa, dist_ln3, dist_pe3, dist_wei)
        min_dist <- min(dist, na.rm = TRUE)
        minimum <- min(which(min_dist == dist))
    }else if (method == 2){
        dist_gev <- cvm_con_probs(model_gev)
        dist_glo <- cvm_con_probs(model_glo)
        dist_gpa <- cvm_con_probs(model_gpa)
        dist_ln3 <- cvm_con_probs(model_ln3)
        dist_pe3 <- cvm_con_probs(model_pe3)
        dist_wei <- cvm_con_probs(model_wei)
        dist <- c(dist_gev, dist_glo, dist_gpa, dist_ln3, dist_pe3, dist_wei)
        min_dist <- max(dist, na.rm = TRUE)
        minimum <- max(which(min_dist == dist))
    }else if (method == 3){
        dist_gev <- ad_con_probs(model_gev)
        dist_glo <- ad_con_probs(model_glo)
        dist_gpa <- ad_con_probs(model_gpa)
        dist_ln3 <- ad_con_probs(model_ln3)
        dist_pe3 <- ad_con_probs(model_pe3)
        dist_wei <- ad_con_probs(model_wei)
        dist <- c(dist_gev, dist_glo, dist_gpa, dist_ln3, dist_pe3, dist_wei)
        min_dist <- max(dist, na.rm = TRUE)
        minimum <- max(which(min_dist == dist))
    }else{
        stop("Method not found")
    }

	distrib <- c("gev", "glo", "gpa", "ln3", "pe3", "wei")
	selected <- distrib[minimum]
	return(list(dist = dist, selected = selected))
}

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
#' @param method 1 = WD (Weighted Distance), 2 = Cramer-von Mises or 3 = Anderson-Darling
#'
#' @return Index serie and parameters
index_data <- function(function_month_data, serie, serie_par, scale, fr, method = 1){

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

        # sha_month <- array(NA, length(ssi_month), dimnames = list(names(ssi_month)))
        # sha_month["gev"] <- calc_sha(ssi_month[["gev"]])
        # sha_month["glo"] <- calc_sha(ssi_month[["glo"]])
        # sha_month["gpa"] <- calc_sha(ssi_month[["gpa"]])
        # sha_month["ln3"] <- calc_sha(ssi_month[["ln3"]])
        # sha_month["pe3"] <- calc_sha(ssi_month[["pe3"]])
        # sha_month["wei"] <- calc_sha(ssi_month[["wei"]])
        # statistics[, imonth] <- sha_month
        # max_month <- names(sha_month)[sha_month == max(sha_month, na.rm = FALSE)][1]

        selec_distrib <- selec_distrib_dist(serie_month_par[!is.na(serie_month_par)], method = method)
        if(!is.na(selec_distrib$selected)){
            max_month <- selec_distrib$selected
            statistics[, imonth] <- selec_distrib$dist
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
SPI <- function(serie, scale, fr, ref.start = NULL, ref.end = NULL, method = 1){

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

    serie_par <- suppressWarnings(window(serie, ref.start, ref.end, frequency = fr))
    spi <- suppressWarnings(index_data(function_month_data = spi_month_data, serie = serie, serie_par = serie_par, scale = scale, fr = fr, method = method))

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
SPEI <- function(serie, scale, fr, ref.start = NULL, ref.end = NULL, method = 1){

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

    serie_par <- suppressWarnings(window(serie, ref.start, ref.end, frequency = fr))
    spei <- suppressWarnings(index_data(function_month_data = spei_month_data, serie = serie, serie_par = serie_par, scale = scale, fr = fr, method = method))

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
SEDI <- function(serie, scale, fr, ref.start = NULL, ref.end = NULL, method = 1){

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

    serie_par <- suppressWarnings(window(serie, ref.start, ref.end, frequency = fr))
    sedi <- suppressWarnings(index_data(function_month_data = sedi_month_data, serie = serie, serie_par = serie_par, scale = scale, fr = fr, method = method))
    return(sedi)
}



