
ProbabSphere <- function( Rad, dim){
    SQRT2 <- 2^(0.5)
	SQRT2DIVPI <- (2/pi)^(0.5)
    return <- 0
	m <- 0
	i <- 0
	coef <- c(0,0)
	factor <- Rad^2
	mult <- 0.0
	temp <- 0.0
	prom <- 0.0
	coef[1] <- exp(-1*factor*0.5)
	if(dim%%2==0){
        m <- dim/2
        coef[2] <- 1.0
        for(i in 1:(m-1)){
            coef[2] <- coef[2]*(2*i)
        }
        temp <- i-coef[1]
        i<-3
        mult <- factor
        while(i < dim){
            temp <- temp*(i-1)-coef[1]*mult
            mult <- mult * factor
            i <- i + 2
        }
        return <- (temp/coef[2])
	}
	else{
        #print("The dimension is not even.")
		m <- dim/2
        coef[2] <- 1.0
        #print ( "m is: ")
        #print( m )
        if(m >= 2)
        for(i in 2:m){
            coef[2] <- coef[2]*(2*i-1)
            #print("Coef[2] is: ")
            #print(coef[2])
        }
        mult <- Rad
        #print("ERF::Rad/SQRT2 is: ")
        #print(erf(Rad/SQRT2))
        temp <- .erf(Rad/SQRT2)/SQRT2DIVPI
        #print("Temp is: ")
        #print(temp)
        prom <- temp
        i <- 1
        while(i<dim){
            temp <- temp*i - coef[1]*mult
            #print("temp is now: ")
            #print(temp)
            mult <- mult * factor
            i <- (i + 2)
        }
        #print("temp is finally: ")
        #print(temp)
        #print("coef is:")
        #print(coef[2])
        #print("The probabsphere return value will be: ")
        #print((temp*SQRT2DIVPI/coef[2]))
        return <- (temp*SQRT2DIVPI/coef[2])
	}
return
}

setGeneric(".OPV",function( transcriptAbundanceValues, ...) {
	standardGeneric(".OPV")
})

setMethod(".OPV", signature(transcriptAbundanceValues="vector"), function(transcriptAbundanceValues, libraryTotals, concentrationLimit=0.0, pVarianceLimit=0.0, messagePrefix="", ...){
	#check to see if there are as many readCounts as there are libraryTotals
	if(length(transcriptAbundanceValues)!=length(libraryTotals)){
		warning(paste(messagePrefix,"The number of transcriptAbundance entries do not match the number of libraryTotals values. p-value of -1 is returned."))
	    return(-1)
	}
	else if(sum(transcriptAbundanceValues)==0){
		warning(paste(messagePrefix,"This transcript has no non-zero abundance values. p-value of -1 is returned."))
	    return(-1)
	}
	else{
		#define the concentration ratio value for each transcriptAbundance across libraries
		concentrations <- transcriptAbundanceValues/libraryTotals
		if(var(concentrations) < pVarianceLimit){
			warning(paste(messagePrefix,"The observed variance between abundance concentrations is not larger than the limit of",pVarianceLimit,"specified. p-value of -1 is returned."))
			return(-1)
		}
		else if(max(concentrations) < concentrationLimit){
			warning(paste(messagePrefix,"No transcript abundance ratio is larger than the limit of",concentrationLimit,"specified. p-value of -1 is returned."))
			return(-1)
		}
		else{
			
			#all input checking complete, calculate the concentration vector length and return a p-value
			maximumLikelihoodConcentration <- sum(transcriptAbundanceValues) / sum(libraryTotals) 
			#subtract the mean, or MLE, concentration
			centeredConcentration <- concentrations - maximumLikelihoodConcentration
			concentrationVectorDistance <- sqrt( sum(centeredConcentration^2 * libraryTotals / (maximumLikelihoodConcentration*(1-maximumLikelihoodConcentration))) )
			return(1.0 - ProbabSphere(concentrationVectorDistance,length(transcriptAbundanceValues)))
		}
	}
})




.erf <- function(x) 2 * pnorm(x * sqrt(2)) - 1

setGeneric("cltDexp",function( object, ...) {
	standardGeneric("cltDexp")
})

setMethod("cltDexp", signature(object="matrix"), function(object, abundanceLimit=0.0, varianceLimit=0.0, ... ){

	data.totals <- colSums(object)
	
	pvals <- apply(object,1,.OPV,libraryTotals=data.totals,concentrationLimit=abundanceLimit,pVarianceLimit=varianceLimit)
	
	return(pvals)
})

setMethod("cltDexp", signature(object="vector"), function(object, treatmentTotals, abundanceLimit=0.0, varianceLimit=0.0,... ){ 
 
	message("vectorFound")
	if(missing(treatmentTotals)) stop("treatmentTotals not defined: must be a numeric vector")
	return(.OPV(object,libraryTotals=treatmentTotals,concentrationLimit=abundanceLimit,pVarianceLimit=varianceLimit) )
 
} )

