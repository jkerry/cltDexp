
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
        #print("The dimension is even.")
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
        #print(paste("The probabsphere return value will be: ",(temp/coef[2])))
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



OlmanPValue <- function(readCounts, libraryTotals, concentrationLimit=0.0, pVarianceLimit=0.0){
    Sum <- c(0,0)
    i <- 0
    Rm <- 0
    Estimation <- 0.0
    Stat <- 0.0
    ratio <- c(rep(0.0,length(readCounts)))
    
    #test for equal dimension
    if(length(readCounts)==length(libraryTotals)){
        for(i in 1:length(readCounts)){
            ratio[i] <- readCounts[i]/libraryTotals[i]
            if(libraryTotals[i]>0 && ratio[i]<concentrationLimit){
                Rm <- Rm + 1
            }
        }
        #test for passing concentration limit
        if(Rm != length(readCounts)){
	    
	    if(var(ratio) >= pVarianceLimit){
            #calculate MLE P
            for(i in 1:length(readCounts)){
                Sum[1] <- Sum[1]+readCounts[i]
                Sum[2] <- Sum[2]+libraryTotals[i]
            }
            if(Sum[1]>0){
                Estimation <- Sum[1]/Sum[2]
                #calc vector distance
                for(i in 1: length(readCounts)){
                    Stat <- Stat + (ratio[i]-Estimation)^2 * libraryTotals[i]
                }
                
		Stat <- (Stat/(Estimation*(1-Estimation)))^0.5
                return(1.0-ProbabSphere(Stat, length(readCounts)))
            }
            else{
                warning("No read information found for this transcript.")
		return(-1)
            }
	    
	    }
	    else{
		message("The variance among the reads is under the minimum threshold in order to calculate a p-value")
		return(-1)
	    }
        }
        else{
            message("All read counts are under the minimum concentration threshold in order to calculate a p-value")
	    return(-1)
        }
    }
    else{
        warning("The number of read counts does not equal the number of library sizes.\n")
	return(-1)
    }
}


.erf <- function(x) 2 * pnorm(x * sqrt(2)) - 1

setGeneric("cltDexp",function( object, ...) {
	standardGeneric("cltDexp")
})

setMethod("cltDexp", signature(object="matrix"), function(object, abundanceLimit=0.0, varianceLimit=0.0, ... ){

	data.totals <- colSums(object)
	
	pvals <- apply(object,1,OlmanPValue,libraryTotals=data.totals,concentrationLimit=abundanceLimit,pVarianceLimit=varianceLimit)
	
	return(pvals)
})

setMethod("cltDexp", signature(object="vector"), function(object, treatmentTotals, abundanceLimit=0.0, varianceLimit=0.0,... ){ 
 
	message("vectorFound")
	if(missing(treatmentTotals)) stop("treatmentTotals not defined: must be a numeric vector")
	return(OlmanPValue(object,libraryTotals=treatmentTotals,concentrationLimit=abundanceLimit,pVarianceLimit=varianceLimit) )
 
} )

