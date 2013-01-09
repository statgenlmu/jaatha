# --------------------------------------------------------------
# ModelFitting.R
# Functions for the machine learning part of Jaatha. 
# 
# Authors:  Lisha Naduvilezhath & Paul R. Staab
# Date:     2012-10-05
# Licence:  GPLv3 or later
# --------------------------------------------------------------

## Function to fit a glm for each summary statistic with the generated
## parameter combinations in the block. The first 'bObject@nPar'
## columns are assumed to be the columns for the parameters, the
## following 'nTotalSumstat' colums contain the results for the summary
## statistic.
glmFitting <- function(bObject,nTotalSumstat,weighting){ 
			#cat("Fitting model ... \n")
			##+3 bc intercept, convergence, sumOfSumstat
			modFeld <- array(-1,dim=c(nTotalSumstat,(bObject@nPar+3)))
			dat <- data.frame(bObject@parNsumstat)
			## Create a formula for a model with nPar variables:
			## they should explain the observed SS
			explanatory <- paste("X", 1:bObject@nPar, sep="",collapse= "+")
			##=X1+..+X(nPar)  these are parameter columns
			response <- paste("X", (1:nTotalSumstat)+bObject@nPar, sep="")
			##=X(nPar+1),...,X(nPar+nTotalSumstat) these are 'nTotalSumstat' ss columns
			
			## suppresses warnings; which ocurr whenever sumstat is 0 "suppressWarnings()"		   
			for(s in 1:nTotalSumstat) {
				## if glm function did not converge, set coefficients & convergence to 0
				tryCatch({
					mod <- suppressWarnings(glm(as.formula(paste(response[s]," ~ ", explanatory)),
								data=dat, family=poisson, weights=weighting,
                                control = list(maxit = 200)))}, 
					error=function(e) {
						print(list("Caught error of GLM function!"))
						cat("sumstat",s," sum=",sum(dat[,s+bObject@nPar]),"\n")
						mod <- list(coef=rep(0,bObject@nPar+1),conv=0) })
				#cat(s,coef(summary(mod))[,4],mod$aic,"\n",file="p-values.txt",append=T,sep="\t")
			
				modFeld[s,] <- c(mod$coef,mod$conv,sum(dat[,s+bObject@nPar])) 
				if (!mod$conv){
					cat("WARNING: sumstat",s," did not converge, sum = ",
							sum(dat[,s+bObject@nPar]),"\n")
				} else{}
			}  
			return (modFeld)
		}


## Function to estimate the best parameters within the block and for
## those the highest composite likelihood for the specific model
## (specified in Simulator.simulateWithinBlock()). Input:
## Block-'object' to be optimized within, 'jObject' for getting the
## parameter ranges for the search, 'modFeld' which holds the
## coefficients, convergence and sum of the fitted glm for each
## summary statistic, 'ssData' the summary statistic of the data for
## which parameters to search. 'boarder' determines how much of the
## around the boundaries of the blocks should not be used for the
## optimization procedure. With 'opt' the optimization strategie can
## be chosen: 1=fast, 2=vote,3=medium.
estimate <- function(bObject,jObject,modFeld,ssData,boarder=0.25,opt=3){
			dimSize <- bObject@upperBound-bObject@lowerBound
			mitte <- dimSize/2 + bObject@lowerBound
			## all dimensions are assumed to be of same length
			#cat("Searching for ML and MLest ... \n")
			if (jObject@externalTheta){ 
				##############
				## FAST METHOD
				if(opt < 3) {
					if (opt==1) cat("Using fast optimization. \n") else{}
					score <-  .calcLikelihoodWithModelfeld(param=mitte,
							modelCoefficients=modFeld,
							observedSS=ssData,
							jObject=jObject,bObject=bObject)  
				}
				###################
				## MEDIUM FAST METHOD
				else if(opt==3) {
					##optimization function
					optfunk <- function(par) {                                 
						if(min(par)<0 || max(par)>1){ return(1e11)} else{}
						score <- .calcLikelihoodWithModelfeld(param=par,
								modelCoefficients=modFeld,
								observedSS=ssData, jObject=jObject,bObject=bObject)  
						#print(par)
						#cat("score:",score,"\n")
						return(-score)  ## '-' enables Maximization
					}                
					##boarder to leave free in the block for optimization
					##(depends on boundaries)
					puffer <- boarder*dimSize   
					##calculate limits for the optimization procedure
					untere <- sapply(1:jObject@nPar,
							function(p) min(bObject@lowerBound[p]+puffer[p], 1))
					obere <-  sapply(1:jObject@nPar,
							function(p) max(bObject@upperBound[p]-puffer[p], 0))
					#cat(boarder,"unt:",untere,"obe:",obere,"puffer",puffer,"\n")
					
					##describes 'boarder'% of values that will be excluded on
					##either side of the block in optimization
					OOO <- optim(mitte, optfunk, lower=untere, upper=obere,
							method="L-BFGS-B")
					##the optimal parameters are contained in the vector mitte
					mitte <- OOO$par  
					score <- -OOO$value  
				}
				####################
				## WEIGHTED VOTE METHOD
				theta <- sum(ssData)/sum(
							exp(modFeld[1:jObject@nTotalSumstat,1:(bObject@nPar+1)]
							%*%c(1,mitte))) * (5*bObject@nLoci)
				if(opt==2) { 
					cat("Using vote optimization. \n")
					parEst <- rep(0,bObject@nPar)
					wsum <- 0
					logThetaEst <- 0
					slL <- sum(score)
					score <- score - max(score,na.rm=T)
					if(!is.na(score)) {
						wsum <- wsum + exp(score)
						parEst <- parEst + exp(score)* mitte
						logThetaEst <- logThetaEst + exp(score)*log(theta)  
					} else{}
					#cat(" opt2: best par:",parEst,exp(LogThetaEst/wsum),"\n")
					return (list(est= parEst,
									theta= exp(logThetaEst/wsum), 
									score= slL))
				}
				else{
					#cat("-> best score:",round(score,2),"param:",
					#    round(.calcAbsParamValue(jObject@dm,mitte),3),"\n")
					return (list(est = mitte,
						     theta = theta/jObject@nLoci,
						     score = score))            
				}
			}
			
			
			###### if not externalTheta  ##########################################
			else{
				###################
				## MEDIUM FAST METHOD
				if(opt==3) {
					#cat("Using medium optimization. \n")
					##optimization function
					optfunk <- function(par) {                                 
						if(min(par)<0 || max(par)>1){ return(1e11)} else{}
						score <- .calcLikelihoodWithModelfeld(param=par,
								modelCoefficients=modFeld,
								observedSS=ssData, jObject=jObject,
								bObject=bObject)
						return(-score)  ## '-' enables Maximization
					}
					##boarder to leave free in the block for optimization
					##(depends on boundaries)
					puffer <- boarder*dimSize   
					##calculate limits for the optimization procedure
					untere <- sapply(1:jObject@nPar,
							function(p) min(bObject@lowerBound[p] + puffer,1))
					obere <-  sapply(1:jObject@nPar,
							function(p) max(bObject@upperBound[p] - puffer,0))
					#cat(boarder,"unt:",untere,"obe:",obere,"puffer",puffer,"\n")
					
					##describes 'boarder'% of values that will be excluded
					##on either side of the block in optimization
					OOO <- optim(mitte,optfunk,lower= untere,upper= obere,
							method="L-BFGS-B")
					mitte <-OOO$par
					##the optimal parameters are contained in the vector mitte
					score <- -OOO$value  
				}
			
				#cat("-> best score:",round(score,2),
				#	"par:",round(.calcAbsParamValue(jObject@dm,mitte),3),"\n")
				return (list(est= mitte, score= score))                   
			}            
		}


.calcLikelihoodWithModelfeld <- function(param, modelCoefficients, observedSS, 
		jObject, bObject){
	##modelCoefficients: +3 bc intercept, convergence, sumOfSumstat
	score <- 0	
	if (jObject@externalTheta){ ## theta gets only calculated for the middle 	
		## middle point of current block
		point <- (bObject@upperBound-bObject@lowerBound)/2 + bObject@lowerBound
		#scale <- 5*bObject@nLoci
		thetaTotal <- sum(observedSS)/sum(exp(
					modelCoefficients[1:jObject@nTotalSumstat,1:(bObject@nPar+1)] 
					%*%	c(1,point)))   

		for (s in 1:jObject@nTotalSumstat) {     
			##if glm did not converge, take sum(SS[s]) or a small number like 0.5
			if (modelCoefficients[s,(bObject@nPar+2)]<0.5) {				  
				loglambda <- max(0.5,modelCoefficients[s,(bObject@nPar+3)]) +
							 log(thetaTotal) 
				## 350 = 10 * 5 *7 = repetitions * no of loci * theta
				##cat("theta=",theta,"\n")
			} else {                        
				## log(exp(modelC))=modelC; if scale is needed: - log(bObject@nLoci/jObject@nLoci)
				loglambda <- (modelCoefficients[s,1:(bObject@nPar+1)]%*%
						c(1,param)) + log(thetaTotal) 
			}
			score <- score + observedSS[s]*loglambda - exp(loglambda)
		}
	} 
	else {  ##if NOT externalTheta		
		point <- param
		#scale <- bObject@nLoci/jObject@nLoci
	
		for (s in 1:jObject@nTotalSumstat) {     
			##if glm did not converge, take sum(SS[s]) or a small number like 0.5
			if (modelCoefficients[s,(bObject@nPar+2)]<0.5) {				  
				loglambda <- max(0.5,modelCoefficients[s,(bObject@nPar+3)])	 
				## 350 = 10 * 5 *7 = repetitions * no of loci * theta
				##cat("theta=",theta,"\n")
			} else {                        
				## log(exp(modelC))=modelC; if scale is needed: - log(bObject@nLoci/jObject@nLoci)
				loglambda <- (modelCoefficients[s,1:(bObject@nPar+1)]%*%c(1,param))
			}
			score <- score + observedSS[s]*loglambda - exp(loglambda)
		}
	}
	##cat("score:",score,"\n")
	return(score)
}


.calcLikelihoodWithModelfeld.old <- function(param, modelCoefficients, observedSS, 
		jObject, bObject){
	##modelCoefficients: +3 bc intercept, convergence, sumOfSumstat
	score <- 0	
	if (jObject@externalTheta){ ## theta gets only calculated for the middle 	
		## middle point of current block
		point <- (bObject@upperBound-bObject@lowerBound)/2 + bObject@lowerBound
		#scale <- 5*bObject@nLoci


	} else{  ##if NOT externalTheta		
		point <- param
		#scale <- bObject@nLoci/jObject@nLoci
	}
	thetaTotal <- sum(observedSS)/sum(exp(
					modelCoefficients[1:jObject@nTotalSumstat,1:(bObject@nPar+1)] 
					%*%	c(1,point)))   
	## likelihood is calculated by going through each SS seperately 
	for (s in 1:jObject@nTotalSumstat) {     
		##if glm did not converge, take sum(SS[s]) or a small number like 0.5
		if (modelCoefficients[s,(bObject@nPar+2)]<0.5) {				  
			loglambda <- max(0.5,modelCoefficients[s,(bObject@nPar+3)]) +
							 log(thetaTotal) 
			## 350 = 10 * 5 *7 = repetitions * no of loci * theta
			##cat("theta=",theta,"\n")
		} else {                        
			## log(exp(modelC))=modelC; if scale is needed: - log(bObject@nLoci/jObject@nLoci)
			loglambda <- (modelCoefficients[s,1:(bObject@nPar+1)]%*%
						c(1,param)) + log(thetaTotal) 
		}
		score <- score + observedSS[s]*loglambda - exp(loglambda)
	}
	##cat("score:",score,"\n") 
	return(score)
}
