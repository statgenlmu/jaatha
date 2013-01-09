# --------------------------------------------------------------
# DataProcessor.R
# Various functions for importing and processing data from external 
# sources.
# 
# Authors:  Lisha Naduvilezhath & Paul R. Staab
# Date:     2012-10-05
# Licence:  GPLv3 or later
# --------------------------------------------------------------

##Function to generate new msoutput file with given parameter values,
##calculate its jsfs and return its summary statistics. Command for ms
##needs to be put in manually (Simulator.simulate()).
#Jaatha.msSimulatedSS <- function(object,dm,sumstatFunc,nTotalSumstat,finiteS=FALSE){
#	cat("just called generic msSimulatedSS \n")}
#setMethod("Jaatha.msSimulatedSS", signature(object = "DataProcessor"),
#		function(object,dm,sumstatFunc,nTotalSumstat,finiteS=FALSE){
#
#			Jaatha.simulate(dm, nSample=object@popSampleSizes, nLoci=object@nLoci,
#					par=c(.5,.5,.5,.5,.5),finite=finiteS, fileName="../data/msoutputToFind")
#			return(sumstatFunc(dm.getJSFS(dm)))
#		})



## load the C++ code platform independently into workspace
#dyn.load(paste("../../R-eclipse/ms2jsfs_C++",.Platform$dynlib.ext, sep=""))
#dyn.load(paste("msFile2jsfs_C++",.Platform$dynlib.ext, sep=""))
## if this gives an error: type "R CMD SHLIB msFile2jsfs_C++.cc" into shell
## R CMD SHLIB seqGenFile2jsfs_C++.cc

## different number of ss
#dyn.load(paste("../../R-eclipse/seqGenFile2jsfs_C++",.Platform$dynlib.ext, sep=""))  ## 23ss
#dyn.load(paste("seqGenFile2jsfs_30ss_C++",.Platform$dynlib.ext, sep=""))
#dyn.load(paste("../../R-eclipse/seqGenFile2jsfs_42ss_C++",.Platform$dynlib.ext, sep=""))

## Reads msout from internal memory and calculates jsfs with use of the C++
## code 'ms2jsfs_C'. The C++ code needs to be read in before. The
## result is returned as a matrix.
#Jaatha.getJSFS <- function(msout,samplesizes,anzLoci){
#	return ( matrix(.C("ms2jsfs",as.character(msout),
#							as.integer(samplesizes[1]),
#							as.integer(samplesizes[2]),
#							as.integer(anzLoci),
#							res=integer((samplesizes[1]+1)*(samplesizes[2]+1)))$res,
#					nrow=samplesizes[1]+1,ncol=samplesizes[2]+1) )
#}

#dyn.load(paste("../../R-eclipse/seqGenFile2jsfs_C++",.Platform$dynlib.ext, sep=""))  ## 23ss
#dyn.load(paste("../../R-eclipse/seqGenFile2jsfs_30ss_C++",.Platform$dynlib.ext, sep=""))
#dyn.load(paste("../../R-eclipse/seqGenFile2jsfs_42ss_C++",.Platform$dynlib.ext, sep=""))

## Reads in file 'filename' and calculates jsfs or any sumstats with use of the 
## C++ code 'msFile2jsfs_C'. The C++ code needs to be read in before. The
## result is returned as an array!!
#Jaatha.calcJSFSfromFile <- function(fileName="output",samplesizes,anzLoci,
#								nTotalSS,nAdditionalSS=0,finiteSites=FALSE){
#	resultSize <- (samplesizes[1]+1)*(samplesizes[2]+1)+ nAdditionalSS
#	if(finiteSites){
#		## 23 ss
#		return (.C( "seqFile2jsfs",as.character(fileName),
#			as.integer(samplesizes[1]),
#			as.integer(samplesizes[2]),
#			as.integer(anzLoci),
#			as.integer(resultSize),
#			res=integer(resultSize))$res )
#		
#		## 30 ss
#		return (.C("seqFile2jsfs_30ss",as.character(fileName),
#						as.integer(samplesizes[1]),
#						as.integer(samplesizes[2]),
#						as.integer(anzLoci),
#						as.integer(resultSize),
#						res=integer(resultSize))$res )	
		##42 ss
#		return (.C( "seqFile2jsfs_42ss",as.character(fileName),
#				as.integer(samplesizes[1]),
#				as.integer(samplesizes[2]),
#				as.integer(anzLoci),
#				as.integer(resultSize),
#				res=integer(resultSize))$res )
#	
#	} else{ ## infinite sites jsfs calculation	
#		return ( .C("msFile2jsfs",as.character(fileName),
#							as.integer(samplesizes[1]),
#							as.integer(samplesizes[2]),
#							as.integer(anzLoci),
#							res=integer(resultSize))$res )
#	}
#}




## Funtion to calculate the likelihood based on simulations with the
## given parameters.  Order of parameters should be the same as needed
## for the simulate-function (in Simulator.R).
Jaatha.calcLikelihood <- function(jObject, nSimulations, par){
    i <- NULL   # To make R CMD check stop complaining

	.log2("Called Jaatha.calcLikelihood()")
    par <- matrix(rep(par, each=nSimulations), nrow=nSimulations)

    sim.packages <- createSimulationPackages(par, jObject@sim.package.size)
    seeds <- generateSeeds(length(sim.packages)+1)

    # Simulate each package, maybe on different cores
    simSS  <- foreach(i = seq(along = sim.packages), .combine='rbind') %dopar% {
      set.seed(seeds[i])
      sim.pars <- .deNormalize(jObject,
                               sim.packages[[i]],
                               withoutTheta=jObject@externalTheta)
      sumStats <- dm.simSumStats(jObject@dm, sim.pars, jObject@sum.stats.func)
      return(sumStats)
    }
    set.seed(seeds[length(seeds)])

    # Average the values of each summary statistic
    simSS <- apply(simSS, 2, mean)

	.log2("Calculating Likelihood...")
	logL <- 0
	simSS[simSS==0] <- 0.5
	for (s in 1:jObject@nTotalSumstat){
		logL <- logL + jObject@sumStats[s] * log(simSS[s]) - simSS[s] - .logfac(jObject@sumStats[s])
	}
	.log2("Finished Jaatha.calcLikelihood(). Return:",logL)
	return(logL)
}

## Reads msout from internal memory and calculates jsfs with use of the C++
## code 'ms2jsfs_C'. The C++ code needs to be read in before. The
## result is returned as a matrix.
#getJSFS <- function(msout,samplesizes,anzLoci){
#	return ( matrix(.C("ms2jsfs",as.character(msout),
#							as.integer(samplesizes[1]),
#							as.integer(samplesizes[2]),
#							as.integer(anzLoci),
#							res=integer((samplesizes[1]+1)*(samplesizes[2]+1)))$res,
#					nrow=samplesizes[1]+1,ncol=samplesizes[2]+1) )
#}


## Returns the logarithm of the factorial of k. Recursively implemented.
.logfac <- function(k) {
	if(k<2) return(0)
	if(k>10) return(lgamma(k+1))
	else return(log(k)+.logfac(k-1))
}


##Calculates the joint site frequency spectrum for the sequence data
##(dnadaten) and the information (index of outgroup sequence[[1]], of pop1[[2]],
##of pop2[[3]]) on the data given in datenInfo. Positions containing
##"-","_", "N", "n", and "?" get excluded.
.jsfsData <- function(dnadaten,datenInfo) {
	daten <- matrix(unlist(dnadaten),byrow=TRUE,nrow=length(dnadaten))
	jsfs <- array(0,dim=c(length(datenInfo[[2]])+1,length(datenInfo[[3]])+1))
	gapfree=apply(daten,2,
			function(x) sum(x=="-" | x=="_"| x=="?"|x=="N"|x=="n")==0 )
	message("Sequence length without special characters is: ",sum(gapfree))
	for( i in 1:ncol(daten)) {
		if(gapfree[i]) {
			a <- sum(daten[datenInfo[[2]],i]!=daten[datenInfo[[1]],i]) 
			b <- sum(daten[datenInfo[[3]],i]!=daten[datenInfo[[1]],i])
			jsfs[a+1,b+1] <- jsfs[a+1,b+1]+1
		}
	}
	return(jsfs)
}



##Reads in all files in fileFolder of fileType and return a list of
##all the sequence data. Each file contains a sepertate locus.
.readIn <- function(object){
	# reads in all the nexus files in the folder (linux sorted read in)
	datei <- system(paste("ls ",object@fileFolder,"*.",object@fileType,sep=""),
			intern=TRUE)
	#print(datei)
	daten <- list()
	
	if(object@fileType=="nex"){
		for (s in 1:object@nLoci){
			#message("Reading in ",(datei[s]))
			daten[[s]] <- Jaatha.readNexus(datei[s])   
		}
	}
	else if (object@fileType=="phy"){
		for (s in 1:object@nLoci){
			daten[[s]] <- Jaatha.readPhylip(datei[s])
		}
	}
	else{
		cat("File format can not be read in. Please use nex, phy or msoutput as 
						file formats!")
	}	
	return(daten)     
	#print(daten)
}


##Function to read in phylip format files and return the sequences.
Jaatha.readPhylip <- function(file){
	lines <- readLines(file)
	#print(lines)
	nInd <- as.numeric(strsplit(lines[1]," ")[[1]][1])
	nPos <- as.numeric(strsplit(lines[1]," ")[[1]][2])
	daten <- list()
	for (i in 1:nInd){
		daten[[i]] <- strsplit(substr(lines[i+1],11,nPos+10),"")[[1]]
	}
	return(daten)  
}


##Function to read in nexus format files and return the sequences as a 
## list[[sample]][pos].
## Assumes the following order of the KEYS in the file: NTAX,NCHAR,FORMAT,MATRIX
Jaatha.readNexus <- function(file, debug=FALSE){
	if (file.exists(file)){
		lines <- readLines(file)
		fileLength <- length(lines)
		#print(lines)
		l <- 1
		## read out ntax
		while(regexpr("NTAX",lines[l],ignore.case=T)<0 & l<=fileLength){
			l<-l+1
		}	
		nInd <- as.numeric(substr(lines[l],start=.getFirstPos("=",lines[l])+1,
						stop=.getFirstPos(";",lines[l])-1))
		cat("File ",file, " contains\n\t\t",nInd, " samples with each ",sep="")
		## read out ntax
		while(regexpr("NCHAR",lines[l],ignore.case=T)<0 & l<=fileLength){
			l<-l+1
		}
		nPos <- as.numeric(substr(lines[l],start=.getFirstPos("=",lines[l])+1,
						stop=.getFirstPos(";",lines[l])-1))
		cat(nPos," positions.\n")
		## read out matchChar
		while(regexpr("FORMAT",lines[l],ignore.case=T)<0 & l<=fileLength){
			l<-l+1
		}
		p <- regexpr("MATCHCHAR",lines[l],ignore.case=T)
		if (p>0){ 
			matchChar <-  gsub("\\s","",substr(lines[l],start=p+10, stop=p+11))
		}else{ 
			matchChar <- ""
		}
		## find line where data starts
		while(regexpr("MATRIX",lines[l],ignore.case=T) < 0 & l<=fileLength){
			l<-l+1
		}
		l <- l+1
		daten <- list()
		## if there is a match character
		if (matchChar!=""){   ## read in first inidi seperately bc it contains names
			p <- .getSeqStart(lines[l])
			daten[[1]] <- strsplit(substr(lines[l],p,p+70),"")[[1]]
			left2read <- nPos -70
			while (left2read>0){
				l<- l+1
				daten[[1]] <- c(daten[[1]],strsplit(lines[l],"")[[1]])
				left2read <- left2read -70
				#print(left2read)
			}		
			if (debug){ print(length(daten[[1]]))}else{}
			l <- l+1
			for (i in 2:nInd){
				if (debug) message("************",i)
				l <- l+1
				left2read <- nPos 
				p <- .getSeqStart(lines[l])
				if (debug) print(substr(lines[l],1,p-1))
				data <- strsplit(substr(lines[l],p,p+70),"")[[1]]
				if (debug) print(data)
				f <- 1
				daten[[i]] <- ""
				while (left2read >0){    
					if(length(data)>0){
						for (p in 1:length(data)){
							if (data[p]==matchChar){
								daten[[i]][f] <- daten[[1]][f]						
							}else{
								daten[[i]][f] <- data[p]			
							}
							f <- f+1
						}
						left2read <- left2read -length(data)
					}else{}
					l <- l+1
					data <- strsplit(substr(lines[l],1,70),"")[[1]]
					if (debug){ cat("left2read:",left2read,"\n")}else{}
				}
				if (debug){ cat("seqLength:",length(daten[[i]]),"\n")}
			}
		}else{    # if there is no matchChar
			for (i in 1:nInd){
				if (debug) message("************",i)
				left2read <- nPos 
				p <- .getSeqStart(lines[l])
				if (debug) print(substr(lines[l],1,p-1))
				data <- strsplit(substr(lines[l],p,p+70),"")[[1]]
				if (debug) print(data)
				daten[[i]]<-""
				while (left2read >0){    
					if(length(data)>0){
						daten[[i]]<- c(daten[[i]],data)	
						left2read <- left2read -length(data)
					}else{}
					l <- l+1
					data <- strsplit(substr(lines[l],1,70),"")[[1]]
					if (debug){ cat("left2read:",left2read,"\n")}
				}
				daten[[i]] <- daten[[i]][-1]  ## to remove the initialization ""
				if (debug) cat("length of sample",i,":",length(daten[[i]]))
				l <- l+1
			}
		}
		return(daten)  
	} else{
		stop("Error: ",file, "not found.\n")
	}
}

## Returns the starting position of the sequence
.getSeqStart <- function(line){
	p <- 2
	s <- strsplit(line,"\\s")[[1]]
	while(p<=length(s) & (s[p]=="")){
		p <- p+1
	}
	if (p<=length(s)){ 
		return(p+nchar(s[1]))
	}else{
		stop("There is no white space character in this line:\n",line)
		return (-1)
	}
}

## Return the index of first occurrence of 'pattern' in 'source'.
.getFirstPos <- function(pattern,source,caseInsensitive=F){
	s <- strsplit(source,"")[[1]]
	pos <- grep(pattern,s,ignore.case=caseInsensitive)
	if(length(pos)>1){
		message(pattern," occurs ",length(pos),
				" times. Only first occurrence reported.")
		return(pos[1])
	}
	else if(length(pos)==1){
		return(pos)
	}else{}
}

## Function to calculate the joint site frequency spectrum from lines
## of a msoutput file. Sample sizes for both populations and the
## number of loci need to be defined.
#jsfsMS <- function(msoutputLines,nSamplePop1,nSamplePop2,nLoci){
#  jsfs <- array(rep(0,(nSamplePop1+1)*(nSamplePop2+1)),
#                c((nSamplePop1+1),(nSamplePop2+1)))
#  for(z in 0:(nLoci-1)*(4+nSamplePop1+nSamplePop2)+7) {
#    locus <- as.matrix(
#                       data.frame(strsplit(msoutputLines[z:(z+nSamplePop1+nSamplePop2-1)],"")))
#    for(i in 1:nrow(locus)) {
#      x <- sum(locus[i,1:nSamplePop1]=="1")+1
#      y <- sum(locus[i,(nSamplePop1+1):(nSamplePop1+nSamplePop2)]=="1")+1
#      jsfs[x,y] <- jsfs[x,y] + 1
#    }              
#  }
#  jsfs
#}

