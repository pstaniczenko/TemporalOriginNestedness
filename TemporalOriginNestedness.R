# Temporal origin of nestedness in interaction networks
# version 1 -- 2023-09-21
# pstaniczenko@brooklyn.cuny.edu
# Please cite:
# Staniczenko, P.P.A. & Panja, D. (2023). Temporal origin of nestedness in interaction networks. PNAS Nexus, 2, pgad412.

library(Rlab)
library(poisbinom)

#############################################
## Calculate TP, FP, TN, FN, sensitivity, specificity, precision, recall, F-score, Informedness
#############################################
## Compares empirical and model matrix and returns
## TP, FP, TN, FN, sensitivity, specificity, precision, recall, F-score
CalculateTPetc <- function(Empirical, Model){
	TP <- sum((Empirical==1)*(Model==1))
	FP <- sum((Empirical==0)*(Model==1))
	TN <- sum((Empirical==0)*(Model==0))
	FN <- sum((Empirical==1)*(Model==0))
	sensitivity <- TP/(TP+FN)
	specificity <- TN/(TN+FP)
	precision <- TP/(TP+FP)
	recall <- sensitivity
	Fscore <- 2 * ((precision*recall)/(precision+recall))
	Informedness <- sensitivity + specificity - 1
	Phi <- ((TP*TN)-(FN*FP))/sqrt(as.double(TP+FP)*as.double(FP+TN)*as.double(TP+FN)*as.double(FN+TN))
	return(list(TP=TP, FP=FP, TN=TN, FN=FN, sensitivity=sensitivity, specificity=specificity, recall=recall, Fscore=Fscore, Informedness=Informedness, Phi=Phi))
}

#############################################
## Calculate spectral radius
#############################################
## Takes a (possibly rectangular) incidence 
## matrix A where A_ij > 0 means that species i (row)
## interacts with species j (column) and return
## the dominant eigenvalue of the corresponding
## adjacency matrix
DominantEigenFromIncidence <- function(A){
  NumP <- dim(A)[1]
  NumA <- dim(A)[2]
  S <- NumP + NumA
  Adj <- matrix(0, S, S)
  Adj[1:NumP, (NumP + 1):S] <- A
  Adj <- Adj + t(Adj)
  return(Re(eigen(Adj)$value[1]))
}

#############################################
# TEST FOR  NESTEDNESS -- Bernoulli probability
#############################################
## Try matrices which have the same number
## of plants and animals of the original
## and with interactions occurring with probability C = L/(P*A)
TestNestednessBernoulli <- function(A, HowMany){
  ## Dominant eigenvalue of the original matrix (after removing redundant rows and columns)
  A_rowSums <- rowSums(A)
  A_colSums <- colSums(A)
  A_extract <- A[which(A_rowSums > 0, arr.ind=TRUE), which(A_colSums > 0, arr.ind=TRUE)]
  EigenA <- DominantEigenFromIncidence(A_extract)
  EigenA_max <- sqrt((2*sum(A_extract))-(dim(A_extract)[1]+dim(A_extract)[2])+1)
  EigenA_normalized <- EigenA/EigenA_max
  ## Vector to store eigen for randomizations
  EigenRandomizations <- rep(-1.0, HowMany)
  ## Vectors to store sensitivity, specificity, Fscore values
  sensitivity_vec <- rep(-1.0, HowMany)
  specificity_vec <- rep(-1.0, HowMany)
  Fscore_vec <- rep(-1.0, HowMany)
  Informedness_vec <- rep(-1.0, HowMany)
  Phi_vec <- rep(-1.0, HowMany)
  #
  NumP <- dim(A_extract)[1]
  NumA <- dim(A_extract)[2]
  NumLinks <- sum(A_extract)
  for (Tried in 1:HowMany){
  	## Bernoulli probability
  	Bernoulli_probability <- NumLinks/(NumP * NumA)
    ## Build a matrix by drawing from Bernoulli distribution
    B <- matrix(rbern(NumP * NumA, Bernoulli_probability), nrow=NumP, ncol=NumA)
    ## Compute its eigenvalue
    EigenRandomizations[Tried] <- DominantEigenFromIncidence(B)/sqrt((2*sum(B))-(dim(B)[1]+dim(B)[2])+1)
    # Store Fscore etc.
    TPetc <- CalculateTPetc(A_extract, B)
    sensitivity_vec[Tried] <- TPetc$sensitivity
    specificity_vec[Tried] <- TPetc$specificity
    Fscore_vec[Tried] <- TPetc$Fscore
    Informedness_vec[Tried] <- TPetc$Informedness
    Phi_vec[Tried] <- TPetc$Phi
  }
  ## Probability of finding a matrix with a larger eigenvalue at random
  ProbabilityLargerEigenvalue <- sum(EigenRandomizations >= EigenA_normalized) / HowMany
  ## z-score
  nullmodel_mean <- mean(EigenRandomizations)
  nullmodel_sd <- sd(EigenRandomizations)
  zscore <- (EigenA_normalized - nullmodel_mean)/nullmodel_sd
  nullmodel_var <- var(EigenRandomizations)
  
  return(list(Success=1,
              A=A, A_extract=A_extract, Eigen_normalized=EigenA_normalized,
              EigenRand=EigenRandomizations,
              Probability=ProbabilityLargerEigenvalue,
              zscore=zscore,
              nullmodel_mean=nullmodel_mean, nullmodel_sd=nullmodel_sd, nullmodel_var=nullmodel_var,
              sensitivity_mean=mean(sensitivity_vec),
              sensitivity_var=var(sensitivity_vec),
              specificity_mean=mean(specificity_vec),
              specificity_var=var(specificity_vec),
              Fscore_mean=mean(Fscore_vec),
              Fscore_var=var(Fscore_vec),
              Informedness_mean=mean(Informedness_vec),
              Informedness_var=var(Informedness_vec),
              Phi_mean=mean(Phi_vec),
              Phi_var=var(Phi_vec),
              Test="Bernoulli"))
}

#############################################
# TEST FOR  NESTEDNESS -- Average degree distribution probability
#############################################
## Try matrices which have the same number
## of plants and animals of the original
## and with interactions occurring with probability (k_i/n + d_i/m)/2
TestNestednessDegreeDist <- function(A, HowMany){
  ## Dominant eigenvalue of the original matrix (after removing redundant rows and columns)
  A_rowSums <- rowSums(A)
  A_colSums <- colSums(A)
  A_extract <- A[which(A_rowSums > 0, arr.ind=TRUE), which(A_colSums > 0, arr.ind=TRUE)]
  EigenA <- DominantEigenFromIncidence(A_extract)
  EigenA_max <- sqrt((2*sum(A_extract))-(dim(A_extract)[1]+dim(A_extract)[2])+1)
  EigenA_normalized <- EigenA/EigenA_max
  ## Vector to store eigen for randomizations
  EigenRandomizations <- rep(-1.0, HowMany)
  ## Vectors to store sensitivity, specificity, Fscore values
  sensitivity_vec <- rep(-1.0, HowMany)
  specificity_vec <- rep(-1.0, HowMany)
  Fscore_vec <- rep(-1.0, HowMany)
  Informedness_vec <- rep(-1.0, HowMany)
  Phi_vec <- rep(-1.0, HowMany)
  #
  NumP <- dim(A_extract)[1]
  NumA <- dim(A_extract)[2]
  NumLinks <- sum(A_extract)
  for (Tried in 1:HowMany){
  	B <- matrix(rep(0, NumP * NumA), nrow=NumP, ncol=NumA)
  	for (i in 1:NumP){
  		for (j in 1:NumA){
  			Bernoulli_probability <- 0.5 * ((rowSums(A_extract)[i]/NumA) + (colSums(A_extract)[j]/NumP))
  			## Build a matrix by drawing from Bernoulli distribution
    		B[i,j] <- rbern(1, Bernoulli_probability)
  		}
  	}
    ## Compute its eigenvalue
    EigenRandomizations[Tried] <- DominantEigenFromIncidence(B)/sqrt((2*sum(B))-(dim(B)[1]+dim(B)[2])+1)
    # Store Fscore etc.
    TPetc <- CalculateTPetc(A_extract, B)
    sensitivity_vec[Tried] <- TPetc$sensitivity
    specificity_vec[Tried] <- TPetc$specificity
    Fscore_vec[Tried] <- TPetc$Fscore
    Informedness_vec[Tried] <- TPetc$Informedness
    Phi_vec[Tried] <- TPetc$Phi
  }
  ## Probability of finding a matrix with a larger eigenvalue at random
  ProbabilityLargerEigenvalue <- sum(EigenRandomizations >= EigenA_normalized) / HowMany
  ## z-score
  nullmodel_mean <- mean(EigenRandomizations)
  nullmodel_sd <- sd(EigenRandomizations)
  zscore <- (EigenA_normalized - nullmodel_mean)/nullmodel_sd
  nullmodel_var <- var(EigenRandomizations)
  
  return(list(Success=1,
              A=A, A_extract=A_extract, Eigen_normalized=EigenA_normalized,
              EigenRand=EigenRandomizations,
              Probability=ProbabilityLargerEigenvalue,
              zscore=zscore,
              nullmodel_mean=nullmodel_mean, nullmodel_sd=nullmodel_sd, nullmodel_var=nullmodel_var,
              sensitivity_mean=mean(sensitivity_vec),
              sensitivity_var=var(sensitivity_vec),
              specificity_mean=mean(specificity_vec),
              specificity_var=var(specificity_vec),
              Fscore_mean=mean(Fscore_vec),
              Fscore_var=var(Fscore_vec),
              Informedness_mean=mean(Informedness_vec),
              Informedness_var=var(Informedness_vec),
              Phi_mean=mean(Phi_vec),
              Phi_var=var(Phi_vec),
              Test="Degree Dist"))
}

#############################################
# TEST FOR  NESTEDNESS -- Poisson binomial with SINGLE probability
#############################################
## A phenological null model for nestedness, only aggregate network available
TestNestednessPoissonbinomialSINGLE <- function(A, cop, samplingfreq, HowMany){
  ## A is empirical matrix and cop is a co-presence matrix with max element <= samplingfreq
  ## Dominant eigenvalue of the original matrix (after removing redundant rows and columns)
  A_rowSums <- rowSums(A)
  A_colSums <- colSums(A)
  A_extract <- A[which(A_rowSums > 0, arr.ind=TRUE), which(A_colSums > 0, arr.ind=TRUE)]
  EigenA <- DominantEigenFromIncidence(A_extract)
  EigenA_max <- sqrt((2*sum(A_extract))-(dim(A_extract)[1]+dim(A_extract)[2])+1)
  EigenA_normalized <- EigenA/EigenA_max
  ## calculate mle single Bernoulli probability
  zeros_empirical <- sum(A_extract==0)
  cop_extract <- cop[which(A_rowSums > 0, arr.ind=TRUE), which(A_colSums > 0, arr.ind=TRUE)]
  poly_coeffs <- rep(-1, samplingfreq + 1)
  for (n in 1:length(poly_coeffs)){
  	poly_coeffs[n] <- sum(cop_extract==(n-1)) # obtain cofficients for polynomial
  }
  if (poly_coeffs[1] >= zeros_empirical){
  	print("Model is trivial: All zeros in empirical matrix are due to species never co-occurring!!")
  }
  poly_coeffs[1] <- poly_coeffs[1] - zeros_empirical
  Bernoulli_probability_q <- polyroot(poly_coeffs)
  # ensure the solution with only positive real part is used
  Bernoulli_probability_q <- Bernoulli_probability_q[which(abs(Im(Bernoulli_probability_q)) < 0.000000001, arr.ind=TRUE)]
  Bernoulli_probability_q <- Bernoulli_probability_q[which(Re(Bernoulli_probability_q) > 0, arr.ind=TRUE)]
  Bernoulli_probability_q <- Re(Bernoulli_probability_q[1])
  Bernoulli_probability <- 1 - Bernoulli_probability_q
  ## Vector to store eigen for randomizations
  EigenRandomizations <- rep(-1.0, HowMany)
  ## Vectors to store sensitivity, specificity, Fscore values
  sensitivity_vec <- rep(-1.0, HowMany)
  specificity_vec <- rep(-1.0, HowMany)
  Fscore_vec <- rep(-1.0, HowMany)
  Informedness_vec <- rep(-1.0, HowMany)
  Phi_vec <- rep(-1.0, HowMany)
  #
  NumP <- dim(A_extract)[1]
  NumA <- dim(A_extract)[2]
  NumLinks <- sum(A_extract)
  for (Tried in 1:HowMany){
  	## Generate synthetic matrix
	B <- matrix(rep(0, NumP * NumA), nrow=NumP, ncol=NumA)
	for (i in 1:NumP){
	  for (j in 1:NumA){
	  	poisbinom_probs <- c(rep(Bernoulli_probability, cop_extract[i,j]), rep(0, samplingfreq-cop_extract[i,j]))
		B[i, j] <- rpoisbinom(1, poisbinom_probs)
	  }
    }
  B[B>0] <- 1
  ## Compute its eigenvalue
  EigenRandomizations[Tried] <- DominantEigenFromIncidence(B)/sqrt((2*sum(B))-(dim(B)[1]+dim(B)[2])+1)
  # Store Fscore etc.
  TPetc <- CalculateTPetc(A_extract, B)
  sensitivity_vec[Tried] <- TPetc$sensitivity
  specificity_vec[Tried] <- TPetc$specificity
  Fscore_vec[Tried] <- TPetc$Fscore
  Informedness_vec[Tried] <- TPetc$Informedness
  Phi_vec[Tried] <- TPetc$Phi
  }
  ## Probability of finding a matrix with a larger eigenvalue at random
  ProbabilityLargerEigenvalue <- sum(EigenRandomizations >= EigenA_normalized) / HowMany
  ## z-score
  nullmodel_mean <- mean(EigenRandomizations)
  nullmodel_sd <- sd(EigenRandomizations)
  zscore <- (EigenA_normalized - nullmodel_mean)/nullmodel_sd
  nullmodel_var <- var(EigenRandomizations)
  
  return(list(Success=1,
              A=A, A_extract=A_extract, Eigen_normalized=EigenA_normalized,
              EigenRand=EigenRandomizations,
              Probability=ProbabilityLargerEigenvalue,
              Bernoulli_probability=Bernoulli_probability,
              zscore=zscore,
              nullmodel_mean=nullmodel_mean, nullmodel_sd=nullmodel_sd, nullmodel_var=nullmodel_var, 
              sensitivity_mean=mean(sensitivity_vec),
              sensitivity_var=var(sensitivity_vec),
              specificity_mean=mean(specificity_vec),
              specificity_var=var(specificity_vec),
              Fscore_mean=mean(Fscore_vec),
              Fscore_var=var(Fscore_vec),
              Informedness_mean=mean(Informedness_vec),
              Informedness_var=var(Informedness_vec),
              Phi_mean=mean(Phi_vec),
              Phi_var=var(Phi_vec),
              Test="Poisson binomial SINGLE"))
}

#############################################
# MAIN

# start timer
start_time <- Sys.time()

# number synthetic matrices for nestedness analysis
HowMany <- 100

# parameters set by data
nrows <- 46 # number of rows in entire aggregate matrix
ncols <- 93 # number of columns in entire aggregate matrix
n_days <- 106 # number of daily matrices
temporal_aggreation <- "weeks" # chose from "weeks", "months", "years"

# load all daily matrices
for (n in 1:n_days){
	mydatafolder <- "ExampleData/daily_matrices/"
	myfilename <- paste("daily_adjacency_", toString(n), ".csv", sep="")
	myfile <- paste(mydatafolder, myfilename, sep="")
	mymatrix <- read.csv(myfile, header=FALSE, sep=",")
	mymatrix_dim <- dim(mymatrix)
	mymatrix <- mapply(mymatrix, FUN=as.numeric)
	mymatrix <- matrix(data=mymatrix, nrow=mymatrix_dim[1], ncol=mymatrix_dim[2])
	mymatrix <- mymatrix[1:nrows, 1:ncols]
	mymatrix[mymatrix>0] <- 1 # CONVERT MATRICES TO BINARY
	assign(paste("daily_adjacency_", toString(n), sep=""), mymatrix)
}

# generate vector of which week (column entry) each day (row) belongs to
if (temporal_aggreation == "weeks"){
	n_weeks <- 24
	myfilename <- "ExampleData/days_in_weeks.csv" # first col is week id and second col is num days in week
} else if (temporal_aggreation == "months"){
	n_weeks <- 6
	myfilename <- "ExampleData/days_in_months.csv" 
} else if (temporal_aggreation == "years"){
	n_weeks <- 3
	myfilename <- "ExampleData/days_in_years.csv" 
}
week_numberdays <- read.csv(myfilename, header=FALSE, sep=",")
daybelongtoweek <- rep(-1, n_days)
day_counter <- 0
for (week_counter in 1:n_weeks){
	daysinweek <- week_numberdays[week_counter, 2]
	for (dayinweek in 1:daysinweek){
		day_counter <- day_counter + 1
		daybelongtoweek[day_counter] <- week_counter
	}
}

# vectors for results
results_empirical_normalized_nestedness <- rep(-1, n_weeks)
# Bernoulli, i.e., Erdos-Renyi random graph
results_Bernoullinullmodel_mean <- rep(-1, n_weeks)
results_Bernoullinullmodel_var <- rep(-1, n_weeks)
results_Bernoullinullmodel_zscore <- rep(-1, n_weeks)
results_Bernoullinullmodel_sensitivity_mean <- rep(-1, n_weeks)
results_Bernoullinullmodel_sensitivity_var <- rep(-1, n_weeks)
results_Bernoullinullmodel_specificity_mean <- rep(-1, n_weeks)
results_Bernoullinullmodel_specificity_var <- rep(-1, n_weeks)
results_Bernoullinullmodel_Fscore_mean <- rep(-1, n_weeks)
results_Bernoullinullmodel_Fscore_var <- rep(-1, n_weeks)
results_Bernoullinullmodel_Informedness_mean <- rep(-1, n_weeks)
results_Bernoullinullmodel_Informedness_var <- rep(-1, n_weeks)
results_Bernoullinullmodel_Phi_mean <- rep(-1, n_weeks)
results_Bernoullinullmodel_Phi_var <- rep(-1, n_weeks)
# Degree distribution model
results_degreedistnullmodel_mean <- rep(-1, n_weeks)
results_degreedistnullmodel_var <- rep(-1, n_weeks)
results_degreedistnullmodel_zscore <- rep(-1, n_weeks)
results_degreedistnullmodel_sensitivity_mean <- rep(-1, n_weeks)
results_degreedistnullmodel_sensitivity_var <- rep(-1, n_weeks)
results_degreedistnullmodel_specificity_mean <- rep(-1, n_weeks)
results_degreedistnullmodel_specificity_var <- rep(-1, n_weeks)
results_degreedistnullmodel_Fscore_mean <- rep(-1, n_weeks)
results_degreedistnullmodel_Fscore_var <- rep(-1, n_weeks)
results_degreedistnullmodel_Informedness_mean <- rep(-1, n_weeks)
results_degreedistnullmodel_Informedness_var <- rep(-1, n_weeks)
results_degreedistnullmodel_Phi_mean <- rep(-1, n_weeks)
results_degreedistnullmodel_Phi_var <- rep(-1, n_weeks)
# Phenology model
results_phenologynullmodelSINGLE_mean <- rep(-1, n_weeks)
results_phenologynullmodelSINGLE_var <- rep(-1, n_weeks)
results_phenologynullmodelSINGLE_zscore <- rep(-1, n_weeks)
results_phenologynullmodelSINGLE_sensitivity_mean <- rep(-1, n_weeks)
results_phenologynullmodelSINGLE_sensitivity_var <- rep(-1, n_weeks)
results_phenologynullmodelSINGLE_specificity_mean <- rep(-1, n_weeks)
results_phenologynullmodelSINGLE_specificity_var <- rep(-1, n_weeks)
results_phenologynullmodelSINGLE_Fscore_mean <- rep(-1, n_weeks)
results_phenologynullmodelSINGLE_Fscore_var <- rep(-1, n_weeks)
results_phenologynullmodelSINGLE_Informedness_mean <- rep(-1, n_weeks)
results_phenologynullmodelSINGLE_Informedness_var <- rep(-1, n_weeks)
results_phenologynullmodelSINGLE_Phi_mean <- rep(-1, n_weeks)
results_phenologynullmodelSINGLE_Phi_var <- rep(-1, n_weeks)
# Misc
results_empirical_connectance <- rep(-1, n_weeks)
results_phenology_prob <- rep(-1, n_weeks)

for (myweek in 1:n_weeks){

start_time_week <- Sys.time()
print(paste("Week/Month/Year ", toString(myweek), sep=""))

# days in week of interest
mydays <- which(daybelongtoweek==myweek, arr.ind=TRUE)

# Get binary aggregate network for week and co-presence matrix for week
aggregate_network_week <- matrix(rep(0, nrows * ncols), nrow=nrows, ncol=ncols)
copresence_network_week <- matrix(rep(0, nrows * ncols), nrow=nrows, ncol=ncols)
# loop over days in week of interest
for (n in 1:length(mydays)){
	myday <- mydays[n]
	assign("mymatrix", get(paste("daily_adjacency_", toString(myday), sep="")))
	aggregate_network_week <- aggregate_network_week + mymatrix
	for (i in 1:nrows){
		for (j in 1:ncols){
			if (rowSums(mymatrix)[i]>0 && colSums(mymatrix)[j]>0){ # if two actors are present on same day
				copresence_network_week[i,j] <- copresence_network_week[i,j] + 1
			}
		}
	}
}
aggregate_network_week_weights <- aggregate_network_week
aggregate_network_week[aggregate_network_week>0] <- 1 # map to binary matrix

# now we have all information to run phenology model
# aggregate_network_week is a binary matrix of all interactions observed during a given week, myweek
# copresence_network_week is a matrix of the number of times each pair of actors is co-present during a given week, myweek

results_Bernoullinullmodel <- TestNestednessBernoulli(aggregate_network_week, HowMany)
#
results_empirical_normalized_nestedness[myweek] <- results_Bernoullinullmodel$Eigen_normalized
#
results_Bernoullinullmodel_mean[myweek] <- results_Bernoullinullmodel$nullmodel_mean
results_Bernoullinullmodel_var[myweek] <- results_Bernoullinullmodel$nullmodel_var
results_Bernoullinullmodel_zscore[myweek] <- results_Bernoullinullmodel$zscore
results_Bernoullinullmodel_sensitivity_mean[myweek] <- results_Bernoullinullmodel$sensitivity_mean
results_Bernoullinullmodel_sensitivity_var[myweek] <- results_Bernoullinullmodel$sensitivity_var
results_Bernoullinullmodel_specificity_mean[myweek] <- results_Bernoullinullmodel$specificity_mean
results_Bernoullinullmodel_specificity_var[myweek] <- results_Bernoullinullmodel$specificity_var
results_Bernoullinullmodel_Fscore_mean[myweek] <- results_Bernoullinullmodel$Fscore_mean
results_Bernoullinullmodel_Fscore_var[myweek] <- results_Bernoullinullmodel$Fscore_var
results_Bernoullinullmodel_Informedness_mean[myweek] <- results_Bernoullinullmodel$Informedness_mean
results_Bernoullinullmodel_Informedness_var[myweek] <- results_Bernoullinullmodel$Informedness_var
results_Bernoullinullmodel_Phi_mean[myweek] <- results_Bernoullinullmodel$Phi_mean
results_Bernoullinullmodel_Phi_var[myweek] <- results_Bernoullinullmodel$Phi_var
#
results_degreedistnullmodel <- TestNestednessDegreeDist(aggregate_network_week, HowMany)
results_degreedistnullmodel_mean[myweek] <- results_degreedistnullmodel$nullmodel_mean
results_degreedistnullmodel_var[myweek] <- results_degreedistnullmodel$nullmodel_var
results_degreedistnullmodel_zscore[myweek] <- results_degreedistnullmodel$zscore
results_degreedistnullmodel_sensitivity_mean[myweek] <- results_degreedistnullmodel$sensitivity_mean
results_degreedistnullmodel_sensitivity_var[myweek] <- results_degreedistnullmodel$sensitivity_var
results_degreedistnullmodel_specificity_mean[myweek] <- results_degreedistnullmodel$specificity_mean
results_degreedistnullmodel_specificity_var[myweek] <- results_degreedistnullmodel$specificity_var
results_degreedistnullmodel_Fscore_mean[myweek] <- results_degreedistnullmodel$Fscore_mean
results_degreedistnullmodel_Fscore_var[myweek] <- results_degreedistnullmodel$Fscore_var
results_degreedistnullmodel_Informedness_mean[myweek] <- results_degreedistnullmodel$Informedness_mean
results_degreedistnullmodel_Informedness_var[myweek] <- results_degreedistnullmodel$Informedness_var
results_degreedistnullmodel_Phi_mean[myweek] <- results_degreedistnullmodel$Phi_mean
results_degreedistnullmodel_Phi_var[myweek] <- results_degreedistnullmodel$Phi_var
#
results_phenologynullmodelSINGLE <- TestNestednessPoissonbinomialSINGLE(aggregate_network_week, copresence_network_week, samplingfreq=length(mydays), HowMany)
results_phenologynullmodelSINGLE_mean[myweek] <- results_phenologynullmodelSINGLE$nullmodel_mean
results_phenologynullmodelSINGLE_var[myweek] <- results_phenologynullmodelSINGLE$nullmodel_var
results_phenologynullmodelSINGLE_zscore[myweek] <- results_phenologynullmodelSINGLE$zscore
results_phenology_prob[myweek] <- results_phenologynullmodelSINGLE$Bernoulli_probability
results_phenologynullmodelSINGLE_sensitivity_mean[myweek] <- results_phenologynullmodelSINGLE$sensitivity_mean
results_phenologynullmodelSINGLE_sensitivity_var[myweek] <- results_phenologynullmodelSINGLE$sensitivity_var
results_phenologynullmodelSINGLE_specificity_mean[myweek] <- results_phenologynullmodelSINGLE$specificity_mean
results_phenologynullmodelSINGLE_specificity_var[myweek] <- results_phenologynullmodelSINGLE$specificity_var
results_phenologynullmodelSINGLE_Fscore_mean[myweek] <- results_phenologynullmodelSINGLE$Fscore_mean
results_phenologynullmodelSINGLE_Fscore_var[myweek] <- results_phenologynullmodelSINGLE$Fscore_var
results_phenologynullmodelSINGLE_Informedness_mean[myweek] <- results_phenologynullmodelSINGLE$Informedness_mean
results_phenologynullmodelSINGLE_Informedness_var[myweek] <- results_phenologynullmodelSINGLE$Informedness_var
results_phenologynullmodelSINGLE_Phi_mean[myweek] <- results_phenologynullmodelSINGLE$Phi_mean
results_phenologynullmodelSINGLE_Phi_var[myweek] <- results_phenologynullmodelSINGLE$Phi_var

print(Sys.time() - start_time_week)

} # end of looping over weeks

####
# Calculate normalized delta-nestedness
# specify variables
RG_mean <- results_Bernoullinullmodel_mean
RG_var <- results_Bernoullinullmodel_var
DD_mean <- results_degreedistnullmodel_mean
DD_var <- results_degreedistnullmodel_var
PhS_mean <- results_phenologynullmodelSINGLE_mean
PhS_var <- results_phenologynullmodelSINGLE_var
Empirical <- results_empirical_normalized_nestedness

remove_outlier_points <- 1 # remove outlier points (>3 standard deviations from mean before error propagation): yes (1) or no (0)
if (remove_outlier_points==1){
	# DD
	delta_nestedness_norm_DD_vec <- (DD_mean - RG_mean)/(Empirical - RG_mean)
	outlier_points <- which(abs((delta_nestedness_norm_DD_vec-mean(delta_nestedness_norm_DD_vec))/(3*sd(delta_nestedness_norm_DD_vec)))>1, arr.ind=TRUE)
	if (length(outlier_points)==0){
		outlier_points <- -1 * c(1:length(delta_nestedness_norm_DD_vec))
	}
	RG_mean <- RG_mean[c(-1 * outlier_points)]
	RG_var <- RG_var[c(-1 * outlier_points)]
	DD_mean <- DD_mean[c(-1 * outlier_points)]
	DD_var <- DD_var[c(-1 * outlier_points)]
	PhS_mean <- PhS_mean[c(-1 * outlier_points)]
	PhS_var <- PhS_var[c(-1 * outlier_points)]
	Empirical <- Empirical[c(-1 * outlier_points)]
	# error propagation (in terms of variance until last step)
	delta_nestedness_norm_DD_vec <- (DD_mean - RG_mean)/(Empirical - RG_mean)
	delta_nestedness_norm_DD <- mean(delta_nestedness_norm_DD_vec)
	numerator_DD <- DD_mean - RG_mean
	numerator_error_DD <- DD_var + RG_var
	denominator_DD <- Empirical - RG_mean
	denominator_error_DD <- RG_var
	overall_error_DD <- delta_nestedness_norm_DD_vec^2 * ((numerator_error_DD/numerator_DD^2) + (denominator_error_DD/denominator_DD^2))
	overall_error_DD <- mean(overall_error_DD)
	overall_error_DD <- sqrt(overall_error_DD)
	# PhS
	RG_mean <- results_Bernoullinullmodel_mean
	RG_var <- results_Bernoullinullmodel_var
	DD_mean <- results_degreedistnullmodel_mean
	DD_var <- results_degreedistnullmodel_var
	PhS_mean <- results_phenologynullmodelSINGLE_mean
	PhS_var <- results_phenologynullmodelSINGLE_var
	Empirical <- results_empirical_normalized_nestedness
	delta_nestedness_norm_PhS_vec <- (PhS_mean - RG_mean)/(Empirical - RG_mean)
	outlier_points <- which(abs((delta_nestedness_norm_PhS_vec-mean(delta_nestedness_norm_PhS_vec))/(3*sd(delta_nestedness_norm_PhS_vec)))>1, arr.ind=TRUE)
	if (length(outlier_points)==0){
		outlier_points <- -1 * c(1:length(delta_nestedness_norm_PhS_vec))
	}
	RG_mean <- RG_mean[c(-1 * outlier_points)]
	RG_var <- RG_var[c(-1 * outlier_points)]
	DD_mean <- DD_mean[c(-1 * outlier_points)]
	DD_var <- DD_var[c(-1 * outlier_points)]
	PhS_mean <- PhS_mean[c(-1 * outlier_points)]
	PhS_var <- PhS_var[c(-1 * outlier_points)]
	Empirical <- Empirical[c(-1 * outlier_points)]
	# error propagation (in terms of variance until last step)
	delta_nestedness_norm_PhS_vec <- (PhS_mean - RG_mean)/(Empirical - RG_mean)
	delta_nestedness_norm_PhS <- mean(delta_nestedness_norm_PhS_vec)
	numerator_PhS <- PhS_mean - RG_mean
	numerator_error_PhS <- PhS_var + RG_var
	denominator_PhS <- Empirical - RG_mean
	denominator_error_PhS <- RG_var
	overall_error_PhS <- delta_nestedness_norm_PhS_vec^2 * ((numerator_error_PhS/numerator_PhS^2) + (denominator_error_PhS/denominator_PhS^2))
	overall_error_PhS <- mean(overall_error_PhS)
	overall_error_PhS <- sqrt(overall_error_PhS)	
} else {
	delta_nestedness_norm_DD_vec <- (DD_mean - RG_mean)/(Empirical - RG_mean)
	delta_nestedness_norm_DD <- mean(delta_nestedness_norm_DD_vec)
	delta_nestedness_norm_PhS_vec <- (PhS_mean - RG_mean)/(Empirical - RG_mean)
	delta_nestedness_norm_PhS <- mean(delta_nestedness_norm_PhS_vec)
	# error propagation (in terms of variance until last step)
	# DD
	numerator_DD <- DD_mean - RG_mean
	numerator_error_DD <- DD_var + RG_var
	denominator_DD <- Empirical - RG_mean
	denominator_error_DD <- RG_var
	overall_error_DD <- delta_nestedness_norm_DD_vec^2 * ((numerator_error_DD/numerator_DD^2) + (denominator_error_DD/denominator_DD^2))
	overall_error_DD <- mean(overall_error_DD)
	overall_error_DD <- sqrt(overall_error_DD)
	# PhS
	numerator_PhS <- PhS_mean - RG_mean
	numerator_error_PhS <- PhS_var + RG_var
	denominator_PhS <- Empirical - RG_mean
	denominator_error_PhS <- RG_var
	overall_error_PhS <- delta_nestedness_norm_PhS_vec^2 * ((numerator_error_PhS/numerator_PhS^2) + (denominator_error_PhS/denominator_PhS^2))
	overall_error_PhS <- mean(overall_error_PhS)
	overall_error_PhS <- sqrt(overall_error_PhS)
}

print("Nestedness results (Delta rho-tilde values):")
print("Degree distribution model -- Mean and stdev:")
print(delta_nestedness_norm_DD)
print(overall_error_DD)
#
print("Phenology model -- Mean and stdev:")
print(delta_nestedness_norm_PhS)
print(overall_error_PhS)
print("See Results.csv for fit metrics.")
####

####
# Generate results table
model_list <- c("RG", "DD", "Phen")
metrics_list <- c("Nestedness","Sen", "Spec", "Fscore", "Inform", "Phi")
headings_models <- c(rep(model_list[1], length(metrics_list)), rep(model_list[2], length(metrics_list)), rep(model_list[3], length(metrics_list)))
headings_metrics <- c(rep(metrics_list, length(model_list)))
metric_means <- rep(-1, length(model_list)*length(metrics_list))
metric_stdevs <- rep(-1, length(model_list)*length(metrics_list))
#
metric_means[1] <- 0
metric_means[2] <- mean(results_Bernoullinullmodel_sensitivity_mean)
metric_means[3] <- mean(results_Bernoullinullmodel_specificity_mean)
metric_means[4] <- mean(results_Bernoullinullmodel_Fscore_mean)
metric_means[5] <- mean(results_Bernoullinullmodel_Informedness_mean)
metric_means[6] <- mean(results_Bernoullinullmodel_Phi_mean)
#
metric_means[7] <- delta_nestedness_norm_DD
metric_means[8] <- mean(results_degreedistnullmodel_sensitivity_mean)
metric_means[9] <- mean(results_degreedistnullmodel_specificity_mean)
metric_means[10] <- mean(results_degreedistnullmodel_Fscore_mean)
metric_means[11] <- mean(results_degreedistnullmodel_Informedness_mean)
metric_means[12] <- mean(results_degreedistnullmodel_Phi_mean)
#
metric_means[13] <- delta_nestedness_norm_PhS
metric_means[14] <- mean(results_phenologynullmodelSINGLE_sensitivity_mean)
metric_means[15] <- mean(results_phenologynullmodelSINGLE_specificity_mean)
metric_means[16] <- mean(results_phenologynullmodelSINGLE_Fscore_mean)
metric_means[17] <- mean(results_phenologynullmodelSINGLE_Informedness_mean)
metric_means[18] <- mean(results_phenologynullmodelSINGLE_Phi_mean)
#
metric_stdevs[1] <- sqrt(mean(RG_var))
metric_stdevs[2] <- sqrt(mean(results_Bernoullinullmodel_sensitivity_var))
metric_stdevs[3] <- sqrt(mean(results_Bernoullinullmodel_specificity_var))
metric_stdevs[4] <- sqrt(mean(results_Bernoullinullmodel_Fscore_var))
metric_stdevs[5] <- sqrt(mean(results_Bernoullinullmodel_Informedness_var))
metric_stdevs[6] <- sqrt(mean(results_Bernoullinullmodel_Phi_var))
#
metric_stdevs[7] <- overall_error_DD
metric_stdevs[8] <- sqrt(mean(results_degreedistnullmodel_sensitivity_var))
metric_stdevs[9] <- sqrt(mean(results_degreedistnullmodel_specificity_var))
metric_stdevs[10] <- sqrt(mean(results_degreedistnullmodel_Fscore_var))
metric_stdevs[11] <- sqrt(mean(results_degreedistnullmodel_Informedness_var))
metric_stdevs[12] <- sqrt(mean(results_degreedistnullmodel_Phi_var))
#
metric_stdevs[13] <- overall_error_PhS
metric_stdevs[14] <- sqrt(mean(results_phenologynullmodelSINGLE_sensitivity_var))
metric_stdevs[15] <- sqrt(mean(results_phenologynullmodelSINGLE_specificity_var))
metric_stdevs[16] <- sqrt(mean(results_phenologynullmodelSINGLE_Fscore_var))
metric_stdevs[17] <- sqrt(mean(results_phenologynullmodelSINGLE_Informedness_var))
metric_stdevs[18] <- sqrt(mean(results_phenologynullmodelSINGLE_Phi_var))
#
myresults_table <- rbind(headings_models, headings_metrics, metric_means, metric_stdevs)
write.table(myresults_table, "Results.csv", sep=",", row.names=FALSE, col.names=FALSE)
####

# print runtime
print("TOTAL RUNTIME:")
print(Sys.time() - start_time)

