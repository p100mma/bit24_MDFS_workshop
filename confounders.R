#' p-value of mutual information of one discretized variable x and discrete y as a function of breaks
#'
#' @param x numeric vector, a continuous variable to be discretized
#' @param y numeric vector, discrete variable
#' @param breaks vector of division points of \code{x} to test
#' @return a chi-squared p-value of the mutual information 
#' @examples
#' x<- runif(1000)
#' pv.1D(x, (x > 0.5)*1., c(0.2,0.5,0.7) )
#' @export
pv.1D <- function(x,
                  y,
                  breaks) {
 
 y.classes <- unique(y)
 breaks <- c(-Inf,breaks,Inf)
 
 kx <- length(breaks)-1
 ky <- length(y.classes)
 
 #counters
 nxy<-matrix(0,kx,ky) 
 for (i in 1:kx) for (j in 1:ky) nxy[i,j] <- sum(x >= breaks[i] & x < breaks[i+1] & y == y.classes[j])
 nx_ <- rowSums(nxy)
 n_y <- colSums(nxy)
 n <- sum(nx_)
 
 #optimal pseudocounts
 px <- (nx_+1)/(n+2)
 py <- (n_y+1)/(n+2)
 pxy <- px%*%t(py)
 
 #auxiliary matrices
 qxy <- (1-px)%*%t(1-py)
 q2xy <- (1-px^2)%*%t(1-py^2)
 
 #coefficients
 c.est <- (sum(q2xy/pxy)-2*(kx-1)*(ky-1))/12/n
 b.est <- sum(qxy/pxy)/2/n
 
 pc <- c.est/b.est
 
 #Mutual information
 nxy <- nxy+pc
 nx_ <- rowSums(nxy)
 n_y <- colSums(nxy)
 n <- sum(nx_)
 
 IG <- sum(nxy*log(nxy))-sum(nx_*log(nx_))-sum(n_y*log(n_y))+n*log(n)
 
 return(pchisq(2*IG, df=(kx-1)*(ky-1), lower.tail=F))  
}


#' Discretize continuous variable in an optimal way according to p-value of mutual information with discrete y
#'
#' @param x numeric vector, a continuous variable to be discretized
#' @param y numeric vector, discrete variable
#' @param divisions number of divisions of \code{x} (n_classes -1)
#' @param range discretization range
#' @param grid number of random tries per each class of discretized \code{x}
#' @param seed random seed for randomized discretization procedure
#' @return a \code{list} with the following fields: 
#'  \itemize{
#'    \item \code{breaks.opt} -- division points of \code{x} in optimal discretization
#'    \item \code{pv.opt} -- chi-squared p-value of optimal discretization
#'    \item \code{x} -- discretized \code{x}
#'  }
#' @examples
#' x<- runif(1000)
#' opt.discretize(x, (x>0.5)*1.)
#' @export
opt.discretize <- function(x, y, divisions=1, range=0.5, grid=10, seed=NULL) {

 if (!is.null(seed)) set.seed(seed)
 npt<-grid^(divisions+1)
 
 pv.opt <- 2
 for (i in 1:npt) {
  
  breaks <- cumsum((1-range)+range*runif(divisions+1))
  breaks <- breaks/max(breaks)*length(x)
  breaks <- (x[order(x)][ceiling(breaks)])[1:divisions]
  pv <- pv.1D(x, y, breaks)
  
  if (pv < pv.opt) {
   breaks.opt <- breaks
   pv.opt <- pv
  }
 }

 breaks <- c(-Inf,breaks.opt,Inf)
 x.discr <- rep(NA, length(x))
 for (i in 1:(divisions+1)) x.discr[x >= breaks[i] & x < breaks[i+1]] <- i

 return(list(breaks.opt=breaks.opt, pv.opt=pv.opt, x=x.discr))
}


#' Pack many confounders into a single value
#'
#' @param confounders object containing discrete variables, gets coerced to \code{data.frame}
#' @param all.combinations boolean stating whether to include all possible combinations in the packed confounder (\code{TRUE}) or only those appearing in data provided (\code{FALSE})
#' @return a vector containing discrete variable, a packed confounder with one value per each combination of \code{confounders}
#' @examples
#' pack.confounders(madelon$data[,1:3], all.combinations=TRUE)
#' @export 
pack.confounders <- function(confounders, all.combinations=F) {
 confounders<-as.data.frame(confounders)
 
 if (length(unique(confounders)) > sqrt(nrow(confounders))) stop('Confounders must be discrete!')
 
 #transform many confounders into a single one
 if (!all.combinations)  combinations <- unique(confounders) else {
  uniques<-list()
  for (i in 1:ncol(confounders)) uniques[[length(uniques)+1]]<-unique(confounders[,i])
  combinations<-expand.grid(uniques)   
 } 
 confounder.classes <- nrow(combinations)
 confounder <- rep(NA, nrow(confounders))
 for (i in 1:confounder.classes) confounder[apply(confounders, 1, function(x){ all(x==combinations[i,]) })] <- i
 
 return(confounder)
} 



############################################################################################
################################### Feature selection ######################################
############################################################################################


##########################################################################
#continuous data
##########################################################################

#' Run end-to-end MDFS
#'
#' @details
#' In case of FDR control it is recommended to use Benjamini-Hochberg-Yekutieli p-value adjustment
#' method (\code{"BY"} in \code{\link[stats]{p.adjust}}) due to unknown dependencies between tests.
#'
#' @param data input data where columns are variables and rows are observations (all numeric)
#' @param decision decision variable as a boolean vector of length equal to number of observations
#' @param confounders object containing discrete variables, gets coerced to \code{data.frame}
#' @param n.contrast number of constrast variables (defaults to max of 1/10 of variables number and 30)
#' @param dimensions number of dimensions (a positive integer; on CUDA limited to 2--5 range)
#' @param divisions number of divisions (from 1 to 15)
#' @param discretizations number of discretizations
#' @param range discretization range (from 0.0 to 1.0; \code{NULL} selects probable optimal number)
#' @param pc.xi parameter xi used to compute pseudocounts (the default is recommended not to be changed)
#' @param p.adjust.method method as accepted by \code{\link[stats]{p.adjust}} (\code{"BY"} is recommended for FDR, see Details)
#' @param level statistical significance level
#' @param seed seed for PRNG used during discretizations (\code{NULL} for random)
#' @param use.CUDA whether to use CUDA acceleration (must be compiled with CUDA; NOTE: the CUDA version might provide a slightly lower sensitivity due to a lack of native support for \code{contrast_data})
#' @return A \code{\link{list}} with the following fields:
#'  \itemize{
#'    \item \code{contrast.indices} -- indices of variables chosen to build contrast variables
#'    \item \code{contrast.variables} -- built contrast variables
#'    \item \code{MIG.Result} -- result of ComputeMaxInfoGains
#'    \item \code{MDFS} -- result of ComputePValue (the MDFS object)
#'    \item \code{statistic} -- vector of statistic's values (IGs) for corresponding variables
#'    \item \code{p.value} -- vector of p-values for corresponding variables
#'    \item \code{adjusted.p.value} -- vector of adjusted p-values for corresponding variables
#'    \item \code{relevant.variables} -- vector of relevant variables indices
#'  }
#' @examples
#' \donttest{
#' MDFS.confounded(madelon$data[,-c(1:3)], madelon$decision, confounders=madelon$data[, 1:3], dimensions = 2, divisions = 1,
#'      range = 0, seed = 0)
#' }
#' @importFrom stats p.adjust
#' @export
MDFS.confounded<-function(data, 
                          decision, 
                          confounders, 
                          n.contrast = max(ncol(data), 30),
                          dimensions = 1,
                          divisions = 1,
                          discretizations = 1,
                          range = NULL,
                          pc.xi = 0.25,
                          p.adjust.method = "holm",
                          level = 0.05,
                          seed = NULL,
                          use.CUDA = FALSE
) {                          
 data <- as.data.frame(data)
 if (!is.null(seed)) set.seed(seed)

 if ((dimensions > 1 || discretizations > 1) && (is.null(n.contrast) || n.contrast < 30)) stop('Specify the number of contrast variables!')

 #transform many confounders into a single one
 confounder <- pack.confounders(confounders)
 confounder.classes <- unique(confounder)

 if (length(confounder) != nrow(data)) stop('Confounders must have the same number of rows as data!')

 if (any(table(data.frame(cbind(decision,confounder)))==0)) stop ('Confounders fully explain the decision!')

 #contrast variables
 if (n.contrast > 0) {
  contrast <- GenContrastVariables(data, n.contrast)
  contrast.indices <- contrast$indices
  contrast_data <- contrast$contrast_data
  contrast.mask <- c(rep.int(F, ncol(data)), rep.int(T, ncol(contrast_data)))
  IG <- rep(0, ncol(data)+ncol(contrast_data))
 }
 else {
  contrast.mask <- contrast.indices <- contrast_data <- NULL
  IG <- rep(0, ncol(data))
 }

 #compute IG
 for (cc in confounder.classes) {
  
  MIG.Result <- ComputeMaxInfoGains(data=data[confounder==cc,],
                                    decision=decision[confounder==cc],
                                    contrast_data=contrast_data[confounder==cc,],
                                    dimensions = dimensions,
                                    divisions = divisions,
                                    discretizations = discretizations,
                                    seed = seed,
                                    range = range,
                                    pc.xi = pc.xi, 
                                    return.tuples = F,
                                    use.CUDA=use.CUDA)

  IG <- IG+c(MIG.Result$IG, attr(MIG.Result, "contrast_igs"))
 }
 
 #p-value
 df <- divisions*(divisions+1)^(dimensions-1)*length(confounder.classes)
 PV <- ComputePValue(IG, dimensions, divisions, discretizations, response.divisions=1, df=df, contrast.mask=contrast.mask)
 
 #result
 result <- list(contrast.indices = contrast.indices,
                statistic = PV$statistic,
                p.value = PV$p.value,
                adjusted.p.value = p.adjust(PV$p.value,method=p.adjust.method),
                relevant.variables = which(p.adjust(PV$p.value,method=p.adjust.method)<level))

 return(result)
}



##########################################################################
#discrete data with the same number of categories
##########################################################################
#' Run end-to-end MDFS for discretized data with same number of categories and discrete confounding factors
#'
#' @details
#' In case of FDR control it is recommended to use Benjamini-Hochberg-Yekutieli p-value adjustment
#' method (\code{"BY"} in \code{\link[stats]{p.adjust}}) due to unknown dependencies between tests.
#'
#' @param data input data.frame where columns are variables and rows are observations (all numeric)
#' @param decision decision variable as a boolean vector of length equal to number of observations
#' @param confounders object containing discrete variables, gets coerced to \code{data.frame}
#' @param n.contrast number of constrast variables (defaults to max of 1/10 of variables number and 30)
#' @param dimensions number of dimensions 
#' @param pc.xi parameter xi used to compute pseudocounts (the default is recommended not to be changed)
#' @param p.adjust.method method as accepted by \code{\link[stats]{p.adjust}} (\code{"BY"} is recommended for FDR, see Details)
#' @param level statistical significance level
#' @param seed seed for PRNG used during discretizations (\code{NULL} for random)
#' @return A \code{\link{list}} with the following fields:
#'  \itemize{
#'    \item \code{contrast.indices} -- indices of variables chosen to build contrast variables
#'    \item \code{statistic} -- vector of statistic's values (IGs) for corresponding variables
#'    \item \code{p.value} -- vector of p-values for corresponding variables
#'    \item \code{adjusted.p.value} -- vector of adjusted p-values for corresponding variables
#'    \item \code{relevant.variables} -- vector of relevant variables indices
#'  }
#' @examples
#' \donttest{
#' MDFS.discrete(madelon$data[,-c(1:3)], madelon$decision,confounders=madelon$data[,1:3], dimensions = 2)
#' }
#' @importFrom stats p.adjust
#' @export
MDFS.discrete.confounded<-function(data, 
                                   decision, 
                                   confounders, 
                                   n.contrast = max(ncol(data), 30),
                                   dimensions = 1,
                                   pc.xi = 0.25,
                                   p.adjust.method = "holm",
                                   level = 0.05,
                                   seed = NULL
) {                          
 data <- as.data.frame(data)
 if (!is.null(seed)) set.seed(seed)
 
 if (dimensions > 1  && (is.null(n.contrast) || n.contrast < 30)) stop('Specify the number of contrast variables!')

 #transform many confounders into a single one
 confounder <- pack.confounders(confounders)
 confounder.classes <- unique(confounder)
 
 if (length(confounder) != nrow(data)) stop('Confounders must have the same number of rows as data!')
 
 if (any(table(data.frame(cbind(decision,confounder)))==0)) stop ('Confounders fully explain the decision!')

 #contrast variables
 if (n.contrast > 0) {
  contrast <- GenContrastVariables(data, n.contrast)
  contrast.indices <- contrast$indices
  contrast_data <- contrast$contrast_data
  contrast.mask <- c(rep.int(F, ncol(data)), rep.int(T, ncol(contrast_data)))
  IG <- rep(0, ncol(data)+ncol(contrast_data))
 }
 else {
  contrast.mask <- contrast.indices <- contrast_data <- NULL
  IG <- rep(0, ncol(data))
 }
 
 #compute IG
 for (cc in confounder.classes) {
  
  MIG.Result <- ComputeMaxInfoGainsDiscrete(data=data[confounder==cc,],
                                    decision=decision[confounder==cc],
                                    contrast_data=contrast_data[confounder==cc,],
                                    dimensions = dimensions,
                                    pc.xi = pc.xi)
  
  IG <- IG+c(MIG.Result$IG, attr(MIG.Result, "contrast_igs"))
 }
 
 #p-value
 divisions <- length(unique(data[,1]))-1
 df <- divisions*(divisions+1)^(dimensions-1)*length(confounder.classes)
 PV <- ComputePValue(IG, dimensions, divisions, discretizations=1, response.divisions=1, df=df, contrast.mask=contrast.mask)
 
 #result
 result <- list(contrast.indices = contrast.indices,
                statistic = PV$statistic,
                p.value = PV$p.value,
                adjusted.p.value = p.adjust(PV$p.value,method=p.adjust.method),
                relevant.variables = which(p.adjust(PV$p.value,method=p.adjust.method)<level))
 
 return(result)
}


############################################################################################
###################################### Similarity ##########################################
############################################################################################


##########################################################################
#discrete data with the same number of categories
##########################################################################

# PS: this function will find indexes of all discrete variables 'V' in data.frame 'data'
#     for which there exists a value of 'confounder' at which some class of 'V' does not appear.

#' Find indexes of variables fully explained by a confounding factor
#'
#' @param data \code{data.frame} of discrete variables
#' @param confounder discrete confounding variable
#' @return a vector of indexes of variables in \code{data} which are fully explained by \code{confounder} 
#' @examples
#' find_fully_explained_vars_idx(madelon$data, madelon$data[,2])
#' @export
find_fully_explained_vars_idx<- function(data,confounder){
which(unlist(lapply(data, function(x) any(table(data.frame(cbind(x,confounder)))==0) )))
}

#' Compute similarity matrix of discrete variables with same number of categories based on their mutual information taking into the account confounding variables
#'
#' @details
#' In case of FDR control it is recommended to use Benjamini-Hochberg-Yekutieli p-value adjustment
#' method (\code{"BY"} in \code{\link[stats]{p.adjust}}) due to unknown dependencies between tests.
#' If any variable in \code{data} is fully explained by any confounder, execution stops.
#' To find fully explained variables one can use \code{find_fully_explained_vars_idx} function.
#'
#' @param data input data.frame where columns are variables and rows are observations (all numeric, discrete variables with same number of categories)
#' @param confounders object containing discrete variables, gets coerced to \code{data.frame}
#' @param divisions number of classes of each variable (same for all variables in data)
#' @param p.adjust.method method as accepted by \code{\link[stats]{p.adjust}} (\code{"BY"} is recommended for FDR, see Details)
#' @param level statistical significance level
#' @return A \code{\link{list}} of symmetric matrices, with the following fields:
#'  \itemize{
#'    \item \code{S} -- similarity matrix of variables where \code{S[i,j]} is equal to \code{MI[i,j]} if \code{pv.adjusted[i,j]} is \code{<level} else 0
#'    \item \code{pv.adjusted} -- \code{pv.adjusted[i,j]} is adjusted p-value of mutual information \code{MI[i,j]}
#'    \item \code{pv} -- \code{pv[i,j]} is p-value of mutual information \code{MI[i,j]} 
#'    \item \code{MI} -- \code{MI[i,j]} is mutual information computed between variables \code{data[,i]} and \code{data[,j]}
#'  }
#' @examples
#' \donttest{
#' pairwiseMIsimilarityDiscrete.confounded(madelon$data[,-1], madelon$data[,1])
#' }
#' @importFrom stats p.adjust
#' @export
pairwiseMIsimilarityDiscrete.confounded <- function(data, 
                                                    confounders,
                                                    divisions=1, 
                                                    p.adjust.method="holm", 
                                                    level=0.05
){
 data <- as.data.frame(data)
 
 #transform many confounders into a single one
 confounder <- pack.confounders(confounders)
 confounder.classes <- length(unique(confounder))
 
 if (length(confounder) != nrow(data)) stop('Confounders must have the same number of rows as data!')
 n_fully_explained<- length(find_fully_explained_vars_idx(data, confounder)	) 
 if ( n_fully_explained ) stop (sprintf('Confounders fully explain %d of variables in data!',n_fully_explained) )

 #compute IG
 MI <- matrix(0, ncol(data), ncol(data))
 for (cc in confounder.classes) {
  MI <- MI+ComputeInterestingTuplesDiscrete(data=data[confounder==cc,],
                                            dimensions=2,
                                            return.matrix=TRUE,
                                            stat_mode="MI")*sum(confounder==cc)
 }
 diag(MI) <- 0

 divisions <- length(unique(data[,1]))-1
 df <- divisions^2*length(confounder.classes)
 pv<- MI
 #p-value
 pv[,] <- ComputePValue(IG=as.vector(MI), 
                     dimensions=1, 
                     divisions=divisions, 
                     response.divisions=divisions,
                     discretizations=1,
                     df=df)$p.value	    
 diag(pv) <- 1

 #adjust p-value
 ltr_selector <- lower.tri(pv)
 adjustedFlatPV <- p.adjust(pv[ltr_selector], method=p.adjust.method)

 pv.adj <- matrix(nrow=ncol(data), ncol=ncol(data))
 pv.adj[ltr_selector] <- adjustedFlatPV
 pv.adj <- t(pv.adj)
 pv.adj[ltr_selector] <- adjustedFlatPV
 pv.adj <- t(pv.adj)

 diag(pv.adj) <- 1
 
 #similarity
 S <- MI
 S[pv.adj >= level] <-0
 
 return (list(S=S,
              pv.adjusted=pv.adj,
              pv=pv,
              MI=MI))
}


##########################################################################
#continuous variables
##########################################################################


#' Compute similarity matrix of continuous data based on their mutual information after discretization, taking into the account the influence of confounding variables
#'
#' @details
#' In case of FDR control it is recommended to use Benjamini-Hochberg-Yekutieli p-value adjustment
#' method (\code{"BY"} in \code{\link[stats]{p.adjust}}) due to unknown dependencies between tests.
#'
#' @param data input data.frame where columns are variables and rows are observations (all numeric, discrete variables with same number of categories)
#' @param confounders object containing discrete variables, gets coerced to \code{data.frame}
#' @param divisions number of divisions
#' @param discretizations number of discretizations
#' @param n.contrast number of contrast variables
#' @param p.adjust.method method as accepted by \code{\link[stats]{p.adjust}} (\code{"BY"} is recommended for FDR, see Details)
#' @param level statistical significance level
#' @param seed random seed for randomized discretization procedure
#' @return A \code{\link{list}} of symmetric matrices, with the following fields:
#'  \itemize{
#'    \item \code{S} -- similarity matrix of variables where \code{S[i,j]} is equal to \code{MI[i,j]} if \code{pv.adjusted[i,j]} is \code{<level} else 0
#'    \item \code{pv.adjusted} -- \code{pv.adjusted[i,j]} is adjusted p-value of mutual information \code{MI[i,j]}. 
#'    \item \code{pv} -- \code{pv[i,j]} is p-value of mutual information \code{MI[i,j]} 
#'    \item \code{MI} -- \code{MI[i,j]} is mutual information computed between variables \code{data[,i]} and \code{data[,j]}
#'  }
#' @examples
#' \donttest{
#' pairwiseMIsimilarity.confounded(madelon$data[,-1], madelon$data[,1])
#' }
#' @importFrom stats p.adjust
#' @export
pairwiseMIsimilarity.confounded <- function(data, 
                                            confounders,
                                            divisions=1,
                                            discretizations=1,
                                            n.contrast = max(ncol(data), 30),
                                            p.adjust.method="holm", 
                                            level=0.05,
                                            seed=NULL){

 if (discretizations > 1 && (is.null(n.contrast) || n.contrast < 30)) stop('Specify the number of contrast variables!')
 
 data <- as.data.frame(data)
 if (!is.null(seed)) set.seed(seed)
 
 
 #transform many confounders into a single one
 confounder <- pack.confounders(confounders)
 confounder.classes <- length(unique(confounder))
 
 if (length(confounder) != nrow(data)) stop('Confounders must have the same number of rows as data!')

 #contrast variables
 if (n.contrast > 0) {
  contrast <- GenContrastVariables(data, n.contrast)
  contrast.indices <- contrast$indices
  contrast_data <- contrast$contrast_data
 }
 else {
  contrast.mask <- contrast.indices <- contrast_data <- NULL
 }

 #compute IG
 MI <-  matrix(0, ncol(data), ncol(data))
 if (n.contrast>0) MI.contrast <- matrix(0, n.contrast, n.contrast) else MI.contrast<-as.matrix(0)
 for (cc in confounder.classes) {
  MI <- MI+ComputeInterestingTuples(data=data[confounder==cc,],
                                    dimensions=2,
                                    divisions=divisions,
                                    discretizations=discretizations,
                                    return.matrix=TRUE,
                                    stat_mode="MI")*sum(confounder==cc)
  
  if (n.contrast > 0)   MI.contrast <- MI.contrast+ComputeInterestingTuples(data=contrast_data[confounder==cc,],
                                                                            dimensions=2,
                                                                            divisions=divisions,
                                                                            discretizations=discretizations,
                                                                            return.matrix=TRUE,
                                                                            stat_mode="MI")*sum(confounder==cc)
                    
    
 }
 diag(MI) <- diag(MI.contrast) <- 0

 df <- divisions^2*length(confounder.classes)
 ltr_selector <- lower.tri(MI)
 ltr_selector.contrast <- lower.tri(MI.contrast) 
 
 flatMI <- c(MI[ltr_selector],MI.contrast[ltr_selector.contrast])
 contrast.mask <- c(rep.int(F, sum(ltr_selector)), rep.int(T, sum(ltr_selector.contrast)))
 
 
 #p-value
 flatPV <- ComputePValue(IG=flatMI, 
                         dimensions=1, 
                         divisions=divisions, 
                         response.divisions=divisions,
                         discretizations=discretizations,
                         df=df,
                         contrast.mask=contrast.mask)$p.value	    
 
 adjustedFlatPV <- p.adjust(flatPV[!contrast.mask], method=p.adjust.method)
 
 pv <- pv.adj <- matrix(nrow=ncol(data), ncol=ncol(data))
 pv[ltr_selector] <- flatPV[!contrast.mask]
 pv <- t(pv)
 pv[ltr_selector] <- flatPV[!contrast.mask]
 pv <- t(pv)
 
 pv.adj[ltr_selector] <- adjustedFlatPV
 pv.adj <- t(pv.adj)
 pv.adj[ltr_selector] <- adjustedFlatPV
 pv.adj <- t(pv.adj)
 
 diag(pv) <- diag(pv.adj) <- 1
 
 #similarity
 S <- MI
 S[pv.adj >= level] <-0
 
 return (list(S=S,
              pv.adjusted=pv.adj,
              pv=pv,
              MI=MI))
}


 
