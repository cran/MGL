#' Module network inference
#' 
#' Takes a high-dimensional data matrix, initial values of the module latent variables, and a penalty parameter,
#' and returns the final assignment of the data points to the modules, the values of the module latent variables, 
#' and the conditional dependency network among the module latent variables.
#' @param data An nxp matrix which contains n samples from p variables, where typically p>>n
#' @param L An nxk matrix which contains the initial latent variable values, a column for each module
#' @param lambda A penalty parameter controlling the sparsity of the conditional dependency network among the modules
#' @param printoutput 1 if the user wants the output from each iteration to be displayed, 0 for silent run
#' @param maxiter Maximum number of iterations to be performed
#' @param threshold Threshold for convergence
#' @return \item{L}{An nxk matrix which contains the final latent variable values, a column for each module}
#' @return \item{theta}{A kxk symmetric positive-semidefinite matrix respresenting the conditional dependency network among the modules}
#' @return \item{Z}{A p-vector containing values between 1 to k, representing the assignment of the p variables to k modules}
#' @useDynLib MGL, .registration = TRUE
#' @examples
#' \dontrun{
#' library(MGL)
#' n = 20 #sample size
#' p = 100 #variable size
#' k = 5 #module size
#' lambda = .01 #penalty parameter to induce sparsity
#' data = matrix(rnorm(n*p), ncol=p)
#' # to start with initial random module latent variables
#' L = matrix(rnorm(n*k), ncol=k)
#' MGL(data, L, lambda)
#' # to start with k-means cluster centroids as module latent variables
#' L = t(kmeans(t(data), k)$centers)
#' MGL(data, L, lambda)
#' }
#' @export

MGL = function(data, L, lambda, printoutput=0, maxiter=100, threshold=1e-2){
	sz = nrow(data)
	p = ncol(data)
	mcnt = ncol(L)
	thetaout = numeric(mcnt*mcnt)
	Zout = numeric(p)
	# below two lines are because of the fact that R passes matrices to C column by column
	data = t(data)
	L = t(L)
	storage.mode(data) = 'double'
	storage.mode(L) = 'double'
	storage.mode(thetaout) = 'double'
	storage.mode(Zout) = 'integer'
	res = .C('MGL', data, Lout=L, as.integer(sz), as.integer(p), as.integer(mcnt), as.double(lambda), as.integer(maxiter), as.double(threshold), as.integer(printoutput), thetaout=thetaout, Zout=Zout)
	return(list('L'=t(matrix(res$Lout, nrow=mcnt, ncol=sz)), 'theta'=matrix(res$thetaout, nrow=mcnt, ncol=mcnt), 'Z'=t(res$Zout+1)))
}
