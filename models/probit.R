# probit model using "pmvnorm" approximation for the multivariate normal CDF
#
# Jean-Michel Tremblay
# University of Maryland
# jeanmi.tremblay@gmail.com

library("mvtnorm")

source("utils/create_X.R")
source("utils/deriv.R")
source("utils/probitUtils.R")

source("models/logit.R")

probit = list(
	LLVec = function(x,args){
		X = args[["X"]]
		Y = args[["Y"]]
		
		beta = x[1:ncol(X)]
		#print(beta)
		L = vec2mat(c(sqrt(2), x[(ncol(X)+1):length(x)]))
		S = L %*% t(L)
	
		nalt = nrow(S) + 1
		n = length(Y)
	
		U = matrix(X %*% as.matrix(beta), nrow = n, ncol = nalt, byrow = T)
		
		p = rep(0,length(Y))
		for(i in 1:n){
			p[i] = probitL(list(U = U[i,], choice = Y[i], SCond = S
						,muCond = rep(0,nalt-1), nalt = nalt, method = args$method))
			}
		p
		},

	computeArgs = function(spec,D){
		comPlusSpec = list(common = spec$common, specific = spec$specific)
		X = getDiscreteMatrix(comPlusSpec, D)
		#print(X)
		Y = D[[spec$Y]]
		Y = Y - min(Y) + 1
		delta = getElem(spec, name = "delta", default = 0.1)
		args = list(X = X, Y = Y, delta = delta, method = spec$method)
		
		
		nalt = length(spec$specific)
		n = nrow(D)
		# if we simulate the probas, generate the error terms here
		if(spec$method == "simulation"){
			args[["Z"]] = list()
			for(i in 1:n)
				args[["Z"]][[i]] = matrix(rnorm(spec$nsim * (nalt-1))
						, (nalt -1), spec$nsim)
			}
		args
		},

	computeStart = function(spec, D){
		# estimate the same specification with logit and use that
		# for starting values. that saves a few iterations
		L = model(logit, spec, D)
		beta = L$results$beta_hat / (pi / sqrt(6))
		
		nalt = length(spec$specific)
		# if errors have covariance identity
		# then the differences have covariance identity + 1
		S = diag(nalt-1) + 1
		# the first element is set, we don't estimate it
		c(beta, mat2vec(S)[-1])
		},

	computeOther = function(spec, D){
		comPlusSpec = list(common = spec$common, specific = spec$specific)
		nalt =  length(spec$specific)
		X = getDiscreteMatrix(comPlusSpec, D)
		var_names = c(colnames(X))
		for(i in 2:(nalt-1))
			for(j in 1:i)
				var_names = c(var_names, paste("L_",i,j,sep=""))
		list(names = var_names)
		}	
	)


















probit_LL = function(beta, X, Y, delta = 0.1){
	betaPr = beta[1:ncol(X)]
	L = vec2mat(c(1,beta[(ncol(X)+1):length(beta)]))
	S = L %*% t(L)
	
	nalt = nrow(S) + 1
	n = length(Y)
	
	U = matrix(X %*% as.matrix(betaPr), nrow = n, ncol = nalt, byrow = T)
	
	LL = 0
	for(i in 1:n){
		M = reparamErrors(nalt-1, currentBase = 1, newBase = Y[i])
		SReparam = M %*% S %*% t(M)
		muReparam = rep(0,nalt-1) # cause M times 0 = 0
		
		UDiff = U[i,Y[i]] - U[i, -Y[i]]
		
		LL = LL + log(pmvnorm(upper = UDiff, mean = muReparam, sigma = SReparam))
		}
	LL
	}

# example use
if(FALSE){
	# D is given somewhere else, just make sure that all the variables are in
	common = list(c("X1","X3","X5","X7","X9"))
	specific = list(c("X2"), c("X4"), c("X6"), c("X8"), c("X10"))
	YPrString = "YPr"
	P = probit(common, specific, YPrString, D, delta = 0.1)	
	}
	
