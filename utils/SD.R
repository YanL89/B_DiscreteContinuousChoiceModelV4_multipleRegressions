source("utils/deriv.R")

NoneSD = function(beta_hat){
	rep(9999999, length(beta_hat))
	}

HessianSD = function(beta_hat, modelFns, args){
	H = LLHessianWrapper(beta_hat, modelFns$LLVec, args)
	sqrt(diag(solve(-H)))	
	}	

# stuff for bootstrap variance estimation
LLWrapperBoot = function(b, f, args, w){
	p = f(b,args)
	p[p < minProba] = minProba
	# w tell how many times we should count each obs
	sum(w * log(p))
	}
	
LLGradWrapperBoot = function(b, f, args, w){
	# print(LLWrapper(b,f,args))
  
	delta = args[["delta"]]
	JMgrad(f = LLWrapperBoot, b, c(list(f = f), list(args = args)
			, list(w = w)), delta)
	}

#BootstrapSD = function(b, args, modelFns, spec, start, lb, ub){
#	n = length(modelFns$LLVec(b,args))
#	nboot = spec$nboot
#	all_estimates = matrix(0,length(b), nboot)
#	for(i in 1:nboot){
#		#if(0 == i %% 10)
#			cat("bootstrap iteration ",i,"\n")
#		bootsamp = sample(n,n,replace=TRUE)
#		w = rep(0,n)
#		for(j in bootsamp)
#			w[j] = w[j] + 1
#		start_boot = b
#		O = optim(fn = LLWrapperBoot, gr = LLGradWrapperBoot
#				, method = "BFGS", hessian =F, par = start_boot
#				, control = list(fnscale = -1, reltol = spec[["reltol"]])
#				, args = args, w = w, f = modelFns$LLVec)
#		all_estimates[,i] = O$par
#		
#		}
#	#print(t(all_estimates))
#	sqrt(diag(var(t(all_estimates))))
#	}


BootstrapSD = function(b, args, modelFns, spec, start, lb, ub){
	n = length(modelFns$LLVec(b,args))
	nboot = spec$nboot
	all_estimates = matrix(0,length(b), nboot)
	for(i in 1:nboot){
		#if(0 == i %% 10)
			cat("bootstrap iteration ",i,"\n")
		bootsamp = sample(n,n,replace=TRUE)
		w = rep(0,n)
		for(j in bootsamp)
			w[j] = w[j] + 1
		start_boot = b
            O = optim(fn = LLWrapperBoot, gr = LLGradWrapperBoot, method = "BFGS"
			, hessian = F, par = start, control = list(fnscale = -1
			, reltol = spec$reltol), args = args, w = w, f = modelFns$LLVec)
		all_estimates[,i] = O$par
		
		}
	#print(t(all_estimates))
	sqrt(diag(var(t(all_estimates))))
	}

	
