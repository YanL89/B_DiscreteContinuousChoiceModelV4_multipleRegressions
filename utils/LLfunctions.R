source("utils/deriv.R")

minProba = 10e-20

LLWrapper = function(b, f, args){
	p = f(b,args)
	p[p < minProba] = minProba
	sum(log(p))
	}
	
LLGradWrapper = function(b, f, args){
	if(args[["verbose"]])
		print(LLWrapper(b,f,args))
	delta = args[["delta"]]
	JMgrad(f = LLWrapper, b, c(list(f = f), list(args = args)), delta)
	}

LLHessianWrapper = function(b, f, args){
	JMhessian(LLWrapper, b, c(list(f = f), list(args = args))
			, args[["delta"]])
	}
