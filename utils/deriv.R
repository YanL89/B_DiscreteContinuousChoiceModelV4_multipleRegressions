JMcx = c(1,-1,-1,1)
JMstep = c(1,-1)

JMmixedDeriv = function(fn, x, pos, method.args = list(), delta = 0.01){
	if(1 == length(delta))
		delta = rep(delta, length(x))
	h = 0
	index = 1
	for(p1 in JMstep)
		for(p2 in JMstep){
			xdiff = x
			xdiff[pos[2]] = xdiff[pos[2]] + p2 * delta[pos[2]]
			xdiff[pos[1]] = xdiff[pos[1]] + p1 * delta[pos[1]]
			h = h + JMcx[index] * do.call(fn,c(list(xdiff),method.args))
			index = index + 1
			}
	h = h / (4 * delta[pos[1]] * delta[pos[2]])		
	h
	}

JMderiv = function(fn, x, pos, method.args = list(), delta = 0.01){
	if(1 == length(delta))
		delta = rep(delta, length(x))	
	xmin = x
	xmax = x
	xmin[pos] = x[pos] - delta[pos]
	xmax[pos] = x[pos] + delta[pos]
	
	as.numeric((do.call(fn, c(list(xmax),method.args)) 
			- do.call(fn, c(list(xmin),method.args))) / (2*delta[pos]))
	}

JMgrad = function(fn, x, method.args = list(), delta = 0.01){
	n = length(x)
	g = rep(0,n)
	for(i in 1:n)
		g[i] = JMderiv(fn,x,i,method.args,delta)
	g
	}
	
JMhessian = function(fn, x, method.args = list(), delta = 0.01){
	n = length(x)
	H = matrix(0,n,n)
	for(i in 1:n)
		for(j in 1:i)					
			H[i,j] = H[j,i] = JMmixedDeriv(fn,x,c(i,j),method.args,delta)
	H
	}
