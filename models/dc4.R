# Dicrete-Continuous model version 4
# fits a probit and multiple regressions with all error terms
# correlated


library(mvtnorm)

source("utils/create_X.R")
source("utils/probitUtils.R")

source("models/logit.R")
#source("models/probit.R")
#source("models/dc4.R")

dc4 = list(
  LLVec = function(b, args){
    params = dc4Paramsnew(b, nparPr = ncol(args$XPr)
                           , nparReg_Primary = ncol(args$XReg_Primary), 
                           nparReg_Secondary = ncol(args$XReg_Secondary),
                           nparReg_Third = ncol(args$XReg_Third))
    L = params$L
    YReg_Primary = args$YReg_Primary
    YReg_Secondary = args$YReg_Secondary
    YReg_Third = args$YReg_Third
    YPr = args$YPr
    S = L %*% t(L)
    
    #observed YReg (specific numbers are used to select corresponding rows in the data structure, need to be generalized for multiple regressions)
    YReg = matrix(0,1289,3)
    YReg[108:1289, 1] = YReg_Primary
    YReg[438:1289, 2] = YReg_Secondary
    YReg[1033:1289, 3] = YReg_Third
    
    # estimated YReg_hat
    YReg_Primary_hat = args$XReg_Primary %*% as.matrix(params$betaReg_Primary)
    YReg_Secondary_hat = args$XReg_Secondary %*% as.matrix(params$betaReg_Secondary)
    YReg_Third_hat = args$XReg_Third %*% as.matrix(params$betaReg_Third)
    YReg_hat = matrix(0,1289,3)
    YReg_hat[108:1289, 1] = YReg_Primary_hat
    YReg_hat[438:1289, 2] = YReg_Secondary_hat
    YReg_hat[1033:1289, 3] = YReg_Third_hat
    
    LL=0
    nalt= 4
    nDiff = nalt-1
    nobs = length(args$YPr)
    LLReg=rep(0,nobs)
    LLPr=rep(0,nobs)
    U = matrix(args$XPr %*% params$betaPr, nrow = args$n
               , ncol = args$nalt, byrow = TRUE)
    
    for(i in 1:nobs){
      # LL for the regression
      numOfY = YPr[i]-1
      SReg = S[(nDiff + 1):(nDiff + numOfY),(nDiff + 1):(nDiff + numOfY)]
      preg = -1
      if(numOfY > 1)
        preg = dmvnorm(YReg[i,1:numOfY],YReg_hat[i,1:numOfY],SReg) 
      if(1 == numOfY)
        preg = dnorm(YReg[i,1],YReg_hat[i,1],SReg)
      if(0 == numOfY)
        preg = 1
      LLReg[i] = preg
      
      # compute conditional dist of probit
      nKeep = nDiff + numOfY
      Schop = S[1:nKeep, 1:nKeep]
      Spr = S[1:nDiff, 1:nDiff]
      S21 = as.matrix(S[(nDiff + 1):(nDiff + numOfY),1:nDiff])
      if(ncol(S21) == 1)
        S21 = t(S21)
      S12 = t(S21)
      
      res = as.matrix(YReg[i,1:numOfY] - YReg_hat[i,1:numOfY])   #epsilon
      mu1 = as.matrix(rep(0,nDiff))
      mu2 = as.matrix(rep(0,length(res)))
      
      muCond = mu1
      Scond = Schop
      
      if(numOfY > 0){  
        muCond = mu1 + S12 %*% solve(SReg) %*% (res - mu2)
        Scond = Spr - S12 %*% solve(SReg) %*% S21
      }
      
      # choice start at 1
      argsLL = list(nalt = 4, choice = numOfY +1, mu = muCond, S = Scond
                    , U = U[i,])
      Ppr = pmvProbitL(argsLL)
      LLPr[i]=Ppr
    }
    LL = LLReg * LLPr

  },


#dc4 = list(
	#LLVec = function(b, args){
		#params = dc4Params(b, nparPr = ncol(args$XPr)
				#, nparReg = ncol(args$XReg))
	#	L = params$L
	#	YReg = args$YReg
		#YPr = args$YPr
		#S = L %*% t(L)
		
		# LL for the regression
		#epsilon = YReg - (args$XReg %*% as.matrix(params$betaReg))
	#	sigma2 = S[nrow(S), ncol(S)]
		#LReg = dnorm(epsilon, 0, sigma2)
		
		# utilities
		#U = matrix(args$XPr %*% params$betaPr, nrow = args$n
		#		, ncol = args$nalt, byrow = TRUE)

		# we simulate for the LL of the probit, 
		# we need a loop because at every line we need to reparametrize
		#LPr = rep(0, length(args$YPr))
	#	for(i in 1:length(args$YPr)){
		#	SigmaDiffCond = condCov(S, args$nalt)
		#	muDiffCond = condMean(S, mu = rep(0, args$nalt)
			#		, posObs = args$nalt, obs = epsilon[i])
		#	p = 0
			#S = SigmaDiffCond
			#muCond = muDiffCond
			#Uline = U[i,]
		#	probitArgs = list(nalt = args$nalt, choice = YPr[i], mu = muDiffCond
			#		, S = SigmaDiffCond, U = U[i,]) 
			#print(args)
			#if(args$method == "simulation")
			#	p = simProbitL(c(list(Z = args$Z[[i]]), probitArgs))
			#if(args$method != "simulation")
				#p = probitL(probitArgs)
		#	LPr[i] = p
		#	}
		#LReg * LPr
	#	},
	
	computeArgs = function(spec, D){
	  XReg_Primary = as.matrix(D[108:1289, spec$reg])
		XReg_Secondary = as.matrix(D[1290:2141, spec$reg])
		XReg_Third = as.matrix(D[2142:2398, spec$reg])
		YReg_Primary = D[108:1289, spec$YReg ]
    YReg_Secondary = D[1290:2141, spec$YReg ]
		YReg_Third = D[2142:2398, spec$YReg ]
    XPr = getDiscreteMatrix(spec$probit, D[1:1289,])   # "long" type data for probit model  1420*5=7100rows  attributes-25columns 
		YPr = D[1:1289, spec$YPr]
		YPr = YPr - min(YPr) + 1            #	transfer 0-4 to 1-5
		nalt = length(spec$probit$specific) # number of alternatives
		n = nrow(D[1:1289,])                         # number of observations
		args = list(XReg_Primary = XReg_Primary, YReg_Primary = YReg_Primary, XReg_Secondary = XReg_Secondary, YReg_Secondary = YReg_Secondary, 
                XReg_Third = XReg_Third, YReg_Third = YReg_Third, XPr = XPr, YPr = YPr, nalt = nalt, n = n, method = spec$method)

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
		nalt = length(spec$probit$specific)   # number of alternatives
		L = model(logit, c(list("Y" = spec$YPr), spec$probit), D[1:1289,])
		formReg = as.formula(paste(spec$YReg, paste(spec$reg,collapse="+")
				, sep = " ~ 0 +"))
    D_Primary = D[108:1289,]
    D_Secondary = D[1290:2141,]
		D_Third = D[2142:2398,]
		R_Primary = lm(formReg, D_Primary)
    R_Secondary = lm(formReg, D_Secondary)
		R_Third = lm(formReg, D_Third)
		betaPr = L$results$beta_hat / (pi / sqrt(6))
		betaReg_Primary = R_Primary$coefficients
		betaReg_Secondary = R_Secondary$coefficients
		betaReg_Third = R_Third$coefficients
		sigmaReg_Primary = mean(R_Primary$residuals ** 2)     # (sigma-reg)^2
		sigmaReg_Secondary = mean(R_Secondary$residuals ** 2)
		sigmaReg_Third = mean(R_Third$residuals ** 2)
    sigma2 = as.matrix(c(sigmaReg_Primary, sigmaReg_Secondary, sigmaReg_Third))      #sigma22
		nReg = length(sigma2)
    
		# if errors of utilities ~ N(0,I) then differences ~ N(0, I + 1)
    ndim = nalt-1+nReg
		S = diag(ndim)
		S[ndim-2, ndim-2] = sigmaReg_Primary 
		S[ndim-1, ndim-1] = sigmaReg_Secondary
		S[ndim, ndim] = sigmaReg_Third
		S[-ndim, -ndim] = S[-ndim, -ndim] + 1     # different covariance matrix
		L = t(chol(S))         #decholesky 
		# we fix the first element of L so we remove it
		c(betaPr, betaReg_Primary, betaReg_Secondary, betaReg_Third, mat2vec(L)[-1])
		},
		
	computeOther = function(spec, D)	{
    nReg = 3
		nalt = length(spec$probit$specific)
    ndim = nalt-1+nReg
	
		# record the names of coefficients
		namesPr = getNames(spec$probit, D)
		namesReg_Primary = spec$reg
		namesReg_Secondary = spec$reg
		namesReg_Third = spec$reg
		namesSigma = c()
		for(i in 2:ndim)
			for(j in 1:i)
				namesSigma = c(namesSigma,paste("L_",i,j,sep=""))
		names = c(namesPr, namesReg_Primary, namesReg_Secondary, namesReg_Third, namesSigma)
	  list(names = names)
		},
	
	findStartL = function(m){
		i = 1
		while(m$results$name[i] != "L_21")
			i = i + 1
		i
		},	
		
	reparam = function(m){
		startL = dc4$findStartL(m)
		endL = length(m$results$beta_hat)
		L = vec2mat(c(sqrt(2), m$results$beta_hat[startL:endL]))
		S = L %*% t(L)
		m2 = m
		m2$results = m2$results[1:(startL -1),]
		m2$SigmaDiffFirst = S
		class(m2) = "dc4"
		m2
		},
		
	comparePmvSim = function(m, spec, D, output){
		# calculate LL for Pmv
		specPmv = spec
		specPmv[["method"]] = "pmvnorm"
		
		specSim = spec
		specSim[["method"]] = "simulation"
		if(! "nsim" %in% names(specSim))
			specSim[["nsim"]] = 500
			
		argsPmv = dc4$computeArgs(specPmv, D)
		argsSim = dc4$computeArgs(specSim, D)
		
		pSim = dc4$LLVec(m$results$beta_hat, argsSim)
		pPmv = dc4$LLVec(m$results$beta_hat, argsPmv)
		
		pdf(output)
		plot(pSim, pPmv, xlab = "probas computed with simulation"
				, ylab = "probas computed with pmv"
				, main = "comparision of invididual probabilities"
				)
		abline(lm(pPmv ~ 0 + pSim))
		dev.off()
		},
		
	outputUtilities = function(m, spec, D, output){
		#b = as.vector(m$results$beta_hat)
		args = dc4$computeArgs(spec,D)
		#print(args$X)
		nparPr = length(spec$probit$common)
		for(i in spec$probit$specific)
			nparPr = nparPr + length(i)
			
		betaPr = m$results$beta_hat[1:nparPr]
		U = getUtilities(args$XPr, betaPr, nrow(D))
		U = round(1000*U) / 1000
		p = dc4$LLVec(m$results$beta_hat, args)
		p = round(1000*p) / 1000
		choice = D[,spec[["YPr"]]]
		
		out = data.frame(U = U, choice = choice, p = p)
		write.csv(out, output, quote = FALSE, row.names = FALSE)
		},
		
	comparePis = function(b1, b2, spec, D, output){
		args = dc4$computeArgs(spec, D)
		p1 = dc4$LLVec(b1, args)
		p2 = dc4$LLVec(b2, args)
		pdf(output)
		plot(p1,p2,xlab="p_i with first set of coefficients",
				ylab="p_i with second set", main="comparision of 2 cx sets")
		abline(lm(p2 ~ p1))
		dev.off()
		},
	
	apply = function(b, spec, D){
		XPr = create_X(spec$probit$common, spec$probit$specific, D)	
		XReg = as.matrix(D[,reg])
		params = dc4Params(b,ncol(XPr),ncol(XReg))	
		nalt = length(spec$probit$specific)
		nerr = nalt -1 +1
		n = nrow(D)

		# error terms
		err = params$L %*% matrix(rnorm(n*nerr),nrow = nerr, ncol = n)
		errsDiff = t(err[1:(nalt-1),])
		errsReg = err[nalt,]

		# calculate dependant varialbe
		b = params$betaPr
		YPr = genChoiceProbit(spec$probit, params$betaPr, D, errsDiff)
		YReg = XReg %*% params$betaReg + errsReg
		list(Ydisc = YPr, Yreg = YReg)
		}
	
	)
