
#' Print output for statistical (fitted) model
#'
#' @param m a model with attributes $results, a data frame with first column
#'	the name of the estimates and last 4 columns the beta_hat, SD, t and p.
print.dcModel = function(m){
	r = m$results
	r[,2:5] = round(r[,2:5] * 1000) / 1000
	print(r)

	printIfElem(m, "maxLL")
	printIfElem(m, "LLAtIntercept")
	printIfElem(m, "LLAtZero")
	}

print.dc4 = function(m){
	print.dcModel(m)
	print("Variance of difference (with first alt) of error terms:")
	print(m$SigmaDiffFirst)
	}
