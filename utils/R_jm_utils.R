#' print the member of a list to the screen
#'
#' @param l a list
#' @param name the name of an element
#'
#' @output name : value if name is defined in l
printIfElem = function(l,name){
	if(name %in% names(l))
		cat(name, " : ", l[[name]], "\n")
	}

getElem = function(l, name, default){
	if(name %in% names(l))
		return(l[[name]])
	return(default)
	}

logicalAnd = function(b){
	length(b) == sum(b)
	}
	
matLogicalAnd = function(b){
	dim(b)[1] * dim(b)[2] == sum(b)
	}
