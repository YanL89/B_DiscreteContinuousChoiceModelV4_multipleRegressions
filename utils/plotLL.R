plotLL = function(plotrange, numpoints, b, output, f, args, names){
	margin = plotrange / 2
	np = numpoints
	
	pdf(output)
	for(i in 1:length(b)){
		b2 = b
		LL = c()
		all_bmods = seq(from = b[i] - margin, to = b[i] + margin, length = np)
		for(bmod in all_bmods){
			b2[i] = bmod
			LL = c(LL, LLWrapper(b2, f, args))
			}
		plot(all_bmods, LL, type = "l", xlab = names[i], ylab = "LL")
		}
	dev.off()	
	}
