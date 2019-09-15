plotLongCART<- function(x, uniform = FALSE, branch = 1, nspace,
     margin = 0, minbranch = 0.3, ...){
	plot(x=x, uniform = uniform, branch = branch, compress = FALSE,
           nspace=nspace, margin = margin, minbranch = minbranch, ...)
}


textLongCART<- function(x, splits = TRUE, all = FALSE, use.n = FALSE, 
                minlength = 1L, ...){
text(x=x, splits = splits , all = all, use.n = use.n, minlength = minlength, ...)
}