args <- commandArgs(trailingOnly = TRUE)

datadir <- args[1] 

files <- system(paste('ls', file.path(datadir, '"*.RData'), intern=T)

for ( file in file ){
	load( file );
	chrSums = Matrix::Matrix( chrSums, sparse=T)

	min=min(chrSums@x)
	if( min > 0 ){
		min = 0
	}
	max = max( chrSums@x )
	max = quantile(chrSums@x, .90)
	
	dim(chrSums)
	chrSums= chrSums[ goodGenes(chrSums, max),]
	dim(chrSums)
	
	
	rm ( chrSums );
	gc()
}
