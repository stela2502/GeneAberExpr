context( 'Class usage')

prefix= './';

load( file.path(prefix, 'data', 'testD.RData'))

## extremely simple test data first 100 cells - cancer next 100 cells healthy
min = min(testD)
max = quantile(testD@x, .90)
table(goodGenes(testD, max)) 
expect_equal( table(goodGenes(testD, max)), table( c(rep( FALSE, 37) , rep(TRUE, 2694) )))

testD= testD[ goodGenes(testD, max),]

model = GetTestModel( testD, seq( min, max,(max - min)/9 ) ,as.integer(1:100), 
	as.integer(101:200), phony=FALSE)


getSimple <- function(x, range, at ){
	x = x[at]
	ret = vector(mode="numeric", length=length(range) )
	for ( i in 2:(length(range)) ){
		rem = which( x < range[i] )
		ret[i-1] = length(rem);
		if ( length(rem) > 0 ){
			x = x[-rem]
		}
	}
	ret[length(range)] = length(x)
	ret
}

getVals <- function(x, range, at){
	ret = getSimple(x, range, at)

	m = min(ret)

	ret = ((ret - m) / sum(ret)) + 1e-9
	log(ret)
}

for( i in  seq(2, nrow(model),2) ) {
	expect_equal ( 
		as.vector( model[i,]), 
		getVals( as.vector(testD[i/2, ]),seq( min, max,(max - min)/9 ),101:200),
		info= paste('bad bg model row', i) )
}


for( i in  seq(1, nrow(model),2) ) {
	expect_equal ( 
		as.vector( model[i,]), 
		getVals( 
			as.vector(testD[ceiling(i/2), ]),
			seq( min, max,(max - min)/9 ),
			1:100
		), info= paste('bad int model row', i) )
}


calcCancerCount2 <- function(x) {
	 length(which(x > .9999))
}
expect_true ( round(sum(calcQuality (model) ),0) ==  1819581)

#[1] 0.005600354 #max
#[1] 0.01155574  #quant 99
#    0.04714544  #quant 80
# 80% quantile of expressed - that is better!

#plot(1:10, res[1,], col='red', t='l')
#lines(res[2,], col='green')

res1 = IdentifyStatesTest( testD, seq( min, max,(max - min)/9 ) ,
	as.integer(1:100), as.integer(101:200), phony=FALSE)

load(file.path( prefix,'data','red.RData'))

expect_equal( res1, res, label = 'HMM real live results')

## now the thing is at least not breaking all the time!
## look for logics problems!






