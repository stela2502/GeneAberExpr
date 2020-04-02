context( 'Models reliable calling' )


prefix= './';

load( file.path(prefix, 'data', 'testD.RData'))

min = min(testD)
max = quantile(testD@x, .90)

## this would almost 100% break of a free(): invalid size fault.
cmpModel =  GetTestModel( testD, seq( min, max,(max - min)/9 ) ,as.integer(1:100), as.integer(101:200), phony=TRUE)
for ( i in 1:100){
	model = GetTestModel( testD, seq( min, max,(max - min)/9 ) ,as.integer(1:100), as.integer(101:200), phony=FALSE)
	expect_equal( cmpModel ,model , tolerance = 1e-6 ) # we are still alive!
	if ( ! all(cmpModel== model )){
		break;
	}
}