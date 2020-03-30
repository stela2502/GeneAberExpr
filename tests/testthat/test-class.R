context( 'Class usage')

prefix= './';

load( file.path(prefix, 'data', 'testD.RData'))

## extremely simple test data first 100 cells - cancer next 100 cells healthy
min = min(testD)
max = quantile(testD@x, .90)
table(goodGenes(testD, max)) 
expect_equal( table(goodGenes(testD, max)), table( c(rep( FALSE, 37) , rep(TRUE, 2694) )))

testD= testD[ goodGenes(testD, max),]

model = GetTestModel( testD, seq( min, max,(max - min)/9 ) ,as.integer(1:100), as.integer(101:200))


calcCancerCount2 <- function(x) {
	 length(which(x > .9999))
}
expect_true ( round(sum(calcQuality (model) ),4) ==  268.7359)

#[1] 0.005600354 #max
#[1] 0.01155574  #quant 99
#    0.04714544  #quant 80
# 80% quantile of expressed - that is better!

#plot(1:10, res[1,], col='red', t='l')
#lines(res[2,], col='green')

res1 = IdentifyStatesTest( testD, seq( min, max,(max - min)/9 ) ,as.integer(1:100), as.integer(101:200))

load(file.path( prefix,'data','red.RData'))

expect_equal( res1, res, label = 'HMM real live results')

## now the thing is at least not breaking all the time!
## look for logics problems!

skip("not test relevant")

res2 = IdentifyStates( testD, seq( min, max,(max - min)/9 ) ,as.integer(1:100), as.integer(101:200))
cancerCount = calcCancerCount(res2 )

boxplot( split( cancerCount, c(rep('cancer', 100), rep('healthy', 100))), main="Full Model")

quality =calcQuality( model )

res3 = IdentifyStates( testD[ sort(order(quality,decreasing=T)[1:800]),], seq( min, max,(max - min)/9 ) ,as.integer(1:100), as.integer(101:200))

cancerCount2 = calcCancerCount( res3 )

x11()
boxplot( split( cancerCount2, c(rep('cancer', 100), rep('healthy', 100))), main="Partial Model")

## the partial model it is then :-(




