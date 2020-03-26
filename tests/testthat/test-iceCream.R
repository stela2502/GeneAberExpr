context( 'IceCream' )

prefix= './';

res = IceCreamTest()

exp=read.delim( file.path( prefix, 'data', 'Eisner_run1.csv'))
expect_equal(as.numeric(exp(res[,1])), as.numeric(as.vector(exp[,1])), label= colnames(exp)[1] )
expect_equal(as.numeric(exp(res[,5])), as.numeric(as.vector(exp[,4])), label= colnames(exp)[4] )
expect_equal(as.numeric(exp(res[,2])), as.numeric(as.vector(exp[,2])), label= colnames(exp)[2] )
expect_equal(as.numeric(exp(res[,6])), as.numeric(as.vector(exp[,5])), label= colnames(exp)[5] )
expect_equal(as.numeric(exp(res[,3])), as.numeric(as.vector(exp[,3])), label= colnames(exp)[3] )
expect_equal(as.numeric(exp(res[,7])), as.numeric(as.vector(exp[,6])), label= colnames(exp)[7] )
expect_equal(round(as.numeric(exp(res[,4])),3), as.numeric(as.vector(exp[,8])), label= colnames(exp)[8] )
expect_equal(round(as.numeric(exp(res[,8])),3), as.numeric(as.vector(exp[,9])), label= colnames(exp)[9] )

## OK that is all that is required here - I am good! 