
context( 'real world example')

if ( ! file.exists( '~/NAS/MattiasMagnusson/Pavan_Prabhala/Seq8/chr3.chrSums.RData' ) ){
	skip('example file not found')
}

skip('too time consuming and outdated')

calcQuality <- function(id){
	 sum((model[id,] - model[id+1,])^2)
}
## this is the better function!
calcCancerCount <- function(x) {
	# ( paste( 'calcCancerCount',length(which(x > .99)), "-", (length(which(x < 0.01))+1)))
	 length(which(x > .99))- length(which(x < 0.01))
}
calcCancerCount2 <- function(x) {
	 length(which(x > .9999))
}

load( '~/NAS/MattiasMagnusson/Pavan_Prabhala/Seq8/chr3.chrSums.RData')
chrSums = Matrix::Matrix( chrSums, sparse=T)
##chrSums@x = log( chrSums@x )
min=min(chrSums@x)
if( min > 0 ){
	min = 0
}
max = max( chrSums@x )
max = quantile(chrSums@x, .90)
dim(chrSums)
chrSums= chrSums[ goodGenes(chrSums, max),]
dim(chrSums)

sname = unlist( stringr::str_replace_all( colnames(chrSums), '_[AGCT]*$', ''))
set.seed(5)
normal = sample( grep( 'NovaSeq6.Pavan1', sname ), 300 )
cancer = sample( grep( 'NovaSeq6.Pavan3_Cancer3', sname ), 300 )

model = GetTestModel( chrSums, seq( min, max,(max - min)/9 ), cancer, normal )
quality = unlist( lapply( seq( 1,nrow(model),2), calcQuality ) )
sum(quality)
#[1] 227.0871

res = IdentifyStates( chrSums[ sort(order(quality,decreasing=T)[1:800]),] , seq( min, max,(max - min)/9 ) , cancer, normal )


cancerCount = apply( res,2, calcCancerCount )
cancerCount.2 = apply( res,2, calcCancerCount2 )


boxplot( split( cancerCount,sname ), border=c(rep('black',13), 'green', 'black','red',rep('black',2)), cex.axis=0.5, las=2)
abline(h=0, col='red')
x11()
boxplot( split( cancerCount.2,sname ), border=c(rep('black',13), 'green', 'black','red',rep('black',2)), cex.axis=0.5, las=2)




#load( '~/NAS/MattiasMagnusson/Pavan_Prabhala/Seq8/chr3.chrSums.RData')
#chrSums = Matrix::Matrix( chrSums, sparse=T)

chrSums@x = log(chrSums@x)
image(as.matrix(Matrix::t( chrSums[,order(cancerCount)])), col=gplots::bluered(40))
image(as.matrix(Matrix::t( chrSums[,order(cancerCount.2)])), col=gplots::bluered(40))

## I think this tool is fooled by the high expression samples - again!

# normal_reaccessed <-  which( cancerCount3 <= quantile(cancerCount3, .2 ))
# if ( length(normal_reaccessed) > 500 ) {
# 	normal_reaccessed = sample( normal_reaccessed, 500)
# }
# cancer_reaccessed <-  which( cancerCount3 >=  quantile(cancerCount3, .8 ))
# if ( length(cancer_reaccessed) > 500 ) {
# 	cancer_reaccessed = sample( cancer_reaccessed, 500)
# }

# load( '~/NAS/MattiasMagnusson/Pavan_Prabhala/Seq8/chr3.chrSums.RData')
# chrSums = Matrix::Matrix( chrSums, sparse=T)
# min=0
# max = quantile(chrSums@x, .80)
# model = GetTestModel( chrSums, seq( min, max,(max - min)/9 ) ,cancer_reaccessed, normal_reaccessed )
# quality = unlist( lapply( seq( 1,nrow(model),2),calcQuality ) )
# sum(quality)
# # [1] 23.05016
# res5 = IdentifyStates( chrSums, seq( min, max,(max - min)/9 ) , cancer_reaccessed, normal_reaccessed )
# cancerCount5 = apply( res5,2, calcCancerCount )

# boxplot( split( cancerCount5,sname ), border=c(rep('black',13), 'green', 'black','red',rep('black',2)), cex.axis=0.5, las=2)



# chrSums =chrSums [sort(order(quality,decreasing=T)[1:800]),]
# min=0
# max = quantile(chrSums@x, .80)

# res4 = IdentifyStates( chrSums, seq( min, max,(max - min)/9 ) , cancer_reaccessed, normal_reaccessed )

# cancerCount4 = apply( res4,2, calcCancerCount )

# boxplot( split( cancerCount4,sname ), border=c(rep('black',13), 'green', 'black','red',rep('black',2)), cex.axis=0.5, las=2)

# boxplot( split( cancerCount4,sname ))
# ## reaccess the cancer and healthy c