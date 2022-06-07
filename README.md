# High-Dimensional-Data-Analysis-HarvardX-

Week 1 :Distance
Many clustering and machine learning techniques rely on being able to define distance, using features or predictors. 
The euclidean distance between A and B is simply:

(Ax−Bx)2+(Ay−By)2−−−−−−−−−−−−−−−−−−−−√

Distance in High Dimensions
We introduce a dataset with gene expression measurements for 22,215 genes from 189 samples. The R objects can be downloaded like this:
   #library(devtools)
   #install_github("genomicsclass/tissuesGeneExpression")
The data represent RNA expression levels for eight tissues, each with several individuals.
    library(tissuesGeneExpression)
    data(tissuesGeneExpression)
    dim(e) ##e contains the expression data   
    #write.csv(e,"tissuesGeneExpression.csv")
    tissue 
    table(tissue) ##tissue[i] tells us what tissue is represented by e[,i]
    
    
Examples
We can now use the formulas above to compute distance. Let’s compute distance between samples 1 and 2, both kidneys, and then to sample 87, a colon.
  x <- e[,1]
  y <- e[,2]
  z <- e[,87]
  tissue[c(1,2,87)]                      ## [1] "kidney" "kidney" "colon"
  sqrt(sum((x-y)^2))                     
  
As expected, the kidneys are closer to each other. A faster way to compute this is using matrix algebra:
  sqrt( crossprod(x-y) )
  
Now to compute all the distances at once, we have the function dist. Because it computes the distance between each row, and here we are interested in the distance between samples, we transpose the matrix
  d <- dist(t(e))
  class(d)

Note that this produces an object of class dist and, to access the entries using row and column, indexes we need to coerce it into a matrix:
  as.matrix(d)[1,2]

Q:How many biological replicates for hippocampus?     
   sum(tissue=="hippocampus")
   
What is the distance between samples 3 and 45?        
   d <- dist(t(e))
   class(d)
   as.matrix(d)[3,45]     or     sqrt( crossprod(e[,3]-e[,45]) )        or        sqrt( sum((e[,3]-e[,45])^2 ))


What is the distance between gene 210486_at and 200805_at
   x<-e["210486_at",]
   y<-e["200805_at",]
   sqrt( crossprod(x-y) )        or    sqrt( crossprod(e["210486_at",]-e["200805_at",]) )     or    sqrt( sum((e["210486_at",]-e["200805_at",])^2 ))

Compute the distance between all pairs of samples: Q:How many distances are stored in d? (Hint: What is the length of d)?
   d = dist(t(e))
   length(d)

We will describe powerful techniques for exploratory data analysis based on dimension reduction.The general idea is to reduce the dataset to have fewer dimensions, yet approximately preserve important properties, such as the distance between samples.

Example: Reducing two dimensions to one

z1 = (y[1,]+y[2,])/2 #the sum 
z2 = (y[1,]-y[2,])   #the difference

z = rbind( z1, z2) #matrix now same dimensions as y

thelim <- c(-3,3)
mypar(1,2)

plot(y[1,],y[2,],xlab="Twin 1 (standardized height)",ylab="Twin 2 (standardized height)",xlim=thelim,ylim=thelim)
points(y[1,1:2],y[2,1:2],col=2,pch=16)

plot(z[1,],z[2,],xlim=thelim,ylim=thelim,xlab="Average height",ylab="Differnece in height")
points(z[1,1:2],z[2,1:2],col=2,pch=16)







A <- 1/sqrt(2)*matrix(c(1,1,1,-1),2,2)
z <- A%*%y
d <- dist(t(y))
d2 <- dist(t(z))
mypar(1,1)
plot(as.numeric(d),as.numeric(d2)) #as.numeric turns distnaces into long vector
abline(0,1,col=2)

We call this particular transformation a rotation of y.

mypar(1,2)

thelim <- c(-3,3)
plot(y[1,],y[2,],xlab="Twin 1 (standardized height)",ylab="Twin 2 (standardized height)",xlim=thelim,ylim=thelim)
points(y[1,1:2],y[2,1:2],col=2,pch=16)

plot(z[1,],z[2,],xlim=thelim,ylim=thelim,xlab="Average height",ylab="Differnece in height")
points(z[1,1:2],z[2,1:2],col=2,pch=16)

So this rotation actually achieves what we originally wanted: we can preserve the distances between points with just one dimension. Let’s remove the second dimension of z and recompute distances:

d3 = dist(z[1,]) ##distance computed using just first dimension
mypar(1,1)
plot(as.numeric(d),as.numeric(d3)) 
abline(0,1)

This first dimension of the transformed data is actually the first principal component

Week 2: Dimension Reduction
Projections

Simple example with N=2
If we let Y=(2 3). We can plot it like this:

mypar (1,1)
plot(c(0,4),c(0,4),xlab="Dimension 1",ylab="Dimension 2",type="n")
arrows(0,0,2,3,lwd=3)
text(2,3," Y",pos=4,cex=3)

Singular Value Decomposition (SVD)

Applying the SVD to the motivating example we have:

library(rafalib)
library(MASS)
n <- 100
y <- t(mvrnorm(n,c(0,0), matrix(c(1,0.95,0.95,1),2,2)))

mypar()
LIM<-c(-3.5,3.5)
plot(y[1,],y[2,],xlim = LIM,ylim = LIM)

s <- svd(y)
round(sqrt(2) * s$u , 3)

The plot we showed after the rotation, was showing what we call the principal components: the second plotted against the first. To obtain the principal components from the SVD, we simply need the columns of the rotation U⊤Y :

PC1 = s$d[1]*s$v[,1]
PC2 = s$d[2]*s$v[,2]
plot(PC1,PC2,xlim=c(-3,3),ylim=c(-3,3))


Compute the SVD of e:                                         s = svd(e)
Now compute the mean of each row:                             m = rowMeans(e)
What is the correlation between the first column of  and m?   cor(s$u[,1],m)

we saw how the first column relates to the mean of the rows of e.Note that if we change these means, the distances between columns do not change. Here is some R code showing how changing the means does not change the distances:

newmeans = rnorm(nrow(e)) ##random values we will add to create new means
newe = e+newmeans ##we change the means
sqrt(crossprod(e[,3]-e[,45]))
sqrt(crossprod(newe[,3]-newe[,45]))

So we might as well make the mean of each row 0 since it does not help us approximate the column distances. We will define y as the detrended e and recompute the SVD:

y = e - rowMeans(e)
s = svd(y)

We showed that  is equal to y up to numerical error:

resid = y - s$u %*% diag(s$d) %*% t(s$v)
max(abs(resid))

The above can be made more efficient in two ways. First, using the crossprod() and second not creating a diagonal matrix. Note that in R we can multiply a matrix x by vector a. The result is a matrix with row i equal to x[i,]*a[i]. Here is an example to illustrate this.

x=matrix(rep(c(1,2),each=5),5,2)
x
x*c(1:5)

Note that the above code is actually equivalent to:

sweep(x,1,1:5,"*")

This means that we don't have to convert s$d into a matrix to obtain .

Which of the following gives us the same as diag(s$d)%*%t(s$v)?
s$d * t(s$v)

Let z = s$d * t(s$v). We showed a derivation demonstrating that because  is orthogonal, the distance between e[,3] and e[,45] is the same as the distance between y[,3] and y[,45], which is the same as z[,3] and z[,45]:

z = s$d * t(s$v)
sqrt(crossprod(e[,3]-e[,45]))
sqrt(crossprod(y[,3]-y[,45]))
sqrt(crossprod(z[,3]-z[,45]))

Note that the columns z have 189 entries, compared to 22,215 for e.

What is the difference (in absolute value) between the actual distance sqrt(crossprod(e[,3]-e[,45])) and the approximation using only two dimensions of z?

realdistance = sqrt(crossprod(e[,3]-e[,45]))
approxdistance = sqrt(crossprod(z[1:2,3]-z[1:2,45]))
abs(realdistance - approxdistance)

What is the minimum number of dimensions we need to use for the approximation in SVD Exercises #4 to be within 10% or less?

ks = 1:189
realdistance = sqrt(crossprod(e[,3]-e[,45]))
approxdistances = sapply(ks,function(k){
    sqrt(crossprod(z[1:k,3,drop=FALSE]-z[1:k,45,drop=FALSE] )) 
  })
percentdiff = 100*abs(approxdistances - realdistance)/realdistance
plot(ks,percentdiff) ##take a look
min(ks[which(percentdiff < 10)])

Compute distances between sample 3 and all other samples:
distances = sqrt(apply(e[,-3]-e[,3],2,crossprod))

Recompute this distance using the 2 dimensional approximation.

What is the Spearman correlation between this approximate distance and the actual distance?

approxdistances = sqrt(apply(z[1:2,-3]-z[1:2,3],2,crossprod))
plot(distances,approxdistances) ##take a look
cor(distances,approxdistances,method="spearman")





