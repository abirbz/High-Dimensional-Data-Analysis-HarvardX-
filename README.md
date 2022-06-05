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
