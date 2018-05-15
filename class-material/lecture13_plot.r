#' ---
#' output: github_document
#' ---

#' # Plotting Expresion Values
#' 
#' This is a complement to the hands on session for lecture 13.
#' First we need to load our expresion data from a file (that you downloaded) or from online.

#expr <- read.table("rs8067378_ENSG00000172057.6.txt")
url <- "https://bioboot.github.io/bggn213_S18/class-material/rs8067378_ENSG00000172057.6.txt"
expr <- read.table(url)

#' Check the genotype sample size and overal summary stats 
summary(expr)

#' Lets break it down to the three genotypes by 
#' 
summary( expr[expr[,2] == "A/A",3] )
summary( expr[expr[,2] == "G/G",3] )
summary( expr[expr[,2] == "A/G",3] )


#' We could of also used `expr$geno == "A/A"` in the above, e.g.
summary( expr$exp[ expr$geno == "A/A" ] )

#' Lets generate some nice plots of our results

# Generate a boxplot for expression data of different genotype groups
p <- boxplot(exp~geno, data=expr, xlab="rs8067378 genotype", ylab="ENSG00000172057.4 (RPKM)", notch=T)

#' Look at the output 'p' and you will see the median values etc. used to generate the plot
p

#' ## Optional: Using ggplot2

#' We could also use the **ggplot2** package to make our plots


library(ggplot2)

## Boxplot
ggplot(expr, aes(geno, exp)) + geom_boxplot()

## Histogram of the exp column with ggplot2
ggplot(expr, aes(exp, fill = geno)) + geom_density(alpha = 0.2)
