#' ---
#' title: Section 4. Plotting Expression values
#' author: 
#' output: github_document
#' ---


#' This is a complement to the hands on session for lecture 13.
#' First we need to load our expression data from a file (that you downloaded):

#expr <- read.table("rs8067378_ENSG00000172057.6.txt")

#' Or we can read from online.
url <- "https://bioboot.github.io/biggn213_S18/class-material/rs8067378_ENSG00000172057.6.txt"
expr <- read.table(url)

#' Check the genotype sample size and overall summary stats 
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

# Boxplot with the data shown
ggplot(expr, aes(geno, exp, fill=geno)) + 
  geom_boxplot(notch=TRUE, outlier.shape = NA) + 
  geom_jitter(shape=16, position=position_jitter(0.2), alpha=0.4)

