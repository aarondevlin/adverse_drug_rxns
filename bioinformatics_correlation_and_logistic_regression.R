# Predicting Adverse Drug Reactions
# The files “gene_expression_n438x978.txt” and “ADRs_HLGT_n438x232.txt” must be in the working directory
# These exercises are adapted from from the NIH LINCS DCIC Crowdsourcing Portal and Ma’ayan Lab @ Mt Sinai, New York.
# http://www.maayanlab.net/crowdsourcing/megatask1.php


# Correlation and data exploration

# Read in the gene expression data and place it in a matrix, M.
mydata <- read.table("gene_expression_n438x978.txt", sep="\t", header=T)
dimnames(mydata)[[1]] <- mydata[,1]
mydata <- mydata[,-1]
t <- t(mydata)
drugcolm <- as.data.frame.matrix(t)
M <- matrix(nrow = 437, ncol = 438, mode(numeric))

# Compute the Correlation Coefficient between all pairs of drugs
for (i in 1:ncol(drugcolm)) {
  cor <- apply (drugcolm[, -i], 2, function(x) {cor(drugcolm[, i], x)})
  M[,i] <- cor
}

# Turn the matrix into a vector in order to plot a histogram of the correlations
nm <- as.vector(M, mode = "numeric")
print("Correlation Coefficients Between All Pairs of Drugs")
print (hist(nm, xlab = "Correlation Coefficient", main = "Correlation of All Drug Pairs"))
print (c("The median correlation coefficient is", median(nm)))


# Name the drug pair that gives the highest correlation coefficient (i.e.closest to 1) and
# produce a scatter plot to show the relationship between this pair of drugs using R

mydata<-read.table("gene_expression_n438x978.txt", sep="\t", header=F)
rnam <- as.vector(mydata[,1])  # set row names
rnam <- rnam[-1]
mydata <- mydata[-1, -1]
mydata <- as.matrix(mydata)
mode(mydata) <- "numeric"

# Create a matrix, df, to store correlation values
df <- matrix(nrow=438, ncol=438)
mode(df)<-"numeric"
colnames(df) <- rnam
rownames(df) <- rnam

# Finds the correlation coefficients between all pairs of drugs
for (i in 1:438){
  thecor <- apply (mydata, 1, function(x) {cor(mydata[i, ], x)})
  df[i,] <- thecor
}

# Set all diagonals (drugs that are correlated with themselves, with values = 1) to 0 so they wont be seen by min and max
for (j in 1:438){
  df[j,j] <- 0
}

print("Question 1B")
print("Drug Pair with Highest Correlation Value:")
print(max(df))
print("Drug Pair:")

# Find the row and column names with the value number
rcnum<-which.max(df)
ro<-floor(rcnum/438)
rmdr<-(rcnum/438)-ro
co<-round(rmdr*438)
ro<-ro+1
print(colnames(df)[co])
print(rownames(df)[ro])

plot(mydata[co,], mydata[ro,], main = "Plot of Drug Pair with Highest Correlation Value", xlab = colnames(df)[co], ylab = rownames(df)[ro])

# Name 10 drug pairs that give the top 10 highest correlation coefficients (i.e. closest to 1)
dfc<- df
mode(dfc)<-"numeric"
dimnames(dfc)[[1]] <- rnam
colnames(dfc) <- rnam

# The matrix is made up of two equal triangles(halves) separated by the diagonal 1's
# here we remove one of those duplicate halves of the matrix
for (j in 1:438) {
  for (i in 1: 438){
    if(dfc[i,j] == dfc[j,i])
      dfc[j,i]<- 0
  }
}

# Turn the matrix into a vector in order to find the top 10 by sorting the vector
tempc <- as.vector (unlist(dfc), mode="numeric")
tempd<-sort(tempc, decreasing = TRUE)
print("Question 1C")
print("The 10 Drug Pairs with Highest Correlation:")
for (i in 1:10) {    # use inds variable to find the row and col names for each of the top 10 correlated drugs pairs
  inds = which(dfc == tempd[i], arr.ind=TRUE)
  inds
  rnames = rownames(dfc)[inds[,1]]
  cnames = colnames(dfc)[inds[,2]]
  print(c(rnames, "+", cnames))
}


# Name the drug pair that gives the lowest correlation (i.e. closest to -1) and make a scatter plot
# to show the relationship

mydata<-read.table("gene_expression_n438x978.txt", sep="\t", header=F)
rnam <- as.vector(mydata[,1])  #row names
rnam <- rnam[-1]
mydata <- mydata[-1, -1]
mydata <- as.matrix(mydata)
mode(mydata) <- "numeric"

df <- matrix(nrow=438, ncol=438)
mode(df)<-"numeric"
colnames(df) <- rnam
rownames(df) <- rnam

for (i in 1:438){
  thecor <- apply (mydata, 1, function(x) {cor(mydata[i, ], x)})
  df[i,] <- thecor
}

for (j in 1:438){
  df[j,j] <- 0
}

print("Minimum Drug Pair Correlation Value:")
print(min(df))
print("The Minimum Correlated Drug Pair:")
rcnum<-which.min(df)
ro<-floor(rcnum/438)
rmdr<-(rcnum/438)-ro
co<-round(rmdr*438)
ro<-ro+1
print(colnames(df)[co])
print(rownames(df)[ro])
plot(mydata[co,], mydata[ro,], main = "Plot of Minimum Correlated Drug Pair", xlab = colnames(df)[co], ylab = rownames(df)[ro])


# Find the drug that is most similar to each of the following drugs:

# CLOFARABINE
mydata<-read.table("gene_expression_n438x978.txt", sep="\t", header=T)
dimnames(mydata)[[1]]<-mydata[,1]
mydata<-mydata[,-1]
t <- t(mydata)
drugcolm <- as.data.frame.matrix(t)
M <- matrix(nrow = 438, ncol = 2, mode(numeric))
dimnames(M)[[1]] <- rnam

# The most similar drug to CLOFARABINE (396)
for (i in 1:438) {
  thecor<- cor(drugcolm[, 396], drugcolm[, i])
  M[i,1] <- thecor
  M[i,2] <- i
}

# Turn drug correlation with itself to 0 so we dont see it in the max
M[396,1]<-0
maxnum<-max(M[,1]) 
maxrow<-which(M == maxnum)
print(c("The most similar drug to CLOFARABINE is:", rownames(M)[maxrow]))

# Most similar drug to DAUNORUBICIN (382)
for (i in 1:438) {
  thecor<- cor(drugcolm[, 382], drugcolm[, i])
  M[i,1] <- thecor
  M[i,2] <- i
}

M[382,1]<-0
maxnum<-max(M[,1]) 
maxrow<-which(M == maxnum)
print(c("The most similar drug to DAUNORUBICIN is:", rownames(M)[maxrow]))

# Most similar drug to FLUDARABINE (211)
for (i in 1:438) {
  thecor<- cor(drugcolm[, 211], drugcolm[, i])
  M[i,1] <- thecor
  M[i,2] <- i
}

M[211,1]<-0
maxnum<-max(M[,1]) 
maxrow<-which(M == maxnum)
print(c("The most similar drug to FLUDARABINE is:", rownames(M)[maxrow]))


# Using both txt files and considering the first 50 genes only, we will use forward 
# stepwise logistic regression to build predictive models for each of the 232 side effects (adr).
# This is done for each of the 438 drugs.
# Which adr can be predicted with the smallest # of errors? Smallest AIC?

# Loading both tables
mydata <- read.table("gene_expression_n438x978.txt", sep="\t", header=T)
rownames(mydata) <- mydata[,1]
mydata <- mydata[,-1]
adrmat <- read.table ("ADRs_HLGT_n438x232.txt", sep="\t", header=T)
rownames(adrmat) <- adrmat[,1]
adrmat <- adrmat[,-1]

# Subset of first 50 genes from the drug and gene data
mydatasub1 <- data.frame(mydata[, 1:50])
rownames(mydatasub1) <- rownames(adrmat)

# Storage for the sum of errors & AIC
numErrors <- matrix(ncol = 2, nrow = 232)
rownames(numErrors) <- colnames(adrmat)
colnames(numErrors) <- c("SumOfErrors", "AIC")

# Loop through all of the reactions for forward stepwise selection
# glmFitFull is the full model, glmFitNull is the null model
for (i in 1:232) {
  glmFitFull<-glm(adrmat[, i] ~ ., data = mydatasub1, family = binomial)
  glmFitNull<-glm(adrmat[, i] ~ 1, data = mydatasub1, family = binomial)
  aDrugStep <- step(glmFitNull, scope = list(upper=glmFitFull, data = mydatasub1, direction = "forward"))
  
  # Storing the sum of errors & AIC
  pred1 <- predict(aDrugStep, type = "response")
  errorTable <- table(adrmat[, i], round(pred1))
  if (sum(round(pred1)) != 0) {
    numErrors[i, 1] <- sum(adrmat[, i]) + errorTable[1,2]
    numErrors[i, 2] <- aDrugStep$aic
  } else {
    numErrors[i, 1] <- sum(adrmat[, i])
    numErrors[i, 2] <- aDrugStep$aic
  }
}

# The name of the side effect with the smallest no. of errors
x <- min(numErrors[,1])
y <- which(numErrors[,1]==x, arr.ind = TRUE)
print(c("The ADR that was predicted with the minimum error FORWARD STEPWISE:", y))

# show the model for the side effect with smallest no. of errors
glmFitFull<-glm(adrmat[, 88] ~ ., data = mydatasub1, family = binomial) #full model (among the first 50 genes)

glmFitNull<-glm(adrmat[, 88] ~ 1, data = mydatasub1, family = binomial) #our null model

minErrorStep <- step(glmFitNull, scope = list(upper=glmFitFull, data = mydatasub1, direction = "forward"))
print("Model for ADR with min errors")
print(minErrorStep$model)
write.csv(minErrorStep$model, file = "FwdModel.csv")

# Shows the name of the side effect with smallest AIC
x <- min(numErrors[,2])
y <- which(numErrors[,2]==x, arr.ind = TRUE)
print(c("The ADR that was predicted with the minimum AIC FORWARD STEPWISE:", y))

# Show the model for the side effect with smallest AIC
glmFitFull2<-glm(adrmat[, 88] ~ ., data = mydatasub1, family = binomial) #full model (among the first 50 genes)

glmFitNull2<-glm(adrmat[, 88] ~ 1, data = mydatasub1, family = binomial) #our null model

minErrorStep2 <- step(glmFitNull2, scope = list(upper=glmFitFull2, data = mydatasub1, direction = "forward"))
print("Model for ADR with min AIC")
print(minErrorStep2$model)


# Backward elimination logistic regression
# start with all 50 variables loop through all of the reactions

for (i in 1:232) {
  glmFitFull<-glm(adrmat[, i] ~ ., data = mydatasub1, family = binomial)
  glmFitNull<-glm(adrmat[, i] ~ 1, data = mydatasub1, family = binomial)
  aDrugStep <- step(glmFitFull, data=mydatasub1, direction = "backward")
  
  # Storing the sum of errors & AIC
  pred1 <- predict(aDrugStep, type = "response")
  errorTable <- table(adrmat[, i], round(pred1))
  if (sum(round(pred1)) != 0) {
    numErrors[i, 1] <- sum(adrmat[, i]) + errorTable[1,2]
    numErrors[i, 2] <- aDrugStep$aic
  } else {
    numErrors[i, 1] <- sum(adrmat[, i])
    numErrors[i, 2] <- aDrugStep$aic
  }
}

# The name of the side effect with the smallest no. of errors from backward elimination
x <- min(numErrors[,1])
y <- which(numErrors[,1]==x, arr.ind = TRUE)
print(c("The ADR that was predicted with the minimum error BACKWARD ELIMINATION:", y))

# Shows the model for the side effect with smallest no. of errors
print("Model for ADR with min errors")
glmFitFull<-glm(adrmat[, 88] ~ ., data = mydatasub1, family = binomial)

minErrorStep <- step(glmFitFull, data=mydatasub1, direction = "backward")
print(minErrorStep$model)
write.csv(minErrorStep$model, file = "BwdModel.csv")

# The name of the side effect with smallest AIC
x <- min(numErrors[,2])
y <- which(numErrors[,2]==x, arr.ind = TRUE)
print(c("The ADR that was predicted with the minimum AIC BACKWARD ELIMINATION:", y))

# Show the model for the side effect with smallest AIC
print("Model for ADR with min AIC")
print(minErrorStep$model)