#Functions for statistical calculations of Simple linear regression models
#Created by Tomi Lehto
#
#This package contains some simple and useful functions for basic statistical calculations,
#mostly related to the Simple linear regression model. These functions could be useful for
#example when studying simple linear regression models. The function printAll is the most usefull
#function in this collection, because it gives fast overview about all calculation
#related to simple linear regression.

#-------------------------------------------------------------------------------------------

#Sums multiplications of pairs of two vectors indices 0...n
SumOfPairs <- function(x2,y2){
  s=0
  i = 1
  for (i in 1:length(x2)) {
    s <- s + x2[i]*y2[i]  
  }
  return(s)
}
#-------------------------------------------------------------------------------------------

#Returns the Sum of multiplied errors from the mean of vector x and y (Sxy)
Sxy <- function(x,y){
  return(SumOfPairs(x,y) - length(x)*mean(x)*mean(y))
}
#-------------------------------------------------------------------------------------------

#Returns the Sum of squared errors from the mean for given vector a. (Saa)
Saa <- function(a){
  return(sum(a^2) - length(a)* (mean(a)^2))
}
#-------------------------------------------------------------------------------------------

#Returns estimated beta for Simple linear regression model calculated from vectors x and y
hatB <- function(x,y){
  return(Sxy(x,y)/Saa(x))
}
#-------------------------------------------------------------------------------------------

#Returns estimated alpha for Simple linear regression model calculated from vectors x and y.
hatA <- function(x,y){
  return(mean(y)-mean(x)*hatB(x,y))  
}
#-------------------------------------------------------------------------------------------

#returns Residual Sum of Squares (RSS) calculated from vectors x and y.
SSr <- function(x,y){
  return((Saa(x)*Saa(y)-Sxy(x,y)^2)/Saa(x))
}
#-------------------------------------------------------------------------------------------

#Returns estimated varians (sigma^2) of Simple regression Model
hatS <- function(x,y){
  return(SSr(x,y)/(length(x)-2))
}
#-------------------------------------------------------------------------------------------

#Returs Coefficient of determination, (R^2- value), calculated from vectors x and y
R2 <- function(x,y){
  return(1-(SSr(x,y)/Saa(y)))
}
#-------------------------------------------------------------------------------------------

#Returns the parameter W for the Prediction interval of Simple linear regression model.
#W is calculated from vectors x, y and from estimation point x0.
#valiW <- function(x,y,x0){
 # return(sqrt((1+(1/length(x)) + (((x0-mean(x))^2)/(Saa(x))))*hatS(x,y)))
#}
valW <- function(x,y,x0){
  return(sqrt((1+(1/length(x)) + (((x0-mean(x))^2)/(Saa(x))))*(SSr(x,y)/(length(x)-2))))
}
#-------------------------------------------------------------------------------------------

#Returns estimated beta values: B0,B1,B2 for two-variable Simple linear regression model.
#beta values are calculated from vectors of Explanatory variables x1, x2 and
#from vector of Dependent variable y.

hatBm <- function(x1,x2,y){
  mx <- matrix(c(rep(1,length(x1)),x1,x2),byrow = FALSE,nrow = length(x1))
  return(solve(t(mx)%*%mx)%*%t(mx)%*%y)
}

#-------------------------------------------------------------------------------------------
#Returns the prediction interval calculated on given parameters x, y, x0 and aplha, where x is
#vector of the explanatory variables and y is vector of the response variables. x0 is the
#estimated explanatory variable. Aplha value corresponds to the value of  1-d where
#d is the desired level of confidence.

preInt <- function(x,y,x0,alpha){
  w <- valW(x,y,x0)
  a <- hatA(x,y)
  b <- hatB(x,y)
  f <- length(x)-2
  t <- qt(1-(alpha/2),f) 
  
  return(c(a+b*x0-t*w,a+b*x0+t*w))
}

#-------------------------------------------------------------------------------------------
#
#printAll prints all functions which relate to the simple linear regression model with given
#explanatory variables vector x and response variable vector y. Only function that is not printed
#is the hatBm because it calculates the hatB values for two variable linear regression. 
#x0 and Alpha are optional and if they are provided the funtion prints also 
#the prediction interval and the W-value for the calculations of the prediction interval. 

printAll <- function(x,y,x0,alpha){
  printPaste("Sum of pairs:",SumOfPairs(x,y))
  printPaste("Sxy:",Sxy(x,y))
  printPaste("Sxx:",Saa(x))
  printPaste("Syy:",Saa(y))
  printPaste("hatB:",hatB(x,y))
  
  printPaste("hatA:",hatA(x,y))
  printPaste("SSr:",SSr(x,y))
  printPaste("hatS:",hatS(x,y))
  printPaste("R^2:",R2(x,y))
  if(!missing(x0)) {
    if(!missing(alpha)){
      print(paste("Prediction_interval:", "[",preInt(x,y,x0,alpha)[1],",",preInt(x,y,x0,alpha)[2],"]"))
    }
    printPaste("W-value:", valW(x,y,x0))
  }
  
} 
  
#---------------------------------------------------------------------------------------------------------
#printPaste is a helper funtion for te printAll function. It transforms x to character and 
#then pastes and prints x with s.
printPaste <- function(s,x){
  print(paste(s,as.character(x)))
}
