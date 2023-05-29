
## 1. LOAD ALL LIBRARIES RELEVANT TO THIS PORTFOLIO ASSIGNMENT. 
install.packages("zoo")
install.packages("data.table")
install.packages("readr")
install.packages("ggplot2")
install.packages("plotly")
install.packages("reshape2")
install.packages("rollRegres")
install.packages("dplyr")

library(zoo)
library(data.table)
library(readr)
library(ggplot2)
library(plotly)
library(reshape2)
library(rollRegres)
library(dplyr)

## 2. LOAD ALL DATA RELEVANT TO THIS PORTFOLIO ASSIGNMENT.

## Test assets

p100 <- read.csv("./testassets2022-23-exam.csv",head=TRUE)
View(p100)
col.var <- expand.grid(1:5,1:5) ##(first number = ME, second number = Book to Market)
col.var <- paste(col.var[,1],col.var[,2],sep=".")
colnames(p100)[3:ncol(p100)] <- col.var
p100$month <- as.yearmon(as.character(p100$month), format = "%Y%m")
month_col <- p100[-which(is.na(p100[,2])), 2]
p100 <- p100[-which(is.na(p100[,2])), ]
p100 <- zoo(p100[,3:ncol(p100)],month_col)
p100 <- window(p100, start = "1963-07", end = "2021-10")

## FF 5 Factors
ff.5f <- read.csv("./factors2022-23-exam.csv",head=TRUE)
ff.5f<-ff.5f[ ,-1]
colnames(ff.5f) <- c("month","rm.rf","smb","hml","rmw","cma","rf")
ff.5f$month <- as.yearmon(as.character(ff.5f$month), format = "%Y%m")
ff.5f <- zoo(ff.5f[, -1], order.by = ff.5f$month) 
ff.5f <- window(ff.5f, start = "1963-07", end = "2021-10")
View(ff.5f)

##------------ Make pretty table --------------
pretty.table <- function(x){
  ## Assume names(x) is 1.1, ... ,5.5
  z <- strsplit(names(x),split="[.]")
  z <- data.frame(do.call("rbind",z))
  z <- apply(z,2,function(x) as.numeric(as.character(x)))
  colnames(z) <- c("column","rows")
  z <- data.frame(z)
  z$values <- x
  z$values <- round(z$values,3)
  res <- acast(z,rows ~ column, value.var="values")
  return(res)
}


## 3. PART A


## 1(a) Average monthly returns
avg <- apply(p100,2,mean,na.rm=TRUE)
pretty.table(avg)


## 2-3 line summary of what you find here
## Here, we know that the rows are the size and columns are the value factors. As we go down the rows, we go from small to large size and on the other hand, as we go right to the column, we go low value from high value.
## Keeping that in mind, we can see that as we move from small to large stock,the return is actually increasing except the last one and on the contrary, as we move from the low value to high value stock, the return is increasing except the third year.
## Now, if we look at the smallest and highest value stock in this matrix, it is, indeed giving me the highest return; however, if we look at the largest and the lowest value stock, it not probably giving the lowest return.

## The matrix of the average return has been stored below.
a_answer <- pretty.table(avg)

## 1(b) Std. Dev and Sharpe Ratio calculations

## Standard Deviation:
sdev <- apply(p100,2,sd,na.rm=TRUE)
pretty.table(sdev)

## Sharpe Ratio
rf <-ff.5f$rf # This is the risk-free rate 
eret <- p100 - rf # Compute excess returns
sharpe <- apply(eret,2,function(x){mean(x,na.rm=TRUE)/sd(x,na.rm=TRUE)})
pretty.table(sharpe)

## 2-3 line summary of what you find here
## In the standard deviation matrix, we can see see interesting findings. We can see that the smallest and the lowest value stock has the highest standard deviation and as we move from small to large, the standard deviation is decreasing. Additionally, as we are moving from low value to high value stock the standard deviation is decreasing except the last one.
## In the sahrpe ratio matrix, we have the monthly sharpe ratio which is a risk-adjusted return. As we are moving from the small to the large stock, the sharpe ratio is increasing and as we are moving from the low value to high value stock the sharpe ratio is incresing significatly.
## Moreover, the smallest and the highest value stock has the highest sharpe ratio.

## The matrix of the Standard Deviation has been stored below.
b_answer1 <- pretty.table(sdev)

## The matrix of the Sharpe Ration has been stored below.
b_answer2 <- pretty.table(sharpe)

## 1(c.i) FMB Regression Stage 1
## Factor model
run.factor.model <- function(lhs,rhs){
alphas <- NULL
betas <- NULL
resids <- zoo(rep(NA,nrow(lhs)),index(lhs))
  
## Check they are of the same window
start.win <- max(start(lhs),start(rhs))
end.win <- min(end(lhs),end(rhs))
lhs <- window(lhs,start=start.win,end=end.win)
rhs <- window(rhs,start=start.win,end=end.win)
## Store residuals and alphas from regressions
for(i in 1:NCOL(lhs)){
  m <- lm(lhs[,i] ~ rhs)
  alphas <- c(alphas,coef(m)[1])
  b <- coef(m)[2:length(coef(m))]
  betas <- rbind(betas,b)
  resids <- merge(resids,residuals(m),all=TRUE)
  colnames(resids) <- NULL
}
resids <- resids[,2:NCOL(resids)]
colnames(resids) <- names(alphas) <- 1:NCOL(lhs)
colnames(betas) <- colnames(rhs)
rownames(betas) <- 1:NCOL(lhs)
  
## Various variables for the GRS statistic
T <- nrow(lhs)
N <- ncol(lhs)
K <- ncol(rhs)
Sigma <- var(resids,na.rm=TRUE)
demeaned <- apply(rhs,2,function(x){x - mean(x,na.rm=TRUE)})
Omega <- (t(demeaned)%*%demeaned)/T
f <- apply(rhs,2,mean,na.rm=TRUE)
GRS <- ((T-N-K)/N)*solve(1+t(f)%*%solve(Omega)%*%f)*t(alphas)%*%solve(Sigma)%*%alphas
GRS <- as.numeric(GRS)
finres <- list(alpha=alphas,beta=betas,residuals=resids,GRS=GRS,T=T,N=N,K=K,Sigma=Sigma,Omega=Omega)
print(paste("The GRS statistic is ", round(GRS,3) ,sep=""))
 return(finres)
}

rhs.m1 <- ff.5f[,1:5]
rhs.m2 <- ff.5f[,1:3]
# Model 1: FF 5 Factors -- This uses a function written for this purpose in the lines of code above. 
m1 <- run.factor.model(lhs=eret,rhs=rhs.m1)
m2 <- run.factor.model(lhs=eret,rhs=rhs.m2)
pv1 <-pf(m1$GRS,694,5,lower.tail = FALSE)
pv2 <-pf(m2$GRS,696,3,lower.tail = FALSE)
## 2-3 line summary of what you find here
## Here, I have tried to create two models and the first model with 5 factors that includes all and the second model is 3 factors including Market,SMB, and HML. The GRS statistics for the for the both models are 3.029 and 3.785 respectively. 
## Then we calculated the P-Values to check whether they are statistically significant or not. The reported P-Values for the both models are .1051 and .1488 respectively. So, with .05 significance level, we fail to reject the null for the both GRS statistics.
##Considering the associated GRSs, we can say that the first model that includes all factors has lower GRS than second model with 3 factors. However, in both cases, we fail to reject the null, so the GRSs are not statistically significant. 

c_i_answer <- list()
c_i_answer[[1]] <- c(GRS1=m1$GRS,P_Value1=pv1)
c_i_answer[[2]] <- c(GRS2=m2$GRS,P_Value2=pv2)

## 1(c.ii) FMB Regression Stage 2
## To observe what range of values are estimated for the various $\beta$s estimated in the regression, we plot an empirical probability density function (PDF): 

dat = data.frame(m1$beta)
p <- ggplot(dat, aes(rm.rf)) + geom_density()
p <- p + labs(x = "The Market Beta from the Fama-French 5 factor model", y = "Density",title="Fama-French model estimates")
p <- p +  theme(text = element_text(size=7),axis.text.x = element_text(angle=90, hjust=1)) 
ggplotly(p)

dat = data.frame(m1$beta)[,-1]
dat = melt(dat)
p <- ggplot(dat, aes(value)) + geom_density() + facet_grid(. ~ variable)
p <- p + labs(x = "The Betas from the Fama-French 5 factor model", y = "Density",title="Fama-French model estimates")
p <- p +  theme(text = element_text(size=7),axis.text.x = element_text(angle=90, hjust=1)) 
ggplotly(p)

## Are these "Priced" Factors? 

second.step <- function(lhs,rhs,mu.m,var.m){
  lambda <- NULL
  for(i in 1:NCOL(lhs)){
    m <- lm(lhs[,i] ~ rhs)
    tmp <- t(summary(m)$coefficients[,c(1,3)])
    lambda <- rbind(lambda,tmp[1,])
  }
  c.alpha <- lambda[,1]
  lambda <- lambda[,2:NCOL(lambda)]
  mean.alpha <- mean(c.alpha)
  mean.lambda <- apply(lambda,2,mean)
  
  var.alpha <- (1/(NCOL(lhs)*(NCOL(lhs)-1)))*sum((c.alpha - mean.alpha)^2)
  var.lambda <- (1/(NCOL(lhs)*(NCOL(lhs)-1)))*apply((lambda - mean.lambda)^2,2,sum)
  
  res <- rbind(c(mean.alpha,mean.lambda),c(mean.alpha,mean.lambda)/c(sqrt(var.alpha),sqrt(var.lambda)))
  rownames(res) <- c("lambda","t(lambda)")
  colnames(res) <- c("alpha",colnames(res)[2:NCOL(res)])
  return(res)
}

## FMB: Model 1
fmb1 <- second.step(lhs=t(eret),rhs= m1$beta)
## FMB: Model 2
fmb2 <- second.step(lhs=t(eret),rhs= m2$beta)


## 2-3 line summary of what you find here
## We can see that in both models, we have positive Alphas and both are statistically significant. One more thing we can see is that in both models, market factors are negative and they are also statistically insignificant.
## In the first model, rest of the coefficients are positive even though the last factor CMA is statistically insignificant. On the other hand, in the second model, rest of the coefficients are positive too, but the SMB factor is statistically insignificant.
## Even after considering the five factor model, we still have significant alpha left in the table. We can also see that CAPM model fails spectacularly.

c_ii_answer <- list(fmb1, fmb2)

## 2(a) Rolling factor betas (3-year).
#Setting up the function to compute 3 year rolling beta 
#Setting up the function to compute 3 year rolling beta 
run.factor.model.roll <- function(lhs,rhs){
  alphas <- NULL
  rf_betas <- NULL
  smb_betas <- NULL
  hml_betas <- NULL
  
  ## Uses rollapply function to the regression of the rolling window
  for(i in 1:ncol(lhs)){
    
    ff_roll <- merge(lhs[,i], rhs, all = FALSE)
    colnames(ff_roll) <- c('returns', "rm.rf", "smb", "hml")
    
    #setting the rolling window size of 36 
    fits <- roll_regres(returns ~ rm.rf + smb + hml, ff_roll, width = 36)
    return (fits$coefs[,-1])
  }
}

roll.model <- run.factor.model.roll(lhs = eret,rhs = rhs.m2)
roll.model <- na.omit(roll.model)

## 2-3 line summary of what you find here
## We can see that the value factor that we have considered is giving us a negative coefficients which means as we increase the value factor by one unit, our excess return is going down.
## On the other hand, we have consistently positive coefficients for the size factor which means as we load up an additional unit in size factor, our excess return goes up.
## Market factor is also generating positive coefficients.

Second_a_answer <-roll.model

#Setting up the function 
model.fit <- function(lhs,rhs){
  predval <- list()
  for(i in 1:NCOL(lhs)){
    m <- lm(lhs[,i] ~ rhs, na.action = na.exclude)
    predval[[i]] <- predict(m)
  }
  predictions <- do.call("rbind",predval)
  return(data.frame(predictions))
}


#Creating an equal length for both dataset  
neret <- tail(eret, nrow(eret)-35)

model.new <- model.fit(lhs=neret,rhs=roll.model)

excdata <- as.data.frame(neret)

#Transposing the dataset of predicted values
model.new.new <- t(model.new)

cor(model.new.new[,1],excdata[,1])

Output_list <- list()

#loop to generate the correlation
for (i in 1:25){
  
  output <- cor(model.new.new[,i], excdata[,i])
  
  Output_list <- append(Output_list, output )
  
}

cor25 <- Output_list

#Estimating cross section regression excess return and rolling betas 
model.fit <- function(lhs,rhs){
  predval <- list()
  for(i in 1:NCOL(lhs)){
    m <- lm(lhs[,i] ~ rhs, na.action = na.exclude)
    predval[[i]] <- predict(m)
  }
  predictions <- do.call("rbind",predval)
  predictions <- apply(predictions,2,mean,na.rm=TRUE)
  actuals <- apply(lhs,1,mean,na.rm=TRUE)
  return(data.frame(actuals,predictions))
}


#plotting the data to check predictibility

model.new <- model.fit(lhs=neret,rhs=roll.model)
model.new$type = "Full Sample"
p <- ggplot(data = model.new, aes(x = predictions, y = actuals)) + geom_point(color='red') +  geom_smooth(method = "lm", se = TRUE) + facet_grid(type ~ .)
p <- p + labs(x = "Predicted Expected Returns", y = "Sample Expected Returns",title="Model Validity: Model 1")
ggplotly(p)

## 2-3 line summary of what you find here
## From the reported 25 coefficients, we can see that there is a poor correlation between the predicted values and the actual observed values. All the correlation results are below 10% which clearly shows the poor relationship between the both variables. 
## To illustrate this, we have also considered the additional step of plotting these two variables to look at the fitted line. From the graph that we have generated above, we can see a straight line which is not upward or downward sloping.
## The straight line that we have generated means slope is zero as it is not upward or downward sloping. If the slope is zero,then it is means there is no predictability in the model and that leads up to conclude that the model is not good at all.

Second_b_answer <-cor25  ## This is how I expect you to store the 25 correlation numbers asked for in the portfolio. 



## PART B


## 3 (a)
pe <- read.csv("./prediction-exercise-2022-23.csv",head=TRUE)
pe<-data.frame(pe)
pe<-pe[,-1]
pe$yyyymm <- as.yearmon(as.character(pe$yyyymm), format = "%Y%m")
#View(pe)
## Checking the percentage of NA value in each column
na_percentages <- apply(pe, 2, function(x) mean(is.na(x)))*100

## b. Create 1 month lag values of all of the predictor variables

## dropping the columns that has excessive NA values. We are dropping the "csp" column as it has nearly 35% NA values.
variables <- c("D12", "E12", "b.m", "tbl", "AAA", "BAA", "lty", "ntis", "Rfree", "infl", "ltr", "corpr", "svar")
for(i in variables) pe[, paste0(i,"_lag")] <- shift(pe[,i],n=1,type="lag")

## The benchmark forecasting model
mean.returns <- rep(NA,1176) 
for(t in 1177:nrow(pe)){
  ## Compute mean until data point t-1: 
  mr <- mean(pe$return[1:(t-1)],na.rm=TRUE)
  mean.returns <- c(mean.returns,mr)
}
pe$mean.returns <- mean.returns

## Compute the Out of Sample Squared Prediction error for the "historical mean return" model:

## It is the difference between "return" in month $t$  and the "mean.returns" until then.
pe$oos_meanret <- (pe$return - pe$mean.returns)^2

## Conver it to a timeseries object for easy work
pe <- zoo(pe[,2:ncol(pe)],order.by=pe$yyyymm)

## The forecasting model 

## Input:
## y = the outcome variable (Returns at t+1)
## x = the predictor variable (Variable at time t) 

## Output:

## The OOS Squared Prediction Error for the model

## Function: 
forecast.model <- function(y,x,initial.sample.end="Jan 2019"){
  fulldat <- merge(y,x)
  months <- index(fulldat)[121:nrow(fulldat)]
  oos.forecasterror <- as.numeric(rep(NA,(nrow(fulldat)-length(months))))
  for(t in 1:length(months)){
    ## create the dataset for regression:
    tmp.d <- window(fulldat,end=(as.yearmon(months[t])-1/12))
    ## Run a linear regression:
    m <- lm(y ~ x,data=tmp.d) 
    
    ## Predict return for time t:
    pt <- data.frame("x"= x[months[t]])
    haty <- predict(m,newdata=pt)
    
    ## Compute out of sample forecast error: 
    oose <- as.numeric(y[months[t]]-haty)
    oos.forecasterror <- c(oos.forecasterror,oose)
  }
  return(oos.forecasterror^2) 
} 

predictors <- paste0(variables,"_lag")

results <- NULL 
for(i in predictors){
  cat("Predicting using:",i,"\n") 
  tmp <- forecast.model(y=pe$return,x=pe[,i]) 
  results <- cbind(results,tmp)
} 

colnames(results) <- paste0(predictors,"_oosse")
pe <- merge(pe,results)

## Evaluate Performance

predictors <- c("D12_lag", "E12_lag", "b.m_lag", "tbl_lag", "AAA_lag", "BAA_lag", "lty_lag", "ntis_lag", "Rfree_lag", "infl_lag", "ltr_lag", "corpr_lag", "svar_lag")

perf <- NULL 
for(i in predictors){
  tmp <- c(rep(NA,1176),cumsum(na.omit(pe$oos_meanret - na.omit(pe[,paste0(i,"_oosse")]))))
  perf <- cbind(perf,tmp)
}
colnames(perf) <- predictors
perf <- zoo(perf,order.by=index(pe))
#View(perf)
##Now we plot each of the variables, and compare the graphs to the annual estimates

toplot <- na.locf(perf)
# Finding which model gives me the lowest error against the benchmark model. 
col_sums <- colSums(toplot)
print(col_sums_sorted <- sort(col_sums, decreasing = TRUE))

#transforming the toplot to data frame

tt<-as.data.frame(toplot)
plot.vars <- c("D12_lag_perf", "E12_lag_perf", "b.m_lag_perf", "tbl_lag_perf", "AAA_lag_perf", "BAA_lag_perf", "lty_lag_perf", "ntis_lag_perf", "Rfree_lag_perf", "infl_lag_perf", "ltr_lag_perf", "corpr_lag_perf", "svar_lag_perf")
for(i in 1:length(plot.vars)){
  plot(tt[,i],
       xlab="month",ylab="Cumulative SSE Difference",
       main=paste0("Variable: ",plot.vars[i]),type = "l")
  grid()
}

## I am storing my best model here.
pdf("model1.pdf")

plot(tt$ntis_lag, xlab="month",ylab="Cummulative SSE Difference", main="ntis_lag_perf", type ="l")

dev.off()






