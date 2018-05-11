#####################################################
##  practical time series analyss using r language ##
#####################################################

#install.packages()
#library()
#old.packages()
#update.packages()

###################################################################
##weekly topic list      
#week1: basic statistic review + measuring linear association
#week2: autogressive processes definition and examples
#week3: partial autocorrelation and the pacf
#week4: measuring the quality of a time series model
#week5: X[t]= noise + AutoRegressive part + Moving Average part 
#week6: forecasting
###################################################################

###################################################################
#week1 example1: atmospheric co2 concentration
###################################################################
help(co2)
class(co2)
summary(co2)
plot(co2,main='atmospheric co2 concentration')
co2.values=as.numeric(co2)
co2.times=as.numeric(time(co2))
ssxx=sum((co2.times-mean(co2.times))*(co2.times-mean(co2.times)))
ssxy=sum((co2.values-mean(co2.values))*(co2.times-mean(co2.times)))
slope=ssxy/ssxx
slope #1.307497
intercept=mean(co2.values)-slope*mean(co2.times)
intercept #-2249.774
co2.linear.model=lm(co2~time(co2))
plot(co2,main='atmosphere co2 concentration with fitted line')
abline(co2.linear.model)
#co2.fitted.values=slope*co2.times+intercept
#co2.residuals=co2.values-co2.fitted.values
co2.residuals=resid(co2.linear.model)
par(mfrow=c(1,3))
hist(co2.residuals,main='histogram of co2 residuals')
qqnorm(co2.residuals,main='normal probability plot')
qqline(co2.residuals)

plot(co2.residuals~time(co2),main='residuals on time')#xlim=c(1960,1990+)
plot(co2.residuals~time(co2),xlim=c(1960,1963),main='residuals on time')

###################################################################
#week1 example2: extra sleep in Gossett data by group
###################################################################
class(sleep)
summary(sleep)
plot(extra~group,data=sleep,main='extra sleep in Gossett data by group')
attach(sleep)#do not use dollar sign
extra.1=extra[group==1]
extra.2=extra[group==2]
t.test(extra.1,extra.2,paired=TRUE,alternative="two.sided")
#t=-4.0621, df=9, p-value=0.002833
diffs=extra.1-extra.2
qqnorm(diffs,main='norm probability plot')
qqline(diffs)

plot(extra.2~extra.1,xlab='extra sleep with drug1', ylab='extra sleep with drug 2',
     main='extra sleep drug2 against drug1')
sleep.linear.model=lm(extra.2~extra.1)
abline(sleep.linear.model)

summary(sleep.linear.model)
#################################################################
#Residuals:
#  Min      1Q    Median    3Q     Max 
#-1.6735 -0.4673 -0.3365  0.3979  2.9375 

#Coefficients:
#  Estimate    Std.      Error   t-value Pr(>|t|)   
#(Intercept)  1.6625     0.4452   3.734  0.00575 **
#  extra.1    0.8899     0.2399   3.709  0.00596 **
#  ---
#  Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

#Residual standard error: 1.288 on 8 degrees of freedom
#Multiple R-squared:  0.6323,	Adjusted R-squared:  0.5863 
#F-statistic: 13.76 on 1 and 8 DF,  p-value: 0.005965
#################################################################

#our model: extra.2=0.8899*extra.1+1.6625
extra.2[3]-(0.8899*extra.1[3]+1.6625) #-0.38452
sleep.residuals=resid(sleep.linear.model)
sleep.residuals[3] #-0.3845478

#################################################################
#week1 example3: girth, height and volume for black cherry trees
#################################################################
summary(trees)
pairs(trees, pch=21, bg=c('red'))
cov(trees)
cor(trees)


#################################################################
#week2: simulating a simple AR(p) process >> AR(1) & AR(2) & AR(3)
#################################################################
set.seed(2018);
N=1000;Z=rnorm(N,0,1);
X=NULL;

X[1]=Z[1];

phi1=0.4;phi2=1.0;
for (t in 2:N){
  X1[t]=Z[t]+phi1*X1[t-1];
  X2[t]=Z[t]+phi2*X2[t-1];
}

par(mfrow=c(4,2))
plot(ts(X1),main='AR(1) time series phi=0.4')
X1.acf=acf(ts(X1),main='Autocorrelation of X[t]=Z[t]+0.4*X[t-1]')

plot(ts(X2),main='AR(1) time series phi=1.0')
X2.acf=acf(ts(X2),main='Autocorrelation of X[t]=Z[t]+X[t-1]')

X.ts<-arima.sim(list(ar=c(.7,.2)),n=1000)
plot(X.ts,main='AR(2) time series phi1=0.7 phi2=0.2')
X.acf=acf(X.ts,main='Autocorrelation of X[t]=Z[t]+0.7*X[t-1]+0.2*X[t-2]')

X.ts<-arima.sim(list(ar=c(.5,-.4)),n=1000)
plot(X.ts,main='AR(2) time series phi1=0.5 phi2=-0.4')
X.acf=acf(X.ts,main='Autocorrelation of X[t]=Z[t]+0.5*X[t-1]-0.4*X[t-2]')

#################################################################
#week2: backshift operator and the acf
#################################################################
set.seed(2018);
N=1000;Z=rnorm(N,0,1);
X1=NULL;X2=NULL;X3=NULL;X4=NULL;
X1[1]=Z[1];X2[1]=Z[1];X3[1]=Z[1];X4[1]=Z[1];

phi1=0.9;phi2=0.2;phi3=-0.3;phi4=-.8;
for (t in 2:N){
  X1[t]=Z[t]+phi1*X1[t-1];
  X2[t]=Z[t]+phi2*X2[t-1];
  X3[t]=Z[t]+phi3*X3[t-1];
  X4[t]=Z[t]+phi4*X4[t-1];
}

par(mfrow=c(2,2))
X1.acf=acf(ts(X1),main='phi=0.9')
X2.acf=acf(ts(X2),main='phi=0.2')
X3.acf=acf(ts(X3),main='phi=-0.3')
X4.acf=acf(ts(X4),main='phi=-0.8')

#################################################################
#week3: partial autocorrelation and the pacf
#################################################################

## a moving average process MA(q) of order q has an ACF that cuts off after q lags ##
## a autogressive process AR(p) of order p has an PACF that cuts off after p lags ##

rm(list = ls(all=TRUE))
par(mfrow=c(3,1))
phi.1=.6; phi.2=.2; 
data.ts=arima.sim(n=500,list(ar=c(phi.1,phi.2)));
plot(data.ts,main = "AR(2) time series phi1=0.6 phi2=0.2")
acf(data.ts,main="Autocorrelation of X[t]=Z[t]+0.6*X[t-1]+0.2*X[t-2]")
acf(data.ts,type="partial",main="partial autocorrelation function")

rm(list = ls(all=TRUE))
par(mfrow=c(3,1))
phi.1=.9; phi.2=-.6; phi.3=.3;
data.ts=arima.sim(n=500,list(ar=c(phi.1,phi.2,phi.3)));
plot(data.ts,main = "AR(3) time series phi1=0.9 phi2=-0.6 phi3=0.3")
acf(data.ts,main="Autocorrelation of X[t]=Z[t]+0.9*X[t-1]-0.6*X[t-2]+0.3*X[t-3]'")
acf(data.ts,type="partial",main="partial autocorrelation function")

#################################################################
#week3 example1: the beveridge wheat price dataset
#################################################################
library(readxl)
beveridge <- read_excel("C:/Users/qren/Downloads/beveridge.xls")
#beveridge <- read_excel("beveridge.xls")
#View(beveridge)
beveridge.ts=ts(beveridge[,2],start = 1500)
plot(beveridge.ts,main='beveridge wheat price data')
beveridge.MA=filter(beveridge.ts,rep(1/31,31),sides = 2)
lines(beveridge.MA,col='red')

par(mfrow=c(3,1))
Y=beveridge.ts/beveridge.MA
plot(Y,main = 'transformed beveridge data')
acf(na.omit(Y),main="acf of transformed data")
acf(na.omit(Y),type="partial",main="Pacf of transformed data")

ar(na.omit(Y),order.max=5)

################################################
#Coefficients:
#  1        2  
#0.7239  -0.2957  

#Order selected 2  sigma^2 estimated as  0.02692
################################################

###################################################################
#week3 example2: the bodyfat dataset
###################################################################
library(faraway)
#library(isdals)
data(bodyfat)
attach(bodyfat)
pairs(cbind(Fat,Triceps,Thigh,Midarm))
cor(cbind(Fat,Triceps,Thigh,Midarm))

################################################
#            Fat    Triceps     Thigh    Midarm
#Fat     1.0000000 0.8432654 0.8780896 0.1424440
#Triceps 0.8432654 1.0000000 0.9238425 0.4577772
#Thigh   0.8780896 0.9238425 1.0000000 0.0846675
#Midarm  0.1424440 0.4577772 0.0846675 1.0000000
################################################

Fat.pred=predict(lm(Fat~Thigh))
Triceps.pred=predict(lm(Triceps~Thigh))
cor((Fat-Fat.pred),(Triceps-Triceps.pred)) #0.1749822

Fat.pred=predict(lm(Fat~Thigh+Midarm))
Triceps.pred=predict(lm(Triceps~Thigh+Midarm))
cor((Fat-Fat.pred),(Triceps-Triceps.pred)) #0.33815

library(ppcor)
pcor(cbind(Fat,Triceps,Thigh,Midarm))


###################################################################
#week4: measuring the quality of a time series model (SSE & AIC)
#We prefer a model with a lower AIC
###################################################################
rm(list=ls(all=TRUE))
set.seed(43)
data=arima.sim(list(order=c(2,0,0),ar=c(0.7,-0.2)), n=2000)
par(mfrow=c(1,2))
acf(data, main='ACF of AR(2)')
acf(data, type='partial', main='PACF of AR(2)')
arima(data, order=c(2,0,0),include.mean=FALSE)

########################################################################
#Coefficients:
#        ar1      ar2
#      0.7111  -0.1912
#s.e.  0.0219   0.0220

#sigma^2 estimated as 0.9985:  log likelihood = -2836.64,  aic = 5679.27
########################################################################

SSE=NULL
AIC=NULL
for(p in 1:5){
  m=arima(data,order=c(p,0,0),include.mean = FALSE)
  SSE[p]=sum(resid(m)^2)
  AIC[p]=m$aic
  print(m$coef)
  print(paste("AIC: ",m$aic,"SSE: ",sum(resid(m)^2)))
}
par(mfrow=c(1,2))
order=c(1,2,3,4,5)
plot(SSE~order,main='sse plot',ylim=c(1800,2100))
plot(AIC~order,main='aic plot',ylim=c(5500,5800))

########################################################################
#ar1 
#0.5969948 
#"AIC:  5751.73196762524 SSE:  2072.83193501059"
#ar1        ar2 
#0.7111457 -0.1911552 
#"AIC:  5679.27375222458 SSE:  1997.00667996082"         >> BEST!!!
#ar1         ar2         ar3 
#0.71359315 -0.20027406  0.01281966 
#"AIC:  5680.94495534325 SSE:  1996.67791506654"
#ar1          ar2          ar3          ar4 
#0.713676747 -0.201599645  0.017553047 -0.006629412 
#"AIC:  5682.85704377107 SSE:  1996.58997811327"
#ar1         ar2         ar3         ar4         ar5 
#0.71410825 -0.20268672  0.03019322 -0.05154692  0.06293048 
#"AIC:  5676.91730818182 SSE:  1988.65973372245"
########################################################################

#OUR MODEL: X[t]=Z[t]+0.7111457*X[t-1]-0.1911552*X[t-2]


###################################################################
#week5: ARMA(p,q) = AR(p) + MA(q) properties & examples
###################################################################
rm(list=ls(all=TRUE))
set.seed(500)
data=arima.sim(list(order=c(1,0,1),ar=0.7,ma=.2),n=1000000)
par(mfcol=c(3,1))
plot(data, main='ARMA(1,1) time series: X[t]= 0.7*X[t-1]+Z[t]+0.2*Z[t-1]', xlim=c(0, 400))
acf(data, main='Autocorrelation of ARMA(1,1)')
acf(data, type='partial', main='Partial Autocorrelation of ARMA(1,1)')

# example dataset: scientific discoveries in a year
par(mfrow=c(4,1))
plot(discoveries, main='Time series of no. of major scientific discoveries in a year')
stripchart(discoveries, method="stack", offset=0.5, at=0.15, pch=19, 
           main='no. of discoveries dotplot', 
           xlab='no. of major scienfic discoveries in a year', ylab = "frequency")

acf(discoveries,main='ACF')
acf(discoveries,type = 'partial',main='PACF')

#HOW TO CHOOSE order=(p, d, q)?
library(forecast)
auto.arima(discoveries, d=0, ic='aic', approximation=TRUE)
###########################################################
#ARIMA(1,0,1) with non-zero mean 

#Coefficients:
#        ar1      ma1    mean
#      0.8353  -0.6243  3.0208
#s.e.  0.1379   0.1948  0.4728

#sigma^2 estimated as 4.538:  log likelihood=-216.1
#AIC=440.2   AICc=440.62   BIC=450.62
###########################################################

#CHOOSE order=(1, 0, 1)


###################################################################
#week6 example1: simple exponential smoothing
#Holt-winters exponential smoothing without trend and without seasonal component
###################################################################

rm(list=ls(all=TRUE))
rain.data <- scan("http://robjhyndman.com/tsdldata/hurst/precip1.dat",skip=1)
rain.ts <- ts(rain.data, start=c(1813))
par( mfrow=c(2,2) )
hist(rain.data, main="Annual London Rainfall 1813-1912", xlab="rainfall in inches")
qqnorm(rain.data,main="Normal Plot of London Rainfall")
qqline(rain.data)
plot.ts(rain.ts, main="Annual London Rainfall 1813-1912", xlab="year", ylab="rainfall in inches")
acf(rain.ts, main="ACF: London Rainfall")

library(forecast)
auto.arima(rain.ts)
############################################################
#ARIMA(0,0,0) with non-zero mean 

#Coefficients:
#        mean
#       24.8239
#s.e.   0.4193

#sigma^2 estimated as 17.76:  log likelihood=-285.25
#AIC=574.49   AICc=574.61   BIC=579.7
############################################################

alpha=0.2 #increase alpha for more rapid decay
n=length(rain.data) #100
forecast.values=NULL
forecast.values[1]=rain.data[1]
for(i in 1:n){
  forecast.values[i+1]=alpha*rain.data[i]+(1-alpha)*forecast.values[i]
}
paste("forecast for time ",n+1, "=", forecast.values[n+1])

#HOW TO CHOOSE alpha?
SSE=NULL
alpha.values=seq(0.001,0.999,by=0.001)
number.alphas=length(alpha.values)
for(k in 1:number.alphas){
  forecast.values=NULL
  alpha=alpha.values[k]
  forecast.values[1]=rain.data[1]
  for(i in 1:n){
    forecast.values[i+1]=alpha*rain.data[i]+(1-alpha)*forecast.values[i]
  }
  SSE[k]=sum((rain.data-forecast.values[1:n])^2)
}
par(mfrow=c(1,2))
plot(SSE~alpha.values, main="optimal alpha value min SSE")
plot(SSE~alpha.values, xlim=c(0, 0.03), main="optimal alpha value min SSE")

index.of.smallest.SSE=which.min(SSE) #24
alpha.values[which.min(SSE)] #0.024
paste("forecast for time ",n+1, "=", forecast.values[which.min(SSE)])
#CHOOSE alpha=0.024

HoltWinters(rain.ts, beta=FALSE, gamma=FALSE)

###################################################################
#Smoothing parameters:
#alpha: 0.02412151
#beta : FALSE
#gamma: FALSE

#Coefficients:
#    [,1]
#a 24.67819
###################################################################
#CHOOSE alpha=0.024


###################################################################
#week6 example2: double exponential smoothing
#Holt-winters exponential smoothing with trend 
###################################################################

#########      Holt-Winters for trend    ##################
#Holt-Winters for trend
#levels: with parameter alpha
#trend: with parameter beta
#seasonal component: with parameter gamma
###########################################################

library(readxl)
money.data= read_excel("C:/Users/qren/Downloads/volume-of-money-abs-definition-m.xls")
#money.data= read_excel("volume-of-money-abs-definition-m.xls")
money.data.ts=ts(money.data[,2],start=c(1960,2),frequency=12)
par(mfrow=c(3,1))
plot(money.data.ts,main='Time plot of volume of money')
acf(money.data.ts,main='ACF of volume of money')
acf(money.data.ts,type="partial",main='PACF of volume of money')

data=money.data[,2]
N=length(data)
alpha=0.7
beta=0.5

forecast=NULL
level=NULL
trend=NULL
level[1]=data[1]
trend[1]=data[2]-data[1]
forecast[1]=data[1]
forecast[2]=data[2]
for(n in 2:N){
  level[n]=alpha*data[n]+(1-alpha)*(level[n-1]+trend[n-1])
  trend[n]=beta*(level[n]-level[n-1])+(1-beta)*trend[n-1]
  forecast[n+1]=level[n]+trend[n]
}

forecast[3:N]
m=HoltWinters(data,alpha=0.7,beta=0.5,gamma=FALSE)
m$fitted[,1]
plot(m,main='Holt winters fitting of money volume with bogus parameters')
m=HoltWinters(data,gamma=FALSE)
plot(m,main='Holt winters fitting of money volume with optimal parameters')


###################################################################
#week6 example3: triple exponential smoothing
#Holt-winters exponential smoothing with trend and seasonality 
###################################################################

par(mfrow=c(4,1))
plot(AirPassengers, main='Number of monthly air passengers(in thousands)')
plot(log10(AirPassengers), main='log10: Number of monthly air passengers(in thousands)')
acf(log10(AirPassengers),main='ACF of airpassenger dataset')
acf(log10(AirPassengers),type="partial",main='PACF of airpassenger dataset')
HoltWinters(x=log10(AirPassengers), beta=FALSE, gamma=FALSE)

############################################################
#Smoothing parameters:
#alpha: 0.9999339
#beta : FALSE
#gamma: FALSE

#Coefficients:
#    [,1]
#a 2.635481(forecasted value)
#############################################################

m1=HoltWinters(x=log10(AirPassengers), beta=FALSE, gamma=FALSE)
m1$SSE #0.3065102
m2=HoltWinters(x=log10(AirPassengers))
m2$SSE #0.0383026

#############################################################
#Call:
#  HoltWinters(x = log10(AirPassengers))

#Smoothing parameters:
#alpha: 0.326612
#beta : 0.005744246
#gamma: 0.8207255

#Coefficients:
#         [,1]
#a    2.680598830
#b    0.003900787
#s1  -0.031790733
#s2  -0.061224237
#s3  -0.015941495
#s4   0.006307818
#s5   0.014138008
#s6   0.067260071
#s7   0.127820295
#s8   0.119893006
#s9   0.038321663
#s10 -0.014181699
#s11 -0.085995400
#s12 -0.044672707

#X[144+1]=2.680598830+1*0.003900787+(-0.031790733)
#X[144+8]=2.680598830+8*0.003900787+0.119893006
#X[144+15]=2.680598830+15*0.003900787+(-0.015941495)
#############################################################

#simplily using 'forecast' method
rm(list=ls(all=TRUE))
library('forecast')
m<-HoltWinters(log10(AirPassengers))
forecast(m)
plot(forecast(m))
#############################################################

#############################################################
#ADDED PART: Yule-Walker Equations
#############################################################

#AR(2) process: X[t]= 1/3*X[t-1]+1/2*Z[t]+Z[t-1]
#polynomial: phi(B)=1-1/3*B-1/2*B^2

sigma=4
phi=NULL
phi[1:2]=c(1/3,1/2)
n=10000
set.seed(2018)
ar.process=arima.sim(n,model=list(ar=c(1/3,1/2)),sd=4)
ar.process[1:5]
r=NULL
r[1:2]=acf(ar.process,plot=F)$acf[2:3]
R=matrix(1,2,2)
R[1,2]=r[1]
R[2,1]=r[1]
b=matrix(r,2,1)
phi.hat=matrix(c(solve(R,b)[1,1],solve(R,b)[2,1]),2,1)
phi.hat 
########################
#phi1 0.3348453 (1/3)
#phi2 0.5115046 (1/2)
########################
c0=acf(ar.process,type='covariance',plot=F)$acf[1]
c0 #41.43329
var.hat=c0*(1-sum(phi.hat*r))
var.hat #16.2185
par(mfrow=c(3,1))
plot(ar.process,main='AR(2) process: X[t]= 1/3*X[t-1]+1/2*Z[t]+Z[t-1]')
acf(ar.process,main='ACF')
pacf(ar.process,main='PACF')


#AR(3) process
sigma=4
phi=NULL
phi[1:3]=c(1/3,1/2,7/100)
n=100000
set.seed(2017)
ar.process=arima.sim(n,model=list(ar=c(1/3,1/2,7/100)),sd=4)
ar.process[1:5]
r=NULL
r[1:3]=acf(ar.process,plot=F)$acf[2:4]
r[1:3]=c(0.8,0.6,0.2)
R=matrix(1,3,3)
R[1,2]=r[1]
R[2,1]=r[1]
R[2,3]=r[1]
R[3,2]=r[1]
R[1,3]=r[2]
R[3,1]=r[2]
b=matrix(,3,1)
b[1,1]=r[1]
b[2,1]=r[2]
b[3,1]=r[3]
phi.hat=(solve(R,b))
phi.hat
c0=acf(ar.process,type='covariance',plot=F)$acf[1]
var.hat=c0*(1-sum(phi.hat*r))
var.hat
par(mfrow=c(3,1))
plot(ar.process,main='AR(3) process')
acf(ar.process,main='ACF')
pacf(ar.process,main='PACF')

# for real case
plot(JohnsonJohnson,main='jj')
jj.log.return=diff(log(JohnsonJohnson))
jj.log.return.mean.zero=jj.log.return-mean(jj.log.return)
par(mfrow=c(3,1))
plot(jj.log.return.mean.zero,main='Log-return (mean zero) of JohnsonJohnson')
acf(jj.log.return.mean.zero,main='ACF')
pacf(jj.log.return.mean.zero,main='PACF')

my.data = jj.log.return
ar.process=jj.log.return.mean.zero
p=4
r=NULL
r[1:p]=acf(ar.process,plot=F)$acf[2:(p+1)]
R=matrix(1,p,p)
for(i in 1:p){
  for(j in 1:p){
    if(i!=j)
      R[i,j]=r[abs(i-j)]
  }
}
b=NULL
b=matrix(r,p,1)
phi.hat=NULL
phi.hat=solve(R,b)[,1]
c0=acf(ar.process,type='covariance',plot=F)$acf[1]
var.hat=c0*(1-sum(phi.hat*r))
phi0.hat=mean(my.data)*(1-sum(phi.hat))
cat("Constant:", phi0.hat," Coeffcinets:", phi.hat, " and Variance:", var.hat, '\n')
#Constant: 0.079781  Coeffcinets: -0.6293492 -0.5171526 -0.4883374 0.2651266  and Variance: 0.01419242 



