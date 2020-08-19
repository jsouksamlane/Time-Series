library(astsa)
library(dplyr)
library(MASS)
library(MuMIn)
library(ggplot2)

#DATA EXPLORATION:

# Reading in the data set
austin <- read.csv("C://Users//Jalen//Desktop//PSTAT174//austin_weather.csv")

HTemp <- austin$Temp.high.f # Recorded highest temperatures per day
LTemp <- austin$Temp.low.f # Recorded lowest temeratures per day 
ATemp <- austin$Temp.avg.f # Avg temperature per day

# Overall avg of each temp
mean(HTemp) # 80.86
mean(LTemp) # 59.9
mean(ATemp) # 70.64

# Median of each temp
median(HTemp) # 83
median(LTemp) # 63
median(ATemp) # 73

# Mode of each temp
get_mode <- function(x) {
  uniqx <- unique(x)
  uniqx[which.max(tabulate(match(x, uniqx)))]
}
get_mode(HTemp) # 86
get_mode(LTemp) # 75
get_mode(ATemp) # 84

# Max of each temp
m1 <- max(HTemp) # 107
m2 <- max(LTemp) # 81
m3 <- max(ATemp) # 93

# Min of each temp
n1 <- min(HTemp) # 32
n2 <- min(LTemp) # 19
n3 <- min(ATemp) # 29

# Range of each temp
(m1-n1) # 75
(m2-n2) # 62
(m3-n3) # 64

LowTemp.pgram <- mvspec(austin$Temp.low.f, log= 'no')
LowTemp.pgram <- mvspec(austin$Temp.low.f, log= 'no', xlim =c(0.0028, 0.003))

#0.002963
#1/0.00296 = 337.5 days that the low temperatures repeat
# This means that it occus slightly less than a year each clycle

#ANALYSIS (Low Temperature):

#from our initial observations, we see that highest temperatures haven't really been growing that much, but low tempertures during the winter periods have been varying (more dense between days 1-700 than the days that come afterwards). Therefore, we'll only look at the austin$Temp.avg.f

#low temperatures plot

plot(austin$Temp.low.f)

LowTemp <- austin$Temp.low.f

#initial ts plot without transformation
tsplot(LowTemp)

#from the plot, it looks seasonal. It also looks as if there are more colder temperature measurements during the winter period (between 1-700 days) than there are during the winter periods past day 700-ish. This implies that temperatures are gradually getting warmer.

#Since the plot is seasonal, it is non-stationary. Let's look at a couple of transformation methods for the data:

#a) log-transform

LowTemp_log <- log(LowTemp)
tsplot(LowTemp_log)

#plot still looks seasonal

#b) square root transform

LowTemp_sqrt <- sqrt(LowTemp)
tsplot(LowTemp_sqrt)

#plot still looks seasonal

#In order to properly remove seasonality from a time series plot, we perform the differencing function on the data.

LowTemp_diff1 <- diff(LowTemp, differences = 1)

tsplot(LowTemp_diff1, main = "De-seasonalized Time Series")

#That's better! The tsplot looks more stationary after the first difference. Most data are stationary after the first difference, aside from financial math data (stocks, for instance) that require 2 or more differences

#Hypothesis: We project that there will be less colder temperature measurements every day (during the winter season). After the end of 2015

#Now let's look at the differenced data's ACF and PACF plots in order to detect an appropriate model

acf2(LowTemp_diff1)

#ACF plot drops to near zero after lag 4. PACF plot appears to have a more geometric decay. Therefore, we deduce MA(4) to be a possible model

#estimate the model using ar()

fit.ar <- ar(LowTemp_diff1, method="yw") #Best model is AR(13)

aiccs <- matrix(NA, nr=11, nc=11)
dimnames(aiccs) = list(p=0:10, q=0:10)
for(p in 0:10){
  for(q in 0:10){
    aiccs[p+1,q+1] =AICc(arima(diff(LowTemp_diff1), order = c(p,0,q), method = 'ML'))
  }
}

aiccs
(aiccs == min(aiccs))

#grabs the lowest aic from the different combinations of ARMA(p,q) models. Outputted ARMA(6,7) as the best model

#candidate models: MA(4), AR(13), ARMA(6,7)

#Perform diagnostics: sarima

#In order to use sarima(), we need to find P, D and Q from the seasonal data's ACF and PACF. 

acf2(LowTemp)

#ACF: gradual decay, PACF: significant until lag 1. Best model: AR(1)


MA_4 <- sarima(LowTemp_diff1, 0, 1, 4, P=1, D=0,Q=0, S=1)
AR_13 <- sarima(LowTemp_diff1, 13, 1, 0, P=1, D=0,Q=0, S=1)
ARMA <- sarima(LowTemp_diff1, 6, 1, 7, P=1, D=0,Q=0, S=1)


#AIC
MA_4$AIC #6.300727
AR_13$AIC #6.436217
ARMA$AIC #6.298166

#BIC
MA_4$BIC #6.328276
AR_13$BIC #6.499186
ARMA$BIC #6.361135


#Conclusion: ARMA(6,7) is the best model due to lowest AIC and BIC scores. MA(4) could be a good alternative

#Normality testing
fit1  <- arima(LowTemp_diff1, order=c(6,1,7), method="ML") #ARMA(6,7)

shapiro.test(residuals(fit1)) #ARMA(6,7)

fit2 <- arima(LowTemp_diff1, order=c(0,1,4), method="ML") #MA(4)

shapiro.test(residuals(fit2)) #MA(4)


#Both shapiro p-values are well below 0.05. Therefore, these models doesn't abide to a normal distribution

#box.test() - independence testing

Box.test(residuals(fit1), type="Ljung")
Box.test(residuals(fit2), type="Ljung")


#Ljung test p-value for  MA(4) is 0.6131. For ARMA(6,7) it's 0.9972. If the p value is greater than 0.05, then the residuals are independent which we want for the model to be correct. Hence ARMA(6,7) has a p-value greater than MA(4), we prefer it over the second model.

par(mfrow=c(1,2),oma=c(0,0,2,0))

# Histogram
hist(residuals(fit2),main = "Histogram")

# q-q plot
qqnorm(residuals(fit3))
qqline(residuals(fit3),col ="blue")

#SARIMA prediction

astsa::sarima.for(LowTemp, n.ahead = 10, p = 6, d = 1, q = 7, P = 1, D = 0, Q = 0, S =1)

#P, D, Q for seasonal plot
#p, d, q for non-seasonal plot

#ANALYSIS: diurnal range 

#Hypothesis: The difference between high temperatures and low temperature is slowly decreasing over the years. Since variation in high temperatures is less prominent than the low temperatures, then we can assume for this hypothesis that colder temperatures during the winter periods are rising (causing the difference to decrease over time).

#For this section, we created a separate column that takes the difference between the high and low temperatures

temp_difference <- austin$Temp.high.f - austin$Temp.low.f
length <- length(names(austin))
austin[length+1] <- temp_difference
names(austin) <- c(names(austin)[1:length], 'Temp_Difference')


write.csv(austin, "C://Users//Jalen//Desktop//PSTAT174//new_austin_weather.csv", row.names = FALSE)

new_austin <- read.csv("C://Users//Jalen//Desktop//PSTAT174//new_austin_weather.csv")

diurnal <- (new_austin$Temp_Difference)

#plot timeseries 

tsplot(diurnal)

#looks non-stationary; has seasonality. Use differnce to get rid of seasonality

diurnal.diff <- diff(diurnal)

#plot timeseries 
tsplot(diurnal.diff)

#look more stationary now

#Look at acf and pacf of difference

acf2(diurnal.diff)

#Looks like an MA(3) Process or an ARMA model

diurnal.fit.MA3 <- arima(diurnal.diff, order=c(0,1,3), method="ML")

#find ARMA model

aiccs1 <- matrix(NA, nr=11, nc=11)
dimnames(aiccs1) = list(p=0:10, q=0:10)
for(p in 0:10){
  for(q in 0:10){
    aiccs1[p+1,q+1] =AICc(arima(diff(diurnal.diff), order = c(p,0,q), method = 'ML'))
  }
}

aiccs1
(aiccs1 == min(aiccs1))

#Fit an ARMA(6,6) (Best model found from AIC by look from 0-10)

diurnal.fit.ARMA66 <- arima(diurnal.diff, order=c(6,1,6), method="ML")

shapiro.test(residuals(diurnal.fit.MA3)) 
shapiro.test(residuals(diurnal.fit.ARMA66))

#both p-values well below 0.05; hence, it does not abide to normality

sarimafitMA <- sarima(diurnal.diff, p =0, d= 1, q =3, P = 3, D= 1, Q= 0, S=1) #AIC 6.732 BIC:6.76
sarimafitARMA <- sarima(diurnal.diff, p =6, d= 1, q =6, P = 3, D= 1, Q= 0, S=1) #AIC:6.637 BIC:6.7

#Pretty close prediction

##predict with sarimaARMA
sarima.for(diurnal, n.ahead = 10, p = 6, d = 1, q = 6, P = 3, D = 0, Q = 0, S =1)

#predicts that difference between high and low temperatures are getting smaller