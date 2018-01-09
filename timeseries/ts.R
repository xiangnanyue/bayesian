#' read from csv the time series data
#' at first we deal with the 5 min data from the last 2 months

datas <- read.csv(file = "~/Documents/Paris_M2/DataScience/Machine_Learning/timeseries/easymoney_fndata_5min.csv")
datas$Date <- as.POSIXct(datas$time, tz="", format="%Y-%m-%d %H:%M")
par(mfrow = c(1,1))
# plot(start_value ~ Date, datas, type = "l")
plot(datas[,2],ylab = "index of HS300", type = "l")

#' draw the graph of white noise 
#' take N = 1400, eg.
par(mfrow = c(2,1))
w = rnorm(1400,0,1)
plot.ts(w, ylab="white noise")

#' draw the moving-average graph
#' side = 2, use symmetric data, centred on the current data point
ma = filter(x = w,sides = 2,rep(1/20, 20))
plot.ts(ma, ylab="moving average")

#' draw the autoregression graph
par(mfrow = c(2,1))
ar = filter(x = w,filter = c(1, -0.9), method = "recursive")
plot.ts(w, ylab="white noise")
plot.ts(ar, ylab="autoregressive")

#' calculate the acf and ccf for datas
par(mfrow = c(1,1))
acf(datas$end_value, 72)
ccf(datas$trade_volume, datas$end_value, 72)

#' ar intro


#' the following specifies the steps for ts analysis arima
# first load the relevant packages
# auto.arima package
install.packages('ggplot2')
install.packages('forecast')
install.packages('tseries')
library('ggplot2')
library('forecast')
library('tseries')

datas$index = 1:length(datas$end_value)
ggplot() +
  geom_line(data = datas, aes(x = index, y = end_value, colour = "end value" ) ) +
  geom_line(data = datas, aes(x = index, y = start_value, colour = "start value")) +
  geom_line(data = datas, aes(x = index, y = highest, colour = "highest")) +
  geom_line(data = datas, aes(x = index, y = lowest, colour = "lowest"))

ggplot() +
  geom_line(data = datas, aes(x = index, y = trade_volume, colour = "trade volume" ) )
#' clean the data 
#' tsclean() as part of its forecast package. 
#' tsclean() identifies and replaces outliers using 
#' series smoothing and decomposition.

#tsclean()
datas$end_value = tsclean(datas$end_value)
#' draw --

datas$index = 1:length(datas$end_value)
datas$ma12_ev = ma(datas$end_value, order=12) # averaging in one hour
datas$ma48_ev = ma(datas$end_value, order=48) # averaging within one day = 4h
ggplot() +
  geom_line(data = datas, aes(x = index, y = end_value, colour = "real value" ) ) +
  geom_line(data = datas, aes(x = index, y = ma12_ev, colour = "one hour averaging")) +
  geom_line(data = datas, aes(x = index, y = ma48_ev, colour = "one day averaging")) 

#' data examination for seasonarity
decomp = stl(na.omit(datas$ma12_ev), s.window="periodic")


#' data examination for stationarity
adf.test(datas$end_value, alternative = "stationary")

#' first order difference
#' xt - xt-1
ev_d1 = diff(datas$end_value, differences = 1)
index_d1 = 1:length(ev_d1)
plot(index_d1, ev_d1, type = "l")
adf.test(ev_d1, alternative = "stationary")

#' find order of p, d q by acf and pacf
par(mfrow = c(1,2))
Acf(ev_d1, lag.max = 200, main='ACF for Differenced Series')
Pacf(ev_d1, lag.max = 200, main='PACF for Differenced Series')

#' Model 1
#' use auto.arima(ev_d1, seasonal=FALSE)
auto.arima(datas$end_value, seasonal=FALSE)
fit<-auto.arima(datas$end_value, seasonal=FALSE)
tsdisplay(residuals(fit), lag.max=100, main='(1,1,0) Model Residuals')

#' Model 1.1
auto.arima(datas$end_value, seasonal=FALSE, max.p=50, max.q=50)
fit1_1<-auto.arima(datas$end_value, seasonal=FALSE)
tsdisplay(residuals(fit), lag.max=160, main='(1,1,0) Model Residuals')


#' Model 2
#' some strong correlation exits at interval = 40, this may be 
#' seasonal because 40 ~ one day
fit2 = arima(datas$end_value, order=c(40,1,0))
tsdisplay(residuals(fit2), lag.max=180, main='(40,1,0) Model Residuals')
#' if we decided to take this model, then 
fcast <- forecast(fit2, h=130)
plot(fcast)

#' Model 3
fit3 = arima(datas$end_value, order=c(1,1,48))
tsdisplay(residuals(fit3), lag.max=80, main='(1,1,48) Model Residuals')
#' if we decided to take this model, then 
fcast <- forecast(fit3, h=230)
plot(fcast)

#' backtest p = 48
l = length(datas$end_value)
hold <- window(ts(datas$end_value), start=1000)
fit_no_holdout = arima(ts(datas$end_value[-c(1000:l)]), order=c(48,1,0))
fcast_no_holdout <- forecast(fit_no_holdout,h=80)
plot(fcast_no_holdout, main=" ")
lines(ts(datas$end_value))

tsdisplay(residuals(fit_no_holdout), lag.max=180, main='(48,1,0) Model Residuals')


#' backtest p = 90
l = length(datas$end_value)
hold <- window(ts(datas$end_value), start=1000)
fit_no_holdout2 = arima(ts(datas$end_value[-c(1000:l)]), order=c(1,1,0), seasonal = list(order = c(0,1,1), period=48))
fcast_no_holdout2 <- forecast(fit_no_holdout2,h=200)
plot(fcast_no_holdout2, main=" ")
lines(ts(datas$end_value))
tsdisplay(residuals(fit_no_holdout2), lag.max=100, main='(1,1,0, 0, 1, 1, 48) Model Residuals')

#' backtest p = 90
l = length(datas$end_value)
hold <- window(ts(datas$end_value), start=1000)
fit_no_holdout3 = arima(ts(datas$end_value[-c(1000:l)]), order=c(1,1,0), seasonal = list(order = c(0,1,1), period=48))
fcast_no_holdout3 <- forecast(fit_no_holdout3,h=200)
plot(fcast_no_holdout3, main=" ")
lines(ts(datas$end_value))
tsdisplay(residuals(fit_no_holdout3), lag.max=100, main='(1,1,0, 0, 1, 1, 48) Model Residuals')

#' now the art: sarima

##################
#' however the trade volume is something easier to be predicted
l = length(datas$trade_volume)
hold <- window(ts(datas$trade_volume), start=1000)
auto.arima(datas$trade_volume, seasonal = TRUE)
par(mfrow = c(1,1))
fit_volume = arima(ts(datas$trade_volume[-c(1000:l)]), order=c(1,0,1),seasonal = list(order = c(1,0,1), period=48) )
fit_volume <- forecast(fit_volume,h=200)
plot(fit_volume, main=" ")
lines(ts(datas$trade_volume))
#' standard method 
#tv_d1 = diff(datas$trade_volume, differences = 1)
#par(mfrow = c(1,2))
#Acf(tv_d1, lag.max = 200, main='ACF for Differenced Series')
#Pacf(tv_d1, lag.max = 200, main='PACF for Differenced Series')

# # trend
# datas$trade_ma = ma(datas$trade_volume, order=144, centre = T)
# par(mfrow = c(1, 1))
# plot(as.ts(datas$trade_volume))
# lines(datas$trade_ma)
# # detrend 
# detrend_volume = datas$trade_volume - datas$trade_ma
# plot(as.ts(detrend_volume))
# 
# # multiply
# detrend_volume2 = datas$trade_volume / datas$trade_ma
# plot(as.ts(detrend_volume2))
# 
# # no season ...
# count_cn = ts(na.omit(detrend_volume2),frequency=200 )
# decomp <- stl(na.omit(detrend_volume2), s.window="periodic")

# #
# par(mfrow = c(2,1))
# Acf(na.omit(detrend_volume2), lag.max = 200, main='ACF for Differenced Series')
# Pacf(na.omit(detrend_volume2), lag.max = 200, main='PACF for Differenced Series')

datas$trade_ma = ma(datas$trade_volume, order=12, centre = T)
count_cn = ts(na.omit(datas$trade_ma),frequency=144 )
plot(count_cn)
decomp <- stl(na.omit(count_cn), s.window="periodic" )
plot(decomp)
deseasonal_cnt <- seasadj(decomp)
plot(deseasonal_cnt)

# result is that is not stationary 
#adf.test(na.omit(datas$trade_ma), alternative = "stationary")
Acf(na.omit(deseasonal_cnt), lag.max = 200, main='ACF for Differenced Series')
Pacf(na.omit(deseasonal_cnt), lag.max = 200, main='PACF for Differenced Series')
adf.test(deseasonal_cnt, alternative = "stationary")

# now do diff
tv_d1 = diff(deseasonal_cnt, differences = 1)
index_d1 = 1:length(tv_d1)
plot(index_d1, tv_d1, type = "l")
adf.test(tv_d1, alternative = "stationary")

# now the new ACF
par(mfrow = c(1,2))
Acf(na.omit(tv_d1), lag.max = 50, main='ACF for Differenced Series')
Pacf(na.omit(tv_d1), lag.max = 50, main='PACF for Differenced Series')


# to be lazy, 112
fit_volume1 <- auto.arima(deseasonal_cnt, seasonal = FALSE)
tsdisplay(residuals(fit_volume1), lag.max=100, main='(1,1,2) Model Residuals')

# to be lazy but with seasonal
fit_w_seasonality = auto.arima(deseasonal_cnt, seasonal=TRUE)
fit_w_seasonality

fit_imp <- arima(deseasonal_cnt, order = c(1,1,0))
tsdisplay(residuals(fit_imp), lag.max=200, main='(1,1,1) Model Residuals')

#' PACF cut off at about 48, try the AR(48), which has been done before
#' PAF cut off at about 12, this means we should add a MA(12)
fit_imp <- arima(deseasonal_cnt, order = c(1,1,12))
tsdisplay(residuals(fit_imp), lag.max=200, main='(1,1,12) Model Residuals')

fit_imp2 <- arima(deseasonal_cnt, order = c(1,2,24))
tsdisplay(residuals(fit_imp2), lag.max=200, main='(1,2,24) Model Residuals')
