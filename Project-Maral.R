#########################
# Step 1
########################

# read the data from the file
data <- read.csv("F:/Tampere Uni/Time Series/Project/Case_study.csv")

time_series <- data$V60 # my own column


###### 1)
# Plot the time series
plot(ts(time_series), ylab='Time series', xlab = 'Lag', type='o')


###### 2 a)

# Calculate differences (∇Yt and ∇2Yt)
diff1_Yt <- diff(time_series,1)
diff2_Yt <- diff(diff1_Yt)

# Plot the autocorrelation functions
acf( time_series , main ="ACF for Yt", lag.max = 40)
# for ∇Yt
acf( diff1_Yt , main ="ACF for ∇Yt", lag.max = 40)
# for ∇Yt
acf( diff2_Yt , main ="ACF for ∇2Yt", lag.max = 40)


###### 2 b)
# Augmented Dickey - Fuller Unit Root Test
#install.packages("urca")
library ('urca')
summary (ur.df( time_series , type = "drift", lags = 1, selectlags = "Fixed"))


#install.packages("forecast")
#library(forecast)
#ndiffs(diff1_Yt)

######## 3
par(mfrow = c(1, 2))
# for ∇Yt
acf(diff1_Yt , main ="ACF for ∇Yt", lag.max = 40) # qmax=4 

pacf(diff1_Yt , main ="PACF for ∇Yt", lag.max = 40) #pmax=2


######################################
#     STEP 2     
######################################

####### 1 & 2)

# Set up pmax and qmax
pmax <- 4
qmax <- 4

# Create a matrix to store AIC and BIC values
AIC_results <- matrix(NA, nrow = pmax , ncol = qmax)
BIC_results <- matrix(NA, nrow = pmax , ncol = qmax)

# Loop over possible (p, q) combinations
for (p in 1:pmax) {
  for (q in 1:qmax) {
    
    # Fit ARIMA model
    model <- arima(time_series, order = c(p, 1, q))
    
    # Compute AIC and BIC
    aic <- AIC(model)
    bic <- BIC(model)
    
    
    # Store AIC and BIC values in the matrix
    AIC_results[p, q] <- aic
    BIC_results[p, q] <- bic
  }
}

# Display the AIC and BIC values for each model
AIC_results   # min AIC p=3, q=1
BIC_results   # min BIC p=3, q=1


# Flatten the AIC results into a vector
AIC_vector <- as.vector(AIC_results)

# Find the indices of the three minimum AIC values
min_aic_indices <- order(AIC_vector)[1:3]
min_aic_indices

# Display the best three models based on AIC
best_models <- matrix(NA, nrow = 3, ncol = 2)
for (i in 1:3) {
  col_index <- floor((min_aic_indices[i] - 1) / ncol(AIC_results)) + 1
  row_index <- (min_aic_indices[i] - 1) %% ncol(AIC_results) + 1
  best_models[i, ] <- c(row_index, col_index)  # Subtract 1 to get the actual p and q values
}

print("Best three models based on AIC:")
print(best_models)


# Flatten the BIC results into a vector
BIC_vector <- as.vector(BIC_results)

# Find the indices of the three minimum BIC values
min_bic_indices <- order(BIC_vector)[1:3]
min_bic_indices

# Display the best three models based on BIC
best_models_bic <- matrix(NA, nrow = 3, ncol = 2)
for (i in 1:3) {
  col_index <- floor((min_bic_indices[i] - 1) / ncol(BIC_results)) + 1
  row_index <- (min_bic_indices[i] - 1) %% ncol(BIC_results) + 1
  best_models_bic[i, ] <- c(row_index, col_index)  # get the actual p and q values
}

print("Best three models based on BIC:")
print(best_models_bic)



######################################
#    STEP 3
######################################

###### 1)

# Fit ARIMA models
model.1 <- arima(time_series, order = c(3, 1, 1))
model.2 <- arima(time_series, order = c(4, 1, 1))
model.3 <- arima(time_series, order = c(3, 1, 2))


# Ljung-Box test
ljung_box_test <- Box.test(residuals(model.1), lag = 10, type = "Ljung-Box")
ljung_box_test  # p-value = 0.9817

ljung_box_test <- Box.test(residuals(model.2), lag = 10, type = "Ljung-Box")
ljung_box_test  # p-value = 0.9826

ljung_box_test <- Box.test(residuals(model.3), lag = 10, type = "Ljung-Box")
ljung_box_test  # p-value = 0.9827
# If the p-value is large (typically > 0.05), it suggests that there
# is no evidence of autocorrelations in the residuals.

windows(width = 8, height = 6)
# Plot ACF and PACF of residuals
par(mfrow=c(3,2)) 

# ACF plot of residuals
acf(residuals(model.1), lag.max = 40, main = "ACF of Residuals-model1") 

# PACF plot of residuals
pacf(residuals(model.1), lag.max = 40, main = "PACF of Residuals-model1")

# ACF plot of residuals
acf(residuals(model.2), lag.max = 40, main = "ACF of Residuals-model2") 

# PACF plot of residuals
pacf(residuals(model.2), lag.max = 40, main = "PACF of Residuals-model2")

# ACF plot of residuals
acf(residuals(model.3), lag.max = 40, main = "ACF of Residuals-model3") 

# PACF plot of residuals
pacf(residuals(model.3), lag.max = 40, main = "PACF of Residuals-model3")

par(mfrow=c(1,1))  # Reset to a single plot
dev.off()

##### 2) 

windows(width = 8, height = 6)
par(mfrow=c(3,2))
# Histogram of residuals
hist(residuals(model.1), main = "Histogram of Residuals 1", xlab = "Residuals")

# QQ plot of residuals
qqnorm(residuals(model.1), main = "QQ Plot of Residuals 1")
qqline(residuals(model.1))

hist(residuals(model.2), main = "Histogram of Residuals 2", xlab = "Residuals")

qqnorm(residuals(model.2), main = "QQ Plot of Residuals 2")
qqline(residuals(model.2))

hist(residuals(model.3), main = "Histogram of Residuals 3", xlab = "Residuals")

qqnorm(residuals(model.3), main = "QQ Plot of Residuals 3")
qqline(residuals(model.3))

# Shapiro-Wilk test for normality
shapiro_test <- shapiro.test(residuals(model.1))
shapiro_test  #p-value = 0.8755

shapiro_test <- shapiro.test(residuals(model.2))
shapiro_test  #p-value = 0.902

shapiro_test <- shapiro.test(residuals(model.3))
shapiro_test  #p-value = 0.9065


#### 4)

# Generate the best-fitted values from the ARIMA model(3,1,1)
library(forecast)
fitted_values1 <- fitted(model.1)

par(mfrow=c(1,1))

# Plot original time series
plot(time_series, type = "l", col = "blue", lwd = 2, main = "Original Time Series vs. Fitted Model",
     xlab = "Lag", ylab = "Value",ylim = range(c(time_series, fitted_values1)))

# Add best-fitted values to the plot
lines(fitted_values1, col = "red", lty = 2, lwd = 2)


# Add legend
legend("topright", legend = c("Original Time Series", "Fitted Model"),
       col = c("blue", "red"), lty = 1:2, cex = 0.7, lwd = 2)



##############################
#   STEP 4
##############################

# Forecast h=10 steps ahead
#forecast_10 <- forecast(model.1, h=25)
forecast_10 <- predict(model.1, n.ahead = 10)

# Forecast h=25 steps ahead
#forecast_25 <- forecast(model.1, h=25)
forecast_25 <- predict(model.1, n.ahead = 25)


windows(width = 8, height = 6)
par(mfrow=c(2,1))

x_axis <- seq_along(time_series)
# Plot original time series and forecast
#plot(time_series, type = "l", xlim = c(1, length(time_series) + 25), ylim = range(time_series, forecast_25$pred))
plot(time_series, type = "l", xlim = c(1, length(time_series) + 25), ylim = range(-20, 2))
lines(x_axis[length(time_series)] + (1:10), forecast_10$pred, lty = 2, col = "blue", lwd = 2)
lines(x_axis[length(time_series)] + (1:10), forecast_10$pred - 1.96 * forecast_10$se, lty = 3, col = "red", lwd = 2)
lines(x_axis[length(time_series)] + (1:10), forecast_10$pred + 1.96 * forecast_10$se, lty = 3, col = "orange", lwd = 2)

# Add a legend to the plot
legend("bottomleft", 
       legend=c("Original", "Lower", "Upper"), 
       col=c("blue", "red", "orange"), lty=3, cex=0.7,lwd = 2)

plot(time_series, type = "l", xlim = c(1, length(time_series) + 25), ylim = range(-20, 4))
lines(x_axis[length(time_series)] + (1:25), forecast_25$pred, lty = 2, col = "blue", lwd = 2)
lines(x_axis[length(time_series)] + (1:25), forecast_25$pred - 1.96 * forecast_25$se, lty = 3, col = "red", lwd = 2)
lines(x_axis[length(time_series)] + (1:25), forecast_25$pred + 1.96 * forecast_25$se, lty = 3, col = "orange", lwd = 2)

# Add a legend to the plot
legend("bottomleft", 
       legend=c("Original", "Lower", "Upper"), 
       col=c("blue", "red", "orange"), lty=3, cex=0.7,lwd = 2)
