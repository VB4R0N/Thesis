################### DAX ANALYSIS ###################
#### 0 Load relevant packages ----
library(tidyverse)
library(rugarch)
library(quantmod)
library(TTR)
library(DescTools)
library(gridExtra)
library(lubridate)
library(scales)
library(sandwich)
library(tseries)
library(lmtest)
library(forecast)
library(plotly)
library(grid)
library(RColorBrewer)
library(patchwork)
library(RColorBrewer)
library(car)
library(imputeTS)
library(Metrics)
library(aTSA)
#### 1 Import and Gather Data ----
rm(list = ls())
set.seed(1892) # Ensure identical results for you 
##### 1.1 Weekly DAX Closing Prices -----
dax<- data.frame(
  getSymbols(
    Symbols = "^GDAXI", 
    src = "yahoo",
    auto.assign = FALSE,
    from="1990-01-02",
    to="2024-12-05",
    periodicity = "daily")) %>%
  mutate(Date=as.Date(rownames(.))) %>% 
  rename(P=GDAXI.Close) %>% 
  select(c(Date, P)) %>% 
  mutate(day=weekdays(Date)) %>% 
  filter(day %in% "Mittwoch") %>% 
  mutate(r=100*ROC(P, type = "continuous")) %>% 
  na.omit()
# Note: I chose Wednesday, as it requires fewest imputations 

##### 1.2 Volatility Index -----
vix<- data.frame(
  getSymbols(
    Symbols = "^VIX", 
    src = "yahoo",
    auto.assign = FALSE,
    from="1990-01-02",
    to="2024-12-05",
    periodicity = "daily")) %>%
  mutate(Date=as.Date(rownames(.))) %>% 
  mutate(day=weekdays(Date)) %>% 
  filter(day %in% "Mittwoch") %>% 
  rename(P=VIX.Close) %>% 
  select(c(Date, P)) %>% 
  na.omit()

##### 1.3 Time series imputation of missing values -----
###### 1.3.1 DAX ------
complete_dates <- data.frame(Date = seq(
  from = min(dax$Date), 
  to = max(dax$Date), 
  by = "7 days"
))
dax_full <- full_join(complete_dates, dax, by = "Date") %>%
  arrange(Date)

dax_full<- dax_full %>%   
  mutate(P=na_interpolation(dax_full$P, option = "linear")) %>% 
  mutate(r=100*ROC(P, type = "continuous")) %>% 
  select(-day) %>% 
  na.omit()

###### 1.3.2 VIX ------
vix_full <- full_join(complete_dates, vix, by = "Date") %>%
  arrange(Date)

vix_full<- vix_full %>%   
  mutate(P=na_interpolation(vix_full$P, option = "linear")) %>% 
  na.omit()

# Create two lags for VIX (will be relevant in training set)
vix_lags<-as.data.frame(sapply(1:2, function(i) lag(vix_full$P, i))) %>% 
  rename(vix_t1=V1, vix_t2=V2) %>% 
  na.omit()
dax_full<-dax_full %>% mutate(vix_t1=vix_lags$vix_t1, vix_t2=vix_lags$vix_t2)
#### 2 Training set ----
training_dax<-dax_full[1:round(0.65*length(dax_full$r), digits = 0),]
##### 2.1. Check presence of GARCH effects -----
###### 2.1.1 Recreate correlograms in Figure 8 -----
# ACF of returns
plot_r <- ggAcf(training_dax$r, lag.max = 26, plot = FALSE) %>%
  autoplot() +
  labs(title = "ACF of DAX Returns", x = "Lags", y = "ACF") +
  theme_minimal(base_size = 14) +
  theme(plot.title = element_text(face = "bold"))

# ACF of squared returns
plot_r_squared <- ggAcf(training_dax$r^2, lag.max = 26, plot = FALSE) %>%
  autoplot() +
  labs(title = "ACF of Squared DAX Returns", x = "Lags", y = "ACF") +
  theme_minimal(base_size = 14) +
  theme(plot.title = element_text(face = "bold")) 

# Combine the two plots
plot_r + plot_r_squared
###### 2.2.2 Perform Lagrange Multiplier Test
est<-estimate(training_dax$r, PDQ = c(0,0,0), intercept = FALSE, output = FALSE)
arch.test(est)
##### 2.3 GARCH-Type Best Fits -----
###### 2.3.1 GARCH Model ------
garchpq <- data.frame(p = integer(), q = integer(), AIC = numeric())
for (i in 1:5) {
  for (j in 1:2) {
    spec_garch <- ugarchspec(
      variance.model = list(model = "sGARCH", garchOrder = c(i, j)),
      mean.model = list(armaOrder = c(0, 0), include.mean=FALSE),
      distribution.model = "norm"
    )
    
    fit_garch <- ugarchfit(spec = spec_garch, data = training_dax$r, silent = TRUE)
    
    aic_value <- infocriteria(fit_garch)[1]
    
    garchpq <- rbind(garchpq, data.frame(p = i, q = j, AIC = aic_value))
  }
}
row_min_garch<-which.min(garchpq$AIC)
garch <- ugarchfit(spec = ugarchspec(
  variance.model = list(model = "sGARCH", garchOrder = c(garchpq[row_min_garch,1], garchpq[row_min_garch,2])),
  mean.model = list(armaOrder = c(0, 0), include.mean=FALSE),
  distribution.model = "norm"),
  data = training_dax$r, silent = TRUE)

###### 2.3.2 EGARCH Model ------
egarchpq <- data.frame(p = integer(), q = integer(), AIC = numeric())

for (i in 1:5) {
  for (j in 1:2) {
    spec_garch <- ugarchspec(
      variance.model = list(model = "eGARCH", garchOrder = c(i, j), variance.targeting=TRUE),
      mean.model = list(armaOrder = c(0, 0), include.mean=FALSE),
      distribution.model = "norm"
    )
    
    
    fit_garch <- ugarchfit(spec = spec_garch, data = training_dax$r, silent = TRUE, solver = "hybrid")
    
    aic_value <- infocriteria(fit_garch)[1]
    
    
    egarchpq <- rbind(egarchpq, data.frame(p = i, q = j, AIC = aic_value))
  }
}
row_min_egarch<-which.min(egarchpq$AIC)

egarch <- ugarchfit(spec = ugarchspec(
  variance.model = list(model = "eGARCH", garchOrder = c(garchpq[row_min_egarch,1], garchpq[row_min_egarch,2]), variance.targeting=TRUE),
  mean.model = list(armaOrder = c(0, 0), include.mean=FALSE),
  distribution.model = "norm"),
  data = training_dax$r, silent = TRUE, solver="hybrid")

###### 2.3.3 GARCHX Model ------
# Correlogram of VIX (Figure 7):
tsdisplay(vix_full$P[1:length(training_dax$r)])
# -> Suggests AR(2) structure -> Include k \in [1,2]
garchxpqr <- data.frame(p = integer(), q = integer(), r = integer(), AIC = numeric())

for (i in 1:5) { 
  for (j in 1:2) { 
    for (k in 1:2) { 
      
      exog_vars <- as.matrix(training_dax[,(4:5)])[, 1:k, drop = FALSE]
      
      spec_garchx <- ugarchspec(
        variance.model = list(
          model = "sGARCH", 
          garchOrder = c(i, j),
          external.regressors = exog_vars
        ),
        mean.model = list(armaOrder = c(0, 0), include.mean=FALSE), 
        distribution.model = "norm"
      )
      
      
      tryCatch({
        fit_garchx <- ugarchfit(spec = spec_garchx, data = training_dax$r, silent = TRUE)
        
        aic_value <- infocriteria(fit_garchx)[1]
        
        garchxpqr <- rbind(garchxpqr, data.frame(p = i, q = j, r = k, AIC = aic_value))
      }, error = function(e) {

        garchxpqr <- rbind(garchxpqr, data.frame(p = i, q = j, r = k, AIC = NA))
      })
      
    }
  }
}

garchx<-ugarchfit(
  spec = ugarchspec(
    variance.model = list(
      model = "sGARCH", 
      garchOrder = c(5, 1),
      external.regressors = as.matrix(training_dax[, 4])
    ),
    mean.model = list(armaOrder = c(0, 0), include.mean=FALSE),
    distribution.model = "norm"
  ),
  data = training_dax$r,  
  silent = TRUE
)

#### 3  RV Proxy Using Daily Returns ----
##### 3.1 Import Daily Returns -----
dax_daily<- data.frame(
  getSymbols(
    Symbols = "^GDAXI", 
    src = "yahoo",
    auto.assign = FALSE,
    from="1990-01-02",
    to="2024-12-05",
    periodicity = "daily")) %>%
  mutate(Date=as.Date(rownames(.))) %>% 
  rename(P=GDAXI.Close) %>% 
  select(c(Date, P)) %>% 
  mutate(day=weekdays(Date)) %>% 
  mutate(r=100*ROC(P, type = "continuous")) %>% 
  na.omit()

##### 3.2 Time series imputation of missing values  -----
# I do this to retain n=5 across time
complete_dates_daily <- data.frame(Date = seq(
  from = min(dax_daily$Date), 
  to = max(dax_daily$Date), 
  by = "1 day"
)) %>% 
  mutate(day=weekdays(Date)) %>% 
  filter(!day %in% c("Samstag", "Sonntag")) %>% 
  select(-day)


dax_full_daily <- full_join(complete_dates_daily, dax_daily, by = "Date") %>%
  arrange(Date)

dax_full_daily<- dax_full_daily %>%   
  mutate(P=na_interpolation(dax_full_daily$P, option = "linear")) %>% 
  mutate(r=100*ROC(P, type = "continuous")) %>% 
  mutate(day=weekdays(Date)) %>% 
  na.omit()


##### 3.3 Compute RV -----
starting_rv <- seq(1, length(dax_full_daily$r) - 1, by = 5)
ending_rv <- seq(5, length(dax_full_daily$r), by=5)
rv <- rep(NA, length(starting_rv))


for (i in 1:length(rv)) {
  rv[i] <- sum(dax_full_daily$r[starting_rv[i]:ending_rv[i]]^2)
}
rv<-rv[-1] # first rv value not in dax_full dataset
tail(dax_full$rv)
# Merge RV with dax_full
dax_full<- dax_full %>% mutate(rv=rv)

##### 3.4 Display in-Sample Volatility (Figure 3) -----
colors<-colors <- c("cyan", "orange", "forestgreen")
data.frame(
  GARCH = garch@fit$sigma,
  GARCHX = garchx@fit$sigma,
  EGARCH = egarch@fit$sigma,
  RV = sqrt(dax_full$rv[1:1184]),
  Date = dax_full$Date[1:1184]
) %>%
  ggplot(aes(x = Date)) +
  geom_line(aes(y = RV, color = "Realized Volatility"), linewidth = 0.5, alpha = 0.6) +
  geom_line(aes(y = GARCH, color = "GARCH Model"), linewidth = 0.75) +
  geom_line(aes(y = EGARCH, color = "EGARCH Model"), linewidth = 0.75) +
  geom_line(aes(y = GARCHX, color = "GARCHX Model"), linewidth = 0.75) +
  scale_color_manual(
    values = c(
      "Realized Volatility" = "grey26",  
      "GARCHX Model" = colors[3],     
      "EGARCH Model" = colors[2],     
      "GARCH Model" = colors[1]    
    )
  ) +
  labs(
    x = "Time",
    y = "Volatility (Standard Deviation)",
    title = "Comparison of Volatility Models: GARCH vs. GARCHX vs. EGARCH",
    color = "Model"
  ) +
  theme_minimal() +
  theme(
    plot.title = element_text(hjust = 0.5, size = 20, face = "bold"),
    axis.text = element_text(size = 12),
    axis.title = element_text(size = 18),
    legend.position = "right",
    legend.title = element_text(face = "bold", size=15),
    legend.text = element_text(size = 18),
    panel.grid.major = element_line(color = "grey", linewidth = 0.1, linetype = "dashed"),
    panel.grid.minor = element_line(color = "grey", linewidth = 0.05, linetype = "dotted")
  ) +
  scale_x_date(date_labels = "%Y", date_breaks = "5 years")
#### 4 Test Set ----
##### 4.1 Rolling Window Forecasting -----
###### 4.1.1 GARCH ------
garch_forecast <- rep(NA, times = length(dax_full$r) - length(training_dax$r))

# Loop for out-of-sample forecasts
for (t in 1:(length(dax_full$r) - length(training_dax$r))) {
  
  
  start_index <- t
  end_index <- round(0.65 * length(dax_full$r)) - 1 + t
  
  if (end_index > length(dax_full$r)) {
    break  
  }
  
  
  current_data <- dax_full$r[start_index:end_index]
  
  
  spec_garch <- ugarchspec(
    variance.model = list(
      model = "sGARCH", 
      garchOrder = c(garchpq[row_min_garch, 1], garchpq[row_min_garch, 2])
    ),
    mean.model = list(armaOrder = c(0, 0), include.mean=FALSE),
    distribution.model = "norm"
  )
  
  
  fit_garch <- ugarchfit(spec = spec_garch, data = current_data, silent = TRUE)
  
  # Forecast
  forecast_garch <- ugarchforecast(fit_garch, n.ahead = 1)
  garch_forecast[t] <- forecast_garch@forecast$sigmaFor[1]^2
}


######  4.1.2 EGARCH ------

egarch_forecast <- rep(NA, times = length(dax_full$r) - length(training_dax$r))

# Loop for out-of-sample forecasts
for (t in 1:(length(dax_full$r) - length(training_dax$r))) {
  
  
  start_index <- t
  end_index <- round(0.65 * length(dax_full$r)) - 1 + t
  
  if (end_index > length(dax_full$r)) {
    break  
  }
  
  
  current_data <- dax_full$r[start_index:end_index]
  
  spec_garch <- ugarchspec(
    variance.model = list(
      model = "eGARCH", 
      garchOrder = c(egarchpq[row_min_egarch, 1], garchpq[row_min_egarch, 2]),
      variance.targeting=TRUE
    ),
    mean.model = list(armaOrder = c(0, 0), include.mean=FALSE),
    distribution.model = "norm"
  )
  
  #Fit
  fit_garch <- ugarchfit(spec = spec_garch, data = current_data, silent = TRUE, solver = "hybrid")
  
  # Forecast 
  forecast_garch <- ugarchforecast(fit_garch, n.ahead = 1)
  egarch_forecast[t] <- forecast_garch@forecast$sigmaFor[1]^2
}

###### 4.1.3 GARCHX ------
garchx_forecast <- rep(NA, times = length(dax_full$r) - length(training_dax$r))

# Loop for out-of-sample forecasts
for (t in 1:(length(dax_full$r) - length(training_dax$r))) {
  
  
  start_index <- t
  end_index <- round(0.65 * length(dax_full$r)) - 1 + t
  
  
  if (end_index > length(dax_full$r)) {
    break  
  }
  
  
  current_data <- dax_full$r[start_index:end_index]
  
  
  exog_vars <- dax_full$vix_t1[start_index:end_index]
  if (!is.matrix(exog_vars)) {
    exog_vars <- matrix(exog_vars, ncol = 1)  
  }
  
  if (nrow(exog_vars) != length(current_data)) {
    stop("Mismatch between external regressors and data length")
  }
  

  spec_garch <- ugarchspec(
    variance.model = list(
      model = "sGARCH", 
      garchOrder = c(5, 1),
      external.regressors = exog_vars
    ),
    mean.model = list(armaOrder = c(0, 0), include.mean = FALSE),
    distribution.model = "norm"
  )
  
  # Fit
  fit_garch <- ugarchfit(spec = spec_garch, data = current_data, silent = TRUE)
  
  # Forecast 
  forecast_garch <- ugarchforecast(fit_garch, n.ahead = 1)
  garchx_forecast[t] <- forecast_garch@forecast$sigmaFor[1]^2
}
##### 4.2 Display Out-of-Sample Forecasts (Figure 4) -----
colors <- c("cyan", "orange", "forestgreen")
data.frame(
  Date = dax_full$Date[1185:1821],
  RV = sqrt(dax_full$rv[1185:1821]),
  GARCH = sqrt(garch_forecast),
  EGARCH = sqrt(egarch_forecast),
  GARCHX = sqrt(garchx_forecast)
) %>%
  ggplot(aes(x = Date)) +
  geom_line(aes(y = RV, color = "Realized Volatility"), linewidth = 0.5, alpha = 0.6) +
  geom_line(aes(y = GARCH, color = "GARCH Model"), linewidth = 0.75) +
  geom_line(aes(y = EGARCH, color = "EGARCH Model"), linewidth = 0.75) +
  geom_line(aes(y = GARCHX, color = "GARCHX Model"), linewidth = 0.75) +
  scale_color_manual(
    values = c(
      "Realized Volatility" = "grey26",  
      "GARCH Model" = colors[1],         
      "EGARCH Model" = colors[2],        
      "GARCHX Model" = colors[3]         
    )
  ) +
  labs(
    x = "Time",
    y = "Volatility (Standard Deviation)",
    title = "Comparison of Out-of-Sample Volatility Forecasts",
    subtitle = "GARCH vs. GARCHX vs. EGARCH",
    color = "Model"
  ) +
  theme_minimal() +
  theme(
    plot.title = element_text(hjust = 0.5, size = 20, face = "bold"),
    plot.subtitle = element_text(hjust = 0.5, size = 18, face = "bold"),
    axis.text = element_text(size = 12),
    axis.title = element_text(size = 18),
    legend.position = "right",
    legend.title = element_text(face = "bold"),
    legend.text = element_text(size = 18),
    panel.grid.major = element_line(color = "grey", linewidth = 0.1, linetype = "dashed"),
    panel.grid.minor = element_line(color = "grey", linewidth = 0.05, linetype = "dotted")
  ) +
  scale_x_date(date_labels = "%Y", date_breaks = "3 years")
##### 4.3 Forecast evaluation -----
###### 4.3.1 Evaluation Metrics ----- 
dax_evaluation<- data.frame(
  Date = dax_full$Date[1185:1821],
  RV = dax_full$rv[1185:1821],
  GARCH = garch_forecast,
  EGARCH = egarch_forecast,
  GARCHX = garchx_forecast
)


# Variance
mae(dax_evaluation$RV, dax_evaluation$GARCH)
mae(dax_evaluation$RV, dax_evaluation$EGARCH)
mae(dax_evaluation$RV, dax_evaluation$GARCHX)
# Volatility
mae(dax_evaluation$RV %>% sqrt(), dax_evaluation$GARCH %>% sqrt())
mae(dax_evaluation$RV %>% sqrt(), dax_evaluation$GARCHX %>% sqrt())
mae(dax_evaluation$RV %>% sqrt() , dax_evaluation$EGARCH %>% sqrt())


###### 4.3.2 Mincer-Zarnowitz-Regression ------
####### 4.3.2.1 GARCH -------
mz_garch<-lm(RV~GARCH, data=dax_evaluation)
coeftest(mz_garch, vcov = vcovHAC(mz_garch))
mz_garch %>% summary()
linearHypothesis(mz_garch, 
                 hypothesis.matrix = c("(Intercept) = 0", "GARCH = 1"), 
                 vcov = vcovHAC(mz_garch))

####### 4.3.2.2 EGARCH -------
mz_egarch<-lm(RV~EGARCH, data=dax_evaluation)
mz_egarch %>% summary()
coeftest(mz_egarch, vcov = vcovHAC(mz_egarch))
linearHypothesis(mz_egarch, 
                 hypothesis.matrix = c("(Intercept) = 0", "EGARCH= 1"), 
                 vcov = vcovHAC(mz_egarch))

####### 4.3.2.3 GARCHX -------
mz_garchx<-lm(RV~GARCHX, data=dax_evaluation)
mz_garchx %>% summary()
coeftest(mz_garchx, vcov = vcovHAC(mz_garchx))
linearHypothesis(mz_garchx, 
                 hypothesis.matrix = c("(Intercept) = 0", "GARCHX = 1"), 
                 vcov = vcovHAC(mz_garchx))


mz_garch_acf <- ggAcf(mz_garch$residuals, lag.max = 26, plot = FALSE) %>%
  autoplot() +
  labs(title = "ACF of MZ GARCH Residuals for DAX", x = "Lags", y = "ACF") +
  theme_minimal(base_size = 14) +
  theme(plot.title = element_text(face = "bold"))


####### 4.3.2.4 Display MZ Residual Correlogram (Figure 13) -------
mz_garch_pacf <- ggPacf(mz_garch$residuals, lag.max = 26, plot = FALSE) %>%
  autoplot() +
  labs(title = "PACF of MZ GARCH Residuals for DAX", x = "Lags", y = "PACF") +
  theme_minimal(base_size = 14) +
  theme(plot.title = element_text(face = "bold")) 

grid.arrange(mz_garch_acf,mz_garch_pacf, ncol=2)

###### 4.3.3 Diebold Mariano Test ------
####### 4.3.3.1 Loss Function -------
dax_evaluation<- dax_evaluation %>% 
  mutate(ae_garch=abs(RV-GARCH)) %>% 
  mutate(ae_egarch=abs(RV-EGARCH)) %>% 
  mutate(ae_garchx=abs(RV-GARCHX))

####### 4.3.3.2 GARCH,EGARCH -------
dm_garch_egarch<-dm.test(e1= dax_evaluation$ae_garch,
                         e2= dax_evaluation$ae_egarch,
                         h=1,
                         alternative = "two.sided",
                         power = 1)
dm_garch_egarch_g<-dm.test(e1= dax_evaluation$ae_garch,
                         e2= dax_evaluation$ae_egarch,
                         h=1,
                         alternative = "greater",
                         power = 1)
dm_garch_egarch_l<-dm.test(e1= dax_evaluation$ae_garch,
                           e2= dax_evaluation$ae_egarch,
                           h=1,
                           alternative = "less",
                           power = 1)

####### 4.3.3.3 GARCH,GARCHX -------
dm_garch_garchx<-dm.test(e1= dax_evaluation$ae_garch,
                         e2= dax_evaluation$ae_garchx,
                         h=1,
                         alternative = "two.sided",
                         power = 1)
dm_garch_garchx_g<-dm.test(e1= dax_evaluation$ae_garch,
                         e2= dax_evaluation$ae_garchx,
                         h=1,
                         alternative = "greater",
                         power = 1)
dm_garch_garchx_l<-dm.test(e1= dax_evaluation$ae_garch,
                           e2= dax_evaluation$ae_garchx,
                           h=1,
                           alternative = "less",
                           power = 1)
####### 4.3.3.4 EGARCH,GARCHX -------
dm_egarch_garchx<-dm.test(e1= dax_evaluation$ae_egarch,
                          e2= dax_evaluation$ae_garchx,
                          h=1,
                          alternative = "two.sided",
                          power = 1)

