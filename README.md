## 📊 Overview  

### 🎯 Objective  
Evaluate and compare the predictive performance of different **GARCH-type models** in forecasting **weekly volatility**, using **out-of-sample evaluation techniques**.  

### 🔍 Key Methodologies  
- **📈 ARCH/GARCH Modeling** – Testing for ARCH effects and fitting **GARCH, GARCHX (with exogenous variables), and EGARCH** models.  
- **⚖ Model Selection** – Using the **Akaike Information Criterion (AIC)** to identify the best in-sample fit.  
- **🔄 Forecasting Approach** – Implementing a **rolling window framework** to generate dynamic volatility predictions.  
- **📊 Forecast Evaluation:**  
  - 📌 **Realized Variance Approach** – Weekly volatility forecasts were evaluated using **realized variance**, computed from **daily returns**.  
  - 📉 **Mean Absolute Error (MAE)** for both **variance and volatility forecasts**.  
  - 📊 **Mincer-Zarnowitz regression** to assess **forecast optimality**.  
  - 🔍 **Diebold-Mariano tests** to compare **predictive accuracy** between models.  
  
### 📂 Data  
📡 **Financial time series data** (DAX, Bitcoin, EUR/RUB) sourced from **Yahoo Finance**, with additional market indicators such as the **Volatility Index (VIX)** used as an **exogenous regressor**.  
