## Goals

* Give a **clear intuition** for Bayesian thinking.
* Teach a **practical calibration workflow** (prior → likelihood → posterior → checks).
* Provide short exercises and a checklist to run real calibrations.

## Audience

Practitioners who want to calibrate models (hydrology, ecology, engineering, etc.) with uncertainty — no deep Bayesian background required.

---

# 1. Quick Bayesian intuition

* **Bayes = update beliefs with data.**
* Start with what you *think* (prior), see data (likelihood), get what you *believe now* (posterior).

### Bayes' rule (compact)

$$p(\theta\mid y) = \frac{p(y\mid \theta)\,p(\theta)}{p(y)}$$

* $p(\theta)$: prior (what you believed before seeing data).
* $p(y\mid\theta)$: likelihood (how likely the data are given parameters).
* $p(\theta\mid y)$: posterior (updated belief).
* $p(y)$: normalizing constant (marginal likelihood).

---

# 2. Key terms — one-line cheatsheet

* **Prior**: beliefs about parameters before data.
* **Likelihood**: model of the data given parameters.
* **Posterior**: updated distribution of parameters.
* **Predictive (posterior predictive)**: distribution of new data given the posterior.
* **MCMC / HMC / NUTS**: algorithms to sample from the posterior.
* **Credible interval**: Bayesian interval for parameter values (e.g., 95% credible interval).
* **PP-checks**: simulate replicated data from posterior; compare to observed.

---

# 3. Bayesian calibration — the practical workflow

1. **Define model (process + noise).** What generates data? Which parameters matter?
2. **Choose priors.** Document choices and reasons.
3. **Specify likelihood.** E.g., Normal errors, log-Normal, Poisson — match data type and error structure.
4. **Fit (compute posterior).** Use sampling (MCMC) or approximation (VI).
5. **Check diagnostics.** Convergence (R̂, ESS), trace plots, divergences.
6. **Posterior predictive checks.** Does the model reproduce key data features?
7. **Summarize & report.** Posterior summaries, credible intervals, predictive uncertainty.
8. **(Optional) Model comparison & robustness.** LOO, WAIC, sensitivity to priors.

---

# 4. Intuitive guide to each component

## 4.1 Priors — minimal words, maximum clarity

* **Purpose**: regularize estimation, incorporate prior knowledge, avoid nonsense.
* **Types**: informative (strong knowledge), weakly informative (mild constraints), non-informative (flat).
* **Rule of thumb**: prefer weakly informative priors that rule out impossible values but don't force a result.
* **Document** your prior choices and test sensitivity by changing them.

## 4.2 Likelihood — what matters

* Choose a likelihood that reflects measurement error and process noise.
* Consider transformations (log for multiplicative errors).
* If residuals show heteroscedasticity or heavy tails, use appropriate error models (e.g., Student-t).

## 4.3 Posterior — what you inspect

* Look at marginal distributions and joint correlations.
* Check credible intervals for practical significance (not just statistical).
* Use posterior predictive distributions to inspect model fit in data space.

---

# 5. Inference methods (brief)

* **MCMC (standard)** — accurate but can be slow. HMC / NUTS are efficient for many parameters.
* **Variational inference (VI)** — approximates posterior fast; can miss tails. Good for large-scale or initial exploration.
* **Optimization / MAP** — single best point estimate; loses uncertainty.

---

# 6. Diagnostics & checks (must-do list)

* **Convergence**: R̂ ≈ 1, no trends in trace plots.
* **Mixing**: effective sample size (ESS) sufficiently large.
* **Divergences** (HMC): fix by reparametrizing or stronger priors.
* **Posterior predictive checks**: overlay simulated vs observed; check summary stats (means, quantiles, autocorrelation).
* **Residuals**: visually inspect residuals for bias, heteroscedasticity, autocorrelation.
* **Sensitivity**: re-run with different priors; check results stable.

---

# 7. Short, practical examples (copy-paste friendly)

### Minimal Stan/brms-style pseudocode (concept)

```r
library(brms)

# Define model
model <- brm(
  y ~ 1,                                    # formula
  prior = c(
    prior(normal(0, 1), class = Intercept), # prior on theta
    prior(normal(0, 1), class = sigma)      # prior on sigma
  ),
  data = data.frame(y = y, x = x),
  chains = 4, 
  iter = 2000, 
  warmup = 1000
)
```

### Posterior predictive (concept)

```r
# Posterior predictive checks
pp_check(model, ndraws = 100)

# Or manually
ppc <- posterior_predict(model)
plot(density(ppc), main = "Posterior Predictive")
points(y, rep(0, length(y)), col = "red")
```

> Adapt to your specific model structure and Stan/brms/rstan workflow.

---

# 8. Reporting & reproducibility — checklist

* Save code, seeds, and data.
* Record priors, likelihood choice, sampler settings, and versions.
* Include diagnostic plots and posterior predictive checks.
* Run sensitivity to priors and alternative noise models.

---

# 9. Common pitfalls (quick)

* Overconfident priors that dominate data.
* Ignoring model misspecification (bad likelihood).
* Accepting convergence diagnostics without visual checks.
* Reporting point estimates without uncertainty.

---

# 10. Priors — expressing what you believe before data

## 10.1 Purpose

The **prior** captures what you know (or assume) about parameters before seeing data. It acts as a regularizer, preventing unrealistic parameter values and stabilizing estimation.

## 10.2 Basic formula

$$p(\theta) = \text{prior distribution of parameters}$$

where $\theta$ represents model parameters (e.g., hydraulic conductivity, reaction rate, etc.).

## 10.3 Common types of priors

| Type                       | Use case                                              | Example                              |
| -------------------------- | ----------------------------------------------------- | ------------------------------------ |
| **Non-informative (flat)** | when you know almost nothing                          | Uniform(a, b)                        |
| **Weakly informative**     | when you want to constrain unreasonable values        | Normal(0, 10) for a regression slope |
| **Informative**            | when expert knowledge or previous data exist          | LogNormal(mean, sd)                  |
| **Hierarchical**           | when parameters vary across groups (e.g., catchments) | θᵢ ~ Normal(μ, σ)                    |

## 10.4 Practical rules

* Use weakly informative priors as a default.
* Check the prior predictive distribution — simulate from priors before using data.
* Avoid overly tight priors unless they reflect strong prior knowledge.

---

# 11. Error Model — describing uncertainty in data

## 11.1 Why it matters

The **error model** defines how observed data deviate from the model output. This is often where misspecification hides.

## 11.2 Simple additive error model

$$y_i = f(x_i, \theta) + \varepsilon_i$$

where $\varepsilon_i \sim \mathcal{N}(0, \sigma^2)$

## 11.3 Multiplicative error model (log-transformed)

$$y_i = f(x_i, \theta) \times (1 + \varepsilon_i), \quad \varepsilon_i \sim \mathcal{N}(0, \sigma^2)$$

Implies:

$$\log(y_i) = \log(f(x_i, \theta)) + \eta_i, \quad \eta_i \sim \mathcal{N}(0, \sigma^2)$$

— useful for skewed or strictly positive data.

## 11.4 Student-t error model (heavy-tailed)

$$\varepsilon_i \sim t_\nu(0, \sigma^2)$$

— robust against outliers; $\nu$ controls tail heaviness.

## 11.5 Heteroscedastic error model

$$\varepsilon_i \sim \mathcal{N}(0, (\sigma_0 + \sigma_1 |f(x_i, \theta)|)^2)$$

— allows larger variance for higher predictions.

## 11.6 Autocorrelated errors (time series)

$$\varepsilon_t = \phi \varepsilon_{t-1} + \eta_t, \quad \eta_t \sim \mathcal{N}(0, \sigma^2)$$

— captures persistent residuals.

---

# 12. Likelihood — connecting data and model

## 12.1 General definition

$$p(y \mid \theta) = \prod_{i=1}^{n} p(y_i \mid \theta)$$

$$\log p(y \mid \theta) = \sum_{i=1}^{n} \log p(y_i \mid \theta)$$

## 12.2 Normal likelihood (additive Gaussian errors)

$$y_i \sim \mathcal{N}(f(x_i, \theta), \sigma^2)$$

$$\log p(y \mid \theta) = -\frac{n}{2}\log(2\pi\sigma^2) - \frac{1}{2\sigma^2}\sum_i (y_i - f(x_i, \theta))^2$$

## 12.3 Log-Normal likelihood (multiplicative errors)

$$\log(y_i) \sim \mathcal{N}(\log(f(x_i, \theta)), \sigma^2)$$

$$\log p(y \mid \theta) = -\frac{n}{2}\log(2\pi\sigma^2) - \sum_i \left( \frac{(\log(y_i) - \log(f_i))^2}{2\sigma^2} + \log(y_i) \right)$$

## 12.4 Student-t likelihood (robust)

$$p(y_i \mid \theta) = \frac{\Gamma((\nu+1)/2)}{\Gamma(\nu/2)\sqrt{\pi\nu}\sigma}\left(1 + \frac{(y_i - f_i)^2}{\nu\sigma^2}\right)^{-(\nu+1)/2}$$

$$\log p(y \mid \theta) = \sum_i \log p(y_i \mid \theta)$$

## 12.5 Autocorrelated likelihood (AR(1))

$$\varepsilon_t = y_t - f_t$$

$$\varepsilon_t = \phi \varepsilon_{t-1} + \eta_t, \quad \eta_t \sim \mathcal{N}(0, \sigma^2)$$

$$\log p(y \mid \theta) = -\frac{1}{2}\sum_t \left( \frac{(\varepsilon_t - \phi\varepsilon_{t-1})^2}{\sigma^2} + \log(2\pi\sigma^2) \right)$$

## 12.6 Heteroscedastic likelihood

$$y_i \sim \mathcal{N}(f_i, (\sigma_0 + \sigma_1|f_i|)^2)$$

$$\log p(y \mid \theta) = -\frac{1}{2}\sum_i \left( \frac{(y_i - f_i)^2}{(\sigma_0 + \sigma_1|f_i|)^2} + \log(2\pi(\sigma_0 + \sigma_1|f_i|)^2) \right)$$

## 12.7 Box–Cox transformation

$$y_i^{(\lambda)} = \begin{cases} \frac{y_i^\lambda - 1}{\lambda}, & \lambda \neq 0 \\ \log(y_i), & \lambda = 0 \end{cases}$$

$$y_i^{(\lambda)} \sim \mathcal{N}(f_i^{(\lambda)}, \sigma^2)$$

---

# 13. Posterior — combining prior and likelihood

## 13.1 Definition

$$p(\theta \mid y) = \frac{p(y \mid \theta) p(\theta)}{p(y)}$$

* Prior: $p(\theta)$
* Likelihood: $p(y \mid \theta)$
* Posterior: $p(\theta \mid y)$
* Evidence (normalization): $p(y) = \int p(y \mid \theta)p(\theta) d\theta$

## 13.2 Summaries

* Mean / median of $\theta$
* 95% credible intervals
* Posterior predictive checks: $y^{rep} \sim p(y^{rep} \mid \theta^{(s)})$

## 13.3 Inference methods

| Method                | Notes                                  |
| --------------------- | -------------------------------------- |
| MCMC (HMC/NUTS)       | Standard, accurate, draws samples      |
| Variational Inference | Fast approximation                     |
| MAP                   | Mode of posterior, no uncertainty info |

---

# 14. Conceptual example

$$y_i \sim \mathcal{N}(f(x_i, \theta), \sigma^2)$$

$$\theta \sim \mathcal{N}(0,1), \quad \sigma \sim \text{HalfNormal}(1)$$

Posterior:

$$p(\theta,\sigma \mid y) \propto \prod_i \mathcal{N}(y_i \mid f(x_i,\theta),\sigma^2) \times \mathcal{N}(\theta\mid 0,1) \times \text{HalfNormal}(\sigma\mid 1)$$

Sample with MCMC → posterior predictive checks.


# Bayesian Linear Regression Example 

## 1. Problem Definition

We want to calibrate a simple linear regression model:

$$y = \beta_0 + \beta_1 x + \varepsilon$$

where:
- $\beta_0$ is the intercept
- $\beta_1$ is the slope
- $\varepsilon \sim \mathcal{N}(0, \sigma^2)$ is Gaussian noise

**Goal**: Estimate the parameters $\theta = (\beta_0, \beta_1, \sigma)$ from observed data using Bayesian inference.

---

## 2. Synthetic Data Generation

We'll create synthetic data using known "true" parameters:

**True parameters:**
- $\beta_0^{true} = 2.5$ (intercept)
- $\beta_1^{true} = 1.8$ (slope)
- $\sigma^{true} = 1.2$ (noise standard deviation)

**Data generation process:**
1. Generate $n = 50$ values of $x$ uniformly between 0 and 10
2. Compute $y = 2.5 + 1.8x + \varepsilon$, where $\varepsilon \sim \mathcal{N}(0, 1.2^2)$

---

## 3. Bayesian Model Specification

### 3.1 Priors

We choose weakly informative priors that allow a wide range of reasonable values:

$$\beta_0 \sim \mathcal{N}(0, 10^2)$$
$$\beta_1 \sim \mathcal{N}(0, 10^2)$$
$$\sigma \sim \text{HalfNormal}(5)$$

The HalfNormal prior for $\sigma$ ensures it stays positive (standard deviations must be > 0).

### 3.2 Likelihood

With Gaussian homoscedastic (constant variance) and non-autocorrelated errors:

$$y_i \sim \mathcal{N}(\beta_0 + \beta_1 x_i, \sigma^2)$$

The log-likelihood is:

$$\log p(y \mid \beta_0, \beta_1, \sigma) = -\frac{n}{2}\log(2\pi\sigma^2) - \frac{1}{2\sigma^2}\sum_{i=1}^{n}(y_i - \beta_0 - \beta_1 x_i)^2$$

### 3.3 Posterior

By Bayes' theorem:

$$p(\beta_0, \beta_1, \sigma \mid y) \propto p(y \mid \beta_0, \beta_1, \sigma) \times p(\beta_0) \times p(\beta_1) \times p(\sigma)$$

The posterior combines:
- The likelihood (how well parameters explain the data)
- The priors (our initial beliefs about parameters)

---

## 4. R Code Implementation with adaptMCMC

### Step 1: Install and load required packages

```r
# Install packages if needed
if (!require("adaptMCMC")) install.packages("adaptMCMC")
if (!require("coda")) install.packages("coda")

# Load libraries
library(adaptMCMC)
library(coda)
```

### Step 2: Generate synthetic data

```r
set.seed(123)  # For reproducibility

# True parameters
beta0_true <- 2.5
beta1_true <- 1.8
sigma_true <- 1.2

# Generate data
n <- 50
x <- runif(n, min = 0, max = 10)
y <- beta0_true + beta1_true * x + rnorm(n, mean = 0, sd = sigma_true)

# Visualize the data
plot(x, y, pch = 19, col = "steelblue", 
     xlab = "x", ylab = "y", 
     main = "Synthetic Data for Linear Regression")
abline(a = beta0_true, b = beta1_true, col = "red", lwd = 2, lty = 2)
legend("topleft", legend = "True relationship", 
       col = "red", lty = 2, lwd = 2)
```

### Step 3: Define the log-posterior function
$$p(\theta\mid y) = \frac{p(y\mid \theta)\,p(\theta)}{p(y)}$$
<img src="https://m-clark.github.io/bayesian-basics/Bayesian-Basics_files/figure-html/prior2post_1-1.svg" alt="Bayesian Workflow">

*Figure: Bayesian inference workflow - combining prior beliefs with data to obtain the posterior distribution.*
```r
# Log-prior function
log_prior <- function(theta) {
  beta0 <- theta[1]
  beta1 <- theta[2]
  sigma <- theta[3]
  
  # Priors: beta0 ~ N(0, 10^2), beta1 ~ N(0, 10^2), sigma ~ HalfNormal(5)
  if (sigma <= 0) return(-Inf)  # Ensure sigma > 0
  
  lp_beta0 <- dnorm(beta0, mean = 0, sd = 10, log = TRUE)
  lp_beta1 <- dnorm(beta1, mean = 0, sd = 10, log = TRUE)
  # HalfNormal(5) is equivalent to |N(0, 5^2)|
  lp_sigma <- dnorm(sigma, mean = 0, sd = 5, log = TRUE) + log(2)
  
  return(lp_beta0 + lp_beta1 + lp_sigma)
}

# Log-likelihood function
log_likelihood <- function(theta, x, y) {
  beta0 <- theta[1]
  beta1 <- theta[2]
  sigma <- theta[3]
  
  # Predicted values
  y_pred <- beta0 + beta1 * x
  
  # Log-likelihood: sum of log-densities
  ll <- sum(dnorm(y, mean = y_pred, sd = sigma, log = TRUE))
  
  return(ll)
}

# Log-posterior = log-prior + log-likelihood
log_posterior <- function(theta, x, y) {
  lp <- log_prior(theta)
  if (is.infinite(lp)) return(lp)  # Don't compute likelihood if prior is -Inf
  
  ll <- log_likelihood(theta, x, y)
  return(lp + ll)
}
```

### Step 4: Run MCMC sampling with adaptMCMC

```r
# Initial values for parameters (starting point for MCMC)
theta_init <- c(0, 0, 1)  # (beta0, beta1, sigma)

# Run adaptive MCMC
# adaptMCMC automatically tunes the proposal distribution
mcmc_result <- MCMC(
  p = log_posterior,           # log-posterior function
  init = theta_init,           # initial values
  n = 10000,                   # number of iterations
  adapt = TRUE,                # enable adaptive sampling
  acc.rate = 0.234,            # target acceptance rate
  scale = c(0.5, 0.5, 0.2),   # initial proposal scales
  x = x,                       # pass data to log_posterior
  y = y
)

# Extract samples (discard burn-in)
burn_in <- 2000
samples <- mcmc_result$samples[(burn_in + 1):nrow(mcmc_result$samples), ]
colnames(samples) <- c("beta0", "beta1", "sigma")

# Convert to coda mcmc object for diagnostics
samples_mcmc <- as.mcmc(samples)
```

### Step 5: Check convergence diagnostics

```r
# Trace plots
par(mfrow = c(3, 1))
plot(samples_mcmc[, "beta0"], main = "Trace plot: beta0", ylab = "beta0")
abline(h = beta0_true, col = "red", lwd = 2, lty = 2)

plot(samples_mcmc[, "beta1"], main = "Trace plot: beta1", ylab = "beta1")
abline(h = beta1_true, col = "red", lwd = 2, lty = 2)

plot(samples_mcmc[, "sigma"], main = "Trace plot: sigma", ylab = "sigma")
abline(h = sigma_true, col = "red", lwd = 2, lty = 2)
par(mfrow = c(1, 1))

# Acceptance rate
cat("Acceptance rate:", mcmc_result$acceptance.rate, "\n")

# Effective sample size
effectiveSize(samples_mcmc)

# Autocorrelation
par(mfrow = c(3, 1))
autocorr.plot(samples_mcmc[, "beta0"], main = "Autocorrelation: beta0")
autocorr.plot(samples_mcmc[, "beta1"], main = "Autocorrelation: beta1")
autocorr.plot(samples_mcmc[, "sigma"], main = "Autocorrelation: sigma")
par(mfrow = c(1, 1))
```

### Step 6: Summarize posterior distributions

```r
# Posterior summaries
summary(samples_mcmc)

# 95% credible intervals
apply(samples, 2, quantile, probs = c(0.025, 0.5, 0.975))

# Posterior means
posterior_means <- colMeans(samples)
cat("\nPosterior means:\n")
print(posterior_means)

cat("\nTrue values:\n")
print(c(beta0 = beta0_true, beta1 = beta1_true, sigma = sigma_true))
```

### Step 7: Posterior density plots

```r
par(mfrow = c(2, 2))

# Beta0
hist(samples[, "beta0"], breaks = 50, freq = FALSE, 
     col = "lightblue", border = "white",
     main = "Posterior: beta0", xlab = "beta0")
abline(v = beta0_true, col = "red", lwd = 2, lty = 2)
abline(v = posterior_means["beta0"], col = "blue", lwd = 2)
legend("topright", legend = c("True", "Posterior mean"), 
       col = c("red", "blue"), lty = c(2, 1), lwd = 2)

# Beta1
hist(samples[, "beta1"], breaks = 50, freq = FALSE, 
     col = "lightblue", border = "white",
     main = "Posterior: beta1", xlab = "beta1")
abline(v = beta1_true, col = "red", lwd = 2, lty = 2)
abline(v = posterior_means["beta1"], col = "blue", lwd = 2)
legend("topright", legend = c("True", "Posterior mean"), 
       col = c("red", "blue"), lty = c(2, 1), lwd = 2)

# Sigma
hist(samples[, "sigma"], breaks = 50, freq = FALSE, 
     col = "lightblue", border = "white",
     main = "Posterior: sigma", xlab = "sigma")
abline(v = sigma_true, col = "red", lwd = 2, lty = 2)
abline(v = posterior_means["sigma"], col = "blue", lwd = 2)
legend("topright", legend = c("True", "Posterior mean"), 
       col = c("red", "blue"), lty = c(2, 1), lwd = 2)

# Joint posterior (beta0 vs beta1)
plot(samples[, "beta0"], samples[, "beta1"], 
     pch = ".", col = rgb(0, 0, 1, 0.1),
     xlab = "beta0", ylab = "beta1",
     main = "Joint Posterior: beta0 vs beta1")
points(beta0_true, beta1_true, col = "red", pch = 19, cex = 2)
points(posterior_means["beta0"], posterior_means["beta1"], 
       col = "blue", pch = 19, cex = 2)
legend("topright", legend = c("True", "Posterior mean"), 
       col = c("red", "blue"), pch = 19)

par(mfrow = c(1, 1))
```

### Step 8: Posterior predictive checks

```r
# Generate posterior predictive samples
n_pred <- 100  # Number of posterior samples to use
pred_indices <- sample(1:nrow(samples), n_pred)

# Create prediction grid
x_pred <- seq(0, 10, length.out = 100)

# Storage for predictions
y_pred_samples <- matrix(NA, nrow = n_pred, ncol = length(x_pred))

for (i in 1:n_pred) {
  idx <- pred_indices[i]
  beta0_s <- samples[idx, "beta0"]
  beta1_s <- samples[idx, "beta1"]
  sigma_s <- samples[idx, "sigma"]
  
  # Generate predictions with uncertainty
  y_mean <- beta0_s + beta1_s * x_pred
  y_pred_samples[i, ] <- rnorm(length(x_pred), mean = y_mean, sd = sigma_s)
}

# Plot posterior predictive distribution
plot(x, y, pch = 19, col = "steelblue",
     xlab = "x", ylab = "y",
     main = "Posterior Predictive Check")

# Add posterior predictive samples (light lines)
for (i in 1:n_pred) {
  lines(x_pred, y_pred_samples[i, ], col = rgb(0, 0, 0, 0.05))
}

# Add posterior mean prediction
y_mean_pred <- posterior_means["beta0"] + posterior_means["beta1"] * x_pred
lines(x_pred, y_mean_pred, col = "blue", lwd = 2)

# Add true relationship
abline(a = beta0_true, b = beta1_true, col = "red", lwd = 2, lty = 2)

# Add data points again (on top)
points(x, y, pch = 19, col = "steelblue")

legend("topleft", 
       legend = c("Data", "True relationship", "Posterior mean", "Posterior samples"),
       col = c("steelblue", "red", "blue", rgb(0, 0, 0, 0.3)),
       pch = c(19, NA, NA, NA),
       lty = c(NA, 2, 1, 1),
       lwd = c(NA, 2, 2, 1))
```

### Step 9: Residual analysis

```r
# Calculate residuals using posterior mean parameters
y_fitted <- posterior_means["beta0"] + posterior_means["beta1"] * x
residuals <- y - y_fitted

# Residual plots
par(mfrow = c(2, 2))

# Residuals vs fitted
plot(y_fitted, residuals, pch = 19, col = "steelblue",
     xlab = "Fitted values", ylab = "Residuals",
     main = "Residuals vs Fitted")
abline(h = 0, col = "red", lty = 2)

# Residuals vs x
plot(x, residuals, pch = 19, col = "steelblue",
     xlab = "x", ylab = "Residuals",
     main = "Residuals vs x")
abline(h = 0, col = "red", lty = 2)

# QQ plot
qqnorm(residuals, pch = 19, col = "steelblue")
qqline(residuals, col = "red", lty = 2)

# Histogram of residuals
hist(residuals, breaks = 20, col = "lightblue", border = "white",
     main = "Histogram of Residuals", xlab = "Residuals")

par(mfrow = c(1, 1))
```

---

## 5. Interpretation

### What we learned:

1. **Convergence**: The trace plots show good mixing with no trends
2. **Parameter recovery**: Posterior means are close to true values
3. **Uncertainty quantification**: 95% credible intervals capture the true parameters
4. **Prediction**: Posterior predictive checks show the model captures the data well
5. **Residuals**: Residuals appear random with no systematic patterns
