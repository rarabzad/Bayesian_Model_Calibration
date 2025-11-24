## 1. Priors — expressing what you believe before data

### Purpose

The **prior** captures what you know (or assume) about parameters before seeing data. It acts as a regularizer, preventing unrealistic parameter values and stabilizing estimation.

### Basic formula

$$p(\theta) = \text{prior distribution of parameters}$$

where $\theta$ represents model parameters (e.g., hydraulic conductivity, reaction rate, etc.).

### Common types of priors

| Type                       | Use case                                              | Example                              |
| -------------------------- | ----------------------------------------------------- | ------------------------------------ |
| **Non-informative (flat)** | when you know almost nothing                          | Uniform(a, b)                        |
| **Weakly informative**     | when you want to constrain unreasonable values        | Normal(0, 10) for a regression slope |
| **Informative**            | when expert knowledge or previous data exist          | LogNormal(mean, sd)                  |
| **Hierarchical**           | when parameters vary across groups (e.g., catchments) | θᵢ ~ Normal(μ, σ)                    |

### Practical rules

* Use weakly informative priors as a default.
* Check the prior predictive distribution — simulate from priors before using data.
* Avoid overly tight priors unless they reflect strong prior knowledge.

---

## 2. Error Model — describing uncertainty in data

### Why it matters

The **error model** defines how observed data deviate from the model output. This is often where misspecification hides.

### Simple additive error model

$$y_i = f(x_i, \theta) + \varepsilon_i$$

where $\varepsilon_i \sim \mathcal{N}(0, \sigma^2)$

### Multiplicative error model (log-transformed)

$$y_i = f(x_i, \theta) \times (1 + \varepsilon_i), \quad \varepsilon_i \sim \mathcal{N}(0, \sigma^2)$$

Implies:

$$\log(y_i) = \log(f(x_i, \theta)) + \eta_i, \quad \eta_i \sim \mathcal{N}(0, \sigma^2)$$

— useful for skewed or strictly positive data.

### Student-t error model (heavy-tailed)

$$\varepsilon_i \sim t_\nu(0, \sigma^2)$$

— robust against outliers; $\nu$ controls tail heaviness.

### Heteroscedastic error model

$$\varepsilon_i \sim \mathcal{N}(0, (\sigma_0 + \sigma_1 |f(x_i, \theta)|)^2)$$

— allows larger variance for higher predictions.

### Autocorrelated errors (time series)

$$\varepsilon_t = \phi \varepsilon_{t-1} + \eta_t, \quad \eta_t \sim \mathcal{N}(0, \sigma^2)$$

— captures persistent residuals.

---

## 3. Likelihood — connecting data and model

### General definition

$$p(y \mid \theta) = \prod_{i=1}^{n} p(y_i \mid \theta)$$

$$\log p(y \mid \theta) = \sum_{i=1}^{n} \log p(y_i \mid \theta)$$

### Normal likelihood (additive Gaussian errors)

$$y_i \sim \mathcal{N}(f(x_i, \theta), \sigma^2)$$

$$\log p(y \mid \theta) = -\frac{n}{2}\log(2\pi\sigma^2) - \frac{1}{2\sigma^2}\sum_i (y_i - f(x_i, \theta))^2$$

### Log-Normal likelihood (multiplicative errors)

$$\log(y_i) \sim \mathcal{N}(\log(f(x_i, \theta)), \sigma^2)$$

$$\log p(y \mid \theta) = -\frac{n}{2}\log(2\pi\sigma^2) - \sum_i \left( \frac{(\log(y_i) - \log(f_i))^2}{2\sigma^2} + \log(y_i) \right)$$

### Student-t likelihood (robust)

$$p(y_i \mid \theta) = \frac{\Gamma((\nu+1)/2)}{\Gamma(\nu/2)\sqrt{\pi\nu}\sigma}\left(1 + \frac{(y_i - f_i)^2}{\nu\sigma^2}\right)^{-(\nu+1)/2}$$

$$\log p(y \mid \theta) = \sum_i \log p(y_i \mid \theta)$$

### Autocorrelated likelihood (AR(1))

$$\varepsilon_t = y_t - f_t$$

$$\varepsilon_t = \phi \varepsilon_{t-1} + \eta_t, \quad \eta_t \sim \mathcal{N}(0, \sigma^2)$$

$$\log p(y \mid \theta) = -\frac{1}{2}\sum_t \left( \frac{(\varepsilon_t - \phi\varepsilon_{t-1})^2}{\sigma^2} + \log(2\pi\sigma^2) \right)$$

### Heteroscedastic likelihood

$$y_i \sim \mathcal{N}(f_i, (\sigma_0 + \sigma_1|f_i|)^2)$$

$$\log p(y \mid \theta) = -\frac{1}{2}\sum_i \left( \frac{(y_i - f_i)^2}{(\sigma_0 + \sigma_1|f_i|)^2} + \log(2\pi(\sigma_0 + \sigma_1|f_i|)^2) \right)$$

### Box–Cox transformation

$$y_i^{(\lambda)} = \begin{cases} \frac{y_i^\lambda - 1}{\lambda}, & \lambda \neq 0 \\ \log(y_i), & \lambda = 0 \end{cases}$$

$$y_i^{(\lambda)} \sim \mathcal{N}(f_i^{(\lambda)}, \sigma^2)$$

---

## 4. Posterior — combining prior and likelihood

### Definition

$$p(\theta \mid y) = \frac{p(y \mid \theta) p(\theta)}{p(y)}$$

* Prior: $p(\theta)$
* Likelihood: $p(y \mid \theta)$
* Posterior: $p(\theta \mid y)$
* Evidence (normalization): $p(y) = \int p(y \mid \theta)p(\theta) d\theta$

### Summaries

* Mean / median of $\theta$
* 95% credible intervals
* Posterior predictive checks: $y^{rep} \sim p(y^{rep} \mid \theta^{(s)})$

### Inference methods

| Method                | Notes                                  |
| --------------------- | -------------------------------------- |
| MCMC (HMC/NUTS)       | Standard, accurate, draws samples      |
| Variational Inference | Fast approximation                     |
| MAP                   | Mode of posterior, no uncertainty info |

---

## 5. Conceptual example

$$y_i \sim \mathcal{N}(f(x_i, \theta), \sigma^2)$$

$$\theta \sim \mathcal{N}(0,1), \quad \sigma \sim \text{HalfNormal}(1)$$

Posterior:

$$p(\theta,\sigma \mid y) \propto \prod_i \mathcal{N}(y_i \mid f(x_i,\theta),\sigma^2) \times \mathcal{N}(\theta\mid 0,1) \times \text{HalfNormal}(\sigma\mid 1)$$

Sample with MCMC → posterior predictive checks.
