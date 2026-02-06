# Bayesian Calibration for Hydrologists

## Theory, Components, and Practical Interpretation

---

## Table of Contents

* [1. Introduction](#1-introduction)
  * [1.1 Sources of Uncertainty in Hydrologic Modeling](#11-sources-of-uncertainty-in-hydrologic-modeling)
  * [1.2 Traditional Calibration Approaches](#12-traditional-calibration-approaches)
  * [1.3 Motivation for Bayesian Calibration](#13-motivation-for-bayesian-calibration)
  * [1.4 Conceptual View of Bayesian Learning](#14-conceptual-view-of-bayesian-learning)
* [2. Foundations of Bayesian Theory](#2-foundations-of-bayesian-theory)
  * [2.1 Bayes' Theorem](#21-bayes-theorem)
  * [2.2 Interpretation of Bayesian Components](#22-interpretation-of-bayesian-components)
  * [2.3 Bayesian Updating as Scientific Reasoning](#23-bayesian-updating-as-scientific-reasoning)
  * [2.4 Bayesian Inference in Hydrologic Context](#24-bayesian-inference-in-hydrologic-context)
  * [2.5 Visualization of Bayesian Updating](#25-visualization-of-bayesian-updating)
  * [2.6 Mathematical Properties of Bayesian Inference](#26-mathematical-properties-of-bayesian-inference)
* [3. Components of Bayesian Calibration](#3-components-of-bayesian-calibration)
  * [3.1 Hydrologic Model Formulation](#31-hydrologic-model-formulation)
  * [3.2 Prior Distribution](#32-prior-distribution)
  * [3.3 Likelihood Function](#33-likelihood-function)
  * [3.4 Homoscedastic and Independent Gaussian Likelihood](#34-homoscedastic-and-independent-gaussian-likelihood)
    * [3.4.1 Likelihood Function](#341-likelihood-function)
    * [3.4.2 Interpretation](#342-interpretation)
  * [3.5 Heteroscedastic and Autocorrelated Gaussian Likelihood](#35-heteroscedastic-and-autocorrelated-gaussian-likelihood)
    * [3.5.1 Heteroscedastic Variance](#351-heteroscedastic-variance)
    * [3.5.2 Autocorrelated Residuals (AR(1))](#352-autocorrelated-residuals-ar1)
    * [3.5.3 Log-Likelihood Function](#353-log-likelihood-function)
    * [3.5.4 Understanding the Jacobian Adjustment](#354-understanding-the-jacobian-adjustment)
      * [3.5.4.1 The Change-of-Variables Problem](#3541-the-change-of-variables-problem)
      * [3.5.4.2 Jacobian for Heteroscedastic Standardization](#3542-jacobian-for-heteroscedastic-standardization)
      * [3.5.4.3 Why the Jacobian Matters](#3543-why-the-jacobian-matters)
      * [3.5.4.4 Practical Hydrologic Interpretation](#3544-practical-hydrologic-interpretation)
      * [3.5.4.5 Implementation Note](#3545-implementation-note)
    * [3.5.5 Summary of Advanced Likelihood Components](#355-summary-of-advanced-likelihood-components)
  * [3.6 Posterior Distribution](#36-posterior-distribution)
  * [3.7 Role of Error Model Selection](#37-role-of-error-model-selection)
---

# 1. Introduction

Hydrologic models are fundamental tools used to represent and simulate the movement, storage, and transformation of water within natural systems. These models are routinely applied to simulate streamflow generation, snow accumulation and melt, groundwater exchange, soil moisture dynamics, evapotranspiration, and river routing processes. Their applications span flood forecasting, drought risk assessment, reservoir management, climate change impact studies, and water resources planning.

Despite their importance, hydrologic models are inherently uncertain representations of natural systems. The hydrologic cycle is governed by complex interactions among atmospheric, terrestrial, and subsurface processes that are difficult to observe directly and challenging to represent mathematically. As a result, model outputs depend heavily on parameters that are often unknown, difficult to measure, or only indirectly observable. These parameters frequently represent conceptual approximations of real physical processes rather than direct measurements.

Model calibration is the process through which unknown parameter values are estimated so that model simulations reproduce observed hydrologic behavior as closely as possible. Calibration is therefore a central component of hydrologic modeling because it directly influences predictive accuracy and reliability.

---

## 1.1 Sources of Uncertainty in Hydrologic Modeling

Hydrologic modeling involves several distinct sources of uncertainty that influence simulation outcomes. Understanding these uncertainty sources is critical before introducing calibration approaches.

### Parameter Uncertainty

Parameter uncertainty arises because many hydrologic parameters cannot be measured directly at watershed scales. For example, soil hydraulic conductivity varies spatially and temporally and is rarely measured at resolutions compatible with model discretization. Similarly, groundwater recession constants or snowmelt factors often represent lumped conceptual approximations rather than physically measurable properties.

### Forcing Uncertainty

Hydrologic models rely on meteorological forcing data such as precipitation, temperature, radiation, and wind speed. These inputs contain measurement errors, spatial interpolation errors, and representativeness uncertainties. For example, precipitation measurements are often sparse and may not capture spatial variability across complex terrain.

### Structural Uncertainty

Structural uncertainty arises because hydrologic models simplify real-world processes. Model equations may omit certain processes, approximate nonlinear interactions, or assume simplified system behavior. For instance, lumped conceptual models often represent spatially distributed processes using aggregated parameters, which introduces structural limitations.

### Observation Uncertainty

Observed hydrologic data, including streamflow measurements, are subject to rating curve errors, sensor inaccuracies, and sampling limitations. These errors influence calibration because the observational dataset itself contains uncertainty.

---

## 1.2 Traditional Calibration Approaches

Traditional hydrologic calibration methods typically rely on optimization techniques that seek parameter values minimizing a selected performance metric. Common objective functions include the Nash–Sutcliffe Efficiency (NSE), Kling–Gupta Efficiency (KGE), or Root Mean Square Error (RMSE). These approaches treat calibration as a deterministic optimization problem and produce a single parameter set that provides the best agreement between simulated and observed data according to the chosen metric.

Although widely used, deterministic calibration methods suffer from important conceptual limitations. Hydrologic systems commonly exhibit a phenomenon known as **equifinality**, where multiple parameter combinations produce similar simulation performance. Deterministic calibration approaches ignore this multiplicity of acceptable solutions by selecting only one parameter set, thereby underrepresenting uncertainty in model predictions.

Another limitation of deterministic calibration is the inability to incorporate prior scientific knowledge in a formal and quantitative manner. Although modelers may restrict parameter ranges manually, traditional calibration methods lack a coherent framework for integrating prior information with observational data.

---

## 1.3 Motivation for Bayesian Calibration

Bayesian calibration addresses the limitations of deterministic calibration by treating model parameters as uncertain quantities described by probability distributions rather than fixed values. Instead of identifying a single optimal parameter set, Bayesian calibration estimates the entire range of plausible parameter values and quantifies their associated uncertainties.

The Bayesian framework provides several important advantages for hydrologic modeling. First, it enables the formal incorporation of prior knowledge derived from field measurements, regional studies, or expert judgment. Second, it allows uncertainty in model parameters to be propagated through the model to generate probabilistic predictions. Third, Bayesian calibration provides a mathematically consistent approach for incorporating realistic residual error models that better represent hydrologic system behavior.

By adopting Bayesian inference, hydrologic calibration becomes a probabilistic learning process in which model parameters are updated as new observational data becomes available.

---

## 1.4 Conceptual View of Bayesian Learning

Bayesian inference can be understood as a process of updating scientific beliefs. Initially, the modeler possesses some level of knowledge about parameter values based on theoretical understanding or empirical evidence. As observational data are introduced, this initial knowledge is updated using mathematical rules that balance prior knowledge with data-driven evidence.

This learning process is iterative and cumulative. As additional observations become available, parameter distributions can be further refined. This adaptive nature of Bayesian inference makes it particularly suitable for hydrologic systems where new data may become available through extended monitoring programs or improved measurement technologies.

---

# 2. Foundations of Bayesian Theory

Bayesian inference provides the theoretical foundation for probabilistic calibration. It offers a rigorous mathematical framework that describes how information from observational data modifies prior knowledge about unknown parameters.

The central principle underlying Bayesian theory is Bayes' theorem, which defines the relationship between prior knowledge, observational data, and updated parameter beliefs.

---

## 2.1 Bayes' Theorem

Bayes' theorem is expressed as:

$$
p(\theta \mid D) = \frac{p(D \mid \theta) \, p(\theta)}{p(D)}
$$

In this expression, $\theta$ represents the set of model parameters and $D$ represents the observed data. The left-hand side of the equation, $p(\theta \mid D)$, is known as the posterior distribution and represents the updated knowledge about parameter values after accounting for observational data.

The posterior distribution is proportional to the product of two key components: the likelihood function and the prior distribution. The denominator, known as the evidence or marginal likelihood, serves as a normalization constant ensuring that the posterior distribution integrates to one.

---

## 2.2 Interpretation of Bayesian Components

Understanding Bayesian inference requires interpreting each component of Bayes' theorem in both mathematical and conceptual terms.

### Prior Distribution

The prior distribution, denoted as $p(\theta)$, represents knowledge about parameter values before considering observational data. Priors can be derived from physical understanding, field experiments, literature values, or regional hydrologic studies. In hydrologic modeling, priors help ensure that parameter estimates remain physically realistic and consistent with known system behavior.

### Likelihood Function

The likelihood function, denoted as $p(D \mid \theta)$, quantifies how well model simulations generated using parameter set $\theta$ reproduce observed data. The likelihood measures the probability of observing the data given a particular parameter combination. Parameter sets that produce simulations closely matching observations receive higher likelihood values.

### Posterior Distribution

The posterior distribution, denoted as $p(\theta \mid D)$, represents updated parameter knowledge obtained by combining prior information with observational evidence. The posterior distribution is the primary result of Bayesian calibration and provides a complete probabilistic description of parameter uncertainty.

### Evidence

The evidence, denoted as $p(D)$, is obtained by integrating the product of the likelihood and prior distributions over the entire parameter space:

$$
p(D) = \int p(D \mid \theta) \, p(\theta) \, d\theta
$$

Although essential for normalization, the evidence is often difficult to compute directly. In practical calibration applications, the evidence is typically not evaluated explicitly because sampling-based algorithms operate using proportional relationships.

---

## 2.3 Bayesian Updating as Scientific Reasoning

Bayesian inference formalizes the intuitive scientific method of learning from evidence. Initially, scientists develop hypotheses or expectations regarding system behavior based on theoretical knowledge or prior observations. As new data are collected, these hypotheses are revised to reflect updated understanding. Bayesian inference translates this reasoning into a mathematical framework that systematically combines prior knowledge with empirical observations.

This perspective highlights an important philosophical distinction between Bayesian and classical statistical approaches. Classical inference typically treats parameters as fixed but unknown constants, whereas Bayesian inference treats parameters as random variables that reflect uncertainty in knowledge about their true values.

---

## 2.4 Bayesian Inference in Hydrologic Context

In hydrologic calibration, Bayesian inference is applied by linking hydrologic model outputs to observational data through probabilistic error models. Consider a hydrologic model that simulates streamflow based on meteorological inputs and parameter values. The difference between simulated and observed streamflow is represented using an error model that accounts for measurement uncertainty, model structural limitations, and natural variability.

The Bayesian framework uses this error model to construct the likelihood function, which quantifies how well different parameter values explain the observed hydrograph. The posterior distribution resulting from Bayesian calibration therefore reflects parameter values that are consistent with both prior scientific understanding and observational evidence.

---

## 2.5 Visualization of Bayesian Updating

Bayesian updating can be visualized as the transformation of an initial probability distribution into an updated distribution after observing data. If observational data are highly informative, the posterior distribution typically becomes narrower than the prior distribution, indicating reduced uncertainty. Conversely, if observational data provide limited information about a parameter, the posterior distribution may remain similar to the prior distribution.

In hydrologic calibration, visualization of prior and posterior distributions provides valuable diagnostic insight. For example, a parameter whose posterior distribution remains nearly identical to its prior distribution may indicate poor parameter identifiability or insufficient observational information.

---

## 2.6 Mathematical Properties of Bayesian Inference

Bayesian inference possesses several mathematical properties that make it particularly useful for hydrologic calibration. Bayesian methods naturally propagate uncertainty through nonlinear models because parameter uncertainty is explicitly represented as probability distributions. Bayesian inference also allows joint estimation of multiple parameters while accounting for parameter correlations. Furthermore, Bayesian updating remains consistent when new observational datasets are incorporated, allowing sequential refinement of parameter estimates.

These properties make Bayesian calibration particularly suitable for hydrologic systems characterized by strong nonlinear interactions and limited observational data availability.

---

# 3. Components of Bayesian Calibration

Bayesian calibration in hydrology combines prior knowledge, model structure, and probabilistic error formulations to infer parameter distributions. The likelihood function plays a central role because it formally links hydrologic simulations to observational data. In realistic hydrologic systems, residual errors are rarely simple: they often exhibit flow-dependent variance (heteroscedasticity) and temporal correlation (autocorrelation).

This section describes the **essential components** of Bayesian calibration, including:

1. The hydrologic model structure
2. Prior distributions for parameters
3. Likelihood functions, including:
   * Homoscedastic and independent Gaussian likelihood
   * Heteroscedastic and autocorrelated Gaussian likelihood with Box-Cox transformation
4. Interpretation of each component in hydrologic terms

---

## 3.1 Hydrologic Model Formulation

A hydrologic model represents the relationship between system inputs, internal states, and simulated outputs:

$$
\hat{Q}_t = f(\theta, X_t, S_0)
$$

where:

* $\hat{Q}_t$ is simulated streamflow at time $t$
* $\theta$ is the set of model parameters (e.g., soil, snow, groundwater parameters)
* $X_t$ represents time-varying forcings, such as precipitation and temperature
* $S_0$ represents initial model states, such as soil moisture or snowpack

Observed flows are denoted as:

$$
Q_t^{obs}
$$

Residual errors represent discrepancies between simulated and observed flows:

$$
\varepsilon_t = Q_t^{obs} - \hat{Q}_t
$$

Residual errors arise from measurement errors, forcing uncertainty, model structural simplifications, and inherent variability in hydrologic systems.

---

## 3.2 Prior Distribution

The prior distribution describes initial knowledge about model parameters:

$$
\theta \sim p(\theta)
$$

In hydrology, priors help constrain parameters within physically reasonable ranges. Examples:

* Soil hydraulic conductivity may be log-normally distributed because it is strictly positive.
* Groundwater recession coefficients may be beta-distributed because they lie in the interval [0,1].
* Snowmelt parameters may follow uniform or truncated normal distributions if prior knowledge is limited.

Priors ensure that posterior parameter estimates remain realistic and physically interpretable. When data are highly informative, the posterior will be dominated by the likelihood; otherwise, the prior may have significant influence.

---

## 3.3 Likelihood Function

The likelihood function describes the probability of observing streamflow data given model parameters and assumptions about residual errors:

$$
p(D \mid \theta)
$$

Residual error behavior is crucial:

1. Residual variance may change with flow magnitude (**heteroscedasticity**)
2. Residuals may exhibit temporal correlation (**autocorrelation**)

We present two likelihood formulations, from simple to advanced.

---

## 3.4 Homoscedastic and Independent Gaussian Likelihood

The simplest assumption is that residuals are independent, identically distributed, and have constant variance:

$$
\varepsilon_t \sim N(0, \sigma^2)
$$

### 3.4.1 Likelihood Function

$$
p(D \mid \theta, \sigma) = \prod_{t=1}^{T} \frac{1}{\sqrt{2 \pi \sigma^2}} \exp\left( -\frac{(Q_t^{obs} - \hat{Q}_t)^2}{2 \sigma^2} \right)
$$

Or equivalently, the log-likelihood:

$$
\log p(D \mid \theta, \sigma) = -\frac{T}{2} \log(2 \pi \sigma^2) - \frac{1}{2\sigma^2} \sum_{t=1}^{T} (Q_t^{obs} - \hat{Q}_t)^2
$$

### 3.4.2 Interpretation

* All residuals are treated equally across flow magnitudes.
* Low-flow errors may dominate calibration because numerical values are smaller.
* This assumption often underestimates high-flow uncertainty and ignores temporal dependence.

---

## 3.5 Heteroscedastic and Autocorrelated Gaussian Likelihood

Real hydrologic residuals often **violate homoscedasticity and independence assumptions**. Observed errors are flow-dependent and temporally correlated. To address this, we define:

1. **Flow-dependent (heteroscedastic) residual variance**
2. **Standardized residuals**
3. **AR(1) autocorrelation** in standardized residuals
4. **Jacobian adjustment** for the standardization transformation

---

### 3.5.1 Heteroscedastic Variance

The flow-dependent residual standard deviation is:

$$
\sigma_{\varepsilon(t)} = a \cdot \hat{Q}_t + b
$$

* $a$ scales error with flow magnitude
* $b$ is the baseline variance (intercept)

Raw residuals are:

$$
\varepsilon_t = Q_t^{obs} - \hat{Q}_t
$$

Standardized residuals are:

$$
\eta_t = \frac{\varepsilon_t}{\sigma_{\varepsilon(t)}} = \frac{Q_t^{obs} - \hat{Q}_t}{\sigma_{\varepsilon(t)}}
$$

This standardization accounts for heteroscedasticity and allows the AR(1) model to assume constant innovation variance.

---

### 3.5.2 Autocorrelated Residuals (AR(1))

Standardized residuals are modeled as an AR(1) process:

$$
\eta_1 \sim N\left(0, \frac{1}{1 - \phi^2}\right)
$$

$$
\eta_t \mid \eta_{t-1} \sim N(\phi \, \eta_{t-1}, 1), \quad t = 2, \dots, T
$$

* $\phi$ is the autocorrelation coefficient
* Innovation variance is 1 due to standardization
* The variance of $\eta_1$ is $1/(1-\phi^2)$ for stationarity

This formulation captures **memory effects in streamflow** and flow-dependent uncertainty.

---

### 3.5.3 Log-Likelihood Function

The log-likelihood function for the heteroscedastic AR(1) error model is:

$$
\begin{align}
\log p(D \mid \theta, a, b, \phi) = &-\frac{1}{2} \log \left( \frac{2\pi}{1 - \phi^2} \right) - \frac{\eta_1^2 (1 - \phi^2)}{2} \\
&- \frac{1}{2} \sum_{t=2}^{T} \left[ \log(2\pi) + (\eta_t - \phi \, \eta_{t-1})^2 \right] \\
&+ \sum_{t=1}^{T} \log \left( \frac{1}{\sigma_{\varepsilon(t)}} \right)
\end{align}
$$

where the last term is the **Jacobian adjustment** for the standardization transformation. The Jacobian adjustment is a critical component of the likelihood function when working with standardized residuals. This section explains **why** the Jacobian is necessary and **how** to derive it correctly for heteroscedastic error models.

This can be rewritten more explicitly as:

$$
\begin{align}
\log p(D \mid \theta, a, b, \phi) = &-\frac{1}{2} \log \left( \frac{2\pi}{1 - \phi^2} \right) - \frac{\eta_1^2 (1 - \phi^2)}{2} \\
&- \frac{T-1}{2} \log(2\pi) - \frac{1}{2} \sum_{t=2}^{T} (\eta_t - \phi \, \eta_{t-1})^2 \\
&- \sum_{t=1}^{T} \log(\sigma_{\varepsilon(t)})
\end{align}
$$

**Components**:
* First line: Initial condition for AR(1) process
* Second line: AR(1) contributions for $t = 2, \ldots, T$
* Third line: Jacobian for heteroscedastic standardization
---

#### 3.5.4.1 The Change-of-Variables Problem

When we transform random variables, probability densities must be adjusted to maintain mathematical validity. Consider a random variable $Y$ with density $p_Y(y)$ and its transformation $Z = g(Y)$ with density $p_Z(z)$. These densities are related by:

$$
p_Y(y) = p_Z(g(y)) \left| \frac{dg(y)}{dy} \right|
$$

The term $\left| \frac{dg(y)}{dy} \right|$ is the **Jacobian** of the transformation. It accounts for how the transformation stretches or compresses probability mass across different regions of the variable space.

**Key principle**: When we write a likelihood in terms of transformed variables, we must include the Jacobian to correctly represent the probability of observing the **original data**.

---

#### 3.5.4.2 Jacobian for Heteroscedastic Standardization

The heteroscedastic error model defines raw residuals:

$$
\varepsilon_t = Q_t^{obs} - \hat{Q}_t
$$

with flow-dependent standard deviation:

$$
\sigma_{\varepsilon(t)} = a \cdot \hat{Q}_t + b
$$

We then standardize residuals to obtain unit-variance innovations:

$$
\eta_t = \frac{\varepsilon_t}{\sigma_{\varepsilon(t)}}
$$

**The transformation**: We are transforming from raw residuals $\boldsymbol{\varepsilon} = (\varepsilon_1, \ldots, \varepsilon_T)$ to standardized residuals $\boldsymbol{\eta} = (\eta_1, \ldots, \eta_T)$.

**Deriving the Jacobian**: The inverse transformation is:

$$
\varepsilon_t = \sigma_{\varepsilon(t)} \cdot \eta_t
$$

The derivative of $\varepsilon_t$ with respect to $\eta_t$ is:

$$
\frac{\partial \varepsilon_t}{\partial \eta_t} = \sigma_{\varepsilon(t)}
$$

For the multivariate transformation, the Jacobian matrix is diagonal:

$$
\frac{\partial \boldsymbol{\varepsilon}}{\partial \boldsymbol{\eta}} = \text{diag}(\sigma_{\varepsilon(1)}, \ldots, \sigma_{\varepsilon(T)})
$$

The determinant is:

$$
\left| \det \left( \frac{\partial \boldsymbol{\varepsilon}}{\partial \boldsymbol{\eta}} \right) \right| = \prod_{t=1}^{T} \sigma_{\varepsilon(t)}
$$

**Likelihood with Jacobian**: The joint density of raw residuals is:

$$
p(\boldsymbol{\varepsilon}) = p(\boldsymbol{\eta}) \left| \det \left( \frac{\partial \boldsymbol{\varepsilon}}{\partial \boldsymbol{\eta}} \right) \right| = p(\boldsymbol{\eta}) \prod_{t=1}^{T} \sigma_{\varepsilon(t)}
$$

In log-likelihood form:

$$
\log p(\boldsymbol{\varepsilon}) = \log p(\boldsymbol{\eta}) + \sum_{t=1}^{T} \log \sigma_{\varepsilon(t)}
$$

**Practical implication**: When we model standardized residuals $\eta_t$ using an AR(1) process with unit innovation variance, the likelihood of the **original residuals** must include the Jacobian term $\sum_{t=1}^{T} \log \sigma_{\varepsilon(t)}$.

Since we want the **negative** log-likelihood (for minimization or to subtract from the log-likelihood), this appears as:

$$
-\sum_{t=1}^{T} \log \sigma_{\varepsilon(t)}
$$

in the complete log-likelihood expression in Section 3.5.3.

---

#### 3.5.4.3 Why the Jacobian Matters

Without the Jacobian adjustment:

1. **Parameter estimates become biased** because the likelihood does not correctly represent the data-generating process
2. **Variance parameters $a$ and $b$ are incorrectly estimated**, leading to wrong uncertainty quantification
3. **Model comparison becomes invalid** because likelihoods are evaluated on different scales
4. **Predictive uncertainty is misrepresented**, with credible intervals not achieving nominal coverage
5. **The optimization problem is mathematically incorrect**

---

#### 3.5.4.4 Practical Hydrologic Interpretation

In hydrologic calibration:

* **Heteroscedastic standardization** ensures that low-flow and high-flow residuals contribute appropriately to the likelihood.

* The Jacobian term $-\sum_{t=1}^{T} \log \sigma_{\varepsilon(t)}$ **penalizes parameter values** that predict unrealistically large error variance. This is a form of automatic regularization.

* **Physical interpretation**: When $\sigma_{\varepsilon(t)}$ is large (high uncertainty), the Jacobian contribution becomes more negative, reducing the likelihood. This encourages the calibration to find parameters that minimize prediction uncertainty.

* **Combined effect**: The standardization removes heteroscedasticity from the AR(1) model (making it have unit innovation variance), while the Jacobian ensures we're still making valid probabilistic statements about the **original streamflow observations**.

---

#### 3.5.4.5 Implementation Note

In practice, the Jacobian term appears naturally in the log-likelihood. Using the relationship $\sigma_{\varepsilon(t)} = a \cdot \hat{Q}_t + b$:

```r
# Compute flow-dependent standard deviations
sigma_eps <- a * Q_sim + b

# Jacobian contribution to log-likelihood
jacobian_term <- -sum(log(sigma_eps))

# Or equivalently, when computing the full likelihood:
log_likelihood <- log_p_eta + jacobian_term
```

where `log_p_eta` is the log-density of the AR(1) process for standardized residuals.

**Critical reminder**: This Jacobian must **always** be included when standardizing residuals. Omitting it invalidates the entire Bayesian inference procedure.

---

### 3.5.5 Summary of Advanced Likelihood Components

The heteroscedastic AR(1) likelihood addresses key limitations of simpler error models:

* **Heteroscedastic variance** ($\sigma_{\varepsilon(t)} = a \cdot \hat{Q}_t + b$) ensures proper weighting of residuals across low and high flows.
* **Standardization** ($\eta_t = \varepsilon_t / \sigma_{\varepsilon(t)}$) transforms residuals to unit variance, enabling standard AR(1) modeling.
* **AR(1) autocorrelation** captures temporal persistence due to storage and routing processes in hydrologic systems.
* **Jacobian adjustment** ($-\sum \log \sigma_{\varepsilon(t)}$) maintains statistical validity when working with standardized residuals.

Calibration using this likelihood leads to **realistic posterior distributions** of model parameters, proper uncertainty quantification across flow regimes, and more reliable predictive uncertainty.

---

## 3.6 Posterior Distribution

Combining prior distributions with the likelihood yields the posterior:

$$
p(\theta \mid D) \propto p(D \mid \theta) \, p(\theta)
$$

Including the heteroscedastic AR(1) error model introduces additional parameters for the residual structure:

$$
p(\theta, a, b, \phi \mid D)
$$

where:
* $\theta$ represents the hydrologic model parameters
* $a, b$ are the heteroscedastic variance parameters
* $\phi$ is the AR(1) autocorrelation coefficient

The posterior distribution describes uncertainty in both **hydrologic model parameters** and **residual error model parameters**, enabling fully probabilistic predictions that account for both parameter uncertainty and realistic error structure.

---

## 3.7 Role of Error Model Selection

* The choice of likelihood directly affects calibration results.
* Simple Gaussian likelihoods may **underestimate uncertainty** in high flows and ignore temporal dependence.
* Complex likelihoods improve realism but increase computational cost and parameter identifiability challenges.
* Diagnostics such as residual analysis, autocorrelation plots, and posterior predictive checks are essential.

---
