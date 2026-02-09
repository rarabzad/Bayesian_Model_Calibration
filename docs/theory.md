# Bayesian Calibration for Hydrologists

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
* [4. Bayesian Linear Regression — A Fully Worked Calibration Example](#4-bayesian-linear-regression--a-fully-worked-calibration-example)
  * [4.1 Problem Definition and Model Structure](#41-problem-definition-and-model-structure)
  * [4.2 Synthetic Data Generation](#42-synthetic-data-generation)
  * [4.3 Bayesian Model Specification](#43-bayesian-model-specification)
    * [4.3.1 Prior Distributions](#431-prior-distributions)
    * [4.3.2 Likelihood Function](#432-likelihood-function)
    * [4.3.3 Posterior Distribution](#433-posterior-distribution)
  * [4.4 MCMC Sampling](#44-mcmc-sampling)
  * [4.5 Convergence Diagnostics](#45-convergence-diagnostics)
  * [4.6 Posterior Summaries](#46-posterior-summaries)
  * [4.7 Joint Posterior Structure](#47-joint-posterior-structure)
  * [4.8 Posterior Predictive Checks](#48-posterior-predictive-checks)
  * [4.9 Residual Analysis](#49-residual-analysis)
  * [4.10 Summary and Learning Outcomes](#410-summary-and-learning-outcomes)
* [5. Bayesian Calibration of GR4J with Homoscedastic Gaussian Errors](#5-bayesian-calibration-of-gr4j-with-homoscedastic-gaussian-errors)
  * [Overview of the GR4J Model](#overview-of-the-gr4j-model)
  * [Why Bayesian Calibration?](#why-bayesian-calibration)
  * [5.1 Model Parameters and Bounds](#51-model-parameters-and-bounds)
    * [Parameter list and interpretation](#parameter-list-and-interpretation)
    * [Detailed Parameter Interpretation](#detailed-parameter-interpretation)
      * [X1: Production Store Capacity (mm)](#x1-production-store-capacity-mm)
      * [X2: Groundwater Exchange Coefficient (dimensionless)](#x2-groundwater-exchange-coefficient-dimensionless)
      * [X3: Routing Store Capacity (mm)](#x3-routing-store-capacity-mm)
      * [X4: Unit Hydrograph Time Base (days)](#x4-unit-hydrograph-time-base-days)
      * [TT: Temperature Threshold (°C)](#tt-temperature-threshold-c)
      * [DDF: Degree-Day Factor (mm °C⁻¹ day⁻¹)](#ddf-degree-day-factor-mm-c1-day1)
      * [σ: Residual Standard Deviation (m³/s)](#σ-residual-standard-deviation-ms)
    * [Recommended example bounds (adjust to local catchment)](#recommended-example-bounds-adjust-to-local-catchment)
    * [Scientific rationale for bounds](#scientific-rationale-for-bounds)
    * [Parameter Transformation Strategies](#parameter-transformation-strategies)
  * [5.2 Log‑Prior Function](#52-logprior-function)
    * [Uniform prior (default)](#uniform-prior-default)
    * [Alternative priors and when to use them](#alternative-priors-and-when-to-use-them)
    * [Prior predictive checks](#prior-predictive-checks)
  * [5.3 Likelihood Function](#53-likelihood-function)
    * [Default: Homoscedastic Gaussian](#default-homoscedastic-gaussian)
    * [Likelihood alternatives](#likelihood-alternatives)
    * [Likelihood selection checklist](#likelihood-selection-checklist)
  * [5.4 Log‑Posterior Function](#54-logposterior-function)
    * [Scientific considerations](#scientific-considerations)
  * [5.5 MCMC Sampling](#55-mcmc-sampling)
    * [Example using adaptMCMC](#example-using-adaptmcmc)
    * [MCMC best practices](#mcmc-best-practices)
    * [Computational tips](#computational-tips)
  * [5.6 Posterior Predictive Checks (PPCs)](#56-posterior-predictive-checks-ppcs)
    * [Theoretical Foundation](#theoretical-foundation)
    * [Procedure](#procedure)
    * [R code (ensemble PPC)](#r-code-ensemble-ppc)
    * [Visual PPC Diagnostics](#visual-ppc-diagnostics)
      * [Time Series Plot with Uncertainty Bands](#time-series-plot-with-uncertainty-bands)
      * [Spaghetti Plot: Overlay of Multiple Realizations](#spaghetti-plot-overlay-of-multiple-realizations)
      * [Flow Duration Curve (FDC) Comparison](#flow-duration-curve-fdc-comparison)
    * [PPC diagnostics and interpretation](#ppc-diagnostics-and-interpretation)
      * [Coverage Analysis](#coverage-analysis)
      * [Conditional Coverage: Coverage by Flow Regime](#conditional-coverage-coverage-by-flow-regime)
      * [Tail Behavior: Extreme Event Capture](#tail-behavior-extreme-event-capture)
    * [Predictive Performance Metrics](#predictive-performance-metrics)
    * [Seasonal and Regime-Specific PPCs](#seasonal-and-regime-specific-ppcs)
    * [Advanced PPC: Test Statistics](#advanced-ppc-test-statistics)
    * [Reporting PPC Results](#reporting-ppc-results)
    * [When PPCs Fail: Next Steps](#when-ppcs-fail-next-steps)
  * [5.7 Residual Diagnostics](#57-residual-diagnostics)
    * [Residual definition](#residual-definition)
    * [5.7.1 Residual vs Fitted](#571-residual-vs-fitted)
    * [5.7.2 Residual vs Time and Autocorrelation](#572-residual-vs-time-and-autocorrelation)
    * [5.7.3 QQ Plot and Histogram](#573-qq-plot-and-histogram)
      * [QQ Plot](#qq-plot)
      * [Histogram](#histogram)
    * [Shapiro-Wilk Test for Normality](#shapiro-wilk-test-for-normality)
    * [5.7.4 Residual Variance by Flow Regime](#574-residual-variance-by-flow-regime)
    * [5.7.5 Seasonal Residual Analysis](#575-seasonal-residual-analysis)
    * [Quantitative diagnostics to report](#quantitative-diagnostics-to-report)
    * [Actionable responses to diagnostics](#actionable-responses-to-diagnostics)
    * [Advanced Diagnostic: Recursive Residuals](#advanced-diagnostic-recursive-residuals)
    * [Residual Diagnostics in Model Selection](#residual-diagnostics-in-model-selection)
    * [Summary: Iterative Refinement](#summary-iterative-refinement)---

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
2. **Prior distributions for parameters**
3. **Likelihood functions**, may include:
   * Homoscedastic and independent Gaussian likelihood
   * Heteroscedastic and autocorrelated Gaussian likelihood with Box-Cox transformation
4. **Posterior distribution**
5. Interpretation of each component in hydrologic terms
![Bayesian Workflow](https://svmiller.com/images/prior-posterior-likelihood.png)

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
![Descriptive alt text](https://ars.els-cdn.com/content/image/1-s2.0-S0022169424011338-gr2.jpg)
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

where the last term is the **Jacobian adjustment** for the standardization transformation.

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

### 3.5.4 Understanding the Jacobian Adjustment

The Jacobian adjustment is a critical component of the likelihood function when working with standardized residuals. This section explains **why** the Jacobian is necessary and **how** to derive it correctly for heteroscedastic error models.

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

**Deriving the Jacobian**: The forward transformation is:

$$
\eta_t = \frac{\varepsilon_t}{\sigma_{\varepsilon(t)}}
$$

The derivative of $\eta_t$ with respect to $\varepsilon_t$ is:

$$
\frac{\partial \eta_t}{\partial \varepsilon_t} = \frac{1}{\sigma_{\varepsilon(t)}}
$$

For the multivariate transformation, the Jacobian matrix is diagonal:

$$
\frac{\partial \boldsymbol{\eta}}{\partial \boldsymbol{\varepsilon}} = \text{diag}\left(\frac{1}{\sigma_{\varepsilon(1)}}, \ldots, \frac{1}{\sigma_{\varepsilon(T)}}\right)
$$

The determinant is:

$$
\left| \det \left( \frac{\partial \boldsymbol{\eta}}{\partial \boldsymbol{\varepsilon}} \right) \right| = \prod_{t=1}^{T} \frac{1}{\sigma_{\varepsilon(t)}}
$$

**Likelihood with Jacobian**: The joint density of standardized residuals is related to the density of raw residuals by:

$$
p(\boldsymbol{\eta}) = p(\boldsymbol{\varepsilon}) \left| \det \left( \frac{\partial \boldsymbol{\eta}}{\partial \boldsymbol{\varepsilon}} \right) \right| = p(\boldsymbol{\varepsilon}) \prod_{t=1}^{T} \frac{1}{\sigma_{\varepsilon(t)}}
$$

Rearranging to express the density of raw residuals (which is what we actually observe):

$$
p(\boldsymbol{\varepsilon}) = p(\boldsymbol{\eta}) \prod_{t=1}^{T} \sigma_{\varepsilon(t)}
$$

In log-likelihood form:

$$
\log p(\boldsymbol{\varepsilon}) = \log p(\boldsymbol{\eta}) + (-\sum_{t=1}^{T} \log \sigma_{\varepsilon(t)})
$$

**Practical implication**: When we model standardized residuals $\eta_t$ using an AR(1) process with unit innovation variance, the likelihood of the **original residuals** must include the Jacobian term $\sum_{t=1}^{T} \log \sigma_{\varepsilon(t)}$.

So we have:

$$
-\sum_{t=1}^{T} \log \sigma_{\varepsilon(t)} = -\sum_{t=1}^{T} \log(a \cdot \hat{Q}_t + b)
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

# 4. Bayesian Linear Regression — A Fully Worked Calibration Example

This section presents a complete Bayesian calibration example using a simple linear regression model. Although linear regression is mathematically simple, it contains every fundamental component required for Bayesian hydrologic calibration, including parameter inference, prior specification, likelihood construction, posterior sampling, convergence assessment, uncertainty quantification, and predictive validation.

The objective of this example is educational rather than application-specific. By working with a simple model where relationships are transparent, participants can focus on understanding Bayesian inference mechanics without being overwhelmed by hydrologic model complexity.

Importantly, linear regression is structurally equivalent to many surrogate hydrologic models and provides an ideal foundation for understanding how Bayesian calibration extends to rainfall–runoff models, snow models, and groundwater models.

---

## 4.1 Problem Definition and Model Structure

We consider the following statistical model:

$$
y = \beta_0 + \beta_1 x + \varepsilon
$$

where:

* $y$ represents observed system responses
* $x$ represents a known predictor or forcing variable
* $\beta_0$ is the intercept parameter
* $\beta_1$ is the slope parameter
* $\varepsilon$ represents residual error

Residual errors are assumed Gaussian:

$$
\varepsilon \sim \mathcal{N}(0, \sigma^2)
$$

This assumption implies two important properties:

### Homoscedasticity

Residual variance is constant across the range of predictions.

### Independence

Residuals are not temporally correlated.

Although these assumptions are often violated in hydrologic modeling, they represent the simplest likelihood formulation and serve as a conceptual baseline.

---

### Conceptual Link to Hydrologic Calibration

The linear model can be interpreted analogously to a rainfall–runoff model:

| Linear Model | Hydrologic Model                     |
|--------------|--------------------------------------|
| $x$          | Precipitation / forcing              |
| $y$          | Streamflow                           |
| $\beta_1$    | Catchment response coefficient       |
| $\beta_0$    | Baseflow or offset term              |
| $\sigma$     | Model structural + measurement error |

Understanding inference for this simple system provides intuition for complex hydrologic calibration workflows.

---

### Calibration Objective

Instead of estimating a single best parameter value, Bayesian calibration aims to estimate the full parameter distribution:

$$
\theta = (\beta_0, \beta_1, \sigma)
$$

These posterior distributions quantify uncertainty in model structure, data, and parameter identifiability.

---

## 4.2 Synthetic Data Generation

To validate Bayesian calibration performance, we generate synthetic observations using known parameter values. This approach allows direct evaluation of whether posterior inference correctly recovers true parameters.

Synthetic data experiments are a cornerstone of statistical method validation. In real-world hydrologic applications, we never know the true parameter values—this is precisely why we calibrate models. However, when learning or testing calibration methods, working with synthetic data provides critical advantages:

1. **Ground truth verification**: We can directly compare posterior estimates to known true values
2. **Algorithm debugging**: If the method fails to recover known parameters, we know something is wrong with our implementation
3. **Uncertainty assessment**: We can evaluate whether credible intervals have correct coverage properties
4. **Sensitivity analysis**: We can systematically vary data quality, sample size, or noise levels to understand method performance

In hydrologic modeling, synthetic experiments help separate three distinct sources of uncertainty: (1) parameter estimation uncertainty, (2) model structural error, and (3) measurement error. By controlling the data-generating process, we can isolate and study each component independently.

---

### True Parameter Values

We define the following true parameter values that will generate our synthetic observations:

$$
\beta_0^{\text{true}} = 2.5
$$

$$
\beta_1^{\text{true}} = 1.8
$$

$$
\sigma^{\text{true}} = 1.2
$$

These values define the true data-generating process. In the context of our linear model $y = \beta_0 + \beta_1 x + \varepsilon$, these parameters have clear interpretations:

* $\beta_0^{\text{true}} = 2.5$ represents the baseline response when the predictor is zero (analogous to baseflow in hydrology)
* $\beta_1^{\text{true}} = 1.8$ represents the sensitivity of the response to changes in the predictor (analogous to a rainfall-runoff coefficient)
* $\sigma^{\text{true}} = 1.2$ represents the magnitude of irreducible uncertainty in predictions (combining measurement error and model structural limitations)

The choice of these specific values is somewhat arbitrary, but they are selected to represent realistic relationships where:
- The intercept is positive (physically meaningful for many hydrologic processes)
- The slope is positive and greater than 1 (representing amplification of the input signal)
- The residual standard deviation is moderate relative to the signal strength (signal-to-noise ratio is reasonable)

---

### Data Generation Procedure

We simulate observations through a two-step stochastic process that mimics real data collection in hydrologic systems.

#### Step 1: Generate Forcing Variable

First, we generate predictor values uniformly distributed across the domain:

$$
x_i \sim \text{Uniform}(0,10), \quad i = 1, \ldots, n
$$

This represents evenly distributed predictor values across the observation space. In hydrologic terms, this could represent precipitation events spanning a range of intensities, or temperature variations across seasons. 

The uniform distribution ensures that our data provide information across the entire range of forcing conditions. This is important for parameter identifiability—if all observations occurred at similar predictor values, the slope parameter would be poorly constrained. In real hydrologic calibration, this corresponds to the importance of having diverse hydrologic conditions (wet and dry periods, high and low flows) in the calibration dataset.

We generate $n = 50$ observations, which represents a moderately sized dataset—large enough to provide reasonable parameter constraints, but small enough that uncertainty remains substantial. This sample size is comparable to, for example, a few years of monthly streamflow observations or one season of daily observations.

#### Step 2: Generate Observations with Noise

Next, we generate the response variable by combining the deterministic model prediction with random noise:

$$
y_i = 2.5 + 1.8 \, x_i + \varepsilon_i
$$

$$
\varepsilon_i \sim \mathcal{N}(0, 1.2^2)
$$

This two-component structure reflects the fundamental nature of hydrologic observations:

**Deterministic component** ($2.5 + 1.8 \, x_i$): Represents the systematic relationship between forcing and response. In a perfect world with no measurement error and a perfect model, this would be the exact observed value.

**Stochastic component** ($\varepsilon_i$): Represents the combined effects of:
- Measurement uncertainty (e.g., rating curve errors in streamflow gauging)
- Model structural error (the linear model is a simplification of reality)
- Natural variability not captured by the predictor (e.g., subsurface processes, spatial heterogeneity)
- Timing errors or temporal aggregation effects

The Gaussian noise assumption ($\varepsilon_i \sim \mathcal{N}(0, 1.2^2)$) implies that errors are symmetric, independent, and homoscedastic (constant variance). While these assumptions are often violated in real hydrologic data—errors may be flow-dependent, autocorrelated, or skewed—they provide a tractable starting point for understanding Bayesian calibration mechanics.

The standard deviation $\sigma = 1.2$ creates a coefficient of variation (noise-to-signal ratio) that varies across the predictor range. For small $x$ values (near 0), the predicted response is approximately 2.5, making the relative error about 48%. For large $x$ values (near 10), the predicted response is approximately 20.5, making the relative error about 6%. This heteroscedasticity would violate our model assumptions in a real application and motivate the heteroscedastic error models discussed in Section 3.

---

### Teaching Insight: The Value of Synthetic Experiments

Synthetic experiments allow researchers to separate algorithmic performance from real-world data complications. If posterior inference cannot recover known parameters, either the likelihood, prior specification, or sampling procedure is incorrect.

This diagnostic capability is invaluable during method development. In real applications, when posterior estimates seem unreasonable, we often cannot determine whether the problem lies in:
- Our calibration algorithm
- Our prior assumptions
- Our likelihood specification
- The data quality
- The model structure itself

Synthetic experiments eliminate this ambiguity. If we generate data from a known model and our calibration fails to recover the true parameters, we have clear evidence of a methodological problem rather than a data or model problem.

Moreover, synthetic experiments allow us to conduct **recovery tests**: we can vary the true parameter values, generate multiple datasets, run calibration, and verify that:
1. Posterior means are unbiased estimators of true values
2. Credible intervals contain true values at the nominal coverage rate (e.g., 95% intervals contain truth in 95% of replications)
3. Posterior widths scale appropriately with data informativeness

In the hydrologic literature, synthetic experiments are often called "twin experiments" or "perfect model experiments" and are standard practice for validating data assimilation systems, forecasting methods, and calibration algorithms.

---

### Code: Synthetic Data Generation

```r
# Set random seed for reproducibility
# This ensures that every run produces identical results
set.seed(123)

# Define true parameter values
beta0_true <- 2.5  # Intercept
beta1_true <- 1.8  # Slope
sigma_true <- 1.2  # Residual standard deviation

# Sample size
n <- 50

# Generate predictor values uniformly between 0 and 10
# This creates a well-distributed set of forcing conditions
x <- runif(n, 0, 10)

# Generate response values using the true model plus Gaussian noise
# This simulates the data-generating process
y <- beta0_true + beta1_true * x + rnorm(n, 0, sigma_true)

# Visualize the synthetic data
plot(x, y, 
     pch = 19, 
     col = "steelblue",
     xlab = "Predictor (x)",
     ylab = "Response (y)",
     main = "Synthetic Data with True Regression Line")

# Add the true regression line for reference
# The red dashed line shows the true relationship without noise
abline(a = beta0_true, b = beta1_true, 
       col = "red", lwd = 2, lty = 2)

# Add legend
legend("topleft", 
       legend = c("Observed data", "True relationship"),
       col = c("steelblue", "red"),
       pch = c(19, NA),
       lty = c(NA, 2),
       lwd = c(NA, 2))
```

**Code explanation**:

The `set.seed(123)` command initializes the random number generator to a specific state, ensuring reproducibility. Without this, each run would produce different random samples, making results impossible to replicate.

The `runif(n, 0, 10)` function generates $n$ uniform random variables between 0 and 10. This creates predictor values that span the entire domain uniformly.

The `rnorm(n, 0, sigma_true)` function generates $n$ Gaussian random errors with mean 0 and standard deviation $\sigma_{\text{true}} = 1.2$. These are added to the deterministic predictions to create realistic noisy observations.

The visualization shows both the noisy observations (blue points) and the true underlying relationship (red dashed line). The scatter of points around the line visually represents the magnitude of residual uncertainty $\sigma$. If we ran Bayesian calibration perfectly, the posterior distribution should concentrate near the true parameter values that generated this data.

This simple data generation procedure embodies the same conceptual structure as complex hydrologic simulators: a deterministic process model (the linear equation) forced by inputs (predictor $x$) and corrupted by uncertainty (Gaussian noise), producing observations (response $y$) that we will use for calibration.

---

## 4.3 Bayesian Model Specification

Bayesian inference combines prior knowledge with observational evidence. The model specification requires defining:

1. Prior distributions
2. Likelihood formulation
3. Posterior inference

---

### 4.3.1 Prior Distributions

Priors represent beliefs about parameter values before analyzing observational data. They formalize what we know (or assume) about parameters from sources other than the current calibration dataset.

In Bayesian inference, prior distributions serve multiple crucial functions:

**1. Encoding existing knowledge**: Priors allow us to incorporate information from previous studies, physical constraints, expert judgment, or theoretical considerations. For example, in hydrology we know that hydraulic conductivity must be positive, or that snow albedo must lie between 0 and 1.

**2. Regularization**: Priors prevent parameter estimates from wandering into unrealistic regions of parameter space, especially when data are limited or parameters are weakly identifiable. This is similar to how ridge regression or LASSO add penalty terms in classical statistics.

**3. Handling identifiability**: When multiple parameter combinations produce similar model outputs (equifinality), priors help break symmetries and produce unique solutions by favoring certain parameter regions.

**4. Propagating uncertainty**: Priors explicitly represent parameter uncertainty, which propagates through to posterior predictions. Wide priors admit more uncertainty; narrow priors are more confident.

For our linear regression example, we specify the following prior distributions:

$$
\beta_0 \sim \mathcal{N}(0, 10^2)
$$

$$
\beta_1 \sim \mathcal{N}(0, 10^2)
$$

$$
\sigma \sim \text{HalfNormal}(5)
$$

Let's examine each prior choice in detail.

---

#### Interpretation of Prior Selection

**Prior for $\beta_0$ (Intercept)**:

The normal distribution $\mathcal{N}(0, 10^2)$ centers the prior on zero with a standard deviation of 10. This is a **weakly informative prior** that expresses the following beliefs:

- We have no strong reason to believe the intercept is positive or negative (centered on 0)
- We believe the intercept is likely between -20 and +20 (approximately 2 standard deviations)
- We assign very low probability to extreme values like $\beta_0 = 100$ or $\beta_0 = -100$

This prior is "weak" because it does not strongly constrain the parameter—the standard deviation of 10 is large relative to our true value of 2.5. The data will easily override this prior if they contain sufficient information. However, the prior is not completely "flat" or "uninformative"—it still expresses skepticism about unreasonably large values.

In hydrologic terms, if $\beta_0$ represented baseflow, a prior centered on zero with moderate variance would express uncertainty about baseflow magnitude while ruling out physically impossible negative flows (through truncation) or implausibly large constant contributions.

**Prior for $\beta_1$ (Slope)**:

Similarly, $\mathcal{N}(0, 10^2)$ for $\beta_1$ expresses:

- No prior preference for positive or negative slopes
- Expectation that slopes are moderate (roughly between -20 and +20)
- Strong discounting of extreme sensitivity coefficients

Our true value of 1.8 falls well within the high-probability region of this prior, so we expect the data to easily refine this prior toward the true value.

In hydrologic calibration, a runoff coefficient might be constrained to (0, 1) or (0, 5) depending on the process. Our current prior allows negative values, which could be inappropriate for some physical parameters. In real applications, we might use truncated normal or log-normal priors to enforce positivity.

**Prior for $\sigma$ (Residual Standard Deviation)**:

The Half-Normal distribution is simply the positive half of a normal distribution:

$$
\text{HalfNormal}(\tau) = \mathcal{N}(0, \tau^2) \text{ truncated to } \sigma > 0
$$

For $\tau = 5$, this prior:

- Enforces $\sigma > 0$ (residual variance must be positive)
- Assigns highest density near $\sigma = 0$ (favoring low noise)
- Decays gradually for larger values
- Has approximately 95% probability mass below $\sigma \approx 10$

This prior is regularizing—it gently pushes $\sigma$ toward smaller values, discouraging the model from attributing all misfit to pure noise. However, it's not so strong that it prevents the data from indicating larger uncertainty if warranted.

The Half-Normal is a common choice for scale parameters (standard deviations, variances) because:
- It respects the positivity constraint naturally
- It's weakly informative (broad but not flat)
- It's computationally stable in MCMC sampling
- It's easier to specify than inverse-gamma or other traditional conjugate priors

Our true value $\sigma^{\text{true}} = 1.2$ falls in the high-probability region of this prior, so calibration should successfully recover it.

---

#### Hydrologic Relevance

In hydrologic calibration, priors commonly incorporate multiple sources of knowledge:

**Physical parameter ranges**: Soil porosity must be between 0 and 1, snow albedo between 0 and 1, temperatures above absolute zero, etc. These hard constraints are implemented through truncated or bounded distributions (uniform, beta, truncated normal).

**Field measurements**: If we have measured soil hydraulic conductivity at a few locations, we might use the mean and variance of these measurements to construct a prior. This is especially useful in distributed models where we calibrate spatially varying parameters.

**Regional transfer**: Catchments in similar climates or geology often share parameter ranges. We can use posterior distributions from previously calibrated nearby watersheds as priors for a new catchment. This is called "hierarchical" or "regional" Bayesian modeling.

**Expert knowledge**: Hydrologists may have intuition about reasonable parameter ranges based on experience. For instance, a modeler familiar with snowmelt processes might specify that a degree-day melt factor should fall between 2 and 8 mm/°C/day based on literature and physical understanding.

**Model structure**: Some parameters are scale-dependent or defined only through their ratio with other parameters. Priors help maintain physical consistency (e.g., ensuring that the sum of soil layer thicknesses equals total soil depth).

**Temporal stability**: If calibrating a model sequentially over time windows, the posterior from period 1 can serve as the prior for period 2, implementing a form of recursive Bayesian updating.

Priors help stabilize calibration when:
- Data are limited (short records, sparse spatial coverage)
- Parameters are weakly identifiable (multiple parameters have similar effects on outputs)
- Model structure is uncertain (compensating errors between parameters)
- Measurements are noisy or biased (prior prevents overfitting to noise)

However, priors must be chosen carefully:

**Overly informative priors** can overwhelm the data, preventing the model from learning from observations. This leads to under-calibration.

**Overly vague priors** (e.g., uniform on $(-\infty, \infty)$) may be improper (don't integrate to 1), cause MCMC sampling problems, or allow unrealistic parameter values.

**Informative but wrong priors** can bias results if the prior conflicts with the data. Always perform prior-posterior comparisons to check whether the data updated the prior meaningfully.

A good rule of thumb: priors should be **weakly informative**—they should regularize and prevent absurdity, but not dominate likelihood information when data are reasonably informative.

---

### Code: Prior Distribution

```r
# Log-prior density function
# Returns the log-probability of parameters under the prior
log_prior <- function(theta) {
  beta0 <- theta[1]  # Extract intercept
  beta1 <- theta[2]  # Extract slope
  sigma <- theta[3]  # Extract residual standard deviation
  
  # Enforce positivity constraint on sigma
  # If sigma <= 0, return -Inf (zero probability)
  if(sigma <= 0) return(-Inf)
  
  # Compute log-prior for beta0: N(0, 10^2)
  log_p_beta0 <- dnorm(beta0, mean = 0, sd = 10, log = TRUE)
  
  # Compute log-prior for beta1: N(0, 10^2)
  log_p_beta1 <- dnorm(beta1, mean = 0, sd = 10, log = TRUE)
  
  # Compute log-prior for sigma: HalfNormal(5)
  # HalfNormal is N(0,5^2) truncated to positive values
  # We compute N(0,5^2) in log space and add log(2) correction
  # because we're only using the positive half
  log_p_sigma <- dnorm(sigma, mean = 0, sd = 5, log = TRUE) + log(2)
  
  # Return sum of log-priors (product of priors in original scale)
  # This assumes prior independence: p(beta0, beta1, sigma) = p(beta0) * p(beta1) * p(sigma)
  return(log_p_beta0 + log_p_beta1 + log_p_sigma)
}
```

**Code explanation**:

The function takes a parameter vector `theta = c(beta0, beta1, sigma)` and returns the log-prior density.

**Why log-probabilities?** We work in log-space for numerical stability. Probabilities can be very small (e.g., $10^{-300}$), causing underflow in floating-point arithmetic. Log-probabilities avoid this problem because:
- Products become sums: $\log(ab) = \log(a) + \log(b)$
- Small probabilities become manageable negative numbers
- MCMC algorithms work naturally in log-space (Metropolis-Hastings uses log-probability ratios)

**Constraint handling**: The check `if(sigma <= 0) return(-Inf)` implements the hard constraint $\sigma > 0$. Returning `-Inf` (log-probability of zero) ensures MCMC never accepts negative $\sigma$ values.

**Half-Normal implementation**: The term `+ log(2)` corrects for truncation. A normal $\mathcal{N}(0, \tau^2)$ has 50% of its mass on the negative side. When we truncate to positive values only, we must renormalize: the positive half needs to integrate to 1 instead of 0.5, which corresponds to multiplying the density by 2, or adding $\log(2)$ in log-space.

**Independence assumption**: We sum log-priors because we assume parameters are independent *a priori*. This is a simplification—in reality, we might believe $\beta_0$ and $\beta_1$ are correlated (e.g., intercept and slope often show negative correlation in regression). For this tutorial, independent priors suffice, but hierarchical or multivariate normal priors can model correlation structure when needed.

**Testing the prior**: Before running calibration, it's good practice to:
1. Sample from the prior distribution
2. Run the model forward with prior samples
3. Visualize prior predictive distributions
4. Check whether prior predictions are reasonable

This is called a **prior predictive check** and helps diagnose overly informative or unrealistic priors before observing data.

---

### 4.3.2 Likelihood Function

The likelihood function defines how probable observed data are given parameter values. It is the mathematical bridge between our model predictions and the observed reality, quantifying the degree of agreement between simulations and measurements.

In Bayesian inference, the likelihood plays a central role:

**1. Data-model comparison**: The likelihood measures how well model predictions match observations for a given parameter set. Higher likelihood means better agreement; lower likelihood means poorer fit.

**2. Weighting evidence**: Different parameter values produce different predictions, which in turn have different likelihoods. The likelihood function creates a landscape over parameter space, with peaks at values that explain the data well.

**3. Balancing with priors**: Through Bayes' theorem, the likelihood combines multiplicatively with the prior. Where they agree, the posterior concentrates. Where they conflict, the posterior represents a compromise weighted by their relative strengths.

**4. Encoding assumptions**: The likelihood embodies our assumptions about error structure—are errors Gaussian? Homoscedastic? Independent? These choices profoundly affect calibration results.

For our linear regression model, we assume that observations are conditionally independent and normally distributed around the model predictions:

$$
y_i \sim \mathcal{N}(\beta_0 + \beta_1 x_i, \sigma^2)
$$

This statement says: "Given parameters $\theta = (\beta_0, \beta_1, \sigma)$ and predictor $x_i$, the observed response $y_i$ is normally distributed with mean $\beta_0 + \beta_1 x_i$ and variance $\sigma^2$."

---

#### Role of the Likelihood

The likelihood formalizes model-data agreement. It quantifies how well the simulation reproduces observations while accounting for residual uncertainty.

**Statistical interpretation**: If we hypothetically could repeat the data collection process many times with the same true parameters and forcings, the distribution of observed values would follow the likelihood distribution. The likelihood thus represents our model of the observational process.

**Calibration perspective**: During calibration, we evaluate the likelihood for many parameter values. Parameters that produce predictions far from observations receive low likelihood; parameters that produce predictions close to observations receive high likelihood. The Bayesian posterior concentrates on high-likelihood regions (weighted by prior plausibility).

**Assumption encoding**: Our likelihood assumes:
- **Normality**: Residuals are Gaussian (symmetric, unbounded)
- **Homoscedasticity**: Variance $\sigma^2$ is constant across all predictions
- **Independence**: Knowing residual $\varepsilon_i$ tells us nothing about $\varepsilon_j$ for $i \neq j$
- **Correct functional form**: The linear relationship $\beta_0 + \beta_1 x$ is the true mean

In real hydrologic applications, these assumptions are often violated:

**Violated normality**: Streamflow errors are often right-skewed (large positive errors more common than large negative errors). This motivates transformed likelihoods (log-space, Box-Cox).

**Violated homoscedasticity**: High flows have larger absolute errors than low flows. This motivates heteroscedastic error models like $\sigma_t = a \cdot \hat{Q}_t + b$ (Section 3.5).

**Violated independence**: Hydrologic processes have memory—today's flow depends on yesterday's flow. This motivates autocorrelated error models like AR(1) (Section 3.5).

**Violated functional form**: Real rainfall-runoff relationships are nonlinear, involve thresholds, and include unmodeled processes. This creates structural model error that cannot be fixed by parameter calibration alone.

Despite these limitations, the simple Gaussian likelihood provides a crucial learning foundation. Once we understand Bayesian calibration with simple likelihoods, we can extend to more realistic error models.

In Bayesian hydrology, **likelihood selection determines which flow features are emphasized during calibration**. Different likelihoods weight different parts of the hydrograph differently:

- **Gaussian in original space**: Emphasizes high flows (large absolute errors)
- **Gaussian in log-space**: Emphasizes low flows (large relative errors)
- **Heteroscedastic models**: Balance emphasis across flow regimes
- **Quantile-based likelihoods**: Target specific exceedance probabilities

The choice of likelihood should align with the intended model use. For flood forecasting, emphasize high flows. For low-flow ecology, emphasize baseflow. For general-purpose simulation, use balanced error models.

---

#### Log-Likelihood Expression

For computational efficiency and numerical stability, we work with the log-likelihood rather than the likelihood itself.

The likelihood for a single observation is:

$$
p(y_i \mid \theta, x_i) = \frac{1}{\sqrt{2\pi\sigma^2}} \exp\left( -\frac{(y_i - \beta_0 - \beta_1 x_i)^2}{2\sigma^2} \right)
$$

For $n$ independent observations, the joint likelihood is the product:

$$
p(y \mid \theta) = \prod_{i=1}^{n} p(y_i \mid \theta, x_i)
$$

Taking logarithms (which is a monotonic transformation, preserving maxima):

$$
\log p(y \mid \theta) = \sum_{i=1}^{n} \log p(y_i \mid \theta, x_i)
$$

Expanding the Gaussian density:

$$
\log p(y \mid \theta) = \sum_{i=1}^{n} \left[ -\frac{1}{2}\log(2\pi) - \frac{1}{2}\log(\sigma^2) - \frac{(y_i - \beta_0 - \beta_1 x_i)^2}{2\sigma^2} \right]
$$

Simplifying:

$$
\log p(y \mid \theta) = -\frac{n}{2}\log(2\pi\sigma^2) - \frac{1}{2\sigma^2} \sum_{i=1}^{n}(y_i - \beta_0 - \beta_1 x_i)^2
$$

This expression reveals several important features:

**First term** ($-\frac{n}{2}\log(2\pi\sigma^2)$): Normalizing constant ensuring the Gaussian density integrates to 1. This term penalizes large $\sigma$—models that claim high uncertainty receive lower likelihood, all else equal. This creates tension: we want small $\sigma$ to maximize likelihood, but if $\sigma$ is too small, the second term (residual sum of squares) will dominate negatively.

**Second term** ($-\frac{1}{2\sigma^2} \sum_{i=1}^{n}(y_i - \beta_0 - \beta_1 x_i)^2$): Sum of squared errors, weighted by $1/\sigma^2$. This is minimized when predictions $\beta_0 + \beta_1 x_i$ closely match observations $y_i$. The weighting by $1/\sigma^2$ means:
- Small $\sigma$: Large penalty for residuals (model claims high precision, so errors are costly)
- Large $\sigma$: Small penalty for residuals (model admits high uncertainty, so errors are expected)

**Relationship to least squares**: When $\sigma$ is fixed, maximizing log-likelihood is equivalent to minimizing sum of squared errors. Classical least squares regression is thus the maximum likelihood estimate under Gaussian errors. Bayesian calibration generalizes this by:
1. Treating $\sigma$ as unknown (estimating uncertainty magnitude)
2. Adding prior information
3. Producing full posterior distributions rather than point estimates

**Numerical advantages of log-likelihood**:
- Products become sums (computationally simpler, more stable)
- Avoids underflow: $p(y|\theta)$ might be $10^{-300}$, but $\log p(y|\theta) = -690$ is manageable
- Derivatives are simpler for gradient-based optimization
- MCMC acceptance ratios involve likelihood ratios, which become differences in log-space

**Connection to information theory**: The log-likelihood relates to Kullback-Leibler divergence, mutual information, and entropy. Maximizing log-likelihood minimizes the divergence between the model's predictive distribution and the data-generating distribution.

---

### Code: Likelihood Function

```r
# Log-likelihood function
# Returns the log-probability of data given parameters
log_likelihood <- function(theta, x, y) {
  beta0 <- theta[1]  # Extract intercept
  beta1 <- theta[2]  # Extract slope
  sigma <- theta[3]  # Extract residual standard deviation
  
  # Compute model predictions
  # These are the mean values E[y|x, theta]
  y_pred <- beta0 + beta1 * x
  
  # Compute residuals
  residuals <- y - y_pred
  
  # Compute log-likelihood
  # dnorm(y, y_pred, sigma, log=TRUE) computes log of Gaussian density
  # Sum over all observations because of independence assumption
  log_lik <- sum(dnorm(y, mean = y_pred, sd = sigma, log = TRUE))
  
  return(log_lik)
}

# Alternative explicit form showing all terms:
log_likelihood_explicit <- function(theta, x, y) {
  beta0 <- theta[1]
  beta1 <- theta[2]
  sigma <- theta[3]
  
  n <- length(y)
  y_pred <- beta0 + beta1 * x
  residuals <- y - y_pred
  
  # First term: normalizing constant
  term1 <- -(n/2) * log(2 * pi * sigma^2)
  
  # Second term: weighted sum of squared errors
  term2 <- -sum(residuals^2) / (2 * sigma^2)
  
  log_lik <- term1 + term2
  
  return(log_lik)
}
```

**Code explanation**:

**Efficient implementation**: The function `dnorm(y, mean = y_pred, sd = sigma, log = TRUE)` is the standard R function for computing Gaussian log-densities. Using built-in functions is:
- More numerically stable (R's implementation handles edge cases)
- More computationally efficient (optimized compiled code)
- Less error-prone (well-tested)
- More readable (clear statistical intent)

**Vectorization**: The expression `dnorm(y, y_pred, sigma, log=TRUE)` is vectorized—it computes the log-density for all $n$ observations simultaneously. This is much faster than a loop in R.

**Input validation**: In production code, we should add checks:
```r
if(sigma <= 0) return(-Inf)  # Invalid variance
if(length(y) != length(x)) stop("x and y must have same length")
if(any(is.na(y)) || any(is.na(x))) stop("Data contain missing values")
```

**Computational cost**: For $n$ observations, this likelihood evaluation requires:
- $n$ multiplications and additions (computing predictions)
- $n$ subtractions (computing residuals)
- $n$ Gaussian density evaluations
- $n$ summation operations

Total: $O(n)$ complexity. This is very fast. Complex hydrologic models may require minutes to hours per likelihood evaluation because each evaluation requires running a simulation. Here, we can evaluate thousands of likelihoods per second.

**Understanding through visualization**: To build intuition, we can evaluate the likelihood over a grid of parameter values:

```r
# Create parameter grid
beta0_grid <- seq(0, 5, length.out = 50)
beta1_grid <- seq(0.5, 3, length.out = 50)
sigma_fixed <- 1.2

# Evaluate log-likelihood over grid
log_lik_surface <- outer(beta0_grid, beta1_grid, function(b0, b1) {
  sapply(1:length(b0), function(i) {
    log_likelihood(c(b0[i], b1[i], sigma_fixed), x, y)
  })
})

# Visualize likelihood surface
contour(beta0_grid, beta1_grid, log_lik_surface,
        nlevels = 20,
        xlab = "beta0", ylab = "beta1",
        main = "Log-Likelihood Surface")
points(beta0_true, beta1_true, pch = 19, col = "red", cex = 2)
```

This visualization reveals:
- The likelihood surface has a unique maximum (convex for linear regression)
- The maximum is near the true parameter values (as expected)
- The curvature indicates parameter uncertainty (flat = uncertain, peaked = certain)
- Parameters may be correlated (elongated contours indicate correlation)

---

#### Connecting to Maximum Likelihood Estimation

In classical (frequentist) statistics, parameters are estimated by maximizing the likelihood:

$$
\hat{\theta}_{\text{MLE}} = \arg\max_\theta \log p(y \mid \theta)
$$

This produces point estimates. In our linear regression example, the MLE solutions have closed forms:

$$
\hat{\beta}_1 = \frac{\sum (x_i - \bar{x})(y_i - \bar{y})}{\sum (x_i - \bar{x})^2}
$$

$$
\hat{\beta}_0 = \bar{y} - \hat{\beta}_1 \bar{x}
$$

$$
\hat{\sigma}^2 = \frac{1}{n}\sum (y_i - \hat{\beta}_0 - \hat{\beta}_1 x_i)^2
$$

These are the familiar least-squares regression formulas.

**Bayesian vs. Frequentist**: Bayesian inference differs by:
1. Adding prior information: $p(\theta|y) \propto p(y|\theta) p(\theta)$ instead of just $p(y|\theta)$
2. Producing distributions: Full posterior $p(\theta|y)$ instead of point estimate $\hat{\theta}$
3. Propagating uncertainty: Predictions integrate over posterior instead of conditioning on $\hat{\theta}$

However, when priors are weak (non-informative) and data are plentiful, the posterior mean often closely approximates the MLE. The likelihood function is thus a shared foundation between Bayesian and frequentist approaches.

---

### 4.3.3 Posterior Distribution

Bayes' theorem combines prior and likelihood to produce the posterior distribution:

$$
p(\theta \mid y) \propto p(y \mid \theta) \, p(\theta)
$$

The posterior distribution represents updated parameter knowledge after observing data. It is the central object of inference in Bayesian statistics and the primary output of Bayesian calibration.

**Reading Bayes' theorem**:

The equation states: "The probability of parameters given data is proportional to the probability of data given parameters, times the probability of parameters."

More intuitively:
- **Left side** $p(\theta|y)$: What we want to learn (parameter knowledge after seeing data)
- **Right side** $p(y|\theta)$: How well parameters explain data (likelihood)
- **Right side** $p(\theta)$: What we knew before seeing data (prior)
- **Proportionality** $\propto$: We ignore the normalizing constant $p(y) = \int p(y|\theta)p(\theta)d\theta$ because MCMC doesn't require it

**The Bayesian update process**:

1. **Start with prior**: Before seeing data, our knowledge is described by $p(\theta)$
2. **Observe data**: We collect observations $y$
3. **Evaluate likelihood**: We compute how probable the data are under different parameter values: $p(y|\theta)$
4. **Combine multiplicatively**: The posterior is proportional to prior times likelihood
5. **Result**: Parameter values that are both a priori plausible AND consistent with data receive high posterior probability

**Why proportionality is sufficient**:

The full Bayes' theorem includes a normalizing constant:

$$
p(\theta \mid y) = \frac{p(y \mid \theta) \, p(\theta)}{p(y)}
$$

where:

$$
p(y) = \int p(y \mid \theta) \, p(\theta) \, d\theta
$$

This integral, called the **marginal likelihood** or **evidence**, is often intractable (impossible to compute analytically). It requires integrating the product of likelihood and prior over the entire parameter space—a high-dimensional integral without closed form.

However, MCMC sampling does not require knowing $p(y)$ because:
- We only need to evaluate $p(\theta|y)$ up to a proportionality constant
- MCMC acceptance ratios involve likelihood ratios, where $p(y)$ cancels out
- The Markov chain still converges to the correct posterior distribution

This is one of the key computational advantages of MCMC: we can sample from complex posteriors without computing intractable normalizing integrals.

---

#### Why Posterior Distributions Matter

Unlike classical point estimates (e.g., maximum likelihood or least squares), the Bayesian posterior provides a complete probabilistic description of parameter uncertainty:

**1. Parameter uncertainty quantification**:

The posterior distribution directly tells us how uncertain we are about each parameter:
- **Narrow posterior**: Parameter is well-constrained by data
- **Wide posterior**: Parameter is poorly identifiable
- **Multi-modal posterior**: Multiple parameter values explain data equally well (equifinality)

We can compute credible intervals: "There is 95% posterior probability that $\beta_1 \in [1.6, 2.0]$" is a probabilistic statement about our uncertainty, not a frequency statement about repeated sampling.

**2. Parameter correlations**:

The joint posterior $p(\beta_0, \beta_1, \sigma | y)$ reveals dependencies between parameters:
- If $\beta_0$ and $\beta_1$ are negatively correlated (common in regression), increasing the intercept requires decreasing the slope to maintain fit
- Correlations reveal compensation effects and identifiability issues
- Marginal posteriors $p(\beta_0|y)$ integrate out (marginalize over) other parameters

**3. Propagation of uncertainty into predictions**:

Instead of making predictions with a single "best" parameter set, we integrate over the posterior:

$$
p(y_{\text{new}} \mid x_{\text{new}}, y) = \int p(y_{\text{new}} \mid x_{\text{new}}, \theta) \, p(\theta \mid y) \, d\theta
$$

This posterior predictive distribution accounts for:
- Parameter uncertainty (different $\theta$ produce different predictions)
- Residual uncertainty (even given $\theta$, predictions are uncertain)

The result is prediction intervals that have correct frequentist coverage: if we make many predictions with 95% intervals, approximately 95% will contain the true value.

**4. Decision-making under uncertainty**:

Posterior distributions enable expected utility calculations and risk-based decisions:
- What is the probability that flow exceeds a threshold?
- What is the expected flood damage given parameter uncertainty?
- Which design alternative minimizes expected cost?

These require integrating over the posterior, which is straightforward with MCMC samples.

**5. Model comparison**:

The marginal likelihood $p(y)$ can be approximated from MCMC output and used for Bayesian model comparison:
- Bayes factors: $\frac{p(y|M_1)}{p(y|M_2)}$ quantify relative support for models
- Information criteria (DIC, WAIC) approximate model evidence
- Cross-validation predictive performance can be computed from posterior samples

---

#### Posterior as a Compromise Between Prior and Likelihood

The posterior represents a Bayesian compromise between prior beliefs and data evidence. Several scenarios illustrate this:

**Scenario 1: Weak prior, informative data**
- Prior is broad and flat
- Data provide strong constraints (large $n$, small $\sigma$)
- Result: Posterior is dominated by likelihood, closely approximates MLE

**Scenario 2: Strong prior, weak data**
- Prior is narrow and peaked
- Data are noisy or sparse (small $n$, large $\sigma$)
- Result: Posterior remains close to prior, data have limited influence

**Scenario 3: Conflicting prior and data**
- Prior centers on one parameter region
- Likelihood centers on a different region
- Result: Posterior is a weighted compromise, potentially bimodal

**Scenario 4: Aligned prior and data**
- Prior and likelihood center on similar parameter values
- Result: Posterior is narrower than either (information is combined)

In hydrologic calibration:
- Early in a study with limited data, priors (from literature or expert judgment) strongly influence posteriors
- As more data accumulate, posteriors become increasingly data-driven
- For well-identified parameters, priors matter little; for poorly identified parameters, priors remain influential

---

### Code: Posterior Function

```r
# Log-posterior density function
# Returns the log-probability of parameters given data
# This is the target distribution for MCMC sampling
log_posterior <- function(theta, x, y) {
  
  # Compute log-prior
  # If prior is zero (e.g., sigma < 0), return -Inf immediately
  lp <- log_prior(theta)
  if(is.infinite(lp)) return(lp)
  
  # Compute log-likelihood
  ll <- log_likelihood(theta, x, y)
  
  # Return log-posterior = log-prior + log-likelihood
  # This is equivalent to prior * likelihood in original scale
  # We don't need the normalizing constant p(y)
  return(lp + ll)
}
```

**Code explanation**:

**Computational efficiency**: The function first evaluates the prior. If the prior is zero (returns `-Inf` in log-space), we immediately return `-Inf` without computing the likelihood. This saves computation for parameter values that are impossible a priori (e.g., negative variance).

**Why check for infinite log-prior?**: 
```r
if(is.infinite(lp)) return(lp)
```
This line is crucial. If we tried to compute `lp + ll` when `lp = -Inf`, we'd get:
- If `ll` is finite: `-Inf + ll = -Inf` ✓ (correct)
- If `ll = Inf`: `-Inf + Inf = NaN` ✗ (undefined, breaks MCMC)

By returning early, we avoid this edge case and save computation.

**Addition in log-space**: 
```r
lp + ll
```
Since $\log(ab) = \log(a) + \log(b)$, adding log-prior and log-likelihood gives us the log-posterior up to the normalizing constant:

$$
\log p(\theta|y) = \log p(y|\theta) + \log p(\theta) - \log p(y)
$$

We omit $\log p(y)$ because it's constant with respect to $\theta$ (doesn't depend on which parameters we're evaluating), so it doesn't affect:
- MCMC acceptance probabilities (they involve ratios)
- Optimization for finding posterior modes
- Relative comparisons between parameter values

**Modularity**: By separating `log_prior` and `log_likelihood` into distinct functions, we gain:
- Code clarity (each component is explicit)
- Reusability (can test each component separately)
- Flexibility (can swap different priors or likelihoods easily)
- Debugging (can check which component causes issues)

**Testing the posterior**: Before running MCMC, we should verify:

```r
# Test at true parameter values
theta_true <- c(beta0_true, beta1_true, sigma_true)
log_post_true <- log_posterior(theta_true, x, y)
print(paste("Log-posterior at true values:", log_post_true))

# Test at random values
theta_random <- c(rnorm(1, 0, 10), rnorm(1, 0, 10), abs(rnorm(1, 0, 5)))
log_post_random <- log_posterior(theta_random, x, y)
print(paste("Log-posterior at random values:", log_post_random))

# The true values should generally have higher log-posterior
# (though not guaranteed with finite data)
```

**Visualizing the posterior** (before MCMC sampling):

For low-dimensional problems, we can evaluate the posterior over a grid and visualize it directly:

```r
# Create parameter grid (fixing sigma at true value for 2D viz)
beta0_seq <- seq(0, 5, length.out = 100)
beta1_seq <- seq(0.5, 3, length.out = 100)

# Compute log-posterior over grid
log_post_grid <- outer(beta0_seq, beta1_seq, function(b0, b1) {
  sapply(1:length(b0), function(i) {
    log_posterior(c(b0[i], b1[i], sigma_true), x, y)
  })
})

# Convert to posterior density (exponentiate and normalize)
post_grid <- exp(log_post_grid - max(log_post_grid))  # Normalize for stability
post_grid <- post_grid / sum(post_grid)  # Normalize to sum to 1

# Visualize
contour(beta0_seq, beta1_seq, post_grid,
        nlevels = 20,
        xlab = "beta0", ylab = "beta1",
        main = "Posterior Distribution (conditional on sigma)")
points(beta0_true, beta1_true, pch = 19, col = "red", cex = 2)
```

This shows the posterior before MCMC—a useful sanity check to ensure our model specification makes sense.

---

#### The Posterior as the Foundation for All Inference

Once we have posterior samples (from MCMC), all downstream inference derives from them:

**Point estimates**:
- Posterior mean: $\mathbb{E}[\theta|y] = \int \theta \, p(\theta|y) \, d\theta \approx \frac{1}{N}\sum_{i=1}^N \theta^{(i)}$
- Posterior median: 50th percentile of marginal posterior
- Posterior mode (MAP): $\arg\max_\theta p(\theta|y)$

**Uncertainty quantification**:
- Posterior standard deviation: $\sqrt{\text{Var}[\theta|y]}$
- Credible intervals: Quantiles of marginal posteriors
- Highest posterior density (HPD) intervals: Shortest intervals containing desired probability mass

**Predictions**:
- Posterior predictive mean: $\mathbb{E}[y_{\text{new}}|x_{\text{new}}, y]$
- Posterior predictive intervals: Quantiles of $p(y_{\text{new}}|x_{\text{new}}, y)$

**Model checking**:
- Posterior predictive p-values
- Residual analysis using posterior mean parameters
- Prior-posterior comparison to check data informativeness

All of these are straightforward to compute from MCMC samples, which is why obtaining good posterior samples is the primary goal of Bayesian calibration.

---

## 4.4 MCMC Sampling

The posterior distribution does not have a closed-form solution for most realistic models, so we approximate it using Markov Chain Monte Carlo (MCMC) sampling.

**The computational challenge**:

For our linear regression model, the posterior is:

$$
p(\theta \mid y) \propto p(y \mid \theta) \, p(\theta)
$$

Even though we can evaluate this function at any point $\theta$, we cannot:
1. Compute the normalizing constant $p(y) = \int p(y|\theta)p(\theta)d\theta$ analytically
2. Find the marginal distributions $p(\beta_0|y)$, $p(\beta_1|y)$, $p(\sigma|y)$ in closed form
3. Compute expectations like $\mathbb{E}[\beta_1|y]$ or $\text{Var}[\beta_1|y]$ directly

The integral is three-dimensional (one for each parameter) and has no analytical solution. Even numerical integration becomes infeasible in higher dimensions—hydrologic models may have 10-50 parameters, making gridding impossible (curse of dimensionality).

**MCMC solution**:

Instead of computing the posterior analytically, MCMC generates a sequence of samples:

$$
\theta^{(1)}, \theta^{(2)}, \theta^{(3)}, \ldots, \theta^{(N)}
$$

that, for large $N$, behave **as if** they were random draws from the true posterior $p(\theta|y)$. We can then approximate any posterior quantity using these samples:

$$
\mathbb{E}[f(\theta)|y] = \int f(\theta) \, p(\theta|y) \, d\theta \approx \frac{1}{N} \sum_{i=1}^{N} f(\theta^{(i)})
$$

For example:
- Posterior mean: $f(\theta) = \theta$
- Posterior variance: $f(\theta) = (\theta - \mathbb{E}[\theta|y])^2$
- Probability of exceedance: $f(\theta) = \mathbb{I}(\theta > c)$
- Predictive distribution: $f(\theta) = p(y_{\text{new}}|x_{\text{new}}, \theta)$

---

### Why MCMC is Necessary

Bayesian calibration often produces complex posterior shapes that require sophisticated sampling methods:

**1. High dimensionality**: Hydrologic models often have many parameters ($d = 10$ to $d = 50$). Direct numerical integration over a $d$-dimensional space is computationally prohibitive. MCMC scales efficiently to high dimensions.

**2. Complex geometry**: Posterior distributions may be:
- **Multimodal**: Multiple parameter combinations fit data equally well (equifinality)
- **Correlated**: Parameters trade off against each other
- **Constrained**: Parameters must satisfy bounds or relationships
- **Heavy-tailed**: Rare parameter combinations still have non-negligible probability

MCMC can explore these complex shapes adaptively.

**3. Non-standard distributions**: Unlike simple distributions (normal, gamma, etc.), posteriors in Bayesian calibration rarely have standard forms. MCMC requires only the ability to **evaluate** the posterior (up to a constant), not to sample from it directly.

**4. Computational efficiency**: For expensive models (like hydrologic simulators that take minutes per run), MCMC focuses samples in high-probability regions rather than wasting evaluations exploring low-probability areas.

---

### The MCMC Algorithm

MCMC constructs a **Markov chain**: a sequence where each state depends only on the previous state:

$$
\theta^{(t+1)} \sim K(\theta^{(t+1)} | \theta^{(t)})
$$

The transition kernel $K$ is designed such that the chain's stationary distribution is the target posterior $p(\theta|y)$.

**Key MCMC properties**:

1. **Markov property**: The next state depends only on the current state, not the full history
2. **Convergence**: For large $t$, $\theta^{(t)}$ is approximately distributed according to $p(\theta|y)$
3. **Ergodicity**: The chain explores all regions of parameter space with non-zero posterior probability
4. **Irreducibility**: Any parameter value can be reached from any other value in finite steps

The most common MCMC algorithm is **Metropolis-Hastings**:

```
1. Initialize: Set theta^(0) to some starting value
2. For t = 1, 2, ..., N:
   a. Propose: theta_prop ~ q(theta | theta^(t-1))
   b. Compute acceptance ratio:
      alpha = min(1, [p(theta_prop|y) * q(theta^(t-1)|theta_prop)] / 
                      [p(theta^(t-1)|y) * q(theta_prop|theta^(t-1))])
   c. Accept or reject:
      u ~ Uniform(0,1)
      if u < alpha:
          theta^(t) = theta_prop  (accept)
      else:
          theta^(t) = theta^(t-1)  (reject, stay in place)
```

**Intuition**: 
- Propose a new parameter value based on the current value
- Accept if it has higher posterior probability
- Sometimes accept even if posterior is lower (to explore the space)
- Rejected proposals mean staying in the current state
- Over time, the chain spends more iterations in high-probability regions

**Symmetric proposals**: When $q$ is symmetric ($q(\theta'|\theta) = q(\theta|\theta')$), the acceptance ratio simplifies to:

$$
\alpha = \min\left(1, \frac{p(\theta_{\text{prop}}|y)}{p(\theta^{(t-1)}|y)}\right)
$$

This is why we only need to evaluate the posterior **up to a constant**—the normalizing constant $p(y)$ cancels in the ratio.

---

### Adaptive MCMC: The `adaptMCMC` Package

Standard Metropolis-Hastings requires careful tuning of the proposal distribution $q$. Poor tuning leads to:
- **Too-small proposals**: Chain moves slowly, high acceptance but poor exploration
- **Too-large proposals**: Most proposals rejected, chain gets stuck
- **Uncorrelated proposals**: Doesn't account for parameter correlations

The optimal acceptance rate for random-walk Metropolis-Hastings is approximately **0.234** (23.4%) in high dimensions, a result from theoretical analysis.

**Adaptive MCMC** automatically tunes the proposal distribution during sampling by:
1. Monitoring acceptance rates
2. Adjusting proposal scale to achieve target acceptance rate
3. Estimating parameter correlations from recent samples
4. Rotating proposal distribution to align with posterior contours

The `adaptMCMC` package implements this using the Adaptive Metropolis algorithm (Haario et al. 2001), which:
- Uses a multivariate normal proposal
- Adapts the covariance matrix based on sample history
- Ensures ergodicity through mixing with a fixed proposal

This relieves the user from manual tuning and significantly improves efficiency, especially for correlated parameters.

---

### Code: MCMC Sampling

```r
# Load required packages
library(adaptMCMC)  # For adaptive MCMC sampling
library(coda)       # For MCMC diagnostics and visualization

# Set initial parameter values
# These should be in a reasonable region (ideally, close to posterior mode)
# We start at the origin for regression parameters and sigma = 1
theta_init <- c(0, 0, 1)

# Run adaptive MCMC
# This automatically tunes the proposal distribution
mcmc_result <- MCMC(
  p = log_posterior,        # Target distribution (log-posterior function)
  init = theta_init,        # Initial parameter values
  n = 10000,                # Number of MCMC iterations
  adapt = TRUE,             # Enable adaptive tuning
  acc.rate = 0.234,         # Target acceptance rate (optimal for random walk)
  scale = c(0.5, 0.5, 0.2), # Initial proposal scales (will be adapted)
  x = x,                    # Additional arguments passed to log_posterior
  y = y
)

# Discard burn-in period
# Early samples are not from the stationary distribution
# We discard the first 2000 iterations (20% of total)
burn_in <- 2000
samples <- mcmc_result$samples[(burn_in + 1):nrow(mcmc_result$samples), ]

# Assign parameter names for clarity
colnames(samples) <- c("beta0", "beta1", "sigma")

# Convert to 'coda' mcmc object for diagnostics
# This enables traceplot(), densplot(), etc.
samples_mcmc <- as.mcmc(samples)

# Print summary
print(paste("Total iterations:", nrow(mcmc_result$samples)))
print(paste("Burn-in discarded:", burn_in))
print(paste("Retained samples:", nrow(samples)))
print(paste("Acceptance rate:", round(mcmc_result$acceptance.rate, 3)))
```

**Code explanation**:

**MCMC function arguments**:

- `p = log_posterior`: The function to sample from. Must return **log** probability. The `adaptMCMC` package evaluates this function at each proposed parameter value.

- `init = theta_init`: Starting values for the chain. Ideally, these should be:
  - In a region of non-zero posterior probability (otherwise, chain may get stuck)
  - Reasonably close to the posterior mode (reduces burn-in)
  - Can be found by quick optimization or using prior mean
  
  Poor initialization can lead to very long burn-in periods or even failure to converge.

- `n = 10000`: Number of MCMC iterations. More iterations give better approximation of the posterior, but take longer. Rules of thumb:
  - For simple problems (2-3 parameters): 5,000-10,000 iterations often sufficient
  - For complex problems (10+ parameters): 50,000-500,000 may be needed
  - Run multiple chains to assess convergence

- `adapt = TRUE`: Enables adaptive tuning. The algorithm adjusts the proposal covariance matrix during sampling to improve efficiency. Adaptation typically stops after a "learning phase" to ensure theoretical convergence guarantees.

- `acc.rate = 0.234`: Target acceptance rate. The algorithm adjusts proposal scales to achieve this rate. 
  - 0.234 is theoretically optimal for high-dimensional Gaussian targets
  - For low dimensions (1-2 parameters), slightly higher rates (0.3-0.5) may be better
  - If your empirical acceptance rate is very different, the chain may not be mixing well

- `scale = c(0.5, 0.5, 0.2)`: Initial proposal standard deviations for each parameter. These are rough guesses that adaptation will refine. Guidelines:
  - Start with approximately 10-20% of the parameter's expected range
  - Smaller for well-identified parameters, larger for uncertain parameters
  - Different scales for different parameters account for different magnitudes

- `x = x, y = y`: Additional arguments passed to `log_posterior()`. This is a convenient feature—we don't need to hard-code data into the function.

**Burn-in**:

```r
burn_in <- 2000
samples <- mcmc_result$samples[(burn_in + 1):nrow(mcmc_result$samples), ]
```

The **burn-in** or **warm-up** period consists of initial iterations before the chain reaches its stationary distribution. During burn-in:
- The chain "forgets" its initial value
- Adaptive tuning is active (adjusting proposal)
- Samples do not accurately represent the posterior

We discard burn-in samples because they bias posterior estimates. How much to discard?
- Visually inspect trace plots—discard until chain appears stationary
- Conservative rule: discard 20-50% of total iterations
- Run diagnostics (e.g., Gelman-Rubin) to verify convergence
- Multiple chains with different initializations help assess when burn-in ends

**Why 2000 iterations here?** With 10,000 total, discarding 2000 leaves 8000 samples—likely more than needed for this simple problem, but ensures safety.

**Acceptance rate**:

```r
mcmc_result$acceptance.rate
```

The acceptance rate tells us what fraction of proposals were accepted:
- **Too low (<10%)**: Proposals too aggressive, chain stuck, poor mixing
- **Optimal (~23%)**: Good balance of exploration and acceptance
- **Too high (>50%)**: Proposals too timid, slow exploration, high autocorrelation

The `adaptMCMC` algorithm automatically adjusts to achieve the target rate, so if it converges to ~23%, adaptation worked well.

**Converting to coda format**:

```r
samples_mcmc <- as.mcmc(samples)
```

The `coda` package provides extensive MCMC diagnostic tools. Converting our samples to a `coda` object enables:
- `traceplot()`: Visualize chain paths
- `densplot()`: Visualize marginal posteriors
- `autocorr.plot()`: Check sample autocorrelation
- `effectiveSize()`: Estimate effective sample size
- `gelman.diag()`: Gelman-Rubin convergence diagnostic (requires multiple chains)

---

### Understanding MCMC Output

The MCMC output is a $8000 \times 3$ matrix where:
- Each row is one parameter sample: $\theta^{(i)} = (\beta_0^{(i)}, \beta_1^{(i)}, \sigma^{(i)})$
- Rows are sequentially dependent (autocorrelated)
- Collectively, they approximate the posterior distribution

We can now compute any posterior quantity:

```r
# Posterior means
post_mean <- colMeans(samples)
print(post_mean)

# Compare to true values
true_vals <- c(beta0_true, beta1_true, sigma_true)
print(rbind(true_vals, post_mean))

# Posterior standard deviations
post_sd <- apply(samples, 2, sd)
print(post_sd)

# 95% credible intervals
post_ci <- apply(samples, 2, quantile, probs = c(0.025, 0.975))
print(post_ci)

# Probability that beta1 > 1.5
prob_beta1_gt_1.5 <- mean(samples[, "beta1"] > 1.5)
print(paste("P(beta1 > 1.5 | y) =", prob_beta1_gt_1.5))
```

These are all straightforward sample statistics—MCMC transforms a complex integration problem into simple Monte Carlo averaging.

---

### Computational Considerations

**Runtime**: For this simple problem, 10,000 MCMC iterations complete in seconds. Each iteration requires:
1. One `log_posterior` evaluation (~microseconds for linear regression)
2. One uniform random number (~nanoseconds)
3. Simple arithmetic for acceptance decision (~nanoseconds)

**Scaling to complex models**: For hydrologic models where each likelihood evaluation requires running a simulation (seconds to minutes), MCMC becomes expensive:
- 10,000 iterations × 10 seconds/iteration = 28 hours
- 100,000 iterations × 1 minute/iteration = 69 days

This motivates:
- **Emulators**: Approximate the model with a fast surrogate (Gaussian process, neural network)
- **Parallel chains**: Run multiple chains on different processors
- **Efficient algorithms**: Use gradient-based methods (Hamiltonian Monte Carlo) if derivatives available
- **Adaptive stopping**: Monitor convergence and stop when sufficient samples obtained

For our tutorial, the fast linear model allows extensive experimentation without computational burden.

---

## 4.5 Convergence Diagnostics

Posterior samples are only meaningful if the Markov chain has converged to its stationary distribution. MCMC convergence diagnostics are essential quality control tools that assess whether our samples reliably represent the posterior.

**The convergence problem**:

MCMC is an asymptotic method—theoretically, as the number of iterations approaches infinity, the sample distribution approaches the true posterior. In practice, we have finite samples and must answer: "Are our samples good enough?"

Convergence diagnostics cannot **prove** convergence (we never know the true posterior), but they can:
1. **Detect problems**: Identify clear failures (chains stuck, not mixing, not stationary)
2. **Build confidence**: Show evidence consistent with convergence
3. **Guide decisions**: Determine if more iterations are needed

**Multiple complementary diagnostics** should always be used—no single test is definitive.

---

### Diagnostic Concepts

**1. Trace plots** reveal sampling stability

Trace plots show parameter values over MCMC iterations: $\theta^{(t)}$ vs. $t$.

**What to look for**:
- **Good mixing**: Chain wanders freely, no long-term trends, looks like "fuzzy caterpillar"
- **Stationarity**: Mean and variance appear constant over time
- **No drift**: Chain doesn't systematically increase or decrease
- **No stuck periods**: No long sequences with identical values

**Warning signs**:
- **Trend**: Chain drifting up or down (not converged, need more burn-in)
- **Stickiness**: Chain stays at same value for many iterations (poor mixing, adjust proposals)
- **Distinct phases**: Abrupt changes in behavior (multiple modes, poor mixing)
- **Extreme autocorrelation**: Very slow movement (inefficient sampling)

Trace plots are usually the first diagnostic to examine—visual inspection often reveals obvious problems.

---

**2. Density plots** show posterior structure

Density plots estimate the marginal posterior distribution $p(\theta_j|y)$ from samples using kernel density estimation.

**What to look for**:
- **Smooth distribution**: Sufficient samples to estimate density reliably
- **Reasonable shape**: Unimodal, roughly symmetric (or at least continuous)
- **Finite support**: Density concentrates in finite region (not spreading indefinitely)
- **Agreement with prior**: Posterior should differ from prior (data informative)

**Warning signs**:
- **Highly irregular**: Spiky, multimodal without physical justification (insufficient mixing or samples)
- **Identical to prior**: Data did not update beliefs (parameter not identifiable or likelihood problem)
- **Extreme skewness**: May indicate poor parameterization or convergence issues
- **Fat tails**: May indicate occasional excursions to extreme values (chain exploring very different modes)

Compare marginal posteriors to priors to assess how much the data learned about each parameter.

---

**3. Autocorrelation** reveals sample dependence

Autocorrelation measures the correlation between samples separated by $k$ iterations (lag $k$):

$$
\rho_k = \text{Corr}(\theta^{(t)}, \theta^{(t+k)})
$$

**What to look for**:
- **Rapid decay**: Autocorrelation drops quickly as lag increases
- **Near zero after ~10-20 lags**: Samples become effectively independent
- **Exponential decay pattern**: Typical for well-mixing chains

**Warning signs**:
- **Slow decay**: Autocorrelation remains high for many lags (inefficient sampling)
- **Periodic patterns**: Oscillation in autocorrelation (strange mixing behavior)
- **Nearly constant**: Chain almost not moving (critical mixing problem)

High autocorrelation means samples are redundant—effective sample size is much smaller than actual sample size.

---

**4. Effective sample size** measures information content

Raw MCMC samples are autocorrelated, so $N$ samples contain less information than $N$ independent samples. The **effective sample size** (ESS) estimates the equivalent number of independent samples:

$$
\text{ESS} = \frac{N}{1 + 2\sum_{k=1}^{\infty} \rho_k}
$$

where $N$ is the actual sample size and $\rho_k$ is the lag-$k$ autocorrelation.

**Guidelines**:
- **ESS > 1000**: Generally sufficient for stable posterior inference
- **ESS > 100**: Minimum for reasonable estimates, but increase samples if possible
- **ESS < 100**: Too few effective samples, need more iterations or better proposals
- **ESS/N ratio**: Measures efficiency (higher is better, perfect mixing gives ESS = N)

Example: If we have 10,000 samples but ESS = 500, we effectively have only 500 independent draws. This is still useful, but we should be cautious about estimating tail probabilities or rare events.

---

### Code: Convergence Diagnostics

```r
# =============================================================================
# 1. TRACE PLOTS - Visual assessment of mixing and stationarity
# =============================================================================

# Create trace plots for all parameters
# These show parameter evolution over iterations
traceplot(samples_mcmc,
          main = "MCMC Trace Plots",
          smooth = TRUE)  # Add smoothed mean line

# What to look for:
# - Good mixing: chain explores parameter space freely
# - Stationarity: no trends or drifts
# - No burn-in artifacts: if trace starts at unusual value and drifts, 
#   may need longer burn-in

# We can also create custom trace plots with more control:
par(mfrow = c(3, 1))
plot(samples[, "beta0"], type = "l", 
     ylab = "beta0", xlab = "Iteration",
     main = "Trace Plot: beta0")
abline(h = beta0_true, col = "red", lty = 2, lwd = 2)

plot(samples[, "beta1"], type = "l",
     ylab = "beta1", xlab = "Iteration", 
     main = "Trace Plot: beta1")
abline(h = beta1_true, col = "red", lty = 2, lwd = 2)

plot(samples[, "sigma"], type = "l",
     ylab = "sigma", xlab = "Iteration",
     main = "Trace Plot: sigma")
abline(h = sigma_true, col = "red", lty = 2, lwd = 2)

# Red dashed lines show true values
# Chain should fluctuate around true values if properly converged

# =============================================================================
# 2. DENSITY PLOTS - Visualize marginal posterior distributions
# =============================================================================

# Create density plots overlaid with prior distributions
densplot(samples_mcmc,
         main = "Posterior Density Plots")

# Custom density plots with prior comparison:
par(mfrow = c(1, 3))

# Beta0: Posterior vs. Prior
hist(samples[, "beta0"], breaks = 30, probability = TRUE,
     col = "lightblue", border = "white",
     main = "beta0: Posterior vs Prior",
     xlab = "beta0")
lines(density(samples[, "beta0"]), col = "blue", lwd = 2)
curve(dnorm(x, 0, 10), add = TRUE, col = "red", lwd = 2, lty = 2)
abline(v = beta0_true, col = "darkgreen", lwd = 2, lty = 3)
legend("topright", 
       legend = c("Posterior", "Prior", "True value"),
       col = c("blue", "red", "darkgreen"),
       lty = c(1, 2, 3), lwd = 2)

# Beta1: Posterior vs. Prior
hist(samples[, "beta1"], breaks = 30, probability = TRUE,
     col = "lightblue", border = "white",
     main = "beta1: Posterior vs Prior",
     xlab = "beta1")
lines(density(samples[, "beta1"]), col = "blue", lwd = 2)
curve(dnorm(x, 0, 10), add = TRUE, col = "red", lwd = 2, lty = 2)
abline(v = beta1_true, col = "darkgreen", lwd = 2, lty = 3)
legend("topright",
       legend = c("Posterior", "Prior", "True value"),
       col = c("blue", "red", "darkgreen"),
       lty = c(1, 2, 3), lwd = 2)

# Sigma: Posterior vs. Prior
hist(samples[, "sigma"], breaks = 30, probability = TRUE,
     col = "lightblue", border = "white",
     main = "sigma: Posterior vs Prior",
     xlab = "sigma")
lines(density(samples[, "sigma"]), col = "blue", lwd = 2)
# Half-normal prior: N(0,5) truncated to positive
sigma_prior_x <- seq(0, 10, length.out = 100)
sigma_prior_y <- dnorm(sigma_prior_x, 0, 5) * 2  # *2 for half-normal
lines(sigma_prior_x, sigma_prior_y, col = "red", lwd = 2, lty = 2)
abline(v = sigma_true, col = "darkgreen", lwd = 2, lty = 3)
legend("topright",
       legend = c("Posterior", "Prior", "True value"),
       col = c("blue", "red", "darkgreen"),
       lty = c(1, 2, 3), lwd = 2)

# Interpretation:
# - Posterior should be narrower than prior (data are informative)
# - Posterior should be centered near true value (calibration successful)
# - Large prior-posterior difference indicates strong data signal

# =============================================================================
# 3. EFFECTIVE SAMPLE SIZE - Quantify information content
# =============================================================================

# Compute effective sample size for each parameter
ess <- effectiveSize(samples_mcmc)
print("Effective Sample Sizes:")
print(ess)

# Compute efficiency: ESS / actual sample size
efficiency <- ess / nrow(samples)
print("Sampling Efficiency (ESS/N):")
print(efficiency)

# Interpretation:
# - ESS > 1000: Excellent, very reliable posterior estimates
# - ESS 100-1000: Good, sufficient for most purposes
# - ESS < 100: Warning, may need more samples
# - Efficiency close to 1: Excellent mixing (samples nearly independent)
# - Efficiency < 0.1: Poor mixing (high autocorrelation)

# If ESS is too low, consider:
# 1. Running longer chains (more iterations)
# 2. Adjusting proposal distribution (though adaptMCMC does this)
# 3. Thinning less aggressively (keep more samples)
# 4. Using more efficient MCMC algorithm (e.g., Hamiltonian Monte Carlo)

# =============================================================================
# 4. ACCEPTANCE RATE - Check proposal tuning
# =============================================================================

# Print acceptance rate
print(paste("Acceptance Rate:", 
            round(mcmc_result$acceptance.rate, 3)))

# Interpretation:
# - 0.15-0.30: Good for random walk Metropolis
# - 0.234: Theoretically optimal for high dimensions
# - < 0.10: Proposals too aggressive, chain stuck
# - > 0.50: Proposals too conservative, slow exploration

# The adaptMCMC algorithm should automatically achieve ~0.23
# If far from this, may indicate problems with:
# - Initialization (started in very low probability region)
# - Posterior shape (extremely narrow or multimodal)
# - Scale parameters (initial scales very wrong)

# =============================================================================
# 5. AUTOCORRELATION PLOTS - Check sample independence
# =============================================================================

# Plot autocorrelation function for each parameter
autocorr.plot(samples_mcmc,
              auto.layout = TRUE,
              main = "Autocorrelation Plots")

# Interpretation:
# - Lag 0: Always 1 (sample correlated with itself)
# - Lag 1-10: Should decay rapidly
# - Lag > 20: Should be near 0
# - Slow decay indicates high autocorrelation (inefficient sampling)

# Custom autocorrelation plot with confidence bands:
par(mfrow = c(1, 3))
for (param in colnames(samples)) {
  acf(samples[, param], 
      main = paste("ACF:", param),
      lag.max = 50)
  # Blue dashed lines show 95% significance bounds
  # If autocorrelation stays outside these bounds at high lags,
  # samples are still significantly correlated
}

# =============================================================================
# 6. POSTERIOR SUMMARY WITH DIAGNOSTICS
# =============================================================================

# Comprehensive summary including diagnostics
summary(samples_mcmc)

# This shows:
# - Mean, SD, Naive SE, Time-series SE
# - Quantiles (2.5%, 25%, 50%, 75%, 97.5%)
# - Time-series SE accounts for autocorrelation
# - If Time-series SE >> Naive SE, high autocorrelation present

# =============================================================================
# 7. GELMAN-RUBIN DIAGNOSTIC (requires multiple chains)
# =============================================================================

# For single chain, we can split it into two halves and compare
# (Not as robust as truly independent chains)
n_half <- floor(nrow(samples) / 2)
chain1 <- samples[1:n_half, ]
chain2 <- samples[(n_half + 1):(2 * n_half), ]

mcmc_list <- mcmc.list(as.mcmc(chain1), as.mcmc(chain2))

# Gelman-Rubin diagnostic
# Values close to 1.0 indicate convergence
# Values > 1.1 suggest lack of convergence
gelman_diag <- gelman.diag(mcmc_list)
print("Gelman-Rubin Diagnostic (Split-Chain):")
print(gelman_diag)

# Note: Splitting a single chain is not ideal
# Better practice: Run multiple chains with different starting values
# and use gelman.diag() on the full mcmc.list

# =============================================================================
# 8. VISUAL SUMMARY
# =============================================================================

# Create comprehensive diagnostic panel
par(mfrow = c(2, 2))

# 1. Trace plot for beta1
plot(samples[, "beta1"], type = "l", 
     ylab = "beta1", xlab = "Iteration",
     main = "Trace: beta1")
abline(h = mean(samples[, "beta1"]), col = "red", lty = 2)

# 2. Density plot for beta1
plot(density(samples[, "beta1"]), 
     main = "Density: beta1", 
     xlab = "beta1")
abline(v = beta1_true, col = "red", lty = 2)

# 3. Autocorrelation for beta1
acf(samples[, "beta1"], main = "ACF: beta1")

# 4. Running mean to check convergence
running_mean <- cumsum(samples[, "beta1"]) / (1:nrow(samples))
plot(running_mean, type = "l",
     ylab = "Running Mean", xlab = "Iteration",
     main = "Running Mean: beta1")
abline(h = beta1_true, col = "red", lty = 2)
# Running mean should stabilize if chain has converged

```

**Interpretation guide**:

After running these diagnostics, ask:

1. **Do trace plots show stationarity and mixing?**
   - ✓ Yes → Good sign
   - ✗ No → Need longer burn-in or better proposals

2. **Are posteriors different from priors?**
   - ✓ Yes → Data are informative
   - ✗ No → Parameters not identifiable or likelihood problem

3. **Is ESS > 1000 for all parameters?**
   - ✓ Yes → Sufficient samples
   - ✗ No → Run longer or improve mixing

4. **Is acceptance rate between 0.15-0.30?**
   - ✓ Yes → Good proposal tuning
   - ✗ No → Check for initialization or posterior problems

5. **Does autocorrelation decay rapidly?**
   - ✓ Yes → Efficient sampling
   - ✗ No → Thin samples or run longer

If all checks pass, we can confidently use the posterior samples for inference. If any fail, diagnose and fix before proceeding.

---

### What to Do If Diagnostics Fail

**Problem**: Trace plots show trend or drift
**Solution**: Increase burn-in period, check initialization

**Problem**: Very low ESS (< 100)
**Solution**: Run more iterations, adjust proposal scale, consider thinning less

**Problem**: Very low acceptance rate (< 10%)
**Solution**: Check initialization (may be in zero-probability region), reduce proposal scale

**Problem**: Very high acceptance rate (> 50%)
**Solution**: Increase proposal scale (chain moving too timidly)

**Problem**: Multimodal posterior, chain stuck in one mode
**Solution**: Run multiple chains with dispersed initializations, use tempering or parallel tempering, consider reparameterization

**Problem**: Posterior identical to prior
**Solution**: Check likelihood implementation, verify data are loading correctly, assess parameter identifiability

**Problem**: Autocorrelation not decaying
**Solution**: Run longer chain, adjust proposals, consider better MCMC algorithm (HMC, NUTS)

Convergence diagnostics are iterative—run diagnostics, identify problems, adjust, re-run, repeat until satisfactory.

---

## 4.6 Posterior Summaries

Posterior samples allow estimation of central tendency, uncertainty, and any function of parameters. Unlike point estimates (maximum likelihood or least squares), Bayesian posteriors provide complete probabilistic descriptions.

**Why summarize the posterior?**

While the full posterior distribution $p(\theta|y)$ contains all information, we often need concise summaries for:
- **Reporting**: Communicating results in tables or text
- **Decision-making**: Choosing parameter values for predictions or design
- **Comparison**: Contrasting different models or datasets
- **Visualization**: Plotting key statistics without showing full distributions

Different summaries emphasize different aspects of the posterior.

---

### Point Estimates

**Posterior mean**: Expected value under the posterior

$$
\mathbb{E}[\theta|y] = \int \theta \, p(\theta|y) \, d\theta \approx \frac{1}{N} \sum_{i=1}^N \theta^{(i)}
$$

- Minimizes squared error loss
- Affected by skewness and outliers
- Most commonly reported point estimate

**Posterior median**: 50th percentile of marginal posterior

- Minimizes absolute error loss  
- Robust to outliers
- Preferred for skewed distributions

**Posterior mode (MAP)**: Maximum a posteriori estimate

$$
\hat{\theta}_{\text{MAP}} = \arg\max_\theta p(\theta|y)
$$

- Corresponds to peak of posterior density
- Can be found by optimization (faster than MCMC for point estimates only)
- Ignores uncertainty structure
- Unstable for multimodal posteriors

---

### Uncertainty Quantification

**Posterior standard deviation**: Spread of marginal posterior

$$
\text{SD}[\theta|y] = \sqrt{\mathbb{E}[(\theta - \mathbb{E}[\theta|y])^2 | y]}
$$

Measures average deviation from the mean.

**Credible intervals**: Bayesian analog of confidence intervals

$$
P(\theta \in [a, b] | y) = \int_a^b p(\theta|y) \, d\theta = \alpha
$$

For $\alpha = 0.95$, we say "there is 95% posterior probability that $\theta$ lies in $[a,b]$."

**Types of credible intervals**:

1. **Equal-tailed interval**: $(q_{0.025}, q_{0.975})$ where $q_p$ is the $p$-th quantile
   - Easy to compute from samples
   - Standard choice
   - May not be shortest interval for skewed posteriors

2. **Highest Posterior Density (HPD) interval**: Shortest interval containing $\alpha$ probability
   - More informative for skewed distributions
   - Requires density estimation
   - Can be discontinuous for multimodal posteriors

**Interpretation**: A 95% credible interval $[a, b]$ means: "Given the data and prior, I believe with 95% probability that the true parameter value lies in $[a, b]$."

This is a direct probability statement about the parameter—distinct from frequentist confidence intervals, which make probability statements about the procedure, not the parameter.

---

### Code: Posterior Summaries

```r
# =============================================================================
# COMPREHENSIVE POSTERIOR SUMMARY
# =============================================================================

# Use coda's built-in summary function
# This provides mean, SD, quantiles, and Monte Carlo error
summary(samples_mcmc)

# Extract key statistics manually for custom reporting

# Posterior means (point estimates minimizing squared error)
posterior_means <- colMeans(samples)
print("Posterior Means:")
print(posterior_means)

# Posterior medians (point estimates minimizing absolute error)
posterior_medians <- apply(samples, 2, median)
print("Posterior Medians:")
print(posterior_medians)

# Posterior standard deviations (uncertainty measure)
posterior_sds <- apply(samples, 2, sd)
print("Posterior Standard Deviations:")
print(posterior_sds)

# 95% equal-tailed credible intervals
posterior_ci_95 <- apply(samples, 2, quantile, probs = c(0.025, 0.975))
print("95% Credible Intervals:")
print(posterior_ci_95)

# 90% equal-tailed credible intervals
posterior_ci_90 <- apply(samples, 2, quantile, probs = c(0.05, 0.95))
print("90% Credible Intervals:")
print(posterior_ci_90)

# Posterior modes (approximate using maximum density)
# For each parameter, find value with highest density
library(hdrcde)  # For highest density region estimation
posterior_modes <- apply(samples, 2, function(x) {
  dens <- density(x)
  dens$x[which.max(dens$y)]
})
print("Posterior Modes (approximate):")
print(posterior_modes)

# =============================================================================
# COMPARISON TABLE: TRUE VS POSTERIOR ESTIMATES
# =============================================================================

# Create comparison table
comparison <- data.frame(
  Parameter = c("beta0", "beta1", "sigma"),
  True_Value = c(beta0_true, beta1_true, sigma_true),
  Post_Mean = posterior_means,
  Post_Median = posterior_medians,
  Post_SD = posterior_sds,
  CI_2.5 = posterior_ci_95[1, ],
  CI_97.5 = posterior_ci_95[2, ],
  row.names = NULL
)

print("Parameter Estimates vs True Values:")
print(comparison)

# Check if true values fall within credible intervals
# This is a basic posterior validation
in_interval <- (comparison$True_Value >= comparison$CI_2.5) & 
               (comparison$True_Value <= comparison$CI_97.5)
print("True values within 95% CI:")
print(in_interval)

# For successful calibration, we expect:
# - Posterior means close to true values
# - True values within credible intervals (about 95% of the time across many experiments)
# - Reasonable posterior standard deviations (not too narrow = overconfident, not too wide = uninformative)

# =============================================================================
# POSTERIOR PROBABILITIES AND DERIVED QUANTITIES
# =============================================================================

# Probability that beta1 > 1.5 (example hypothesis test)
prob_beta1_gt_1.5 <- mean(samples[, "beta1"] > 1.5)
print(paste("P(beta1 > 1.5 | y) =", round(prob_beta1_gt_1.5, 3)))

# Probability that beta1 is in [1.7, 1.9]
prob_beta1_in_interval <- mean((samples[, "beta1"] >= 1.7) & (samples[, "beta1"] <= 1.9))
print(paste("P(1.7 <= beta1 <= 1.9 | y) =", round(prob_beta1_in_interval, 3)))

# Expected value of a function of parameters
# Example: E[beta0 + beta1 * 5 | y] (prediction at x=5)
pred_at_x5 <- mean(samples[, "beta0"] + samples[, "beta1"] * 5)
print(paste("E[y | x=5, data] =", round(pred_at_x5, 3)))

# Coefficient of variation for sigma
cv_sigma <- sd(samples[, "sigma"]) / mean(samples[, "sigma"])
print(paste("Coefficient of variation for sigma:", round(cv_sigma, 3)))

# Posterior probability that parameters are within 10% of true values
# This is a recovery criterion for synthetic experiments
within_10pct <- colMeans(abs(samples - matrix(rep(c(beta0_true, beta1_true, sigma_true), each = nrow(samples)), ncol = 3, byrow = TRUE)) / matrix(rep(c(beta0_true, beta1_true, sigma_true), each = nrow(samples)), ncol = 3, byrow = TRUE) < 0.10)
print("P(parameter within 10% of true value):")
print(within_10pct)

```

**Interpretation**:

From these summaries, we can assess calibration quality:

1. **Accuracy**: Are posterior means/medians close to true values?
2. **Coverage**: Do 95% intervals contain true values?
3. **Precision**: Are posterior SDs appropriately sized (not too narrow/wide)?
4. **Identifiability**: Are posteriors much narrower than priors?

For this example, we expect:
- $\beta_0$ posterior mean ≈ 2.5
- $\beta_1$ posterior mean ≈ 1.8
- $\sigma$ posterior mean ≈ 1.2
- All true values within 95% credible intervals

If these conditions hold, calibration successfully recovered the data-generating parameters.

---

## 4.7 Joint Posterior Structure

Joint posterior visualization reveals parameter interactions and identifiability beyond what marginal distributions show.

**Why examine joint posteriors?**

Marginal posteriors $p(\theta_j|y)$ show individual parameter uncertainty, but they hide:
- **Correlations**: How parameters covary
- **Trade-offs**: Whether increasing one parameter compensates for decreasing another
- **Identifiability**: Whether parameters can be estimated independently
- **Multimodality**: Whether multiple parameter combinations fit data equally well

These structures critically affect:
- Uncertainty propagation in predictions
- Parameter interpretation
- Model selection and comparison
- Experimental design for reducing uncertainty

---

### Parameter Correlation in Regression

For linear regression, $\beta_0$ and $\beta_1$ are typically **negatively correlated**:

- Increasing intercept + decreasing slope can produce similar fits
- This correlation is stronger when $\bar{x}$ (mean of predictors) is far from 0
- If predictors were centered ($\bar{x} = 0$), correlation would be reduced

This is a fundamental identifiability phenomenon in regression and appears in many hydrologic models (e.g., storage coefficients trading off against initial conditions).

---

### Code: Joint Posterior Plots

```r
# =============================================================================
# PAIRS PLOT: Joint posterior visualization
# =============================================================================

# Create scatterplot matrix showing all pairwise posteriors
pairs(samples, 
      pch = 19,                  # Solid points
      col = rgb(0, 0, 1, 0.1),   # Semi-transparent blue
      main = "Joint Posterior Distributions",
      labels = c("beta0", "beta1", "sigma"))

# Add true values as reference
# We can overlay points or lines showing true values
points(beta0_true, beta1_true, pch = 4, col = "red", cex = 2, lwd = 3)

# Interpretation:
# - Diagonal: Marginal densities (histograms)
# - Off-diagonal: Joint distributions (scatterplots)
# - Elliptical patterns indicate correlation
# - Slope of ellipse shows correlation direction
# - Tightness shows strength of correlation

# =============================================================================
# CUSTOM JOINT PLOTS with contours
# =============================================================================

library(MASS)  # For 2D density estimation

# Function to create nice joint posterior plot
plot_joint_posterior <- function(samples, param1, param2, true1, true2) {
  # Extract samples
  x <- samples[, param1]
  y <- samples[, param2]
  
  # Compute 2D kernel density
  dens <- kde2d(x, y, n = 100)
  
  # Create plot
  par(mfrow = c(1, 1))
  contour(dens, 
          xlab = param1, ylab = param2,
          main = paste("Joint Posterior:", param1, "vs", param2),
          nlevels = 15,
          col = "blue")
  
  # Add points (subsample for clarity)
  points(x[seq(1, length(x), by = 10)], 
         y[seq(1, length(y), by = 10)],
         pch = 19, col = rgb(0, 0, 1, 0.1), cex = 0.5)
  
  # Add true values
  points(true1, true2, pch = 4, col = "red", cex = 2, lwd = 3)
  
  # Add posterior means
  points(mean(x), mean(y), pch = 19, col = "darkgreen", cex = 2)
  
  # Add legend
  legend("topright",
         legend = c("True value", "Posterior mean"),
         pch = c(4, 19),
         col = c("red", "darkgreen"),
         cex = 1.2)
}

# Plot all pairwise combinations
plot_joint_posterior(samples, "beta0", "beta1", beta0_true, beta1_true)
plot_joint_posterior(samples, "beta0", "sigma", beta0_true, sigma_true)
plot_joint_posterior(samples, "beta1", "sigma", beta1_true, sigma_true)

# =============================================================================
# CORRELATION MATRIX
# =============================================================================

# Compute posterior correlation matrix
post_cor <- cor(samples)
print("Posterior Correlation Matrix:")
print(round(post_cor, 3))

# Visualize correlation matrix as heatmap
library(corrplot)
corrplot(post_cor, method = "color", type = "upper",
         addCoef.col = "black", number.cex = 1.2,
         tl.col = "black", tl.srt = 45,
         title = "Posterior Parameter Correlations",
         mar = c(0,0,2,0))

# Interpretation:
# - Diagonal = 1 (perfect self-correlation)
# - Off-diagonal shows parameter correlations
# - Strong correlation (|ρ| > 0.7): Parameters not independently identifiable
# - Weak correlation (|ρ| < 0.3): Parameters largely independent
# - Negative correlation: Trade-off relationships

# For linear regression:
# - beta0 and beta1 typically negatively correlated
# - sigma often independent of beta0, beta1
# - Strength depends on data structure

```

**What to look for in joint posteriors**:

1. **Elliptical patterns**: Indicate bivariate normal-like structure, common for well-behaved posteriors
2. **Diagonal elongation**: Strong positive or negative correlation
3. **Circular patterns**: Parameters independent (low correlation)
4. **Curved or nonlinear patterns**: Nonlinear parameter relationships
5. **Multiple clusters**: Multimodality (equifinality)

**Hydrologic relevance**:

In hydrologic models, parameter correlations reveal:
- **Compensation effects**: Increasing soil storage can compensate for decreasing infiltration capacity
- **Identifiability limits**: Which parameters can be estimated independently
- **Process interactions**: How different model components influence each other
- **Experimental design**: Which data types would best reduce correlations

---

## 4.8 Posterior Predictive Checks

Posterior predictive simulation evaluates whether the calibrated model can reproduce observed behavior and generate realistic predictions.

**The posterior predictive distribution** integrates over parameter uncertainty:

$$
p(y_{\text{new}} | x_{\text{new}}, y) = \int p(y_{\text{new}} | x_{\text{new}}, \theta) \, p(\theta | y) \, d\theta
$$

This distribution includes both:
1. **Parameter uncertainty**: Different $\theta$ values produce different predictions
2. **Residual uncertainty**: Even given $\theta$, predictions are stochastic ($\sigma$ adds noise)

**Why posterior predictive checks matter**:

1. **Model validation**: Do simulations look like real data?
2. **Uncertainty calibration**: Are prediction intervals appropriate?
3. **Assumption checking**: Does residual structure match assumptions?
4. **Forecast evaluation**: For out-of-sample data, how well do predictions perform?

**Posterior predictive p-values**: Compare observed statistics to their predictive distribution to identify discrepancies.

---

### Code: Posterior Predictive Simulation

```r
# =============================================================================
# POSTERIOR PREDICTIVE DISTRIBUTIONS
# =============================================================================

# Create prediction grid
x_pred <- seq(0, 10, length.out = 100)

# Storage for predictive samples
# Each row is one draw from posterior predictive distribution
n_pred_samples <- 500
y_pred_samples <- matrix(NA, nrow = n_pred_samples, ncol = length(x_pred))

# Generate predictive samples by:
# 1. Sample parameter set from posterior
# 2. Generate predictions using those parameters
# 3. Add residual noise
for (i in 1:n_pred_samples) {
  # Randomly select one posterior sample
  idx <- sample(1:nrow(samples), 1)
  b0 <- samples[idx, "beta0"]
  b1 <- samples[idx, "beta1"]
  s  <- samples[idx, "sigma"]
  
  # Generate predictions with noise
  # This is one realization from p(y_new | x_new, y)
  y_pred_samples[i, ] <- rnorm(length(x_pred),
                                mean = b0 + b1 * x_pred,
                                sd = s)
}

# =============================================================================
# VISUALIZE POSTERIOR PREDICTIVE DISTRIBUTION
# =============================================================================

# Compute quantiles of predictive distribution
y_pred_mean <- colMeans(y_pred_samples)
y_pred_median <- apply(y_pred_samples, 2, median)
y_pred_lower <- apply(y_pred_samples, 2, quantile, 0.025)
y_pred_upper <- apply(y_pred_samples, 2, quantile, 0.975)
y_pred_lower_50 <- apply(y_pred_samples, 2, quantile, 0.25)
y_pred_upper_50 <- apply(y_pred_samples, 2, quantile, 0.75)

# Create visualization
plot(x, y, 
     pch = 19, col = "black", cex = 0.8,
     xlab = "Predictor (x)", ylab = "Response (y)",
     main = "Posterior Predictive Distribution",
     ylim = range(c(y, y_pred_lower, y_pred_upper)))

# Add 95% prediction interval
polygon(c(x_pred, rev(x_pred)),
        c(y_pred_lower, rev(y_pred_upper)),
        col = rgb(0, 0, 1, 0.2), border = NA)

# Add 50% prediction interval (darker)
polygon(c(x_pred, rev(x_pred)),
        c(y_pred_lower_50, rev(y_pred_upper_50)),
        col = rgb(0, 0, 1, 0.4), border = NA)

# Add posterior mean prediction
lines(x_pred, y_pred_mean, col = "blue", lwd = 2)

# Add true regression line
abline(a = beta0_true, b = beta1_true, col = "red", lwd = 2, lty = 2)

# Add legend
legend("topleft",
       legend = c("Observed data", "Posterior mean prediction",
                  "True relationship", "50% interval", "95% interval"),
       col = c("black", "blue", "red", rgb(0,0,1,0.4), rgb(0,0,1,0.2)),
       pch = c(19, NA, NA, 15, 15),
       lty = c(NA, 1, 2, NA, NA),
       lwd = c(NA, 2, 2, NA, NA),
       pt.cex = c(0.8, NA, NA, 2, 2))

# =============================================================================
# POSTERIOR PREDICTIVE CHECK: Do observations look typical?
# =============================================================================

# For each observed data point, compute how likely it is
# under the posterior predictive distribution

# This is useful for identifying outliers or model misspecification

# For each observation, generate posterior predictive samples at that x
posterior_predictive_check <- function(x_obs, y_obs) {
  y_pred_at_x <- numeric(nrow(samples))
  
  for (i in 1:nrow(samples)) {
    b0 <- samples[i, "beta0"]
    b1 <- samples[i, "beta1"]
    s  <- samples[i, "sigma"]
    
    # Generate one prediction
    y_pred_at_x[i] <- rnorm(1, b0 + b1 * x_obs, s)
  }
  
  # Compute posterior predictive p-value:
  # P(y_pred >= y_obs | x_obs, data)
  p_value <- mean(y_pred_at_x >= y_obs)
  
  # Two-sided p-value (extreme on either tail)
  p_value_2sided <- 2 * min(p_value, 1 - p_value)
  
  return(list(p_value = p_value, 
              p_value_2sided = p_value_2sided,
              predictive_samples = y_pred_at_x))
}

# Check all observations
ppc_results <- lapply(1:length(x), function(i) {
  posterior_predictive_check(x[i], y[i])
})

# Extract p-values
pvalues_2sided <- sapply(ppc_results, function(r) r$p_value_2sided)

# Identify potential outliers (p < 0.05)
outliers <- which(pvalues_2sided < 0.05)

print(paste("Number of potential outliers (p < 0.05):", length(outliers)))
if (length(outliers) > 0) {
  print("Outlier indices:")
  print(outliers)
}

# In a well-calibrated model with correct assumptions:
# - About 5% of observations should have p < 0.05
# - If many more, model may be mis-specified

# Plot observations colored by p-value
plot(x, y, 
     pch = 19, 
     col = ifelse(pvalues_2sided < 0.05, "red", "blue"),
     cex = ifelse(pvalues_2sided < 0.05, 1.5, 0.8),
     xlab = "Predictor (x)", ylab = "Response (y)",
     main = "Posterior Predictive Check (Red = Potential Outliers)")
abline(a = beta0_true, b = beta1_true, lwd = 2, lty = 2)
legend("topleft", 
       legend = c("Typical (p >= 0.05)", "Outlier (p < 0.05)"),
       col = c("blue", "red"),
       pch = 19,
       cex = c(0.8, 1.5))
```

**Interpretation**:

Good posterior predictive checks show:
- Most observations fall within prediction intervals
- Posterior mean prediction close to true relationship
- About 95% of points within 95% interval
- No systematic patterns in outliers

Poor posterior predictive checks indicate:
- Model structural error (wrong functional form)
- Incorrect residual assumptions
- Unmodeled covariates or processes
- Need for heteroscedastic or autocorrelated error models

---

## 4.9 Residual Analysis

Residual diagnostics test likelihood assumptions and guide model improvement. Even though we calibrated with Bayesian methods, classical residual analysis remains valuable for checking model adequacy.

**Why analyze residuals?**

The likelihood assumed:
1. Normality of errors
2. Homoscedasticity (constant variance)
3. Independence (no autocorrelation)
4. Correct mean function (linear relationship)

If any assumption is violated, posterior inference may be misleading, prediction intervals under-cover, and parameter estimates biased.

Residual plots reveal these violations visually and suggest improvements.

---

### Code: Residual Diagnostics

```r
# =============================================================================
# COMPUTE RESIDUALS using posterior mean parameters
# =============================================================================

# Use posterior mean as point estimate for residual analysis
y_fitted <- posterior_means["beta0"] + posterior_means["beta1"] * x
residuals <- y - y_fitted

# Alternatively, could use posterior median or MAP estimate

# =============================================================================
# DIAGNOSTIC PLOT 1: Residuals vs Fitted Values
# =============================================================================

# Tests homoscedasticity assumption
# If variance is constant, should see random scatter around 0
# Patterns indicate heteroscedasticity

par(mfrow = c(2, 2))

plot(y_fitted, residuals, 
     pch = 19, col = "steelblue",
     xlab = "Fitted Values", ylab = "Residuals",
     main = "Residuals vs Fitted")
abline(h = 0, col = "red", lwd = 2)
smooth_line <- lowess(y_fitted, residuals)
lines(smooth_line, col = "blue", lwd = 2)

# What to look for:
# - Random scatter: Homoscedasticity holds ✓
# - Funnel shape: Variance increases with fitted values (heteroscedasticity) ✗
# - Curved pattern: Nonlinear relationship (wrong mean function) ✗
# - Outliers: Points far from 0 (data anomalies or model failure) ✗

# =============================================================================
# DIAGNOSTIC PLOT 2: Residuals vs Predictor
# =============================================================================

# Tests for patterns related to predictor
# Should see random scatter if model is correct

plot(x, residuals,
     pch = 19, col = "steelblue",
     xlab = "Predictor (x)", ylab = "Residuals",
     main = "Residuals vs Predictor")
abline(h = 0, col = "red", lwd = 2)
smooth_line <- lowess(x, residuals)
lines(smooth_line, col = "blue", lwd = 2)

# What to look for:
# - Random scatter: Linear relationship appropriate ✓
# - Systematic curve: Need nonlinear model (polynomial, spline, etc.) ✗
# - Increasing/decreasing variance: Heteroscedasticity ✗

# =============================================================================
# DIAGNOSTIC PLOT 3: Q-Q Plot (Normality Check)
# =============================================================================

# Tests normality assumption
# If residuals are normal, points should follow straight line

qqnorm(residuals, 
       pch = 19, col = "steelblue",
       main = "Normal Q-Q Plot")
qqline(residuals, col = "red", lwd = 2)

# What to look for:
# - Points on line: Normality holds ✓
# - S-curve: Heavy tails (more extreme values than normal) ✗
# - Reverse S-curve: Light tails ✗
# - Points above line at right: Right skewness ✗
# - Points below line at right: Left skewness ✗

# Formal test (though p-values can be misleading for large samples):
shapiro_test <- shapiro.test(residuals)
print(paste("Shapiro-Wilk normality test p-value:", 
            round(shapiro_test$p.value, 4)))
# p > 0.05 suggests normality cannot be rejected

# =============================================================================
# DIAGNOSTIC PLOT 4: Histogram of Residuals
# =============================================================================

# Visual normality check

hist(residuals, 
     breaks = 20,
     col = "lightblue", border = "white",
     xlab = "Residuals", main = "Histogram of Residuals",
     probability = TRUE)

# Overlay normal curve with estimated parameters
curve(dnorm(x, mean = mean(residuals), sd = sd(residuals)),
      add = TRUE, col = "red", lwd = 2)

# Add vertical line at zero
abline(v = 0, col = "darkgreen", lwd = 2, lty = 2)

# What to look for:
# - Bell-shaped: Approximate normality ✓
# - Skewed: Asymmetric errors ✗
# - Bimodal: Mixture of error types or missing covariate ✗

# =============================================================================
# ADDITIONAL DIAGNOSTICS
# =============================================================================

# Autocorrelation check (for time series data)
# Our data are not time-ordered, but if they were:
acf(residuals, main = "ACF of Residuals")
# Should see no significant autocorrelation (all bars within blue lines)
# If time-ordered data show autocorrelation, need AR error model

# Durbin-Watson test for autocorrelation
library(lmtest)
dw_test <- dwtest(y ~ x)
print("Durbin-Watson test (independence):")
print(dw_test)
# p > 0.05 suggests no autocorrelation

# Breusch-Pagan test for heteroscedasticity
bp_test <- bptest(y ~ x)
print("Breusch-Pagan test (homoscedasticity):")
print(bp_test)
# p > 0.05 suggests homoscedasticity holds

# =============================================================================
# INFLUENCE DIAGNOSTICS
# =============================================================================

# Cook's distance: Identify influential observations
# High Cook's distance means observation strongly influences parameter estimates

lm_fit <- lm(y ~ x)  # Fit for diagnostic purposes
cooks_d <- cooks.distance(lm_fit)

plot(cooks_d, type = "h", 
     ylab = "Cook's Distance", xlab = "Observation Index",
     main = "Cook's Distance (Influence)")
abline(h = 4 / length(x), col = "red", lty = 2)
# Points above red line are potentially influential

influential_points <- which(cooks_d > 4 / length(x))
if (length(influential_points) > 0) {
  print(paste("Influential observations:", 
              paste(influential_points, collapse = ", ")))
}
```

**Interpretation and Next Steps**:

**If diagnostics look good**:
- Random residual scatter
- Normal Q-Q plot follows line
- No autocorrelation
- Homoscedastic variance

→ Model assumptions valid, proceed with confidence in posterior inference

**If heteroscedasticity detected**:
- Funnel pattern in residuals
- Breusch-Pagan test significant

→ Use heteroscedastic error model (Section 3.5): $\sigma_t = a \cdot \hat{Q}_t + b$

**If non-normality detected**:
- Q-Q plot deviates from line
- Histogram skewed

→ Consider log-transformation, robust error models, or Student-t errors

**If autocorrelation detected**:
- ACF shows significant lags
- Durbin-Watson test significant

→ Use AR(1) error model (Section 3.5) or more complex time series models

**If nonlinearity detected**:
- Systematic curve in residuals vs predictor

→ Add polynomial terms, use splines, or switch to nonlinear model

Residual analysis is iterative—diagnose issues, improve model, re-calibrate, re-diagnose until assumptions approximately hold.

---

## 4.10 Summary and Learning Outcomes

This comprehensive example demonstrated that Bayesian calibration:

### 1. Produces Full Parameter Uncertainty Distributions

Unlike point estimates, Bayesian methods yield complete posterior distributions describing parameter uncertainty, correlations, and trade-offs.

**Key learning**: Every parameter estimate comes with quantified uncertainty propagated from data limitations and prior knowledge.

### 2. Integrates Prior Knowledge and Observational Data

Bayesian calibration formally combines prior information (physical constraints, expert judgment, regional data) with observational evidence through Bayes' theorem.

**Key learning**: Prior specification is a modeling choice that should be transparent, justified, and tested through sensitivity analysis and prior-posterior comparison.

### 3. Provides Probabilistic Prediction Envelopes

Posterior predictive distributions integrate both parameter and residual uncertainty, producing prediction intervals with correct coverage properties.

**Key learning**: Predictions account for all sources of uncertainty, enabling risk-based decision-making and honest communication of forecast reliability.

### 4. Requires Diagnostic Validation

MCMC convergence diagnostics, posterior predictive checks, and residual analysis are essential quality control procedures.

**Key learning**: Always validate assumptions, assess convergence, and check model adequacy before trusting posterior inference or predictions.

### 5. Forms the Foundation for Hydrologic Model Calibration

Every component demonstrated here—prior specification, likelihood construction, MCMC sampling, convergence assessment, posterior summarization—transfers directly to complex hydrologic models.

**Key learning**: Linear regression is pedagogically simple but methodologically complete. The same workflow applies to rainfall-runoff models, snow models, groundwater models, and integrated Earth system models.

---

# 5. Bayesian Calibration of GR4J with Homoscedastic Gaussian Errors

This section introduces a **full Bayesian calibration workflow** for the **GR4J hydrological model**, assuming **homoscedastic Gaussian residuals**. The aim is to estimate model parameters while **quantifying uncertainty** in both the hydrological model and observational errors.

## Overview of the GR4J Model

The GR4J (Génie Rural à 4 paramètres Journalier) model is a parsimonious rainfall-runoff model developed by Cemagref (now INRAE) in France. It operates at a daily time step and uses only four core hydrological parameters to transform precipitation and potential evapotranspiration into streamflow. The model's structure consists of:

1. **Production store** (governed by X1): Represents the soil moisture accounting component where rainfall is partitioned between actual evapotranspiration and runoff generation based on current storage levels.

2. **Unit hydrographs**: Two unit hydrographs (UH1 and UH2) controlled by parameter X4 route the generated runoff. UH1 handles 90% of the net rainfall while UH2 routes the remaining 10%, introducing different timing characteristics.

3. **Routing store** (governed by X3): A non-linear reservoir that simulates baseflow and slow groundwater contributions to streamflow.

4. **Groundwater exchange** (governed by X2): Represents potential gains or losses to deep aquifers or neighboring catchments that are not directly observed.

The **extended GR4J** version used here incorporates a **degree-day snow module** with two additional parameters (TT and DDF), making it suitable for catchments with seasonal snow accumulation and melt. This extension is critical for accurately simulating streamflow in mountainous or high-latitude regions where snowmelt dominates the hydrological regime during spring and early summer.

The degree-day approach assumes that snowmelt is proportional to the temperature excess above a threshold (TT). When air temperature falls below TT, precipitation accumulates as snow; when temperature exceeds TT, snow melts at a rate determined by DDF (degree-day factor). This simple yet effective parameterization captures first-order snowmelt dynamics without requiring complex energy balance calculations.

## Why Bayesian Calibration?

Bayesian calibration offers several advantages over traditional optimization-based calibration methods:

- **Uncertainty quantification**: Rather than obtaining single "best" parameter values, we obtain full posterior distributions that reflect parameter uncertainty given the available data and prior knowledge.

- **Incorporation of prior information**: Expert knowledge, regional studies, or physical constraints can be formally encoded as prior distributions, preventing unrealistic parameter combinations.

- **Predictive uncertainty**: Posterior distributions propagate through the model to generate probabilistic streamflow predictions with credible intervals, essential for risk-based water resources management.

- **Diagnosis of identifiability issues**: Posterior correlations and marginal distributions reveal which parameters are well-constrained by data and which remain uncertain, guiding model refinement or data collection priorities.

- **Flexible error modeling**: Bayesian frameworks naturally accommodate various error structures (heteroscedastic, autocorrelated, heavy-tailed) through likelihood specification.

Bayesian calibration allows us to combine **prior knowledge about parameters** with **observed streamflow data** to obtain a **posterior distribution** representing updated knowledge about parameter values through Bayes' theorem:

$$p(\theta \mid Q_{\text{obs}}) \propto p(Q_{\text{obs}} \mid \theta) \cdot p(\theta)$$

where $p(\theta \mid Q_{\text{obs}})$ is the posterior distribution, $p(Q_{\text{obs}} \mid \theta)$ is the likelihood, and $p(\theta)$ is the prior distribution.

**`params = c(X1, X2, X3, X4, TT, DDF)`**

**Units:** precipitation and PET in **mm/day**, temperature in **°C**, area in **km²**, simulated discharge in **m³/s**, residual standard deviation **σ** in **m³/s**.

---

## 5.1 Model Parameters and Bounds

**Purpose.** Define parameters and physically plausible ranges to constrain the posterior, improve sampler efficiency, and encode domain knowledge.

### Parameter list and interpretation

| **Parameter** | **Symbol** | **Meaning** | **Units** |
|---------------|------------|-------------|-----------|
| Production store capacity | **X1** | Controls partitioning between evapotranspiration and runoff generation | mm |
| Groundwater exchange coefficient | **X2** | Exchange between fast and slow components (model scaling) | dimensionless |
| Routing store capacity | **X3** | Controls slow flow magnitude and baseflow | mm |
| Unit hydrograph time base | **X4** | Controls timing/spread of quickflow | days |
| Temperature split | **TT** | Threshold for rain/snow split in degree‑day module | °C |
| Degree‑day factor | **DDF** | Melt rate per °C above TT | mm °C⁻¹ day⁻¹ |
| Residual SD | **σ** | Standard deviation of observation‑model residuals | m³/s |

### Detailed Parameter Interpretation

#### X1: Production Store Capacity (mm)

The production store represents the catchment's soil moisture storage capacity. This parameter fundamentally controls the water balance partitioning:

- **Low X1 values** (< 100 mm): The soil saturates quickly, leading to high runoff coefficients and flashy responses. This is typical of shallow soils, urban areas, or catchments with limited infiltration capacity.

- **High X1 values** (> 500 mm): Deep soils with high water-holding capacity. More precipitation is retained and released as evapotranspiration, reducing runoff generation. Common in forested catchments with deep, permeable soils.

- **Physical interpretation**: X1 approximates the product of soil depth, porosity, and the fraction of catchment actively contributing to runoff generation. However, it is an effective parameter that also absorbs model structural errors and spatial heterogeneity.

- **Seasonal dynamics**: X1 implicitly affects seasonal flow patterns. High X1 catchments show delayed responses to precipitation and more sustained baseflows during dry periods due to soil moisture carryover.

#### X2: Groundwater Exchange Coefficient (dimensionless)

X2 represents inter-catchment groundwater exchange or losses to deep aquifers:

- **X2 > 0**: Net gain from neighboring catchments or deep groundwater upwelling. This can occur in karst systems or where regional groundwater gradients favor inflow.

- **X2 < 0**: Net loss to adjacent basins or deep percolation beyond the routing store. Common in catchments overlying permeable geology or fractured bedrock.

- **X2 ≈ 0**: Minimal inter-catchment exchange, typical for hydrologically isolated basins with impermeable boundaries.

- **Calibration challenges**: X2 is often poorly identifiable because its effects can be confounded with other water balance components. Consider fixing X2 = 0 if regional hydrogeology suggests minimal exchange, or use tracer studies or water balance analyses to inform priors.

#### X3: Routing Store Capacity (mm)

The routing store governs baseflow dynamics and recession behavior:

- **Low X3 values** (< 50 mm): Rapid drainage, steep recession curves, and limited baseflow sustenance. Characteristic of impermeable catchments or shallow aquifers with high hydraulic conductivity.

- **High X3 values** (> 200 mm): Slow drainage, gentle recessions, and sustained low flows during dry seasons. Typical of deep aquifer systems, permeable formations, or wetland-dominated basins.

- **Recession analysis**: X3 directly influences the recession coefficient. Observed recession curves (plotting $\log Q$ vs. time during dry spells) can inform X3 priors. The recession slope relates to the routing store's time constant.

- **Baseflow index**: Catchments with high baseflow indices (ratio of baseflow to total flow) generally require larger X3 values to sustain groundwater contributions.

#### X4: Unit Hydrograph Time Base (days)

X4 controls the timing and spread of the quickflow response:

- **Low X4 values** (< 2 days): Sharp, peaked hydrograph responses with minimal lag. Appropriate for small, steep catchments with rapid concentration times.

- **High X4 values** (> 10 days): Delayed, attenuated hydrograph peaks with broad distribution. Suitable for large, flat catchments or those with significant channel storage and floodplain attenuation.

- **Time of concentration**: X4 approximates the catchment's time of concentration, though it also reflects channel routing and hillslope travel times. Standard empirical formulas (e.g., Kirpich, SCS) can provide initial estimates.

- **Spatial scale**: X4 typically scales with catchment area, though terrain slope and drainage density also matter. Regionally calibrated relationships between X4 and catchment characteristics (area, slope, stream length) can inform priors.

#### TT: Temperature Threshold (°C)

TT determines when precipitation falls as snow versus rain:

- **Typical range**: -2°C to +2°C for most temperate and alpine catchments. Maritime climates may have TT closer to 2°C, while continental climates favor values near 0°C or below.

- **Mixed precipitation**: Around TT, precipitation is often mixed (rain and snow coexisting). The degree-day model treats this threshold as sharp, which is a simplification but generally adequate for daily modeling.

- **Elevation effects**: TT should reflect the catchment's hypsometry (elevation distribution). Catchments with wide elevation ranges may require more sophisticated precipitation-phase partitioning (e.g., elevation-dependent TT), but the single-value approach works well for many applications.

- **Observational constraints**: If snow depth or snow-covered area data are available, these can help constrain TT through posterior predictive checks on snow accumulation patterns.

#### DDF: Degree-Day Factor (mm °C⁻¹ day⁻¹)

DDF quantifies the snowmelt rate per degree above TT:

- **Typical range**: 2–6 mm °C⁻¹ day⁻¹ for open sites; lower values (1–3) for forested areas where canopy reduces energy input. Alpine snow with high albedo may also show lower DDF.

- **Energy balance approximation**: DDF is an empirical proxy for the energy balance. It implicitly integrates solar radiation, air temperature, humidity, and wind effects. As such, it varies with season, aspect, and forest cover.

- **Literature priors**: Regional snowmelt studies often report DDF distributions. For example, European alpine studies suggest DDF ~ N(3.5, 1.0). Tailor priors to local conditions where possible.

- **Seasonal variability**: Some studies estimate separate DDF values for different melt periods (early spring vs. late spring) to capture changes in snow albedo and energy input. The single-parameter approach here assumes an effective seasonal average.

#### σ: Residual Standard Deviation (m³/s)

σ quantifies the magnitude of model-data mismatch:

- **Sources of residuals**: Includes measurement errors (rating curve uncertainty, sensor noise), model structural errors (missing processes, parameter aggregation), and input data errors (precipitation gauge network gaps, PET estimation).

- **Scaling with discharge**: Homoscedastic σ assumes constant variance across flow magnitudes. In reality, errors often scale with discharge (heteroscedasticity). If residual diagnostics show funnel-shaped patterns, consider heteroscedastic likelihoods (Section 5.3).

- **Informative priors on σ**: Upper bounds should reflect plausible total error. For example, if rating curve uncertainty is ±10% and typical flows are 20 m³/s, an upper bound of σ < 5 m³/s is reasonable. Avoid overly large σ priors that allow the model to "give up" on fitting data.

### Recommended example bounds (adjust to local catchment)

```r
# Hydrologic parameter bounds (X1,X2,X3,X4,TT,DDF)
lb_hydro <- c(1e-3, -5, 1e-3, 0.1, -15, 0.0)
ub_hydro <- c(3000, 10, 2000, 60, 10, 20.0)

# Residual parameter bounds (sigma in m3/s)
lb_resid <- c(1e-3)
ub_resid <- c(200)

# Combined bounds
lb <- c(lb_hydro, lb_resid)
ub <- c(ub_hydro, ub_resid)
```

### Scientific rationale for bounds

- **X1, X3 (mm):** Reflect catchment storage capacity (soil depth × porosity × contributing area). Upper bounds should be large enough to avoid truncating plausible values but can be narrowed using local soil and geology data. For instance, soil surveys providing field capacity and wilting point estimates can inform X1 ranges. Deep aquifer systems revealed by borehole data or groundwater models can guide X3 bounds.

- **X2 (dimensionless):** Exchange coefficient sign and magnitude depend on model formulation; allow a wide range if sign convention is uncertain. If hydrogeological evidence (e.g., groundwater contour maps) indicates minimal exchange, consider narrowing bounds around zero or even fixing X2 = 0 to reduce dimensionality.

- **X4 (days):** Typical values often fall in 1–20 days for many catchments; slow systems may require larger values. Empirical time-of-concentration formulas or analysis of observed storm hydrographs can provide initial guidance. For large catchments (> 1000 km²), X4 may exceed 20 days due to channel routing delays.

- **TT (°C):** Temperature split should cover the local climatology (negative values for cold climates). If the catchment rarely experiences freezing temperatures, narrow the range (e.g., 0–5°C). Conversely, high-elevation or high-latitude sites may require TT as low as -10°C to -15°C.

- **DDF (mm °C⁻¹ day⁻¹):** Degree‑day factors vary with snowpack properties; literature values often lie in a few mm °C⁻¹ day⁻¹. Forest cover, slope aspect, and snow albedo all influence DDF. If regional studies exist, use those ranges. Otherwise, 0–10 mm °C⁻¹ day⁻¹ is a conservative default, though values above 6 are rare.

- **σ (m³/s):** Expressed in discharge units; upper bound should reflect plausible measurement and model error (e.g., a fraction of peak flows). Avoid setting σ upper bounds so high that the model can fit noise. A practical rule: σ_max ~ 0.2 × (95th percentile of observed discharge).

### Parameter Transformation Strategies

For some parameters, working in transformed space improves sampler efficiency:

- **Log-transform X1, X3**: These span orders of magnitude. Sampling $\log(X1)$ and $\log(X3)$ with uniform or normal priors in log-space often yields better exploration than uniform priors in linear space.

- **Logit-transform for bounded parameters**: If you want to constrain X2 to (-5, 5) but sample more efficiently, transform to logit scale.

- **Standardization**: If using Hamiltonian Monte Carlo (e.g., Stan), standardize parameters to have comparable scales (mean 0, SD 1 in prior space) to improve sampler geometry.

Example: sampling X1 in log-space:

```r
# Prior in log-space: log(X1) ~ Uniform(log(1e-3), log(3000))
log_prior_transformed <- function(log_theta, log_lb, log_ub) {
  if (any(log_theta < log_lb) || any(log_theta > log_ub)) return(-Inf)
  # Jacobian adjustment for log-transform
  return(sum(dunif(log_theta, log_lb, log_ub, log = TRUE)) - sum(log_theta))
}
```

**Practical guidance**
- Document bounds and units in the repository README and function docstrings.
- Use prior predictive checks (simulate from priors) to verify that bounds produce physically plausible hydrographs.
- Tighten bounds when identifiability issues arise or when local studies provide reliable constraints.
- When calibrating multiple catchments, consider hierarchical priors where regional parameters inform individual catchment priors, borrowing strength across sites.

---

## 5.2 Log‑Prior Function

**Purpose.** Encode prior knowledge and enforce parameter bounds. Uniform priors inside bounds are a default non‑informative choice; alternatives are recommended when prior information exists.

### Uniform prior (default)

For each parameter $\theta_i$:

$$
p(\theta_i) =
\begin{cases}
\dfrac{1}{\text{ub}_i - \text{lb}_i} & \theta_i \in [\text{lb}_i,\text{ub}_i] \\\\
0 & \text{otherwise}
\end{cases}
$$

**R implementation**

```r
log_prior <- function(theta, lb, ub) {
  if (any(theta < lb) || any(theta > ub)) return(-Inf)
  return(sum(dunif(theta, lb, ub, log = TRUE)))
}
```

### Alternative priors and when to use them

- **Truncated normal**: use when literature suggests a plausible center and uncertainty (e.g., DDF from snow studies).
- **Log‑normal**: use for strictly positive parameters spanning orders of magnitude (X1, X3).
- **Hierarchical priors**: use when calibrating multiple catchments jointly to share information across basins.

### Prior predictive checks

1. Sample $\theta \sim p(\theta)$.
2. Simulate $Q_{\text{sim}}(\theta)$ with the simulator.
3. Inspect simulated hydrographs for physical realism (seasonality, magnitudes).
4. Revise priors/bounds if prior predictive simulations are unrealistic.

---

## 5.3 Likelihood Function

**Purpose.** Quantify the probability of observed discharge given model simulations and residual assumptions. The likelihood choice affects uncertainty quantification and parameter inference.

### Default: Homoscedastic Gaussian

Assume independent, identically distributed Gaussian residuals:

$$Q_{\text{obs},t} \sim \mathcal{N}\big(Q_{\text{sim},t}(\theta),\sigma^2\big)$$

Log‑likelihood:

$$\log L(\theta) = \sum_{t=1}^T \left[ -\frac{1}{2}\log(2\pi\sigma^2) - \frac{(Q_{\text{obs},t}-Q_{\text{sim},t}(\theta))^2}{2\sigma^2} \right]$$

**R implementation (homoscedastic)**

```r
log_likelihood_homo <- function(theta, P_xts, T_xts, PET_xts = NULL,
                                area_km2, Q_obs_xts) {
  hydro_params <- theta[1:6]
  sigma <- theta[7]
  if (sigma <= 0) return(-Inf)

  Q_sim_xts <- gr4j_sim(P_xts, T_xts, PET_xts, hydro_params, area_km2)
  merged <- merge(Q_sim_xts, Q_obs_xts, all = FALSE)
  Q_sim <- as.numeric(merged[,1]); Q_obs <- as.numeric(merged[,2])
  if (length(Q_obs) == 0) return(-Inf)

  ll <- sum(dnorm(Q_obs, mean = Q_sim, sd = sigma, log = TRUE))
  return(ll)
}
```

### Likelihood alternatives

**Heteroscedastic Gaussian** — variance grows with flow:

$$\sigma_t = \sigma_0 (1 + \alpha Q_{\text{sim},t})$$

**AR(1) residuals** — temporal correlation:

$$r_t = Q_{\text{obs},t} - Q_{\text{sim},t},\quad r_t = \phi r_{t-1} + \varepsilon_t,\quad \varepsilon_t\sim\mathcal{N}(0,\sigma^2)$$

**Student‑t** — robust to outliers; parameterized by degrees of freedom $\nu$.

**Log‑space likelihood** — model $\log Q$ when multiplicative errors dominate.

**R examples (heteroscedastic and AR(1))**

```r
# Heteroscedastic example (alpha fixed or estimated)
log_likelihood_hetero <- function(theta, P_xts, T_xts, PET_xts = NULL,
                                  area_km2, Q_obs_xts, alpha = 0.01) {
  hydro_params <- theta[1:6]
  sigma0 <- theta[7]
  if (sigma0 <= 0) return(-Inf)

  Q_sim_xts <- gr4j_sim(P_xts, T_xts, PET_xts, hydro_params, area_km2)
  merged <- merge(Q_sim_xts, Q_obs_xts, all = FALSE)
  Q_sim <- as.numeric(merged[,1]); Q_obs <- as.numeric(merged[,2])

  sd_t <- pmax(1e-6, sigma0 * (1 + alpha * Q_sim))
  ll <- sum(dnorm(Q_obs, mean = Q_sim, sd = sd_t, log = TRUE))
  return(ll)
}

# AR(1) residuals example (phi included in theta as theta[8])
log_likelihood_ar1 <- function(theta, P_xts, T_xts, PET_xts = NULL,
                               area_km2, Q_obs_xts) {
  hydro_params <- theta[1:6]
  sigma <- theta[7]
  phi <- theta[8]
  if (sigma <= 0 || abs(phi) >= 1) return(-Inf)

  Q_sim_xts <- gr4j_sim(P_xts, T_xts, PET_xts, hydro_params, area_km2)
  merged <- merge(Q_sim_xts, Q_obs_xts, all = FALSE)
  Q_sim <- as.numeric(merged[,1]); Q_obs <- as.numeric(merged[,2])
  res <- Q_obs - Q_sim
  ll <- dnorm(res[1], mean = 0, sd = sigma / sqrt(1 - phi^2), log = TRUE)
  for (t in 2:length(res)) {
    mu_t <- phi * res[t-1]
    ll <- ll + dnorm(res[t], mean = mu_t, sd = sigma, log = TRUE)
  }
  return(ll)
}
```

### Likelihood selection checklist

1. Fit with homoscedastic Gaussian.
2. Compute residual diagnostics (Section 5.7).
3. If residual variance increases with flow → heteroscedastic or log‑space likelihood.
4. If residuals show temporal correlation → AR(1) residuals.
5. If heavy tails/outliers → Student‑t likelihood.

---

## 5.4 Log‑Posterior Function

**Definition.** The log‑posterior is the sum of the log‑prior and the log‑likelihood:

$$\log p(\theta \mid Q_{\text{obs}}) = \log L(\theta) + \log p(\theta)$$

**R wrapper (pluggable likelihood)**

```r
log_posterior <- function(theta, lb, ub, P_xts, T_xts, PET_xts,
                          area_km2, Q_obs_xts, likelihood_fn) {
  lp <- log_prior(theta, lb, ub)
  if (is.infinite(lp)) return(-Inf)
  ll <- likelihood_fn(theta, P_xts, T_xts, PET_xts, area_km2, Q_obs_xts)
  return(lp + ll)
}
```

### Scientific considerations

- **Posterior geometry**: multimodality, ridges, and narrow valleys affect sampler choice and tuning.
- **Identifiability**: strong posterior correlations (e.g., between X1 and DDF) indicate structural non‑identifiability; consider reparameterization, fixing parameters, or informative priors.
- **Transformations**: log‑transform positive parameters to stabilize sampling and reduce skewness.

---

## 5.5 MCMC Sampling

**Sampler selection.** Adaptive Metropolis (e.g., `adaptMCMC`) is a practical default for moderate dimensionality. For complex posteriors, consider DE‑MCMC or Hamiltonian Monte Carlo (`rstan`, `cmdstanr`).

### Example using `adaptMCMC`

```r
library(adaptMCMC)

# initial values: X1,X2,X3,X4,TT,DDF,sigma
theta_init <- c(500, 0.5, 200, 5, 0, 3, 1)

# scale initialization proportional to parameter ranges
scale_init <- (ub - lb) * 0.02
scale_init[scale_init <= 0] <- 0.1

mcmc_result <- MCMC(
  p = function(theta, ...) log_posterior(theta, lb, ub, ...,
                                         likelihood_fn = log_likelihood_homo),
  init = theta_init,
  n = 30000,
  adapt = TRUE,
  acc.rate = 0.234,
  scale = scale_init,
  P_xts = P_xts,
  T_xts = T_xts,
  PET_xts = PET_xts,
  area_km2 = area_km2,
  Q_obs_xts = Q_obs_xts
)

burn_in <- 5000
samples <- mcmc_result$samples[(burn_in + 1):nrow(mcmc_result$samples), ]
colnames(samples) <- c("X1","X2","X3","X4","TT","DDF","sigma")
```

### MCMC best practices

- **Multiple chains**: run at least 3 chains with dispersed initial values.
- **Convergence diagnostics**: compute R̂ (Gelman–Rubin), inspect traceplots, and check effective sample size (ESS).
- **Acceptance rate**: target ~0.2–0.3 for random‑walk Metropolis; adaptMCMC tunes proposals automatically.
- **Posterior exploration**: visualize marginal densities and pairwise scatterplots to detect multimodality or ridges.
- **Reproducibility**: set random seeds, save raw chains, and record `sessionInfo()`.

### Computational tips

- **Parallelize chains** to reduce wall time.
- **Profile** `gr4j_sim` to identify bottlenecks; vectorize or optimize inner loops where possible.
- **Use short test runs** to tune sampler settings before long runs.

---

## 5.6 Posterior Predictive Checks (PPCs)

**Purpose.** Evaluate whether posterior samples reproduce observed flows and quantify predictive uncertainty.

Posterior predictive checks are the cornerstone of Bayesian model validation. They answer the fundamental question: **"If the model and fitted parameters are correct, could they have generated data like what we actually observed?"** PPCs go beyond traditional goodness-of-fit metrics by:

1. **Assessing model adequacy**: Do simulations capture the full range of observed behaviors (extremes, seasonality, autocorrelation)?
2. **Quantifying predictive uncertainty**: Are prediction intervals well-calibrated (coverage matches nominal levels)?
3. **Diagnosing model deficiencies**: Where does the model systematically fail, and what processes might be missing?

### Theoretical Foundation

The posterior predictive distribution is:

$$p(Q_{\text{pred}} \mid Q_{\text{obs}}) = \int p(Q_{\text{pred}} \mid \theta) \, p(\theta \mid Q_{\text{obs}}) \, d\theta$$

This integrates over parameter uncertainty: for each posterior sample $\theta^{(i)}$, we simulate $Q_{\text{sim}}^{(i)}$ and then draw $Q_{\text{pred}}^{(i)}$ from the likelihood. The ensemble of $Q_{\text{pred}}^{(i)}$ represents our predictive uncertainty, combining:

- **Parameter uncertainty**: Different $\theta$ values produce different simulated flows.
- **Residual variability**: The likelihood adds noise (e.g., Gaussian with $\sigma^{(i)}$) around each simulation.

### Procedure

1. Draw posterior samples $\theta^{(i)}$ from the MCMC chain (thin if necessary to reduce autocorrelation).
2. For each $\theta^{(i)}$, run the GR4J simulator to obtain $Q_{\text{sim}}^{(i)}$.
3. Draw predictive replicates $Q_{\text{pred},t}^{(i)} \sim p(Q_t \mid \theta^{(i)})$ using the likelihood. For homoscedastic Gaussian errors:

$$Q_{\text{pred},t}^{(i)} = Q_{\text{sim},t}^{(i)} + \epsilon_t^{(i)}, \quad \epsilon_t^{(i)} \sim \mathcal{N}(0, [\sigma^{(i)}]^2)$$

4. Compare observed $Q_{\text{obs}}$ to the ensemble of $Q_{\text{pred}}^{(i)}$ visually and with summary statistics.

### R code (ensemble PPC)

```r
n_pred <- 200  # Number of posterior samples to use for predictions
pred_idx <- sample(1:nrow(samples), n_pred)
merged_obs <- merge(gr4j_sim(P_xts, T_xts, PET_xts, samples[1,1:6], area_km2),
                    Q_obs_xts, all = FALSE)
time_steps <- seq_len(nrow(merged_obs))
Q_pred_samples <- matrix(NA, nrow = n_pred, ncol = length(time_steps))

for (i in seq_len(n_pred)) {
  theta_i <- samples[pred_idx[i], ]
  hydro_params <- theta_i[1:6]
  sigma_i <- theta_i[7]

  Q_sim_xts <- gr4j_sim(P_xts, T_xts, PET_xts, hydro_params, area_km2)
  merged <- merge(Q_sim_xts, Q_obs_xts, all = FALSE)
  Q_sim <- as.numeric(merged[,1])
  
  # Add residual noise to create predictive replicates
  Q_pred_samples[i, ] <- rnorm(length(Q_sim), mean = Q_sim, sd = sigma_i)
}

# Compute summary statistics
Q_mean_pred <- colMeans(Q_pred_samples)
Q_median_pred <- apply(Q_pred_samples, 2, median)
Q_lower_50 <- apply(Q_pred_samples, 2, quantile, probs = 0.25)
Q_upper_50 <- apply(Q_pred_samples, 2, quantile, probs = 0.75)
Q_lower_90 <- apply(Q_pred_samples, 2, quantile, probs = 0.05)
Q_upper_90 <- apply(Q_pred_samples, 2, quantile, probs = 0.95)
```

### Visual PPC Diagnostics

#### Time Series Plot with Uncertainty Bands

This is the most informative PPC plot: observed flows overlaid on predictive intervals.

```r
# Extract time index from merged data
time_index <- index(merged_obs)
Q_obs <- as.numeric(merged_obs[,2])

# Plot setup
par(mfrow = c(1,1), mar = c(4,4,2,1))
plot(time_index, Q_obs, type = "n", ylim = range(c(Q_obs, Q_pred_samples), na.rm = TRUE),
     xlab = "Time", ylab = "Discharge (m3/s)", 
     main = "Posterior Predictive Check: Observed vs Predicted")

# 90% credible interval (light shading)
polygon(c(time_index, rev(time_index)), 
        c(Q_lower_90, rev(Q_upper_90)),
        col = rgb(0.7, 0.7, 0.7, 0.3), border = NA)

# 50% credible interval (darker shading)
polygon(c(time_index, rev(time_index)), 
        c(Q_lower_50, rev(Q_upper_50)),
        col = rgb(0.5, 0.5, 0.5, 0.4), border = NA)

# Posterior mean prediction
lines(time_index, Q_mean_pred, col = "blue", lwd = 2)

# Observed data
points(time_index, Q_obs, pch = 19, cex = 0.5, col = "black")

legend("topright", 
       c("Observed", "Posterior Mean", "50% Interval", "90% Interval"),
       col = c("black", "blue", rgb(0.5,0.5,0.5,0.4), rgb(0.7,0.7,0.7,0.3)),
       pch = c(19, NA, 15, 15), lty = c(NA, 1, NA, NA), lwd = c(NA, 2, NA, NA))
```

**Interpretation:**
- Observations should be scattered throughout the predictive bands, not systematically above or below.
- Narrower bands indicate lower uncertainty (well-constrained parameters and small $\sigma$).
- If observations frequently fall outside 90% intervals, the model is overconfident (underestimating uncertainty).

#### Spaghetti Plot: Overlay of Multiple Realizations

Instead of summary intervals, plot a subset of individual predictive trajectories.

```r
plot(time_index, Q_obs, type = "l", lwd = 2, col = "black",
     ylim = range(c(Q_obs, Q_pred_samples[1:50,]), na.rm = TRUE),
     xlab = "Time", ylab = "Discharge (m3/s)",
     main = "Spaghetti Plot: 50 Posterior Predictive Realizations")

# Overlay 50 random realizations
for (i in 1:50) {
  lines(time_index, Q_pred_samples[i,], col = rgb(0, 0, 1, 0.1))
}

# Re-draw observed on top
lines(time_index, Q_obs, lwd = 2, col = "black")
```

**Interpretation:**
- The "cloud" of blue lines represents model spread.
- Observed trajectory should weave through the cloud naturally.
- Dense regions indicate where model is confident; sparse regions indicate high uncertainty.

#### Flow Duration Curve (FDC) Comparison

FDCs reveal whether the model captures the full distribution of flows, especially extremes.

```r
# Compute empirical FDC for observed data
Q_obs_sorted <- sort(Q_obs, decreasing = TRUE)
exceedance_prob <- seq_along(Q_obs_sorted) / length(Q_obs_sorted)

# Compute FDCs for each predictive realization
Q_pred_fdc <- apply(Q_pred_samples, 1, function(x) sort(x, decreasing = TRUE))

# Plot observed FDC
plot(exceedance_prob, Q_obs_sorted, type = "l", lwd = 3, col = "black",
     log = "y", xlab = "Exceedance Probability", ylab = "Discharge (m3/s)",
     main = "Flow Duration Curve: Observed vs Predictive Ensemble")

# Overlay predictive FDCs (thinned for clarity)
for (i in seq(1, n_pred, by = 10)) {
  lines(exceedance_prob, Q_pred_fdc[,i], col = rgb(0, 0, 1, 0.2))
}

# Re-draw observed on top
lines(exceedance_prob, Q_obs_sorted, lwd = 3, col = "black")

legend("topright", c("Observed", "Predictive Ensemble"),
       col = c("black", "blue"), lwd = c(3, 1))
```

**Interpretation:**
- **High flows** (left side, low exceedance): Model should capture peak flow magnitudes. Systematic underprediction indicates missing flood processes or parameter issues.
- **Low flows** (right side, high exceedance): Model should reproduce baseflow levels. Overprediction suggests excessive groundwater contribution (X3 too large).
- **Mid-range flows**: Check for bias in the central tendency.

### PPC diagnostics and interpretation

#### Coverage Analysis

Quantify the fraction of observations falling within predictive intervals:

```r
# Compute coverage for 50% and 90% intervals
coverage_50 <- mean(Q_obs >= Q_lower_50 & Q_obs <= Q_upper_50, na.rm = TRUE)
coverage_90 <- mean(Q_obs >= Q_lower_90 & Q_obs <= Q_upper_90, na.rm = TRUE)

cat("50% Credible Interval Coverage:", round(coverage_50, 3), "\n")
cat("90% Credible Interval Coverage:", round(coverage_90, 3), "\n")
```

**Expected coverage:** 50% interval should contain ~50% of observations; 90% interval should contain ~90%.

**Deviations:**
- **Undercoverage** (e.g., 30% in 50% interval, 70% in 90% interval): Model is overconfident. Possible causes:
  - $\sigma$ is underestimated (too small).
  - Model structure is inadequate (missing processes).
  - Likelihood assumption is violated (e.g., heteroscedasticity not accounted for).
  
- **Overcoverage** (e.g., 70% in 50% interval, 98% in 90% interval): Model is too uncertain. Possible causes:
  - $\sigma$ is overestimated.
  - Priors are too diffuse, allowing unrealistic parameter combinations.
  - Data are more informative than model structure allows (overparameterization).

#### Conditional Coverage: Coverage by Flow Regime

Check if coverage varies across flow magnitudes:

```r
# Divide observations into terciles (low, medium, high flows)
Q_terciles <- quantile(Q_obs, probs = c(1/3, 2/3), na.rm = TRUE)
low_flow_idx <- Q_obs <= Q_terciles[1]
mid_flow_idx <- Q_obs > Q_terciles[1] & Q_obs <= Q_terciles[2]
high_flow_idx <- Q_obs > Q_terciles[2]

# Coverage by regime
coverage_50_low <- mean(Q_obs[low_flow_idx] >= Q_lower_50[low_flow_idx] & 
                        Q_obs[low_flow_idx] <= Q_upper_50[low_flow_idx], na.rm = TRUE)
coverage_50_high <- mean(Q_obs[high_flow_idx] >= Q_lower_50[high_flow_idx] & 
                         Q_obs[high_flow_idx] <= Q_upper_50[high_flow_idx], na.rm = TRUE)

cat("50% Coverage (Low Flows):", round(coverage_50_low, 3), "\n")
cat("50% Coverage (High Flows):", round(coverage_50_high, 3), "\n")
```

**Interpretation:**
- **Poor high-flow coverage**: Model struggles with extremes. Consider heteroscedastic likelihood or additional flood processes.
- **Poor low-flow coverage**: Baseflow simulation is biased or uncertain. Check X3, consider recession analysis.

#### Tail Behavior: Extreme Event Capture

Extract maximum observed flow and compare to predictive distribution of maxima:

```r
Q_obs_max <- max(Q_obs, na.rm = TRUE)
Q_pred_max <- apply(Q_pred_samples, 1, max, na.rm = TRUE)

hist(Q_pred_max, breaks = 30, col = "lightblue", border = "white",
     xlab = "Maximum Discharge (m3/s)", main = "Predictive Distribution of Peak Flow")
abline(v = Q_obs_max, col = "red", lwd = 3, lty = 2)
legend("topright", "Observed Max", col = "red", lwd = 3, lty = 2)

# Compute p-value: what fraction of predictive maxima exceed observed max?
p_val_max <- mean(Q_pred_max >= Q_obs_max)
cat("P-value (observed max vs predictive):", round(p_val_max, 3), "\n")
```

**Interpretation:**
- If $p < 0.05$, observed maximum is unusually high given the model—suggests missing flood processes or underestimated parameter uncertainty.
- If $p > 0.95$, observed maximum is unusually low—model may be overestimating flood risk or parameters are biased high.

### Predictive Performance Metrics

Compute NSE and KGE on predictive ensembles:

**Nash–Sutcliffe Efficiency (NSE)**

$$\text{NSE} = 1 - \frac{\sum_t (Q_{\text{obs},t} - Q_{\text{sim},t})^2}{\sum_t (Q_{\text{obs},t} - \overline{Q}_{\text{obs}})^2}$$

For Bayesian ensembles, compute NSE for each posterior sample and report the distribution:

```r
NSE_ensemble <- numeric(n_pred)
Q_obs_mean <- mean(Q_obs, na.rm = TRUE)

for (i in 1:n_pred) {
  Q_sim_i <- Q_pred_samples[i, ]  # or use deterministic Q_sim without noise for NSE
  NSE_ensemble[i] <- 1 - sum((Q_obs - Q_sim_i)^2, na.rm = TRUE) / 
                         sum((Q_obs - Q_obs_mean)^2, na.rm = TRUE)
}

cat("NSE Median:", round(median(NSE_ensemble), 3), "\n")
cat("NSE 90% CI:", round(quantile(NSE_ensemble, c(0.05, 0.95)), 3), "\n")
```

**Kling–Gupta Efficiency (KGE)** (one common formulation)

$$\text{KGE} = 1 - \sqrt{(r-1)^2 + (\alpha-1)^2 + (\beta-1)^2}$$

where $r$ is correlation, $\alpha = \frac{\sigma_{\text{sim}}}{\sigma_{\text{obs}}}$, and $\beta = \frac{\mu_{\text{sim}}}{\mu_{\text{obs}}}$.

```r
KGE_ensemble <- numeric(n_pred)

for (i in 1:n_pred) {
  Q_sim_i <- Q_pred_samples[i, ]
  r <- cor(Q_obs, Q_sim_i, use = "complete.obs")
  alpha <- sd(Q_sim_i, na.rm = TRUE) / sd(Q_obs, na.rm = TRUE)
  beta <- mean(Q_sim_i, na.rm = TRUE) / mean(Q_obs, na.rm = TRUE)
  KGE_ensemble[i] <- 1 - sqrt((r - 1)^2 + (alpha - 1)^2 + (beta - 1)^2)
}

cat("KGE Median:", round(median(KGE_ensemble), 3), "\n")
cat("KGE 90% CI:", round(quantile(KGE_ensemble, c(0.05, 0.95)), 3), "\n")
```

**Interpretation:**
- NSE and KGE both range from $-\infty$ to 1 (perfect fit).
- NSE > 0.5 and KGE > 0.5 are often considered acceptable for hydrological models.
- Report median and credible intervals to show uncertainty in performance metrics themselves.

### Seasonal and Regime-Specific PPCs

Hydrological models often fail in specific seasons or flow regimes. Subset PPCs by:

- **Season**: Calibrate on full year, but check winter vs. summer performance separately.
- **Snowmelt period**: For snow-dominated catchments, isolate April-June and check if melt timing and magnitude are captured.
- **Recession periods**: Subset dry spells and check baseflow recession behavior.

Example: PPC for snowmelt season only:

```r
# Assume time_index is a Date or POSIXct object
snowmelt_months <- format(time_index, "%m") %in% c("04", "05", "06")

Q_obs_snowmelt <- Q_obs[snowmelt_months]
Q_pred_snowmelt <- Q_pred_samples[, snowmelt_months]

# Repeat plots and coverage analysis for snowmelt subset
```

### Advanced PPC: Test Statistics

Instead of raw time series, compute summary statistics on observed vs. predictive data:

- **Autocorrelation at lag 1**: Does model capture persistence?
- **Coefficient of variation**: Does model reproduce flow variability?
- **Skewness**: Does model match the asymmetry of the flow distribution?

```r
# Compute test statistic: lag-1 autocorrelation
acf_obs <- acf(Q_obs, lag.max = 1, plot = FALSE)$acf[2]
acf_pred <- apply(Q_pred_samples, 1, function(x) acf(x, lag.max = 1, plot = FALSE)$acf[2])

hist(acf_pred, breaks = 30, col = "lightblue", border = "white",
     xlab = "Lag-1 Autocorrelation", main = "PPC: Autocorrelation Structure")
abline(v = acf_obs, col = "red", lwd = 3, lty = 2)
legend("topright", "Observed ACF(1)", col = "red", lwd = 3, lty = 2)

# p-value: fraction of predictive ACF more extreme than observed
p_val_acf <- mean(abs(acf_pred - mean(acf_pred)) >= abs(acf_obs - mean(acf_pred)))
cat("P-value (autocorrelation):", round(p_val_acf, 3), "\n")
```

**Interpretation:**
- If observed ACF falls in the tails of the predictive distribution, the model may not capture temporal dependencies well.
- Consider AR(1) residuals or additional routing dynamics if autocorrelation is systematically mismatched.

### Reporting PPC Results

A comprehensive PPC report should include:

1. **Time series plot** with 50% and 90% credible intervals.
2. **Coverage statistics** (overall and by flow regime).
3. **Flow duration curve** comparison.
4. **NSE and KGE** distributions (median and 90% CI).
5. **Extreme event analysis** (maxima, minima).
6. **Seasonal subset** PPCs if relevant.
7. **Test statistics** (autocorrelation, variance, skewness) if model assumptions are questionable.

Example summary table:

| **Metric** | **Observed** | **Predictive Median** | **Predictive 90% CI** |
|---|---:|---:|---:|
| Mean Q (m³/s) | 15.2 | 15.1 | [14.5, 15.8] |
| Max Q (m³/s) | 89.3 | 85.7 | [72.1, 102.4] |
| 50% Coverage | — | 52% | — |
| 90% Coverage | — | 88% | — |
| NSE | — | 0.78 | [0.72, 0.83] |
| KGE | — | 0.81 | [0.76, 0.85] |

### When PPCs Fail: Next Steps

If PPCs reveal systematic deficiencies:

1. **Heteroscedasticity**: Adopt heteroscedastic or log-space likelihood (Section 5.3).
2. **Autocorrelation**: Include AR(1) residuals.
3. **Missing processes**: Add snow module parameters, improve PET estimation, or consider distributed routing.
4. **Parameter identifiability**: Strong correlations or wide posteriors suggest fixing some parameters or using informative priors.
5. **Structural inadequacy**: GR4J may be too simple for the catchment. Consider more complex models (e.g., SWAT, VIC) or hybrid approaches.

---

## 5.7 Residual Diagnostics

**Purpose.** Validate likelihood assumptions and guide model or likelihood refinement.

Residual diagnostics are essential for assessing whether the statistical assumptions underlying the Bayesian inference are met. Violations of these assumptions (e.g., non-constant variance, autocorrelation, non-normality) can lead to:

- **Biased parameter estimates**: Parameters may be pulled toward values that compensate for unmodeled error structures.
- **Incorrect uncertainty quantification**: Credible intervals and predictive intervals may be too narrow or too wide.
- **Poor out-of-sample performance**: Models calibrated with incorrect error assumptions often fail when applied to new data or future periods.

### Residual definition

$$r_t = Q_{\text{obs},t} - \bar{Q}_{\text{pred},t}$$

where $\bar{Q}_{\text{pred},t}$ is the posterior predictive mean (or median; choice depends on reporting preference).

**Compute residuals (aligned index)**

```r
residuals <- as.numeric(merged_obs[,2]) - Q_mean_pred
```

Note: For heteroscedastic or AR(1) likelihoods, compute **standardized residuals** to account for varying variances or temporal structure.

### 5.7.1 Residual vs Fitted

This plot reveals heteroscedasticity (variance changing with discharge magnitude).

```r
plot(Q_mean_pred, residuals, pch=19, col="steelblue", cex = 0.7,
     xlab="Posterior mean prediction (m3/s)", ylab="Residuals (m3/s)",
     main="Residuals vs Posterior Mean Prediction")
abline(h=0, col="red", lty=2, lwd=2)

# Add loess smoother to detect systematic bias
loess_fit <- loess(residuals ~ Q_mean_pred)
lines(sort(Q_mean_pred), predict(loess_fit, sort(Q_mean_pred)), 
      col = "orange", lwd = 2)
```

**Interpretation:**

- **Centered residuals** around zero (red line) indicate unbiased predictions on average.
- **Loess curve deviating from zero**: Systematic bias. For example:
  - Curve above zero at high flows → model underpredicts peaks.
  - Curve below zero at low flows → model overpredicts baseflow.
  
- **Funnel pattern** (variance increasing with $Q_{\text{pred}}$): Classic sign of heteroscedasticity. Residual variance grows with discharge magnitude, violating the homoscedastic Gaussian assumption. 
  - **Action**: Adopt heteroscedastic likelihood where $\sigma_t = \sigma_0 (1 + \alpha Q_{\text{sim},t})$, or work in log-space.
  
- **Reverse funnel** (variance decreasing with $Q_{\text{pred}}$): Rare but possible; may indicate measurement errors are proportionally larger at low flows.

- **Outliers**: Points far from the zero line. Investigate specific events:
  - Extreme rainfall events not captured by precipitation data?
  - Anthropogenic influences (dam releases, abstractions)?
  - Rating curve errors at extreme stages?

### 5.7.2 Residual vs Time and Autocorrelation

Temporal plots reveal trends, seasonality, and autocorrelation.

```r
# Residual time series
plot(1:length(residuals), residuals, type="b", pch=19, col="darkgreen", cex=0.5,
     xlab="Time step", ylab="Residuals (m3/s)", main="Residuals vs Time")
abline(h=0, col="red", lty=2, lwd=2)

# Add seasonal loess smoother
loess_time <- loess(residuals ~ seq_along(residuals), span = 0.1)
lines(seq_along(residuals), predict(loess_time), col = "blue", lwd = 2)
```

**Interpretation:**

- **Random scatter around zero**: No systematic temporal bias—good!
- **Seasonal pattern**: Loess curve shows periodic swings (e.g., consistent overprediction in summer, underprediction in winter). Possible causes:
  - Seasonal bias in PET estimation (e.g., overestimating summer ET).
  - Missing snow dynamics (if TT and DDF are poorly constrained).
  - Seasonal human influences (irrigation, reservoir operations).
  
- **Drift or trend**: Residuals trending upward or downward over time. May indicate:
  - Non-stationarity in catchment behavior (land use change, climate trends).
  - Model warm-up issues (initial conditions not equilibrated).

**Autocorrelation Function (ACF)**

```r
acf(residuals, main="Autocorrelation of Residuals", lag.max = 30)
```

**Interpretation:**

- **No significant lags beyond lag 0**: Residuals are approximately independent—homoscedastic Gaussian assumption is reasonable.
- **Significant ACF at lag 1 or higher** (bars exceeding blue dashed confidence bands): Residuals are temporally correlated. The model is not capturing all persistence in the system.
  - **Action**: Include AR(1) residuals in the likelihood (Section 5.3), or add missing process memory (e.g., deeper groundwater store, multi-store routing).

**Partial Autocorrelation Function (PACF)**

```r
pacf(residuals, main="Partial Autocorrelation of Residuals", lag.max = 30)
```

PACF shows direct correlation at each lag, controlling for shorter lags. If PACF is significant only at lag 1, an AR(1) model is appropriate. If multiple lags are significant, consider higher-order AR models.

**Ljung-Box Test**

Formal test for autocorrelation:

```r
Box.test(residuals, lag = 10, type = "Ljung-Box")
```

- **p-value < 0.05**: Significant autocorrelation detected. Reject the hypothesis of independent residuals.
- **p-value ≥ 0.05**: No strong evidence of autocorrelation.

### 5.7.3 QQ Plot and Histogram

These assess normality of residuals.

#### QQ Plot

```r
qqnorm(residuals, pch=19, col="steelblue", main="QQ Plot of Residuals", cex=0.7)
qqline(residuals, col="red", lty=2, lwd=2)
```

**Interpretation:**

- **Points along the diagonal line**: Residuals follow a normal distribution—Gaussian likelihood is appropriate.
- **Heavy upper tail** (points curve above line at right end): More extreme positive residuals than normal distribution predicts. Model underestimates large flows more severely than Gaussian errors allow.
  - **Action**: Consider Student-t likelihood (robust to outliers) or heteroscedastic likelihood.
  
- **Heavy lower tail** (points curve below line at left end): More extreme negative residuals. Model overestimates low flows more than expected.
  
- **S-shaped curve**: Skewed distribution. If upper tail is heavier, consider log-transformation or heteroscedastic errors. If lower tail is heavier, model may systematically overpredict in some regime.

- **Systematic deviation from line**: Non-normal residuals. May indicate:
  - Wrong likelihood family (not Gaussian).
  - Structural model error producing systematic patterns.

#### Histogram

```r
hist(residuals, breaks=30, col="lightblue", border="white", freq=FALSE,
     main="Histogram of Residuals", xlab="Residuals (m3/s)")

# Overlay normal distribution with mean=0, sd=empirical sd
curve(dnorm(x, mean = 0, sd = sd(residuals)), add = TRUE, col = "red", lwd = 2)
```

**Interpretation:**

- **Symmetric, bell-shaped**: Consistent with normality.
- **Long tails**: Indicates outliers; Student-t may be better.
- **Skewness**: Asymmetric distribution. Positive skew (long right tail) is common in hydrology (extreme floods). Log-transformation or heteroscedastic likelihood can help.
- **Bimodal**: Two peaks suggest different error regimes (e.g., snowmelt vs. rainfall-dominated periods). Consider regime-specific calibration or mixture models.

### Shapiro-Wilk Test for Normality

```r
shapiro.test(residuals)
```

- **p-value < 0.05**: Residuals significantly deviate from normality. Consider alternative likelihoods.
- **p-value ≥ 0.05**: No strong evidence against normality.

**Caution**: Shapiro-Wilk is sensitive to sample size. With thousands of data points, minor deviations from normality may be statistically significant but practically unimportant. Always combine formal tests with visual diagnostics.

### 5.7.4 Residual Variance by Flow Regime

Stratify residuals by flow magnitude to detect heteroscedasticity more explicitly:

```r
# Define flow regimes based on quantiles
Q_quantiles <- quantile(Q_mean_pred, probs = c(0.33, 0.67), na.rm = TRUE)
regime <- cut(Q_mean_pred, breaks = c(-Inf, Q_quantiles, Inf), 
              labels = c("Low", "Medium", "High"))

boxplot(residuals ~ regime, col = c("lightblue", "lightgreen", "lightcoral"),
        xlab = "Flow Regime", ylab = "Residuals (m3/s)",
        main = "Residual Variance by Flow Regime")
abline(h = 0, col = "red", lty = 2, lwd = 2)
```

**Interpretation:**

- **Similar box widths**: Variance is roughly constant across flow regimes—homoscedastic assumption holds.
- **Widening boxes from low to high flows**: Clear heteroscedasticity. Adopt heteroscedastic or log-space likelihood.
- **Median lines off-center**: Systematic bias in specific regimes. For example, median above zero in "High" regime means model underpredicts peaks.

**Levene's Test for Homogeneity of Variance**

```r
library(car)
leveneTest(residuals ~ regime)
```

- **p-value < 0.05**: Variance differs significantly across regimes. Heteroscedasticity present.
- **p-value ≥ 0.05**: No strong evidence of variance differences.

### 5.7.5 Seasonal Residual Analysis

For snow-dominated or strongly seasonal catchments, examine residuals by month or season:

```r
# Assume time_index is available as Date or POSIXct
month <- as.numeric(format(time_index, "%m"))

boxplot(residuals ~ month, col = rainbow(12), 
        xlab = "Month", ylab = "Residuals (m3/s)",
        main = "Residual Distribution by Month", names = month.abb)
abline(h = 0, col = "red", lty = 2, lwd = 2)
```

**Interpretation:**

- **Consistent medians near zero across months**: No seasonal bias.
- **Positive residuals in spring**: Model underpredicts snowmelt. Check TT and DDF; consider increasing DDF or adjusting TT.
- **Negative residuals in summer**: Model overpredicts summer flows, possibly due to overestimated PET or underestimated ET efficiency.
- **Varying box widths by season**: Seasonal heteroscedasticity. May require season-specific $\sigma$ or complex error models.

### Quantitative diagnostics to report

A comprehensive residual diagnostic summary should include:

1. **Mean residual** (should be near zero): $\bar{r} = \frac{1}{T}\sum_t r_t$
2. **Standard deviation of residuals**: Compare to posterior mean of $\sigma$. If residual SD $\gg \hat{\sigma}$, model is overconfident.
3. **NSE** and **KGE** computed on posterior predictive ensembles (report median and credible intervals).
4. **Coverage** of predictive intervals (50%, 90%).
5. **Autocorrelation** statistics:
   - Lag-1 ACF value
   - Ljung-Box test p-value
6. **Normality tests**:
   - Shapiro-Wilk p-value (with caution on large samples)
   - QQ plot deviations (qualitative)
7. **Heteroscedasticity tests**:
   - Levene's test p-value
   - Visual funnel pattern in residual vs. fitted plot

Example summary table:

| **Diagnostic** | **Value** | **Interpretation** |
|---|---:|---|
| Mean residual (m³/s) | -0.12 | Near-zero, minimal bias |
| SD of residuals (m³/s) | 2.34 | Slightly higher than median $\hat{\sigma}$ (2.10) |
| Lag-1 ACF | 0.28 | Moderate autocorrelation detected |
| Ljung-Box p-value | 0.003 | Significant autocorrelation (p < 0.05) |
| Shapiro-Wilk p-value | 0.08 | No strong deviation from normality |
| Levene's test p-value | 0.02 | Significant heteroscedasticity (p < 0.05) |
| 50% Coverage | 52% | Well-calibrated |
| 90% Coverage | 86% | Slight undercoverage |

### Actionable responses to diagnostics

Based on residual diagnostics, take these corrective actions:

- **Heteroscedasticity** (funnel pattern, Levene's test p < 0.05): 
  - Adopt heteroscedastic Gaussian likelihood: $\sigma_t = \sigma_0 (1 + \alpha Q_{\text{sim},t})$
  - Or use log-space likelihood: model $\log(Q)$ instead of $Q$
  - Or use multiplicative errors: $Q_{\text{obs},t} = Q_{\text{sim},t} \cdot \exp(\epsilon_t)$ where $\epsilon_t \sim \mathcal{N}(0, \sigma^2)$

- **Autocorrelation** (ACF significant, Ljung-Box p < 0.05):
  - Include AR(1) residual structure: $r_t = \phi r_{t-1} + \varepsilon_t$
  - Or add missing process memory (e.g., additional groundwater store, multi-timescale routing)
  - Check if temporal aggregation (e.g., monthly instead of daily) reduces autocorrelation

- **Poor extremes** (heavy tails in QQ plot, high flows poorly captured):
  - Adopt Student-t likelihood with estimated degrees of freedom $\nu$
  - Or use Box-Cox transformation to normalize residuals
  - Or implement regime-specific calibration (separate parameters for flood vs. baseflow periods)

- **Strong parameter correlations** (check pairwise posterior scatterplots):
  - Reparameterize model (e.g., combine X1 and DDF into a composite "effective storage" parameter)
  - Fix less identifiable parameters based on literature or regional studies
  - Use informative priors to break degeneracies

- **Seasonal bias** (residuals systematically positive/negative in certain months):
  - Improve PET estimation (use alternative methods like Penman-Monteith instead of Oudin)
  - Check snow module: TT and DDF may be season-dependent
  - Investigate anthropogenic influences (irrigation, reservoir operations) and include them explicitly if data are available

- **Systematic bias in flow regimes** (consistent over/underprediction at high or low flows):
  - High flow bias: Check precipitation input quality, consider flood routing enhancements
  - Low flow bias: Examine baseflow separation, adjust X3 priors, check for groundwater abstractions
  - Consider dual-objective calibration (e.g., weight NSE for high flows and log-NSE for low flows)

### Advanced Diagnostic: Recursive Residuals

For detecting non-stationarity, compute recursive residuals (one-step-ahead prediction errors):

```r
# Placeholder concept: use rolling window to check if residual properties change over time
window_size <- 365  # one year
n_windows <- floor(length(residuals) / window_size)

mean_by_window <- sapply(1:n_windows, function(i) {
  idx <- ((i-1)*window_size + 1):(i*window_size)
  mean(residuals[idx], na.rm = TRUE)
})

plot(1:n_windows, mean_by_window, type = "b", pch = 19, col = "purple",
     xlab = "Year", ylab = "Mean Residual (m3/s)",
     main = "Temporal Evolution of Mean Residual")
abline(h = 0, col = "red", lty = 2)
```

**Interpretation:**
- Drift in mean residuals over time suggests non-stationarity (climate change, land use change).
- Action: Consider time-varying parameters or split calibration into sub-periods.

### Residual Diagnostics in Model Selection

Use residual diagnostics to compare alternative model structures or likelihood assumptions:

1. Calibrate with homoscedastic Gaussian likelihood.
2. Compute residual diagnostics.
3. If diagnostics suggest heteroscedasticity, recalibrate with heteroscedastic likelihood.
4. Compare:
   - Posterior predictive coverage (should improve)
   - Residual plots (funnel pattern should disappear)
   - Information criteria (DIC, WAIC) if available

Similarly, if autocorrelation is detected:

1. Recalibrate with AR(1) likelihood.
2. Check ACF of residuals from AR(1) model—should be white noise.
3. Compare coverage and predictive performance.

### Summary: Iterative Refinement

Residual diagnostics are not a one-time check but part of an **iterative model development cycle**:

1. **Calibrate** with simplest assumptions (homoscedastic Gaussian).
2. **Diagnose** residuals to identify violations.
3. **Refine** likelihood or model structure.
4. **Re-calibrate** and re-diagnose.
5. **Repeat** until residuals approximate i.i.d. Gaussian (or other assumed distribution).

This process ensures that Bayesian inference is based on valid statistical foundations, leading to reliable parameter estimates and credible predictive intervals.
