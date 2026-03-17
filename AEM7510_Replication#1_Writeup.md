# AEM 7510 Replication Project

## Burke, Hsiang & Miguel (2015) — Assignment Write-Up

**Name:** Yifan Luo (yl3699@cornell.edu)
**Date:** March 9, 2026

---

## Question 1: Does the Replication Package Exactly Replicate Main-Text Results?

### Results That Do Replicate

| Main-text claim | Replication result | Exact? |
|:---|---:|:---:|
| Temperature optimum at 13°C | 13.1°C | Yes |
| Sample size N = 6,584 | 6,584 | Yes |
| 166 countries in sample | 166 | Yes |
| β₁ (temp) positive, *p* < 0.001 | +0.01272, SE = 0.00374, *t* = 3.40 | Yes |
| β₂ (temp²) negative, *p* < 0.001 | −0.000487, SE = 0.000117, *t* = −4.17 | Yes |
| 77% of countries poorer in 2100 | 77.3% | Yes |
| 5% of countries poorer than today (SSP5) | 5.52% | Yes |
| 43% of countries poorer than today (SSP3) | 43.56% | Yes |

### Results That Do Not Replicate

| Main-text claim | Replication result | Difference |
|:---|---:|:---|
| Global output loss by 2100: −23% (SSP5, benchmark) | −20.1% | 3 pp |
| Probability of positive global impact: 0.29 | 0.317 | 0.027 |
| Probability losses exceed 20%: 0.44–0.87 | 0.39–0.84 | ~5 pp |
| Most flexible model gives 2.2× larger losses | 2.48× | 0.28× |
| Poorest 40% decline by 75% by 2100 | −80.6% | 6 pp |

### Conclusion

No. The regression coefficients and sample statistics replicate exactly, but several headline projection numbers in the main text do not match the package output exactly.

---

## Question 2: Weighting in the Regressions and Alternative Weighting Scheme

### Where Weighting Plays a Role

BHM uses weighting at three distinct stages:

1. **Data construction**: Country climate variables are already population-weighted (`UDel_temp_popweight`, `UDel_precip_popweight`). Each grid cell is weighted by its resident population before being averaged to the national level.
2. **Main regressions**: Baseline panel regressions are unweighted country-year OLS (clustered by country). This identifies the effect for the "average country-year observation," not the average person:
$$g_{it} = \beta_1 T_{it} + \beta_2 T_{it}^2 + \gamma_1 P_{it} + \gamma_2 P_{it}^2 + \mu_i + \lambda_t + \text{country trends} + u_{it}$$
where $g_{it}$ is GDP per-capita growth.
3. **Projection aggregation**: `ComputeMainProjections.R` computes global average GDP/capita impacts using `weighted.mean(..., population)`.

### Alternative Weighting Scheme: Population Weights on Regression Observations

**Rationale**: If the goal is to estimate the temperature–growth relationship as experienced by the average person globally, each country-year observation should be weighted by population. This is internally consistent: BHM already uses population weights at Stage 1 (data construction) and Stage 3 (projection aggregation), but not at Stage 2. Population weighting at Stage 2 also prevents micro-nations at climate extremes from driving the response function. The weight used is `Pop` — contemporaneous national population from `GrowthClimateDataset.csv`, applied via `feols(..., weights = ~Pop)`.

### Results

| | Unweighted (BHM baseline) | Population-weighted |
|:---|---:|---:|
| β₁ (temp) | 0.01272 | 0.01361 |
| SE(β₁) | 0.00374 | 0.00858 |
| *t*(β₁) | 3.40\*\*\* | 1.59 (n.s.) |
| β₂ (temp²) | −0.000487 | −0.000499 |
| SE(β₂) | 0.000117 | 0.000237 |
| *t*(β₂) | −4.17\*\*\* | −2.11\*\* |
| Optimal T | 13.1°C | 13.6°C |
| ME at 20°C | −0.677%/°C | −0.637%/°C |

**Conclusion**: The shape remains an inverted-U; weighted estimates are close in magnitude, with a slightly higher optimum and a slightly less negative high-temperature effect.

---

## Question 3: Inference and Conley Spatial HAC (200–2,000 km)

### How BHM Handles Inference

The baseline regression is:

$$\Delta y_{it} = \beta_1 T_{it} + \beta_2 T^2_{it} + \beta_3 P_{it} + \beta_4 P^2_{it} + \alpha_i + \gamma_t + \delta_{it} + \varepsilon_{it}$$

where $\alpha_i$ are country fixed effects, $\gamma_t$ year fixed effects, and $\delta_{it}$ country-specific quadratic time trends.

Standard errors are clustered by country (`cluster = ~iso_id`), allowing arbitrary serial correlation within each country but assuming zero cross-country correlation conditional on the fixed effects. This assumption is plausible for within-country dynamics but misses spatial correlation. Year fixed effects absorb the common global component, but regional co-movement (e.g., ENSO) creates cross-country residual correlation that clustering ignores.

### Conley (2008) Spatial HAC Estimator

Conley (1999, 2008) extends HAC variance estimation to two-dimensional (space × time) dependence:

$$\hat{V} = (X'X)^{-1} \left( \sum_{i,t} \sum_{j,s} k\!\left(\frac{d_{ij}}{D}\right) k\!\left(\frac{|t-s|}{L}\right) x_{it}\hat{u}_{it}\, x_{js}'\hat{u}_{js} \right) (X'X)^{-1}$$

where $k(\cdot)$ is the Bartlett kernel, $d_{ij}$ is great-circle distance, $D$ is the spatial cutoff, and $L$ is the temporal lag cutoff.

### Results

Baseline (country-clustered) SEs for comparison: β₁ SE = 0.003737 (*t* = 3.40), β₂ SE = 0.000117 (*t* = −4.17).

| Cutoff (km) | SE(β₁) | \|*t*(β₁)\| | SE(β₂) | \|*t*(β₂)\| |
|---:|---:|---:|---:|---:|
| Baseline | 0.003737 | 3.40 | 0.000117 | 4.17 |
| 200 | 0.003620 | 3.51 | 0.000113 | 4.30 |
| 300 | 0.003683 | 3.45 | 0.000114 | 4.26 |
| 400 | 0.003748 | 3.39 | 0.000115 | 4.22 |
| 500 | 0.003796 | 3.35 | 0.000116 | 4.20 |
| 600 | 0.003837 | 3.31 | 0.000117 | 4.17 |
| 700 | 0.003871 | 3.29 | 0.000117 | 4.15 |
| 800 | 0.003895 | 3.27 | 0.000118 | 4.14 |
| 900 | 0.003913 | 3.25 | 0.000118 | 4.12 |
| 1,000 | 0.003925 | 3.24 | 0.000119 | 4.11 |
| 1,100 | 0.003934 | 3.23 | 0.000119 | 4.09 |
| 1,200 | 0.003941 | 3.23 | 0.000119 | 4.08 |
| 1,300 | 0.003946 | 3.22 | 0.000120 | 4.07 |
| 1,400 | 0.003947 | 3.22 | 0.000120 | 4.06 |
| 1,500 | 0.003947 | 3.22 | 0.000120 | 4.05 |
| 1,600 | 0.003946 | 3.22 | 0.000121 | 4.04 |
| 1,700 | 0.003946 | 3.22 | 0.000121 | 4.03 |
| 1,800 | 0.003945 | 3.22 | 0.000121 | 4.02 |
| 1,900 | 0.003945 | 3.22 | 0.000122 | 4.01 |
| 2,000 | 0.003945 | 3.22 | 0.000122 | 4.00 |

### Do the Estimates Ever Become Insignificant?

No. Neither β₁ nor β₂ ever falls below the 5% significance threshold (\|*t*\| = 1.96) or even the 10% threshold (\|*t*\| = 1.645) across the entire range from 200 km to 2,000 km.

- **β₁ (temp, linear)**: \|*t*\| decreases from 3.51 (200 km) to a floor of 3.22 (1,400–2,000 km) and then stabilizes. The minimum \|*t*\| = 3.22 leaves a margin of 1.26 above the 5% critical value.
- **β₂ (temp²)**: \|*t*\| decreases from 4.30 (200 km) to 4.00 (2,000 km). The minimum \|*t*\| = 4.00 leaves an even larger margin above 1.96.

---

## Question 4: Econometric Issues with BHM's Quadratic Specification and Three Alternatives

### Econometric Issues with BHM's Quadratic Specification

#### McIntosh & Schlenker (2006): Within vs. Global Non-linearity

In a fixed-effects model with a quadratic term, the standard formulation squares the covariate first, then absorbs fixed effects (demeaning). This is not the same as measuring within-group non-linearity:

$$\overline{T^2_{it}} \neq (\bar{T}_{it})^2$$

MS (2006) distinguish two conceptually different non-linearities:

1. **Global non-linearity** (BHM's $T^2_{it}$): the marginal effect depends on the absolute temperature level.
2. **Within non-linearity** (MS's $(T_{it} - \frac{1}{T}\sum_{t=1}^{T} T_{i,t})^2$): the marginal effect depends on deviation from the country's own mean.

#### Nath, Ramey & Klenow (2024): Three Additional Concerns

**1. Serial correlation of temperature as a treatment variable**: Decomposing $T_{it} = \bar{T}_i + \bar{T}_t + \tau_{it}$ (country mean + common year effect + shock), BHM's regression conflates the contemporaneous effect of $\tau_{it}$ with the distributed-lag effects of past shocks that persist in $T_{it}$. When temperature is persistent over time, regressing GDP growth on contemporaneous temperature mixes the impact of a current temperature shock with the lingering effects of earlier shocks.

**2. The quadratic identification problem**: Substituting into the quadratic shows that β₂ in BHM is identified partly not only from within-country variation in temperature, but also from variation in the average temperature across countries $\bar{T}_i$.

**3. Omitted lags**: If temperature shocks persist and GDP growth responds with lags, omitting lagged temperature can make transient level effects appear as permanent growth-rate effects.

### Three Alternative Specifications

#### (a1) NRK Simple Shock Estimator

**Specification:**

$$\Delta y_{it} = \theta_1 \tau^{\text{simple}}_{it} + \theta_2 \left(\tau^{\text{simple}}_{it} \times \bar{T}^{30\text{yr}}_i\right) + X_{it} + \eta_{it}$$

where $\tau^{\text{simple}}_{it} = T_{it} - \bar{T}^{30\text{yr}}_i$ is the deviation from the country's average temperature over the first 30 years of the sample, and $\bar{T}^{30\text{yr}}_i$ is the country mean temperature.

| Parameter | Estimate | SE | *t* |
|:---|---:|---:|---:|
| θ₁ (shock) | 0.01232 | 0.00361 | 3.42\*\*\* |
| θ₂ (shock × T̄) | −0.000961 | 0.000230 | −4.17\*\*\* |
| Bliss T̄\* | 12.8°C | | |
| ME at T̄ = 20°C | −0.0069%/°C | | |
| N | 6,584 | | |

#### (a2) NRK Sophisticated Shock Estimator (AR Innovation)

**First-stage AR(3) model:**

$$T_{it} = \sum_{j=1}^{3} \gamma_j T_{i,t-j} + \sum_{j=1}^{3} \theta_j (T_{i,t-j} \times \bar{T}_i) + \mu_i + \mu_t + \tau_{it}$$

fitted with country and year FEs. The lagged interactions allow the AR dynamics to vary with mean temperature. The innovation to temperature $\tau_{it}$ is then used in the same state-dependent GDP regression.

**Second-stage GDP regression:** Same state-dependent form as above.

| Parameter | Estimate | SE | *t* |
|:---|---:|---:|---:|
| γ₁ (shock) | 0.00867 | 0.00304 | 2.85\*\*\* |
| θ₂ (shock × T̄) | −0.000793 | 0.000219 | −3.63\*\*\* |
| Bliss T̄\* | 10.9°C | | |
| ME at T̄ = 20°C | −0.0072%/°C | | |
| N | 6,277 | | |

#### (b) McIntosh–Schlenker Within Estimator

**Specification:** Replace $T^2_{it}$ with $(T_{it} - \frac{1}{T}\sum_{t=1}^{T} T_{i,t})^2$ so non-linearity is identified purely from within-country variation:

$$\Delta y_{it} = \beta_1 T_{it} + \beta_3 \left(T_{it} - \frac{1}{T}\sum_{t=1}^{T} T_{i,t}\right)^2 + X_{it} + \eta_{it}$$

| Parameter | Estimate | SE | *t* |
|:---|---:|---:|---:|
| β₁ (temp) | −0.000580 | 0.002029 | −0.29 (n.s.) |
| β₃ (within²) | −0.003909 | 0.002145 | −1.82\* (10%) |
| N | 6,584 | | |

#### (c) McIntosh–Schlenker Hybrid Estimator

**Specification:** Include both the global and within quadratics:

$$\Delta y_{it} = \beta_1 T_{it} + \beta_2 T^2_{it} + \beta_3 \left(T_{it} - \frac{1}{T}\sum_{t=1}^{T} T_{i,t}\right)^2 + X_{it} + \eta_{it}$$

This is consistent under any DGP regardless of which form of non-linearity is true.

| Parameter | Estimate | SE | *t* |
|:---|---:|---:|---:|
| β₁ (temp) | 0.01178 | 0.003481 | 3.38\*\*\* |
| β₂ (temp²) | −0.000456 | 0.000113 | −4.05\*\*\* |
| β₃ (within²) | −0.002529 | 0.002104 | −1.20 (n.s.) |
| Optimal T | 12.9°C | | |
| N | 6,584 | | |

### Summary and Discussion

#### Summary Table

| Model | β₁/θ₁/γ₁ | *t* | β₂/β₃/θ₂ | *t* | Bliss/Opt T | ME @ 20°C |
|:---|---:|---:|---:|---:|---:|---:|
| BHM baseline | 0.01272 | 3.40\*\*\* | −0.000487 | −4.17\*\*\* | 13.1°C | −0.68%/°C |
| (a1) NRK simple shock | 0.01232 | 3.42\*\*\* | −0.000961 | −4.17\*\*\* | 12.8°C | −0.69%/°C |
| (a2) NRK AR shock | 0.00867 | 2.85\*\*\* | −0.000793 | −3.63\*\*\* | 10.9°C | −0.72%/°C |
| (b) MS within | −0.00058 | −0.29 | −0.003909 (β₃) | −1.82\* | — | — |
| (c) MS hybrid | 0.01178 | 3.38\*\*\* | −0.000456 (β₂), −0.002529 (β₃) | −4.05\*\*\*, −1.20 | 12.9°C | −0.64%/°C |

#### Discussion

- **BHM baseline**: Strong inverted-U with β₁ = 0.01272 (*t* = 3.40) and β₂ = −0.000487 (*t* = −4.17), implying an optimum of 13.1°C and ME at 20°C of −0.68%/°C.

- **(a1) NRK simple shock**: Essentially identical to BHM, with θ₁ = 0.01232 (*t* = 3.42) and θ₂ = −0.000961 (*t* = −4.17), bliss 12.8°C and ME at 20°C of −0.69%/°C, confirming that the NRK reframing is principally conceptual rather than empirical.

- **(a2) NRK AR shock**: Persistence-adjusted estimates attenuate but remain significant, with γ₁ = 0.00867 (*t* = 2.85) and θ₂ = −0.000793 (*t* = −3.63), shifting bliss to 10.9°C and ME at 20°C to −0.72%/°C, consistent with BHM partly capturing lagged effects when temperature is serially correlated.

- **(b) MS within**: The linear term collapses (β₁ = −0.00058, *t* = −0.29) and within-country curvature is only marginal (β₃ = −0.003909, *t* = −1.82), reflecting the algebraic limitation that within-only curvature cannot locate countries on the global temperature axis.

- **(c) MS hybrid**: BHM's global quadratic is stable (β₁ = 0.01178, *t* = 3.38; β₂ = −0.000456, *t* = −4.05) with optimal 12.9°C, while the added within-squared term is not supported (β₃ = −0.002529, *t* = −1.20), confirming that the data support BHM's global nonlinear form with no additional within-country curvature.

---

## Question 5: NRK-Inspired Re-estimation on BHM Data

In order to address the critiques of Nath, Ramey & Klenow (2024), three upgrades have been implemented: (i) AR(3) temperature innovations rather than levels, (ii) distributed lags of the shock (*k* = 0–5), and (iii) a direct comparison of growth-rate and income-level outcomes.

### Model 1: Distributed Lag Shock Regression

$$\Delta y_{it} = \sum_{k=0}^{5}\bigl[\delta_k\,\tau_{i,t-k} + \delta'_k\,(\tau_{i,t-k}\times\bar{T}_i)\bigr] + \beta_P P_{it} + \beta_{P^2} P_{it}^2 + \mu_i(t,t^2) + \gamma_t + \varepsilon_{it}$$

where $\mu_i(t,t^2)$ denotes country-specific quadratic time trends and $\gamma_t$ are year fixed effects. Marginal effect at $\bar{T}$: $\text{ME}_k(\bar{T}) = \delta_k + \delta'_k\,\bar{T}$.

| Lag *k* | $\hat{\delta}_k$ | *t* | $\hat{\delta}'_k$ | *t* | ME cold (5°C) | ME mod (15°C) | ME hot (25°C) |
|:---:|---:|---:|---:|---:|---:|---:|---:|
| 0 | +0.00730 | 2.01\* | −0.000759 | −3.12\*\* | +0.00350 | −0.00409 | −0.01168 |
| 1 | −0.00276 | −0.65 | +0.000274 | 1.04 | −0.00140 | +0.00134 | +0.00408 |
| 2 | −0.00726 | −1.55 | +0.000304 | 1.16 | −0.00574 | −0.00270 | +0.00034 |
| 3 | −0.00541 | −1.39 | −0.000000 | 0.00 | −0.00541 | −0.00541 | −0.00541 |
| 4 | −0.00686 | −1.97\* | +0.000525 | 2.19\* | −0.00423 | +0.00102 | +0.00627 |
| 5 | −0.00411 | −0.85 | +0.000259 | 0.90 | −0.00282 | −0.00023 | +0.00236 |
| Σ | −0.01911 | | | | −0.01611 | −0.01007 | −0.00404 |

#### Key Findings

- At *k* = 0: cold countries ($\bar{T}$ = 5°C) gain +0.0035 pp from a warm year (below optimum → benefits), while hot countries ($\bar{T}$ = 25°C) lose −0.0117 pp, directly replicating BHM's inverted-U from innovations alone.
- Lags 1–5 show reversal: the contemporaneous effect is partially undone in subsequent years. The 6-year cumulative ME shrinks to −0.0040 pp for hot countries (from −0.0117 at *k* = 0) and −0.0101 pp for moderate climates — strong evidence of mean-reversion, not compounding.

### Model 2: Local Projections

$$\sum_{s=0}^{h}\Delta y_{i,t+s} = \beta_h\,\tau_{it} + \beta'_h\,(\tau_{it}\times\bar{T}_i) + \sum_{k=1}^{3}\Delta y_{i,t-k} + P_{it} + P_{it}^2 + \mu_i(t,t^2) + \gamma_t + \varepsilon_{iht}$$

The left-hand side is the cumulative growth sum over $h+1$ periods. State-dependent Impulse Response Function (IRF): $\text{IRF}_h(\bar{T}) = \beta_h + \beta'_h\,\bar{T}$, with SE via delta method. Under a permanent growth effect, $\text{IRF}_h$ should grow proportionally with $h$; under a transient level effect, it plateaus or reverts.

| Horizon *h* | IRF cold (5°C) | *t* | IRF mod (15°C) | *t* | IRF hot (25°C) | *t* |
|:---:|---:|---:|---:|---:|---:|---:|
| 0 | +0.00463 | 2.11\* | −0.00362 | −1.92 | −0.01187 | −3.43\*\*\* |
| 1 | +0.00272 | 0.70 | −0.00189 | −0.82 | −0.00650 | −1.44 |
| 2 | −0.00225 | −0.44 | −0.00344 | −1.11 | −0.00463 | −0.85 |
| 3 | −0.00580 | −0.94 | −0.00829 | −2.34\* | −0.01079 | −1.70 |
| 5 | −0.00727 | −1.27 | −0.00338 | −0.80 | +0.00052 | 0.07 |
| 10 | +0.00111 | 0.19 | −0.00424 | −0.81 | −0.00960 | −0.98 |

#### Key Findings

- **Hot countries** ($\bar{T}$ = 25°C): IRF = −0.0119 pp at *h* = 0 (\*\*\*), decays to insignificance by *h* = 2 and oscillates near zero.
- **Moderate countries** ($\bar{T}$ = 15°C): small negative impact at *h* = 0 (−0.0036, marginal), significant at *h* = 3 (−0.0083\*), then reverts.
- **Cold countries** ($\bar{T}$ = 5°C): significant positive impact at *h* = 0 (+0.0046, *p* < 0.05), consistent with cold countries benefiting from warming; effect fades and eventually turns slightly negative at *h* = 3–5.
- **Level vs. growth test**: Under a permanent growth effect, the *h* = 10 cumulative IRF should be ≈11× the *h* = 0 value. For hot countries the ratio is 0.96/1.19 ≈ 0.81 < 1: after 10 years the total loss is smaller than the first-year impact. This directly mirrors NRK's finding and rules out a compounding growth mechanism.

### Model 3: Levels Regression

$$\ln(\text{GDPpc})_{it} = \delta_1\,\tau_{it} + \delta_2\,(\tau_{it}\times\bar{T}_i) + P_{it} + P_{it}^2 + \alpha_i + \gamma_t + \varepsilon_{it}$$

Country fixed effects $\alpha_i$ absorb permanent income differences; year effects $\gamma_t$ absorb global trends. If temperature shocks permanently affect the *stock* of income, $\delta_1$ and $\delta_2$ should be significant.

| | $\hat{\delta}_1$ (shock) | $\hat{\delta}_2$ (shock × $\bar{T}$) | Implied bliss $\bar{T}$ |
|:---|---:|---:|---:|
| Estimate | 0.01065 | −0.000776 | 13.7°C |
| *t* | 0.64 | −0.46 | — |

#### Key Findings

Neither coefficient is significant, which means no detectable permanent footprint on income stocks. The implied bliss point of 13.7°C aligns with BHM (13°C) and NRK (13.2°C), but the imprecision confirms that any level effect is too transient to show up in income stocks.

### Conclusion

All three models support a transient level effect, consistent with NRK (2024). BHM's inverted-U shape and the state-dependent nonlinearity are confirmed: cold countries benefit and hot countries lose from a 1°C innovation at *h* = 0. But the effect does not compound because lags 1–5 reverse most of the initial impact: LP IRFs become insignificant within 2 years, and income levels are unresponsive. BHM's −23% projection by 2100 likely overstates long-run damages by treating this transient shock as a permanent annual reduction in the growth rate. This aligns with NRK's central critique: damages are real and state-dependent, but their long-run magnitude is closer to a level effect than a permanent growth effect.

---

## Question 6: Extending the Panel to 2022 and Testing Climate Data Sensitivity

### Panel Specifications

| Specification | Climate | GDP | Years | N |
|:---|:---|:---|---:|---:|
| BHM original | UDel | BHM | 1960–2013 | 6,584 |
| Extended CRU | CRU TS | BHM + WB | 1960–2022 | 8,038 |
| Extended ERA5 | ERA5 | BHM + WB | 1960–2022 | 8,038 |
| Hybrid UDel+CRU | UDel (1960–2013) + CRU (2014–2022) | BHM + WB | 1960–2022 | 8,038 |

### Results

| | BHM original | Extended CRU | Extended ERA5 | Hybrid |
|:---|---:|---:|---:|---:|
| β₁ (temp) | 0.01272\*\*\* | 0.01487\*\*\* | 0.01315\*\*\* | 0.01174\*\*\* |
| SE(β₁) | 0.00374 | 0.00395 | 0.00349 | 0.00351 |
| *t*(β₁) | 3.40 | 3.77 | 3.77 | 3.35 |
| β₂ (temp²) | −0.000487\*\*\* | −0.000504\*\*\* | −0.000523\*\*\* | −0.000401\*\*\* |
| SE(β₂) | 0.000117 | 0.000107 | 0.0000994 | 0.000109 |
| *t*(β₂) | −4.17 | −4.70 | −5.26 | −3.67 |
| Optimal T | 13.1°C | 14.8°C | 12.6°C | 14.7°C |
| ME @ 20°C | −0.0068 | −0.0053 | −0.0078 | −0.0043 |

### Are the Extended Results Different from BHM?

The core finding is robust. Extending the sample through 2022 does not overturn BHM's main result. The inverted-U relationship is significant in all panels. Quantitatively, the cross-source range of β₁ (0.0117–0.0149) and β₂ (−0.000401 to −0.000523) is modest relative to the standard errors. The optimal temperature ranges from 12.6°C (ERA5) to 14.8°C (CRU).

### Are the Results Sensitive to the Climate Data Source?

Yes, results are meaningfully sensitive in magnitude across data sources. The main source of variation is in marginal effects at hot temperatures:

- **ERA5** produces the largest-magnitude damage (ME @ 20°C = −0.0078) and the closest optimal T to BHM (12.6°C).
- **CRU TS** produces the weakest damage (ME @ 20°C = −0.0053) and a higher optimal T (14.8°C).
- **Hybrid UDel+CRU** shows the weakest overall non-linearity.
