# Predictive Regressions with Persistent Regressors: IVX and RA-IVX

**Undergraduate thesis — BSc Economics, Universidad Carlos III de Madrid (2025/2026)**
Enrique Ruiz de Almirón Da Silva · Thesis PDF: [English](thesis/TFG_RA-IVX_Phillips_EN.pdf) · [Español](thesis/TFG_RA-IVX_Phillips_ES.pdf) · R

This repository contains the code and empirical analysis for my undergraduate thesis on predictive
regressions with highly persistent regressors.

The project studies the finite-sample behavior of OLS, IVX and Residual-Augmented IVX (RA-IVX)
inference when the predictor is persistent and its innovations are correlated with the regression
error. Its main contribution is to document an **asymmetry**: RA-IVX's advantage over IVX is
substantially larger when the true predictive coefficient is **negative**. This region had not been
examined for RA-IVX. Demetrescu & Rodrigues (2022) study only positive alternatives, and
Hosseinkouchack & Demetrescu (2021) consider both signs only for IVX-type statistics. The
predictability literature was built around positive-sign predictors such as the dividend yield.

The analysis combines Monte Carlo simulations with an empirical application to the Phillips curve,
where the coefficient is naturally negative, for the Euro Area and the United States.

---

## Research Question

Predictive regressions test whether a variable contains information about the future evolution of
another. Conventional OLS inference becomes unreliable when the predictor is highly persistent and
its innovations are correlated with the regression error. The OLS estimator inherits the bias of the
autoregressive coefficient (Stambaugh bias), and the *t*-statistic follows a non-standard
distribution that depends on an unknown persistence parameter.

The analysis focuses on three questions:

1. How do OLS, IVX and RA-IVX behave in finite samples under persistent and endogenous predictors?
2. Does RA-IVX improve the finite-sample size and power of IVX?
3. Is the relative performance of RA-IVX asymmetric for positive and negative predictive coefficients?

---

## Methodology

The baseline predictive regression is

$$
y_t = \alpha + \beta x_{t-1} + u_t,
$$

where the predictor follows a highly persistent AR(1) process,

$$
x_t = \rho x_{t-1} + v_t, \qquad \rho = 1 + c/T,
$$

and the innovations $u_t$ and $v_t$ are correlated.

The project compares three inference procedures:

- **OLS:** conventional benchmark. It is subject to substantial size distortions when the predictor
  is persistent and endogenous.
- **IVX** (Phillips & Magdalinos, 2009; Kostakis, Magdalinos & Stamatogiannis, 2015): builds a
  mildly integrated instrument from the changes in the predictor, with
  $\rho_z = 1 - C_z/T^{\eta}$ and $\eta < 1$. The resulting Wald statistic has a standard χ² limit
  whatever the true degree of persistence, without having to estimate $c$.
- **RA-IVX** (Demetrescu & Rodrigues, 2022): augments the IVX regression with the estimated
  innovations of the predictor process, $\hat v_t$. This removes the component of the error that is
  correlated with $v_t$, reducing the asymptotic variance and correcting the finite-sample bias IVX
  inherits.

**Why the sign matters.** RA-IVX works through two channels. The variance reduction is symmetric in
β, but the bias correction has a sign. When β < 0, both channels push the *t*-statistic into the
rejection region; when β > 0, they partly offset each other.

---

## Monte Carlo Design

The data-generating process is

$$
y_t = \beta x_{t-1} + u_t, \qquad x_t = \rho x_{t-1} + v_t, \qquad \operatorname{Corr}(u_t, v_t) = -0.9.
$$

### Empirical Size

Under $H_0: \beta = 0$, rejection frequencies at the 5% level are computed across:

- Sample sizes: $n \in \{50, 100, 200, 500\}$
- Persistence levels: $\rho \in \{0.70, 0.90, 0.95, 0.99\}$
- 5,000 Monte Carlo replications

### Empirical Power

- $\rho = 0.99$, $\operatorname{Corr}(u_t, v_t) = -0.9$
- Sample sizes from 50 to 500
- Predictive coefficients ranging from negative to positive values
- 2,000 Monte Carlo replications

Additional simulations with 5,000 replications examine the distribution of the *t*-statistics under
the null and under positive and negative alternatives.

---

## Monte Carlo Results

### Size

Empirical size at the 5% nominal level, $\rho = 0.99$:

| $n$ | OLS | IVX | RA-IVX |
|---:|---:|---:|---:|
| 50 | 0.219 | 0.066 | 0.034 |
| 100 | 0.199 | 0.059 | 0.035 |
| 200 | 0.158 | 0.058 | 0.028 |
| 500 | 0.102 | 0.059 | 0.031 |

OLS over-rejects severely: about one in five tests rejects a true null in small samples. IVX stays
close to the nominal level. RA-IVX is slightly conservative, consistent with the behavior reported
by Demetrescu & Rodrigues (2022).

### Power and the Sign Asymmetry

When the true coefficient is positive, the gains from residual augmentation are partly offset by a
shift in the location of the *t*-statistic distribution, and IVX and RA-IVX display similar power.
For negative coefficients, the RA-IVX *t*-statistic shifts further into the rejection region while
also benefiting from the lower variance. The result is a much larger power gain over IVX.

Mean simulated *t*-statistics with $\rho = 0.99$, $\operatorname{Corr}(u_t, v_t) = -0.9$ and
$|\beta| = 0.075$:

| True coefficient | IVX mean *t*-stat | RA-IVX mean *t*-stat |
|---|---:|---:|
| $\beta = +0.075$ | 2.65 | 2.52 |
| $\beta = -0.075$ | −1.15 | −1.92 |

The advantage is largest in small samples and fades as $n$ grows, with the power curves overlapping
around $n \approx 300$.

<!-- If these are in the final thesis, keep them; otherwise delete:
- **Fixed-c robustness check.** With ρ fixed, c grows with n, so sample size and effective persistence
  move together. Repeating the power exercise with c = 1 fixed shows the convergence persists, so it
  is a genuine sample-size effect.
- **Decomposition of the t-statistic.** E[t] = E[β̂]·E[1/SE] + Cov(β̂, 1/SE). For β < 0 the covariance
  term accounts for ~68% of the total under IVX but at most ~24% under RA-IVX.
-->

---

## Empirical Application: Phillips Curve

The predictive specification relates future inflation to current unemployment,

$$
\pi_{t+1} = \alpha + \beta u_t + \varepsilon_{t+1},
$$

where unemployment is treated as a persistent predictor. The analysis covers the Euro Area and the
United States, at monthly and quarterly frequencies, across several subsamples. Before estimation,
persistence and endogeneity diagnostics ($\hat\rho$, ADF / ADF-GLS, $\hat\delta$) locate each
sample within the Monte Carlo design.

**Main results (monthly, *t*-statistics)**

| Panel | Sample | β̂ IVX | *t* OLS | *t* IVX | *t* RA-IVX |
|---|---|---:|---:|---:|---:|
| Euro Area | Full sample | −0.042 | −5.05 | −3.72 | −3.71 |
| Euro Area | Pre-2008 | −0.036 | −1.36 | −0.94 | −0.92 |
| Euro Area | Post-2008 | −0.042 | −4.23 | −4.58 | −4.58 |
| United States | Full sample | 0.008 | 0.49 | 1.41 | 1.05 |
| United States | Pre-Volcker | 0.045 | 1.50 | 3.95 | 4.12 |
| United States | Volcker-Greenspan | 0.030 | 2.54 | 2.53 | 2.58 |
| United States | Great Moderation | −0.024 | −0.16 | −0.76 | −0.75 |
| United States | Bernanke/Yellen/Powell | −0.012 | −1.31 | −1.44 | −1.97 |

Quarterly results, which point the same way, are in `Results/Tables/Main/phillips_curve_results.csv`.

- **Euro Area:** unemployment has a negative and significant effect on future inflation under all
  three estimators. The predictability is a **post-2008** phenomenon; before 2008 there is no
  evidence.
- **United States:** the full sample shows no predictability. It mixes regimes with opposite
  slopes. The pre-Volcker period shows a **positive** slope, a supply-side regime, that OLS misses
  but IVX and RA-IVX detect ($t \approx 4$). The Great Moderation shows a flat curve. In the
  Bernanke/Yellen/Powell period, the negative relationship reappears. Only RA-IVX, the estimator
  that gains most power against negative alternatives in the Monte Carlo, reaches the 5% threshold
  ($t = -1.97$, against $-1.44$ for IVX).
- Empirical endogeneity is much lower than in the Monte Carlo design. When OLS and RA-IVX give
  similar estimates in the data, this reflects low-bias conditions, not OLS reliability. The
  simulations serve as a worst-case benchmark that supports this reading.

---

## Repository Structure

```text
predictive-regressions/
│
├── R/
│   └── ivx_functions.R
│
├── Scripts/
│   ├── 01_build_ea_data.R
│   ├── 02_build_ea_data_quarterly.R
│   ├── 03_build_us_data.R
│   ├── 04_build_us_data_quarterly.R
│   ├── 05_run_phillips_ivx.R
│   ├── 06_diagnostics.R
│   ├── 07_mc_size.R
│   ├── 08_mc_power.R
│   └── 09_run_xdlasso.R        ← high-dimensional extension (see XDLASSO.md)
│
├── Data/
│   ├── EA/
│   └── US/
│
├── Results/
│   ├── Tables/
│   │   ├── Main/
│   │   ├── Diagnostics/
│   │   └── Monte_Carlo/
│   │
│   └── Figures/
│       ├── Diagnostics/
│       └── Monte_Carlo/
│
├── thesis/
│   ├── TFG_RA-IVX_Phillips_EN.pdf
│   └── TFG_RA-IVX_Phillips_ES.pdf
│
├── .gitignore
├── LICENSE
├── predictive-regressions.Rproj
├── README.md
└── XDLASSO.md
```

## Code Organization

- `R/ivx_functions.R`: core OLS, IVX and RA-IVX estimation and inference functions.
- `01–04`: build the Euro Area and US datasets at monthly and quarterly frequencies.
- `05_run_phillips_ivx.R`: main Phillips-curve estimations.
- `06_diagnostics.R`: persistence and endogeneity diagnostics.
- `07_mc_size.R`: Monte Carlo size experiment.
- `08_mc_power.R`: power and *t*-statistic distribution experiments.
- `09_run_xdlasso.R`: high-dimensional extension with XDlasso (see [`XDLASSO.md`](XDLASSO.md)).

In the code, the predictive-regression error $u_t$ is named `eps` and the predictor innovation
$v_t$ is named `nu`.

## Main Outputs

```text
Results/Tables/Main/phillips_curve_results.csv
Results/Tables/Diagnostics/ivx_diagnostics.csv
Results/Tables/Diagnostics/persistence_diagnostics.csv
Results/Tables/Monte_Carlo/mc_size_results.csv
Results/Tables/Monte_Carlo/mc_power_results.csv
Results/Figures/Monte_Carlo/mc_size_plot.pdf
Results/Figures/Monte_Carlo/mc_power_rho99.pdf
Results/Figures/Monte_Carlo/mc_tdist_rho99.pdf
```

---

## Reproducing the Analysis

Open `predictive-regressions.Rproj` in RStudio. Scripts are numbered in execution order: `01`–`06`
reproduce the empirical analysis, and `07`–`08` run the Monte Carlo experiments, which take longer.
All paths are relative to the project root.

## Data

| Region | Source | Inflation | Unemployment |
|---|---|---|---|
| United States | FRED-MD (McCracken & Ng, 2016) | `100·Δlog(CPIAUCSL)` | `UNRATE` |
| Euro Area | EA-MD-QD (Barigozzi, Lissona & Tonni, 2024) | `100·Δlog(HICP)` | EA unemployment rate |

The data belong to their original sources and are included only to ease replication.

## References

- Demetrescu, M. & Rodrigues, P. M. M. (2022). Residual-augmented IVX predictive regression. *Journal of Econometrics*, 227(2).
- Hosseinkouchack, M. & Demetrescu, M. (2021). Finite-sample size control of IVX-based tests in predictive regressions. *Econometric Theory*, 37(4).
- Kostakis, A., Magdalinos, T. & Stamatogiannis, M. P. (2015). Robust econometric inference for stock return predictability. *Review of Financial Studies*, 28(5).
- Phillips, P. C. B. & Magdalinos, T. (2009). Econometric inference in the vicinity of unity. Working paper, Singapore Management University.
- Stambaugh, R. F. (1999). Predictive regressions. *Journal of Financial Economics*, 54(3).
- McCracken, M. W. & Ng, S. (2016). FRED-MD: A monthly database for macroeconomic research. *Journal of Business & Economic Statistics*, 34(4).
- Barigozzi, M., Lissona, C. & Tonni, L. (2024). Large datasets for the Euro Area and its member countries and the dynamic effects of the common monetary policy. arXiv:2410.05082.

## Extension: XDlasso

[`XDLASSO.md`](XDLASSO.md) tests whether the Euro Area result survives controlling for the broader
macroeconomic environment. It uses the XDlasso estimator of Gao, Lee, Mei & Shi (2024). Unemployment
keeps its negative sign, with a magnitude close to the univariate estimate, but is no longer
significant once the macroeconomic panel is included.

## License

Code released under the MIT License.
