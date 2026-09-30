# Phillips Curve in High Dimensions: XDlasso for the Euro Area

High-dimensional extension of my undergraduate thesis on predictive regressions with persistent
regressors (IVX and RA-IVX; BSc Economics, UC3M, 2025/2026).

The thesis finds that unemployment predicts inflation in the Euro Area using a **single** predictor.
This extension asks whether that result survives once we control for the broader macroeconomic
environment. It applies the **XDlasso** estimator of Gao, Lee, Mei & Shi (2024) to the Euro Area
using the EA-MD-QD dataset. Gao et al. apply the method to the US; as a validation step, their US
results are first replicated.

Estimation uses the authors' own R implementation, the `debias_ivx()` function from their package
[`LasForecast`](https://github.com/zhan-gao/LasForecast). The FRED-MD data construction follows
their replication code in [`lasso_inference`](https://github.com/zhan-gao/lasso_inference).

---

## Main Result

> **In the Euro Area, unemployment predicts inflation on its own, but its predictive content is not
> distinguishable from that of the broader macroeconomic panel.** Once 38–42 macroeconomic controls
> are included, XDlasso finds no significant effect in any subsample. Unlike in the US, however, the
> effect does not vanish: every XDlasso estimate remains negative, and in the full sample it is
> almost identical in magnitude to the univariate estimate of the thesis (−0.041 vs −0.042). What
> changes is the precision. The Euro Area Phillips curve is not absent but **not detectable** once
> the macroeconomic environment is controlled for.

---

## Why a High-Dimensional Test

With a single predictor, the coefficient on unemployment may pick up the effect of omitted variables
correlated with it, such as wages, interest rates, industrial production or commodity prices.
Including them all is not feasible with OLS: there are dozens to over a hundred controls and only a
few hundred monthly observations. Standard high-dimensional methods such as the Lasso are not
designed for a nearly integrated regressor like unemployment.

$$
\pi_{t+1} = \alpha + \beta\, u_t + \mathbf{w}_t'\boldsymbol{\gamma} + \varepsilon_{t+1}
$$

Here $u_t$ is unemployment, the persistent regressor of interest, and $\mathbf{w}_t$ is a
high-dimensional vector of macroeconomic controls.

---

## How XDlasso Works

**1. Standardized Lasso (Slasso).** The Lasso handles more regressors than observations by shrinking
less relevant coefficients toward zero. Persistent regressors, whose sample variance grows with $T$,
coexist here with stationary ones. Under a single penalty the Lasso would favor the persistent ones,
so each regressor is rescaled by its own standard deviation before penalization.

**2. Debiasing.** Shrinkage toward zero is a known bias of the Lasso. It can set to zero coefficients
that are not zero, and the resulting estimator has no standard distribution on which to base a
*t*-test. The debiased Lasso (Zhang & Zhang, 2014) corrects this for the coefficient of interest only,
with a correction term built from an auxiliary residual (the *score*). The score comes from
regressing the variable of interest on the controls.

**3. The IVX instrument in the score.** When the variable of interest is nearly integrated, that
auxiliary regression becomes **spurious** (Granger & Newbold, 1974). The residual inherits the
persistence and the Stambaugh-type bias comes back. Gao et al. build the score from the **IVX
instrument** of unemployment instead:

$$
\tilde z_t = \sum_{j=1}^{t} \rho_z^{\,t-j}\,\Delta u_j, \qquad \rho_z = 1 - C_z/T^{\tau},
$$

with $C_z = 5$ and $\tau = 0.5$ as in the paper. The instrument remains correlated with unemployment
but is less persistent, so the auxiliary residual is no longer spurious.

The resulting estimator handles many controls, allows standard inference on β, and remains valid
under near-integration without knowing the true degree of persistence.

### Three Estimators Side by Side

| Estimator | Controls | Persistence correction | What it isolates |
|---|---|---|---|
| **IVX** | No | Yes | Predictability of unemployment on its own |
| **Dlasso** | Yes (score built from $u_t$) | No | Controlling without correcting for persistence |
| **XDlasso** | Yes (score built from the IVX instrument) | Yes | The full Gao et al. estimator |

All three are computed with `LasForecast::debias_ivx()`.

---

## Data

| Region | Source | Sample | Inflation | Unemployment | Controls |
|---|---|---|---|---|---|
| Euro Area | EA-MD-QD (Barigozzi, Lissona & Tonni, 2024) | 2000:03–2025:02, monthly | `100·Δlog(HICP)` | total EA unemployment rate | 42 / 38 |
| United States | FRED-MD (McCracken & Ng, 2016) | 1960:01–2025:04, monthly | `100·Δlog(CPIAUCSL)` | `UNRATE` | 110 / 85 |

The Euro Area series are aggregates for the countries that use the euro, not for the whole EU.

Two versions of the control panel are used for each region:

- **Transformed:** each series made stationary according to its transformation code (*tcode*).
- **Untransformed:** series in levels or log-levels, excluding I(2) variables and series without
  enough history.

---

## Step 1 — Validation: Replicating Gao et al. for the US

Following the authors' replication code, FRED-MD is loaded with `fbi::fredmd(transform = TRUE)`,
and every series with missing values up to April 2025 is dropped automatically. This yields exactly
the **110 controls** used in the paper.

**US monthly, transformed panel — this repository vs. Gao et al. (2024), Table 7(a)**

| Period | IVX (ours) | IVX (paper) | Dlasso (ours) | Dlasso (paper) | XDlasso (ours) | XDlasso (paper) |
|---|---:|---:|---:|---:|---:|---:|
| Full sample | −0.014 (0.018) | −0.014 (0.018) | 0.019 (0.006)*** | 0.018 (0.006)*** | −0.021 (0.073) | −0.024 (0.077) |
| Pre-Volcker | 0.077 (0.069) | 0.080 (0.069) | 0.074 (0.018)*** | 0.074 (0.017)*** | 0.002 (0.223) | 0.013 (0.224) |

The replication matches up to small differences attributable to the FRED-MD vintage. With the
untransformed panel the point estimates differ from the paper, but the conclusion is the same:
neither IVX nor XDlasso finds significant predictability.

---

## Step 2 — Results: Euro Area

**Table 1 — Euro Area, monthly, transformed panel (p = 42 controls)**

| Period | T | IVX | Dlasso | XDlasso |
|---|---:|---:|---:|---:|
| Full sample | 299 | −0.120 (0.044)*** | −0.017 (0.007)** | −0.041 (0.034) |
| Pre-2008 | 93 | −0.058 (0.059) | −0.026 (0.026) | −0.055 (0.054) |
| Post-2008 | 205 | −0.154 (0.065)** | −0.016 (0.008)** | −0.055 (0.047) |

**Table 2 — Euro Area, monthly, untransformed panel (p = 38 controls)**

| Period | T | IVX | Dlasso | XDlasso |
|---|---:|---:|---:|---:|
| Full sample | 284 | −0.130 (0.048)*** | +0.069 (0.035)** | −0.100 (0.097) |
| Pre-2008 | 93 | −0.058 (0.059) | −0.097 (0.196) | −0.256 (0.158) |
| Post-2008 | 190 | −0.181 (0.082)** | +0.077 (0.054) | −0.221 (0.238) |

Standard errors in parentheses. \*\*\* p < 0.01, \*\* p < 0.05, \* p < 0.10.

---

## What the Results Say

**1. On its own, unemployment predicts inflation.** IVX is negative and significant in the full
sample and after 2008, and not significant before 2008. This is the same pattern found in the thesis,
where the time pattern is identical. The result survives a change of instrument ($\tau = 0.5$ here,
against the Demetrescu–Rodrigues construction used in the thesis).

**2. Controlling for the macroeconomic panel, the effect is no longer significant.** XDlasso is not
significant in any subsample or panel version. The Euro Area therefore reaches the same conclusion
Gao et al. document for the US: the predictability found with a single regressor does not survive a
high-dimensional control.

**3. But the effect does not disappear: it loses precision.** All six XDlasso estimates are negative,
consistent with the Phillips curve. In the full sample, XDlasso gives −0.041, almost exactly the
thesis estimate (−0.042), but with a standard error three times larger (0.034 vs 0.011). In the US,
by contrast, the XDlasso estimates hover around zero and change sign. The difference between the two
economies is a matter of degree:

| | United States | Euro Area |
|---|---|---|
| Univariate predictability (thesis) | None in the full sample | Negative and significant |
| XDlasso estimate | Around zero, changing sign | Negative in every subsample |
| Reading | No effect | Effect in the expected direction, not detectable |

**4. Without the IVX correction, the debiased Lasso is unreliable.** In the untransformed panel,
Dlasso turns **positive and significant** in the full sample (+0.069), against negative estimates
from IVX (−0.130) and XDlasso (−0.100). The only difference between Dlasso and XDlasso is the IVX
instrument in the debiasing step. With persistent controls in levels, the auxiliary regression
behaves spuriously, as the theory predicts and as Gao et al. find in the US.

---

## Conclusions

1. In the Euro Area, unemployment carries information about future inflation **on its own**, and the
   result holds with a different IVX instrument.
2. That information is **not distinguishable** from the information in the broader macroeconomic
   panel: XDlasso finds no significant effect in any subsample.
3. Unlike in the US, the Euro Area effect keeps the expected negative sign and a magnitude close to
   the univariate estimate. The Phillips curve there is weak and imprecise rather than absent.
4. Correcting for persistence remains essential in high dimensions: the standard debiased Lasso can
   reverse the sign of the effect.

**Limitations.**

- **Sample size.** With T ≈ 300 monthly observations (T = 93 before 2008), the data cannot separate
  two explanations: other variables carrying the same information as unemployment, or insufficient
  power to detect a small effect.
- **Labor-market controls.** The Euro Area panel includes other labor-market variables closely
  related to unemployment. Gao et al. also report results excluding them. For the Euro Area this
  robustness check is only available at quarterly frequency.
- **Quarterly results.** They are computed but not used for inference. The Euro Area quarterly panel
  has more controls than observations (T = 94, p = 114; T = 26 before 2008), a setting where the
  Lasso cannot reliably identify the model. The same happens for the US pre-Volcker quarterly
  sample (T = 78, p = 111).

---

## Code

```text
Scripts/
├── 01_build_ea_data.R             ← EA monthly: IVX objects + transformed / untransformed panels
├── 02_build_ea_data_quarterly.R   ← EA quarterly
├── 03_build_us_data.R             ← FRED-MD via fbi::fredmd(): 110 controls
├── 04_build_us_data_quarterly.R   ← US quarterly
└── 09_run_xdlasso.R               ← main_xdlasso(): IVX / Dlasso / XDlasso by subsample and panel
```

## Reproducing

```r
install.packages(c("dplyr", "tibble", "glmnet", "remotes"))
remotes::install_github("zhan-gao/LasForecast")   # debias_ivx()
remotes::install_github("cykbennie/fbi")          # fredmd()

source("Scripts/09_run_xdlasso.R")                # sources the data-building scripts 01–04
```

## References

- Gao, Z., Lee, J. H., Mei, Z. & Shi, Z. (2024). LASSO inference for high dimensional predictive regressions. arXiv:2409.10030. Code: [`LasForecast`](https://github.com/zhan-gao/LasForecast), [`lasso_inference`](https://github.com/zhan-gao/lasso_inference).
- Zhang, C.-H. & Zhang, S. S. (2014). Confidence intervals for low dimensional parameters in high dimensional linear models. *Journal of the Royal Statistical Society: Series B*, 76(1).
- Granger, C. W. J. & Newbold, P. (1974). Spurious regressions in econometrics. *Journal of Econometrics*, 2(2).
- Kostakis, A., Magdalinos, T. & Stamatogiannis, M. P. (2015). Robust econometric inference for stock return predictability. *Review of Financial Studies*, 28(5).
- Demetrescu, M. & Rodrigues, P. M. M. (2022). Residual-augmented IVX predictive regression. *Journal of Econometrics*, 227(2).
- McCracken, M. W. & Ng, S. (2016). FRED-MD: A monthly database for macroeconomic research. *Journal of Business & Economic Statistics*, 34(4).
- Barigozzi, M., Lissona, C. & Tonni, L. (2024). Large datasets for the Euro Area and its member countries and the dynamic effects of the common monetary policy. arXiv:2410.05082.

## License

Code in this repository is released under the MIT License. `LasForecast` and `lasso_inference` are
the work of their authors and are used under their own licenses.
