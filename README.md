[![](https://img.shields.io/badge/DOI-10.1021%2Facssynbio.0c00288-blue.svg)](https://doi.org/10.1021/acssynbio.0c00288)

Pre-print: [![](https://img.shields.io/badge/DOI-10.1101%2F2020.05.30.122077-blue.svg)](https://doi.org/10.1101/2020.05.30.122077)

# GZM_InducibleTF
Modeling inducible transcriptional regulators, GEM &amp; ZPM, dose-response.

Code asociated to [Gómez-Schiavon, Dods_et al. (2020; _ACS Synthetic Biology_)](https://doi.org/10.1021/acssynbio.0c00288) and [Dods, Gómez-Schiavon _et al._ (2020; _bioRxiv_)](https://doi.org/10.1101/2020.05.30.122077).

## Main files:
- `FN_SS_MODELNAME.m`: Finding steady state solution the output (![formula](https://render.githubusercontent.com/render/math?math=Y_{ss})) for the specific model (e.g. `FN_SS_SimpleHill.m` for simple hill model; see options below) given a set of biophysical paramers (e.g. ![formula](https://render.githubusercontent.com/render/math?math=\mu_Y,\gamma)), inducer concentration (![formula](https://render.githubusercontent.com/render/math?math=H)), and total transcription factor concentration (![formula](https://render.githubusercontent.com/render/math?math=X_{ss})).
- `FN_FitError.m`: Calculate the sum of square errors between model's steady state (![formula](https://render.githubusercontent.com/render/math?math=Y_{ss})) and observed data (![formula](https://render.githubusercontent.com/render/math?math=D)), given a set of biophysical paramers (e.g. ![formula](https://render.githubusercontent.com/render/math?math=\mu_Y,\gamma)), inducer concentration (![formula](https://render.githubusercontent.com/render/math?math=H)), and total transcription factor concentration (![formula](https://render.githubusercontent.com/render/math?math=X_{ss})).
- `FN_FitMRW.m`: Find the set of parameters that best fit the data.
- `RUN_FitMRW.m`: Determine model, conditions, and parameters to be fitted, and call `FN_FitMRW.m`.


## MODEL: Simple Hill model

`FN_SS_SimpleHill.m`

![formula](https://render.githubusercontent.com/render/math?math=Y_{ss}\=f_{SH}(X_{ss},H)/\gamma)

![formula](https://render.githubusercontent.com/render/math?math=f_{SH}(X_{ss},H)\=\mu_Y(\alpha%2B(1-\alpha)\frac{(X_{ss}H)^n}{(X_{ss}H)^n%2BK_D^n}))

## MODEL: Hill + Basal model

`FN_SS_HillxBasal.m`

![formula](https://render.githubusercontent.com/render/math?math=Y_{ss}\=f_{HB}(X_{ss},H)/\gamma)

![formula](https://render.githubusercontent.com/render/math?math=f_{HB}(X_{ss},H)\=\mu_Y(\alpha%2B(1-\alpha)\frac{(X_{ss}H)^n}{(X_{ss}H)^n%2BK_D^n})%2B{\beta}X_{ss})

## MODEL: Mechanistic model

`FN_SS_Mechanistic.m`

![formula](https://render.githubusercontent.com/render/math?math=Y_{ss}\=f_{M}(X_{ss},H)/\gamma)

![formula](https://render.githubusercontent.com/render/math?math=f_{M}(X_{ss},H)\=\mu_Y(1-(1-\alpha)e^{(-{\beta}X_O(X_{ss},H)\/\mu_Y)}))

![formula](https://render.githubusercontent.com/render/math?math=X_O(X_{ss},H)\=\frac{(X_a(X_{ss},H)%2B{\beta}(X_{ss}-X_a(X_{ss},H)))^n}{(X_a(X_{ss},H)%2B{\beta}(X_{ss}-X_a(X_{ss},H)))^n{%2B}K^n})

![formula](https://render.githubusercontent.com/render/math?math=0\=X_a(X_{ss},H)^2-(H%2BX_{ss}%2BK_X)X_a(X_{ss},H)%2B(X_{ss}H))

![formula](https://render.githubusercontent.com/render/math?math=X_a(X_{ss},H){&lt;}X_{ss})

## MODEL: Expanded Hill model

\* NOTE: Also refered as _allosteric_

`FN_SS_Allosteric.m`

![formula](https://render.githubusercontent.com/render/math?math=Y_{ss}\=f_{EH}(X_{ss},H)/\gamma)

![formula](https://render.githubusercontent.com/render/math?math=f_{EH}(X_{ss},H)\=\mu_Y(\alpha%2B(1-\alpha)\frac{(X_a(X_{ss},H)%2B{\beta}(X_{ss}-X_a(X_{ss},H)))^n}{(X_a(X_{ss},H)%2B{\beta}(X_{ss}-X_a(X_{ss},H)))^n%2BK_D^n}))

![formula](https://render.githubusercontent.com/render/math?math=0\=X_a(X_{ss},H)^2-(H%2BX_{ss}%2BK_X)X_a(X_{ss},H)%2B(X_{ss}H))

![formula](https://render.githubusercontent.com/render/math?math=X_a(X_{ss},H){&lt;}X_{ss})
