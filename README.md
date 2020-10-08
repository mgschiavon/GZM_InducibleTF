[![](https://img.shields.io/badge/DOI-10.1101%2F2020.05.30.122077-blue.svg)](https://doi.org/10.1101/2020.05.30.122077)

# GZM_InducibleTF
Modeling inducible transcriptional regulators, GEM &amp; ZPM, dose-response.

Code asociated to [Dods, Gómez-Schiavon et al. (2020)](https://doi.org/10.1101/2020.05.30.122077).

## MODEL: Simple Hill model

'FN_SS_SimpleHill.m'

![formula](https://render.githubusercontent.com/render/math?math=Y_{ss}\=f_{SH}(X_{ss},H)/\gamma)

![formula](https://render.githubusercontent.com/render/math?math=f_{SH}(X_{ss},H)\=\mu_Y(\alpha+(1-\alpha)\frac{(X_{ss}H)^n}{(X_{ss}H)^n%2BK_D^n}))

## MODEL: Hill + Basal model

'FN_SS_HillxBasal.m'

## MODEL: Mechanistic model

'FN_SS_Mechanistic.m'

## MODEL: Expanded Hill model

\* NOTE: Also refered as _allosteric_

'FN_SS_Allosteric.m'
