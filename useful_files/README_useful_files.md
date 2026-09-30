## GWJULIA (a.k.a. GWInference.jl)

This folder contains useful files for the code, they are divided in:
- Waveform files (in the folder WFfiles)
- Power Spectral Densities (in the folder PSDs)

In writing this folder we acknowledge the help of the GWFAST documentation.


### Waveform Files (WFfiles)

These files are used to interpolate the PhenomD quasinormal modes, this can be useful to speed up a bit the calculation. This procedure is used by LAL and GWFAST so we kept it to reach the same results


### Power Spectral Densities (PSDs)

Note that power spectral densities are the square of the amplitude spectral densities(ASD) $PSD = ASD^2$. The code offers both _readASD() and _readPSD() accordingly.

The files can have 2 or 4 columns, the first one is always the freuencies and the last one is always the full PSD. In the case there are 4 columns there are the PSD for high and low frequncies (HF and LF). 
To select different columns from the function _readPSD(..., cols = [1,4]) where the first index is the frequency and the second one is the PSD.

Inside this folder there are the Power Spectral Densities (or ASD) divided in 5 sub-folder:
- curves\_Jan\_2020
- ET\_curves
- CE\_curves
- LVC\_O1O2O3
- observing\_scenarios\_paper


#### curves\_Jan\_2020
ASD curves from [LIGO](https://dcc.ligo.org/LIGO-T1500293/public) 

#### ET\_curves
ET curves (updated 23/3) from ET collaboration [link](https://apps.et-gw.eu/tds/?content=3&r=18213) plus ETD from [link](https://arxiv.org/abs/1012.0908). Note that there is a pdf document inside the folder explaining the curves updated at 23/3. 


#### CE\_curves

Cosmic Explorer ASDs from [*Science-Driven Tunable Design of Cosmic Explorer Detectors*](https://arxiv.org/abs/2201.10668), available at [https://dcc.cosmicexplorer.org/cgi-bin/DocDB/ShowDocument?.submit=Identifier&docid=T2000017&version=](https://dcc.cosmicexplorer.org/cgi-bin/DocDB/ShowDocument?.submit=Identifier&docid=T2000017&version=).

The folder contains ASDs for:

* the baseline 40km detector (```cosmic_explorer```)
* the baseline 20 km detector compact binary tuned (```cosmic_explorer_20km```)
* the 20 km detector tuned for post-merger signals (```cosmic_explorer_20km_pm```)
* the 40 km detector tuned for low-freqency signals (```cosmic_explorer_40km_lf```)


#### LVK_O1O2O3/

The folder contains ASDs for the LIGO and Virgo detectors during their O1, O2 and O3 observing runs, extracted in specific moment from actual data.

Available at [https://dcc.ligo.org/P1800374/public/](https://dcc.ligo.org/P1800374/public/) for O1 and O2, [https://dcc.ligo.org/LIGO-P2000251/public](https://dcc.ligo.org/LIGO-P2000251/public) for O3a, and computed using [PyCBC](https://pycbc.org) around the times indicated in the caption of Fig. 2 of [https://arxiv.org/abs/2111.03606](https://arxiv.org/abs/2111.03606).

#### observing\_scenarios\_paper/

ASDs used for the paper [*Prospects for observing and localizing gravitational-wave transients with Advanced LIGO, Advanced Virgo and KAGRA*, KAGRA Collaboration, LIGO Scientific Collaboration and Virgo Collaboration](https://link.springer.com/article/10.1007/s41114-020-00026-9).

Available at [https://dcc.ligo.org/LIGO-T2000012/public](https://dcc.ligo.org/LIGO-T2000012/public). 

The folder contains ASDs for the Advanced LIGO, Advanced Virgo and KAGRA detectors during the O3, O4 and O5 observing runs.


#### EOS\_table

Tables of masses and tidal deformabilities for different Equation of states (EOS), used in catalog.jl. The tidal deformabilities are calculated using the LALsuite tool and we stick with its EOS name convention.