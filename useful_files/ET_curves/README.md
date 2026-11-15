# ET Sensitivity Curves

**Important!** These curves represent a snapshot of a work-in-progress, as of January 2026. These are not meant to reflect the final sensitivity of the detector.

This repository contains public data files for different noise contributions to current ET sensitivity, computed using pygwinc (https://git.ligo.org/gwinc/pygwinc). Note that custom models were used for some noises and basic pygwinc will not give 1-to-1 match with curves presented here. See more details below.

The accompanying document describing the noise budget can be found here: TBD

## Short description of the noises

The files in the "Data" folder descibe each sub-noise and the total noise of the LF and HF detector in power spectral density. There are four colum in each file. 
- 1 column is 'frequency (Hz)', 
- 2 column is the 'displacement ASD', not scaled by length and with no detector response included;
- 3 column is the default pyGWINC sensitivity: it includes the high-frequency correction for the signal coming 15 degrees off zenith;
- 4 column is the sensitivity with high-frequency correction for the signal coming 0 degrees off zenith.

**Note** that the sensitivity is computed for the 90 degrees opening angle between the arms. 