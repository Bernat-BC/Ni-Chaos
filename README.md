# Noise-induced chaos: a conditioned random dynamics perspective

Bernat Bassols-Cornudella<sup>1</sup> and Jeroen S.W. Lamb<sup>1,2,3</sup>

## Affiliations
<sup>1</sup> Department of Mathematics, Imperial College London, London SW7 2AZ, UK

<sup>2</sup> International Research Center for Neurointelligence, The University of Tokyo, Tokyo 113-0033, Japan

<sup>3</sup> Centre for Applied Mathematics and Bioinformatics, Department of Mathematics and Natural Sciences, Gulf University for Science and Technology, Halwally 32093, Kuwait


## Overview

This repository contains the code and data necessary to reproduce the results of the paper titled "Noise-induced chaos: a conditioned random dynamics perspective." The paper's data was obtained using MATLAB 2023b, and this README provides step-by-step instructions to replicate the study.

## Paper Citation

If you use this code or the data for your research, please consider citing the original paper: 

[![DOI](https://zenodo.org/badge/703998819.svg)](https://zenodo.org/doi/10.5281/zenodo.10209879)

https://arxiv.org/abs/2308.07116


## Table of Contents

1. [Requirements](#requirements)
2. [Installation](#installation)
3. [Usage](#usage)
4. [License](#license)
5. [Contact](#contact)
6. [Acknowledgements](#acks)

## Requirements

To reproduce the study, you will need the following:

- MATLAB R2023b (or later)
- GAIO for matlab, see https://github.com/gaioguy/GAIO.git (this is included in `Data-scripts` when cloning this repository)


## Installation

1. Clone this repository to your local machine:

   ```bash
   git clone https://github.com/Bernat-BC/Ni-Chaos.git

## Usage

This repository is structured as follows:

```
.
├── Data-scripts
│   ├── GAIO-master 
│   │   ├── ...
│   │   └── tpmatrix1d.m
│   ├── all_data_eD_eL.m
│   └── measures.m
├── LICENSE
├── Plot-scripts
│   ├── plots_eD_eL.m
│   └── plots_measures.m
├── README.md
└── Workspaces
    ├── data-eD-eL_depth14_ws.mat
    └── data_00125_depth14_ws.mat
```

### Data-scripts

This folder contains the relevant scripts to generate data.

1. GAIO-master was cloned from https://github.com/gaioguy/GAIO.git and the file `tpmatrix1d.m` was incldued afterwards. (This is effectively the same as the file `tpmatrix.m` in the `GAIO-master' folder but outputs a progress bar to keep track of run time.)

2. `all_data_eD_eL.m` can be run to generate all data shown in Figures 3, 4 and 6, when $\varepsilon$ is used as a varying parameter.

3. `measures.m` can be run to produce the data in Figure 5.

### Plot-scripts

This folder contains the relevant scripts to produce plots.

1. `plots_eD_eL.m` produces Figure 3, 4 and 6.
2. `plots_measures.m` produces Figure 5.

These can be run after importing their respective workspace (see below) or after generating the data by running the scripts in `Data-scripts`.

### Workspaces

This folder contains the worksapces with the exact data plotted in the paper. 

1. `data-eD-eL_depth14_ws.mat` contains the data for Figures 3, 4 and 6.
2. `data_00125_depth14_ws.mat` contains the data for Figure 5.

Run the scripts in `Plot-scripts` to visualise these.

## License

This project is licensed under the GNU GENERAL PUBLIC LICENSE version 3 - see the LICENSE.md file for details.

## Contact

Feel free to send me an [email](bernat.bassols-cornudella20@imperial.ac.uk) with any questions or comments.

## Acknowledgements

The authors gratefully acknowledge support from the EPSRC Centre for Doctoral Training in Mathematics of Random Systems: Analysis, Modelling and Simulation (EP/S023925/1). JSWL thanks IRCN (Tokyo) and GUST (Kuwait) for their research support. We would like to thank M. M. Castro, H. Chu, E. Gibson, T. Pereira, M. Rasmussen, Y. Sato, M. Tabaro, G. Tenaglia and D. Turaev for useful discussions.
