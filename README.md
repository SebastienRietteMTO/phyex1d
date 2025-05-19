# phyex1d: 1D model using PHYEX

## Description
The phyex1D python module runs 1D atmospheric simulations using ideal
cases described in a "common format" file. The physics is provided by
the PHYEX package, also used by the Meso-NH and AROME models.

## Usage

### Installation

The installation procedure is as follows:
```
# phyex1d and DEPHY-SCM installation
pip install phyex1d
git clone https://github.com/GdR-DEPHY/DEPHY-SCM.git

# PHYEX compilation
git clone https://github.com/UMR-CNRM/PHYEX.git
. PHYEX/tools/env.sh
cd PHYEX/build/with_ecbuild
./make_ecbuild.sh
cd ../../..
```

### Experiment

The ARMCU case can be run with the following command:
```
. PHYEX/tools/env.sh
phyex1d ./DEPHY-SCM/ARMCU/REF/ARMCU_REF_SCM_driver.nc --exp dt=60 --plot rc.png rc y_var=P
```
And a rc plot is available in the ```rc.png``` file.

## Main references
  - PHYEX:
    - [code](https://github.com/UMR-CNRM/PHYEX)
  - MesoNH:
    - [reference](https://doi.org/10.5194%2Fgmd-11-1929-2018)
    - [documentation](http://mesonh.aero.obs-mip.fr/)
    - [code](https://src.koda.cnrs.fr/mesonh/mesonh-code)
  - AROME:
    - [reference](https://doi.org/10.1175/2010MWR3425.1)
  - Common format:
    - [documentation and cases](https://github.com/GdR-DEPHY/DEPHY-SCM)
