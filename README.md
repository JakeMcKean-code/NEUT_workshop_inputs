# NEUT_workshop_inputs
Scripts and macros for the NEUT workshop (第一回NEUT講習会)

## Setup scripts 
There are two setup scripts, one for Sukap and one for KEKCC. Please use the respective script for your cluster.
```
source setup_<cluster>.sh
```

## Run NEUT
```
neutroot2 NEUT_6.1.4_CC_numu.toml <file.root>
```

## Analysis scripts
There are analysis scripts prepared in both C++ and Python. Run by 
```
root -l -q 'quickNEUTanalysis.cpp("<file.root>")'
```
or 
```
python quickNEUTanalysis.py <file.root>
```

## NuHepMC example
There is an example script using NuHepMC output. It must be compiled first, do so by executing the `build_NuHepMC_example.sh` script as 
```
./build_NuHepMC_example.sh
./NuHepMC_example <file.hepmc3>
```
