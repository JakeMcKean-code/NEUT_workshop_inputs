# source ROOT on KEKCC (required for NEUT to run)
source /cvmfs/sft.cern.ch/lcg/app/releases/ROOT/6.36.04/x86_64-almalinux9.6-gcc115-opt/bin/thisroot.sh;
# source the setup scripts to get neut-config and includes for NEUT analysis scripts
export NEUT_ROOT="/group/t2k/beam/work/jmckean/neut/VSoN_NEUT_inputs/neut/build/Linux"
export ROOT_INCLUDE_PATH="${NEUT_ROOT}/include${ROOT_INCLUDE_PATH:+:${ROOT_INCLUDE_PATH}}"
export CPLUS_INCLUDE_PATH="${NEUT_ROOT}/include${CPLUS_INCLUDE_PATH:+:${CPLUS_INCLUDE_PATH}}"
source ${NEUT_ROOT}/bin/setup.NEUT.sh
# source and setup NUISANCE 
source /group/t2k/beam/work/jmckean/neut/VSoN_NEUT_inputs/nuisance/build/Linux/setup.sh




