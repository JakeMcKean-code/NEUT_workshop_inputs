# source the setup scripts to get neut-config and includes for analysis scripts
export NEUT_ROOT="/usr/local/root6_cxx17/neut_6.1.4/build/Linux"
export ROOT_INCLUDE_PATH="${NEUT_ROOT}/include${ROOT_INCLUDE_PATH:+:${ROOT_INCLUDE_PATH}}"
export CPLUS_INCLUDE_PATH="${NEUT_ROOT}/include${CPLUS_INCLUDE_PATH:+:${CPLUS_INCLUDE_PATH}}"
source /usr/local/root6_cxx17/setup_neut_6.1.4.sh



