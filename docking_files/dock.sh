#!/bin/bash
#
# please remove the log directory before you being running the docking
# rm -rf log/log20141202/
#
EXE=${1}
LOGNAME=${2}
BRUN=${3}
ERUN=${4}
rootpa=${5}

#### predefined working directories ####
LOGDIR=${LOGNAME}/DOCK_LOG

if [ -e ${LOGDIR} ] ; then
  rm -r ${LOGDIR}
fi

if [ ! -e ${LOGDIR} ] ; then
   mkdir -p ${LOGDIR}
fi


for i in `seq ${BRUN} ${ERUN}`   # run how many repeats of docking
do
  ${rootpa}/docking_files/run_vina.sh ${EXE} ${i} ${LOGNAME} ${rootpa}
  sh ${rootpa}/docking_files/collect_score.sh ${i} ${LOGNAME}
#  sh docking_file/collect_rmsd.sh ${i} ${LOGDIR}
#  sh ./collect_runtime.sh ${i} ${LOGDIR} ${DATASETSUFFIX}
#  sh ./collect_gbest.sh ${i} ${LOGDIR} ${DATASETSUFFIX}
#  sh ./collect_pbest.sh ${i} ${LOGDIR} ${DATASETSUFFIX}
# sh ./collect_conditionVector.sh $i ${LOGDIR} $DATASETSUFFIX
#  sh ./collect_convergestep.sh ${i} ${LOGDIR} ${DATASETSUFFIX}
#  sh ./collect_localcall.sh ${i} ${LOGDIR} ${DATASETSUFFIX}
#  sh ./collect_weight.sh ${i} ${LOGDIR} ${DATASETSUFFIX}
done
