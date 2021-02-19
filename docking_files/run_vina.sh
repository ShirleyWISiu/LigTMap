#!/bin/bash
# to run:
#    run_vina.sh [pso|vina] [1-30] [log directory suffix]
   
##EXE_PSOVINA=${PSOVINA}"/psovina"

VINA=psovina ##${EXE_PSOVINA}


RUN=${2}
LOGDIR=${3}/DOCK_LOG
testfull=${4}"/docking_files"
#infile=${testfull}/file_list
infile=${3}/file_list
count=0

for nameList in `grep -v "\#" $infile | cut -f1 -d\ ` ; do  
   count=`expr ${count} + 1`

   mkdir -p ${LOGDIR}/${nameList}
  
   #cd ${testfull}/${nameList}

   if [ -e ${LOGDIR}/${nameList}/${nameList}_ligand_${RUN}.pdbqt ] ; then
	rm ${LOGDIR}/${nameList}/${nameList}_ligand_${RUN}.pdbqt
   fi

   countdock=0
   maxdock=100
   # repeat if the docking expt doesn't produce result
   while [ ! -e ${LOGDIR}/${nameList}/${nameList}_ligand_${RUN}.pdbqt ] && [ "$countdock" -lt "$maxdock" ]; do
      (time -p $VINA --receptor ${testfull}/pdbqt/${nameList}_protein.pdbqt --ligand ${3}/input.pdbqt --config ${testfull}/config/${nameList}_config.txt --num_modes 1  --log ${LOGDIR}/${nameList}/vina_${RUN}.log --cpu 8 --exhaustiveness 8 --out ${LOGDIR}/${nameList}/${nameList}_ligand_${RUN}.pdbqt) >& ${LOGDIR}/${nameList}/run_${RUN}.log ## 2>&1 
      countdock=`expr ${countdock} + 1`
   done 

   if [ "$countdock" -eq "$maxdock" ]; then
      sed -i '/'${nameList}'/d' $infile
   else
      # compare to expt structure
      #utako cut -c-66 ${LOGDIR}/${nameList}/${nameList}_ligand_${RUN}.pdbqt > ${LOGDIR}/${nameList}/${nameList}_ligand_${RUN}.pdb
      obabel ${LOGDIR}/${nameList}/${nameList}_ligand_${RUN}.pdbqt -opdb -O ${LOGDIR}/${nameList}/${nameList}_ligand_${RUN}.pdb 2> /dev/null
   fi
   #echo -e "0\n0" |g_confrms -f1 ${nameList}_ligand.pdb -f2 ${LOGDIR}/${nameList}/${nameList}_ligand_${RUN}.pdb  -nofit -o ${LOGDIR}/${nameList}/${nameList}_ligand_fit_${RUN}.pdb &>   ${LOGDIR}/${nameList}/${nameList}_ligand_fit_${RUN}.log 

   #cd - >& /dev/null
done

exit 0
