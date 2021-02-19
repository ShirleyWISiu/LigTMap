#!/bin/bash

infile=${2}/file_list
LOGDIR=${2}/DOCK_LOG
#refer=Psovina_dataset/file_list
count=0

if [ -e ${LOGDIR}/score_${1}.dat ]; then
   rm ${LOGDIR}/score_${1}.dat
fi

# convert binding affinity in kcal/mol (from vina) to binding constant (expt)
# RT=2.494339
# ki=exp(vinascore*4.184/RT)
# pki=-1.0*alog10(ki)

for nameList in `grep -v "\#" $infile | cut -f1 -d\ ` ; do  
   count=`expr $count + 1`

   score=`cat ${LOGDIR}/${nameList}/vina_$1.log | grep -n2 '| (kcal/mol) | rmsd l.b.| rmsd u.b.' | grep '0.000' | awk '{ print $3 }'`
   echo $score >> ${LOGDIR}/score_${1}.dat
done

exit 0
