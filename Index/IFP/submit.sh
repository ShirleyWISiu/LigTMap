for i in {00..42}
do
stdbuf -oL ~/anaconda2/envs/rdkit/bin/python2.7 convert_final.py ${i} &> log/${i} &
done

