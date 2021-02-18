# LigTMap Target Prediction Method for Small Molecules
Computational Biology and Bioinformatics Lab (CBBIO) (https://cbbio.online/LigTMap)
University of Macau

*Developed by*
*Faraz Shaikh, Giotto Tai, Shirley Siu*

shaikh.faraz78@gmail.com; shirleysiu@um.edu.mo

Welcome to the LigTMap target and activity prediction for
small molecules. This method currently support prediction for 
17 target classes including 6000+ protein targets.

Follow the installation instructions for the standalone program setup.
For online access, please visit https://cbbio.online/LigTMap/.

## Software requirements:
Python 2.7, RDKit, Openbabel, MOPAC2016, ODDT

Specifically, our method has been tested with these versions:
- rdkit-2016.03.4
- numpy-1.11.3
- openbabel-3.0.0
- pychem-1.0
- pybel-0.12.2
- scikit-learn-0.19.2
- scipy-1.1.0
- pandas-0.23.4

## INSTALLATION

### 1. MOPAC2016 
MOPAC2016 can be downloaded from http://openmopac.net/
You need to request for a license to use, please go to the homepage
to obtain a license. The license key will be emailed to you within
one day.

In essence, the installation steps are:

Create the directory:
```
% sudo mkdir -p /opt/mopac
% sudo chmod 777 /opt/mopac
```
Copy over the MOPAC executable that obtained after unpacking the
downloaded package:
```
% cp <somewhere>/MOPAC2016.exe /opt/mopac
% chmod +x /opt/mopac/MOPAC2016.exe
```
Add the following line to your .bashrc start-up script:
```  
alias mopac='/opt/mopac/MOPAC2016.exe'
```
Source the start-up script, e.g.
```
% source ~/.bashrc
```  
Install the license key that you have received in your email:
``` 
% /opt/mopac/MOPAC2016.exe <license-key>
```
Test the installation using the given example:
```
% mopac Example_data_set.mop
```
If the run is completed with output at Example_data_set.out,
then your installation is successful!


### 2. Anaconda 

Download and install Anaconda for Python 2.7 from
https://www.anaconda.com/download/.


### 3. Setup environment in your anaconda
```
% conda create -n ligtmap -c rmg rdkit python=2.7 
% conda activate ligtmap
```
From now on, your python should be the one from the ligtmap env.
Check to confirm:
```
% which python
e.g. /<path>/anaconda3/envs/ligtmap/bin/python
```

### 4. Openbabel
Follow http://openbabel.org/wiki/Category:Installation
to install Openbabel that suits your platform.

For example, for MacOS:
```
% conda install -c conda-forge openbabel  
```

### 5. PyChem

Download and install PyChem from
https://code.google.com/archive/p/pychem/downloads
```
% tar cvfz pychem-1.0.tar.gz
% cd pychem-1.0
% python setup.py install
```

### 6. PyBel
```
% python -m pip install pybel
```

### 7. Scikit-learn + Pandas
```
% conda install -c conda-forge scikit-learn=0.19.2
% conda install -c conda-forge pandas=0.23.4 
```

### 8. ODDT

Make sure you have all previous libraries installed with the
correct version before running this:
```
% python -m pip install oddt
```


> #### Troubleshooting 
> 
> In case you meet errors in between and want to remove and reinstall from step 3:
```
% conda env remove --name ligtmap
```

### 9. PSOVina
Download and install `psovina-2.0.tar.gz` from
https://sourceforge.net/projects/psovina/
```
% tar xfz psovina-2.0.tar.gz
% cd psovina-2.0/build/<your-platform>/release
```
modify `Makefile` to suit your system setting, specifically
give the location of the boost, e.g.:
`BASE=/usr/local/opt/boost-1.59.0`
```
% make
```

### 10. MGLTools
Download and install MGLTools of your platform from
http://mgltools.scripps.edu/downloads


### 11. gsplit
Install some GNU utilities via Homebrew:
```
% brew install gsplit
```


### 12. LigTMap 
Download and unpack `ligtmap-0.1.tar.gz`. You can move the 
program directory to anywhere.
```
% tar xfz ligtmap-0.1
% mv ligtmap-0.1 <your-installed-path>
```

### 13. Setting environment variables
Define necessary environment variables in the `.bashrc`
start-up script file:
```
export LIGTMAP=/<your-installed-path>/ligtmap-0.1
export MGLTools=/<your-installed-path>/mgltools-1.5.7rc1/
export PATH=/<your-installed-path>/psovina-2.0/build/mac/release:$PATH
```
Make sure replace each <your-installed-path> with your real paths.

Finally, source the script file.
```
% source ~/.bashrc
```

## HOW TO RUN TARGET PREDICTION

1. Prepare your molecule(s) to be predicted in `input.smi`

```
E.g., our HIV benchmark molecules:
c1ccccc1Oc(ccc2)c(c23)n(c(=O)[nH]3)CC
c1c(C)cc(C)cc1Oc(ccc2)c(c23)n(c(=O)[nH]3)CC
N#Cc(c1)cc(Cl)cc1Oc(ccc2)c(c23)n(c(=O)[nH]3)CC
N#Cc(c1)cc(Cl)cc1Oc(ccc2)c(c23)n(C)c(=O)[nH]3
```

2. Prepare the list of targets in `target.lst`
```
HIV
HCV
```
For a complete list of targets, refer to `$LIGTMAP/target.lst`.

3. Activate the condo environment
```
% condo activate ligtmap 
```

4. Run the prediction
```
% $LIGTMAP/predict
```

The run will generate two directories `Input` and `Output`.
Input stores each molecule SMILES in a separate file: 
>   input_00001, input_00002, ...

Output stores prediction results for each molecule separately
in directories. 

In case you have a previous run, the Input and Output directories
will be renamed to `Input.xxx` and `Output.xxx`.

5. Examine prediction results

In the summary section, target class in which target proteins 
have been identified for the query molecule is marked `Complete`,
Otherwise `Fail`.

For a molecule `Input_xxxxx`, the top-ranked targets sorted by 
the `LigTMapScore` can be found in `Output/Input_xxxxx/IFP_result.csv`.

This file contain 9 columns of data of the relevant targets:
1. PDB   
2. Class   
3. TargetName   
4. LigandName 
5. LigandSimilarityScore
6. BindingSimilarityScore 
7. LigTMapScore  
8. PredictedAffinity 
9. DockingScore

The binding mode PDB of the molecule at the target protein can be 
found in the corresponding directory
`Output/Input_xxxxx/TargetName/Complex`
