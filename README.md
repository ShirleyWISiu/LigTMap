# LigTMap Target Prediction Method for Small Molecules

Welcome to the LigTMap target and activity prediction for
small molecules. This method currently support prediction for 
17 target classes including 6000+ protein targets. This code includes the main prediction workflow and all data/models so it can be run offline in your working computer. However, for better visualization of the prediction result, our web server is recommended. 

> Visit our online server at https://cbbio.online/LigTMap/

The code is still in its developmental stage. You are welcome to feedback or join us to contribute to making target prediction a truly powerful method for novel drug discovery! -- Shirley Siu

## Software requirements:
Anaconda, RDKit, Openbabel, MOPAC2016, ODDT, PSOVina, MGLTools, and Python libraries.

Specifically, our method has been tested with these versions:
- python 2.7 (from anaconda)
- rdkit-2016.03.4
- numpy-1.11.3
- openbabel-3.0.0
- pychem-1.0
- pybel-0.12.2
- scikit-learn-0.19.2
- scipy-1.1.0
- pandas-0.23.4
- boost-1.59.0

## INSTALLATION

This code was tested on MacOS X 11.2, CentOS 7.6 and 7.8. We will be glad to know if it works also on your platform!

### 1. MOPAC2016 
MOPAC2016 can be downloaded from http://openmopac.net/MOPAC2016.html
You need a license to use, please go to the homepage
to obtain a license. The license key will be emailed to you.

In essence, the installation steps are:

Create the directory:
```
% sudo mkdir -p /opt/mopac
% sudo chmod 777 /opt/mopac
```
Copy over the MOPAC executable and library that are obtained 
after unpacking the downloaded package:
```
% cp <source-path>/MOPAC2016.exe /opt/mopac
% cp <source-path>/libiomp5.so /opt/mopac
% chmod +x /opt/mopac/MOPAC2016.exe
```
Add the following lines to your .bashrc start-up script:
```  
alias mopac='/opt/mopac/MOPAC2016.exe'
export LD_LIBRARY_PATH=/opt/mopac:$LD_LIBRARY_PATH
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

Download and install the latest version of Anaconda from https://www.anaconda.com/download/.
Simply run the `Anaconda3-xxx.sh` file and provide an installation directory, e.g. 
```
% ./Anaconda3-2020.11-Linux-x86_64.sh
...
/home/user/opt/anaconda3
```
> It's good to organize your program files in one central place like `/home/user/opt/`

### 3. Setup a Python 2.7 environment in your anaconda 
```
% conda create -n ligtmap -c rmg rdkit python=2.7 
% conda activate ligtmap
```
After activation, your default `python` interpreter should be the one from the `ligtmap` env.
Check to confirm:
```
% which python
e.g. /home/user/opt/anaconda3/envs/ligtmap/bin/python
```

### 4. Openbabel
Follow http://openbabel.org/wiki/Category:Installation
to install Openbabel that suits your platform.

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
> In case you meet errors in between and want to remove and reinstall from Step 3:
```
% conda env remove --name ligtmap
```

### 9. PSOVina
Download and install `boost-1.59.0.tar.gz` from https://sourceforge.net/projects/boost/files/boost/1.59.0/ if `boost` is not yet in your system.
```
% tar xfz boost_1_59_0.tar.gz
% cd boost_1_59_0
% ./bootstrap.sh --prefix=/home/user/opt/boost-1.59.0
% ./b2 -j 4
% ./b2 install
```
Add boost to the library path in `.bashrc`
```
export LD_LIBRARY_PATH=$HOME/opt/boost-1.59.0/lib:$LD_LIBRARY_PATH
```

Once your boost is in place, download and install `psovina-2.0.tar.gz` from
https://sourceforge.net/projects/psovina/
```
% tar xfz psovina-2.0.tar.gz
% cd psovina-2.0/build/<your-platform>/release
```
Modify `Makefile` to suit your system setting, specifically
give the location of the boost, e.g.:
`BASE=/home/user/opt/boost-1.59.0`
```
% make
% mkdir /home/user/opt/psovina-2.0
% cp psovina psovina_split /home/user/opt/psovina-2.0
```
Make it accessible by adding the location of the compiled `psovina` to the `PATH` in `.bashrc`
```
export PATH=/home/user/opt/psovina-2.0:$PATH

```
### 10. MGLTools
Download and install MGLTools of your platform from
http://mgltools.scripps.edu/downloads

```
% tar xfz mgltools_x86_64Linux2_1.5.6.tar.gz
% mv mgltools_x86_64Linux2_1.5.6 /home/user/opt
% cd /home/user/opt/mgltools_x86_64Linux2_1.5.6
% ./install.sh
```
Following the instructions at the end of the installation to include some variables in your `.bashrc` file.


### 11. gsplit (for MacOS X only)
Install some GNU utilities via Homebrew, especially, we need `gsplit` as an alternative to the darwin `split`.
```
% brew install coreutils  
```


### 12. LigTMap 
Download and unpack `ligtmap-0.1.tar.gz`. You can move the 
program directory to anywhere.
```
% tar xfz ligtmap-0.1
% mv ligtmap-0.1 /home/user/opt
```

### 13. Setting environment variables
Define necessary environment variables in the `.bashrc`
start-up script file:
```
export LIGTMAP=/home/user/opt/ligtmap-0.1
export MGLTools=/home/user/opt/mgltools_x86_64Linux2_1.5.6/
```

Finally, source the script file.
```
% source ~/.bashrc
```

## HOW TO RUN TARGET PREDICTION

1. Prepare your molecule(s) to be predicted in `input.smi`, e.g. our benchmark molecules for HIV. Make sure you don't leave any empty lines in the file:

```
c1ccccc1Oc(ccc2)c(c23)n(c(=O)[nH]3)CC
c1c(C)cc(C)cc1Oc(ccc2)c(c23)n(c(=O)[nH]3)CC
N#Cc(c1)cc(Cl)cc1Oc(ccc2)c(c23)n(c(=O)[nH]3)CC
N#Cc(c1)cc(Cl)cc1Oc(ccc2)c(c23)n(C)c(=O)[nH]3
```

2. Prepare the list of targets in `target.lst`. For a complete list of supported targets, refer to [$LIGTMAP/target.lst](https://github.com/siuwengin/LigTMap/blob/master/target.lst).
```
HIV
HCV
```


3. Activate the condo environment
```
% condo activate ligtmap 
```

4. Run the prediction
```
% $LIGTMAP/predict
```

The run will generate two directories `Input` and `Output`.
Input stores each molecule `SMILES` in a separate file: 
>   input_00001, input_00002, ...

Output stores prediction results for each molecule separately
in directories. 

> In case you have a previous run, the Input and Output directories
will be backuped to `Input.xxx` and `Output.xxx`.

5. Examine prediction results

In the summary section, the target class for which target proteins 
have been identified for the query molecule is marked `Complete`,
Otherwise `Fail`.

For a molecule `Input_xxxxx`, the top-ranked targets sorted by 
the `LigTMapScore` can be found in `Output/Input_xxxxx/IFP_result.csv`.

This file contain 9 columns of data of the identified targets:
1. PDB   
2. Class   
3. TargetName   
4. LigandName 
5. LigandSimilarityScore
6. BindingSimilarityScore 
7. LigTMapScore  
8. PredictedAffinity 
9. DockingScore

The binding mode (PDB) of the molecule at the target protein can be 
found in the corresponding directory
`Output/Input_xxxxx/TargetName/Complex`

## Citation
Our method paper is currently under review:

>Shaikh, Faraz; Tai, Hio Kuan; Desai, Nirali; Siu, Shirley (2020): LigTMap: Ligand and Structure-Based Target Identification and Activity Prediction for Small Molecular Compounds. ChemRxiv. Preprint. https://doi.org/10.26434/chemrxiv.12923474.v2 

## Contact
Developer:
Faraz Shaikh (shaikh.faraz78@gmail.com), Giotto Tai (hiokuantai@gmail.com)

Project PI:
Shirley Siu (siuwengin@gmail.com / shirleysiu@um.edu.mo)


[Computational Biology and Bioinformatics Lab (CBBIO)](https://cbbio.online/LigTMap)

University of Macau

