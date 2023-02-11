# kineticsform  

## 1. Overvies  
"kineticsform" is a Python library to set up reaction kinetics equations from reaction formulas and the rate constats provided as csv file.  

## 2. Current version and requirements  
- current version = 0.0.1   
- pyhon 3.7, 3.8, 3.9, 3.10  
- numpy >= 1.20.2  
- pandas >= 1.2.4  

## 3. Copyright and license  
Copyright (c) 2023 Mitsuru Ohno  
Released under the BSD-3 license, license that can be found in the LICENSE file.  

## 4. Installation  
Download the directry "src" to a suitable directry on your computer.  
To start quickly, put the sample script "kineticsform_sample_script.ipynb"  on the same directry described above. In this case, Jupyter Notebook was required.  

## 5. Usage  
### 5-1. Prepare csv file 
1. Describe raction formula in csv format. THe first row shoud be header. The first column and teh second column should be reaction ID (named RID) and the rate constant (named k). 
2. The reactants were putted in alternately with the coeficient and the chemical species.  
3. To separate the reactants and the products, two ">" were setted with two columns. 
4. The products were putted in alternately with the coeficient and the chemical species.  
Note the reaction formula should be inputted by the left filling. The example of a csv format was as follows. The coefficient '1'　can be substituted for blank.  

example of a csv format [1]  

    RID, k,,,,,,,,,,  
    1, 0.054,,AcOEt,,OH-,>,>,,AcO-,,EtOH  
    2, 1.4,,AciPr,,OH-,>,>,,AcO-,,iPrOH  
    3, 0.031,,EGAc2,2,OH-,>,>,2,Ac,,EG  

### 5-2. Run kineticsform   
1. Import kineticsform.  
```py
from src import kineticsform as kf
```
2. Read the csv file as Pandas DataFrame.  

3. Run the the function "react2kinetics". The argument of the function is the name of the DataFrame.  
```py
kf.react2kinetic(df)
```
If the function ran successfully, the number of the unique chemical species, the list of unique chemidcal species and the kinetic equations as text form were returned.  


## References
1) https://doi.org/10.1246/nikkashi.1973.1831
