# Code documentation for the Bioprospected_LCA repository

## Overview of the repository

This repository contains the code used to reproduce the results of the manuscript: XXXX 
 
### Overview of folders and files:

**Data**

+ **_elemental_contents.csv_** includes average elemental compositions (N, P, C) of microalgal macromolecules (Lipid, Phospholipids, proteins, Carbohydrates) as given 
by Geider et al. 2011.  


+ **_Feed_composition.csv_** includes the fish feed composition (wheat, oil etc.) and the calculated composition in term of biochemical classes (Lipid, proteins, carbohydrates, ash, water).  

+ **_Microalgae_foreground_2.json_** includes the foreground database.


**Climatic_data**

+ After running the code, this folder will contain the csv files with temperature and irradiance data (W.m<sup>-2</sup>) over an average day of each month of the cultivation period in the  locations of the grid modeled in the simulations.

**Plot**

All plots from the scripts are saved in this folder. The plot files process the results saved in the folders Outputs-Mono and Outputs_Multi.

+ **_Plot_Multi-dimensional sampling-Figure3 and SI..R_** Script to process results and produce Figure 3 and Figures S14 to S21.  

+ **_Plot_Multi-dimensional sampling-Figure4_** Script to generate the production mixes and Figure 4.

+ **_Plot-Mono-dimensional sampling_** Script to produce the figures associated to mono-dimensional Sampling : Figure 2 and Figures S5 to S13 in SI I.

+ **_Sensitivity_Plot_** Script to produce Figure5 and Figure S22 and S23.

**Environment**

+ **environment_microalgae_windows.yml** File needed to create the virtual environment on WINDOWS.
+ **environment_microalgae_ubuntu.yml** File needed to create the virtual environment on UBUNTU.


**Outputs_Mono**

All outputs of csv and xlsx types from Simulate_Mono_2nd are saved in this folder.

**Outputs_Multi**

All outputs of csv and xlsx types from Simulate_Multi_2nd are saved in this folder.

**Scripts**

+ Ten **.py** files: python scripts including the model itself and needed to run the simulations. 

Files, scripts, and their functions'interconnections are mapped below.  
<br>  

<img src="Code map_2nd.jpg"
     alt="Markdown Monster icon"
     style="float: left; margin-right: 10px;" />  
<br>  




**Retrieving_solar_and_climatic_data_2nd**

Contains functions which download solar and temperature data from the European Photovoltaic Geographical System and estimate the solar power received by a given PBR geometry during the day.


**Functions_for_physical_and_biological_calculations_2nd**

Contains functions to calculate values related to :

+ Strain biology and metabolism
+ PBR geometry
+ Culture physical properties and centrifugation
+ Optimization for fish feed substitution
+ Anaerobic digestion


**Cultivation_simul_2nd**

Contains functions which simulate the cultivation over a day with solving of differential equations for temperature evolution and associated thermoregulation.

**Main_simulations_functions_Mono_2nd**

Contains the main functions calling the other scripts to calculate the LCAs in mono-dimensional sampling. 


**Main_simulations_functions_Multi_2nd**

Contains the main functions calling the other scripts to calculate the LCAs in multi-dimensional sampling. 


**Map_2nd**

Contains the functions generating and handling the geographic files/locations grid/geodataframes.



**Simulate_Mono_2nd** 

To be executed to conduct the mono-dimensional sampling of the research article and the associated LCAs. Generates the results from mono-dimensional sampling.

Exports results in folder "Output-Mono".

 See *Reproducing results from the article.*

**Simulate_Multi_2nd** 

To be executed to conduct the multi-dimensional sampling and the associated calculations. Generates the results from multi-dimensional sampling with the exception of the production mixes.

Exports results in folder "Output-Multi". 

See *Reproducing results from the article.*



**Prepare_project_2nd** 

Creates the Brightway2 project and loads the foreground database in it. Imports your local version of ecoinvent 3.6 consequential in the new project and loads bioshpere3.

See *Reproducing results from the article.*

<br>

### Reproducing results from the article

*Requirements*

+ Miniconda or Anaconda
https://docs.conda.io/projects/conda/en/latest/user-guide/install/index.html

+ A local version of the ecoinvent 3.6 consequential database

+ A python interface (e.g., Spyder) and a R interface (e.g., Rstudio)

+ An internet connection as the code downloads climatic data from an API.

*Step by step procedure:*

1. **Download or clone a local copy of the repository. Keep the folders structure.**

2. **Prepare a conda environment with all needed packages**

+ From terminal or Anaconda/Miniconda terminal access the "Environment" folder. The path depends on where you saved the repository:

```
 cd <yourpathtothefolder/Environment>
```

+ Create the conda environment with all necessary packages using the .yml file corresponding to your OS.

**For Windows:**


```
conda env create --file environment_microalgae_windows.yml
```

+ Activate the newly created environment:

```
conda activate environment_microalgae_windows
```

+ Install the ray package from pip:
```
pip install ray
```

**For Ubuntu:**


```
conda env create --file environment_microalgae_ubuntu.yml
```

+ Activate the newly created environment:

```
conda activate environment_microalgae_ubuntu
```

+ Install the geopandas package:
```
conda install geopandas 
```

+ Install the shapefile package:
```
pip install pyshp 
```


+ Install the ray package from pip:
```
pip install ray
```

For MACOS, you should try to install the environment from the ubuntu file and, in case of issues, complete the environment by installing the problematic packages manually. 




3. **Set up the Brigtway2 project**

+ In the Scripts directory, open the file **Prepare_project_2nd.py** in a text editor or python interface and change the value of the variable ```ei36dir``` by specifying the directory where the ecoinvent files are on your drive: ```ei36dir = <yourpathtoecoinventfiles>```. 

+ From the python interface or from command, execute the whole script to prepare the Brightway2 project (```python Prepare_project_2nd.py```) .

4. **Run the simulations using the model** 

Mono-dimensional and multi-dimensional samplings are done separately with two different scripts.

+ In the Scripts directory, open the file **Simulate_Mono_2nd** or/and **Simulate_Multi_2nd**. 
  Read the instructions at the top of the file, and, as indicated, change the value of the sample size to quickly test the model .
<br>
  &#128680;&#128680;&#128680;**WARNING 
  These scripts were run on a remote server with 64 cores, 256 GB RAM and it took respectively 3 days and 7 days to produce the results with the chosen   samples'sizes for mono and multi-dimensional samplings.
  A laptop with fewer cores will take much more time.
  To test the code, one can use lower sample sizes.**

+ Run the simulation by executing the whole script from the python interface or from command line (```python Simulate.py```). 

+ Wait for all the simulations to be finished (Takes a few minutes to several days depending on the sample size and your computer). The script will export excel and csv files in the folder "Outputs-Mono" or "Outputs-Multi".

5. **Plot the figures based on the generated excel files **

+ Open the R script **Plot-Mono-dimensional sampling**. Execute to plot and save the figures associated with Mono-dimensional sampling (Article and SI). 
+ Open the R script **Plot_Multi-dimensional sampling-Figure4**. Execute to process production mixes and plot and save Figure 4. 
+ Open the R script **Plot_Multi-dimensional sampling-Figure3 and SI**. Execute to plot and save Figure 3 and associated Figures in SI. 
+ Open the R script **Sensitivity_Plot**. Read the instructions and change the ID of the sensitivity file to match the one created in your folder *Outputs-Multi*. Execute to plot and save Figure 5 and associated figures in SI.
+ Open the R script **Plot_geo_2nd**. Read the instructions and change the ID of the sensitivity file to match the one created in your folder *Outputs-Multi*. Execute to plot and save Figure 5 and associated figures in SI.



<br>  


REFERENCES

Geider, R. & La Roche, J. Redfield revisited: variability of C:N:P in marine microalgae and its biochemical basis, European Journal of Phycology, 37:1, 1-17, (2002)

Perez-Lopez, P. et al. Comparative life cycle assessment of real pilot reactors for microalgae cultivation in different seasons. Appl. Energy 205, 1151-1164 (2017)
