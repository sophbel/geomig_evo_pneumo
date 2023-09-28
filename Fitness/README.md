# Analysis for fitness of *Streptococcus pneumoniae* including vaccine status, individual serotypes, and penicillin resistance
## Fitness model
*This analysis was run using R version 4.0.5* <br />
<br />
Install packages ```install.packages(c("RColorBrewer","rstan","binom","stringr","data.table","doParallel","loo","ggplot2"))```<br />
### Model Types
Each top folders is subset by the three different types of models we run: <br />
1)  NVT_PCV7_PCV13: Fitness of all serotypes grouped by whether they are in the vaccine.
2) Serotypes: Fitness of each serotype across provinces in South Africa.
3) AMR_VaxStat: Fitness of penicillin resistance within non-vaccine type and vaccine type serotypes.

### 1: Raw data
```./1_raw_data/```<br />
In this folder, you will find the raw data used, consisting of the metadat for each sequence in this study.  

### 2: Processed data
```./2_processed_data/```<br />
#### NVT_PCV_PCV13 & Serotypes
In these folders you will find data_creator script which will take the raw data and create the input files for the model. 
#### AMR_VaxStat
1) Fitness Formatting script ```1_fitness.formatting.VaccineStatus_AMR_AMRVaxStat.R``` which will create the input files for the datacreator. These will go into the respective folders for penicillin resistance ```AMR``, Vaccine Type and Non-Vaccine type ```VT_NVT```, AMR within each of those groups ```AMR_VTNVT``` <br />
2) Data creator script ```2_data_creator_fitness_model_AMRVT_together.R``` will create the data input file for this model.

### 3: Model
```./3_model/```<br />
You will find the various models used in this study in their respective folders. These will be ```*.stan``` files and once these compile from the next step there will be ```*.rds``` files here. 

### 4: Run Model
```./4_run_model/```<br />
In this folder you will find the scripts to run each model for each group respectively as well as an output folder where your chains will save. <br>
Link to chains to download https://figshare.com/s/6ce649eb497ee9455ca2

### 5: Figures
```./4_figures/```<br />
In this folder you will find the folders for generating the fitness for each group respectively as well as for plotting the fits. 
