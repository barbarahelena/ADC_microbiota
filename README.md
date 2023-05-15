# ADC_microbiota

## Paper
The results of these analyses were published [here](https://doi.org/10.3389/fimmu.2021.794519) in January 2021.

## Methods
We included 170 patients from the Amsterdam Dementia Cohort, comprising 33 with AD dementia (66 ± 8 years, 46%F, mini-mental state examination (MMSE) 21[19-24]), 21 with MCI (64 ± 8 years, 43%F, MMSE 27[25-29]) and 116 with SCD (62 ± 8 years, 44%F, MMSE 29[28-30]). Fecal samples were collected and gut microbiome composition was determined using 16S rRNA sequencing. Biomarkers of AD included cerebrospinal fluid (CSF) amyloid-beta 1-42 (amyloid) and phosphorylated tau (p-tau), and MRI visual scores (medial temporal atrophy, global cortical atrophy, white matter hyperintensities). Associations between gut microbiota composition and dichotomized AD biomarkers were assessed with machine learning classification models. The two models with the highest area under the curve (AUC) were selected for logistic regression, to assess associations between the 20 best predicting microbes and the outcome measures from these machine learning models while adjusting for age, sex, BMI, diabetes, medication use, and MMSE.
For further details on the design of the study, please check out the paper (see above).

## Processing 16S data
16S rRNA primers were removed from the sequencing reads using seqtk (v. 1.3). The reads were subsequently processed using dada2 (v.1.18) as follows. After examining read quality profiles, 50 bases were trimmed from the 5’ end of the forward reads, and 60 bases from the 5’ of the reverse reads, respectively. The reads were truncated at the first base with a Q score lower than 4, then quality filtered using 2 maximum expected error for the forward reads and 4 maximum expected errors for the reverse reads, allowing for no ambiguous bases. The filtered reads were used to learn the error rates and to infer Amplicon Sequence Variants (ASVs) separately for the forward and the reverse reads. Forward and reverse ASVs were merged allowing no mismatching bases and requiring a minimum overlap of 20 bases. ASVs shorter than 350 bp, longer than 500bp, and chimeric ASVs were removed. An ASV table was constructed for the remaining ASVs. ASV taxonomy was then assigned using the dada2 ‘assignTaxonomy’ function and the SILVA database (release 138) allowing up to 3 multiple species-level assignments (Quast et al., 2013; Callahan et al., 2016). The ASV table and taxonomy were integrated using the phyloseq R package (v.1.34.0). The ASV table was rarefied to 20000 counts per sample (McMurdie and Holmes, 2013). Of 175 sequenced samples, 5 had insufficient counts (<20000 counts per sample) and were excluded at the rarefaction stage.
See also the Methods section of the paper and Supplement 1.

## Analyses
The analysis consisted of several steps:
- Data cleaning, create 1 phyloseq object with clinical data and count table (`phyloseq_prep.R`). Create dichotomous variables for the classification models (amy+/-, p-tau+/-, etcetera)
- Table 1 (`table.R`) and plots of distribution CSF biomarker levels (presented in Supplement 3)
- Descriptive plots (`ordinationplots.R` and `compositionalplots.R`) which is Figure 2 in the paper, including a compositional plot at genus level, bray-curtis distance and alpha diversity between diagnosis groups. The ggarrange at the end of the `ordinationplots.R` script uses the genus level composition plot of the `compositionalplots.R` script.
- Create input data for 6 XGBoost models: create-input-data-classification scripts. Input data (x and y) are stored in a folder named after the variable.
- Run XGBoost classification models (`XGBoost.py`). You need to set up the conda environment for this, see below for a separate explanation on the machine learning models.
- Process output of classification models (`process_model_results_class.R`): feature importance plots and violin plots of differences. These can be found in supplement 4. (However the AUC plots in Supplement 4 are python plots, resulting from the XGBoost script above.) This script relies on the `function.R` script.
- Explained variance plot (`explainedvariance.R`; Figure 3 in the paper) showing the explained variance distribution from the XGBoost models.
- Logistic regression models (`logregression-new.R`; Figure 4). This script runs regression models for the best predicting microbes and either amy or p-tau positive status. The results of these models are plotted in a forest plot.

## Renv file
Most analyses (except for the XGBoost models) were performed in RStudio (v.2022.7.2.576) using R (v.4.2.1). We used renv and uploaded a lockfile in this repository to reconstruct the renv.

## Machine learning models
### Introduction
We used a machine learning algorithm to assess which gut microbes (ASVs) were most predictive for amyloid status, p-tau status, MTA, GCA, Fazekas and microbleeds on MRI (all dichotomized). All machine learning models used the XGBoost algorithm in a nested cross-validation design. In each iteration, the dataset was randomly split into a test set containing 20% of the subjects and a training set with the remaining 80%. Within the train set, 5-fold cross-validation was performed in order to optimize the model hyperparameters. Two random variables were added to the determinants in each iteration to serve as a benchmark. The resulting model was evaluated on the test set which yielded an area under the receiver-operator curve (AUC) as main model quality metric. In addition, each iteration resulted in a ranked list of metabolites with their relative importance to the prediction, with the first ranked ASV set at 100% and the other ASV’s importance calculated relative to the first. These were recorded for each iteration and were averaged across 200 iterations.

<img src="https://user-images.githubusercontent.com/34349946/220138280-b58d5408-0fcc-4fbd-812c-df4e402925d3.png" width="500" alt = "flowchart machine learning model">

### Installing conda env
For installation of the conda environment, you will need the yaml file in this folder (xgb_mac_env.yaml). First change the name / path of the environment from the yaml file, if needed. You can install the conda environment using:
`conda env create --file xgb_mac_env.yaml`

Then activate the environment before running a model with:
`conda activate xgb`

### Running the models
As soon as the conda environment is activate, you will need the the following items for running the model:
- Input data (x and y) in the right formats as created by the create_input_data scripts
- XGBoost model python script (`XGBoost.py`)
- Parameter grid (`param_grid_tuned.json` for this project)
- Bash script with the commands to run the model (`xgboost-models-commands.sh`)

The XGBoost.py has a --help section with the possible arguments for the function.

### Output of the models
An output folder is created in the same folder as the input folder when running the model.The XGBoost python script outputs the following the results in the output folder:
- aggregated_metrics_classification.txt: aggregated main model metrics (AUC)
- all_model_parameters.txt: all parameters used, partly defined by the parameter grid
- conda_environment.txt: file that can be used to create conda environment
- feature_importance.txt: file with the ranked microbes and their relative importance for the model
- model_results_per_iteration.txt: main model metrics per iteration
- system_info.txt

