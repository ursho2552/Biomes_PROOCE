%% Analysis of the Environmental parameters of our biomes

% =========================================================================
% This script is a stand-alone script that uses the monthly, seasonal,
% annual biomes, and sometimes the raw presence/absence projections. Thus
% we first make sure that our workspace is empty, and everything is setup
% correctly.
% =========================================================================

clear all
folder_main = '/net/kryo/work/ursho/Damiano_Presence_data/presence_absence_tables_ensemble_averages/Group_specific_background_approach';
addpath(genpath(folder_main))
cd(folder_main)

% =========================================================================
% Load all variables needed
% =========================================================================

%load trained SOM
cd('/net/kryo/work/ursho/Damiano_Presence_data/presence_absence_tables_ensemble_averages/Group_specific_background_approach/Data/01NeuronsError/')
load('Single_run_11.mat')

%load help variables
cd('/net/kryo/work/ursho/Damiano_Presence_data/presence_absence_tables_ensemble_averages/Group_specific_background_approach/Data/')
load('HelpVariables.mat')
load('Area_map.mat')

%construct or load simplified version of raw data
cd('/net/kryo/work/ursho/Damiano_Presence_data/presence_absence_tables_ensemble_averages/Group_specific_background_approach/Data/00Probabilities/')
load('Simple_sort_Data.mat')
load('Transformed_CompleteSuitePhyto.mat')
load('Seasonally_corrected_data.mat')
load('Names_species')

%load biomes
cd('/net/kryo/work/ursho/Damiano_Presence_data/presence_absence_tables_ensemble_averages/Group_specific_background_approach/Data/05Biomes/')
load('No_mean_PCA_biomes_9_v2_5_perc.mat')
load('No_mean_PCA_biomes_annual_9_v2_5_perc.mat')
load('No_mean_PCA_biomes_seasonal_9_v2_5_perc.mat')

%% Get dendrogram of environmental data

GetEnvironmentalData

%% Get box plots

GetBoxPlotsEnv

%% Multicomparison

EnvMulticomparison

%% Get best combination of env variables

GetEnvCombinations

%% Koeppen analysis

Koeppen

%% Calculate niches of species

Niches


