#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Oct 15 15:26:33 2020

@author: ursho
This script is used to read in csv files and transform then into .mat file.
The routine used in MATLAB seems to be slower for large datasets
"""
import numpy as np
import scipy.io as sio

names = ['Jan','Feb','Mar', 'Apr', 'May', 'Jun', 'Jul', 'Aug', 'Sep', 'Oct', 'Nov', 'Dec']
directory = "/net/kryo/work/ursho/PhD/Projects/Biomes/Scripts/Biomes_PROOCE/Data/00Probabilities/"

for n in range(12):
    #read csv in
    data_gam = np.genfromtxt(directory+'GAM/'+'dfs_month_'+str(n+1)+'_trsp_gam_pa_gridded_gr_bg(sit_ove).csv',delimiter=",", dtype=float, autostrip=True)
    data_glm = np.genfromtxt(directory+'GLM/'+'dfs_month_'+str(n+1)+'_trsp_glm_pa_gridded_gr_bg(sit_ove).csv',delimiter=",", dtype=float, autostrip=True)
    data_rf = np.genfromtxt(directory+'RF/'+'dfs_month_'+str(n+1)+'_trsp_rf_pa_gridded_gr_bg(sit_ove).csv',delimiter=",", dtype=float, autostrip=True)
    
    #save to mat
    sio.savemat(directory+'GAM/'+names[n]+'_GAM.mat', {names[n]+'_GAM': data_gam})
    sio.savemat(directory+'GLM/'+names[n]+'_GLM.mat', {names[n]+'_GLM': data_glm})
    sio.savemat(directory+'RF/'+names[n]+'_RF.mat', {names[n]+'_RF': data_rf})
