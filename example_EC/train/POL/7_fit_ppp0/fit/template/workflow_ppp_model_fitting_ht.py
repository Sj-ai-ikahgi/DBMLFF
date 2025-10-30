


import os
from os.path import join
import sys

whome="/data/home/xiongrui1/sabaisheng/Shenjie/fitting/EC"
sys.path.append(whome + '/POL/7_fit_ppp0/script')

import ppp_model_fitting
import utils




# Please use absolute path for path_main
path_main = join(os.getcwd())

ind_model_0 = (os.path.basename(os.getcwd()).split(sep='_')[0], int(os.path.basename(os.getcwd()).split(sep='_')[1]) - 1)

path_saving = join(path_main, 'temp_folder')
path_polar_param_input  = whome + '/POL/7_fit_ppp0/parameters/polar_param.input_ppp0'
path_in_neighbore_dir   = whome + '/POL/7_fit_ppp0/parameters/IN.NEIGHBORE_3'
path_movement_type_dir  = whome + '/POL/7_fit_ppp0/parameters'
path_fit_linearmm_input = whome + '/POL/7_fit_ppp0/parameters/fit_linearMM.input'
path_pwmat_mlff         = whome + '/POL/7_fit_ppp0/pwmlff_ppp0'
path_parameters         = join(path_pwmat_mlff, 'parameters_EC_ppp.py')
path_ppp_model          = whome + '/POL/7_fit_ppp0/results/ppp_model'
path_ppp_neighbore      = whome + '/POL/7_fit_ppp0/results/ppp_neighbore'
path_job_array          = whome + '/POL/6_generate_ppp0/job_array.txt'

path_lmod_init = '/opt/ohpc/admin/lmod/lmod/init'
path_log = join(path_main, 'log')
path_python3 = '/share/app/anaconda3/bin/python3'
arg_load_intel = ['load', 'intel/2020']
remove_temp_after_fitting = False
movement_exist = False
# path_existed_movement_dir = '/data/home/sqjiang/fitting/md_src_pol/test/fit/existed_movement'
path_existed_movement_dir = None
weight_exist = False
path_existed_weight_dir = None

# def generate_movement_for_ppp_model_fitting(path_job_array, path_saving, ind_ppp_0, sign):
    # path_1 = '/data/home/xueyefengxiang/project_2/1_20230828_ppp_model/1_tools/auto_ppp_model_fitting/test/1_model_1/temp_folder/movement'
    # path_2 = '/data/home/xueyefengxiang/project_2/1_20230828_ppp_model/1_tools/auto_ppp_model_fitting/test/movement'
    # utils.rm_fr(path_1)
    # utils.cp_r(path_2, path_1)

# ppp_model_fitting.generate_movement_for_ppp_model_fitting = generate_movement_for_ppp_model_fitting

ppp_model_fitting.ppp_model_fitting(
    ind_model_0, 
    path_saving, 
    path_polar_param_input, 
    path_job_array, 
    path_in_neighbore_dir, 
    path_movement_type_dir, 
    path_fit_linearmm_input, 
    path_pwmat_mlff, 
    path_parameters, 
    path_ppp_model, 
    path_ppp_neighbore, 
    path_lmod_init, 
    path_log, 
    path_python3, 
    arg_load_intel, 
    remove_temp_after_fitting=remove_temp_after_fitting, 
    movement_exist=movement_exist, 
    path_existed_movement_dir=path_existed_movement_dir, 
    weight_exist=weight_exist, 
    path_existed_weight_dir=path_existed_weight_dir, 
)