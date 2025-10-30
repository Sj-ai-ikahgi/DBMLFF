import os
from os.path import join

import sys

import importlib

import get_energy_and_force
import generate_movement
import utils


def generate_ind_ppp_0_according_to_polar_param_input(path_polar_param_input, ind_model_0):
    with open(join(path_polar_param_input)) as f:
        natom = int(f.readline().split()[0])
        
        ipol = []
        indp = []
        ind_ppp_model_atom = []
        sign_ppp_model_atom = []
        for i in range(natom):
            content = f.readline().split()
            if (int(content[0]) - 1) != i:
                print('ppp_model_fitting.generate_ind_ppp_0_according_to_polar_param_input: ERROR!!! (int(content[0]) - 1) != i, EXIT!')
                exit()
            else:
                ipol.append([int(j) for j in content[1:1+4]])
            
            indp.append([int(j) for j in f.readline().split()[0:0+6]])

            content = f.readline().split()
            ind_ppp_model_atom.append([int(j) for j in content[0:0+6]])
            sign_ppp_model_atom.append([int(j) for j in content[0+6:0+6+6]])
            
        numb = int(f.readline().split()[0])
        
        indb = []
        ind_ppp_model_bond = []
        sign_ppp_model_bond = []
        for i in range(numb):
            content = f.readline().split()
            indb.append([int(j) for j in content[0:0+3]])
            ind_ppp_model_bond.append([int(j) for j in content[0+3:0+3+4]])
            sign_ppp_model_bond.append([int(j) for j in content[0+3+4:0+3+4+4]])
    
    list_ind_ppp_0 = []
    list_sign = []
    if ind_model_0[0] == 'atom':
        ind_ppp_model = ind_ppp_model_atom
        sign_ppp_model = sign_ppp_model_atom
    elif ind_model_0[0] == 'bond':
        ind_ppp_model = ind_ppp_model_bond
        sign_ppp_model = sign_ppp_model_bond
    else:
        raise ValueError("ind_model_0[0] should be either 'atom' or 'bond'. However, ind_model_0[0] == %s" % (ind_model_0[0]))
    
    for i in range(len(ind_ppp_model)):
        for j in range(len(ind_ppp_model[i])):
            if (ind_ppp_model[i][j] - 1) == ind_model_0[1]:
                list_ind_ppp_0.append((ind_model_0[0], i, j//2, j%2))
                list_sign.append(sign_ppp_model[i][j])
    
    return list_ind_ppp_0, list_sign, ipol, indb


def generate_movement_for_ppp_model_fitting(
    path_job_array, 
    path_saving, 
    ind_ppp_0, 
    sign, 
    movement_exist=False, 
    path_existed_movement_dir=None, 
):
    movement_name = 'MOVEMENT_%s_%d_%d_%d' % tuple(ind_ppp_0)
    
    if (movement_exist == True) and (path_existed_movement_dir != None) and (os.path.isfile(join(path_existed_movement_dir, movement_name))):
        os.symlink(join(path_existed_movement_dir, movement_name), join(path_saving))
    else:
        with open(join(path_job_array)) as f:
            job_array = f.read().splitlines()
        
        with open(join(path_saving), mode='w') as f_write:
            for i in job_array:
                with open(join(i, 'atom.config')) as f:
                    atom_config = f.readlines()
                
                num_atom = int(atom_config[0].split()[0])
                
                lattice = atom_config[2:2+3]
                
                position = atom_config[6:6+num_atom]
                element = [int(j.split()[0]) for j in position]
                
                force = [[element[j], 0.5, 0.5, 0.5] for j in range(num_atom)]
                
                ppp = get_energy_and_force.get_ppp_0(join(i, 'OUT.PPP'), ind_ppp_0)
                ppp_per_atom = ppp / num_atom
                atomic_energy = [[element[j], ppp_per_atom*sign, 0.5, 0.5] for j in range(num_atom)]
                
                f_write.write(generate_movement.generate_movement_piece(
                    num_atom, 
                    lattice, 
                    position, 
                    nonperiodic_position=None, 
                    force=force, 
                    velocity=None, 
                    atomic_energy=atomic_energy, 
                    de=None, 
                ))
        
        if (path_existed_movement_dir != None) and (not os.path.exists(join(path_existed_movement_dir, movement_name))):
            utils.cp_r(join(path_saving), join(path_existed_movement_dir, movement_name))


def generate_weight_for_ppp_model_fitting(
    path_job_array, 
    path_saving, 
    ind_atom_0, 
    weight_exist=False, 
    path_existed_weight_dir=None, 
):
    weight_name = 'linear_fitting_with_weight_%d' % (ind_atom_0)
    
    if (weight_exist == True) and (path_existed_weight_dir != None) and (os.path.isfile(join(path_existed_weight_dir, weight_name))):
        os.symlink(join(path_existed_weight_dir, weight_name), join(path_saving))
    else:
        with open(join(path_job_array)) as f:
            job_array = f.read().splitlines()
        
        num_image = len(job_array)
        
        with open(join(path_saving), mode='w') as f_write:
            f_write.write('%d\n' % num_image)
            for i in job_array:
                with open(join(i, 'OUT.PROJECTION_PPP')) as f:
                    projection_ppp = f.read().splitlines()
                
                num_atom = int(projection_ppp[0])
                
                if num_atom < (ind_atom_0 + 1):
                    print('ppp_model_fitting.generate_weight_for_ppp_model_fitting: ERROR!!! num_atom < (ind_atom_0 + 1), num_atom = %d, ind_atom_0 = %d, EXIT!' % (num_atom, ind_atom_0))
                    exit()
                
                f_write.write('%s\n' % (projection_ppp[ind_atom_0+1]))
        
        if (path_existed_weight_dir != None) and (not os.path.exists(join(path_existed_weight_dir, weight_name))):
            utils.cp_r(join(path_saving), join(path_existed_weight_dir, weight_name))


def _get_input_name(input_type, ind_ppp_0, indb, path_input_dir):
    name_input = None
    name_1 = '%s_%s_%d_%d' % (input_type, ind_ppp_0[0], ind_ppp_0[1]+1, ind_ppp_0[2]+1)
    if os.path.isfile(join(path_input_dir, name_1)):
        name_input = name_1
    else:
        ind_atom = None
        if ind_ppp_0[0] == 'atom':
            ind_atom = ind_ppp_0[1] + 1
        elif ind_ppp_0[0] == 'bond':
            ind_atom = indb[ind_ppp_0[1]][ind_ppp_0[2]]
        else:
            raise ValueError("ind_ppp_0[0] should be either 'atom' or 'bond'. However, ind_ppp_0[0] == %s" % (ind_ppp_0[0]))
        
        name_2 = '%s_%d' % (input_type, ind_atom)
        name_3 = input_type
        if os.path.isfile(join(path_input_dir, name_2)):
            name_input = name_2
        elif os.path.isfile(join(path_input_dir, name_3)):
            name_input = name_3
        else:
            raise FileNotFoundError('No %s correspond to %s' % (input_type, str(ind_ppp_0)))
    
    return name_input


def generate_feature(
    path_saving, 
    path_movement, 
    path_movement_type, 
    path_in_neighbore, 
    path_pwmat_mlff, 
    path_parameters, 
    path_log, 
    path_python3, 
    iat_and_idir, 
):
    utils.rm_fr(join(path_saving))
    os.mkdir(join(path_saving))
    
    # file preparation
    os.mkdir(join(path_saving, 'movement'))
    os.symlink(join(path_movement), join(path_saving, 'movement', 'MOVEMENT'))
    os.symlink(join(path_movement_type), join(path_saving, 'movement', 'MOVEMENT.type'))
    
    os.mkdir(join(path_saving, 'input'))
    with open(join(path_saving, 'input', 'find_neighbore.in'), mode='w') as f:
        f.write('3\n')
    os.symlink(join(path_in_neighbore), join(path_saving, 'input', 'IN.NEIGHBORE'))
    
    os.mkdir(join(path_saving, 'output'))
    
    # utils.rm_fr(join(path_pwmat_mlff, 'parameters.py'))
    # utils.cp_r(join(path_parameters), join(path_pwmat_mlff, 'parameters.py'))
    
    cwd = os.getcwd()
    os.chdir(join(path_saving))
    os.system('%s %s True False False 4 "%s" > %s 2> %s' % (
        join(path_python3), 
        join(path_pwmat_mlff, 'main_for_auto.py'), 
        iat_and_idir, 
        join(path_log, 'generate_feature_%s.out' % os.path.basename(path_saving)), 
        join(path_log, 'generate_feature_%s.err' % os.path.basename(path_saving))
    )
    )
    os.chdir(join(cwd))


def combine_feature(path_saving, path_feature_dir):
    utils.rm_fr(path_saving)
    os.mkdir(path_saving)
    
    with open(join(path_saving, 'location'), mode='w') as f_location:
        with open(join(path_saving, 'trainData.txt.Ftype1'), mode='w') as f_type1:
            with open(join(path_saving, 'trainData.txt.Ftype2'), mode='w') as f_type2:
                file_feature_dir = os.listdir(path_feature_dir)
                
                num_feature = len(file_feature_dir)
                f_location.write('%d\n' % (num_feature))
                f_location.write('%s\n' % (join(path_saving)))
                
                ind = 0
                for i in file_feature_dir:
                    os.symlink(join(path_feature_dir, i, 'movement'), join(path_saving, str(ind)))
                    f_location.write('%s\n' % (join(path_saving, str(ind))))
                    
                    if not os.path.isfile(join(path_feature_dir, i, 'trainData.txt.Ftype1')):
                        utils.rm_fr(join(path_feature_dir, i, 'trainData.txt.Ftype1'))
                        f = open(join(path_feature_dir, i, 'trainData.txt.Ftype1'), mode='w')
                        f.close()
                    with open(join(path_feature_dir, i, 'trainData.txt.Ftype1')) as f:
                        f_type1.write(f.read())
                    
                    if not os.path.isfile(join(path_feature_dir, i, 'trainData.txt.Ftype2')):
                        utils.rm_fr(join(path_feature_dir, i, 'trainData.txt.Ftype2'))
                        f = open(join(path_feature_dir, i, 'trainData.txt.Ftype2'), mode='w')
                        f.close()
                    with open(join(path_feature_dir, i, 'trainData.txt.Ftype2')) as f:
                        f_type2.write(f.read())
                    
                    ind += 1
                
    utils.cp_r(join(path_feature_dir, os.listdir(path_feature_dir)[0], 'input'), join(path_saving, 'input'))
    for i in ['find_neighbore.in', 'IN.NEIGHBORE', 'location']:
        os.remove(join(path_saving, 'input', i))
    utils.cp_r(join(path_feature_dir, os.listdir(path_feature_dir)[0], 'output'), join(path_saving, 'output'))


def pca(path_fitting, arg_load_intel, path_pwmat_mlff, path_lmod_init, path_python3, path_log, iflag_weight):
    sys.path.insert(0, path_lmod_init)
    import env_modules_python
    env_modules_python.module(*arg_load_intel)
    
    # pca
    cwd = os.getcwd()
    os.chdir(join(path_fitting))
    print("@@@ here 1 ", path_python3, path_pwmat_mlff)
    os.system('%s %s False True %s 0 > %s 2> %s' % (
        join(path_python3), 
        join(path_pwmat_mlff, 'main_for_auto.py'), 
        iflag_weight, 
        join(path_log, 'pca.out'), 
        join(path_log, 'pca.err')
    )
    )

    print("@@@ here 2 ")

    os.chdir(join(cwd))
    
    sys.path = sys.path[1:]
    

def fit_and_evaluate(path_fitting, path_fit_linearmm_input, path_lmod_init, arg_load_intel, path_log, iflag_weight):
    num_feature = None
    with open(join(path_fitting, 'location')) as f:
        num_feature = int(f.readline().split()[0])
    
    with open(join(path_fitting, 'fread_dfeat', 'fit_linearMM.input'), mode='w') as f_write:
        with open(join(path_fit_linearmm_input)) as f_read:
            fit_linearMM_input = f_read.read().splitlines()
        
        if iflag_weight:
            fit_linearMM_input[0] += ' 1'
        else:
            fit_linearMM_input[0] += ' 0'
        
        print("@@@ fit_linearMM_input[0] ", fit_linearMM_input[0])

        fit_linearMM_input.append('%d' % (num_feature))
        
        for i in range(num_feature):
            fit_linearMM_input.append('%.2f' % (1.0))
        
        for i in fit_linearMM_input:
            f_write.write(i+'\n')
    
    sys.path.insert(0, path_lmod_init)
    import env_modules_python
    env_modules_python.module(*arg_load_intel)
    
    cwd = os.getcwd()
    os.chdir(join(path_fitting, 'fread_dfeat'))
    os.system('make linear_fitB.ntype > %s 2> %s' % (join(path_log, 'fitting.out'), join(path_log, 'fitting.err')))
    os.system('make energyL.pred.1 > %s 2> %s' % (join(path_log, 'pred.out'), join(path_log, 'pred.err')))
    with open(join(path_fitting, 'RMSE'), mode='w') as f_write:
        with open(join(path_log, 'pred.out')) as f_read:
            rmse = f_read.readlines()[-4:]
            f_write.writelines(rmse)
    os.chdir(join(cwd))
    
    sys.path = sys.path[1:]
    
    return rmse


def wrap(
    path_fitting, 
    path_ppp_model, 
    path_ppp_neighbore, 
    path_in_neighbore_dir, 
    path_movement_type_dir, 
    ind_model_0, 
    list_ind_ppp_0, 
    indb, 
):
    if not os.path.isdir(join(path_ppp_model)):
        os.mkdir(join(path_ppp_model))
    if not os.path.isdir(join(path_ppp_neighbore)):
        os.mkdir(join(path_ppp_neighbore))    
    
    path_model_saving = join(path_ppp_model, '%s_%d' % (ind_model_0[0], ind_model_0[1]+1))
    utils.rm_fr(join(path_model_saving))
    os.mkdir(join(path_model_saving))
    for i in ['feat.info', 'fit_linearMM.input', 'linear_fitB.ntype']:
        utils.cp_r(join(path_fitting, 'fread_dfeat', i), join(path_model_saving, i))
    for i in ['feat_PV.']:
        for j in os.listdir(join(path_fitting, 'fread_dfeat')):
            if j.startswith(i):
                utils.cp_r(join(path_fitting, 'fread_dfeat', j), join(path_model_saving, j))
    for i in ['input', 'output', 'RMSE']:
        utils.cp_r(join(path_fitting, i), join(path_model_saving, i))
    
    for i in list_ind_ppp_0:
        name_in_neighbore = _get_input_name('IN.NEIGHBORE', i, indb, path_in_neighbore_dir)
        name_movement_type = _get_input_name('MOVEMENT.type', i, indb, path_movement_type_dir)
        
        path_neighbore_saving = join(path_ppp_neighbore, '%s_%d_%d' % (i[0], i[1]+1, i[2]+1))
        utils.rm_fr(join(path_neighbore_saving))
        os.mkdir(join(path_neighbore_saving))
    
        os.mkdir(join(path_neighbore_saving, 'input'))
        utils.cp_r(join(path_in_neighbore_dir, name_in_neighbore), join(path_neighbore_saving, 'input', 'IN.NEIGHBORE'))
        with open(join(path_neighbore_saving, 'input', 'find_neighbore.in'), mode='w') as f:
            f.write('3\n')
        
        utils.cp_r(join(path_movement_type_dir, name_movement_type), join(path_neighbore_saving, 'input', 'MOVEMENT.type'))


def ppp_model_fitting(
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
    remove_temp_after_fitting=True, 
    movement_exist=False, 
    path_existed_movement_dir=None, 
    weight_exist=False, 
    path_existed_weight_dir=None, 
):
    utils.rm_fr(join(path_saving))
    os.mkdir(join(path_saving))
    
    utils.rm_fr(join(path_log))
    os.mkdir(join(path_log))
    
    # step 1: read polar_param.input, get the list of ind_ppp_0 according to ind_model_0 (get all 
    # qualified ppp according to the specified model)
    print('-------------------------------------------------------------')
    print('ppp_model_fitting.ppp_model_fitting: STEP_1_PPP_DETERMINATION')
    print('-------------------------------------------------------------')
    list_ind_ppp_0, list_sign, ipol, indb = generate_ind_ppp_0_according_to_polar_param_input(path_polar_param_input, ind_model_0)
    for i in range(len(list_ind_ppp_0)):
        print('%s %d %d %d, sign: %d' % tuple(list(list_ind_ppp_0[i])+[list_sign[i]]))
    
    # step 2: generate MOVEMENT files according to the list of ind_ppp_0 obtained in step 1
    # (also generate linear_fitting_with_weight if ind_ppp_0[0]=='atom' and ind_ppp_0[2]==2)
    print('---------------------------------------------------------------')
    print('ppp_model_fitting.ppp_model_fitting: STEP_2_MOVEMENT_GENERATION')
    print('---------------------------------------------------------------')
    os.mkdir(join(path_saving, 'movement'))
    for i in range(len(list_ind_ppp_0)):
        print('MOVEMENT_%s_%d_%d_%d' % tuple(list_ind_ppp_0[i]))
        generate_movement_for_ppp_model_fitting(
            path_job_array=path_job_array, 
            path_saving=join(path_saving, 'movement', 'MOVEMENT_%s_%d_%d_%d' % tuple(list_ind_ppp_0[i])), 
            ind_ppp_0=list_ind_ppp_0[i], 
            sign=list_sign[i], 
            movement_exist=movement_exist, 
            path_existed_movement_dir=path_existed_movement_dir, 
        )
    # if (list_ind_ppp_0[0][0]=='atom') and (list_ind_ppp_0[0][2]==2):
        # os.mkdir(join(path_saving, 'linear_fitting_with_weight'))
        # for i in range(len(list_ind_ppp_0)):
            # print('linear_fitting_with_weight_%d' % (list_ind_ppp_0[i][1]))
            # generate_weight_for_ppp_model_fitting(
                # path_job_array=path_job_array, 
                # path_saving=join(path_saving, 'linear_fitting_with_weight', 'linear_fitting_with_weight_%d' % (list_ind_ppp_0[i][1])), 
                # ind_atom_0=list_ind_ppp_0[i][1], 
                # weight_exist=weight_exist, 
                # path_existed_weight_dir=path_existed_weight_dir, 
            # )
    
    # step 3: generate feature using MOVEMENT files obtained in step 2 and corresponding 
    # IN.NEIGHBORE files
    print('--------------------------------------------------------------')
    print('ppp_model_fitting.ppp_model_fitting: STEP_3_FEATURE_GENERATION')
    print('--------------------------------------------------------------')
    os.mkdir(join(path_saving, 'feature'))
    for i in list_ind_ppp_0:        
        name_in_neighbore = _get_input_name('IN.NEIGHBORE', i, indb, path_in_neighbore_dir)
        name_movement_type = _get_input_name('MOVEMENT.type', i, indb, path_movement_type_dir)
                
        iat_and_idir = None
        if (i[0]=='atom') and (i[2]==0):
            iat_and_idir = '%d %d %d %d 1' % tuple(ipol[i[1]])
        elif (i[0]=='atom') and (i[2]==1):
            iat_and_idir = '%d %d %d %d 2' % tuple(ipol[i[1]])
        elif (i[0]=='atom') and (i[2]==2):
            iat_and_idir = '%d %d %d %d 3' % tuple(ipol[i[1]])
        elif (i[0]=='bond'):
            iat_and_idir = '%d %d 0 0 1' % tuple(indb[i[1]][0:0+2])
        
        print('%s %d %d %d, IN.NEIGHBORE: %s, MOVEMENT.type: %s, iat_and_idir: %s' % tuple(list(i)+[name_in_neighbore, name_movement_type, iat_and_idir]))
        
        generate_feature(
            path_saving=join(path_saving, 'feature', '%s_%d_%d_%d' % tuple(i)), 
            path_movement=join(path_saving, 'movement', 'MOVEMENT_%s_%d_%d_%d' % tuple(i)), 
            path_movement_type=join(path_movement_type_dir, name_movement_type), 
            path_in_neighbore=join(path_in_neighbore_dir, name_in_neighbore), 
            path_pwmat_mlff=path_pwmat_mlff, 
            path_parameters=path_parameters, 
            path_log=path_log, 
            path_python3=path_python3, 
            iat_and_idir=iat_and_idir, 
        )
    
    # step 4: combine generated features, do pca, and fit
    print('-------------------------------------------------------------------')
    print('ppp_model_fitting.ppp_model_fitting: STEP_4_COMBINATION_PCA_FITTING')
    print('-------------------------------------------------------------------')
    print('STEP_4_COMBINATION')
    combine_feature(join(path_saving, 'fitting'), join(path_saving, 'feature'))
    # if (list_ind_ppp_0[0][0]=='atom') and (list_ind_ppp_0[0][2]==2):
        # print('STEP_4_LINK_WEIGHT')
        # for i in range(len(list_ind_ppp_0)):
            # os.symlink(
                # join(path_saving, 'linear_fitting_with_weight', 'linear_fitting_with_weight_%d' % (list_ind_ppp_0[i][1])), 
                # join(path_saving, 'fitting', str(i), 'linear_fitting_with_weight')
            # )
    print('STEP_4_PCA')
    # if (list_ind_ppp_0[0][0]=='atom') and (list_ind_ppp_0[0][2]==2):
        # iflag_weight = True
    # else:
        # iflag_weight = False
    iflag_weight = False
    pca(join(path_saving, 'fitting'), arg_load_intel, path_pwmat_mlff, path_lmod_init, path_python3, path_log, iflag_weight)
    print('STEP_4_FITTING')
    rmse = fit_and_evaluate(join(path_saving, 'fitting'), path_fit_linearmm_input, path_lmod_init, arg_load_intel, path_log, iflag_weight)
    for i in rmse:
        print(i, end='')
    
    # step 5: wrap
    print('------------------------------------------------')
    print('ppp_model_fitting.ppp_model_fitting: STEP_5_WRAP')
    print('------------------------------------------------')
    wrap(
        path_fitting=join(path_saving, 'fitting'), 
        path_ppp_model=path_ppp_model, 
        path_ppp_neighbore=path_ppp_neighbore, 
        path_in_neighbore_dir=path_in_neighbore_dir, 
        path_movement_type_dir=path_movement_type_dir, 
        ind_model_0=ind_model_0, 
        list_ind_ppp_0=list_ind_ppp_0, 
        indb=indb, 
    )
    
    # step 6: remove temporary directory
    print('------------------------------------------------------------------')
    print('ppp_model_fitting.ppp_model_fitting: STEP_6_TEMP_DIRECTORY_REMOVAL')
    print('------------------------------------------------------------------')
    if remove_temp_after_fitting:
        utils.rm_fr(join(path_saving))
        print('temp directory removed')
    else:
        print('temp directory retained')
