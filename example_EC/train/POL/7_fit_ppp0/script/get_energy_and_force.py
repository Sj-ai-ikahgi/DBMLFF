import os
from os.path import join

import sys

def get_dft_energy(path, energy_type):
    if energy_type == 'Etot':
        if os.path.isfile(join(path, 'REPORT')):
            with open(join(path, 'REPORT')) as f:
                data = f.read().splitlines()

            for i in data:
                if i.startswith(' E_tot(eV)    ='):
                    return float(i.split()[2])
        else:
            return None
    elif energy_type == 'Etot,ENDIV':
        if os.path.isfile(join(path, 'OUT.ENDIV')):
            with open(join(path, 'OUT.ENDIV')) as f:
                return float(f.readline().split()[2])
        else:
            return None


def get_dbmlff_energy(path, energy_type):
    if os.path.isfile(join(path, 'OUT.ENERGY')):
        with open(join(path, 'OUT.ENERGY')) as f:
            data = f.read().splitlines()

        if energy_type == 'Etot':
            for i in data:
                if i.startswith(' Etot(eV)='):
                    return float(i.split()[1])
        elif energy_type == 'Eintra,MLFF':
            for i in data:
                if i.startswith(' Eintra,MLFF(eV)='):
                    return float(i.split()[1])                
        elif energy_type == 'Einter,DB':
            for i in data:
                if i.startswith(' Einter,DB(eV)='):
                    return float(i.split()[1])                
        elif energy_type == 'Einter,DB_tot':
            for i in data:
                if i.startswith(' Einter,DB_tot(eV)='):
                    return float(i.split()[1])                
        elif energy_type == 'Einter,DB_M':
            for i in data:
                if i.startswith(' Einter,DB_M(eV)='):
                    return float(i.split()[1])                
        elif energy_type == 'Einter,POL':
            for i in data:
                if i.startswith(' Einter,POL(eV)='):
                    return float(i.split()[1])                
        elif energy_type == 'Einter,PC':
            for i in data:
                if i.startswith(' Einter,PC(eV)='):
                    return float(i.split()[1])                
        elif energy_type == 'Einter,DB+POL':
            for i in data:
                if i.startswith(' Einter,DB+POL(eV)='):
                    return float(i.split()[1])                
        elif energy_type == 'Einter,DB+POL+PC':
            for i in data:
                if i.startswith(' Einter,DB+POL+PC(eV)='):
                    return float(i.split()[1])                
    else:
        return None


def get_dbmlff_force(path_dir, force_type):
    num_atom = None
    with open(join(path_dir, 'atom.config')) as f:
        content = f.readline()
        while 'atoms,Iteration' not in content:
            content = f.readline()
        num_atom = int(content.split()[0])
    
    if os.path.isfile(join(path_dir, 'OUT.FORCE')):
        with open(join(path_dir, 'OUT.FORCE')) as f:
            data = f.read().splitlines()
        
        force = []
        if force_type == 'fatom_tot':
            ind = data.index(' fatom_tot')
            for i in range(ind+1, ind+1+num_atom):
                force.append([float(j) for j in data[i].split()])
        elif force_type == 'fatom_MLFF':
            ind = data.index(' fatom_MLFF')
            for i in range(ind+1, ind+1+num_atom):
                force.append([float(j) for j in data[i].split()])
        elif force_type == 'fatom_DB':
            ind = data.index(' fatom_DB')
            for i in range(ind+1, ind+1+num_atom):
                force.append([float(j) for j in data[i].split()])
        elif force_type == 'fatom_DB_tot':
            ind = data.index(' fatom_DB_tot')
            for i in range(ind+1, ind+1+num_atom):
                force.append([float(j) for j in data[i].split()])
        elif force_type == 'fatom_DB_M':
            ind = data.index(' fatom_DB_M')
            for i in range(ind+1, ind+1+num_atom):
                force.append([float(j) for j in data[i].split()])
        elif force_type == 'fatom_POL':
            ind = data.index(' fatom_POL')
            for i in range(ind+1, ind+1+num_atom):
                force.append([float(j) for j in data[i].split()])
        elif force_type == 'fatom_PC':
            ind = data.index(' fatom_PC')
            for i in range(ind+1, ind+1+num_atom):
                force.append([float(j) for j in data[i].split()])
        elif force_type == 'fatom_DB+fatom_POL':
            ind = data.index(' fatom_DB+fatom_POL')
            for i in range(ind+1, ind+1+num_atom):
                force.append([float(j) for j in data[i].split()])
        elif force_type == 'fatom_DB+fatom_POL+fatom_PC':
            ind = data.index(' fatom_DB+fatom_POL+fatom_PC')
            for i in range(ind+1, ind+1+num_atom):
                force.append([float(j) for j in data[i].split()])
        
        return force
    else:
        return None


def get_energy_from_dft_md_0(path, energy_type, start, end=None, step=None):
    energy = []
    
    if end == None:
        end = start + 1
    if step == None:
        step = 1
    
    if (energy_type=='Etot') or (energy_type=='Ep') or (energy_type=='Ek'):
        if os.path.isfile(join(path, 'MDSTEPS')):
            with open(join(path, 'MDSTEPS')) as f:
                for i in range(start):
                    f.readline()
                
                for i in range(start, end):
                    temp = f.readline()
                    if ((i-start)%step) == 0:
                        if energy_type == 'Etot':
                            energy.append(float(temp.split()[3]))
                        elif energy_type == 'Ep':
                            energy.append(float(temp.split()[4]))
                        elif energy_type == 'Ek':
                            energy.append(float(temp.split()[5]))
            
            return energy
        elif os.path.isfile(join(path, 'MOVEMENT')):
            num_atom = -1
            with open(join(path, 'MOVEMENT')) as f:
                num_atom = int(f.readline().split()[0])
            
            num_line_one_step = get_num_line_one_step(join(path, 'MOVEMENT'))
            
            with open(join(path, 'MOVEMENT')) as f:
                for i in range(start):
                    for j in range(num_line_one_step):
                        f.readline()
                
                for i in range(start, end):
                    for j in range(num_line_one_step):
                        temp = f.readline()
                        
                        if (((i-start)%step)==0) and (j==0):
                            if energy_type == 'Etot':
                                energy.append(float(temp.split()[8]))
                            elif energy_type == 'Ep':
                                energy.append(float(temp.split()[9]))
                            elif energy_type == 'Ek':
                                energy.append(float(temp.split()[10].split(sep=',')[0]))
            
            return energy
        else:
            return None
    elif energy_type == 'Ep,ENDIV':
        if os.path.isfile(join(path, 'MOVEMENT')):
            num_atom = -1
            with open(join(path, 'MOVEMENT')) as f:
                num_atom = int(f.readline().split()[0])
            
            num_line_one_step = get_num_line_one_step(join(path, 'MOVEMENT'))
            ind_line_atomic_energy = get_ind_line_atomic_energy_0(join(path, 'MOVEMENT'))
            
            with open(join(path, 'MOVEMENT')) as f:
                for i in range(start):
                    for j in range(num_line_one_step):
                        f.readline()
                
                for i in range(start, end):
                    ep_endiv = 0.0
                    
                    for j in range(num_line_one_step):
                        temp = f.readline()
                        
                        if (((i-start)%step)==0) and (j in list(range(ind_line_atomic_energy+1, ind_line_atomic_energy+1+num_atom))):
                            ep_endiv += float(temp.split()[1])
                    
                    if ((i-start)%step) == 0:
                        energy.append(ep_endiv)
            
            return energy
        else:
            return None


def get_num_line_one_step(path):
    with open(join(path)) as f:
        f.readline()
        
        counter = 0
        temp = ''
        while 'atoms,Iteration' not in temp:
            temp = f.readline()
            counter += 1
        
    return counter


def get_ind_line_atomic_energy_0(path):
    with open(join(path)) as f:
        counter = 0
        temp = f.readline()
        while 'Atomic-Energy' not in temp:
            temp = f.readline()
            counter += 1
        
    return counter


def get_ppp_0(path_out_ppp, ind_ppp_0):
    with open(join(path_out_ppp)) as f:
        out_ppp = f.read().splitlines()
    
    num_atom = int(out_ppp[0])
    num_bond = int(out_ppp[1])
    
    ind_line = None
    if ind_ppp_0[0] == 'atom':
        ind_line = ind_ppp_0[1] * 6 + ind_ppp_0[2] * 2 + ind_ppp_0[3] + 3
    elif ind_ppp_0[0] == 'bond':
        ind_line = num_atom * 6 + ind_ppp_0[1] * 4 + ind_ppp_0[2] * 2 + ind_ppp_0[3] + 3
    else:
        raise ValueError("ind_ppp_0[0] should be either 'atom' or 'bond'. However, ind_ppp_0[0] == %s" % (ind_ppp_0[0]))
    
    return float(out_ppp[ind_line])


if __name__ == '__main__':
    if sys.argv[1] == 'get_dft_energy':
        print(get_dft_energy(sys.argv[2], sys.argv[3]))
    elif sys.argv[1] == 'get_dbmlff_energy':
        print(get_dbmlff_energy(sys.argv[2], sys.argv[3]))
    elif sys.argv[1] == 'get_energy_from_dft_md_0':
        if len(sys.argv) == 5:
            print(get_energy_from_dft_md_0(sys.argv[2], sys.argv[3], int(sys.argv[4])))
        elif len(sys.argv) == 6:
            print(get_energy_from_dft_md_0(sys.argv[2], sys.argv[3], int(sys.argv[4]), int(sys.argv[5])))
        elif len(sys.argv) == 7:
            print(get_energy_from_dft_md_0(sys.argv[2], sys.argv[3], int(sys.argv[4]), int(sys.argv[5]), int(sys.argv[6])))
    elif sys.argv[1] == 'get_dbmlff_force':
        print(get_dbmlff_force(sys.argv[2], sys.argv[3]))
