import numpy as np

# atom number of 1 mol
Avogadro = 6.02214076*10**23

# g/mol
mH = 1.00794
mLi = 6.9410
mC = 12.0107
mO = 15.9994
mF = 18.998403
mP = 30.973761

def mass(nH, nC, nO):
  return nH*mH + nC*mC + nO*mO

mPF6 = mP + 6*mF
mEC  = mass(4,3,3)
mDEC = mass(10,5,3)
mDMC = mass(6,3,3)
mEMC = mass(8,4,3)
mPC  = mass(6,4,3)
mDME = mass(6,2,1)
mEA = mass(8,4,2)

mTTE = 8*mF + 1*mO + 5*mC + 4*mH



mol_mass = {"EC":mEC, "DEC":mDEC, "DMC":mDMC, "EMC":mEMC, "PC":mPC, \
            "PF6":mPF6, "Li":mLi, "TTE":mTTE, "DME":mDME,"EA":mEA}

mol_Num = {"EC":0, "DEC":0, "DMC":0, "EMC":0, "PC":0, \
           "PF6":0, "Li":0, "TTE":0, "DME": 0, "EA": 15}
atom_num_of_mol = {"EC":10, "DEC":18, "DMC":12, "EMC":15, "PC":13, \
           "PF6":7, "Li":1, "TTE":18, "DME":9, "EA":14}
           
density = 0.902



total_mass = 0.0
total_atom_num = 0
total_mol_num = 0
for key in mol_mass:
  mass = mol_mass[key]
  num = mol_Num[key]
  total_mass += mass*num
  total_atom_num += atom_num_of_mol[key]*num
  total_mol_num += num


lx = pow(total_mass/Avogadro/10**-24/density, 1.0/3)
lz = ly = lx

print("lx %12.6f volume %12.6f " % (lx, pow(lx*0.1, 3)))

# solution concentration -- molarity
# molarity = number of moles of solute / number of liters of solution
mole =  mol_Num["PF6"] / (Avogadro * (lx*ly*lz*10**-27))

print("ionNum %7d  mole %12.6f" % (mol_Num["PF6"], mole))
print("atom num %7d"            % (total_atom_num))
print("mol num %7d"             % (total_mol_num))