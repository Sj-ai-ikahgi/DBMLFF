
import sys
import numpy as np
from copy import deepcopy
from units_dbmlff import *



MolEindexList = {
  "EC":  [8,8,8,6,6,\
          6,1,1,1,1],
  "Li":  [3],
  "PF6": [15,9,9,9,9,\
          9,9],
  "DEC": [8,8,8,6,6,\
          6,6,6,1,1,\
          1,1,1,1,1,\
          1,1,1],
  "DMC": [8,8,8,6,6,\
          6,1,1,1,1,\
          1,1],
  "PC":  [8,8,8,6,6,\
          6,6,1,1,1,\
          1,1,1],
  "EMC": [8,8,8,6,6,\
          6,6,1,1,1,\
          1,1,1,1,1],
  "VC":  [8,8,8,6,6,\
          6,1,1],
  "TTE": [9,9,9,9,9,\
          9,9,9,8,6,\
          6,6,6,6,1,\
          1,1,1],
  "DME": [8,6,6,1,1,\
          1,1,1,1],
  "EA": [6,6,6,6,8,\
          8,1,1,1,1,\
          1,1,1,1]
}


MolAtomSysList = {
  "EC": [1,1,2,3,3,\
         4,5,5,5,5],
  "Li": [1],
  "PF6": [1,2,2,2,2,\
          2,2],
  "DEC": [1,2,2,3,4,\
          4,5,5,6,6,\
          7,7,7,6,6,\
          7,7,7],
  "DMC": [1,2,2,3,4,\
          4,5,5,5,5,\
          5,5],
  "PC":  [1,2,3,4,5,\
          6,7,8,8,9,\
          9,9,10],
  "EMC": [1,2,3,4,5,\
          6,7,8,8,9,\
          9,9,10,10,10],
  "VC":  [1,2,2,3,4,\
          4,5,5],
  "TTE": [1,1,2,2,3,\
          3,4,4,5,6,\
          7,8,9,10,11,\
          11,12,13],
  "DME": [1,2,2,3,3,\
          3,3,3,3],
   "EA": [1,2,3,4,5,\
          6,7,7,7,8,\
          8,9,9,9]
}


zValenceList   = {1:1.0, 6:4.0,  8:6.0,  9:7.0,  15:5.0}
atomicMassList = {1:1.0, 6:12.0, 8:16.0, 9:19.0, 15:31.0}
moveStartAtomId = {"EC":3, "Li":1, "PF6":7, "DEC":18, "DMC":12, \
                   "PC":3, "EMC":15, "VC":1, "TTE":5, "DME":6, "EA":7}

# -------------------------------------------------
#  change real coordinates to fraction coordinates
# -------------------------------------------------
def Cart2Lamda(lattice, coords):

  matrix = np.array(lattice).reshape((3,3))
  inv_matrix = np.linalg.inv(matrix)

  frac_coord = []
  for r in range(len(coords)):
    frac_coord.append(np.dot(coords[r],inv_matrix))

  return frac_coord

# -------------------------------------------------
#  change fraction coordinates to real coordinates
# -------------------------------------------------
def Lamda2Cart(lattice, xs):

  xhi = lattice[0][0]
  yhi = lattice[1][1]
  zhi = lattice[2][2]
  xy  = lattice[1][0]
  xz  = lattice[2][0]
  yz  = lattice[2][1]

  Cart = []
  for i in range(len(xs)):
    tx = xs[i][0]*xhi + xs[i][1]*xy  + xs[i][2]*xz
    ty =                xs[i][1]*yhi + xs[i][2]*yz
    tz =                               xs[i][2]*zhi
    Cart.append([tx,ty,tz])

  return Cart



def Composition(elem):

  natoms = len(elem)
  atype = np.zeros(natoms, dtype=int)
  charge = np.zeros(natoms, dtype=int)
  eindex = np.zeros(natoms, dtype=int)
  elem_list = list(dict.fromkeys(elem))
  num_atype = len(elem_list)  
  mass = np.zeros(num_atype, dtype=int)

  for i in range(natoms):
    eindex[i] = element_index[elem[i]]
    atype[i] = elem_list.index(elem[i]) + 1
    charge[i] = charge_list[elem[i]]

  for i in range(num_atype):
    mass[i] = mass_list[elem_list[i]]

  return num_atype, atype, mass, charge, eindex, elem_list


def eindex2atype(eindex):

  natoms = len(eindex)
  atype = np.zeros(natoms, dtype=int)

  eindex_list = list(dict.fromkeys(eindex))
  eindex_list.sort()
  elem_list = []
  emass_list = []
  natype = len(eindex_list)

  for i in range(natype): 
    elem = index_element[eindex_list[i]]
    elem_list.append(elem)
    emass_list.append(mass_list[elem])

  for i in range(natoms):
    atype[i] = eindex_list.index(eindex[i]) + 1

  return elem_list, eindex_list, emass_list, natype, atype



def New_MOVEMENT_Conf(isfile):

  headline = ''

  line = isfile.readline()
  while line:
    if "Iteration" in line:
      headline = line
      break
    line = isfile.readline()

  return headline


def Get_One_MOVEMENT_Conf(headline, isfile):

  tokens = headline.split()

  natoms = int(tokens[0])

  #ntimestep = float(tokens[4][:-1])
  ntimestep = -1 

  Ep = 0.0
  lattice = []
  eindex = []
  coords = []
  forces = []
  velocity = []
  atomic_eng = []

  for i,t in enumerate(tokens):
    if t == "Etot,Ep,Ek":
      Ep = float(tokens[i+4])

  while True:

    line = isfile.readline()

    if "Lattice vector (Angstrom)" in line:
      for j in range(3):
        tokens = isfile.readline().split()[:3]
        lattice.append([float(c) for c in tokens])

    elif ("Position (normalized), move_x, move_y, move_z" in line) and \
      (not ("nonperiodic_Position" in line)): 
      for j in range(natoms):
        tokens = isfile.readline().split()
        eindex.append(int(tokens[0]))
        coords.append([float(t) for t in tokens[1:4]])

    elif "Force (-force, eV/Angstrom)" in line:
      for j in range(natoms):
        tokens = isfile.readline().split()
        forces.append([float(f) for f in tokens[1:4]])

    elif "Velocity (bohr/fs)" in line:
      for j in range(natoms):
         tokens = isfile.readline().split()
         velocity.append([float(v) for v in tokens[1:4]])

    elif "Atomic-Energy, Etot(eV),E_nonloc(eV),Q_atom:dE(eV)= " in line:
      for j in range(natoms):
         tokens = isfile.readline().split()
         atomic_eng.append([float(t) for t in tokens[1:4]])
    elif ("----------------" in line) or\
         (line == ''):
      break

  if (natoms == 0) or \
     (abs(Ep) < 1.0e-9) or\
     (np.shape(lattice)[1] != 3) or \
     (len(eindex) != natoms) or\
     (np.shape(coords)[0] != natoms): 

    print("@@@ READ MOVEMENT CONF ERROR")
    print("natoms, Etot, np.shape(lattice) ", end=' ')
    print("len(eindex), np.shape(coords), \
           np.shape(forces), np.shape(velocity)")
    print(natoms, abs(Ep), np.shape(lattice), \
          len(eindex), np.shape(coords), \
          np.shape(forces), np.shape(velocity))
    sys.exit()

  lattice = np.array(lattice)
  eindex = np.array(eindex)
  coords = np.array(coords)

  if len(forces) != natoms:
    forces = np.zeros([natoms,3], dtype=float)

  if len(velocity) != natoms:
    velocity = np.zeros([natoms,3], dtype=float)

  if len(atomic_eng) != natoms:
    atomic_eng = np.zeros([natoms,3], dtype=float)

  forces = np.array(forces)
  velocity = np.array(velocity)
  atomic_eng = np.array(atomic_eng)

  return ntimestep, Ep, natoms, lattice, eindex, \
         coords, forces, velocity, atomic_eng 



def Read_One_Config(isfile):
  
  line = isfile.readline()
  natoms = int(line.split()[0])
  lattice = []
  eindex = []
  coords = []

  while True:
    line = isfile.readline()
    if "Lattice vector (Angstrom)" in line:
      for j in range(3):
        tokens = isfile.readline().split()
        lattice.append([float(c) for c in tokens[:3]])

    elif "Position (normalized), move_x, move_y, move_z" in line:
      for j in range(natoms):
        tokens = isfile.readline().split()
        eindex.append(int(tokens[0]))
        coords.append([float(j) for j in tokens[1:4]])
    elif line == '':
      break

  lattice = np.array(lattice)
  eindex = np.array(eindex)
  coords = np.array(coords)
 
  return natoms, lattice, eindex, coords



def Write_LAMMPS_Dump(\
  cid, \
  ofilename, \
  natoms, \
  lattice,
  atype, \
  coords):

  xlo = ylo = zlo = 0.0
  xhi = lattice[0][0]
  yhi = lattice[1][1]
  zhi = lattice[2][2]
  xy = lattice[1][0]
  xz = lattice[2][0]
  yz = lattice[2][1]

  tril = False
  if not( (abs(xy) < 1.0e-6) and \
     (abs(xz) < 1.0e-6) and \
     (abs(yz) < 1.0e-6)):
    tril = True
 
  xlo_bound = xlo + min([0.0,xy,xz,xy+xz])
  xhi_bound = xhi + max([0.0,xy,xz,xy+xz])
  ylo_bound = ylo + min([0.0,yz])
  yhi_bound = yhi + max([0.0,yz])
  zlo_bound = zlo
  zhi_bound = zhi

  osfile = open(ofilename, 'a')

  osfile.write("ITEM: TIMESTEP\n")
  osfile.write("%-9d\n" % (cid))
  osfile.write("ITEM: NUMBER OF ATOMS\n")
  osfile.write("%-12d\n" % (natoms))

  if tril:
    osfile.write("ITEM: BOX BOUNDS xy xz yz pp pp pp\n")
    osfile.write("%16.12f %16.12f %16.12f xlo xhi\n" %  (xlo_bound, xhi_bound, xy))
    osfile.write("%16.12f %16.12f %16.12f ylo yhi\n" %  (ylo_bound, yhi_bound, xz))
    osfile.write("%16.12f %16.12f %16.12f zlo zhi\n" %  (zlo_bound, zhi_bound, yz))
  else:
    osfile.write("ITEM: BOX BOUNDS pp pp pp\n")
    osfile.write("%16.12f %16.12f xlo xhi\n" %  (xlo_bound, xhi_bound))
    osfile.write("%16.12f %16.12f ylo yhi\n" %  (ylo_bound, yhi_bound))
    osfile.write("%16.12f %16.12f zlo zhi\n" %  (zlo_bound, zhi_bound))

  osfile.write("ITEM: ATOMS id type x y z\n")

  for i in range(natoms):
    if type(atype[i]) == type(1):
      osfile.write("%9d %5d %16.12f %16.12f %16.12f\n" %
                  (i+1, atype[i], \
                   coords[i][0], coords[i][1], coords[i][2]))
    else:
      osfile.write("%9d %5s %16.12f %16.12f %16.12f\n" % \
                  (i+1, atype[i], \
                  coords[i][0], coords[i][1], coords[i][2]))

  osfile.close()


def Write_Atom_Config(osfile,lattice,eindex,coords,MD=''):
  natoms = len(coords)
  if MD == '':
    osfile.write(" %7d \n" % (natoms))
  elif MD == "DBMD":
    osfile.write(" %7d 4 10 1 7 18 0 2 2 0 \n" % (natoms))

  osfile.write("Lattice vector (Angstrom)\n")
  for r in range(3):
    for c in range(3):
      osfile.write("%12.9f " % (lattice[r][c]))
    osfile.write("\n")

  osfile.write("Position (normalized), move_x, move_y, move_z\n")
  for i in range(natoms):
    osfile.write("%7d %12.9f %12.9f %12.9f 1 1 1\n" \
      % (eindex[i], coords[i][0], coords[i][1], coords[i][2]))


def Write_One_MOVEMENT(
  osfile, \
  oheadline,\
  lattice,\
  eindex,\
  coords,\
  forces = '',\
  velocity = '',\
  atomic_eng = 0.0):

  natoms = len(eindex)

  osfile.write(oheadline)

  osfile.write("Lattice vector (Angstrom)\n")
  for r in range(3):
    for c in range(3):
      osfile.write("%12.9f " % (lattice[r][c]))
    osfile.write("\n")

  osfile.write("Position (normalized), move_x, move_y, move_z\n")
  for aid in range(natoms):
    osfile.write("%   -7d %18.16f %18.16f %18.16f 1 1 1\n" \
      % (eindex[aid], coords[aid][0], coords[aid][1], coords[aid][2]))

  if len(forces) > 0:
    osfile.write("Force (-force, eV/Angstrom)\n")
    for aid in range(natoms):
      osfile.write("%   -7d %20.16f %20.16f %20.16f\n" \
        % (eindex[aid], forces[aid][0], forces[aid][1], forces[aid][2]))

  if len(velocity) > 0:
    osfile.write("Velocity (bohr/fs)\n")
    for aid in range(natoms):
      osfile.write("%   -7d %12.9f %12.9f %12.9f\n" \
        % (eindex[aid], velocity[aid][0], velocity[aid][1], velocity[aid][2]))

  if abs(atomic_eng) > 0:
    osfile.write("Atomic-Energy, Etot(eV),E_nonloc(eV),Q_atom:dE(eV)=  %15.6f\n" % atomic_eng)
    aeng = atomic_eng/natoms
    for aid in range(natoms):
      osfile.write("%   -7d %12.9f %12.9f %12.9f\n" \
        % (eindex[aid], aeng, 0.0, 0.0))

  osfile.write("-----------------------------------------\n")

  return



def Write_MLFF_Init_Config(
  osfile, \
  oheadline,\
  lattice,\
  atomId,\
  eindex,\
  coords,\
  forces,\
  velocity):

  natoms = len(eindex)

  osfile.write(oheadline)

  osfile.write("Lattice vector (Angstrom)\n")
  for r in range(3):
    for c in range(3):
      osfile.write("%12.9f " % (lattice[r][c]))
    osfile.write("\n")

  osfile.write("Position (normalized), move_x, move_y, move_z\n")
  for i in range(natoms):
    aid = atomId[i]
    osfile.write("%-7d %18.16f %18.16f %18.16f 1 1 1\n" \
      % (eindex[aid], coords[aid][0], coords[aid][1], coords[aid][2]))

  osfile.write("Force (-force, eV/Angstrom)\n")
  for i in range(natoms):
    aid = atomId[i]
    osfile.write("%-7d %18.16f %18.16f %18.16f\n" \
      % (eindex[aid], forces[aid][0], forces[aid][1], forces[aid][2]))

  osfile.write("Velocity (bohr/fs)\n")
  for i in range(natoms):
    aid = atomId[i]
    osfile.write("%-7d %12.9f %12.9f %12.9f\n" \
      % (eindex[aid], velocity[aid][0], velocity[aid][1], velocity[aid][2]))

  return

def Write_Pwmat_Init(\
  osfile,\
  lattice,\
  atomId,\
  eindex,\
  coords,\
  Ep=0.0):

  natoms = len(coords)
  osfile.write("%9d %20.9f\n" % (natoms, Ep))

  osfile.write("Lattice vector\n")
  for r in range(3):
    for c in range(3):
      osfile.write("%12.9f " % (lattice[r][c]))
    osfile.write("\n")

  osfile.write("Position, move_x, move_y, move_z\n")
  for i in range(natoms):
    aid = atomId[i]
    osfile.write("%-7d %18.16f %18.16f %18.16f 1 1 1\n" \
      % (eindex[aid], coords[aid][0], coords[aid][1], coords[aid][2]))

 
def Write_XYZ(\
  cid, \
  ofilename, \
  natoms, \
  eindex, \
  coords, \
  lattice = ''):

  osfile = open(ofilename, 'a')

  if lattice == '':
    osfile.write("%12d\n" % natoms)
  else:
    osfile.write("%12d %20.9f %20.9f %20.9f\n" % \
      (natoms, lattice[0][0], lattice[1][1], lattice[2][2]))
  
  osfile.write("Iteration %12d \n" % (cid))

  for i in range(natoms):
    osfile.write("%5s %16.9f %16.9f %16.9f\n" % \
                 (index_element[eindex[i]], \
                  coords[i][0], coords[i][1], coords[i][2]))

  osfile.close()



def two_point_distance(coord1, coord2, AL):

  dr = np.zeros(3, dtype=float)

  for j in range(3):
    dx = coord2[j] - coord1[j]
    if abs(dx+1.0) < abs(dx): dx += 1.0
    if abs(dx-1.0) < abs(dx): dx -= 1.0
    dr[j] = dx*AL[j][j]

  d = np.sqrt(np.dot(dr,dr))

  return dr, d



def two_atom_distance(xatom, AL, iat, jat):

  dr = np.zeros(3, dtype=float)

  for j in range(3):
    dx = xatom[jat][j] - xatom[iat][j]
    if abs(dx+1.0) < dx: dx += 1.0
    if abs(dx-1.0) < dx: dx -= 1.0
    dr[j] = dx*AL[j][j]

  d = np.sqrt(np.dot(dr,dr))

  return dr, d


def Move_Atom_Together(natoms, coords, caId):

  new_coords = deepcopy(coords)
  center_coord = np.zeros(3, dtype=float)

  for aid in range(natoms):
    for j in range(3):
      dr = new_coords[aid][j] - new_coords[caId][j]
      if np.fabs(dr) > 0.5:
        new_coords[aid][j] -= np.sign(dr)
      center_coord[j] += new_coords[aid][j]

  return center_coord/natoms, new_coords


def Move_Enlarge_Box(new_box_len, old_lattice, natoms, coords):

  old_coords = Lamda2Cart(old_lattice, coords)

  center_coord = np.zeros(3, dtype=float)
  for i in range(natoms):
    for j in range(3):
      center_coord[j] += old_coords[i][j]
  center_coord /= natoms

  for i in range(natoms):
    for j in range(3):
      old_coords[i][j] = old_coords[i][j] - center_coord[j] + 0.5*new_box_len

  new_lattice = np.zeros([3,3], dtype=float)
  for i in range(3):
    new_lattice[i][i] = new_box_len

  new_coords = Cart2Lamda(new_lattice, old_coords)

  return new_lattice, new_coords