
import sys
import numpy as np
from copy import deepcopy
#from units import *


# -------------------------------------------------
#  change real coordinates to fraction coordinates
# -------------------------------------------------
def Cart2Lamda(lattice, coords):

  matrix = np.array(lattice).reshape((3,3))
  inv_matrix = np.linalg.inv(matrix)

  frac_coord = []
  for r in range(len(coords)):
    frac_coord.append(np.dot(coords[r],inv_matrix))

  # frac_coord = np.array(frac_coord)

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

  Cart = np.array(Cart)

  return Cart


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




def New_DUMP_Conf(isfile):

  line = isfile.readline()
  while(line):
    if "TIMESTEP" in line:
      return True 
    line = isfile.readline()
    print(line)
  return False



def Get_One_DUMP_Conf(isfile):
  
  lattice = np.zeros([3,3], dtype=float)
  natoms = 0
  forces = []

  line = isfile.readline()
  ntimestep = int(line.split()[0])
  line = isfile.readline()
  line = isfile.readline()
  natoms = int(line.split()[0])

  line = isfile.readline()

  for i in range(3):
    line = isfile.readline().split()
    lattice[i][i] = float(line[1]) - float(line[0])

  line = isfile.readline()

  atype = np.zeros(natoms, dtype=int)
  coords = np.zeros([natoms,3], dtype=float)

  for i in range(natoms):
    line = isfile.readline()
    tokens = line.split()
    aid = int(tokens[0])-1 
    atype[aid] = int(tokens[1])
    for j in range(3):
      coords[aid][j] = float(tokens[2+j])

  return ntimestep, natoms, lattice, atype, coords, forces




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



def Write_LAMMPS_Dump(\
  cid, \
  ofilename, \
  natoms, \
  lattice,
  atype, \
  coords,\
  forces = []):

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

  if forces == []:
   osfile.write("ITEM: ATOMS id type x y z\n")
  else:
   osfile.write("ITEM: ATOMS id type x y z fx fy fz\n")

  for i in range(natoms):
    if type(atype[i]) == type(1):
      if forces == []:
        osfile.write("%9d %5d %16.12f %16.12f %16.12f\n" %
                    (i+1, atype[i], \
                    coords[i][0], coords[i][1], coords[i][2]))
      else:
        osfile.write("%9d %5d %16.12f %16.12f %16.12f %16.12f %16.12f %16.12f\n" %
                    (i+1, atype[i], \
                    coords[i][0], coords[i][1], coords[i][2], \
                    forces[i][0], forces[i][1], forces[i][2]))
    else:
      if forces == []:
        osfile.write("%9d %5s %16.12f %16.12f %16.12f\n" % \
                    (i+1, atype[i], \
                    coords[i][0], coords[i][1], coords[i][2]))
      else:
        osfile.write("%9d %5s %16.12f %16.12f %16.12f %16.12f %16.12f %16.12f\n" % \
                    (i+1, atype[i], \
                    coords[i][0], coords[i][1], coords[i][2], \
                    forces[i][0], forces[i][1], forces[i][2]))

  osfile.close()


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






