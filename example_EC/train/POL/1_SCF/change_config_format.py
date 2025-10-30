import numpy as np 


def Read_Config(ifilename):

  isfile = open(ifilename, 'r')

  token = isfile.readline().split()
  natoms = int(token[0])

  isfile.readline()
  AL = []
  for i in range(3):
    token = isfile.readline()
    AL.append(token)

  isfile.readline()
  coords = []
  for i in range(natoms):
    token = isfile.readline()
    coords.append(token)

  isfile.close()

  return natoms, AL, coords


def Write_Config(ofilename, natoms, AL, coords):

  osfile = open(ofilename, 'w')

  osfile.write("  %7d\n" % natoms)
  osfile.write(" Lattice vector (Angstrom)\n")
  for i in range(3):
    osfile.write(AL[i])

  osfile.write(" Position (normalized), move_x, move_y, move_z\n")
  for i in range(natoms):
    osfile.write(coords[i])

  osfile.close()


natoms, AL, coords = Read_Config("./xatom0.config")
Write_Config("./xatom0.config", natoms, AL, coords)

