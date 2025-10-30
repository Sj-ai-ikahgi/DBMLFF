

def Get_REPOPRT_E_tot(isfilename):

  E_tot = -1.0

  isfile = open(isfilename, 'r')

  while True:
    token = isfile.readline()
    if token == "":
      break
    elif "E_tot(eV)    =" in token:
      E_tot = float(token.split()[2])

  return E_tot



E_tot = Get_REPOPRT_E_tot("./REPORT")

isfile = open("./OUT.ENDIV", 'r')

token = isfile.readline()
token = isfile.readline().split()
energy_shift = float(token[2])

aeng = []
while True:
  token = isfile.readline()
  if token == "":
    break
  else:
    token = token.split()
    aeng.append(float(token[1]))
isfile.close()

atomic_eng_sum = sum(aeng)

print("@@@ energy_shift                           %12.6f " % (energy_shift))
print("@@@ atomic_eng_sum                         %12.6f" % (atomic_eng_sum))
print("@@@ E_shift (energy_shift+atomic_eng_sum)  %12.6f" % (energy_shift+atomic_eng_sum))
print("@@@ E_tot                                  %12.6f" % (E_tot))
print("@@@ E_tot - E_shift =                      %12.6f" % (E_tot - (energy_shift+atomic_eng_sum)))


osfile = open("./energy_shift.dat", 'w')
osfile.write("%20.9f" % (energy_shift))
osfile.close()


