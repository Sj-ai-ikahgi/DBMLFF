import sys



def Get_Point_All(ifilename):

  sstring = len("0.50000E+00   0.34000E+01    ")

  isfile = open(ifilename, 'r')

  pnum_list = []
  cpoint_list = []

  while True:
    token = isfile.readline()
    if token == "":
      break
    else:
      token = token.split()
      pnum = int(token[0])
      for p in range(pnum):
        token = isfile.readline()

        if len(cpoint_list) > 2 and \
          token[sstring:] == cpoint_list[-1][sstring:] and \
          token[sstring:] == cpoint_list[-2][sstring:]:
          continue
        cpoint_list.append(token)
        pnum_list.append(pnum)

  isfile.close()

  return pnum_list, cpoint_list



def Write_Point(ofilename, pnum_list, cpoint_list):

  osfile = open(ofilename, 'w')

  for i in range(len(pnum_list)):
    osfile.write("%3d %5d \n" % (pnum_list[i], i+1))
    for p in range(pnum_list[i]):
      osfile.write(cpoint_list[i+p])



ifilename = sys.argv[1]

pnum_list, cpoint_list = Get_Point_All(ifilename)
Write_Point(ifilename + "_remove", pnum_list, cpoint_list)




