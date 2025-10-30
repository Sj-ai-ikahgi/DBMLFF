#!/bin/bash

  pDir=../MLFF_dataset/parameters

  cp -p ${pDir}/bond.molecule ./template
  cp -p ${pDir}/MOVEMENT.type ./template

  cp -p ${pDir}/find_neighbore.in ./template/input
  cp -p ${pDir}/IN.NEIGHBORE ./template/input

