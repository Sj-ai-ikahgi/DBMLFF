#!/bin/bash

  rm -rf ./data_set

  export PYTHONPATH="../../dataset":${PYTHONPATH}

  python dft_dbmlff_diff.py
  python dft_dbmlff_diff_single.py
  
