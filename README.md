# DBMLFF Tutorial: Building Force Field for EC Molecule

## Overview
This tutorial provides a comprehensive guide to constructing the Density-Based Machine Learning Force Field (DBMLFF) for EC molecules. The DBMLFF framework consists of three key components: **MLFF**, **DB**, and **POL**, which collectively describe molecular interactions with high accuracy and computational efficiency.

## Prerequisites
- Access to PWmat for DFT calculations
- Basic familiarity with molecular dynamics simulations
- Computational resources for quantum chemistry calculations

## Step-by-Step Procedure

### 1. Code Compilation
**Working Directory:** `train/dataset/md_src`

```bash
# Load compilation environment
source env.sh

# Compile programs
make
```

### 2. Parameter File Generation
**Working Directory:** `train/dataset/generate_file`

#### 2.1 Parameter Setup
Modify `aux_dbmlff.py` with molecular-specific parameters:

```python
# Example parameter settings
MolEindexList = [element_indices]
MolAtomSysList = [atomic_symmetry_info]  
atomicMassList = [atomic_masses]
zValenceList = [valence_charges]
moveStartAtomId = [reference_atom_index]
```

#### 2.2 Structure Optimization
- Optimize target molecular structure in `reference_conf/`
- Save optimized structure as `atom_EC.config`

#### 2.3 Bond Parameterization
```bash
# Update molType in qsub scripts
sed -i 's/molType=.*/molType=EC/g' qsub_1.sh qsub_2.sh

# Execute parameter generation
bash qsub_1.sh

# Modify bond file and rerun
vim ./result/bond.molecule_EC_need_modify  # Add bond types and counts
mv ./result/bond.molecule_EC_need_modify ./result/bond.molecule_EC_bt
bash qsub_2.sh
```

**Output Location:** `train/dataset/generate_file/results`

### 3. DB Component: Monomer Charge Distribution Fitting
**Working Directory:** `train/DB` (Requires MLFF Step 2 completion)

```bash
# Step 3.1: Configuration Selection
cd DB_dataset/
bash choose_conf.sh  # Randomly selects 50 EC configurations

# Step 3.2: DFT Charge Density Calculation  
cd ../DB/
bash qsub_1.sh

# Step 3.3: Charge Density Fitting
# Edit Z-value in cp_rho_and_prepare_input.sh
vim 2_fit_DB/cp_rho_and_prepare_input.sh  # Set total valence electrons
bash qsub_2.sh

# Step 3.4: Validation
bash qsub_3.sh  # Tests DB convergence
```

**Quality Checks:**
- Check `DB/2_fit_DB/Smoothen_funcr/data_out` for smooth curves
- Verify convergence in `DB/3_check/egg_box/sample/data_out` at n123 ≈ 100

### 4. MLFF Component: Linear Model for DFT-DB Differences
**Working Directory:** `train/MLFF`

#### 4.1 AIMD Data Preparation
**Multi-molecule System:**
- 15 EMC molecules at experimental density
- AIMD: 400K, 700K, 1200K (5000 steps, 0.5 fs, NVT)

**Monomer Simulations:**
- 3 selected molecules
- Single-molecule AIMD: 400K, 600K, 1900K

#### 4.2 Configuration Collection
```bash
cd MLFF_dataset/
# Configure parameters in get_conf.sh
vim get_conf.sh  # Set temperature points and sampling parameters
sbatch get_conf.sh  # Generates Conf_* directories
```

#### 4.3 SCF Calculations
```bash
cd ../1_SCF/

# SCF calculations
bash qsub_1.1_scf.sh

# Extract energies and forces  
bash qsub_1.2_get_dft_eng_force.sh

# Energy correction
bash qsub_2_energy_shift.sh
```

#### 4.4 DB Calculations and Model Fitting
```bash
# DB energies/forces
cd ../2_dbmlff/
bash qsub.sh

# Model fitting
cd ../3_fit_mlff/
bash get_dataset.sh  # Prepares dataset
sbatch run.sh        # Fits MLFF model
```

### 5. POL Component: Intermolecular Polarization
**Working Directory:** `train/POL`

```bash
# Complete polarization workflow
bash qsub_0_1.sh    # Structure relaxation
bash qsub_2.sh      # Point charge generation
bash qsub_3.1.sh    # Task splitting
bash qsub_3.2.sh    # DFT calculations  
bash qsub_3.3.sh    # Result integration
bash qsub_4_5.sh    # Format conversion and fitting
bash qsub_6_7.sh    # Monomer polarization
```

### 6. Molecular Dynamics Simulation
**Working Directory:** `MD`

```bash
# Automated DBMLFF construction and MD
bash qsub.sh
```

## Required Files for MD Simulation

### Force Field File
- `mol.1`: Automatically generated from previous steps

### Control File (`MD.input`)
```
[Basic Parameters]
md_type = 1
md_step = 10000
...

[Force Field Settings]
ff_type = dbmlff
...
```

### Structure File (`atom.config`)
```
Total_Atoms  Molecule_Types  Atoms_per_Molecule  Molecule_Counts
[Atomic coordinates...]
```

### Supporting Files
- Pseudopotentials: `*.nc.pbe.UPF.ionrhoR` (provided)
- Submission script: `run_mpich.sh` (provided)

## Output Structure
```
MD/md_test/
├── mol.1                    # Complete force field
├── MD.input                 # Simulation parameters
├── atom.config              # Initial structure
├── *.nc.pbe.UPF.ionrhoR     # Pseudopotentials
├── run_mpich.sh             # Submission script
└── output/                  # Simulation results
```

## Validation Metrics
- Charge density fitting convergence (DB)
- MLFF model accuracy vs DFT reference
- Polarization energy convergence
- MD simulation stability

## Troubleshooting
1. **Compilation errors**: Check environment variables in `env.sh`
2. **SCF convergence**: Adjust DFT parameters in PWmat input
3. **Force field instability**: Verify training data quality and coverage
4. **MD crashes**: Check structure file format and force field parameters

This tutorial provides a complete workflow for constructing accurate machine learning force fields. For additional support, consult the DBMLFF documentation and user forums.
