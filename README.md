DBMLFF: Linear scaling machine learning force fields via electron density decomposition for molecular electrolytes
=
We introduce a density-based machine learning force field (DBMLFF) that calculates intermolecular interactions based on fitted electron charge densities, while intramolecular energies are modeled using traditional atomic structure-based MLFF models. In addition, a local potential polarization model based on charge density calculations is introduced to describe polarization interactions between molecules. The main advantage of DBMLFF is that the form of its intermolecular interaction energies is physical based and only separate models and parameters need to be developed for each molecule, thus simplifying the potential field development process. For simulations of electrolyte systems, DBMLFF can accurately characterize the intermolecular interactions and is computationally three orders of magnitude faster than density functional theory in molecular dynamics.
Files
=
**AIMD** Contains AIMD simulation raw data of the LiFSI/EA/TTE system. It includes key simulation results such as trajectory files, atomic configurations, and control files.

**DBMLFF_LiFSI_EGDME_TTE** Stores MD trajectory files of the LiFSI/EGDME/TTE system calculated based on the DBMLFF. These files can be used to analyze the microscopic dynamic behaviors of the electrolyte system, such as ion migration and solvation structure.

**DBMLFF_LiPF6_EC_DMC** Stores MD trajectory files of the LiPF6/EC/DMC system calculated based on the DBMLFF. 

**model** Holds core module files of the DBMLFF force field, including force field parameters, feature calculation configurations, and fitting models for multi-component molecules (e.g., solvent molecules, lithium salt ions). It serves as the fundamental dependency for implementing DBMLFF force field simulations.

