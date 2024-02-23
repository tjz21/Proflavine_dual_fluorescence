## Proflavine_dual_fluorescence


```bash
├── EOM_CCSD_input_data
├── MD_input_data
│   ├── acetonitrile
│   ├── acetonitrile_mm
│   ├── acetonitrile_qmmm
│   ├── methanol
│   ├── methanol_mm
│   ├── methanol_qmmm
│   └── vacuum
├── README.md
├── Results
│   ├── Acetonitrile
│   │   ├── QM_correction
│   │   ├── TDDFT
│   │   └── shifted
│   ├── MeOH
│   │   ├── EOM-CCSD
│   │   ├── QM_correction
│   │   ├── TDDFT
│   │   ├── shifted
│   │   └── stim_emission
│   ├── MeOH_stripped
│   │   └── QM_correction
│   ├── Vacuum
│   │   ├── EOM_CCSD
│   │   ├── TDDFT
│   │   ├── gap_scan
│   │   │   ├── 0.16
│   │   │   ├── 0.18
│   │   │   ├── 0.20
│   │   │   ├── 0.21
│   │   │   ├── 0.22
│   │   │   ├── 0.23
│   │   │   ├── 0.24
│   │   │   └── 0.26
│   │   ├── model_SDs
│   │   │   ├── 0.0005
│   │   │   ├── 0.0025
│   │   │   ├── 0.0075
│   │   │   ├── Debye
│   │   │   └── Debye_gaussian
│   │   ├── shifted
│   │   │   ├── fully_correlated
│   │   │   └── uncorrelated
│   │   └── stim_emission
│   └── vacuum_large_gap
├── TDDFT_input_data
│   ├── methanol
│   └── vacuum
├── chain_coefficients
│   ├── chain_coefficients
│   │   ├── chain_coeffs_MeOH.hdf5
│   │   ├── chain_coeffs_MeOH_stripped.hdf5
│   │   ├── chain_coeffs_acetonitrile.hdf5
│   │   ├── chain_coeffs_vacuum.hdf5
│   │   ├── chain_coeffs_vacuum_model_coupling_debye.hdf5
│   │   ├── chain_coeffs_vacuum_model_coupling_debye_and_gaussian.hdf5
│   │   ├── chain_coeffs_vacuum_model_coupling_lorentzian_0.0005.hdf5
│   │   ├── chain_coeffs_vacuum_model_coupling_lorentzian_0.0025.hdf5
│   │   └── chain_coeffs_vacuum_model_coupling_lorentzian_0.0075.hdf5
│   ├── python_scripts
│   │   ├── compute_diabatic_sds.py
│   │   └── spectral_dens_to_chain
│   ├── raw_diabatic_data
│   │   ├── full_energy_coupling_dipole_proflavine_camb3lyp_RPA_MeOH_diabatic_coupling.dat
│   │   ├── full_energy_coupling_dipole_proflavine_camb3lyp_RPA_MeOH_stripped_diabatic_coupling.dat
│   │   ├── full_energy_coupling_dipole_proflavine_camb3lyp_RPA_acetonitrile_diabatic_coupling.dat
│   │   └── full_energy_coupling_dipole_proflavine_camb3lyp_RPA_vacuum_diabatic_coupling.dat
│   └── spectral_densities
│       ├── acetonitrile
│       ├── methanol
│       ├── methanol_stripped
│       ├── vacuum
│       └── vacuum_model_sds
└── quantum_dynamics_input_data
    ├── example_fluorescence_input
    ├── hamiltonian_parameters.txt
    └── process_output
        ├── access_file.py
        └── dat_u45pW.jld

```
