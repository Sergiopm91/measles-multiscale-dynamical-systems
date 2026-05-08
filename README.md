# Measles multiscale dynamical-systems model

This repository contains the reproducible MATLAB implementation supporting the manuscript:

**A multiscale dynamical-systems model of measles immuno-epidemiology with ODE-to-cellular-automaton coupling**

The project implements a multiscale mathematical-biology framework for measles immuno-epidemiology. The framework couples within-host ordinary differential equation (ODE) dynamics with a stochastic spatial cellular automaton (CA), allowing host-level viral and immune trajectories to inform population-level transmission, severity and vaccination scenarios.

## Overview

The model integrates:

- a Morris-type four-variable within-host measles ODE core;
- a seven-variable ODE extension including:
  - susceptible target cells;
  - infected cells;
  - IFN-gamma-dominant cellular immune response;
  - IL-17-associated cellular immune response;
  - infectious virus;
  - persistent viral RNA;
  - neutralizing antibodies;
- six host archetypes:
  - healthy adult;
  - child;
  - elderly;
  - immunocompromised;
  - partially vaccinated;
  - strongly seroprotected breakthrough profile;
- in silico immune perturbation modules;
- an explicit ODE-to-cellular-automaton coupling map;
- an ODE-informed emergent mortality mapping;
- a stochastic spatial cellular automaton with vaccination and reactive-campaign scenarios;
- global Sobol sensitivity analysis and convergence diagnostics.

The repository is intended for scientific reproducibility and methodological transparency. It is not intended for clinical prediction, medical decision-making or therapeutic guidance.

## Repository structure

```text
.
├── README.md
├── LICENSE
├── CITATION.cff
├── .gitignore
├── manuscript/
│   └── main_measles_aims_v9.tex
├── src/
│   ├── measles_params.m
│   ├── measles_ode_rhs.m
│   ├── phase1_validate_morris.m
│   ├── phase2_extended_model.m
│   ├── phase3_immunotherapy.m
│   ├── phase4_cellular_automaton.m
│   ├── phase5_sensitivity_analysis.m
│   ├── emergent_mortality.m
│   └── run_all_phases.m
├── scripts/
│   ├── reproduce_all.m
│   ├── reproduce_figures.m
│   ├── reproduce_tables.m
│   └── check_outputs.m
├── export/
│   ├── figures/
│   ├── tables/
│   ├── metadata/
│   └── mat/
├── docs/
│   ├── model_overview.md
│   ├── ode_to_ca_coupling.md
│   ├── mortality_mapping.md
│   ├── sensitivity_analysis.md
│   └── reproducibility_notes.md
└── data/
    └── README.md
Main model components
Phase 1: Four-variable ODE core

The first phase validates a Morris-type within-host measles immunodynamics model with four state variables:

susceptible target cells;
infected cells;
activated T cells;
infectious virus.

This phase is used as the acute-infection reference model.

Phase 2: Seven-variable ODE extension

The second phase extends the within-host system to seven variables by adding:

IFN-gamma-dominant cellular immune activity;
IL-17-associated immune activity;
persistent viral RNA;
neutralizing antibodies.

The extension is designed to preserve the acute infection-clearance behavior of the four-variable core while adding post-acute and humoral-immunity descriptors.

Phase 3: In silico immune perturbations

The third phase evaluates mathematical immune-perturbation modules, including:

dendritic-cell-like stimulation;
CTL-like immune input;
combined DC-like and CTL-like perturbation;
diagnostic additive DC-like control.

These modules are exploratory mathematical perturbations and should not be interpreted as clinically validated measles therapies.

Phase 4: ODE-to-CA coupling and cellular automaton

The fourth phase maps ODE-derived descriptors into cellular-automaton parameters:

rash-centered infectious timing;
daily infectivity profiles;
cumulative infectivity burden;
emergent mortality probabilities.

The CA uses a spatial lattice with local Moore-neighborhood transmission, explicit vaccinated states, host archetypes and stochastic progression.

Phase 5: Sensitivity analysis

The fifth phase performs:

Tier 1 Sobol sensitivity analysis for within-host and multiscale outputs;
convergence assessment up to N_base = 4096;
Tier 2 coupling-map surrogate diagnostics;
functional sensitivity of antibody stimulation;
demographic sensitivity to the immunocompromised fraction.
Reproducibility

The full workflow can be reproduced from MATLAB by running:

run('src/run_all_phases.m')

or, if using the helper script:

run('scripts/reproduce_all.m')

The expected outputs are written to the export/ directory.

Main outputs

The workflow generates:

ODE trajectory figures;
host-archetype comparison figures;
in silico immune-perturbation figures;
cellular-automaton vaccination and campaign figures;
Sobol sensitivity figures;
mortality-map validation figures;
CSV summary tables;
simulation metadata;
MATLAB .mat files for selected phases.

Representative output files include:

export/figures/fig1_healthy_adult_dynamics.png
export/figures/fig7_archetypes_7var.png
export/figures/fig11_vax_sweep.png
export/figures/fig16_tier1_sobol_indices.png
export/figures/fig21_antibody_stimulation_sensitivity.png
export/tables/phase4_table_coverage_sweep.csv
export/tables/phase4_table_therapy_strategies.csv
export/tables/phase5_tier1_sobol_N4096.csv
export/tables/phase5_mortality_validation.csv
Interpretation and scope

This repository provides a computational framework for mathematical-biology research. The model is intended to study qualitative and semi-mechanistic links among measles within-host immunodynamics, host heterogeneity, vaccination coverage, spatial transmission and severity.

The framework does not provide:

clinical diagnosis;
clinical prognosis;
treatment recommendations;
vaccine-effectiveness estimates for a specific population;
outbreak forecasts for a specific region.

Any clinical or public-health interpretation would require additional calibration, validation and domain-specific review.

Software requirements

The code was developed for MATLAB.

Recommended MATLAB components:

MATLAB R2023a or later;
Statistics and Machine Learning Toolbox;
Parallel Computing Toolbox, optional but useful for large sensitivity analyses.

The scripts use standard MATLAB numerical solvers, including ode15s.

How to cite

If this repository is used, please cite the associated manuscript:

Pérez Montes, S., and Chimal Eguía, J. C.
A multiscale dynamical-systems model of measles immuno-epidemiology with ODE-to-cellular-automaton coupling.
Manuscript submitted to AIMS Mathematics.

After publication, this section should be updated with the DOI and final journal citation.

Authors
Sergio Pérez Montes
Laboratorio de Ciencias Matemáticas y Computacionales, Centro de Investigación en Computación, Instituto Politécnico Nacional, Ciudad de México, Mexico.
Juan Carlos Chimal Eguía
Laboratorio de Ciencias Matemáticas y Computacionales, Centro de Investigación en Computación, Instituto Politécnico Nacional, Ciudad de México, Mexico.
License

The source code is released under the MIT License unless otherwise specified.

Data availability

The repository contains generated CSV tables, figure files and metadata used to support the manuscript. No human-subject, animal, clinical-trial or cell-line data are included.

Disclaimer

This repository is provided for academic research in mathematical biology. It is not intended for clinical, diagnostic, therapeutic or public-health decision-making.
