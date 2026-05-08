% RUN_ALL_PHASES  Master script for the multiscale measles model.
%
%   Executes Phases 1-5 sequentially, generating all figures (1-17)
%   and tables (A-F) needed for the paper.
%
%   Required files in working directory:
%     measles_params.m            - Parameter sets for 6 archetypes
%     measles_ode_rhs.m           - ODE right-hand side (4/7 var)
%     emergent_mortality.m        - Shared mortality function
%     phase1_validate_morris.m    - Core model validation
%     phase2_extended_model.m     - 7-variable extended model
%     phase3_immunotherapy.m      - DC/CTL therapy analysis
%     phase4_cellular_automaton.m - Population CA (uses emergent_mortality)
%     phase5_sensitivity_analysis.m - Sobol SA (uses emergent_mortality)
%
%   Outputs (all saved to ./export/):
%     Figures:  fig1-fig17  (.png, 300 dpi)
%     Data:     phase4_data.mat, phase5_sensitivity_data.mat
%     Tables:   Printed to console (copy for LaTeX)
%
%   Usage:
%     >> run_all_phases          % run everything
%     >> run_all_phases(3)       % run from Phase 3 onward
%     >> run_all_phases(4, 5)    % run only Phases 4 and 5
%
%   Estimated runtime: ~15-30 min total (Phase 5 is the slowest)

function run_all_phases(varargin)

if nargin == 0
    phases = 1:5;
elseif nargin == 1
    phases = varargin{1}:5;
else
    phases = [varargin{:}];
end

fprintf('\n');
fprintf('================================================================\n');
fprintf('  MULTISCALE MEASLES MODEL — MASTER RUNNER\n');
fprintf('  Phases to execute: %s\n', mat2str(phases));
fprintf('================================================================\n\n');

total_tic = tic;

if ~exist('export', 'dir'), mkdir('export'); end

% Verify all files exist
required = {'measles_params.m', 'measles_ode_rhs.m', 'emergent_mortality.m'};
for i = 1:numel(required)
    if ~exist(required{i}, 'file')
        error('Missing required file: %s', required{i});
    end
end

% ---- Phase 1: Core model validation ----
if ismember(1, phases)
    fprintf('\n====== PHASE 1 ======\n');
    phase_tic = tic;
    phase1_validate_morris();
    fprintf('  Phase 1 completed in %.1f s\n', toc(phase_tic));
end

% ---- Phase 2: Extended 7-variable model ----
if ismember(2, phases)
    fprintf('\n====== PHASE 2 ======\n');
    phase_tic = tic;
    phase2_extended_model();
    fprintf('  Phase 2 completed in %.1f s\n', toc(phase_tic));
end

% ---- Phase 3: Immunotherapy ----
if ismember(3, phases)
    fprintf('\n====== PHASE 3 ======\n');
    phase_tic = tic;
    phase3_immunotherapy();
    fprintf('  Phase 3 completed in %.1f s\n', toc(phase_tic));
end

% ---- Phase 4: Cellular automaton ----
if ismember(4, phases)
    fprintf('\n====== PHASE 4 ======\n');
    phase_tic = tic;
    phase4_cellular_automaton();
    fprintf('  Phase 4 completed in %.1f s\n', toc(phase_tic));
end

% ---- Phase 5: Sensitivity analysis ----
if ismember(5, phases)
    fprintf('\n====== PHASE 5 ======\n');
    phase_tic = tic;
    phase5_sensitivity_analysis('paper');
    phase5_antibody_stimulation_sensitivity();
    phase4_demographic_ic_sensitivity('paper');
    phase4_outbreak_plausibility_check();
    fprintf('  Phase 5 and v9 diagnostics completed in %.1f s\n', toc(phase_tic));
end

% ---- Summary ----
fprintf('\n================================================================\n');
fprintf('  ALL PHASES COMPLETE\n');
fprintf('  Total runtime: %.1f s (%.1f min)\n', toc(total_tic), toc(total_tic)/60);
fprintf('  Output directory: ./export/\n');
fprintf('================================================================\n');

% List generated files
files = dir(fullfile('export', '*'));
fprintf('\n  Generated files:\n');
for i = 1:numel(files)
    if ~files(i).isdir
        fprintf('    %-40s  %6.0f KB\n', files(i).name, files(i).bytes/1024);
    end
end
fprintf('\n');

end
