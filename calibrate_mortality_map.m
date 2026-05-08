function T = calibrate_mortality_map(profiles, map_names)
% CALIBRATE_MORTALITY_MAP  Export mortality-map calibration diagnostics.
%
%   T = calibrate_mortality_map(profiles)
%
% If profiles is omitted, Phase 4 profiles are recomputed through the local
% ODE-to-CA extraction logic when available. This function compares named
% maps against internal benchmark anchors and exports absolute, relative,
% and log-scale errors.

if nargin < 2 || isempty(map_names)
    map_names = {'conservative','benchmark_aligned','high_risk_IC'};
end
if nargin < 1 || isempty(profiles)
    if exist('extract_coupled_profiles','file') == 2
        profiles = extract_coupled_profiles();
    else
        error('Provide profiles from phase4_cellular_automaton or run this after Phase 4.');
    end
end

bench = benchmark_cfr_table();
rows = {};
for mi = 1:numel(map_names)
    cfg = mortality_config(map_names{mi});
    for w = 1:6
        p = profiles(w);
        [M, d] = emergent_mortality(p.peak_V, p.T_infectious_duration, p.peak_Ag, p.auc_beta, cfg);
        b = bench.Benchmark_CFR(w);
        rows(end+1,:) = {string(cfg.map_name), string(p.name), p.peak_V, p.T_infectious_duration, ...
            p.peak_Ag, p.auc_beta, d.Z, M, b, M-b, (M-b)/max(b,eps), ...
            log10(max(M,1e-12))-log10(max(b,1e-12)), d.viral_burden, ...
            d.infectious_duration, d.immune_deficit, d.infectivity_burden, ...
            d.immune_infectivity_interaction}; %#ok<AGROW>
    end
end

T = cell2table(rows, 'VariableNames', {'Map','Archetype','Peak_V','T_inf','Peak_Ag','AUC_beta', ...
    'Z_score','M_predicted','Benchmark_CFR','Absolute_error','Relative_error','Log10_error', ...
    'Z_viral','Z_duration','Z_immune_deficit','Z_AUC','Z_deficit_AUC_interaction'});
if ~exist('export','dir'), mkdir('export'); end
writetable(T, fullfile('export','phase5_mortality_map_calibration_scenarios.csv'));

fprintf('\nMortality-map calibration scenarios exported: export/phase5_mortality_map_calibration_scenarios.csv\n');
end

function T = benchmark_cfr_table()
Archetype = {'Healthy adult';'Child';'Elderly';'Immunocomp.';'Partial vacc';'Full vacc'};
Benchmark_CFR = [0.002; 0.010; 0.015; 0.300; 0.0005; 0.0001];
T = table(Archetype, Benchmark_CFR);
end
