function export_dashboard_json()
% EXPORT_DASHBOARD_JSON  Export MATLAB results as JSON for the React dashboard.
%
%   Run after run_all_phases (or at least Phases 4-5).
%   Reads: export/phase4_data.mat, export/phase5_sensitivity_data.mat
%   Writes: export/dashboard_data.json
%
%   >> export_dashboard_json

fprintf('  Exporting dashboard JSON...\n');

d4 = load(fullfile('export','phase4_data.mat'));

out = struct();

% Profiles (all 9: 6 baseline + 3 therapy)
for w = 1:numel(d4.profiles)
    p = d4.profiles(w);
    out.profiles(w).name = p.name;
    out.profiles(w).T_lat = p.T_lat;
    out.profiles(w).T_pre_rash = p.T_pre_rash;
    out.profiles(w).T_inf = p.T_inf;
    out.profiles(w).T_rash = p.T_rash;
    out.profiles(w).M_emergent = p.M_emergent;
    out.profiles(w).peak_V = p.peak_V;
    out.profiles(w).peak_Ag = p.peak_Ag;
    out.profiles(w).auc_beta = p.auc_beta;
end

% Vax sweep
out.vax_levels = d4.vax_levels;
out.ar_mean = d4.ar_mean;
out.ar_std = d4.ar_std;
out.d_mean = d4.d_mean;
out.d_std = d4.d_std;

% Per-archetype outcomes at each vax level
for vi = 1:numel(d4.vax_levels)
    out.arch_outcomes(vi).vax = d4.arch(vi).vax;
    out.arch_outcomes(vi).stats = d4.arch(vi).stats;
end

% Therapy results
for ti = 1:numel(d4.therapy_results)
    out.therapy(ti).arch_stats = d4.therapy_results(ti).arch_stats;
    out.therapy(ti).deaths_mean = mean(d4.therapy_results(ti).deaths);
    out.therapy(ti).ar_mean = mean(d4.therapy_results(ti).ar);
end

% Campaign results
for ci = 1:numel(d4.campaign_results)
    out.campaigns(ci).ar_mean = mean(d4.campaign_results(ci).ar);
    out.campaigns(ci).deaths_mean = mean(d4.campaign_results(ci).deaths);
end

% Sensitivity data (if available)
if exist(fullfile('export','phase5_sensitivity_data.mat'),'file')
    d5 = load(fullfile('export','phase5_sensitivity_data.mat'));
    out.sobol.S1 = d5.S1;
    out.sobol.ST = d5.ST;
    out.sobol.M_values = d5.M_values;
    out.sobol.severity_data = d5.severity_data;
    out.sobol.S1_coupling = d5.S1_coupling;
    out.sobol.ST_coupling = d5.ST_coupling;
    out.sobol.coupling_param_names = d5.coupling_param_names;
    out.sobol.coupling_qoi_names = d5.coupling_qoi_names;
    out.sobol.param_names = d5.param_names;
    out.sobol.qoi_names = d5.qoi_names;
end

json_str = jsonencode(out);
fid = fopen(fullfile('export','dashboard_data.json'), 'w');
fprintf(fid, '%s', json_str);
fclose(fid);
fprintf('  Saved: export/dashboard_data.json\n');
end
