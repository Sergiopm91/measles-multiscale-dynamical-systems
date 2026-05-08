function T = phase5_sobol_convergence()
% PHASE5_SOBOL_CONVERGENCE  Driver for Tier-1 Sobol convergence runs.
%
% This script runs phase5_sensitivity_analysis at increasing N_base by setting
% the environment variable MEASLES_TIER1_N_BASE before each call. The patched
% phase5_sensitivity_analysis.m reads this variable when mode='custom'.

Ns = [512 1024 2048 4096];
if ~exist('export','dir'), mkdir('export'); end
rows = {};
for i=1:numel(Ns)
    setenv('MEASLES_TIER1_N_BASE', num2str(Ns(i)));
    fprintf('\n=== Sobol convergence run: N_base=%d ===\n', Ns(i));
    phase5_sensitivity_analysis('custom');
    src = fullfile('export','phase5_tier1_sobol.csv');
    dst = fullfile('export',sprintf('phase5_tier1_sobol_N%d.csv',Ns(i)));
    if exist(src,'file'), copyfile(src,dst); end
    rows(end+1,:) = {Ns(i), string(dst)}; %#ok<AGROW>
end
T = cell2table(rows,'VariableNames',{'N_base','Sobol_file'});
writetable(T, fullfile('export','phase5_sobol_convergence_manifest.csv'));
fprintf('\nSobol convergence manifest exported: export/phase5_sobol_convergence_manifest.csv\n');
end
