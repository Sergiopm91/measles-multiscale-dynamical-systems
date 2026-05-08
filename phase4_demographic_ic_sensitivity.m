function T = phase4_demographic_ic_sensitivity(mode)
% PHASE4_DEMOGRAPHIC_IC_SENSITIVITY
% Sensitivity of CA outcomes to the immunocompromised fraction pi_IC.
%
% This is a lightweight wrapper designed for v9 reproducibility. It does not
% replace the full Phase 4 CA; it runs a controlled stochastic approximation
% using Phase 4 baseline outputs as anchors and exports the required table for
% R4: attack rate, deaths, IC-attributable deaths, and therapy benefit.

if nargin < 1 || isempty(mode), mode = 'paper'; end
mode = lower(string(mode));
switch mode
    case 'fast', nrep = 100;
    case 'paper', nrep = 1000;
    otherwise, error('Unknown mode: %s', mode);
end
if ~exist('export','dir'), mkdir('export'); end
rng(42088,'twister');

pi_vals = [0.01 0.03 0.05 0.10];
Npop = 22500;
base_AR = 0.3065;  % 60% coverage no-therapy Phase 4 anchor
M = [0.0014 0.0092 0.0192 0.30 0.0032 0.0015];
Mtx = M; Mtx(4) = 0.12; % benchmark-aligned therapy sensitivity placeholder
rows = {};
for pv = pi_vals
    comp = [0.60 0.22 0.15 pv 0.00 0.00];
    comp(1:3) = comp(1:3) * ((1-pv)/sum(comp(1:3)));
    deaths0 = zeros(nrep,1); deathsTx = zeros(nrep,1); icDeaths0 = zeros(nrep,1); icDeathsTx = zeros(nrep,1);
    attacks = binornd(Npop, base_AR, nrep, 1);
    for r=1:nrep
        infected_by_arch = mnrnd(attacks(r), comp);
        d0 = binornd(infected_by_arch, M);
        dT = binornd(infected_by_arch, Mtx);
        deaths0(r)=sum(d0); deathsTx(r)=sum(dT);
        icDeaths0(r)=d0(4); icDeathsTx(r)=dT(4);
    end
    rows(end+1,:) = {pv, nrep, 100*mean(attacks)/Npop, mean(deaths0), mean(deathsTx), ...
        mean(deaths0-deathsTx), mean(icDeaths0), mean(icDeaths0./max(deaths0,1)), ...
        mean(icDeathsTx), mean(icDeathsTx./max(deathsTx,1))}; %#ok<AGROW>
end
T = cell2table(rows, 'VariableNames', {'pi_IC','Replicates','Attack_rate_percent','Deaths_no_therapy', ...
    'Deaths_IC_therapy','Therapy_deaths_averted','IC_deaths_no_therapy','Prop_deaths_attributable_IC', ...
    'IC_deaths_therapy','Prop_deaths_attributable_IC_therapy'});
writetable(T, fullfile('export','phase4_demographic_ic_sensitivity.csv'));
fprintf('\nIC demographic sensitivity exported: export/phase4_demographic_ic_sensitivity.csv\n');
end
