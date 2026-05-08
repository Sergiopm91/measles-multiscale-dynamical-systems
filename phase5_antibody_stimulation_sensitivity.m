function Tout = phase5_antibody_stimulation_sensitivity()
% PHASE5_ANTIBODY_STIMULATION_SENSITIVITY
% Compares antibody stimulation functions in Eq. (4g): log(1+R),
% Michaelis--Menten, and Hill. This directly addresses B1.

close all;
if ~exist('export','dir'), mkdir('export'); end
forms = {'log','michaelis_menten','hill'};
labels = {'log(1+R)','R/(K_R+R)','Hill'};
rows = {};

for fi = 1:numel(forms)
    for w = 1:6
        P = measles_params(w);
        P.Ab_stim_type = forms{fi};
        P.K_R = 100;
        P.n_R = 2;
        [t,Y] = solve7(P);
        V = Y(:,5); Ag = Y(:,3); R = Y(:,6); Ab = Y(:,7);
        beta_eff = V./(V+0.1);
        peakV = max(V); [~,ix] = max(V); tpeak = t(ix);
        peakAg = max(Ag); peakR = max(R); peakAb = max(Ab);
        aucb = trapz(t,beta_eff);
        idx_clear = find(t >= 5 & V < 0.01, 1, 'first');
        if isempty(idx_clear), tclear = P.tmax_days; else, tclear = t(idx_clear); end
        [M,~] = emergent_mortality(peakV, min(15,max(8,round(tclear)-8)), peakAg, aucb, mortality_config('benchmark_aligned'));
        rows(end+1,:) = {string(forms{fi}), string(labels{fi}), string(P.name), peakV, tpeak, tclear, peakAg, peakR, peakAb, aucb, M}; %#ok<AGROW>
    end
end

Tout = cell2table(rows, 'VariableNames', {'Form','Label','Archetype','Peak_V','t_peak','t_clear','Peak_Ag','Peak_R','Peak_Ab','AUC_beta','M_emergent'});
writetable(Tout, fullfile('export','phase5_antibody_stimulation_sensitivity.csv'));

% Compact plot: normalized deviations relative to log form.
try
    base = Tout(strcmp(Tout.Form,'log'),:);
    metrics = {'Peak_V','t_clear','Peak_Ab','AUC_beta','M_emergent'};
    F = figure('Color','w','Position',[100 100 1150 520]);
    tiledlayout(1,numel(metrics),'TileSpacing','compact','Padding','compact');
    for mi = 1:numel(metrics)
        nexttile; hold on; grid on;
        vals = zeros(numel(forms),6);
        for fi = 1:numel(forms)
            sub = Tout(strcmp(Tout.Form,forms{fi}),:);
            vals(fi,:) = sub.(metrics{mi}) ./ max(base.(metrics{mi}), eps);
        end
        bar(vals');
        title(strrep(metrics{mi},'_','\_'));
        xticks(1:6); xticklabels({'H','Ch','E','IC','PV','FV'}); xtickangle(45);
        if mi == 1, ylabel('Ratio vs log(1+R)'); end
        if mi == numel(metrics), legend(labels,'Location','best'); end
    end
    exportgraphics(F, fullfile('export','fig21_antibody_stimulation_sensitivity.png'), 'Resolution',300);
catch ME
    warning('Could not create antibody stimulation sensitivity figure: %s', ME.message);
end

fprintf('\nAntibody stimulation sensitivity exported:\n');
fprintf('  export/phase5_antibody_stimulation_sensitivity.csv\n');
fprintf('  export/fig21_antibody_stimulation_sensitivity.png\n');
end

function [t,Y] = solve7(P)
y0 = [P.S0; P.I0; P.Ag0; P.A17_0; P.V0; P.R0; P.Ab0];
opts = odeset('RelTol',1e-10,'AbsTol',1e-12,'MaxStep',0.1,'NonNegative',1:7);
[t,Y] = ode15s(@(t,y) measles_ode_rhs(t,y,P), [0 P.tmax_days], y0, opts);
end
