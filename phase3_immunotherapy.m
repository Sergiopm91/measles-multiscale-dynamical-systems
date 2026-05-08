function phase3_immunotherapy()
% PHASE3_IMMUNOTHERAPY  In silico immunotherapy on immunocompromised host.
%
% Main scenarios:
%   1) No therapy
%   2) Multiplicative DC-like stimulation
%   3) CTL transfer
%   4) Combined DC + CTL
%
% A fifth diagnostic scenario tests an additive DC-control architecture.
% It is exported separately and should be described as an architecture
% control, not as the main therapy result and not propagated to the CA.
%
% Scientific interpretation:
%   - DC-only multiplicative stimulation is expected to overlap almost
%     completely with the no-therapy baseline in the immunocompromised
%     archetype, because the baseline effector substrate is nearly absent.
%   - To avoid visual ambiguity, the baseline curve is plotted last as a
%     black dashed line with larger line width.
%   - CTL-containing regimens are interpreted as in silico immune
%     perturbations, not as clinically validated measles therapies.

close all; clc;

fprintf('=== Phase 3: Immunotherapy on Immunocompromised Host ===\n\n');

if ~exist('export','dir')
    mkdir('export');
end

%% ------------------------------------------------------------------------
%  Scenario definitions
% -------------------------------------------------------------------------

scenarios = cell(5,1);

% Scenario 1: No therapy
P1 = measles_params(4);
scenarios{1} = make_scenario(P1, ...
    'No therapy', ...
    [0.00 0.00 0.00], ...
    '--', ...
    2.8);

% Scenario 2: Multiplicative DC-like stimulation
P2 = measles_params(4);
P2.DC_active   = true;
P2.DC_mode     = 'multiplicative';
P2.eta_D       = 0.5;
P2.eta_AbD     = 0.3;
P2.DC_doses    = 4;
P2.DC_start    = 12;
P2.DC_interval = 3;
P2.DC_d0       = 1.0;
P2.DC_ef       = 0.8;
P2.DC_muD      = 0.3;
scenarios{2} = make_scenario(P2, ...
    'DC therapy', ...
    [0.00 0.45 0.00], ...
    '-', ...
    2.0);

% Scenario 3: CTL transfer
P3 = measles_params(4);
P3.CTL_active  = true;
P3.CTL_delta   = 200;
P3.CTL_time    = 10;
P3.CTL_tau     = 1.5;
scenarios{3} = make_scenario(P3, ...
    'CTL transfer', ...
    [0.00 0.45 0.74], ...
    '-', ...
    2.0);

% Scenario 4: Combined DC + CTL
P4 = measles_params(4);
P4.DC_active   = true;
P4.DC_mode     = 'multiplicative';
P4.eta_D       = 0.5;
P4.eta_AbD     = 0.3;
P4.DC_doses    = 4;
P4.DC_start    = 12;
P4.DC_interval = 3;
P4.DC_d0       = 1.0;
P4.DC_ef       = 0.8;
P4.DC_muD      = 0.3;
P4.CTL_active  = true;
P4.CTL_delta   = 200;
P4.CTL_time    = 10;
P4.CTL_tau     = 1.5;
scenarios{4} = make_scenario(P4, ...
    'DC + CTL', ...
    [0.85 0.33 0.10], ...
    '-', ...
    2.0);

% Scenario 5: Architecture diagnostic only
% This intentionally asks whether the DC-only result is due to multiplicative
% scaling. Do not use as a main therapy scenario in the manuscript.
P5 = measles_params(4);
P5.DC_active    = true;
P5.DC_mode      = 'additive_control';
P5.DC_doses     = 4;
P5.DC_start     = 12;
P5.DC_interval  = 3;
P5.DC_d0        = 1.0;
P5.DC_ef        = 0.8;
P5.DC_muD       = 0.3;
P5.DC_add_q     = 0.15;
P5.DC_add_rhoAb = 0.006;
P5.DC_add_Ag    = 8.0;
scenarios{5} = make_scenario(P5, ...
    'DC additive control', ...
    [0.49 0.18 0.56], ...
    '--', ...
    2.0);

%% ------------------------------------------------------------------------
%  Solve all scenarios
% -------------------------------------------------------------------------

sols = cell(numel(scenarios),1);

for i = 1:numel(scenarios)

    P = scenarios{i}.P;
    y0 = [P.S0; P.I0; P.Ag0; P.A17_0; P.V0; P.R0; P.Ab0];

    opts = odeset( ...
        'RelTol', 1e-10, ...
        'AbsTol', 1e-12, ...
        'MaxStep', 0.1, ...
        'NonNegative', 1:7);

    fprintf('  Solving: %s... ', scenarios{i}.name);
    tic;
    [t, Y] = ode15s(@(t,y) measles_ode_rhs(t,y,P), ...
        [0 P.tmax_days], y0, opts);
    fprintf('%.2fs\n', toc);

    Ag  = Y(:,3);
    A17 = Y(:,4);
    V   = Y(:,5);
    R   = Y(:,6);
    Ab  = Y(:,7);

    Atot = Ag + A17;

    rash = Ag ./ (Ag + P.hR + eps);

    clearance_rate = P.c + P.theta_Ab .* Ab + P.k .* Atot;

    Kv = 0.1;
    beta_eff = V ./ (V + Kv);

    sols{i} = struct( ...
        't', t, ...
        'S', Y(:,1), ...
        'I', Y(:,2), ...
        'Ag', Ag, ...
        'A17', A17, ...
        'V', V, ...
        'R', R, ...
        'Ab', Ab, ...
        'Atot', Atot, ...
        'rash', rash, ...
        'clearance', clearance_rate, ...
        'beta_eff', beta_eff, ...
        'P', P);
end

%% ------------------------------------------------------------------------
%  Figure 9: Main therapy scenarios only
%  Plot order intentionally places no-therapy last so it remains visible.
% -------------------------------------------------------------------------

main_idx   = 1:4;
plot_order = [2 3 4 1];

fig9 = figure('Color','w','Position',[50 30 1200 800]);

panels = {'V','Ag','Ab','R','rash','beta_eff'};
titles = {'Viremia', ...
          'IFN\gamma T cells', ...
          'Antibodies', ...
          'Viral RNA', ...
          'Rash proxy', ...
          'Infectivity \beta_{eff}'};

ylabs  = {'V', ...
          'A_\gamma', ...
          'Ab', ...
          'R', ...
          'Rash', ...
          '\beta_{eff}'};

use_log = [true false false false false false];

for p = 1:6

    subplot(2,3,p);
    hold on;

    for kk = 1:numel(plot_order)

        i = plot_order(kk);
        yy = sols{i}.(panels{p});

        if use_log(p)
            semilogy(sols{i}.t, max(yy,1e-10), ...
                'Color', scenarios{i}.color, ...
                'LineStyle', scenarios{i}.style, ...
                'LineWidth', scenarios{i}.line_width);
        else
            plot(sols{i}.t, yy, ...
                'Color', scenarios{i}.color, ...
                'LineStyle', scenarios{i}.style, ...
                'LineWidth', scenarios{i}.line_width);
        end
    end

    xlabel('Days');
    ylabel(ylabs{p});
    title(titles{p});
    grid on;
    xlim([0 sols{1}.P.tmax_days]);

    if use_log(p)
        set(gca,'YScale','log');
    end

    if p == 1
        legend({'DC therapy', 'CTL transfer', 'DC + CTL', 'No therapy'}, ...
            'Location','best', ...
            'FontSize',7);
    end
end

sgtitle('Immunocompromised Host: Therapy Comparison', ...
    'FontWeight','bold', ...
    'FontSize',14);

exportgraphics(fig9, fullfile('export','fig9_therapy_comparison.png'), ...
    'Resolution',300);

savefig(fig9, fullfile('export','fig9_therapy_comparison.fig'));

%% ------------------------------------------------------------------------
%  Figure 10: Intervention timeline and clearance rate
% -------------------------------------------------------------------------

fig10 = figure('Color','w','Position',[100 50 1000 500]);

subplot(1,2,1);
hold on;

for kk = 1:numel(plot_order)

    i = plot_order(kk);

    semilogy(sols{i}.t, max(sols{i}.V,1e-10), ...
        'Color', scenarios{i}.color, ...
        'LineStyle', scenarios{i}.style, ...
        'LineWidth', scenarios{i}.line_width);
end

xline(10, ':', 'CTL pulse', ...
    'Color', [0.00 0.45 0.74], ...
    'LineWidth', 1.5, ...
    'LabelHorizontalAlignment','left', ...
    'LabelVerticalAlignment','bottom');

for d = 0:3
    xline(12 + d*3, ':', '', ...
        'Color', [0.00 0.45 0.00], ...
        'LineWidth', 1.0);
end

text(12, 1e-8, 'DC doses', ...
    'Color', [0.00 0.45 0.00], ...
    'FontSize', 8);

xlabel('Days');
ylabel('Viremia (log)');
title('Viremia + Intervention Timeline');
set(gca,'YScale','log');
grid on;
xlim([0 80]);

legend({'DC therapy', 'CTL transfer', 'DC + CTL', 'No therapy'}, ...
    'Location','best', ...
    'FontSize',7);

subplot(1,2,2);
hold on;

for kk = 1:numel(plot_order)

    i = plot_order(kk);

    plot(sols{i}.t, sols{i}.clearance, ...
        'Color', scenarios{i}.color, ...
        'LineStyle', scenarios{i}.style, ...
        'LineWidth', scenarios{i}.line_width);
end

xlabel('Days');
ylabel('Clearance rate (day^{-1})');
title('Effective Viral Clearance Rate');
grid on;
xlim([0 80]);

legend({'DC therapy', 'CTL transfer', 'DC + CTL', 'No therapy'}, ...
    'Location','best', ...
    'FontSize',7);

sgtitle('Immunotherapy Impact', ...
    'FontWeight','bold', ...
    'FontSize',13);

exportgraphics(fig10, fullfile('export','fig10_therapy_timeline.png'), ...
    'Resolution',300);

savefig(fig10, fullfile('export','fig10_therapy_timeline.fig'));

%% ------------------------------------------------------------------------
%  Figure 10b: Diagnostic DC architecture-control figure
% -------------------------------------------------------------------------

fig10b = figure('Color','w','Position',[150 60 1000 420]);

dc_diag_idx   = [1 2 5];
dc_plot_order = [2 5 1];

subplot(1,2,1);
hold on;

for kk = 1:numel(dc_plot_order)

    i = dc_plot_order(kk);

    semilogy(sols{i}.t, max(sols{i}.V,1e-10), ...
        'Color', scenarios{i}.color, ...
        'LineStyle', scenarios{i}.style, ...
        'LineWidth', scenarios{i}.line_width);
end

xlabel('Days');
ylabel('Viremia (log)');
title('DC architecture-control viremia');
set(gca,'YScale','log');
grid on;
xlim([0 100]);

legend({'DC therapy', 'DC additive control', 'No therapy'}, ...
    'Location','best', ...
    'FontSize',8);

subplot(1,2,2);
hold on;

for kk = 1:numel(dc_plot_order)

    i = dc_plot_order(kk);

    plot(sols{i}.t, sols{i}.Ag, ...
        'Color', scenarios{i}.color, ...
        'LineStyle', scenarios{i}.style, ...
        'LineWidth', scenarios{i}.line_width);
end

xlabel('Days');
ylabel('A_\gamma');
title('DC architecture-control IFN\gamma response');
grid on;
xlim([0 100]);

legend({'DC therapy', 'DC additive control', 'No therapy'}, ...
    'Location','best', ...
    'FontSize',8);

sgtitle('Diagnostic control: multiplicative versus additive DC-like stimulation', ...
    'FontWeight','bold');

exportgraphics(fig10b, fullfile('export','fig10b_dc_architecture_control.png'), ...
    'Resolution',300);

savefig(fig10b, fullfile('export','fig10b_dc_architecture_control.fig'));

%% ------------------------------------------------------------------------
%  Tables and exports
% -------------------------------------------------------------------------

rows = cell(numel(scenarios),1);

for i = 1:numel(scenarios)

    s = sols{i};

    YY = [s.S s.I s.Ag s.A17 s.V s.R s.Ab];

    M = ode_common_metrics(s.t, YY, 'extended7', s.P);

    rows{i} = { ...
        scenarios{i}.name, ...
        scenarios{i}.P.DC_mode, ...
        logical(scenarios{i}.P.DC_active), ...
        logical(scenarios{i}.P.CTL_active), ...
        M.Peak_V, ...
        M.t_peak, ...
        M.T_clear, ...
        M.Peak_Ag, ...
        M.Peak_A17, ...
        M.Peak_Ab, ...
        M.Peak_R, ...
        M.AUC_beta, ...
        M.AUC_Ag, ...
        M.AUC_A17, ...
        M.A17_Ag_AUC_ratio, ...
        M.t50_sigma17, ...
        M.t90_sigma17};
end

T = cell2table(vertcat(rows{:}), ...
    'VariableNames', { ...
    'Scenario', ...
    'DC_mode', ...
    'DC_active', ...
    'CTL_active', ...
    'Peak_V', ...
    't_peak_days', ...
    'T_clear_days', ...
    'Peak_Ag', ...
    'Peak_A17', ...
    'Peak_Ab', ...
    'Peak_R', ...
    'AUC_beta', ...
    'AUC_Ag', ...
    'AUC_A17', ...
    'AUC_A17_to_Ag_ratio', ...
    't50_sigma17_days', ...
    't90_sigma17_days'});

% Relative changes versus no-therapy baseline.
base = T(1,:);

T.Peak_V_change_percent = ...
    (T.Peak_V - base.Peak_V) ./ base.Peak_V * 100;

T.T_clear_change_days = ...
    T.T_clear_days - base.T_clear_days;

T.AUC_beta_change_percent = ...
    (T.AUC_beta - base.AUC_beta) ./ base.AUC_beta * 100;

T.Architecture_role = repmat( ...
    "main multiplicative therapy scenario", ...
    height(T), 1);

T.Architecture_role(strcmp(T.Scenario,'DC additive control')) = ...
    "diagnostic additive-control scenario; not propagated to CA";

writetable(T, fullfile('export','phase3_therapy_metrics.csv'));

%% ------------------------------------------------------------------------
%  Console report
% -------------------------------------------------------------------------

fprintf('\n  ============================================================\n');
fprintf('  TABLE: Therapy Impact on Immunocompromised Host\n');
fprintf('  ============================================================\n');

fprintf('  %-22s %8s %8s %8s %8s %10s %10s\n', ...
    'Scenario', 'peakV', 't_clear', 'peakAg', 'peakAb', 'peakR', 'AUC_beta');

fprintf('  %s\n', repmat('-',1,84));

for i = 1:height(T)
    fprintf('  %-22s %8.4f %8.1f %8.1f %8.3f %10.2f %10.2f\n', ...
        T.Scenario{i}, ...
        T.Peak_V(i), ...
        T.T_clear_days(i), ...
        T.Peak_Ag(i), ...
        T.Peak_Ab(i), ...
        T.Peak_R(i), ...
        T.AUC_beta(i));
end

fprintf('\n  Relative to no-therapy baseline:\n');

for i = 2:height(T)

    fprintf('  %s:\n', T.Scenario{i});
    fprintf('    Peak viremia:   %+.1f%%\n', T.Peak_V_change_percent(i));
    fprintf('    Clearance time: %+.1f days\n', T.T_clear_change_days(i));
    fprintf('    AUC(beta_eff):  %+.1f%%\n', T.AUC_beta_change_percent(i));
end

fprintf('\n  Saved: export/phase3_therapy_metrics.csv\n');
fprintf('  Figures -> ./export/fig9, fig10, fig10b\n');
fprintf('  NOTE: DC-only multiplicative results should be interpreted as a consequence\n');
fprintf('        of the implemented architecture, not as a validated clinical law.\n');
fprintf('  NOTE: The no-therapy baseline is plotted last as a black dashed curve so\n');
fprintf('        overlap with DC-only remains visually explicit.\n');
fprintf('  ============================================================\n');

end

%% =========================================================================
%  Local helper
% =========================================================================

function S = make_scenario(P, name, color, style, line_width)

if nargin < 5
    line_width = 2.0;
end

S = struct( ...
    'P', P, ...
    'name', name, ...
    'color', color, ...
    'style', style, ...
    'line_width', line_width);

end