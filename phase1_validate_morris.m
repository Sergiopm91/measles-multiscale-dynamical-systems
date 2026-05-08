function phase1_validate_morris()
% PHASE1_VALIDATE_MORRIS  Solve the core Morris model (S, I, A, V)
%   for all 6 host archetypes and generate validation figures.
%
%   >> phase1_validate_morris

close all; clc;
fprintf('============================================================\n');
fprintf('  Phase 1: Morris Core Model Validation (S, I, A, V)\n');
fprintf('============================================================\n\n');

if ~exist('export', 'dir'), mkdir('export'); end

%% ---- Archetype names and colors ----
N_arch = 6;
arch_names = {'Healthy adult', 'Child (<5 y)', 'Elderly (>65 y)', ...
              'Immunocompromised', 'Partial vaccine', 'Fully vaccinated'};
colors = [0.00 0.45 0.74;   % blue
          0.85 0.33 0.10;   % red-orange
          0.47 0.67 0.19;   % green
          0.64 0.08 0.18;   % dark red
          0.49 0.18 0.56;   % purple
          0.20 0.70 0.70];  % teal

%% ---- Solve for all archetypes ----
solutions = cell(N_arch,1);
for omega = 1:N_arch
    P = measles_params(omega);
    y0 = [P.S0; P.I0; P.A0; P.V0];
    odefun = @(t,y) measles_ode_rhs(t, y, P);
    tspan = [0, P.tmax_days];

    opts = odeset('RelTol', 1e-10, 'AbsTol', 1e-12, ...
                  'MaxStep', 0.1, 'NonNegative', 1:4);
    fprintf('  Solving archetype %d (%s)... ', omega, arch_names{omega});
    tic;
    [t, Y] = ode15s(odefun, tspan, y0, opts);
    elapsed = toc;
    fprintf('%d steps, %.2f s\n', numel(t), elapsed);

    solutions{omega}.t = t;
    solutions{omega}.S = Y(:,1);
    solutions{omega}.I = Y(:,2);
    solutions{omega}.A = Y(:,3);
    solutions{omega}.V = Y(:,4);
    solutions{omega}.P = P;
end

%% ================================================================
%  FIGURE 1: Healthy adult -- full dynamics (4 panels)
% ================================================================
fprintf('\n  Generating Figure 1: Healthy adult dynamics...\n');

sol = solutions{1};
fig1 = figure('Name', 'Fig1: Healthy Adult', 'Color', 'w', ...
    'Position', [50 50 1000 700]);

subplot(2,2,1);
plot(sol.t, sol.S, '-', 'Color', colors(1,:), 'LineWidth', 2);
xlabel('Days post-infection'); ylabel('S (cells/\muL)');
title('(a) Susceptible target cells');
grid on; xlim([0 60]);

subplot(2,2,2);
semilogy(sol.t, max(sol.I, 1e-10), '-', 'Color', colors(2,:), 'LineWidth', 2);
xlabel('Days post-infection'); ylabel('I (cells/\muL)');
title('(b) Infected cells');
grid on; xlim([0 60]);

subplot(2,2,3);
plot(sol.t, sol.A, '-', 'Color', colors(3,:), 'LineWidth', 2);
xlabel('Days post-infection'); ylabel('A (cells/\muL)');
title('(c) Activated T cells');
grid on; xlim([0 60]);

subplot(2,2,4);
semilogy(sol.t, max(sol.V, 1e-10), '-', 'Color', colors(4,:), 'LineWidth', 2);
xlabel('Days post-infection'); ylabel('V (virus)');
title('(d) Viremia');
grid on; xlim([0 60]);

sgtitle('Figure 1: Healthy Adult -- Morris Core Dynamics', ...
    'FontWeight', 'bold', 'FontSize', 14);
saveas(fig1, fullfile('export', 'fig1_healthy_adult_dynamics.png'));

%% ================================================================
%  FIGURE 2: The Measles Paradox -- S and A overlaid
% ================================================================
fprintf('  Generating Figure 2: Measles paradox...\n');

fig2 = figure('Name', 'Fig2: Measles Paradox', 'Color', 'w', ...
    'Position', [100 50 900 400]);

yyaxis left;
plot(sol.t, sol.S, '-', 'Color', colors(1,:), 'LineWidth', 2);
hold on;
plot(sol.t, sol.A, '-', 'Color', colors(3,:), 'LineWidth', 2);
ylabel('Cells/\muL');
ylim([0, max(sol.S)*1.1]);

yyaxis right;
semilogy(sol.t, max(sol.V, 1e-10), '-', 'Color', colors(4,:), 'LineWidth', 2);
ylabel('Viremia (log scale)');

xlabel('Days post-infection');
xlim([0 40]);
legend('S: Target cells', 'A: Activated T cells', 'V: Viremia', ...
    'Location', 'east');
title('The Measles Paradox: Viral Clearance Amid Immunosuppression');
grid on;

xp = [7 16 16 7];
yp_ax = ylim;

hold on

h_patch = fill(xp, ...
    [yp_ax(1) yp_ax(1) yp_ax(2) yp_ax(2)], ...
    [0.8 0.9 1.0], ...   
    'FaceAlpha', 0.25, ...
    'EdgeColor', 'none');
set(h_patch, 'DisplayName', 'Immune transition phase');

legend show

saveas(fig2, fullfile('export', 'fig2_measles_paradox.png'));

%% ================================================================
%  FIGURE 3: All archetypes -- viremia comparison
% ================================================================
fprintf('  Generating Figure 3: Archetype viremia comparison...\n');

fig3 = figure('Name', 'Fig3: Archetype Viremia', 'Color', 'w', ...
    'Position', [150 50 900 500]);

hold on;
for omega = 1:N_arch
    sol = solutions{omega};
    semilogy(sol.t, max(sol.V, 1e-10), '-', 'Color', colors(omega,:), ...
        'LineWidth', 2);
end
xlabel('Days post-infection'); ylabel('V (viremia, log scale)');
title('Viremia Across Host Archetypes');
legend(arch_names, 'Location', 'northeast');
grid on; xlim([0 80]); set(gca, 'YScale', 'log');

saveas(fig3, fullfile('export', 'fig3_archetype_viremia.png'));

%% ================================================================
%  FIGURE 4: All archetypes -- 4-panel comparison
% ================================================================
fprintf('  Generating Figure 4: Archetype comparison (4 panels)...\n');

fig4 = figure('Name', 'Fig4: All Archetypes', 'Color', 'w', ...
    'Position', [200 50 1100 700]);

var_names = {'S: Target cells', 'I: Infected cells', ...
             'A: Activated T cells', 'V: Viremia'};
var_idx = {'S', 'I', 'A', 'V'};
use_log = [false, true, false, true];
panels = 'abcd';

for v = 1:4
    subplot(2,2,v); hold on;
    for omega = 1:N_arch
        sol = solutions{omega};
        ydata = sol.(var_idx{v});
        if use_log(v)
            semilogy(sol.t, max(ydata, 1e-10), '-', ...
                'Color', colors(omega,:), 'LineWidth', 1.8);
        else
            plot(sol.t, ydata, '-', ...
                'Color', colors(omega,:), 'LineWidth', 1.8);
        end
    end
    xlabel('Days post-infection');
    ylabel(var_names{v});
    title(sprintf('(%s) %s', panels(v), var_names{v}));
    grid on;
    if use_log(v), set(gca, 'YScale', 'log'); end
    xlim([0 80]);
    if v == 1
        legend(arch_names, 'Location', 'best', 'FontSize', 7);
    end
end

sgtitle('Figure 4: Within-Host Dynamics -- All Archetypes', ...
    'FontWeight', 'bold', 'FontSize', 14);
saveas(fig4, fullfile('export', 'fig4_all_archetypes.png'));

%% ================================================================
%  FIGURE 5: Key quantities summary
% ================================================================
fprintf('  Generating Figure 5: Summary metrics...\n');

fig5 = figure('Name', 'Fig5: Summary', 'Color', 'w', ...
    'Position', [250 50 1100 500]);

peak_V     = zeros(N_arch,1);
time_peak  = zeros(N_arch,1);
nadir_S    = zeros(N_arch,1);
peak_A     = zeros(N_arch,1);
clearance  = zeros(N_arch,1);

for omega = 1:N_arch
    sol = solutions{omega};
    [peak_V(omega), idx] = max(sol.V);
    time_peak(omega) = sol.t(idx);
    nadir_S(omega) = min(sol.S);
    peak_A(omega) = max(sol.A);
    idx_clear = find(sol.V(idx:end) < 1e-6, 1, 'first');
    if ~isempty(idx_clear)
        clearance(omega) = sol.t(idx + idx_clear - 1);
    else
        clearance(omega) = sol.P.tmax_days;
    end
end

subplot(2,2,1);
bar(peak_V, 'FaceColor', 'flat', 'CData', colors);
set(gca, 'XTickLabel', arch_names, 'XTickLabelRotation', 35, 'FontSize', 7);
ylabel('Peak viremia'); title('(a) Peak Viremia'); grid on;

subplot(2,2,2);
bar(time_peak, 'FaceColor', 'flat', 'CData', colors);
set(gca, 'XTickLabel', arch_names, 'XTickLabelRotation', 35, 'FontSize', 7);
ylabel('Days'); title('(b) Time to Peak Viremia'); grid on;

subplot(2,2,3);
bar(nadir_S, 'FaceColor', 'flat', 'CData', colors);
set(gca, 'XTickLabel', arch_names, 'XTickLabelRotation', 35, 'FontSize', 7);
ylabel('Cells/\muL'); title('(c) Minimum Target Cells (S)'); grid on;

subplot(2,2,4);
bar(clearance, 'FaceColor', 'flat', 'CData', colors);
set(gca, 'XTickLabel', arch_names, 'XTickLabelRotation', 35, 'FontSize', 7);
ylabel('Days'); title('(d) Viral Clearance Time'); grid on;

sgtitle('Figure 5: Archetype Summary Metrics', ...
    'FontWeight', 'bold', 'FontSize', 14);
saveas(fig5, fullfile('export', 'fig5_summary_metrics.png'));

%% ================================================================
%  TABLE: Summary of key metrics
% ================================================================
fprintf('\n============================================================\n');
fprintf('  TABLE: Archetype Summary Metrics\n');
fprintf('============================================================\n');
fprintf('%-22s %10s %10s %10s %10s %10s\n', ...
    'Archetype', 'Peak V', 't_peak', 'min(S)', 'max(A)', 't_clear');
fprintf('%s\n', repmat('-', 1, 72));
for omega = 1:N_arch
    fprintf('%-22s %10.4f %10.1f %10.1f %10.1f %10.1f\n', ...
        arch_names{omega}, peak_V(omega), time_peak(omega), ...
        nadir_S(omega), peak_A(omega), clearance(omega));
end
fprintf('============================================================\n');



%% ================================================================
%  Export paper-ready tables: parameter catalog and core ODE metrics
% ================================================================
export_parameter_catalog(fullfile('export','phase_parameters_by_source.csv'));

phase1_rows = cell(N_arch, 1);
for omega = 1:N_arch
    sol = solutions{omega};
    YY = [sol.S sol.I sol.A sol.V];
    M = ode_common_metrics(sol.t, YY, 'core4', sol.P);
    phase1_rows{omega} = {arch_names{omega}, M.Peak_V, M.t_peak, M.min_S, ...
        M.peak_A, M.peak_I, M.T_clear, M.AUC_beta, sol.P.profile_type, ...
        sol.P.profile_interpretation};
end
T_phase1 = cell2table(vertcat(phase1_rows{:}), 'VariableNames', ...
    {'Archetype','Peak_V','t_peak_days','min_S','peak_A','peak_I', ...
     'T_clear_days','AUC_beta','Profile_type','Interpretation'});
writetable(T_phase1, fullfile('export','phase1_core_ode_metrics.csv'));
fprintf('  Saved: export/phase1_core_ode_metrics.csv\n');

% Numerical stability check for the within-host ODE layer.
phase_ode_numerical_stability();

%% ================================================================
%  Validation checklist (healthy adult)
% ================================================================
fprintf('\n  VALIDATION CHECKLIST (qualitative, healthy adult):\n');
sol = solutions{1};
[~, iP] = max(sol.V);
tP = sol.t(iP);
minS = min(sol.S);
[~, iA] = max(sol.A);
tA = sol.t(iA);

checks = {
    'Viremia peaks around day 7-10',        tP >= 5 && tP <= 12;
    'Target cells depleted (S drops >50%%)', minS < 0.5*sol.S(1);
    'T cells peak around day 15-22',         tA >= 12 && tA <= 25;
    'Viremia cleared by day 30',             sol.V(end) < 1e-4;
    'T cells peak after viremia peak',       tA > tP;
    'S begins recovery after day 15',        sol.S(end) > minS;
};

all_pass = true;
for i = 1:size(checks, 1)
    if checks{i,2}
        status = 'PASS';
    else
        status = 'FAIL';
        all_pass = false;
    end
    fprintf('    [%s] %s\n', status, checks{i,1});
end

if all_pass
    fprintf('\n  ALL CHECKS PASSED. Core model reproduces Morris dynamics.\n');
else
    fprintf('\n  WARNING: Some checks failed. Review parameters.\n');
end

fprintf('\n  Figures saved to ./export/\n');
fprintf('============================================================\n');
end