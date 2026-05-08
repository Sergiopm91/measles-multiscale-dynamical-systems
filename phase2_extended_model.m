function phase2_extended_model()
% PHASE2_EXTENDED_MODEL  7-variable model: S, I, Ag, A17, V, R, Ab
%   with shift gate sigma17, rash proxy, clearance rate.
%   Generates figures and paper-ready tables for all 6 archetypes.
%
%   >> phase2_extended_model

close all; clc;
fprintf('=== Phase 2: Extended Model (7 variables) ===\n\n');
if ~exist('export','dir'), mkdir('export'); end

names = {'Healthy adult','Child (<5 y)','Elderly (>65 y)', ...
         'Immunocompromised','Partial vaccine', 'Fully vaccinated'};
col = [0 .45 .74; .85 .33 .10; .47 .67 .19; .64 .08 .18; .49 .18 .56; 0.2, 0.2, 0.2];

sols = cell(6,1);
for w = 1:6
    P = measles_params(w);
    y0 = [P.S0; P.I0; P.Ag0; P.A17_0; P.V0; P.R0; P.Ab0];
    opts = odeset('RelTol',1e-10,'AbsTol',1e-12,'MaxStep',0.1,'NonNegative',1:7);
    fprintf('  Archetype %d (%s)... ', w, names{w});
    tic;
    [t, Y] = ode15s(@(t,y) measles_ode_rhs(t,y,P), [0 P.tmax_days], y0, opts);
    fprintf('%d steps, %.2fs\n', numel(t), toc);

    % Compute derived quantities
    Ag = Y(:,3); A17 = Y(:,4); V = Y(:,5); R = Y(:,6); Ab = Y(:,7);
    Atot = Ag + A17;

    % Shift gate sigma17
    sig17 = zeros(size(t));
    for i = 1:numel(t)
        z = P.alpha_t*(t(i)-P.t17) + P.alpha_R*log(1+R(i));
        sig17(i) = 1/(1+exp(-z));
    end

    % Rash proxy
    rash = Ag ./ (Ag + P.hR);

    % Clearance rate
    clearance_rate = P.c + P.theta_Ab*Ab + P.k*Atot;

    sols{w} = struct('t',t, 'S',Y(:,1), 'I',Y(:,2), 'Ag',Ag, ...
        'A17',A17, 'V',V, 'R',R, 'Ab',Ab, 'Atot',Atot, ...
        'sig17',sig17, 'rash',rash, 'clearance',clearance_rate, 'P',P);
end

%% ================================================================
%  FIGURE 6: Healthy adult — complete 7-variable dynamics
% ================================================================
fprintf('\n  Generating figures...\n');
s = sols{1};

fig6 = figure('Color','w','Position',[50 30 1200 900]);

subplot(3,3,1);
plot(s.t, s.S, '-', 'Color', col(1,:), 'LineWidth', 2);
xlabel('Days'); ylabel('cells/\muL'); title('S: Target cells');
grid on; xlim([0 120]);

subplot(3,3,2);
semilogy(s.t, max(s.I,1e-10), '-', 'Color', col(2,:), 'LineWidth', 2);
xlabel('Days'); ylabel('cells/\muL'); title('I: Infected cells');
grid on; xlim([0 120]);

subplot(3,3,3);
semilogy(s.t, max(s.V,1e-10), '-', 'Color', [.64 .08 .18], 'LineWidth', 2);
xlabel('Days'); ylabel('virus'); title('V: Viremia');
grid on; xlim([0 120]);

subplot(3,3,4);
plot(s.t, s.Ag, '-', 'Color', [0 .6 .3], 'LineWidth', 2); hold on;
plot(s.t, s.A17, '--', 'Color', [.85 .33 .10], 'LineWidth', 2);
xlabel('Days'); ylabel('cells/\muL'); title('T cells: IFN\gamma vs IL-17');
legend('A_\gamma (IFN\gamma)', 'A_{17} (IL-17)', 'Location', 'best');
grid on; xlim([0 120]);

subplot(3,3,5);
plot(s.t, s.sig17, '-', 'Color', [.49 .18 .56], 'LineWidth', 2);
xlabel('Days'); ylabel('\sigma_{17}'); title('Shift gate IFN\gamma \rightarrow IL-17');
ylim([0 1]); grid on; xlim([0 120]);

subplot(3,3,6);
plot(s.t, s.R, '-', 'Color', [.64 .08 .18], 'LineWidth', 2);
xlabel('Days'); ylabel('RNA'); title('R: Persistent viral RNA');
grid on; xlim([0 120]);

subplot(3,3,7);
plot(s.t, s.Ab, '-', 'Color', [0 .45 .74], 'LineWidth', 2);
xlabel('Days'); ylabel('Ab titer'); title('Ab: Neutralizing antibodies');
grid on; xlim([0 120]);

subplot(3,3,8);
plot(s.t, s.rash, '-', 'Color', [.85 .2 .2], 'LineWidth', 2);
xlabel('Days'); ylabel('Rash proxy'); title('Rash proxy');
ylim([0 1]); grid on; xlim([0 120]);

subplot(3,3,9);
plot(s.t, s.clearance, '-', 'Color', [.3 .3 .3], 'LineWidth', 2);
xlabel('Days'); ylabel('day^{-1}'); title('Effective clearance rate');
grid on; xlim([0 120]);

sgtitle('Healthy Adult — Complete 7-Variable Dynamics', ...
    'FontWeight','bold','FontSize',14);
saveas(fig6, fullfile('export','fig6_healthy_7var.png'));

%% ================================================================
%  FIGURE 7: Key variables across all archetypes (6-panel)
% ================================================================
fig7 = figure('Color','w','Position',[100 30 1200 800]);

% Viremia
subplot(2,3,1); hold on;
for w=1:6, semilogy(sols{w}.t, max(sols{w}.V,1e-10),'-','Color',col(w,:),'LineWidth',1.8); end
xlabel('Days'); ylabel('V'); title('Viremia'); grid on; set(gca,'YScale','log'); xlim([0 100]);
legend(names,'Location','best','FontSize',6);

% IFN-gamma T cells
subplot(2,3,2); hold on;
for w=1:6, plot(sols{w}.t, sols{w}.Ag,'-','Color',col(w,:),'LineWidth',1.8); end
xlabel('Days'); ylabel('A_\gamma'); title('IFN\gamma T cells'); grid on; xlim([0 100]);

% IL-17 T cells
subplot(2,3,3); hold on;
for w=1:6, plot(sols{w}.t, sols{w}.A17,'-','Color',col(w,:),'LineWidth',1.8); end
xlabel('Days'); ylabel('A_{17}'); title('IL-17 T cells'); grid on; xlim([0 100]);

% Persistent RNA
subplot(2,3,4); hold on;
for w=1:6, plot(sols{w}.t, sols{w}.R,'-','Color',col(w,:),'LineWidth',1.8); end
xlabel('Days'); ylabel('R'); title('Persistent viral RNA'); grid on; xlim([0 120]);

% Antibodies
subplot(2,3,5); hold on;
for w=1:6, plot(sols{w}.t, sols{w}.Ab,'-','Color',col(w,:),'LineWidth',1.8); end
xlabel('Days'); ylabel('Ab'); title('Neutralizing antibodies'); grid on; xlim([0 120]);

% Rash proxy
subplot(2,3,6); hold on;
for w=1:6, plot(sols{w}.t, sols{w}.rash,'-','Color',col(w,:),'LineWidth',1.8); end
xlabel('Days'); ylabel('Rash'); title('Rash proxy'); grid on; xlim([0 60]); ylim([0 1]);

sgtitle('Extended Model — All Archetypes','FontWeight','bold','FontSize',14);
saveas(fig7, fullfile('export','fig7_archetypes_7var.png'));

%% ================================================================
%  FIGURE 8: Shift gate sigma17 for all archetypes
% ================================================================
fig8 = figure('Color','w','Position',[150 50 700 400]);
hold on;
for w=1:6, plot(sols{w}.t, sols{w}.sig17,'-','Color',col(w,:),'LineWidth',2); end
xlabel('Days post-infection'); ylabel('\sigma_{17}(t)');
title('IFN\gamma \rightarrow IL-17 Shift Gate'); ylim([0 1]);
legend(names,'Location','best'); grid on; xlim([0 120]);
saveas(fig8, fullfile('export','fig8_shift_gate.png'));

%% ================================================================
%  Backward compatibility check
% ================================================================
fprintf('\n  Backward compatibility check (7-var vs 4-var, healthy adult):\n');
P = measles_params(1);
% Solve 4-var
y0_4 = [P.S0; P.I0; P.A0; P.V0];
opts = odeset('RelTol',1e-10,'AbsTol',1e-12,'MaxStep',0.1,'NonNegative',1:4);
[t4, Y4] = ode15s(@(t,y) measles_ode_rhs(t,y,P), [0 P.tmax_days], y0_4, opts);

% Solve 7-var with extensions zeroed out
P0 = P; P0.kappa17 = 0; P0.theta_Ab = 0; P0.rho_R = 0; P0.rho_Ab = 0;
P0.alpha_t = 0; P0.alpha_R = 0; P0.eta_R = 0;
y0_7 = [P0.S0; P0.I0; P0.Ag0; 0; P0.V0; 0; 0];
opts7 = odeset('RelTol',1e-10,'AbsTol',1e-12,'MaxStep',0.1,'NonNegative',1:7);
[t7, Y7] = ode15s(@(t,y) measles_ode_rhs(t,y,P0), [0 P0.tmax_days], y0_7, opts7);

% Compare S, I, V at shared time points
tc = linspace(0, P.tmax_days, 500)';
Y4i = interp1(t4, Y4, tc, 'pchip');
Y7i = interp1(t7, Y7, tc, 'pchip');

err_S = max(abs(Y4i(:,1) - Y7i(:,1)));
err_I = max(abs(Y4i(:,2) - Y7i(:,2)));
% In 7-var, A = Ag + A17 = Y7(:,3) + Y7(:,4)
err_A = max(abs(Y4i(:,3) - (Y7i(:,3)+Y7i(:,4))));
err_V = max(abs(Y4i(:,4) - Y7i(:,5)));

fprintf('    max|S4-S7| = %.2e\n', err_S);
fprintf('    max|I4-I7| = %.2e\n', err_I);
fprintf('    max|A4-(Ag+A17)| = %.2e\n', err_A);
fprintf('    max|V4-V7| = %.2e\n', err_V);
if max([err_S err_I err_A err_V]) < 1e-2
    fprintf('    [PASS] 7-var collapses to Morris core.\n');
else
    fprintf('    [WARN] Discrepancy detected. Check extension params.\n');
end

%% ================================================================
%  Summary table
% ================================================================
fprintf('\n  Summary Table (7-var model):\n');
fprintf('  %-20s %8s %8s %8s %8s %8s %8s\n', ...
    'Archetype','peakV','t_peak','peakAg','peakA17','peakR','peakAb');
fprintf('  %s\n', repmat('-',1,72));
for w = 1:6
    s = sols{w};
    [pV,iV] = max(s.V);
    fprintf('  %-20s %8.4f %8.1f %8.1f %8.1f %8.2f %8.3f\n', ...
        names{w}, pV, s.t(iV), max(s.Ag), max(s.A17), max(s.R), max(s.Ab));
end

fprintf('\n  Figures -> ./export/\n');
end
