function T = phase_ode_numerical_stability()
% PHASE_ODE_NUMERICAL_STABILITY  Numerical robustness check for ODE settings.
%
% Varies RelTol, AbsTol, and MaxStep one at a time around the manuscript
% baseline and exports the stability of Peak V, t_peak, T_clear, AUC_beta,
% and key immune variables for the healthy-adult extended model.

fprintf('\n=== ODE Numerical Stability Check ===\n');
if ~exist('export','dir'), mkdir('export'); end

P = measles_params(1);
y0 = [P.S0; P.I0; P.Ag0; P.A17_0; P.V0; P.R0; P.Ab0];
base = struct('RelTol',1e-10,'AbsTol',1e-12,'MaxStep',0.1,'Label','baseline');
settings = { ...
    base, ...
    struct('RelTol',1e-8,  'AbsTol',1e-12,'MaxStep',0.1,'Label','RelTol_1e-8'), ...
    struct('RelTol',1e-11, 'AbsTol',1e-12,'MaxStep',0.1,'Label','RelTol_1e-11'), ...
    struct('RelTol',1e-10, 'AbsTol',1e-10,'MaxStep',0.1,'Label','AbsTol_1e-10'), ...
    struct('RelTol',1e-10, 'AbsTol',1e-13,'MaxStep',0.1,'Label','AbsTol_1e-13'), ...
    struct('RelTol',1e-10, 'AbsTol',1e-12,'MaxStep',0.05,'Label','MaxStep_0p05'), ...
    struct('RelTol',1e-10, 'AbsTol',1e-12,'MaxStep',0.2,'Label','MaxStep_0p2')};

rows = {};
base_metrics = [];
for i = 1:numel(settings)
    s = settings{i};
    opts = odeset('RelTol',s.RelTol,'AbsTol',s.AbsTol,'MaxStep',s.MaxStep,'NonNegative',1:7);
    tic;
    [t,Y] = ode15s(@(t,y) measles_ode_rhs(t,y,P), [0 P.tmax_days], y0, opts);
    elapsed = toc;
    M = ode_common_metrics(t,Y,'extended7',P);
    if i == 1, base_metrics = M; end

    rows(end+1,:) = {s.Label, s.RelTol, s.AbsTol, s.MaxStep, M.n_steps, elapsed, ...
        M.Peak_V, M.t_peak, M.T_clear, M.AUC_beta, M.Peak_Ag, M.Peak_A17, ...
        M.Peak_R, M.Peak_Ab, ...
        rel_diff(M.Peak_V, base_metrics.Peak_V), ...
        abs(M.t_peak - base_metrics.t_peak), ...
        abs(M.T_clear - base_metrics.T_clear), ...
        rel_diff(M.AUC_beta, base_metrics.AUC_beta), ...
        rel_diff(M.Peak_Ag, base_metrics.Peak_Ag), ...
        rel_diff(M.Peak_A17, base_metrics.Peak_A17)}; %#ok<AGROW>
end

T = cell2table(rows, 'VariableNames', {'Setting','RelTol','AbsTol','MaxStep','Steps','Runtime_s', ...
    'Peak_V','t_peak','T_clear','AUC_beta','Peak_Ag','Peak_A17','Peak_R','Peak_Ab', ...
    'relDiff_Peak_V','absDiff_t_peak','absDiff_T_clear','relDiff_AUC_beta', ...
    'relDiff_Peak_Ag','relDiff_Peak_A17'});

writetable(T, fullfile('export','phase_ode_numerical_stability.csv'));

fprintf('  Saved: export/phase_ode_numerical_stability.csv\n');
fprintf('  Max rel diff Peak V:   %.3e\n', max(T.relDiff_Peak_V));
fprintf('  Max abs diff t_peak:   %.3e days\n', max(T.absDiff_t_peak));
fprintf('  Max abs diff T_clear:  %.3e days\n', max(T.absDiff_T_clear));
fprintf('  Max rel diff AUC_beta: %.3e\n', max(T.relDiff_AUC_beta));

pass = max(T.relDiff_Peak_V) < 1e-3 && max(T.absDiff_t_peak) < 0.25 && ...
       max(T.absDiff_T_clear) < 0.5 && max(T.relDiff_AUC_beta) < 1e-3;
if pass
    fprintf('  [PASS] ODE outputs are stable under tolerance/MaxStep perturbations.\n');
else
    fprintf('  [WARN] Stability differences exceed conservative thresholds; inspect CSV.\n');
end
end

function d = rel_diff(x, ref)
d = abs(x-ref) / max(abs(ref), eps);
end
