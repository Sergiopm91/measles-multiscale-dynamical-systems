function M = ode_common_metrics(t, Y, mode, P)
% ODE_COMMON_METRICS  Shared metric extraction for 4-var and 7-var measles ODEs.
%
%   M = ode_common_metrics(t,Y,'core4',P)
%   M = ode_common_metrics(t,Y,'extended7',P)

if nargin < 4, P = struct(); end
Vth_clear = 1e-6;
Kv = 0.1;

M = struct();
M.n_steps = numel(t);
M.t_end = t(end);

if strcmpi(mode, 'core4')
    S = Y(:,1); I = Y(:,2); A = Y(:,3); V = Y(:,4);
    [M.Peak_V, iPk] = max(V);
    M.t_peak = t(iPk);
    M.min_S = min(S);
    M.peak_A = max(A);
    M.peak_I = max(I);
    idx_clear = find(V(iPk:end) < Vth_clear, 1, 'first');
    if isempty(idx_clear), M.T_clear = t(end); else, M.T_clear = t(iPk + idx_clear - 1); end
    beta_eff = V ./ (V + Kv);
    M.AUC_beta = trapz(t, beta_eff);
elseif strcmpi(mode, 'extended7')
    S=Y(:,1); I=Y(:,2); Ag=Y(:,3); A17=Y(:,4); V=Y(:,5); R=Y(:,6); Ab=Y(:,7);
    Atot = Ag + A17;
    [M.Peak_V, iPk] = max(V);
    M.t_peak = t(iPk);
    M.min_S = min(S);
    M.peak_I = max(I);
    M.Peak_Ag = max(Ag);
    M.Peak_A17 = max(A17);
    M.Peak_R = max(R);
    M.Peak_Ab = max(Ab);
    idx_clear = find(V(iPk:end) < Vth_clear, 1, 'first');
    if isempty(idx_clear), M.T_clear = t(end); else, M.T_clear = t(iPk + idx_clear - 1); end
    beta_eff = V ./ (V + Kv);
    M.AUC_beta = trapz(t, beta_eff);
    M.AUC_Ag = trapz(t, Ag);
    M.AUC_A17 = trapz(t, A17);
    M.A17_Ag_AUC_ratio = M.AUC_A17 / max(M.AUC_Ag, eps);
    M.A17_Ag_peak_ratio = M.Peak_A17 / max(M.Peak_Ag, eps);

    sig17 = compute_sigma17(t, R, P);
    M.t50_sigma17 = crossing_time(t, sig17, 0.50);
    M.t90_sigma17 = crossing_time(t, sig17, 0.90);
    M.final_sigma17 = sig17(end);
    M.rash_peak = max(Ag ./ (Ag + P.hR + eps));
    M.clearance_rate_peak = max(P.c + P.theta_Ab*Ab + P.k*Atot);
else
    error('Unknown mode: %s', mode);
end
end

function sig17 = compute_sigma17(t, R, P)
z = P.alpha_t.*(t - P.t17) + P.alpha_R.*log(1+R);
z = min(max(z, -60), 60);
sig17 = 1 ./ (1 + exp(-z));
end

function tc = crossing_time(t, y, threshold)
idx = find(y >= threshold, 1, 'first');
if isempty(idx)
    tc = NaN;
elseif idx == 1
    tc = t(1);
else
    t0 = t(idx-1); t1 = t(idx);
    y0 = y(idx-1); y1 = y(idx);
    if abs(y1-y0) < eps
        tc = t1;
    else
        tc = t0 + (threshold-y0)*(t1-t0)/(y1-y0);
    end
end
end
