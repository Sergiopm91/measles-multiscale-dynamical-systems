function [M, details] = emergent_mortality(peakV, T_inf, peakAg, auc_beta, cfg)
% EMERGENT_MORTALITY
% ODE-informed severity-to-mortality mapping for the measles multiscale model.
%
%   M = emergent_mortality(peakV, T_inf, peakAg, auc_beta)
%   [M, details] = emergent_mortality(...)
%
% v9 correction:
% The main benchmark-aligned map keeps the original four-term score and adds
% a biologically interpretable high-risk interaction:
%
%   a5 * D_gamma * B_beta,
%
% where
%   D_gamma = max(0, 1 - peakAg/Ag_ref)
%   B_beta  = max(0, auc_beta/AUC_ref - 1).
%
% This term only activates when IFN-gamma effector deficit coexists with
% above-reference cumulative infectivity burden. It therefore lifts the
% immunocompromised high-risk profile without globally inflating healthy or
% vaccinated mortality.
%
% The mapping is semi-mechanistic and phenomenological. It is intended for
% transparent ODE-to-CA coupling, not as a clinically validated CFR model.

if nargin < 5 || isempty(cfg)
    cfg = mortality_config('benchmark_aligned');
elseif ischar(cfg) || isstring(cfg)
    cfg = mortality_config(cfg);
end
if ~isfield(cfg,'a5'), cfg.a5 = 0.0; end

% Numerical safety
peakV    = max(double(peakV),    cfg.eps_val);
T_inf    = max(double(T_inf),    cfg.eps_val);
peakAg   = max(double(peakAg),   0);
auc_beta = max(double(auc_beta), cfg.eps_val);

% Individual severity contributions
viral_burden = cfg.a1 * log(peakV / cfg.V_ref);
infectious_duration = cfg.a2 * (T_inf / cfg.T_ref - 1);
immune_deficit_fraction = max(0, 1 - peakAg / cfg.Ag_ref);
infectivity_excess_fraction = max(0, auc_beta / cfg.AUC_ref - 1);
immune_deficit = cfg.a3 * immune_deficit_fraction;
infectivity_burden = cfg.a4 * (auc_beta / cfg.AUC_ref - 1);
immune_infectivity_interaction = cfg.a5 * immune_deficit_fraction * infectivity_excess_fraction;

Z = viral_burden + infectious_duration + immune_deficit + ...
    infectivity_burden + immune_infectivity_interaction;

% Stable logistic transform
x = Z - cfg.z0;
if x >= 0
    M = cfg.M_max / (1 + exp(-x));
else
    ex = exp(x);
    M = cfg.M_max * ex / (1 + ex);
end

% Floor only for stochastic CA numerical stability. This is not a biological
% lower bound and should not be interpreted as a minimum clinical CFR.
M_raw = M;
M = max(M, cfg.M_floor);

if nargout > 1
    details = struct();
    details.label = 'ODE-informed severity-to-mortality mapping';
    if isfield(cfg,'map_name'), details.map_name = cfg.map_name; end
    if isfield(cfg,'description'), details.map_description = cfg.description; end
    details.interpretation = ['Semi-mechanistic phenomenological map; ', ...
        'not an independently validated clinical CFR model'];
    details.peakV = peakV;
    details.T_inf = T_inf;
    details.peakAg = peakAg;
    details.auc_beta = auc_beta;
    details.V_ref = cfg.V_ref;
    details.T_ref = cfg.T_ref;
    details.Ag_ref = cfg.Ag_ref;
    details.AUC_ref = cfg.AUC_ref;
    details.a1 = cfg.a1;
    details.a2 = cfg.a2;
    details.a3 = cfg.a3;
    details.a4 = cfg.a4;
    details.a5 = cfg.a5;
    details.z0 = cfg.z0;
    details.M_max = cfg.M_max;
    details.M_floor = cfg.M_floor;
    details.Z = Z;
    details.M_raw = M_raw;
    details.M = M;
    details.viral_burden = viral_burden;
    details.infectious_duration = infectious_duration;
    details.immune_deficit = immune_deficit;
    details.infectivity_burden = infectivity_burden;
    details.immune_infectivity_interaction = immune_infectivity_interaction;
    details.immune_deficit_fraction = immune_deficit_fraction;
    details.infectivity_excess_fraction = infectivity_excess_fraction;
    details.contribution_names = { ...
        'viral_burden', ...
        'infectious_duration', ...
        'immune_deficit', ...
        'infectivity_burden', ...
        'immune_infectivity_interaction'};
    details.contribution_values = [ ...
        viral_burden, infectious_duration, immune_deficit, ...
        infectivity_burden, immune_infectivity_interaction];
    details.calibration_note = ['Benchmarks are internal consistency anchors ', ...
        'for a bounded severity-to-risk mapping. They are not external ', ...
        'clinical validation data.'];
end
end
