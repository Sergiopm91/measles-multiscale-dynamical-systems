function gR = antibody_stimulation_function(R, P)
% ANTIBODY_STIMULATION_FUNCTION  RNA-driven antibody stimulation term.
%
%   gR = antibody_stimulation_function(R, P)
%
% Supported forms selected by P.Ab_stim_type:
%   'log'              : log(1+R)
%   'michaelis_menten' : R/(K_R+R)
%   'hill'             : R^n/(K_R^n+R^n)
%
% The default remains 'log' to preserve backward compatibility. The other
% forms are used for the v9 functional-sensitivity analysis requested by
% the technical review.

R = max(R, 0);
if ~isfield(P, 'Ab_stim_type') || isempty(P.Ab_stim_type)
    stim_type = 'log';
else
    stim_type = lower(string(P.Ab_stim_type));
end

if ~isfield(P, 'K_R') || isempty(P.K_R), P.K_R = 100; end
if ~isfield(P, 'n_R') || isempty(P.n_R), P.n_R = 2; end

switch stim_type
    case {"log", "log1p", "logarithmic"}
        gR = log1p(R);
    case {"michaelis_menten", "michaelis", "mm"}
        gR = R ./ (P.K_R + R + eps);
    case "hill"
        n = P.n_R;
        gR = R.^n ./ (P.K_R.^n + R.^n + eps);
    otherwise
        error('Unknown Ab_stim_type: %s', stim_type);
end

% Guard against numerical artifacts.
gR = max(gR, 0);
end
