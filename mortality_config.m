function cfg = mortality_config(map_name)
% MORTALITY_CONFIG  Named mortality-map calibrations for v9.
%
%   cfg = mortality_config('benchmark_aligned')
%   cfg = mortality_config('conservative')
%   cfg = mortality_config('high_risk_IC')
%
% v9 correction:
% The benchmark-aligned and high-risk maps include an interaction penalty
% between IFN-gamma effector deficit and cumulative infectivity burden:
%
%   a5 * max(0,1 - Ag_peak/Ag_ref) * max(0,AUC_beta/AUC_ref - 1)
%
% This term is biologically interpretable: mortality risk should increase
% sharply when prolonged/cumulative infectivity coexists with failure of the
% antiviral IFN-gamma-dominant effector branch. It resolves the previous
% problem in which lowering z0 globally increased mortality for healthy and
% vaccinated profiles while still underestimating the immunocompromised host.

if nargin < 1 || isempty(map_name), map_name = 'benchmark_aligned'; end
map_name = char(lower(string(map_name)));

cfg = struct();
cfg.map_name = map_name;
cfg.V_ref   = 57.42;
cfg.T_ref   = 15.0;
cfg.Ag_ref  = 119.8;
cfg.AUC_ref = 13.0;
cfg.M_max   = 0.45;
cfg.M_floor = 1e-5;
cfg.eps_val = 1e-9;

switch map_name
    case 'conservative'
        % v8-compatible behavior. Kept only as sensitivity map.
        cfg.a1 = 1.8; cfg.a2 = 2.0; cfg.a3 = 2.0; cfg.a4 = 1.2;
        cfg.a5 = 0.0;
        cfg.z0 = 4.8;
        cfg.description = 'Conservative v8-compatible mortality map without interaction term.';

    case 'benchmark_aligned'
        % Main v9 map. Keeps the original four interpretable terms and adds
        % an immune-deficit x infectivity-burden interaction. The threshold is
        % not globally lowered; therefore healthy and vaccinated profiles remain
        % below the therapy-trigger threshold, while the IC profile is lifted.
        cfg.a1 = 1.8; cfg.a2 = 2.0; cfg.a3 = 2.0; cfg.a4 = 1.2;
        cfg.a5 = 2.0;
        cfg.z0 = 5.2;
        cfg.description = ['Benchmark-aligned v9 mortality map with ', ...
            'immune-deficit/infectivity-burden interaction.'];

    case 'high_risk_ic'
        % Sensitivity stress-test close to the high-risk IC benchmark.
        cfg.a1 = 1.8; cfg.a2 = 2.0; cfg.a3 = 2.0; cfg.a4 = 1.2;
        cfg.a5 = 3.0;
        cfg.z0 = 5.4;
        cfg.description = ['High-risk immunocompromised sensitivity map with ', ...
            'stronger immune-deficit/infectivity-burden interaction.'];

    otherwise
        error('Unknown mortality map: %s', map_name);
end
end
