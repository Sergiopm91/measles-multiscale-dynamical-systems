function phase5_sensitivity_analysis(mode)
% PHASE5_SENSITIVITY_ANALYSIS
% Global sensitivity analysis for the measles multiscale manuscript.
%
% This version strengthens the defensibility of Phase 5 by:
%   1) treating mortality as an ODE-informed severity-to-mortality mapping;
%   2) exporting mortality score contributions per archetype;
%   3) comparing mapped mortality against benchmark CFR anchors;
%   4) adding pseudo leave-one-out calibration diagnostics;
%   5) adding Sobol metadata, raw indices, clipped indices, and bootstrap CIs;
%   6) making Tier 2 explicitly a coupling-map surrogate analysis;
%   7) avoiding claims that AR_surrogate represents the full CA;
%   8) documenting non-identifiability or near-zero-variance cases.
%
% Usage:
%   phase5_sensitivity_analysis
%   phase5_sensitivity_analysis('fast')
%   phase5_sensitivity_analysis('paper')
%
% Outputs:
%   export/phase5_sensitivity_data.mat
%   export/phase5_metadata.csv
%   export/phase5_tier1_sobol.csv
%   export/phase5_tier2_sobol.csv
%   export/phase5_mortality_validation.csv
%   export/phase5_mortality_contributions.csv
%   export/phase5_mortality_pseudo_loo.csv
%   export/phase5_surrogate_vs_phase4.csv
%   export/fig16_tier1_sobol_indices.png
%   export/fig17_mortality_validation_linear.png
%   export/fig18_mortality_validation_log.png
%   export/fig19_mortality_contributions.png
%   export/fig20_tier2_surrogate_sobol.png

if nargin < 1 || isempty(mode)
    mode = 'paper';
end
mode = char(mode);


fprintf('================================================================\n');
fprintf('  Phase 5: Global Sensitivity and Mortality-Mapping Diagnostics\n');
fprintf('================================================================\n\n');

if ~exist('export','dir'), mkdir('export'); end

switch lower(mode)
    case 'fast'
        N_tier1 = 256;
        N_tier2 = 128;
        B_boot  = 80;
    case 'paper'
        N_tier1 = 4096;
        N_tier2 = 512;
        B_boot  = 300;
    case 'custom'
        envN = str2double(getenv('MEASLES_TIER1_N_BASE'));
        if isnan(envN) || envN <= 0, envN = 1024; end
        N_tier1 = envN;
        N_tier2 = 256;
        B_boot  = 120;
    otherwise
        error('Unknown mode "%s". Use fast, paper, or custom.', mode);
end

seed = 42;
rng(seed,'twister');

%% ------------------------------------------------------------------------
% Tier 1 parameter space
% -------------------------------------------------------------------------
param_display = {'\beta','k','p','q','c','\delta', ...
    '\rho_R','\mu_R','\rho_{Ab}','\theta_{Ab}', ...
    '\alpha_t','\alpha_R','\kappa_{17}','t_{17}'};

param_fields = {'beta','k','p','q','c','delta', ...
    'rho_R','mu_R','rho_Ab','theta_Ab', ...
    'alpha_t','alpha_R','kappa17','t17'};

n_params = numel(param_fields);
P0 = measles_params(1);

base_vals = zeros(n_params,1);
for i = 1:n_params
    base_vals(i) = P0.(param_fields{i});
end

ranges = zeros(n_params,2);
for i = 1:6
    ranges(i,:) = base_vals(i) * [0.70, 1.30];
end
for i = 7:14
    ranges(i,:) = base_vals(i) * [0.50, 1.50];
end
ranges(:,1) = max(ranges(:,1), 1e-9);
ranges(14,:) = [15, 55];

fprintf('  Tier 1 parameters: %d\n', n_params);
fprintf('  Tier 1 N_base=%d, total evaluations=%d\n', ...
    N_tier1, N_tier1*(2+n_params));
fprintf('  Random seed=%d\n\n', seed);

for i = 1:n_params
    fprintf('    %-12s: [%.6g, %.6g] base=%.6g\n', ...
        param_fields{i}, ranges(i,1), ranges(i,2), base_vals(i));
end

%% ------------------------------------------------------------------------
% Tier 1 Saltelli sampling
% -------------------------------------------------------------------------
[A01, B01, sampling_method] = make_saltelli_unit_samples(N_tier1, n_params, 1000, 5000);
A_phys = scale_unit_to_ranges(A01, ranges);
B_phys = scale_unit_to_ranges(B01, ranges);

qoi_names = {'peak_V','t_clear','peak_Ag','AUC_beta','M_emergent'};
qoi_display = {'Peak V','t_{clear}','Peak A_\gamma','AUC(\beta_{eff})','M_{emergent}'};
n_qoi = numel(qoi_names);

fprintf('\n  Evaluating Tier 1 matrix A... '); tic;
Y_A = eval_samples(A_phys, param_fields, n_qoi);
fprintf('%.1f s\n', toc);

fprintf('  Evaluating Tier 1 matrix B... '); tic;
Y_B = eval_samples(B_phys, param_fields, n_qoi);
fprintf('%.1f s\n', toc);

Y_AB = zeros(N_tier1, n_qoi, n_params);
for i = 1:n_params
    fprintf('  Evaluating Tier 1 AB_%02d (%s)... ', i, param_fields{i}); tic;
    ABi = A_phys;
    ABi(:,i) = B_phys(:,i);
    Y_AB(:,:,i) = eval_samples(ABi, param_fields, n_qoi);
    fprintf('%.1f s\n', toc);
end

fprintf('\n  Computing Tier 1 Sobol indices and bootstrap CIs...\n');
tier1 = compute_sobol_with_bootstrap(Y_A, Y_B, Y_AB, B_boot, seed+100);

%% ------------------------------------------------------------------------
% Mortality validation and contribution diagnostics
% -------------------------------------------------------------------------
fprintf('\n  Computing mortality-mapping diagnostics...\n');
[mort_table, contrib_table, loo_table, profiles] = mortality_diagnostics();

%% ------------------------------------------------------------------------
% Tier 2 coupling-map surrogate sensitivity
% -------------------------------------------------------------------------
fprintf('\n================================================================\n');
fprintf('  TIER 2: Coupling-Map Surrogate Sensitivity\n');
fprintf('================================================================\n');
fprintf(['  This tier varies the mortality-map coefficients and a transmission\n', ...
         '  scaling parameter in a deterministic surrogate. It is not the full CA.\n']);

cp_names_display = {'a_1','a_2','a_3','a_4','z_0','\lambda'};
cp_fields = {'a1','a2','a3','a4','z0','lambda'};
cp_base = [1.8, 2.0, 2.0, 1.2, 4.8, 0.45];
cp_ranges = cp_base(:) .* [0.70, 1.30];
n_cp = numel(cp_base);

fprintf('  Tier 2 N_base=%d, total evaluations=%d\n', ...
    N_tier2, N_tier2*(2+n_cp));

[Acp01, Bcp01, sampling_method_t2] = make_saltelli_unit_samples(N_tier2, n_cp, 2000, 8000);
Acp = scale_unit_to_ranges(Acp01, cp_ranges);
Bcp = scale_unit_to_ranges(Bcp01, cp_ranges);

cp_qoi_names = {'M_bar_pop_weighted','M_IC','AR_surrogate_60pct'};
cp_qoi_display = {'Population-weighted M','M_{IC}','AR surrogate at 60% coverage'};
n_cp_qoi = numel(cp_qoi_names);

Y_Acp = eval_coupling_surrogate(Acp, profiles);
Y_Bcp = eval_coupling_surrogate(Bcp, profiles);
Y_ABcp = zeros(N_tier2, n_cp_qoi, n_cp);

for i = 1:n_cp
    ABi = Acp;
    ABi(:,i) = Bcp(:,i);
    Y_ABcp(:,:,i) = eval_coupling_surrogate(ABi, profiles);
end

tier2 = compute_sobol_with_bootstrap(Y_Acp, Y_Bcp, Y_ABcp, B_boot, seed+200);

%% ------------------------------------------------------------------------
% Surrogate readback against Phase 4 exported CA outputs, if available
% -------------------------------------------------------------------------
surrogate_phase4_table = compare_surrogate_with_phase4(profiles);

%% ------------------------------------------------------------------------
% Export tables
% -------------------------------------------------------------------------
fprintf('\n  Exporting Phase 5 tables and figures...\n');

metadata = build_metadata_table(mode, seed, sampling_method, sampling_method_t2, ...
    N_tier1, N_tier2, B_boot, n_params, n_cp, param_fields, ranges, cp_fields, cp_ranges);
writetable(metadata, fullfile('export','phase5_metadata.csv'));

tier1_table = sobol_to_table(tier1, param_fields, param_display, qoi_names, 'Tier1_within_host');
writetable(tier1_table, fullfile('export','phase5_tier1_sobol.csv'));

tier2_table = sobol_to_table(tier2, cp_fields, cp_names_display, cp_qoi_names, 'Tier2_coupling_surrogate');
writetable(tier2_table, fullfile('export','phase5_tier2_sobol.csv'));

writetable(mort_table, fullfile('export','phase5_mortality_validation.csv'));
writetable(contrib_table, fullfile('export','phase5_mortality_contributions.csv'));
writetable(loo_table, fullfile('export','phase5_mortality_pseudo_loo.csv'));
writetable(surrogate_phase4_table, fullfile('export','phase5_surrogate_vs_phase4.csv'));

%% ------------------------------------------------------------------------
% Figures
% -------------------------------------------------------------------------
make_tier1_figure(tier1, param_display, qoi_display);
make_mortality_figures(mort_table, contrib_table);
make_tier2_figure(tier2, cp_names_display, cp_qoi_display);

%% ------------------------------------------------------------------------
% MAT export
% -------------------------------------------------------------------------
sa = struct();
sa.mode = mode;
sa.seed = seed;
sa.N_tier1 = N_tier1;
sa.N_tier2 = N_tier2;
sa.B_boot = B_boot;
sa.param_fields = param_fields;
sa.param_display = param_display;
sa.param_ranges = ranges;
sa.qoi_names = qoi_names;
sa.tier1 = tier1;
sa.coupling_fields = cp_fields;
sa.coupling_display = cp_names_display;
sa.coupling_ranges = cp_ranges;
sa.coupling_qoi_names = cp_qoi_names;
sa.tier2 = tier2;
sa.mortality_validation = mort_table;
sa.mortality_contributions = contrib_table;
sa.mortality_pseudo_loo = loo_table;
sa.surrogate_vs_phase4 = surrogate_phase4_table;

save(fullfile('export','phase5_sensitivity_data.mat'),'-struct','sa');

%% ------------------------------------------------------------------------
% Console tables
% -------------------------------------------------------------------------
print_sobol_summary('TABLE E1: Tier 1 Sobol Indices', tier1, param_fields, qoi_names);
print_mortality_tables(mort_table, contrib_table, loo_table);
print_sobol_summary('TABLE G1: Tier 2 Coupling-Surrogate Sobol Indices', tier2, cp_fields, cp_qoi_names);

fprintf('\n================================================================\n');
fprintf('  Phase 5 complete.\n');
fprintf('  Saved MAT: export/phase5_sensitivity_data.mat\n');
fprintf('  Saved CSV: export/phase5_*.csv\n');
fprintf('  Saved PNG: export/fig16-fig20_*.png\n');
fprintf('================================================================\n');
end

%% =========================================================================
% Helper functions
% =========================================================================

function [A01, B01, method] = make_saltelli_unit_samples(N, d, skipA, skipB)
try
    ss = sobolset(d,'Skip',skipA,'Leap',37);
    ss = scramble(ss,'MatousekAffineOwen');
    A01 = net(ss,N);

    ss2 = sobolset(d,'Skip',skipB,'Leap',53);
    ss2 = scramble(ss2,'MatousekAffineOwen');
    B01 = net(ss2,N);

    method = 'sobolset_scrambled_MatousekAffineOwen';
catch
    warning('sobolset unavailable. Falling back to Latin-hypercube-like random stratification.');
    A01 = zeros(N,d);
    B01 = zeros(N,d);
    for j = 1:d
        pm = randperm(N);
        A01(:,j) = (pm(:) - rand(N,1)) / N;
        pm = randperm(N);
        B01(:,j) = (pm(:) - rand(N,1)) / N;
    end
    method = 'fallback_stratified_random';
end
end

function X = scale_unit_to_ranges(U, ranges)
N = size(U,1);
d = size(U,2);
X = zeros(N,d);
for i = 1:d
    X(:,i) = ranges(i,1) + U(:,i) .* (ranges(i,2)-ranges(i,1));
end
end

function Y = eval_samples(X, param_fields, n_qoi)
N = size(X,1);
Y = NaN(N,n_qoi);

for j = 1:N
    P = measles_params(1);
    for i = 1:numel(param_fields)
        P.(param_fields{i}) = X(j,i);
    end

    y0 = [P.S0; P.I0; P.Ag0; P.A17_0; P.V0; P.R0; P.Ab0];
    opts = odeset('RelTol',1e-8,'AbsTol',1e-10,'MaxStep',0.5,'NonNegative',1:7);

    try
        [t,Ys] = ode15s(@(t,y) measles_ode_rhs(t,y,P), [0 P.tmax_days], y0, opts);

        V  = Ys(:,5);
        Ag = Ys(:,3);

        peakV = max(V);
        peakAg = max(Ag);

        [~,iPk] = max(V);
        idx_clear = find(V(iPk:end) < 0.01, 1, 'first');
        if isempty(idx_clear)
            t_clear = P.tmax_days;
        else
            t_clear = t(iPk + idx_clear - 1);
        end

        beta_eff = V ./ (V + 0.1);
        auc_beta = trapz(t,beta_eff);

        % CA-compatible infectious duration.
        % Detectable viremia onset is an ODE descriptor, but it is not used as the
        % CA infectious window. This helper keeps Phase 5 aligned with the
        % rash-centered timing logic used in Phase 4.
        T_inf = estimate_ca_infectious_duration(t, V, Ag, P);
        
        M = emergent_mortality(peakV, T_inf, peakAg, auc_beta, mortality_config('benchmark_aligned'));

        Y(j,:) = [peakV, t_clear, peakAg, auc_beta, M];

    catch ME
        warning('Tier 1 sample %d failed: %s', j, ME.message);
        Y(j,:) = NaN(1,n_qoi);
    end
end
end

function T_inf = estimate_ca_infectious_duration(t, V, Ag, P)
% ESTIMATE_CA_INFECTIOUS_DURATION
% CA-compatible infectious duration used for Tier 1 sensitivity.
%
% The ODE viral-threshold crossing is retained as a within-host descriptor,
% but the population CA uses a rash-centered infectious window. For the
% healthy-adult-like Tier 1 perturbations, the default infectious duration
% is 8 days. Conservative prolongation is allowed only when ODE clearance is
% substantially delayed.
%
% This keeps Phase 5 consistent with the Phase 4 coupling philosophy:
%   - viremia onset is not the CA latent period;
%   - infectiousness is centered around rash timing;
%   - prolonged ODE clearance can extend the infectious window.

Vth = 0.01;

if isempty(t) || isempty(V) || all(~isfinite(V))
    T_inf = 8;
    return;
end

V = max(V(:), 0);
t = t(:);

[~, i_peak] = max(V);

idx_clear = find(V(i_peak:end) < Vth, 1, 'first');
if isempty(idx_clear)
    t_ode_clear = max(t);
else
    t_ode_clear = t(i_peak + idx_clear - 1);
end

% Default measles-compatible infectious window: approximately 4 days before
% rash to 4 days after rash.
T_inf = 8;

% Conservative prolongation for unusually delayed viral clearance.
% This is intentionally bounded to avoid reverting to the old
% threshold-to-threshold ODE duration.
if t_ode_clear > 21
    T_inf = 8 + round(0.35 * (t_ode_clear - 21));
end

T_inf = max(8, min(T_inf, 21));

% Severe immune dysfunction can be represented by longer windows only when
% the parameter structure clearly resembles the immunocompromised setting.
if isstruct(P) && isfield(P,'q') && isfield(P,'k')
    if P.q < 0.35 || P.k < 0.006
        T_inf = max(T_inf, 12);
        if t_ode_clear > 35
            T_inf = min(21, T_inf + round(0.20 * (t_ode_clear - 35)));
        end
    end
end
end

function out = compute_sobol_with_bootstrap(Y_A, Y_B, Y_AB, B_boot, seed)
rng(seed,'twister');

[N, n_qoi] = size(Y_A);
n_params = size(Y_AB,3);

S1_raw = NaN(n_params,n_qoi);
ST_raw = NaN(n_params,n_qoi);
S1 = NaN(n_params,n_qoi);
ST = NaN(n_params,n_qoi);
VarY = NaN(1,n_qoi);
n_valid = zeros(1,n_qoi);
low_variance = false(1,n_qoi);

for q = 1:n_qoi
    yA = Y_A(:,q);
    yB = Y_B(:,q);

    valid = isfinite(yA) & isfinite(yB);
    for i = 1:n_params
        valid = valid & isfinite(Y_AB(:,q,i));
    end

    n_valid(q) = sum(valid);
    if n_valid(q) < max(20, round(0.20*N))
        continue;
    end

    yA = yA(valid);
    yB = yB(valid);
    yAB = Y_AB(valid,q,:);

    VarY(q) = var([yA; yB], 1);
    if VarY(q) < 1e-14
        low_variance(q) = true;
        continue;
    end

    for i = 1:n_params
        yi = yAB(:,:,i);
        S1_raw(i,q) = mean(yB .* (yi - yA)) / VarY(q);
        ST_raw(i,q) = mean((yA - yi).^2) / (2*VarY(q));
    end
end

S1 = min(max(S1_raw,0),1);
ST = min(max(ST_raw,0),1);

S1_lo = NaN(n_params,n_qoi);
S1_hi = NaN(n_params,n_qoi);
ST_lo = NaN(n_params,n_qoi);
ST_hi = NaN(n_params,n_qoi);

for q = 1:n_qoi
    yA0 = Y_A(:,q);
    yB0 = Y_B(:,q);

    valid = isfinite(yA0) & isfinite(yB0);
    for i = 1:n_params
        valid = valid & isfinite(Y_AB(:,q,i));
    end

    idx_valid = find(valid);
    if numel(idx_valid) < max(20, round(0.20*N)) || low_variance(q)
        continue;
    end

    bootS1 = NaN(B_boot,n_params);
    bootST = NaN(B_boot,n_params);

    for b = 1:B_boot
        idx = idx_valid(randi(numel(idx_valid), numel(idx_valid), 1));
        yA = Y_A(idx,q);
        yB = Y_B(idx,q);
        VarB = var([yA; yB],1);
        if VarB < 1e-14
            continue;
        end
        for i = 1:n_params
            yi = Y_AB(idx,q,i);
            s1b = mean(yB .* (yi - yA)) / VarB;
            stb = mean((yA - yi).^2) / (2*VarB);
            bootS1(b,i) = min(max(s1b,0),1);
            bootST(b,i) = min(max(stb,0),1);
        end
    end

    for i = 1:n_params
        S1_lo(i,q) = prctile(bootS1(:,i),2.5);
        S1_hi(i,q) = prctile(bootS1(:,i),97.5);
        ST_lo(i,q) = prctile(bootST(:,i),2.5);
        ST_hi(i,q) = prctile(bootST(:,i),97.5);
    end
end

out = struct();
out.S1_raw = S1_raw;
out.ST_raw = ST_raw;
out.S1 = S1;
out.ST = ST;
out.S1_ci_low = S1_lo;
out.S1_ci_high = S1_hi;
out.ST_ci_low = ST_lo;
out.ST_ci_high = ST_hi;
out.VarY = VarY;
out.n_valid = n_valid;
out.low_variance = low_variance;
out.note = ['Raw indices are exported. Negative estimates can occur from ', ...
    'finite Monte Carlo error. Clipped indices are used only for ranking ', ...
    'and plotting, with raw values preserved for transparency.'];
end

function [mort_table, contrib_table, loo_table, profiles] = mortality_diagnostics()
% MORTALITY_DIAGNOSTICS
% Mortality diagnostics aligned with Phase 4.
%
% Priority:
%   1) If export/phase4_coupling_profiles.csv exists, use it as the
%      authoritative coupling source for T_inf and M_emergent.
%   2) If the Phase 4 file is unavailable, fall back to ODE recomputation
%      with CA-compatible infectious-duration estimation.
%
% This prevents Phase 5 from recomputing mortality using the old
% threshold-to-threshold ODE infectious duration.

arch_names = {'Healthy adult','Child','Elderly','Immunocomp.','Partial vacc','Full vacc'};

% Internal consistency benchmarks only.
% These are not independent clinical validation data.
benchmark_cfr = [0.0020; 0.0100; 0.0150; 0.3000; 0.0005; 0.0001];

n = numel(arch_names);

peakV = zeros(n,1);
T_inf = zeros(n,1);
peakAg = zeros(n,1);
auc_beta = zeros(n,1);
M = zeros(n,1);
Z = zeros(n,1);

viral_burden = zeros(n,1);
infectious_duration = zeros(n,1);
immune_deficit = zeros(n,1);
infectivity_burden = zeros(n,1);
immune_infectivity_interaction = zeros(n,1);

profiles = repmat(struct(), n, 1);

phase4_file = fullfile('export','phase4_coupling_profiles.csv');
used_phase4 = false;

if exist(phase4_file,'file')
    try
        T4 = readtable(phase4_file);

        base_arch = get_table_var(T4, {'Base_archetype','BaseArchetype','base_archetype'});
        is_tx     = get_table_var(T4, {'Is_therapy','IsTherapy','is_therapy'});

        name_col  = get_table_var(T4, {'Name','name'});
        tinf_col  = get_table_var(T4, {'T_infectious_duration','Tinf','T_inf','T_infectious'});
        m_col     = get_table_var(T4, {'M_emergent','M','M_emerg'});
        pv_col    = get_table_var(T4, {'Peak_V','PeakV','peakV'});
        pag_col   = get_table_var(T4, {'Peak_Ag','PeakAg','peakAg'});
        auc_col   = get_table_var(T4, {'AUC_beta','AUC_b','AUCbeta','auc_beta'});

        for w = 1:n
            idx = find(double(base_arch) == w & double(is_tx) == 0, 1, 'first');

            if isempty(idx)
                error('No baseline Phase 4 coupling profile found for archetype %d.', w);
            end

            peakV(w)    = double(pv_col(idx));
            T_inf(w)    = double(tinf_col(idx));
            peakAg(w)   = double(pag_col(idx));
            auc_beta(w) = double(auc_col(idx));
            M(w)        = double(m_col(idx));

            profiles(w).name = char(string(name_col(idx)));
            profiles(w).source = 'phase4_coupling_profiles.csv';
        end

        used_phase4 = true;

    catch ME
        warning('Could not use Phase 4 coupling profiles: %s', ME.message);
        used_phase4 = false;
    end
end

if ~used_phase4
    warning(['Phase 4 coupling profile file was not available or could not be read. ', ...
        'Falling back to ODE recomputation with CA-compatible infectious duration.']);

    for w = 1:n
        P = measles_params(w);
        y0 = [P.S0; P.I0; P.Ag0; P.A17_0; P.V0; P.R0; P.Ab0];
        opts = odeset('RelTol',1e-10,'AbsTol',1e-12,'MaxStep',0.1,'NonNegative',1:7);

        [t,Y] = ode15s(@(t,y) measles_ode_rhs(t,y,P), [0 P.tmax_days], y0, opts);

        V = Y(:,5);
        Ag = Y(:,3);
        beta_eff = V ./ (V + 0.1);

        peakV(w)    = max(V);
        peakAg(w)   = max(Ag);
        auc_beta(w) = trapz(t,beta_eff);
        T_inf(w)    = estimate_ca_infectious_duration(t, V, Ag, P);

        [M(w), ~] = emergent_mortality(peakV(w), T_inf(w), peakAg(w), auc_beta(w), mortality_config('benchmark_aligned'));

        profiles(w).name = arch_names{w};
        profiles(w).source = 'ODE fallback with CA-compatible duration';
    end
end

% Compute score contributions using the same formula used by the mortality map.
% M itself is taken from Phase 4 when available, because Phase 4 is the
% authoritative ODE-to-CA coupling source.
for w = 1:n
    [M_formula, details] = emergent_mortality(peakV(w), T_inf(w), peakAg(w), auc_beta(w), mortality_config('benchmark_aligned'));

    Z(w) = details.Z;
    viral_burden(w) = details.viral_burden;
    infectious_duration(w) = details.infectious_duration;
    immune_deficit(w) = details.immune_deficit;
    infectivity_burden(w) = details.infectivity_burden;
    if isfield(details,'immune_infectivity_interaction')
        immune_infectivity_interaction(w) = details.immune_infectivity_interaction;
    else
        immune_infectivity_interaction(w) = 0;
    end

    profiles(w).peakV = peakV(w);
    profiles(w).T_inf = T_inf(w);
    profiles(w).peakAg = peakAg(w);
    profiles(w).auc_beta = auc_beta(w);
    profiles(w).M = M(w);
    profiles(w).M_formula_check = M_formula;
    profiles(w).Z = Z(w);
    profiles(w).benchmark_cfr = benchmark_cfr(w);
end

abs_diff = M - benchmark_cfr;
abs_diff_percent_points = 100 * abs_diff;
rel_error = abs_diff ./ max(benchmark_cfr, 1e-9);
log10_error = log10(max(M,1e-9)) - log10(max(benchmark_cfr,1e-9));

Archetype = string(arch_names(:));

if used_phase4
    Mortality_source = repmat("Phase 4 ODE-to-CA coupling profiles", n, 1);
else
    Mortality_source = repmat("ODE fallback with CA-compatible infectious duration", n, 1);
end

Benchmark_role = repmat("internal consistency benchmark; not independent clinical validation", n, 1);

mort_table = table(Archetype, peakV, T_inf, peakAg, auc_beta, Z, M, benchmark_cfr, ...
    abs_diff, abs_diff_percent_points, rel_error, log10_error, Mortality_source, Benchmark_role, ...
    'VariableNames', {'Archetype','Peak_V','T_inf','Peak_Ag','AUC_beta','Z_score', ...
    'M_emergent','Benchmark_CFR','Absolute_difference','Difference_percent_points', ...
    'Relative_error','Log10_error','Mortality_source','Benchmark_role'});

contrib_table = table(Archetype, viral_burden, infectious_duration, immune_deficit, ...
    infectivity_burden, immune_infectivity_interaction, Z, M, ...
    'VariableNames', {'Archetype','Z_viral_burden','Z_infectious_duration', ...
    'Z_immune_deficit','Z_infectivity_burden','Z_deficit_AUC_interaction', ...
    'Z_total','M_emergent'});

loo_table = pseudo_leave_one_out(Archetype, peakV, T_inf, peakAg, auc_beta, benchmark_cfr);
end

function loo_table = pseudo_leave_one_out(Archetype, peakV, T_inf, peakAg, auc_beta, benchmark_cfr)
% Pseudo-LOO diagnostic:
% With fixed weights and Mmax, fit only z0 using all benchmarks except one,
% then predict the held-out benchmark. This tests intercept robustness only.
% It is not a clinical validation procedure.

n = numel(benchmark_cfr);
cfg = local_mortality_cfg();

Z = zeros(n,1);
for i = 1:n
    dgamma = max(0, 1 - peakAg(i)/cfg.Ag_ref);
    bexcess = max(0, max(auc_beta(i),cfg.eps_val)/cfg.AUC_ref - 1);
    z1 = cfg.a1 * log(max(peakV(i),cfg.eps_val)/cfg.V_ref);
    z2 = cfg.a2 * (max(T_inf(i),cfg.eps_val)/cfg.T_ref - 1);
    z3 = cfg.a3 * dgamma;
    z4 = cfg.a4 * (max(auc_beta(i),cfg.eps_val)/cfg.AUC_ref - 1);
    z5 = cfg.a5 * dgamma * bexcess;
    Z(i) = z1 + z2 + z3 + z4 + z5;
end

z0_loo = zeros(n,1);
M_pred = zeros(n,1);
abs_diff = zeros(n,1);
rel_error = zeros(n,1);

for hold = 1:n
    train = true(n,1);
    train(hold) = false;

    y = min(max(benchmark_cfr(train), 1e-6), cfg.M_max-1e-6);
    logit_scaled = log(y ./ (cfg.M_max - y));

    % Since logit(M/Mmax) = Z - z0, z0 = Z - logit_scaled.
    z0_candidates = Z(train) - logit_scaled;
    z0_loo(hold) = median(z0_candidates);

    x = Z(hold) - z0_loo(hold);
    if x >= 0
        M_pred(hold) = cfg.M_max / (1 + exp(-x));
    else
        ex = exp(x);
        M_pred(hold) = cfg.M_max * ex / (1 + ex);
    end
    M_pred(hold) = max(M_pred(hold), cfg.M_floor);

    abs_diff(hold) = M_pred(hold) - benchmark_cfr(hold);
    rel_error(hold) = abs_diff(hold) / max(benchmark_cfr(hold), 1e-9);
end

Diagnostic_scope = repmat("pseudo leave-one-out; refits z0 only; not external validation", n, 1);

loo_table = table(Archetype, Z, benchmark_cfr, z0_loo, M_pred, abs_diff, rel_error, Diagnostic_scope, ...
    'VariableNames', {'Held_out_archetype','Z_score','Benchmark_CFR','z0_fit_without_held_out', ...
    'Predicted_M','Absolute_difference','Relative_error','Diagnostic_scope'});
end

function cfg = local_mortality_cfg()
cfg.V_ref   = 57.42;
cfg.T_ref   = 15.0;
cfg.Ag_ref  = 119.8;
cfg.AUC_ref = 13.0;
cfg.a1 = 1.8;
cfg.a2 = 2.0;
cfg.a3 = 2.0;
cfg.a4 = 1.2;
cfg.a5 = 2.0;
cfg.M_max = 0.45;
cfg.z0 = 5.2;
cfg.M_floor = 1e-5;
cfg.eps_val = 1e-9;
end

function Y = eval_coupling_surrogate(X, profiles)
% Deterministic Tier-2 surrogate.
%
% This is deliberately NOT the full cellular automaton. It measures how
% mortality-map coefficients and a transmission scaling parameter affect:
%   1) population-weighted mortality among baseline non-vaccinated hosts,
%   2) immunocompromised mortality,
%   3) a non-saturated final-size-inspired attack-rate surrogate at 60%
%      baseline vaccination.
%
% The surrogate is used only to verify coupling-map robustness.

N = size(X,1);
Y = NaN(N,3);

cfg = local_mortality_cfg();

% Demographic fractions among baseline non-vaccinated archetypes.
demo = [0.58, 0.22, 0.15, 0.05];

% Susceptibility-weighted fraction at 60% coverage:
% non-vaccinated 40%; vaccinated are 15% partial and 85% full.
S_eff_60 = 0.40 + 0.60*(0.15*0.35 + 0.85*0.08);

% A modest spatial-clustering damping factor calibrated only to keep the
% deterministic surrogate in the same order as Phase 4, not to replace it.
cluster_damping = 0.65;

beta_sum_ref = max(profiles(1).auc_beta, 1e-9);

for j = 1:N
    a1 = X(j,1);
    a2 = X(j,2);
    a3 = X(j,3);
    a4 = X(j,4);
    z0 = X(j,5);
    lambda = X(j,6);

    M_pop = 0;
    M_IC = NaN;

    for w = 1:4
        dgamma = max(0, 1 - profiles(w).peakAg/cfg.Ag_ref);
        bexcess = max(0, max(profiles(w).auc_beta,cfg.eps_val)/cfg.AUC_ref - 1);
        z = a1*log(max(profiles(w).peakV,cfg.eps_val)/cfg.V_ref) ...
          + a2*(max(profiles(w).T_inf,cfg.eps_val)/cfg.T_ref - 1) ...
          + a3*dgamma ...
          + a4*(max(profiles(w).auc_beta,cfg.eps_val)/cfg.AUC_ref - 1) ...
          + cfg.a5*dgamma*bexcess;

        Mw = cfg.M_max / (1 + exp(-(z - z0)));
        Mw = max(Mw, cfg.M_floor);

        M_pop = M_pop + demo(w)*Mw;
        if w == 4
            M_IC = Mw;
        end
    end

    % Non-saturated final-size-inspired surrogate.
    % R_contact is scaled down from local-contact CA pressure so that
    % AR_surrogate remains responsive over the lambda range.
    R_contact = 0.55 * lambda * beta_sum_ref;
    Reff = R_contact * S_eff_60;

    final_size_sus = solve_final_size(Reff);
    AR_surrogate = cluster_damping * S_eff_60 * final_size_sus;
    AR_surrogate = min(max(AR_surrogate,0),1);

    Y(j,:) = [M_pop, M_IC, AR_surrogate];
end
end

function z = solve_final_size(R)
% Solve z = 1 - exp(-R z), z in [0,1].
if R <= 1e-12
    z = 0;
    return;
end

z = max(1e-6, 1 - 1/max(R,1.000001));
for k = 1:100
    f = z - 1 + exp(-R*z);
    df = 1 - R*exp(-R*z);
    if abs(df) < 1e-10
        break;
    end
    zn = z - f/df;
    zn = min(max(zn,0),1);
    if abs(zn-z) < 1e-10
        z = zn;
        return;
    end
    z = zn;
end
end

function T = compare_surrogate_with_phase4(profiles)
coverage = [0; 60; 70; 80; 85; 90; 92; 95];
phase4_AR_mean = NaN(size(coverage));
phase4_AR_sd = NaN(size(coverage));

csvfile = fullfile('export','phase4_table_coverage_sweep.csv');
if exist(csvfile,'file')
    try
        T4 = readtable(csvfile);
        for i = 1:numel(coverage)
            idx = find(T4.Vax_percent == coverage(i),1,'first');
            if ~isempty(idx)
                phase4_AR_mean(i) = T4.AR_mean_percent(idx)/100;
                phase4_AR_sd(i) = T4.AR_sd_percent(idx)/100;
            end
        end
    catch ME
        warning('Could not read Phase 4 coverage table: %s', ME.message);
    end
end

cfg_lambda = 0.45;
beta_sum_ref = max(profiles(1).auc_beta,1e-9);
cluster_damping = 0.65;

surrogate_AR = zeros(size(coverage));
for i = 1:numel(coverage)
    v = coverage(i)/100;
    S_eff = (1-v) + v*(0.15*0.35 + 0.85*0.08);
    R_contact = 0.55 * cfg_lambda * beta_sum_ref;
    Reff = R_contact * S_eff;
    surrogate_AR(i) = min(max(cluster_damping * S_eff * solve_final_size(Reff),0),1);
end

absolute_difference = surrogate_AR - phase4_AR_mean;

Interpretation = repmat("surrogate readback only; full CA results remain authoritative", numel(coverage), 1);

T = table(coverage, surrogate_AR, phase4_AR_mean, phase4_AR_sd, absolute_difference, Interpretation, ...
    'VariableNames', {'Coverage_percent','AR_surrogate','Phase4_AR_mean','Phase4_AR_sd', ...
    'Surrogate_minus_Phase4','Interpretation'});
end

function metadata = build_metadata_table(mode, seed, sampling1, sampling2, N1, N2, B, np, ncp, pf, ranges, cpf, cpranges)
Key = strings(0,1);
Value = strings(0,1);

add('mode', mode);
add('seed', string(seed));
add('tier1_sampling', sampling1);
add('tier2_sampling', sampling2);
add('tier1_N_base', string(N1));
add('tier1_total_evaluations', string(N1*(2+np)));
add('tier2_N_base', string(N2));
add('tier2_total_evaluations', string(N2*(2+ncp)));
add('bootstrap_replicates', string(B));
add('sobol_negative_indices_policy', 'raw exported; clipped to [0,1] only for ranking/figures');
add('mortality_map_scope', 'ODE-informed severity-to-mortality mapping; not independent clinical CFR validation');
add('tier2_scope', 'coupling-map deterministic surrogate; not full cellular automaton');
add('tier3_CA_real_status', 'not executed here because Phase 4 CA engine is nested; surrogate is compared against exported Phase 4 outputs when available');

for i = 1:numel(pf)
    add("tier1_range_" + string(pf{i}), sprintf('[%.8g, %.8g]', ranges(i,1), ranges(i,2)));
end
for i = 1:numel(cpf)
    add("tier2_range_" + string(cpf{i}), sprintf('[%.8g, %.8g]', cpranges(i,1), cpranges(i,2)));
end

metadata = table(Key,Value);

    function add(k,v)
        Key(end+1,1) = string(k);
        Value(end+1,1) = string(v);
    end
end

function T = sobol_to_table(sobol, fields, display_names, qoi_names, tier_label)
rows = {};
for q = 1:numel(qoi_names)
    for i = 1:numel(fields)
        rows(end+1,:) = { ...
            string(tier_label), ...
            string(qoi_names{q}), ...
            string(fields{i}), ...
            string(display_names{i}), ...
            sobol.S1_raw(i,q), ...
            sobol.ST_raw(i,q), ...
            sobol.S1(i,q), ...
            sobol.ST(i,q), ...
            sobol.S1_ci_low(i,q), ...
            sobol.S1_ci_high(i,q), ...
            sobol.ST_ci_low(i,q), ...
            sobol.ST_ci_high(i,q), ...
            sobol.VarY(q), ...
            sobol.n_valid(q), ...
            sobol.low_variance(q)};
    end
end

T = cell2table(rows, 'VariableNames', {'Tier','QoI','Parameter','Parameter_display', ...
    'S1_raw','ST_raw','S1_clipped','ST_clipped','S1_ci95_low','S1_ci95_high', ...
    'ST_ci95_low','ST_ci95_high','QoI_variance','N_valid','Low_variance_flag'});
end

function make_tier1_figure(tier1, param_display, qoi_display)
fig = figure('Color','w','Position',[50 50 1250 650]);
n_qoi = numel(qoi_display);
n_params = numel(param_display);

for q = 1:n_qoi
    subplot(2,3,q); hold on;
    x = 1:n_params;
    bar(x-0.16, tier1.S1(:,q), 0.30, 'EdgeColor','none');
    bar(x+0.16, tier1.ST(:,q), 0.30, 'EdgeColor','none');
    set(gca,'XTick',x,'XTickLabel',param_display,'XTickLabelRotation',45,'FontSize',7);
    ylabel('Sobol index');
    title(qoi_display{q});
    ymax = max([tier1.ST(:,q); tier1.S1(:,q); 0.1]);
    ylim([0, min(1, ymax*1.25)]);
    grid on;
    if q == 1
        legend('S_i','S_{Ti}','Location','best','FontSize',7);
    end
end

subplot(2,3,6);
imagesc(tier1.ST');
colorbar;
set(gca,'XTick',1:n_params,'XTickLabel',param_display,'XTickLabelRotation',45, ...
    'YTick',1:n_qoi,'YTickLabel',qoi_display,'FontSize',7);
title('Total-effect heatmap');

sgtitle('Figure 16: Tier 1 within-host Sobol sensitivity','FontWeight','bold');
hide_tb(fig);
exportgraphics(fig, fullfile('export','fig16_tier1_sobol_indices.png'), 'Resolution', 300);
end

function make_mortality_figures(mort_table, contrib_table)
fig1 = figure('Color','w','Position',[100 80 720 460]);
hold on;
plot(mort_table.Benchmark_CFR*100, mort_table.M_emergent*100, 'o', ...
    'MarkerSize',9,'MarkerFaceColor',[0.35 0.35 0.35],'MarkerEdgeColor','k');
plot([1e-4 45],[1e-4 45],'k--','LineWidth',1);
xlabel('Benchmark CFR (%)');
ylabel('ODE-informed M (%)');
title('Figure 17: Mortality mapping versus benchmark CFR anchors');
grid on;
axis square;
xlim([0 35]);
ylim([0 35]);
text(mort_table.Benchmark_CFR*100, mort_table.M_emergent*100, mort_table.Archetype, ...
    'VerticalAlignment','bottom','HorizontalAlignment','left','FontSize',8);
hide_tb(fig1);
exportgraphics(fig1, fullfile('export','fig17_mortality_validation_linear.png'), 'Resolution', 300);

fig2 = figure('Color','w','Position',[120 100 720 460]);
loglog(mort_table.Benchmark_CFR, mort_table.M_emergent, 'o', ...
    'MarkerSize',9,'MarkerFaceColor',[0.35 0.35 0.35],'MarkerEdgeColor','k');
hold on;
loglog([1e-5 0.5],[1e-5 0.5],'k--','LineWidth',1);
xlabel('Benchmark CFR');
ylabel('ODE-informed M');
title('Figure 18: Mortality mapping on logarithmic scale');
grid on;
xlim([1e-5 0.5]);
ylim([1e-5 0.5]);
text(mort_table.Benchmark_CFR, mort_table.M_emergent, mort_table.Archetype, ...
    'VerticalAlignment','bottom','HorizontalAlignment','left','FontSize',8);
hide_tb(fig2);
exportgraphics(fig2, fullfile('export','fig18_mortality_validation_log.png'), 'Resolution', 300);

fig3 = figure('Color','w','Position',[80 80 1000 470]);
C = [contrib_table.Z_viral_burden, contrib_table.Z_infectious_duration, ...
     contrib_table.Z_immune_deficit, contrib_table.Z_infectivity_burden, ...
     contrib_table.Z_deficit_AUC_interaction];
bar(C,'stacked');
set(gca,'XTick',1:height(contrib_table),'XTickLabel',contrib_table.Archetype, ...
    'XTickLabelRotation',30,'FontSize',8);
ylabel('Contribution to Z');
title('Figure 19: ODE-derived severity-score contributions by archetype');
legend({'Viral burden','Infectious duration','Immune deficit', ...
    'AUC/infectivity burden','Deficit x AUC interaction'}, 'Location','best');
grid on;
hide_tb(fig3);
exportgraphics(fig3, fullfile('export','fig19_mortality_contributions.png'), 'Resolution', 300);
end

function make_tier2_figure(tier2, cp_names_display, cp_qoi_display)
fig = figure('Color','w','Position',[100 100 980 360]);
for q = 1:numel(cp_qoi_display)
    subplot(1,3,q); hold on;
    x = 1:numel(cp_names_display);
    bar(x-0.16, tier2.S1(:,q), 0.30, 'EdgeColor','none');
    bar(x+0.16, tier2.ST(:,q), 0.30, 'EdgeColor','none');
    set(gca,'XTick',x,'XTickLabel',cp_names_display,'XTickLabelRotation',45,'FontSize',8);
    ylabel('Sobol index');
    title(cp_qoi_display{q});
    ymax = max([tier2.ST(:,q); tier2.S1(:,q); 0.1]);
    ylim([0, min(1, ymax*1.25)]);
    grid on;
    if q == 1
        legend('S_i','S_{Ti}','Location','best','FontSize',8);
    end
end
sgtitle('Figure 20: Tier 2 coupling-map surrogate sensitivity','FontWeight','bold');
hide_tb(fig);
exportgraphics(fig, fullfile('export','fig20_tier2_surrogate_sobol.png'), 'Resolution', 300);
end

function print_sobol_summary(title_text, sobol, fields, qoi_names)
fprintf('\n================================================================\n');
fprintf('  %s\n', title_text);
fprintf('================================================================\n');

for q = 1:numel(qoi_names)
    fprintf('\n  %s:\n', qoi_names{q});
    fprintf('  %-14s %10s %10s %18s %18s\n', 'Parameter','S1','ST','S1 95% CI','ST 95% CI');
    fprintf('  %s\n', repmat('-',1,78));

    [~,ix] = sort(sobol.ST(:,q),'descend','MissingPlacement','last');
    for r = 1:min(6,numel(fields))
        i = ix(r);
        fprintf('  %-14s %10.4f %10.4f   [%6.3f,%6.3f]   [%6.3f,%6.3f]\n', ...
            fields{i}, sobol.S1(i,q), sobol.ST(i,q), ...
            sobol.S1_ci_low(i,q), sobol.S1_ci_high(i,q), ...
            sobol.ST_ci_low(i,q), sobol.ST_ci_high(i,q));
    end
end

fprintf('\n  Note: raw and clipped indices were exported. Negative raw estimates\n');
fprintf('  may occur because of finite Monte Carlo error and are not hidden.\n');
end

function print_mortality_tables(mort_table, contrib_table, loo_table)
fprintf('\n================================================================\n');
fprintf('  TABLE F1: ODE-Informed Mortality Mapping Diagnostic\n');
fprintf('================================================================\n');
fprintf('  %-16s %8s %8s %8s %8s %9s %9s %9s\n', ...
    'Archetype','peakV','T_inf','peakAg','AUC_b','M','Bench','AbsDiff');
fprintf('  %s\n', repmat('-',1,88));
for i = 1:height(mort_table)
    fprintf('  %-16s %8.2f %8.2f %8.2f %8.2f %9.4f %9.4f %9.4f\n', ...
        mort_table.Archetype(i), mort_table.Peak_V(i), mort_table.T_inf(i), ...
        mort_table.Peak_Ag(i), mort_table.AUC_beta(i), mort_table.M_emergent(i), ...
        mort_table.Benchmark_CFR(i), mort_table.Absolute_difference(i));
end

fprintf('\n================================================================\n');
fprintf('  TABLE F2: Severity-Score Contribution Breakdown\n');
fprintf('================================================================\n');
fprintf('  %-16s %10s %10s %10s %10s %10s %10s\n', ...
    'Archetype','Viral','Duration','Immune','AUC','Interact','Z_total');
fprintf('  %s\n', repmat('-',1,90));
for i = 1:height(contrib_table)
    fprintf('  %-16s %10.3f %10.3f %10.3f %10.3f %10.3f %10.3f\n', ...
        contrib_table.Archetype(i), contrib_table.Z_viral_burden(i), ...
        contrib_table.Z_infectious_duration(i), contrib_table.Z_immune_deficit(i), ...
        contrib_table.Z_infectivity_burden(i), ...
        contrib_table.Z_deficit_AUC_interaction(i), contrib_table.Z_total(i));
end

fprintf('\n================================================================\n');
fprintf('  TABLE F3: Pseudo Leave-One-Out Diagnostic\n');
fprintf('================================================================\n');
fprintf('  %-16s %9s %9s %9s %9s\n', ...
    'Held-out','Bench','Pred','z0_LOO','RelErr');
fprintf('  %s\n', repmat('-',1,62));
for i = 1:height(loo_table)
    fprintf('  %-16s %9.4f %9.4f %9.3f %9.2f\n', ...
        loo_table.Held_out_archetype(i), loo_table.Benchmark_CFR(i), ...
        loo_table.Predicted_M(i), loo_table.z0_fit_without_held_out(i), ...
        loo_table.Relative_error(i));
end

fprintf('\n  Benchmark comparisons and pseudo-LOO diagnostics are internal\n');
fprintf('  consistency checks, not independent clinical validation.\n');
end

function hide_tb(fig)
ax = findall(fig,'Type','axes');
for k = 1:numel(ax)
    try
        ax(k).Toolbar.Visible = 'off';
    catch
    end
end
end

function v = get_table_var(T, candidates)
% GET_TABLE_VAR
% Robustly retrieve a table variable using several possible names.

names = T.Properties.VariableNames;

for i = 1:numel(candidates)
    c = candidates{i};

    idx = find(strcmp(names, c), 1, 'first');
    if ~isempty(idx)
        v = T.(names{idx});
        return;
    end

    c_valid = matlab.lang.makeValidName(c);
    idx = find(strcmp(names, c_valid), 1, 'first');
    if ~isempty(idx)
        v = T.(names{idx});
        return;
    end
end

% More permissive fallback: ignore case and underscores.
norm_names = lower(regexprep(names, '_', ''));
for i = 1:numel(candidates)
    c_norm = lower(regexprep(candidates{i}, '_', ''));
    idx = find(strcmp(norm_names, c_norm), 1, 'first');
    if ~isempty(idx)
        v = T.(names{idx});
        return;
    end
end

error('None of the requested table variables were found: %s', strjoin(candidates, ', '));
end
