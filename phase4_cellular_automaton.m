function phase4_cellular_automaton(mode)
% PHASE4_CELLULAR_AUTOMATON  Corrected multiscale ODE->CA measles phase.
%
%   phase4_cellular_automaton()          runs PAPER mode by default.
%   phase4_cellular_automaton('fast')    quick smoke test.
%   phase4_cellular_automaton('paper')   publication tables/figures.
%
% Main corrections relative to v5:
%   1) Viremia onset from the ODE is no longer used as the CA latent period.
%      It is exported as T_viremia_onset, an intra-host descriptor only.
%   2) The CA uses epidemiological measles timings:
%        T_infectious_onset = T_rash - 4 days,
%        T_infectious_duration = approximately 8 days for baseline hosts,
%      with conservative prolongation for immunocompromised profiles.
%   3) Coverage sweep, therapy strategies, and campaigns use explicit scenario
%      metadata, matched seeds, and identical baseline definitions.
%   4) All reported quantities from stochastic ensembles are printed as
%      mean +/- SD, with 95% normal-approximation CIs exported to CSV.
%   5) Population balance is checked every simulated day:
%      S + E + Ipre + Irash + R + D + Vac = N^2.
%
% Required files:
%   measles_params.m, measles_ode_rhs.m, emergent_mortality.m

if nargin < 1 || isempty(mode), mode = 'paper'; end
mode = lower(string(mode));
if ~ismember(mode, ["fast","paper"])
    error('Unknown mode "%s". Use ''fast'' or ''paper''.', mode);
end

close all; clc;
fprintf('================================================================\n');
fprintf('  Phase 4: Corrected ODE-to-CA population model (%s mode)\n', upper(mode));
fprintf('================================================================\n\n');

if ~exist('export','dir'), mkdir('export'); end
rng(42, 'twister');

%% ------------------------------------------------------------------------
%  1. ODE-derived profile descriptors and epidemiological CA timings
% -------------------------------------------------------------------------
fprintf('--- Section 1: Coupling map C: ODE -> CA profiles ---\n');
profiles = extract_coupled_profiles();
try
    estimate_ca_reproduction_number(char(mode));
catch ME
    warning('CA reproduction estimate failed: %s', ME.message);
end

%% ------------------------------------------------------------------------
%  2. CA configuration
% -------------------------------------------------------------------------
CA = default_ca_config(mode, profiles);
print_ca_config(CA);

vax_levels = [0.00, 0.60, 0.70, 0.80, 0.85, 0.90, 0.92, 0.95];
arch_names = {'Healthy adult','Child (<5y)','Elderly (>65y)', ...
              'Immunocomp.','Partial vacc','Full vacc'};

%% ------------------------------------------------------------------------
%  3. Vaccination coverage sweep
% -------------------------------------------------------------------------
fprintf('\n--- Section 3: Vaccination coverage sweep ---\n');
results = run_vax_sweep(CA, vax_levels, CA.nruns_main);

%% ------------------------------------------------------------------------
%  4. Therapy strategies at 60% baseline coverage
% -------------------------------------------------------------------------
fprintf('\n--- Section 4: Therapy comparison (60%% baseline coverage) ---\n');
[therapy_results, therapy_names] = run_therapy_comparison(CA, 0.60, CA.nruns_main);

%% ------------------------------------------------------------------------
%  5. Vaccination campaigns
% -------------------------------------------------------------------------
fprintf('\n--- Section 5: Vaccination campaigns ---\n');
[campaign_results, campaign_names] = run_campaigns(CA, CA.nruns_main);

%% ------------------------------------------------------------------------
%  6. Finite-size check
% -------------------------------------------------------------------------
fprintf('\n--- Section 6: Finite-size check ---\n');
CAfs = CA;
CAfs.N = CA.N_finite;
CAfs.nruns_main = CA.nruns_finite;
% Keep the same imported-case density in the finite-size comparison.
CAfs.n_seeds = max(1, round(CA.n_seeds * (CAfs.N/CA.N)^2));
vax_check = [0.00, 0.80, 0.92, 0.95];
results_fs = run_vax_sweep(CAfs, vax_check, CA.nruns_finite);
print_finite_size_table(vax_levels, results, vax_check, results_fs, CA, CAfs);

%% ------------------------------------------------------------------------
%  7. Lattice snapshots for one illustrative trajectory
% -------------------------------------------------------------------------
fprintf('\n--- Section 7: Lattice snapshots ---\n');
snap_days = [0, 10, 20, 40, 60, 100];
cspec = struct('scenario','Combined campaign snapshot', ...
               'base_vax',0.60, 'use_therapy',true, ...
               'events',{{{10,'targeted',0.80,2}, {15,'ring',0.90}}});
rng(CA.base_seed + 900001, 'twister');
snapshots = run_ca_snapshots(CA, cspec, snap_days);

%% ------------------------------------------------------------------------
%  8. Figures and tables
% -------------------------------------------------------------------------
fprintf('\n--- Generating figures and tables ---\n');
make_figures(CA, vax_levels, results, therapy_results, therapy_names, ...
             campaign_results, campaign_names, snapshots, snap_days);

Tcoverage  = build_coverage_table(vax_levels, results);
Ttherapy   = build_named_table(therapy_names, therapy_results, 'Strategy');
Tcampaign  = build_named_table(campaign_names, campaign_results, 'Campaign');
Tarch80    = build_archetype_table(arch_names, results(vax_levels==0.80).arch_stats, profiles);
Tmetadata  = build_metadata_table(CA, mode);
Tprofiles  = build_profiles_table(profiles);

print_profiles_table(Tprofiles);
print_coverage_table(Tcoverage);
print_named_outcomes(Ttherapy, 'TABLE B: Therapy Impact (60% coverage)');
print_named_outcomes(Tcampaign, 'TABLE C: Campaign Outcomes');
print_archetype_table(Tarch80);

writetable(Tcoverage, fullfile('export','phase4_table_coverage_sweep.csv'));
writetable(Ttherapy,  fullfile('export','phase4_table_therapy_strategies.csv'));
writetable(Tcampaign, fullfile('export','phase4_table_campaigns.csv'));
writetable(Tarch80,   fullfile('export','phase4_table_archetype_outcomes_80pct.csv'));
writetable(Tmetadata, fullfile('export','phase4_metadata.csv'));
writetable(Tprofiles, fullfile('export','phase4_coupling_profiles.csv'));

phase4 = struct();
phase4.mode = char(mode);
phase4.CA = CA;
phase4.profiles = profiles;
phase4.vax_levels = vax_levels;
phase4.coverage_results = results;
phase4.therapy_results = therapy_results;
phase4.therapy_names = therapy_names;
phase4.campaign_results = campaign_results;
phase4.campaign_names = campaign_names;
phase4.finite_size_results = results_fs;
phase4.tables.coverage = Tcoverage;
phase4.tables.therapy = Ttherapy;
phase4.tables.campaign = Tcampaign;
phase4.tables.archetype80 = Tarch80;
phase4.tables.metadata = Tmetadata;
phase4.tables.profiles = Tprofiles;
save(fullfile('export','phase4_data.mat'), 'phase4');

fprintf('\n--- Export complete ---\n');
fprintf('  Saved MAT: export/phase4_data.mat\n');
fprintf('  Saved CSV: export/phase4_table_*.csv, phase4_metadata.csv, phase4_coupling_profiles.csv\n');
fprintf('  Saved PNG: export/fig11-fig15_*.png\n');
fprintf('================================================================\n');
end

%% ========================================================================
%  Configuration
% ========================================================================
function CA = default_ca_config(mode, profiles)
CA = struct();
CA.mode = char(mode);
CA.N = 150;
CA.N_finite = 300;
CA.dt = 1;
CA.T_max = 365;
CA.n_seeds = 20;
CA.seed_initial_stage = 'staggered_infectious';
CA.base_seed = 42000;
CA.seed_stride = 100000;
CA.profiles = profiles;

% Local Moore-neighborhood force-of-infection multiplier. The value is kept
% fixed across all scenarios. Calibration is reported as metadata rather than
% re-estimated inside scenarios.
CA.lambda_inf = 0.45;

% Demographic composition among non-vaccinated host archetypes.
CA.host_demo = [0.58, 0.22, 0.15, 0.05];
CA.partial_frac = 0.15;
CA.full_frac = 0.85;
CA.cluster_sigma = 6.0;
CA.cluster_blend = 0.35;

% State codes.
CA.S = 0;
CA.E = 1;
CA.I_pre = 2;
CA.I_rash = 3;
CA.R = 4;
CA.D = 5;
CA.Vac = 6;

% ODE-derived vulnerability threshold for therapy eligibility.
CA.M_therapy_threshold = 0.005;

switch char(mode)
    case 'fast'
        CA.nruns_main = 10;
        CA.nruns_finite = 3;
    case 'paper'
        CA.nruns_main = 100;
        CA.nruns_finite = 30;
end
end

function print_ca_config(CA)
fprintf('\n--- Section 2: CA parameters and stochastic design ---\n');
fprintf('  Grid: N=%d (%d individuals), dt=%g day, T_max=%d days\n', CA.N, CA.N^2, CA.dt, CA.T_max);
fprintf('  Mode: %s | main nruns=%d | finite-size nruns=%d | base_seed=%d\n', ...
    CA.mode, CA.nruns_main, CA.nruns_finite, CA.base_seed);
fprintf('  Therapy threshold M_th = %.4f\n', CA.M_therapy_threshold);
fprintf('  Eligible baseline archetypes: ');
for w = 1:6
    if CA.profiles(w).M_emergent > CA.M_therapy_threshold
        fprintf('%s (M=%.3f) ', CA.profiles(w).name, CA.profiles(w).M_emergent);
    end
end
fprintf('\n  Balance check: S+E+Ipre+Irash+R+D+Vac=N^2 at every simulated day\n');
end

%% ========================================================================
%  Coupling map
% ========================================================================
function profiles = extract_coupled_profiles()
arch_names = {'Healthy adult','Child','Elderly','Immunocomp.','Partial vacc','Full vacc'};
Kv = 0.1;
V_th = 0.01;
rash_proxy_th = 0.30;
profiles = [];

for w = 1:6
    P = measles_params(w);
    [t,Y] = solve_ode(P);
    profiles = [profiles, build_profile(t,Y,P,arch_names{w},Kv,V_th,rash_proxy_th,w,false)]; %#ok<AGROW>
    p = profiles(w);
    fprintf(['  w=%d %-15s: T_viremia=%4.1f, T_inf_onset=%2d, ', ...
             'T_rash=%2d, T_inf_dur=%2d, T_clear=%2d, M=%.4f, peakV=%.2f\n'], ...
        w, arch_names{w}, p.T_viremia_onset, p.T_infectious_onset, ...
        p.T_rash, p.T_infectious_duration, p.T_clearance, p.M_emergent, p.peak_V);
end

tx_configs = {
    2,'Child+CTL',   struct('CTL_active',true,'CTL_delta',120,'CTL_time',8,'CTL_tau',1.3);
    4,'IC+DC+CTL',   struct('CTL_active',true,'CTL_delta',200,'CTL_time',10,'CTL_tau',1.5, ...
                            'DC_active',true,'eta_D',0.5,'eta_AbD',0.3, ...
                            'DC_doses',4,'DC_start',12,'DC_interval',3, ...
                            'DC_d0',1.0,'DC_ef',0.8,'DC_muD',0.3);
    3,'Elderly+CTL', struct('CTL_active',true,'CTL_delta',150,'CTL_time',10,'CTL_tau',1.5)
};

for ti = 1:size(tx_configs,1)
    w = tx_configs{ti,1};
    P = measles_params(w);
    flds = fieldnames(tx_configs{ti,3});
    for fi = 1:numel(flds)
        P.(flds{fi}) = tx_configs{ti,3}.(flds{fi});
    end
    [t,Y] = solve_ode(P);
    profiles = [profiles, build_profile(t,Y,P,tx_configs{ti,2},Kv,V_th,rash_proxy_th,w,true)]; %#ok<AGROW>
    p = profiles(end);
    fprintf('  w=%d+Tx %-12s: T_inf_dur=%2d, T_clear=%2d, M=%.4f (baseline M=%.4f)\n', ...
        w, tx_configs{ti,2}, p.T_infectious_duration, p.T_clearance, ...
        p.M_emergent, profiles(w).M_emergent);
end

p1 = profiles(1);
fprintf('\n  Infectivity burden diagnostic: sum(beta_daily)=%.2f days (not R0)\n', sum(p1.beta_daily));
end

function [t,Y] = solve_ode(P)
y0 = [P.S0; P.I0; P.Ag0; P.A17_0; P.V0; P.R0; P.Ab0];
opts = odeset('RelTol',1e-10,'AbsTol',1e-12,'MaxStep',0.1,'NonNegative',1:7);
[t,Y] = ode15s(@(t,y) measles_ode_rhs(t,y,P), [0 P.tmax_days], y0, opts);
end

function prof = build_profile(t,Y,P,name,Kv,V_th,rash_proxy_th,base_w,is_therapy)
V = Y(:,5);
Ag = Y(:,3);
beta_raw = V ./ (V + Kv);

idx_v = find(V > V_th, 1, 'first');
if isempty(idx_v)
    T_viremia_onset = NaN;
else
    T_viremia_onset = t(idx_v);
end

idx_clear = find(t >= 5 & V < V_th, 1, 'first');
if isempty(idx_clear)
    T_ode_clearance = P.tmax_days;
else
    T_ode_clearance = t(idx_clear);
end

rash_proxy = Ag ./ (Ag + P.hR);
idx_r = find(rash_proxy > rash_proxy_th, 1, 'first');
if isempty(idx_r)
    T_ode_rash_proxy = NaN;
else
    T_ode_rash_proxy = t(idx_r);
end

% Epidemiological timing used by the CA. The ODE viremia threshold is kept as
% an intra-host descriptor and is deliberately not mapped to incubation.
T_rash = epidemiological_rash_day(base_w, T_ode_rash_proxy);
T_infectious_onset = max(1, T_rash - 4);
T_infectious_end = T_rash + 4;

% Conservatively prolong shedding for immunocompromised hosts and partially
% reduce duration for treatment profiles when the ODE indicates earlier control.
if base_w == 4
    T_infectious_end = max(T_infectious_end, min(21, round(T_ode_clearance)));
end
if is_therapy
    baseline_end = T_infectious_end;
    ode_based_end = max(T_rash + 4, min(baseline_end, round(T_ode_clearance)));
    T_infectious_end = max(T_rash + 3, ode_based_end);
end

T_clearance = max(T_infectious_end, T_rash + 1);
T_infectious_duration = max(1, T_clearance - T_infectious_onset);

% Resample the ODE-derived infectivity shape onto the epidemiological
% infectious window. The temporal location is epidemiological; the shape and
% relative magnitude remain ODE-derived.
tau = (0:T_infectious_duration)';
if all(beta_raw == 0)
    beta_daily = zeros(size(tau));
else
    q = linspace(0,1,numel(beta_raw));
    qq = linspace(0,1,numel(tau));
    beta_daily = interp1(q, beta_raw, qq, 'pchip', 'extrap')';
    beta_daily = max(beta_daily, 0);
end

peakV = max(V);
peakAg = max(Ag);
auc_b = trapz(t, beta_raw);
M = emergent_mortality(peakV, T_infectious_duration, peakAg, auc_b);

prof = struct('name',name, ...
    'base_archetype',base_w, ...
    'is_therapy',logical(is_therapy), ...
    'T_viremia_onset',T_viremia_onset, ...
    'T_ode_rash_proxy',T_ode_rash_proxy, ...
    'T_ode_clearance',T_ode_clearance, ...
    'T_infectious_onset',round(T_infectious_onset), ...
    'T_rash',round(T_rash), ...
    'T_infectious_duration',round(T_infectious_duration), ...
    'T_clearance',round(T_clearance), ...
    'M_emergent',M, ...
    'beta_daily',beta_daily, ...
    'peak_V',peakV, ...
    'peak_Ag',peakAg, ...
    'auc_beta',auc_b);
end

function T_rash = epidemiological_rash_day(w, ode_proxy)
% Typical rash onset occurs about 10-14 days after exposure. The ODE proxy is
% allowed to fine-tune within this defensible range, but not outside it.
defaults = [11, 12, 12, 13, 11, 10];
lo = 10; hi = 14;
if w >= 1 && w <= numel(defaults)
    T_rash = defaults(w);
else
    T_rash = 12;
end
if ~isnan(ode_proxy)
    T_rash = round(0.70*T_rash + 0.30*min(max(ode_proxy, lo), hi));
end
T_rash = min(max(T_rash, lo), hi);
end

%% ========================================================================
%  Scenario runners
% ========================================================================
function results = run_vax_sweep(CA, vax_levels, n_runs)
results = struct([]);
for vi = 1:numel(vax_levels)
    vax = vax_levels(vi);
    scen = struct('scenario',sprintf('coverage_%02.0f',100*vax), ...
                  'scenario_type','coverage_sweep', ...
                  'base_vax',vax, 'use_therapy',false, 'events',{{}});
    tmp = run_ensemble(CA, scen, n_runs, vi);
    if vi == 1
        results = repmat(tmp, 1, numel(vax_levels));
    else
        results(vi) = tmp;
    end
    fprintf('  Vax=%3.0f%%: AR=%.2f +/- %.2f%%, deaths=%.2f +/- %.2f, peakI=%.1f +/- %.1f (n=%d)\n', ...
        100*vax, 100*mean(results(vi).attack_rate), 100*std(results(vi).attack_rate), ...
        mean(results(vi).deaths), std(results(vi).deaths), ...
        mean(results(vi).peak_I), std(results(vi).peak_I), n_runs);
end
end

function [tr, names] = run_therapy_comparison(CA, vax_cov, n_runs)
configs = {
    'No therapy',     false, false, false;
    'CTL->Child',     true,  false, false;
    'DC+CTL->IC',     false, true,  false;
    'CTL->Elderly',   false, false, true;
    'All vulnerable', true,  true,  true
};
names = configs(:,1)';
tr = struct([]);
for ti = 1:size(configs,1)
    scen = struct('scenario',configs{ti,1}, 'scenario_type','therapy', ...
        'base_vax',vax_cov, 'use_therapy',any([configs{ti,2:4}]), 'events',{{}}, ...
        'treat_child',configs{ti,2}, 'treat_IC',configs{ti,3}, 'treat_elderly',configs{ti,4});
    tmp = run_ensemble(CA, scen, n_runs, 100 + ti);
    if ti == 1
        tr = repmat(tmp, 1, size(configs,1));
    else
        tr(ti) = tmp;
    end
    fprintf('  %-15s: AR=%.2f +/- %.2f%%, deaths=%.2f +/- %.2f, peakI=%.1f +/- %.1f (n=%d)\n', ...
        names{ti}, 100*mean(tr(ti).attack_rate), 100*std(tr(ti).attack_rate), ...
        mean(tr(ti).deaths), std(tr(ti).deaths), mean(tr(ti).peak_I), std(tr(ti).peak_I), n_runs);
end
end

function [cr, names] = run_campaigns(CA, n_runs)
campaigns = {
    'Base 60%', struct('scenario','Base 60%', 'scenario_type','campaign', 'base_vax',0.60, 'use_therapy',false, 'events',{{}});
    'Ring d15', struct('scenario','Ring d15', 'scenario_type','campaign', 'base_vax',0.60, 'use_therapy',false, 'events',{{{15,'ring',0.90}}});
    'SIA d20 +20%', struct('scenario','SIA d20 +20%', 'scenario_type','campaign', 'base_vax',0.60, 'use_therapy',false, 'events',{{{20,'mass',0.20}}});
    'Child catch-up d10', struct('scenario','Child catch-up d10', 'scenario_type','campaign', 'base_vax',0.60, 'use_therapy',false, 'events',{{{10,'targeted',0.80,2}}});
    'Combined', struct('scenario','Combined', 'scenario_type','campaign', 'base_vax',0.60, 'use_therapy',false, 'events',{{{10,'targeted',0.80,2},{15,'ring',0.90}}});
    'Routine 92%', struct('scenario','Routine 92%', 'scenario_type','campaign', 'base_vax',0.92, 'use_therapy',false, 'events',{{}})
};
names = campaigns(:,1)';
cr = struct([]);
for ci = 1:size(campaigns,1)
    tmp = run_ensemble(CA, campaigns{ci,2}, n_runs, 200 + ci);
    if ci == 1
        cr = repmat(tmp, 1, size(campaigns,1));
    else
        cr(ci) = tmp;
    end
    fprintf('  %-20s: AR=%.2f +/- %.2f%%, deaths=%.2f +/- %.2f, peakI=%.1f +/- %.1f (n=%d)\n', ...
        names{ci}, 100*mean(cr(ci).attack_rate), 100*std(cr(ci).attack_rate), ...
        mean(cr(ci).deaths), std(cr(ci).deaths), mean(cr(ci).peak_I), std(cr(ci).peak_I), n_runs);
end
end

function out = run_ensemble(CA, scenario, n_runs, scenario_id)
T = CA.T_max;
ar = zeros(n_runs,1);
deaths = zeros(n_runs,1);
peakI = zeros(n_runs,1);
tpeak = zeros(n_runs,1);
curves = zeros(T,n_runs);
arch = zeros(6,5,n_runs);
balance_ok = true(n_runs,1);
seed_used = zeros(n_runs,1);

for r = 1:n_runs
    seed = CA.base_seed + scenario_id*CA.seed_stride + r;
    seed_used(r) = seed;
    rng(seed, 'twister');
    [epi,~,c,as,diagnostics] = run_ca_scenario(CA, scenario);
    ar(r) = epi.attack_rate;
    deaths(r) = epi.deaths;
    peakI(r) = epi.peak_I;
    tpeak(r) = epi.t_peak;
    curves(:,r) = c.I_total(:);
    balance_ok(r) = diagnostics.balance_ok;
    for w = 1:6
        arch(w,:,r) = [as(w).initial, as(w).infected, as(w).recovered, as(w).dead, as(w).treated];
    end
end

if ~all(balance_ok)
    error('Population balance failed in scenario %s.', scenario.scenario);
end

out = struct();
out.scenario = scenario.scenario;
out.scenario_type = scenario.scenario_type;
out.base_vax = scenario.base_vax;
out.use_therapy = scenario.use_therapy;
out.n_runs = n_runs;
out.seed_policy = sprintf('base_seed + scenario_id*%d + run_index', CA.seed_stride);
out.seeds = seed_used;
out.attack_rate = ar;
out.deaths = deaths;
out.peak_I = peakI;
out.t_peak = tpeak;
out.curves = curves;
out.arch_stats = mean(arch,3);
out.arch_stats_sd = std(arch,0,3);
out.balance_ok = all(balance_ok);
end

%% ========================================================================
%  Core CA engine
% ========================================================================
function [epi,lattice,curves,arch_stats,diagnostics] = run_ca_scenario(CA, scenario)
N = CA.N;
N2 = N*N;
profiles = CA.profiles;
[omega,state,timer] = init_pop(CA, scenario.base_vax);
therapy_applied = false(N,N);

initial_by_arch = accumarray(omega(:), 1, [6 1], @sum, 0);
infected_by_arch = zeros(6,1);
recovered_by_arch = zeros(6,1);
dead_by_arch = zeros(6,1);
treated_by_arch = zeros(6,1);

% Seed imported index infections among truly susceptible or vaccinated-but-susceptible hosts.
% Index cases are placed directly within the epidemiological infectious window
% with staggered infection ages. This avoids a synchronized artificial cohort in
% which all imported exposed cases become infectious on the same day.
eligible = find(state==CA.S | state==CA.Vac);
ns = min(CA.n_seeds, numel(eligible));
seed_idx = eligible(randperm(numel(eligible), ns));
for k = 1:numel(seed_idx)
    idx = seed_idx(k);
    w_seed = omega(idx);
    prof_seed = profiles(w_seed);
    lo_age = max(0, round(prof_seed.T_infectious_onset));
    hi_age = max(lo_age, round(prof_seed.T_clearance) - 1);
    age = lo_age + floor(rand() * (hi_age - lo_age + 1));
    timer(idx) = age;
    if age < prof_seed.T_rash
        state(idx) = CA.I_pre;
    else
        state(idx) = CA.I_rash;
    end
    infected_by_arch(w_seed) = infected_by_arch(w_seed) + 1;
end

It = zeros(CA.T_max,1);
Dd = zeros(CA.T_max,1);
Ni = zeros(CA.T_max,1);
balance_ok = true;

for day = 1:CA.T_max
    [state, omega] = apply_events_if_due(state, omega, day, scenario, CA);
    [sn,tn,therapy_applied,day_deaths,day_new_inf,infected_by_arch,recovered_by_arch,dead_by_arch,treated_by_arch] = ...
        update_one_day(CA,state,timer,omega,therapy_applied,profiles,scenario, ...
                       infected_by_arch,recovered_by_arch,dead_by_arch,treated_by_arch);
    state = sn;
    timer = tn;

    It(day) = nnz(state==CA.I_pre) + nnz(state==CA.I_rash);
    Dd(day) = day_deaths;
    Ni(day) = day_new_inf;

    if ~check_balance(state, CA, N2)
        balance_ok = false;
        error('Population balance failed on day %d in scenario %s.', day, scenario.scenario);
    end

    if day > 20 && (nnz(state==CA.E) + It(day)) == 0
        It(day+1:end) = 0;
        Dd(day+1:end) = 0;
        Ni(day+1:end) = 0;
        break;
    end
end

[pk,tp] = max(It);
epi = struct('attack_rate', sum(infected_by_arch)/N2, ...
             'deaths', sum(Dd), 'peak_I', pk, 't_peak', tp);
lattice = state;
curves = struct('I_total',It,'D_daily',Dd,'new_inf',Ni);
arch_stats = repmat(struct('initial',0,'infected',0,'recovered',0,'dead',0,'treated',0), 1, 6);
for w = 1:6
    arch_stats(w).initial = initial_by_arch(w);
    arch_stats(w).infected = infected_by_arch(w);
    arch_stats(w).recovered = recovered_by_arch(w);
    arch_stats(w).dead = dead_by_arch(w);
    arch_stats(w).treated = treated_by_arch(w);
end
diagnostics = struct('balance_ok',balance_ok);
end

function [sn,tn,therapy_applied,day_deaths,day_new_inf,infected_by_arch,recovered_by_arch,dead_by_arch,treated_by_arch] = ...
    update_one_day(CA,state,timer,omega,therapy_applied,profiles,scenario, ...
                   infected_by_arch,recovered_by_arch,dead_by_arch,treated_by_arch)
N = CA.N;
sn = state;
tn = timer;
day_deaths = 0;
day_new_inf = 0;
ifield = calc_ifield(N,state,timer,omega,therapy_applied,profiles,CA);

for i = 1:N
    for j = 1:N
        s = state(i,j);
        w = omega(i,j);
        tau = timer(i,j);
        switch s
            case {CA.S, CA.Vac}
                bs = ifield(i,j);
                if bs > 0
                    p_inf = 1 - exp(-CA.lambda_inf * bs * susceptibility_factor(w));
                    if rand() < p_inf
                        sn(i,j) = CA.E;
                        tn(i,j) = 0;
                        day_new_inf = day_new_inf + 1;
                        infected_by_arch(w) = infected_by_arch(w) + 1;
                    end
                end

            case CA.E
                tn(i,j) = tau + CA.dt;
                if tn(i,j) >= profiles(w).T_infectious_onset
                    sn(i,j) = CA.I_pre;
                end

            case CA.I_pre
                tn(i,j) = tau + CA.dt;
                if tn(i,j) >= profiles(w).T_rash
                    sn(i,j) = CA.I_rash;
                    if scenario.use_therapy && is_therapy_eligible(w, CA, scenario)
                        therapy_applied(i,j) = true;
                        treated_by_arch(w) = treated_by_arch(w) + 1;
                    end
                end

            case CA.I_rash
                tn(i,j) = tau + CA.dt;
                prof = get_prof(w, therapy_applied(i,j), profiles);
                if tn(i,j) >= prof.T_clearance
                    if rand() < prof.M_emergent
                        sn(i,j) = CA.D;
                        day_deaths = day_deaths + 1;
                        dead_by_arch(w) = dead_by_arch(w) + 1;
                    else
                        sn(i,j) = CA.R;
                        recovered_by_arch(w) = recovered_by_arch(w) + 1;
                    end
                end
        end
    end
end
end

function tf = is_therapy_eligible(w, CA, scenario)
tf = false;
if CA.profiles(w).M_emergent <= CA.M_therapy_threshold
    return;
end
if isfield(scenario,'treat_child')
    tf = (w==2 && scenario.treat_child) || ...
         (w==4 && scenario.treat_IC) || ...
         (w==3 && scenario.treat_elderly);
else
    tf = (w==2 || w==3 || w==4);
end
end

function sf = susceptibility_factor(w)
switch w
    case 5
        sf = 0.35;
    case 6
        sf = 0.08;
    otherwise
        sf = 1.00;
end
end

function prof = get_prof(w, treated, profiles)
if treated
    if w == 2
        prof = profiles(7); return;
    elseif w == 4
        prof = profiles(8); return;
    elseif w == 3
        prof = profiles(9); return;
    end
end
prof = profiles(w);
end

function ifield = calc_ifield(N,state,timer,omega,therapy_applied,profiles,CA)
ifield = zeros(N,N);
for di = -1:1
    for dj = -1:1
        if di == 0 && dj == 0, continue; end
        src_i = mod((1:N) + di - 1, N) + 1;
        src_j = mod((1:N) + dj - 1, N) + 1;
        for ii = 1:N
            for jj = 1:N
                ns = state(src_i(ii), src_j(jj));
                if ns == CA.I_pre || ns == CA.I_rash
                    w = omega(src_i(ii), src_j(jj));
                    tau = timer(src_i(ii), src_j(jj));
                    prof = get_prof(w, therapy_applied(src_i(ii), src_j(jj)), profiles);
                    bd = prof.beta_daily;
                    k = min(max(round(tau - prof.T_infectious_onset) + 1, 1), numel(bd));
                    ifield(ii,jj) = ifield(ii,jj) + bd(k);
                end
            end
        end
    end
end
end

%% ========================================================================
%  Initial population and interventions
% ========================================================================
function [omega,state,timer] = init_pop(CA,vax_cov)
N = CA.N;
N2 = N*N;
demo = CA.host_demo(:) / sum(CA.host_demo);
cumdemo = cumsum(demo);
r = rand(N2,1);
ov = ones(N2,1);
for k = 2:numel(demo)
    ov(r > cumdemo(k-1)) = k;
end
omega = reshape(ov,N,N);

if vax_cov <= 0
    vaccinated = false(N,N);
elseif vax_cov >= 1
    vaccinated = true(N,N);
else
    r1 = randn(N,N);
    r2 = rand(N,N);
    rad = max(2, ceil(3*CA.cluster_sigma));
    [x,y] = meshgrid(-rad:rad, -rad:rad);
    K = exp(-(x.^2 + y.^2)/(2*CA.cluster_sigma^2));
    K = K / sum(K(:));
    spatial_field = conv2(r1, K, 'same');
    spatial_field = (spatial_field - mean(spatial_field(:))) / (std(spatial_field(:)) + eps);
    random_field = (r2 - mean(r2(:))) / (std(r2(:)) + eps);
    f = CA.cluster_blend*spatial_field + (1-CA.cluster_blend)*random_field;
    vaccinated = f >= quantile(f(:), 1-vax_cov);
end

vi = find(vaccinated);
nv = numel(vi);
np = round(CA.partial_frac * nv);
if nv > 0
    pm = randperm(nv);
    partial_idx = vi(pm(1:np));
    full_idx = vi(pm(np+1:end));
    omega(partial_idx) = 5;
    omega(full_idx) = 6;
end

state = zeros(N,N);
state(vaccinated) = CA.Vac;
timer = zeros(N,N);
end

function [state,omega] = apply_events_if_due(state, omega, day, scenario, CA)
if ~isfield(scenario,'events') || isempty(scenario.events)
    return;
end
for ei = 1:numel(scenario.events)
    ev = scenario.events{ei};
    if day == ev{1}
        [state, omega] = apply_event(ev, state, omega, CA.N, CA); %#ok<ASGLU>
    end
end
end

function [state,omega] = apply_event(ev,state,omega,N,CA)
switch ev{2}
    case 'mass'
        eligible = find(state==CA.S);
        nv = round(ev{3} * numel(eligible));
        if nv > 0
            pick = eligible(randperm(numel(eligible), min(nv,numel(eligible))));
            state(pick) = CA.Vac;
            omega(pick) = 6;
        end

    case 'ring'
        ring_target = false(N,N);
        for i = 1:N
            for j = 1:N
                if state(i,j)==CA.I_rash || state(i,j)==CA.I_pre
                    for di = -2:2
                        for dj = -2:2
                            ni = mod(i+di-1,N)+1;
                            nj = mod(j+dj-1,N)+1;
                            if state(ni,nj)==CA.S
                                ring_target(ni,nj) = true;
                            end
                        end
                    end
                end
            end
        end
        ri = find(ring_target);
        nv = round(ev{3} * numel(ri));
        if nv > 0
            pick = ri(randperm(numel(ri), min(nv,numel(ri))));
            state(pick) = CA.Vac;
            omega(pick) = 6;
        end

    case 'targeted'
        target_w = ev{4};
        eligible = find(state==CA.S & omega==target_w);
        nv = round(ev{3} * numel(eligible));
        if nv > 0
            pick = eligible(randperm(numel(eligible), min(nv,numel(eligible))));
            state(pick) = CA.Vac;
            omega(pick) = 6;
        end

    otherwise
        error('Unknown intervention event type: %s', ev{2});
end
end

function ok = check_balance(state, CA, N2)
total = nnz(state==CA.S) + nnz(state==CA.E) + nnz(state==CA.I_pre) + ...
        nnz(state==CA.I_rash) + nnz(state==CA.R) + nnz(state==CA.D) + nnz(state==CA.Vac);
ok = (total == N2);
end

%% ========================================================================
%  Snapshots and figures
% ========================================================================
function snapshots = run_ca_snapshots(CA, cspec, snap_days)
N = CA.N;
profiles = CA.profiles;
[omega,state,timer] = init_pop(CA, cspec.base_vax);
therapy_applied = false(N,N);
eligible = find(state==CA.S | state==CA.Vac);
seed_idx = eligible(randperm(numel(eligible), min(CA.n_seeds,numel(eligible))));
state(seed_idx) = CA.E;
timer(seed_idx) = 0;
snapshots = cell(numel(snap_days),1);
if ismember(0,snap_days)
    snapshots{snap_days==0} = state;
    fprintf('    Snapshot day 0\n');
end
for day = 1:max(snap_days)
    [state, omega] = apply_events_if_due(state, omega, day, cspec, CA);
    scenario = cspec;
    dummy = zeros(6,1);
    [state,timer,therapy_applied,~,~,~,~,~,~] = update_one_day(CA,state,timer,omega,therapy_applied,profiles,scenario,dummy,dummy,dummy,dummy);
    if ismember(day,snap_days)
        snapshots{snap_days==day} = state;
        fprintf('    Snapshot day %d\n', day);
    end
end
end

function make_figures(CA, vax_levels, results, therapy_results, therapy_names, campaign_results, campaign_names, snapshots, snap_days)
coverage = build_coverage_table(vax_levels, results);

fig11 = figure('Color','w','Position',[50 50 1000 420]);
subplot(1,2,1); hold on;
errorbar(coverage.Vax_percent, coverage.AR_mean_percent, coverage.AR_sd_percent, 'o-', 'LineWidth', 1.8, 'MarkerSize', 6);
yline(5,'--','LineWidth',1.0);
xlabel('Vaccination coverage (%)'); ylabel('Attack rate (% of total population)');
title('(a) Attack rate'); grid on; xlim([-2 97]);
subplot(1,2,2); hold on;
errorbar(coverage.Vax_percent, coverage.Deaths_mean, coverage.Deaths_sd, 's-', 'LineWidth', 1.8, 'MarkerSize', 6);
xlabel('Vaccination coverage (%)'); ylabel('Deaths, mean +/- SD');
title('(b) Deaths'); grid on; xlim([-2 97]);
sgtitle('Figure 11. Vaccination coverage sweep','FontWeight','bold');
hide_tb(fig11);
exportgraphics(fig11, fullfile('export','fig11_vax_sweep.png'), 'Resolution', 300);

fig12 = figure('Color','w','Position',[90 50 1100 500]);
sel = [1,2,4,6,8];
subplot(1,2,1); hold on; leg = cell(numel(sel),1);
for j = 1:numel(sel)
    mc = mean(results(sel(j)).curves,2);
    plot(1:CA.T_max, mc, '-', 'LineWidth', 1.8);
    leg{j} = sprintf('%.0f%%', vax_levels(sel(j))*100);
end
xlabel('Days'); ylabel('Active infections'); title('(a) Linear scale');
legend(leg,'Location','northeast'); grid on; xlim([0 CA.T_max]);
subplot(1,2,2); hold on;
for j = 1:numel(sel)
    mc = mean(results(sel(j)).curves,2);
    semilogy(1:CA.T_max, max(mc,0.1), '-', 'LineWidth', 1.8);
end
xlabel('Days'); ylabel('Active infections (log)'); title('(b) Log scale');
legend(leg,'Location','northeast'); grid on; xlim([0 CA.T_max]);
sgtitle('Figure 12. Epidemic curves by coverage','FontWeight','bold');
hide_tb(fig12);
exportgraphics(fig12, fullfile('export','fig12_epidemic_curves.png'), 'Resolution', 300);

fig13 = figure('Color','w','Position',[140 50 1100 450]);
subplot(1,2,1); hold on;
for ti = 1:numel(therapy_names)
    plot(1:CA.T_max, mean(therapy_results(ti).curves,2), '-', 'LineWidth', 1.7);
end
xlabel('Days'); ylabel('Active infections'); title('(a) Epidemic curves');
legend(therapy_names,'Location','northeast','FontSize',7); grid on; xlim([0 CA.T_max]);
subplot(1,2,2);
bar(arrayfun(@(r) mean(r.deaths), therapy_results));
set(gca,'XTickLabel',therapy_names,'XTickLabelRotation',30,'FontSize',7);
ylabel('Deaths, ensemble mean'); title('(b) Mortality'); grid on;
sgtitle('Figure 13. Therapy strategies at 60% baseline coverage','FontWeight','bold');
hide_tb(fig13);
exportgraphics(fig13, fullfile('export','fig13_therapy_comparison.png'), 'Resolution', 300);

fig14 = figure('Color','w','Position',[180 50 1100 450]);
subplot(1,2,1); hold on;
for ci = 1:numel(campaign_names)
    plot(1:CA.T_max, mean(campaign_results(ci).curves,2), '-', 'LineWidth', 1.7);
end
xlabel('Days'); ylabel('Active infections'); title('(a) Campaign trajectories');
legend(campaign_names,'Location','northeast','FontSize',7); grid on; xlim([0 CA.T_max]);
subplot(1,2,2); hold on;
camp_ar = arrayfun(@(r) mean(r.attack_rate), campaign_results) * 100;
camp_d = arrayfun(@(r) mean(r.deaths), campaign_results);
yyaxis left; bar(1:numel(campaign_names), camp_ar, 0.45); ylabel('Attack rate (%)');
yyaxis right; plot(1:numel(campaign_names), camp_d, 's-', 'LineWidth', 1.8); ylabel('Deaths, mean');
set(gca,'XTick',1:numel(campaign_names),'XTickLabel',campaign_names,'XTickLabelRotation',25,'FontSize',7);
title('(b) Outcomes'); grid on;
sgtitle('Figure 14. Vaccination campaign comparison','FontWeight','bold');
hide_tb(fig14);
exportgraphics(fig14, fullfile('export','fig14_campaigns.png'), 'Resolution', 300);

fig15 = figure('Color','w','Position',[200 30 1200 700]);
cmap = [0.90 0.95 1.00; 1.00 0.85 0.60; 1.00 0.55 0.20; 0.85 0.15 0.15; 0.30 0.70 0.30; 0.10 0.10 0.10; 0.60 0.75 0.95];
for si = 1:numel(snap_days)
    subplot(2,3,si);
    imagesc(snapshots{si});
    colormap(cmap); caxis([-0.5 6.5]);
    title(sprintf('Day %d', snap_days(si)));
    axis square; set(gca,'XTick',[],'YTick',[]);
end
cb = colorbar('Ticks',0:6,'TickLabels',{'S','E','I_{pre}','I_{rash}','R','D','Vac'});
cb.Position = [0.93 0.11 0.02 0.82];
sgtitle('Figure 15. Lattice snapshots','FontWeight','bold');
hide_tb(fig15);
exportgraphics(fig15, fullfile('export','fig15_lattice_snapshots.png'), 'Resolution', 300);
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

%% ========================================================================
%  Tables and summaries
% ========================================================================
function T = build_coverage_table(vax_levels, results)
n = numel(vax_levels);
Vax_percent = 100*vax_levels(:);
N_runs = zeros(n,1);
AR_mean_percent = zeros(n,1); AR_sd_percent = zeros(n,1); AR_ci95_low_percent = zeros(n,1); AR_ci95_high_percent = zeros(n,1);
Deaths_mean = zeros(n,1); Deaths_sd = zeros(n,1); Deaths_ci95_low = zeros(n,1); Deaths_ci95_high = zeros(n,1);
PeakI_mean = zeros(n,1); PeakI_sd = zeros(n,1); PeakI_ci95_low = zeros(n,1); PeakI_ci95_high = zeros(n,1);
for i = 1:n
    N_runs(i) = results(i).n_runs;
    [m,s,lo,hi] = mean_sd_ci(100*results(i).attack_rate);
    AR_mean_percent(i)=m; AR_sd_percent(i)=s; AR_ci95_low_percent(i)=lo; AR_ci95_high_percent(i)=hi;
    [m,s,lo,hi] = mean_sd_ci(results(i).deaths);
    Deaths_mean(i)=m; Deaths_sd(i)=s; Deaths_ci95_low(i)=lo; Deaths_ci95_high(i)=hi;
    [m,s,lo,hi] = mean_sd_ci(results(i).peak_I);
    PeakI_mean(i)=m; PeakI_sd(i)=s; PeakI_ci95_low(i)=lo; PeakI_ci95_high(i)=hi;
end
T = table(Vax_percent,N_runs,AR_mean_percent,AR_sd_percent,AR_ci95_low_percent,AR_ci95_high_percent, ...
    Deaths_mean,Deaths_sd,Deaths_ci95_low,Deaths_ci95_high, ...
    PeakI_mean,PeakI_sd,PeakI_ci95_low,PeakI_ci95_high);
end

function T = build_named_table(names, results, first_name)
n = numel(names);
Name = string(names(:));
N_runs = zeros(n,1);
Base_coverage_percent = zeros(n,1);
AR_mean_percent = zeros(n,1); AR_sd_percent = zeros(n,1); AR_ci95_low_percent = zeros(n,1); AR_ci95_high_percent = zeros(n,1);
Deaths_mean = zeros(n,1); Deaths_sd = zeros(n,1); Deaths_ci95_low = zeros(n,1); Deaths_ci95_high = zeros(n,1);
PeakI_mean = zeros(n,1); PeakI_sd = zeros(n,1); PeakI_ci95_low = zeros(n,1); PeakI_ci95_high = zeros(n,1);
for i = 1:n
    N_runs(i) = results(i).n_runs;
    Base_coverage_percent(i) = 100*results(i).base_vax;
    [m,s,lo,hi] = mean_sd_ci(100*results(i).attack_rate);
    AR_mean_percent(i)=m; AR_sd_percent(i)=s; AR_ci95_low_percent(i)=lo; AR_ci95_high_percent(i)=hi;
    [m,s,lo,hi] = mean_sd_ci(results(i).deaths);
    Deaths_mean(i)=m; Deaths_sd(i)=s; Deaths_ci95_low(i)=lo; Deaths_ci95_high(i)=hi;
    [m,s,lo,hi] = mean_sd_ci(results(i).peak_I);
    PeakI_mean(i)=m; PeakI_sd(i)=s; PeakI_ci95_low(i)=lo; PeakI_ci95_high(i)=hi;
end
T = table(Name,N_runs,Base_coverage_percent,AR_mean_percent,AR_sd_percent,AR_ci95_low_percent,AR_ci95_high_percent, ...
    Deaths_mean,Deaths_sd,Deaths_ci95_low,Deaths_ci95_high,PeakI_mean,PeakI_sd,PeakI_ci95_low,PeakI_ci95_high);
T.Properties.VariableNames{1} = first_name;
end

function T = build_archetype_table(arch_names, arch_mean, profiles)
Archetype = string(arch_names(:));
Pop_mean = arch_mean(:,1);
Infected_mean = arch_mean(:,2);
Recovered_mean = arch_mean(:,3);
Deaths_mean = arch_mean(:,4);
Treated_mean = arch_mean(:,5);
AR_mean_percent = 100 * Infected_mean ./ max(Pop_mean, 1);
CFR_mean_percent = 100 * Deaths_mean ./ max(Infected_mean, 1);
M_emergent = arrayfun(@(p) p.M_emergent, profiles(1:6))';
T = table(Archetype,Pop_mean,Infected_mean,Recovered_mean,Deaths_mean,Treated_mean,AR_mean_percent,CFR_mean_percent,M_emergent);
end

function T = build_profiles_table(profiles)
n = numel(profiles);
Name = strings(n,1); Base_archetype = zeros(n,1); Is_therapy = false(n,1);
T_viremia_onset = zeros(n,1); T_infectious_onset = zeros(n,1); T_rash = zeros(n,1);
T_infectious_duration = zeros(n,1); T_clearance = zeros(n,1); T_ode_clearance = zeros(n,1);
M_emergent = zeros(n,1); Peak_V = zeros(n,1); Peak_Ag = zeros(n,1); AUC_beta = zeros(n,1);
for i = 1:n
    Name(i) = string(profiles(i).name);
    Base_archetype(i) = profiles(i).base_archetype;
    Is_therapy(i) = profiles(i).is_therapy;
    T_viremia_onset(i) = profiles(i).T_viremia_onset;
    T_infectious_onset(i) = profiles(i).T_infectious_onset;
    T_rash(i) = profiles(i).T_rash;
    T_infectious_duration(i) = profiles(i).T_infectious_duration;
    T_clearance(i) = profiles(i).T_clearance;
    T_ode_clearance(i) = profiles(i).T_ode_clearance;
    M_emergent(i) = profiles(i).M_emergent;
    Peak_V(i) = profiles(i).peak_V;
    Peak_Ag(i) = profiles(i).peak_Ag;
    AUC_beta(i) = profiles(i).auc_beta;
end
T = table(Name,Base_archetype,Is_therapy,T_viremia_onset,T_infectious_onset,T_rash,T_infectious_duration,T_clearance,T_ode_clearance,M_emergent,Peak_V,Peak_Ag,AUC_beta);
end

function T = build_metadata_table(CA, mode)
Key = strings(12,1);
Value = strings(12,1);
Key(1) = 'mode'; Value(1) = string(mode);
Key(2) = 'N_main'; Value(2) = string(CA.N);
Key(3) = 'N_finite'; Value(3) = string(CA.N_finite);
Key(4) = 'nruns_main'; Value(4) = string(CA.nruns_main);
Key(5) = 'nruns_finite'; Value(5) = string(CA.nruns_finite);
Key(6) = 'n_seeds'; Value(6) = string(CA.n_seeds);
Key(7) = 'seed_initial_stage'; Value(7) = string(CA.seed_initial_stage);
Key(8) = 'lambda_inf'; Value(8) = string(CA.lambda_inf);
Key(6) = 'base_seed'; Value(6) = string(CA.base_seed);
Key(7) = 'seed_policy'; Value(7) = sprintf('base_seed + scenario_id*%d + run_index', CA.seed_stride);
Key(8) = 'attack_rate_denominator'; Value(8) = 'total population N^2';
Key(9) = 'balance_equation'; Value(9) = 'S+E+Ipre+Irash+R+D+Vac=N^2';
Key(10) = 'viremia_threshold_role'; Value(10) = 'descriptor only; not CA latent period';
Key(11) = 'infectious_window'; Value(11) = 'T_rash-4 to T_rash+4, with conservative IC prolongation';
Key(12) = 'therapy_trigger'; Value(12) = sprintf('M_emergent > %.4f', CA.M_therapy_threshold);
T = table(Key,Value);
end

function [m,s,lo,hi] = mean_sd_ci(x)
x = x(:);
n = numel(x);
m = mean(x);
s = std(x);
if n > 1
    half = 1.96 * s / sqrt(n);
else
    half = NaN;
end
lo = m - half;
hi = m + half;
end

function print_profiles_table(T)
fprintf('\n================================================================\n');
fprintf('  TABLE 4.1: Corrected Coupling Map Timings\n');
fprintf('================================================================\n');
fprintf('  %-15s %7s %7s %7s %7s %7s %9s\n', 'Profile','Viremia','InfOn','Rash','InfDur','Clear','M');
fprintf('  %s\n', repmat('-',1,75));
for i = 1:height(T)
    fprintf('  %-15s %7.1f %7.0f %7.0f %7.0f %7.0f %9.4f\n', ...
        T.Name(i), T.T_viremia_onset(i), T.T_infectious_onset(i), T.T_rash(i), ...
        T.T_infectious_duration(i), T.T_clearance(i), T.M_emergent(i));
end
end

function print_coverage_table(T)
fprintf('\n================================================================\n');
fprintf('  TABLE D: Coverage Sweep\n');
fprintf('================================================================\n');
fprintf('  %-8s %7s %17s %17s %17s\n','Vax(%)','n','AR mean+/-SD','Deaths mean+/-SD','PeakI mean+/-SD');
fprintf('  %s\n', repmat('-',1,80));
for i = 1:height(T)
    fprintf('  %-8.0f %7d %8.2f +/- %-6.2f %8.2f +/- %-6.2f %8.1f +/- %-6.1f\n', ...
        T.Vax_percent(i), T.N_runs(i), T.AR_mean_percent(i), T.AR_sd_percent(i), ...
        T.Deaths_mean(i), T.Deaths_sd(i), T.PeakI_mean(i), T.PeakI_sd(i));
end
end

function print_named_outcomes(T, title_text)
fprintf('\n================================================================\n');
fprintf('  %s\n', title_text);
fprintf('================================================================\n');
first = T.Properties.VariableNames{1};
fprintf('  %-22s %7s %17s %17s %17s\n', first, 'n', 'AR mean+/-SD', 'Deaths mean+/-SD', 'PeakI mean+/-SD');
fprintf('  %s\n', repmat('-',1,90));
for i = 1:height(T)
    nameval = string(T{i,1});
    fprintf('  %-22s %7d %8.2f +/- %-6.2f %8.2f +/- %-6.2f %8.1f +/- %-6.1f\n', ...
        nameval, T.N_runs(i), T.AR_mean_percent(i), T.AR_sd_percent(i), ...
        T.Deaths_mean(i), T.Deaths_sd(i), T.PeakI_mean(i), T.PeakI_sd(i));
end
end

function print_archetype_table(T)
fprintf('\n================================================================\n');
fprintf('  TABLE A: Per-Archetype Outcomes (80%% coverage, no therapy)\n');
fprintf('================================================================\n');
fprintf('  %-15s %9s %10s %10s %10s %9s %9s %9s\n', ...
    'Archetype','Pop','Inf','Rec','Deaths','Treated','AR(%)','CFR(%)');
fprintf('  %s\n', repmat('-',1,95));
for i = 1:height(T)
    fprintf('  %-15s %9.1f %10.1f %10.1f %10.2f %9.1f %9.2f %9.2f\n', ...
        T.Archetype(i), T.Pop_mean(i), T.Infected_mean(i), T.Recovered_mean(i), ...
        T.Deaths_mean(i), T.Treated_mean(i), T.AR_mean_percent(i), T.CFR_mean_percent(i));
end
end

function print_finite_size_table(vax_levels, results, vax_check, results_fs, CA, CAfs)
fprintf('\n  Finite-size comparison (main n=%d, finite n=%d):\n', CA.nruns_main, CAfs.nruns_main);
fprintf('  %-8s | %18s %18s | %18s %18s\n', 'Vax(%)', 'AR main mean', 'AR finite mean', 'Deaths main', 'Deaths finite');
fprintf('  %s\n', repmat('-',1,95));
for vi = 1:numel(vax_check)
    vax = vax_check(vi);
    idx = find(abs(vax_levels - vax) < 1e-12, 1);
    fprintf('  %-8.0f | %16.2f%% %16.2f%% | %16.2f %16.2f\n', ...
        100*vax, 100*mean(results(idx).attack_rate), 100*mean(results_fs(vi).attack_rate), ...
        mean(results(idx).deaths), mean(results_fs(vi).deaths));
end
end
