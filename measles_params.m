function P = measles_params(archetype)
% MEASLES_PARAMS  Parameter sets for within-host measles model.
%
%   P = measles_params(archetype)
%
%   archetype: 1=Healthy adult, 2=Child, 3=Elderly,
%              4=Immunocompromised, 5=Partially vaccinated,
%              6=Fully vaccinated / seroprotected breakthrough profile
%
% The fields P.param_class and P.param_source_label classify every parameter
% into one of the manuscript categories:
%   inherited_morris
%   calibrated
%   phenomenological
%   exploratory_archetype_modifier
%   in_silico_intervention_parameter
%
% This classification is exported by export_parameter_catalog.m and by the
% phase scripts. Archetype-specific parameter changes are exploratory
% mechanistic profiles, not clinical cohort calibrations.

if nargin < 1, archetype = 1; end

%% ---- Core parameters: inherited from Morris et al. (2018) where applicable ----
P.beta  = 0.040;   P.k     = 0.017;   P.p     = 0.089;
P.q     = 0.99;    P.qs    = 0.028;   P.r     = 0.058;
P.s     = 0.0104;  P.c     = 3.0;     P.d     = 1/40;
P.delta = 1/2;     P.td    = 4.2;

%% ---- Extension parameters: calibrated/phenomenological additions ----
P.mu_R     = log(2)/50;  P.rho_R    = 0.08;     P.eta_R   = 2e-3;
P.t_Ab     = 10;         P.rho_Ab   = 0.015;    P.mu_Ab   = 1/120;
P.theta_Ab = 1e-3;       P.t17      = 35;       P.alpha_t = 0.15;
P.alpha_R  = 0.6;        P.kappa17  = 0.02;     P.hR      = 50;

% Antibody stimulation function in Eq. (4g). Default preserves v8, while
% alternative functions are used in the v9 sensitivity analysis.
P.Ab_stim_type = 'log';       % 'log', 'michaelis_menten', or 'hill'
P.K_R = 100;                  % half-saturation for RNA-driven antibody stimulation
P.n_R = 2;                    % Hill coefficient

%% ---- Therapy parameters: inactive by default ----
P.DC_active = false;
P.DC_mode   = 'multiplicative';  % 'multiplicative' or 'additive_control'
P.eta_D = 0;       P.eta_AbD = 0;
P.DC_doses  = 0;   P.DC_start = 0;    P.DC_interval = 0;
P.DC_d0     = 0;   P.DC_ef = 0;       P.DC_muD = 0;
% Additive-control fields are only used for the architecture diagnostic.
P.DC_add_q      = 0;     % absolute addition to q_eff per unit D(t)
P.DC_add_rhoAb  = 0;     % absolute addition to rho_Ab_eff per unit D(t)
P.DC_add_Ag     = 0;     % external Ag recruitment/stimulation per unit D(t)

P.CTL_active = false; P.CTL_delta = 0;   P.CTL_time = 0;
P.CTL_tau    = 1;

%% ---- Initial conditions ----
P.S0 = 3800;  P.I0 = 0;  P.V0 = 1e-4;
P.A0 = 1;          % 4-variable mode
P.Ag0 = 1;         % 7-variable mode: IFN-gamma-dominant T cells
P.A17_0 = 0;       % 7-variable mode: IL-17-dominant T cells
P.R0 = 0;          % persistent RNA
P.Ab0 = 0;         % normalized antibody proxy

%% ---- Simulation ----
P.tmax_days = 120;
P.archetype = archetype;
P.profile_type = 'baseline mechanistic archetype';
P.profile_interpretation = 'Exploratory mechanistic profile, not a calibrated clinical cohort.';

%% ---- Archetype modifications ----
switch archetype
    case 1
        P.name = 'Healthy adult';
        P.profile_type = 'reference archetype';
        P.profile_interpretation = 'Reference host using Morris-core baseline parameters with extension defaults.';

    case 2
        P.name = 'Child (<5 y)';
        P.S0     = 3000;
        P.q      = 0.68;
        P.k      = 0.011;
        P.rho_Ab = 0.008;
        P.delta  = 0.45;
        P.td     = 4.8;

    case 3
        P.name = 'Elderly (>65 y)';
        P.q      = 0.75;
        P.k      = 0.013;
        P.rho_Ab = 0.008;
        P.mu_Ab  = 1/80;
        P.delta  = 0.45;
        P.td     = 4.8;

    case 4
        P.name = 'Immunocompromised';
        P.S0       = 3000;
        P.q        = 0.20;
        P.k        = 0.003;
        P.rho_Ab   = 0.005;
        P.delta    = 0.30;
        P.tmax_days = 180;
        P.profile_interpretation = ['Exploratory severe immune-dysfunction profile; ', ...
            'target cells remain available but effector proliferation/killing are strongly impaired.'];

    case 5
        P.name = 'Partial vaccine';
        P.V0     = 1e-7;
        P.Ab0    = 3.0;
        P.Ag0    = 1;
        P.A0     = 2.5;
        P.k      = 0.020;
        P.rho_Ab = 0.025;
        P.td     = 3.7;
        P.profile_interpretation = ['Exploratory incomplete-protection breakthrough profile ', ...
            'with pre-existing antibody and faster recall.'];

    case 6
        P.name = 'Fully vaccinated';
        P.profile_type = 'seroprotected breakthrough profile';
        P.profile_interpretation = ['Strongly seroprotected breakthrough profile, not sterilizing immunity; ', ...
            'used to represent attenuated infection under high pre-existing humoral protection.'];
        P.V0     = 1e-8;
        P.Ab0    = 30;
        P.Ag0    = 3;
        P.A17_0  = 0;
        P.A0     = 3;
        P.k      = 0.055;
        P.rho_Ab = 0.080;
        P.theta_Ab = 0.030;
        P.t17     = 10;
        P.alpha_t = 0.10;
        P.alpha_R = 0.030;
        P.rho_R   = 0.005;
        P.eta_R   = 0.050;
        P.td      = 2.2;
        P.q       = 0.95;

    otherwise
        error('Unknown archetype: %d. Use 1-6.', archetype);
end

P = attach_parameter_metadata(P);
end

function P = attach_parameter_metadata(P)
% Attach source/type metadata as structs for lightweight backward compatibility.

names = fieldnames(P);
for i = 1:numel(names)
    P.param_class.(names{i}) = 'not_parameter_or_metadata';
    P.param_source_label.(names{i}) = 'internal/model-control field';
end

inherited = {'beta','k','p','q','qs','r','s','c','d','delta','td','S0','A0','I0','V0'};
calibrated = {'rho_R','mu_R','rho_Ab','mu_Ab','theta_Ab','kappa17','t17'};
phenom = {'eta_R','t_Ab','alpha_t','alpha_R','hR','Ag0','A17_0','R0','Ab0','Ab_stim_type','K_R','n_R'};
intervention = {'DC_active','DC_mode','eta_D','eta_AbD','DC_doses','DC_start', ...
    'DC_interval','DC_d0','DC_ef','DC_muD','DC_add_q','DC_add_rhoAb', ...
    'DC_add_Ag','CTL_active','CTL_delta','CTL_time','CTL_tau'};

P = set_class(P, inherited, 'inherited_morris', 'Inherited from Morris et al. core or required core initial condition.');
P = set_class(P, calibrated, 'calibrated', 'Calibrated or literature-anchored extension parameter.');
P = set_class(P, phenom, 'phenomenological', 'Phenomenological proxy or extension parameter.');
P = set_class(P, intervention, 'in_silico_intervention_parameter', 'In silico intervention parameter.');

% Any parameter that differs from the healthy-adult default in archetypes 2-6
% is treated in exported tables as an exploratory archetype modifier.
if isfield(P, 'archetype') && P.archetype ~= 1
    P.archetype_modifier_note = 'Values changed from the healthy-adult reference are exploratory archetype modifiers.';
else
    P.archetype_modifier_note = 'Reference archetype; no exploratory archetype modifier applied.';
end
end

function P = set_class(P, names, cls, source_label)
for i = 1:numel(names)
    nm = names{i};
    if isfield(P, nm)
        P.param_class.(nm) = cls;
        P.param_source_label.(nm) = source_label;
    end
end
end
