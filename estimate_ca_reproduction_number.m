function T = estimate_ca_reproduction_number(mode)
% ESTIMATE_CA_REPRODUCTION_NUMBER  Direct CA reproduction-scale estimate.
%
%   T = estimate_ca_reproduction_number('fast')
%   T = estimate_ca_reproduction_number('paper')
%
% This replaces the weak interpretation of sum(beta_daily) as R0. It simulates
% one index case in an otherwise susceptible lattice and counts secondary
% infections directly caused during the index infectious window.

if nargin < 1 || isempty(mode), mode = 'paper'; end
mode = lower(string(mode));
switch mode
    case 'fast'
        nrep = 100; N = 75;
    case 'paper'
        nrep = 1000; N = 150;
    otherwise
        error('Unknown mode: %s', mode);
end
if ~exist('export','dir'), mkdir('export'); end

% Healthy-adult profile approximation kept independent from private Phase 4
% nested functions. It uses the same daily window and infectivity transform.
P = measles_params(1);
[t,Y] = solve7(P);
V = Y(:,5); beta_raw = V./(V+0.1);
T_rash = 12; T_on = T_rash-4; T_clear = T_rash+4; Tdur = T_clear-T_on;
q = linspace(0,1,numel(beta_raw)); qq = linspace(0,1,Tdur+1);
beta_daily = max(interp1(q,beta_raw,qq,'pchip','extrap'),0);
lambda = 0.45;

secondary = zeros(nrep,1);
rng(42077,'twister');
for r = 1:nrep
    infected = false(N,N);
    newly_by_index = false(N,N);
    c = ceil(N/2);
    infected(c,c) = true;
    for tau = 1:numel(beta_daily)
        neigh = moore_neighbors(c,c,N);
        for k = 1:size(neigh,1)
            ii = neigh(k,1); jj = neigh(k,2);
            if ~infected(ii,jj)
                p = 1 - exp(-lambda*beta_daily(tau));
                if rand < p
                    infected(ii,jj) = true;
                    newly_by_index(ii,jj) = true;
                end
            end
        end
    end
    secondary(r) = nnz(newly_by_index);
end

Rhat = mean(secondary);
sd = std(secondary);
ci = Rhat + [-1 1]*1.96*sd/sqrt(nrep);
T = table(string(mode), N, nrep, lambda, sum(beta_daily), Rhat, sd, ci(1), ci(2), ...
    'VariableNames', {'Mode','N','Replicates','lambda_inf','AUC_beta_daily','R0_CA_hat','SD','CI95_low','CI95_high'});
writetable(T, fullfile('export','phase4_ca_reproduction_estimate.csv'));

fprintf('\nDirect CA reproduction-scale estimate exported: export/phase4_ca_reproduction_estimate.csv\n');
fprintf('  R0_CA_hat = %.3f +/- %.3f, 95%% CI [%.3f, %.3f], n=%d\n', Rhat, sd, ci(1), ci(2), nrep);
end

function [t,Y] = solve7(P)
y0 = [P.S0; P.I0; P.Ag0; P.A17_0; P.V0; P.R0; P.Ab0];
opts = odeset('RelTol',1e-10,'AbsTol',1e-12,'MaxStep',0.1,'NonNegative',1:7);
[t,Y] = ode15s(@(t,y) measles_ode_rhs(t,y,P), [0 P.tmax_days], y0, opts);
end

function nb = moore_neighbors(i,j,N)
nb = zeros(8,2); k=0;
for di=-1:1
    for dj=-1:1
        if di==0 && dj==0, continue; end
        k=k+1;
        nb(k,:) = [mod(i+di-1,N)+1, mod(j+dj-1,N)+1];
    end
end
end
