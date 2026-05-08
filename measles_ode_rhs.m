function dydt = measles_ode_rhs(t, y, P)
% MEASLES_ODE_RHS  Right-hand side of the measles ODE system.
%   Auto-detects 4-variable Morris-core mode or 7-variable extended mode.
%
% Therapy interpretation:
%   - DC_mode='multiplicative' reproduces the manuscript architecture:
%       q_eff = q*(1 + eta_D*D(t)); rho_eff = rho_Ab*(1 + eta_AbD*D(t)).
%     This explains why DC-only stimulation has negligible effect when the
%     immunocompromised effector substrate is nearly absent.
%   - DC_mode='additive_control' is a diagnostic architecture control used
%     only in Phase 3. It adds absolute stimulation/recruitment terms and is
%     not propagated to the CA as the main therapy profile.

n = numel(y);

if n == 4
    S = max(y(1),0); I = max(y(2),0); A = max(y(3),0); V = max(y(4),0);
    fV = V/(P.s + V + eps);
    dh = double(t < P.td);

    dS = -P.beta*S*V + P.qs*dh*S + P.r*(1-fV)*A;
    dI =  P.beta*(S+A)*V - P.delta*I - P.k*I*A;
    dA =  P.q*fV*A - P.beta*A*V - (1-fV)*(P.d+P.r)*A;
    dV =  P.p*I - P.c*V;
    dydt = [dS; dI; dA; dV];

elseif n == 7
    S=max(y(1),0); I=max(y(2),0); Ag=max(y(3),0); A17=max(y(4),0);
    V=max(y(5),0); R=max(y(6),0); Ab=max(y(7),0);
    Atot = Ag + A17;
    fV = V/(P.s + V + eps);
    dh = double(t < P.td);

    % Shift gate; clipping prevents overflow without changing biological scale.
    z = P.alpha_t*(t - P.t17) + P.alpha_R*log(1+R);
    z = min(max(z, -60), 60);
    sig17 = 1/(1+exp(-z));

    q_eff = P.q;
    rho_eff = P.rho_Ab;
    u_DC_Ag = 0;

    if isfield(P,'DC_active') && P.DC_active
        Dt = dc_input(t, P);
        if ~isfield(P,'DC_mode') || strcmpi(P.DC_mode, 'multiplicative')
            q_eff   = P.q      * (1 + P.eta_D   * Dt);
            rho_eff = P.rho_Ab * (1 + P.eta_AbD * Dt);
        elseif strcmpi(P.DC_mode, 'additive_control')
            q_eff   = P.q      + P.DC_add_q     * Dt;
            rho_eff = P.rho_Ab + P.DC_add_rhoAb * Dt;
            u_DC_Ag = P.DC_add_Ag * Dt;
        else
            error('Unknown DC_mode: %s', P.DC_mode);
        end
    end

    u_CTL = 0;
    if isfield(P,'CTL_active') && P.CTL_active
        u_CTL = (P.CTL_delta/(P.CTL_tau*sqrt(2*pi))) * ...
                exp(-0.5*((t-P.CTL_time)/P.CTL_tau)^2);
    end
    H_Ab = double(t >= P.t_Ab);

    dS   = -P.beta*S*V + P.qs*dh*S + P.r*(1-fV)*Atot;
    dI   =  P.beta*(S+Atot)*V - P.delta*I - P.k*I*Atot;
    dAg  =  q_eff*fV*Atot*(1-sig17) - P.beta*Ag*V ...
            -(1-fV)*(P.d+P.r)*Ag - P.kappa17*Ag + u_CTL + u_DC_Ag;
    dA17 =  q_eff*fV*Atot*sig17 + P.kappa17*Ag ...
            -P.beta*A17*V - (1-fV)*(P.d+P.r)*A17;
    dV   =  P.p*I - P.c*V - P.theta_Ab*Ab*V;
    dR   =  P.rho_R*I - P.mu_R*R - P.eta_R*Ab*R;
    gR   =  antibody_stimulation_function(R, P);
    dAb  =  rho_eff*H_Ab*gR - P.mu_Ab*Ab;
    dydt = [dS; dI; dAg; dA17; dV; dR; dAb];
else
    error('Expected 4 or 7 variables, got %d.', n);
end
end

function Dt = dc_input(t, P)
Dt = 0;
for n = 0:(P.DC_doses-1)
    tn = P.DC_start + n*P.DC_interval;
    if t >= tn
        Dt = Dt + P.DC_ef*P.DC_d0*exp(-P.DC_muD*(t-tn));
    end
end
end
