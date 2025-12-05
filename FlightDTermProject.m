function [cg, cp, sigma_cg, sigma_cp] = centerlocation(mass, location, mass_tol, location_tol, geom_tol, alpha)
    % 1) First calculate the center of gravity 
    total_mass = sum(mass);
    moment = location .* mass;
    cg = sum(moment) / total_mass;

    % 2) Calculate the error propogation for the center of gravity
    % partial derivatives
    dcdm = (location * total_mass - sum(moment)) / (total_mass^2);
    dcdx = mass / total_mass;

    % variance from each term
    cg_variance = sum((dcdm(:) .* mass_tol(:)).^2) + sum((dcdx(:) .* location_tol(:)).^2);

    % standard deviation 
    sigma_cg = sqrt(cg_variance);

    % 3) Determine the center of pressure 
    %Given dimensions
    Ln = 1.2;
    d = 0.4;
    df = 0.4; 
    dr = 0.25;
    Lt = 0.8;
    Xp = 7.2;
    Cr = 0.8;
    Ct = 0.7;
    S = 0.2;
    Lf = 0.2;
    R = 0.2;
    Xr = 0.1;
    Xb = 6;
    N = 2;
    
    % Calculate the components of the Borrowman equation 
    CN_nose = 2;
    X_nose = 0.666*Ln;

    CN_trans = 2*((dr/d)^2 - (df/d)^2);
    X_trans = Xp + (Lt/3)*(1 + ((1-(df/dr))/(1-(df/dr)^2)));

    CN_fin = (1 + (R/(S+R)))*((4*N*(S/d)^2)/(1+sqrt(1+(2*Lf/(Cr+Ct))^2)));
    X_fin = Xb + (Xr/3)*(Cr + 2*Ct)/(Cr+Ct) + (1/6)*((Cr+Ct)-((Cr*Ct)/(Cr+Ct)));

    % Plug into the center of pressure function 
    cp = ((CN_nose*X_nose) + (CN_trans*X_trans) + (CN_fin*X_fin))/(CN_nose + CN_trans + CN_fin);

    % 4) Propogate the uncertainty of the center of pressure
    % Use numerical partial derivatives
    geom_tol = geom_tol*ones(13,1);
    delta = 1e-6;

    % List of all geometric variables
    vars = {Ln, d, df, dr, Lt, Xp, Cr, Ct, S, Lf, R, Xr, Xb};

    % Helper function: recompute CP with perturbed geometry
    function out = CP_eval(Ln_, d_, df_, dr_, Lt_, Xp_, Cr_, Ct_, S_, Lf_, R_, Xr_, Xb_)
        CN_t = 2*((dr_/d_)^2 - (df_/d_)^2);
        X_t = Xp_ + (Lt_/3)*(1 + ((1-(df_/dr_))/(1-(df_/dr_)^2)));
        CN_f = (1 + (R_/(S_+R_))) * ((4*N*(S_/d_)^2)/(1+sqrt(1+(2*Lf_/(Cr_+Ct_))^2)));
        X_f = Xb_ + (Xr_/3)*(Cr_ + 2*Ct_)/(Cr_ + Ct_) + ...
              (1/6)*((Cr_ + Ct_) - (Cr_*Ct_)/(Cr_ + Ct_));
        num = CN_nose*(0.666*Ln_) + CN_t*X_t + CN_f*X_f;
        den = CN_nose + CN_t + CN_f;
        out = num / den;
    end

    % Compute partial derivatives for each geometric variable
    nVars = numel(vars);
    dcp = zeros(nVars,1);

    for k = 1:nVars
        perturbed = vars;
        perturbed{k} = perturbed{k} + delta;

        dcp(k) = (CP_eval(perturbed{:}) - cp) / delta;
    end

    % Combine uncertainty: σ^2 = Σ (∂cp/∂xi * σ)^2
    cp_variance = sum((dcp .* geom_tol).^2 );
    sigma_cp = sqrt(cp_variance);

end 

function cg = center_of_gravity(mass, location, mass_tol, location_tol)
    total_mass = sum(mass);
    moment = location .* mass;
    cg = sum(moment) / total_mass;

    % uncertainty propagation
    % partial derivatives
    dcdm = (location * total_mass - sum(moment)) / (total_mass^2);
    dcdx = mass / total_mass;

    % variance from each term
    cg_variance = sum((dcdm(:) .* mass_tol(:)).^2) + sum((dcdx(:) .* location_tol(:)).^2);

    % standard deviation 
    cg_uncertainty = sqrt(cg_variance);
end

function cp = center_of_pressure(location_tol, a)
    Ln = 1.2;
    d = 0.4;
    df = 0.4; 
    dr = 0.25;
    Lt = 0.8;
    Xp = 7.2;
    Cr = 0.8;
    Ct = 0.7;
    S = 0.2;
    Lf = 0.2;
    R = 0.2;
    Xr = 0.1;
    Xb = 6;
    N = 2;

    CN_nose = 2;
    X_nose = 0.666*Ln;

    CN_trans = 2*((dr/d)^2 - (df/d)^2);
    X_trans = Xp + (Lt/3)*(1 + ((1-(df/dr))/(1-(df/dr)^2)));

    CN_fin = (1 + (R/(S+R)))*((4*N*(S/d)^2)/(1+sqrt(1+(2*Lf/(Cr+Ct))^2)));
    X_fin = Xb + (Xr/3)*(Cr + 2*Ct)/(Cr+Ct) + (1/6)*((Cr+Ct)-((Cr*Ct)/(Cr+Ct)));

    cp = ((CN_nose*X_nose) + (CN_trans*X_trans) + (CN_fin*X_fin))/(CN_nose + CN_trans + CN_fin);

    % All geometric distances have tolerance location_tol
    % Use numerical partial derivatives
    delta = 1e-6;

    % List of all geometric variables
    vars = {Ln, d, df, dr, Lt, Xp, Cr, Ct, S, Lf, R, Xr, Xb};

    % Helper function: recompute CP with perturbed geometry
    function out = CP_eval(Ln_, d_, df_, dr_, Lt_, Xp_, Cr_, Ct_, S_, Lf_, R_, Xr_, Xb_)
        CN_t = 2*((dr_/d_)^2 - (df_/d_)^2);
        X_t = Xp_ + (Lt_/3)*(1 + ((1-(df_/dr_))/(1-(df_/dr_)^2)));
        CN_f = (1 + (R_/(S_+R_))) * ((4*N*(S_/d_)^2)/(1+sqrt(1+(2*Lf_/(Cr_+Ct_))^2)));
        X_f = Xb_ + (Xr_/3)*(Cr_ + 2*Ct_)/(Cr_ + Ct_) + ...
              (1/6)*((Cr_ + Ct_) - (Cr_*Ct_)/(Cr_ + Ct_));
        num = CN_nose*(0.666*Ln_) + CN_t*X_t + CN_f*X_f;
        den = CN_nose + CN_t + CN_f;
        out = num / den;
    end

    % Compute partial derivatives for each geometric variable
    nVars = numel(vars);
    dcp = zeros(nVars,1);

    for k = 1:nVars
        perturbed = vars;
        perturbed{k} = perturbed{k} + delta;

        dcp(k) = (CP_eval(perturbed{:}) - cp) / delta;
    end

    % Combine uncertainty: σ^2 = Σ (∂cp/∂xi * σ)^2
    cp_variance = sum( (dcp * location_tol).^2 );
    cp_uncertainty = sqrt(cp_variance);

end
    


location = [0, 0.3, 0.7, 1.1, 1.2, 2.2, 2.3, 2.9, 3.1, 3.3, 3.4, 4, 4.2, 4.3, 4.9, 6, 6.5, 6.7, 6.9, 7.1, 7.6]
mass = [0, 20, 25, 25, 20, 25, 70, 1, 10, 15, 1, 20, 120, 1, 50, 5, 1, 5, 5, 10, 10]

mass_tol = 0.1 
location_tol = 0.002 
geom_tol = 0.002

total_mass = sum(mass);
moment = location .* mass;
cg = sum(moment) / total_mass

% uncertainty propagation
    % partial derivatives
dcdm = (location * total_mass - sum(moment)) / (total_mass^2);
dcdx = mass / total_mass;

    % variance from each term
cg_variance = sum((dcdm(:) .* mass_tol(:)).^2) + sum((dcdx(:) .* location_tol(:)).^2);

    % standard deviation 
cg_uncertainty = sqrt(cg_variance)

Ln = 1.2 
d = 0.4 
df = 0.4 
dr = 0.25
Lt = 0.8
Xp = 7.2
Cr = 0.8 
Ct = 0.7
S = 0.2 
Lf = 0.2
R = 0.2
Xr = 0.1
Xb = 6
N = 2
