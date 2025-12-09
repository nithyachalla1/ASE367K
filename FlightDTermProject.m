function [cg, cg_uncertainty] = center_of_gravity(mass, location, mass_tol, location_tol, a)
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

function [cp, cp_uncertainty] = center_of_pressure(location_tol, a)
   % Since there is no transition in the main body of the rocket, Lt, Xp,
   % df, and dr
   % are excluded from the data set
    Ln = 1.2;      
    d = 0.4;
    Xb = 6.0;       
    Xr = 0.0;       
    Cr = 0.8;       
    Ct = 0.7;      
    S  = 0.2;       
    Lf = 0.2;       
    R  = 0.2;       
    N  = 2;    

    % First determine the coefficient of the normal force and location for
    % the nose
    CN_nose = 2;
    X_nose = 0.666*Ln;

    % No transition normal force coefficient 

    % Determine the coefficient of the normal force and location for the
    % fins
    CN_fin = (1 + (R/(S+R)))*((4*N*(S/d)^2)/(1+sqrt(1+(2*Lf/(Cr+Ct))^2)));
    X_fin = Xb + (Xr/3)*(Cr + 2*Ct)/(Cr+Ct) + (1/6)*((Cr+Ct)-((Cr*Ct)/(Cr+Ct)));

    % Calculate cp using a weighted sum
    cp = ((CN_nose*X_nose) + (CN_fin*X_fin))/(CN_nose + CN_fin);

    % Uncertainty propogation
    vars = {Ln, d, Cr, Ct, S, Lf, R, Xr, Xb};
    delta = 1e-6;
    nVars = numel(vars);
    dcp = zeros(nVars,1);

    function out = CP_eval(Ln_, d_, Cr_, Ct_, S_, Lf_, R_, Xr_, Xb_)
        CNn = 2;
        Xn_ = 0.666 * Ln_;
        CNf_ = (1 + R_/(S_+R_)) * ...
               ( (4*N*(S_/d_)^2) / (1 + sqrt(1 + (2*Lf_/(Cr_+Ct_))^2)) );
        Xf_ = Xb_ + (Xr_/3) * (Cr_ + 2*Ct_)/(Cr_ + Ct_) + ...
               (1/6) * ((Cr_ + Ct_) - (Cr_*Ct_)/(Cr_ + Ct_));
        out = (CNn*Xn_ + CNf_*Xf_) / (CNn + CNf_);
    end

    base_cp = cp;
    for k = 1:nVars
        pert = vars;
        pert{k} = pert{k} + delta;
        dcp(k) = ( CP_eval(pert{:}) - base_cp ) / delta;
    end

    % Combine (σ_cp)^2 = sum( (∂cp/∂xi * σ_xi)^2 ), σ_xi = location_tol
    cp_variance = sum( (dcp * location_tol).^2 );
    cp_uncertainty = sqrt(cp_variance);

end

% Mass and location distribution
mass = [0, 20, 25, 25, 20, 25, 70, 1, 10, 15, 1, 20, 120, 1, 50, 5, 1, 5, 5, 10, 10];
location = [0, 0.3, 0.7, 1.1, 1.2, 2.2, 2.3, 2.9, 3.1, 3.3, 3.4, 4, 4.2, 4.3, 4.9, 6, 6.5, 6.7, 6.9, 7.1, 7.6];

% Tolerances 
mass_tol = 0.1;
location_tol = 0.002;

% Print cg and cp values alongside their uncertainty for aoa = 0
[cg, cg_uncertainty] = center_of_gravity(mass, location, mass_tol, location_tol, 0);
[cp, cp_uncertainty] = center_of_pressure(location_tol, 0);

fprintf('Center of Gravity (CG): %.4f m\n', cg);
fprintf('CG Uncertainty: ±%.4f m \n\n', cg_uncertainty);

fprintf('Center of Pressure (CP): %.4f m\n', cp);
fprintf('CP Uncertainty: ±%.4f m \n\n', cp_uncertainty);

alphas = linspace(0, 5, 50);   % AoA range (deg)

cg_vals = zeros(size(alphas));
cg_err  = zeros(size(alphas));

cp_vals = zeros(size(alphas));
cp_err  = zeros(size(alphas));

% For each value of alpha, graph the center of gravity and center of
% pressure
for i = 1:length(alphas)
    [cg_vals(i), cg_err(i)] = center_of_gravity(mass, location, mass_tol, location_tol, alphas(i));
    [cp_vals(i), cp_err(i)] = center_of_pressure(location_tol, alphas(i));
end

figure;

subplot(2,1,1)
errorbar(alphas, cg_vals, cg_err, 'LineWidth', 1.5);
ylabel('CG (m)');
title('Center of Gravity vs AoA');
grid on;

subplot(2,1,2)
errorbar(alphas, cp_vals, cp_err, 'LineWidth', 1.5);
xlabel('Angle of Attack (deg)');
ylabel('CP (m)');
title('Center of Pressure vs AoA');
grid on;


