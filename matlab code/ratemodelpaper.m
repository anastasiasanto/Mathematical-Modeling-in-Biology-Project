
params_dict = struct(...
    'tau_I', 1, ...                  % timescale of neuroinflammatory reaction [days]
    'tau_B', 10, ...                 % timescale of BBB recovery [days]
    'tau_D', 10, ...                 % timescale of neuronal death process [days]
    'tau_R', 10, ...                 % timescale of circuit remodeling [days]
    'k_IB', 0.1, ...                 % scaling parameter for effect of neuroinflammation on BBB permeability [-]
    'k_BI', 1, ...                   % scaling parameter for proinflammatory effect of BBB leakage [-]
    'k_IS', 2, ...                   % scaling parameter for strength of seizure-promoting effects of neuroinflammation [-]
    'k_RS', 2, ...                   % scaling parameter for strength of seizure-promoting effects of circuit remodeling [-]
    'k_ID', 8, ...                   % scaling parameter for neurotoxic effect of overactivated glia [-]
    'k_BR', 1, ...                   % scaling parameter of BBB leakage on circuit remodeling [-]
    'k_DR', 0.0005, ...              % scaling parameter of neuronal loss on circuit remodeling [-]
    'K_SB', 0.875, ...                % scaling parameter for seizure burden on BBB integrity [-]
    'D_m', 1, ...                    % maximum possible extent of neuronal loss [-]
    'Theta', 0.25, ...                % Neurotoxicity threshold of overactivated glia [-]
    'IBDR_E_duration', [0, 2, 2, 0], ...    % Integration time step [days]
    'IBDR_E_amplitude', [0, 1.65, 1, 0], ... % Integration time step [days]
    'Complex_input', 'no', ...       % Flag for complex input
    'amount_simulations', 1, ...     % Amount of simulations
    'number_simulation', 1 ,...       % Simulation num;ber
    'IC',[0, 0, 0 ,0]...
 );

% Adding simulation-specific parameters
t_end = 500;  % Sostituire con il valore desiderato per t_end
dt = 0.1;  % Sostituire con il valore desiderato per dt
IC=[0 0 0 0]; % initial conditions
t_vec = linspace(0, t_end, int16(t_end/dt)+1);
% I_vec = zeros(size(t_vec));
% B_vec = zeros(size(t_vec));
% D_vec = zeros(size(t_vec));
% R_vec = zeros(size(t_vec));
% Risoluzione con ode45
[t, Y] = ode45(@(t, y) dIBDRdt_Rate(t, y, params_dict) , t_vec, IC);
I_vec = Y(:, 1);
B_vec = Y(:, 2);
D_vec = Y(:, 3);
R_vec = Y(:, 4);
% Creazione di un dizionario dei risultati
results_dict = struct('t_vec', t, 'I_vec', I_vec, 'B_vec', B_vec, 'D_vec', D_vec, 'R_vec', R_vec);

% Salvataggio dei risultati
filename = strcat('C:\Users\santo\OneDrive\Desktop\mathematical modeling\presentazione pugliese\POLIRATE');
save(filename, 'params_dict', 'results_dict');
clear ii filename;
function dxdt = dIBDRdt_Rate(t, y, params_dict)
    I0 = y(1); B0 = y(2); D0 = y(3);R0 = y(4);
    inp = zeros(size(params_dict.IC));

    for ww = 1:length(inp)
        if t <= params_dict.IBDR_E_duration(ww)
            inp(ww) = params_dict.IBDR_E_amplitude(ww);       
        end
    end

    S = SeizureBurden(I0, R0, params_dict);
    dIdt = (1 / (params_dict.tau_I)) * (-I0 + params_dict.k_BI*B0 + inp(1));
    dBdt= (1 / (params_dict.tau_B)) * (-B0 + params_dict.k_IB*I0 + S + inp(2));
    dDdt = (1 / (params_dict.tau_D)) * ((1 - D0 / params_dict.D_m) * params_dict.k_ID * max(0, I0 - params_dict.Theta) + inp(3));
    dRdt = (1 / (params_dict.tau_R)) * (-R0 + params_dict.k_BR*B0 + params_dict.k_DR*D0 + inp(4));
    dxdt=[dIdt; dBdt; dDdt; dRdt];
    
    
end

function SB = SeizureBurden(I, R, par_dict)
    SB = (par_dict.K_SB * (exp(par_dict.k_IS*I^2 + par_dict.k_RS*R) - 1) )/ (exp(par_dict.k_IS*I^2 + par_dict.k_RS*R) + 1);
end

%function result = fun_dbdt(x, params_dict) 
% function db/dt, where D=D_const; R(t) is expressed as R=kappa_BR*B+kappa_DR*D_const;
% I is expressed as I=kappa_BI*B
    % result = -x + params_dict.k_IB*params_dict.k_BI*x + params_dict.K_SB * ...
    %          (exp(params_dict.k_IS*params_dict.k_BI*params_dict.k_BI*x*x + params_dict.k_RS*params_dict.k_BR*x + ...
    %          params_dict.k_RS*params_dict.k_DR*params_dict.D_const) - 1) / ...
    %          (exp(params_dict.k_IS*params_dict.k_BI*params_dict.k_BI*x*x + params_dict.k_RS*params_dict.k_BR*x + ...
    %          params_dict.k_RS*params_dict.k_DR*params_dict.D_const) + 1);
%end