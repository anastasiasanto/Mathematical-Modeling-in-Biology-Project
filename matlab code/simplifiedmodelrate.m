%SIMPLIFICATION OF THE MODEL :
%Assumptions : 
%the process of neuroinflammatory reaction evolves faster than blood-brain barrier (BBB)
%disruption and recovery of its integrity, neuronal loss, and circuit remodeling: tItB,tD,tR. Under the condition of an absent neuroinflammatory external input IE, we can thus perform a time scale separation. At
%equilibrium, this will result in I=kbiB
%we can obtain a system of equations for fixed values of neuronal loss extent D = Dconst,
%where 0 <=Dconst <= Dmax. The resulting system of equations describes the system in the dynamical regimes
%characterized by the absence of neurotoxicity I<theta

%we obtain the system described in B-R dimensions. It is used for 
% analysis and visualization of dynamics with state-space plots for 
% variables B and R:

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
    'k_DR', 0.0001, ...              % scaling parameter of neuronal loss on circuit remodeling [-]
    'K_SB', 0.875, ...                % scaling parameter for seizure burden on BBB integrity [-]
    'D_m', 1, ...                    % maximum possible extent of neuronal loss [-]
    'Theta', 0.25, ...                % Neurotoxicity threshold of overactivated glia [-]
    'IBDR_E_duration', [0, 2, 2, 0], ...    % Integration time step [days]
    'IBDR_E_amplitude', [0, 1.65, 1, 0], ... % Integration time step [days]
    'Complex_input', 'no', ...       % Flag for complex input
    'amount_simulations', 1, ...     % Amount of simulations
    'number_simulation', 1 ,...       % Simulation num;ber
    'IC',[0, 0, 0 ,0],...
    'Dconst',0.41...
 );

% Original color
originalColor =  [.5 .5 .5];
% Make it more transparent (e.g., set alpha to 0.5 for half transparency)
transparentColor = [originalColor(1:3), 0.5];
%where I = kB/IB<Q and D = Dconst
t_end = 500;  % Sostituire con il valore desiderato per t_end
dt = 0.1;  % Sostituire con il valore desiderato per dt
IC=[-0.1 -0.1 ]; % initial conditions
t_vec = linspace(0, t_end, int16(t_end/dt)+1);
% I_vec = zeros(size(t_vec));
% B_vec = zeros(size(t_vec));
% D_vec = zeros(size(t_vec));
% R_vec = zeros(size(t_vec));
% Risoluzione con ode45
[t, Y] = ode45(@(t, y) Simplified_Rate(t, y, params_dict) , t_vec, IC);
Bs_vec = Y(:, 1);
Rs_vec = Y(:, 2);
% Creazione di un dizionario dei risultati
results_dictsimp = struct('t_vec', t, 'Bs_vec', Bs_vec, 'Rs_vec', Rs_vec);
figure(1);
plot(Bs_vec,Rs_vec, 'Color', "#D95319", "LineWidth",1) 
xlabel('Extent of BBB distruption')
ylabel('Degree of circuit Remodelling')
title('B vs R');
grid on;
% Definisci le equazioni differenziali del sistema
dVdt=@(t,B0) (1 / (params_dict.tau_B)) * (-B0 + params_dict.k_IB*params_dict.k_BI*B0 +(params_dict.K_SB * (exp(params_dict.k_IS*I^2 + params_dict.k_RS*R0) - 1) )/ (exp(params_dict.k_IS*I^2 + params_dict.k_RS*R0) + 1))
dndt=@(t,R0) (1 / (params_dict.tau_R)) * (-R0 + params_dict.k_BR*B0 + params_dict.k_DR*params_dict.Dconst) 
V_values=linspace(-0.1,1.1,20); %values for V
n_values=linspace(-0.1,1.1,20);% values for n
[V_grid,n_grid]=meshgrid(V_values,n_values); %grid of (V,n) values
dV_dt=zeros(size(V_grid));
dn_dt=zeros(size(n_grid));
%evaluate derivatives at each point in the grid
for i =1:numel(V_grid)
    x=[V_grid(i);n_grid(i)];
    dx=Simpli(0,x,params_dict);
    dV_dt(i)=dx(1);
    dn_dt(i)=dx(2);
end
%plot isoclines of simplified model
figure(3);
hold on; 
quiver(V_grid, n_grid, dV_dt, dn_dt, 'AutoScaleFactor',1.9, 'Color',transparentColor);
contour(V_grid, n_grid, dV_dt,[2 0], "Color",[0.8,0.13, 0.6], "LineWidth",1.5); 
contour(V_grid, n_grid, dn_dt,[2 0], "Color",[0.2,0.13, 0.8], "LineWidth",1.5);
plot([0.41 0.41],ylim, 'r--', 'LineWidth', 1.5, 'DisplayName','Neurotoxicity threshold');  
%title('Isoclines and Direction Field in the (B, R): K_S_B = 0.875');
title('Isoclines and Direction Field in the (B, R): D_cost = 0.41 ');
xlabel('Extent of BBB distruption');
ylabel('Degree of circuit Remodelling');
legend('dB/dt Isoclines','dR/dt Isoclines' ,'Direction Field');
grid on;  % k-- rappresenta una linea tratteggiata nera
% Calcola i punti d'intersezione tra le due rette
% Definisci le funzioni delle isocline
f1 = @(B0, R0) (1 / (params_dict.tau_B)) * (-B0 + params_dict.k_IB*params_dict.k_BI*B0 +(params_dict.K_SB * (exp(params_dict.k_IS*(params_dict.k_BI*B0)^2 + params_dict.k_RS*R0) - 1) )/ (exp(params_dict.k_IS*(params_dict.k_BI*B0)^2 + params_dict.k_RS*R0) + 1))%... inserisci qui l'espressione della prima isocline;
f2 = @(B0, R0) (1 / (params_dict.tau_R)) * (-R0 + params_dict.k_BR*B0 + params_dict.k_DR*params_dict.Dconst) %... inserisci qui l'espressione della seconda isocline;

% Supponiamo che tu voglia trovare i punti d'intersezione in una certa regione
% Definisci una funzione contenente il sistema di equazioni
equations = @(x) [f1(x(1), x(2)); f2(x(1), x(2))];

% Specifica un punto di partenza per la ricerca della soluzione
initial_guess = [0.3, 0.4]; % Sostituisci con i valori iniziali appropriati

% Risolvi il sistema di equazioni
intersection_point = fsolve(equations, initial_guess);

% Estrai i valori di V e n dall'output
B_intersection = intersection_point(1);
R_intersection = intersection_point(2);

% Ora V_intersection e n_intersection contengono le coordinate del punto d'intersezione
disp(['Punto d''intersezione: B = ', num2str(B_intersection), ', R = ', num2str(R_intersection)]);


% % Sostituisci con i tuoi calcoli effettivi
% intersection_points_x = [value1, value2, value3]; % Sostituisci con i valori effettivi sull'asse x
% intersection_points_y = [value4, value5, value6]; % Sostituisci con i valori effettivi sull'asse y
% % Traccia i cerchi nei punti d'intersezione
% for i = 1:length(intersection_points_x)
%     plot(intersection_points_x(i), intersection_points_y(i), 'ro', 'MarkerSize', 10, 'LineWidth', 2);
% end


function dxdt = Simplified_Rate(t, y, params_dict)
    B0 = y(1);R0 = y(2);
    S = SeizureBurden((params_dict.k_BI*B0), R0, params_dict);
    dBdt= (1 / (params_dict.tau_B)) * (-B0 + params_dict.k_IB*params_dict.k_BI*B0 + S);
    dRdt = (1 / (params_dict.tau_R)) * (-R0 + params_dict.k_BR*B0 + params_dict.k_DR*params_dict.Dconst);
    dxdt=[ dBdt; dRdt];
end
function SB = SeizureBurden(I, R, par_dict)
    SB = (par_dict.K_SB * (exp(par_dict.k_IS*I^2 + par_dict.k_RS*R) - 1) )/ (exp(par_dict.k_IS*I^2 + par_dict.k_RS*R) + 1);
end

function dxdt = Simpli(t, y, params_dict)
    B0 = y(1);R0 = y(2);
    dBdt= (1 / (params_dict.tau_B)) * (-B0 + params_dict.k_IB*params_dict.k_BI*B0 +(params_dict.K_SB * (exp(params_dict.k_IS*(params_dict.k_BI*B0)^2 + params_dict.k_RS*R0) - 1) )/ (exp(params_dict.k_IS*(params_dict.k_BI*B0)^2 + params_dict.k_RS*R0) + 1));
    dRdt = (1 / (params_dict.tau_R)) * (-R0 + params_dict.k_BR*B0 + params_dict.k_DR*params_dict.Dconst);  
    dxdt=[ dBdt; dRdt];
end
%prova D=0, kdr=0.001
%prova D=0.35, kdr=0.001 fig S9