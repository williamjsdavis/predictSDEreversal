function [t,X,eta] = OUcorr(dt,dt_int_max,N,f,g,theta,eta0,X0)
%Stochastic process with correlated noise, Euler-Maruyama method
%   William Davis, 10/01/19
%
%   Notes:
%   Realize a stochastic process driven by Ornstein-Uhlenbeck noise using 
%   the Euler-Maruyama method. Similar to integrating the Langevin 
%   equation:
%
%   X_{n+1} = X_n + f(X_n)\Delta t + g(X_n)\Delta t \eta_n,
%   \eta_{n+1} = \eta_n - \Delta t/\theta +1/\theta \Delta\xi_n,
%
%   where \Delta\xi_n are iid Weiner increments. These are normal variates
%   with zero mean and variance \Delta t. Thus W_t_{n+1}-W_t_n = \Delta W_n
%   which is is distributed N(-,\Delta t)=\sqrt{\delta t}N(0,1).
%   
%
%   Inputs:
%   - "dt"                      Time-step, dimensional
%   - "dt_int_max"              Maximum fraction of time-step for interior
%                               iteration
%   - "N"                       Number of time-steps
%   - "f"                       Drift function
%   - "g"                       Noise function
%   - "theta"                   Correlation time
%
%   Problems:
%   - 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Processing
t = 0:dt:(dt*(N-1)); % Time

% Time variables
min_t = min(dt,theta); % Minimum time
dt_int_aim = dt_int_max*min_t; % Aim for interior time-step
N_int = ceil(dt/dt_int_aim); % Number of time-steps per interior iteration
dt_int = dt/N_int; % True interior time-step

% Random noise (xi) for full internal scheme
mu = zeros(1,N_int);
stdev = ones(1,N_int);

% Preallocating sizes
eta = zeros(size(t)); % Stochastic force
X = zeros(size(t)); % Observed variable (functional eta)
X_int = zeros(1,N_int+1); % Internal vector (functional eta)

% Time checking
tt = zeros(size(t)); % External time
tt_int = zeros(1,N_int+1); % Full internal time

% Initial conditions
eta(1) = eta0; % Stochastic force
X(1) = X0; % Observed variable
tt(1) = 0; % Time

% Functional eta
t_int_in = 0:dt_int:dt; % Time to integrate over
ex = exp(-t_int_in/theta); % Exponential
expl = exp(t_int_in/theta); % Positive exponential

for ii = 2:N
    
    xi = random('norm',mu,stdev); % Create noise (for each internal loop)
    
    % Starting internal values
    X_int(1) = X(ii-1);
    tt_int(1) = tt(ii-1);
    
    % Functional eta
    eta_int_fun = eta(ii-1)*ex + sqrt(dt_int)*ex/theta.*cumsum(expl.*[0,xi]);
    
    % Interior loop
    for jj = 2:N_int+1
        % Observed variable with functional stochastic force
        X_int(jj) = X_int(jj-1) + ...
            f(X_int(jj-1))*dt_int + g(X_int(jj-1))*eta_int_fun(jj)*dt_int;
    end
    
    % Setting external values
    eta(ii) = eta_int_fun(end);
    X(ii) = X_int(jj);
end

end
