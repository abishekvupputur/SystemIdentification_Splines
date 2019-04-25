%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Demonstration of the Linear Kalman Filter
%
%   Author: C.C. de Visser, Delft University of Technology, 2013
%   email: c.c.devisser@tudelft.nl
%   Version: 1.0
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
close all;
clear all;

% Use these commands to initialize the randomizer to a fixed (reproducable) state.
% rng('default'); % init randomizer (default, fixed)-> version 2014a,b
% RandStream.setDefaultStream(RandStream('mt19937ar','seed', 300));-> version 2013a,b

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Set simulation parameters
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
dt              = 0.01;
N               = 10000;
epsilon         = 1e-10;
doIEKF          = 1;
maxIterations   = 150;


printfigs = 0;
figpath = '';
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Set initial values for states and statistics
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Ex_0    = [10; 1; 0 ; 0.01]; % initial estimate of optimal value of x_k_1k_1
x_0     = [0; 0.45; 0 ; 0.01]; % initial state
m       = 3; % number of input dimensions

% Initial estimate for covariance matrix
stdx_0  = [10 10 10 10000];
P_0     = diag(stdx_0.^2);

% System noise statistics:
Ew = [0 0 0 0]; % bias
stdw = [1e-3 1e-3 1e-3 0]; % noise variance
Q = diag(stdw.^2);
n  = length(stdw);
w_k = diag(stdw) * randn(n, N)  + diag(Ew) * ones(n, N);

% Measurement noise statistics:
Ev = [0 0 0]; % bias
stdv = [0.035 0.013 0.11]; % noise variance
R = diag(stdv.^2);
nm = length(stdv);
v_k = diag(stdv) * randn(nm, N)  + diag(Ev) * ones(nm, N);

G = eye(n); % noise input matrix
B = [eye(m);0,0,0]; % input matrix


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Calculate batch with measurement data
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

tic;

% Real simulated state-variable and measurements data:
x = x_0;
X_k = zeros(n, N);
%Z_k = zeros(nm, N);
%U_k = zeros(m, N);
dataname = 'F16traindata_CMabV_2019';
load(dataname, 'Cm', 'Z_k', 'U_k')
Z_k = Z_k';
U_k = U_k';
XX_k1k1 = zeros(n, N);
PP_k1k1 = zeros(n, N);
STDx_cor = zeros(n, N);
z_pred = zeros(nm, N);
IEKFitcount = zeros(N, 1);

x_k_1k_1 = Ex_0; % x(0|0)=E{x_0}
P_k_1k_1 = P_0; % P(0|0)=P(0)

time1 = toc;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Run the Extended Kalman filter
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
tic;
% Extended Kalman Filter (EKF)
ti = 0; 
tf = dt;

% Run the filter through all N samples
for k = 1:N
    % Prediction x(k+1|k) 
    [t, x_kk_1] = rk4(@kf_calc_f, x_k_1k_1,U_k(:,k), [ti tf]); 

    % z(k+1|k) (predicted output)
    z_kk_1 = kf_calc_h(0, x_kk_1, U_k(:,k)); %x_kk_1.^3; 
    z_pred(:,k) = z_kk_1;

    % Calc Phi(k+1,k) and Gamma(k+1, k)
    Fx = zeros(4); % perturbation of f(x,u,t)
    % the continuous to discrete time transformation of Df(x,u,t) and G
    [dummy, Psi] = c2d(Fx, B, dt);   
    [Phi, Gamma] = c2d(Fx, G, dt);   
    
    % P(k+1|k) (prediction covariance matrix)
    P_kk_1 = Phi*P_k_1k_1*Phi' + Gamma*Q*Gamma'; 
    P_pred = diag(P_kk_1);
    stdx_pred = sqrt(diag(P_kk_1));

    
    % Run the Iterated Extended Kalman filter (if doIEKF = 1), else run standard EKF
    if (doIEKF)

        % do the iterative part
        eta2    = x_kk_1;
        err     = 2*epsilon;

        itts    = 0;
        while (err > epsilon)
            if (itts >= maxIterations)
                fprintf('Terminating IEKF: exceeded max iterations (%d)\n', maxIterations);
                break
            end
            itts    = itts + 1;
            eta1    = eta2;

            % Construct the Jacobian H = d/dx(h(x))) with h(x) the observation model transition matrix 
            Hx       = kf_calc_Hx(0, eta1, U_k(:,k)); 
            
            % Check observability of state
            if (k == 1 && itts == 1)
                rankHF = kf_calcObsRank(Hx, Fx);
                if (rankHF < n)
                    warning('The current state is not observable; rank of Observability Matrix is %d, should be %d', rankHF, n);
                end
            end
            
            % The innovation matrix
            Ve  = (Hx*P_kk_1*Hx' + R);

            % calculate the Kalman gain matrix
            K       = P_kk_1 * Hx' / Ve;
            % new observation state
            z_p     = kf_calc_h(0, eta1, U_k(:,k)) ;%fpr_calcYm(eta1, u);

            eta2    = x_kk_1 + K * (Z_k(:,k) - z_p' - Hx*(x_kk_1 - eta1));
            err     = norm((eta2 - eta1), inf) / norm(eta1, inf);
        end

        IEKFitcount(k)    = itts;
        x_k_1k_1          = eta2;

    else
        % Correction
        Hx = kf_calc_Hx(0, x_kk_1, U_k(:,k)); % perturbation of h(x,u,t)
        % Pz(k+1|k) (covariance matrix of innovation)
        Ve = (Hx*P_kk_1 * Hx' + R); 

        % K(k+1) (gain)
        K = P_kk_1 * Hx' / Ve;
        % Calculate optimal state x(k+1|k+1) 
        x_k_1k_1 = x_kk_1 + K * (Z_k(:,k) - z_kk_1); 

    end    
    
    P_k_1k_1 = (eye(n) - K*Hx) * P_kk_1 * (eye(n) - K*Hx)' + K*R*K';  
    P_cor = diag(P_k_1k_1);
    stdx_cor = sqrt(diag(P_k_1k_1));

    % Next step
    ti = tf; 
    tf = tf + dt;
    
    % store results
    XX_k1k1(:,k) = x_k_1k_1;
%     PP_k1k1(k,:) = P_k_1k_1;
    STDx_cor(:,k) = stdx_cor;
end

time2 = toc;
Z_k1k1=zeros(nm,N+1);
dump=XX_k1k1(4,:);
XX_k1k1(4,:)=0;
%% Estimate Output parameters
for i=2:N+1
   Z_k1k1(1:3,i)=kf_calc_h(0, XX_k1k1(:,i-1), zeros(4,1));
end
XX_k1k1(4,:)=dump;
errrr=Z_k-Z_k1k1;

