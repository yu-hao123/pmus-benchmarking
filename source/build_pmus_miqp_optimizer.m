function miqp_optimizer = build_pmus_miqp_optimizer(waveforms, initial_delay, l2_reg, delay_length, epsilon)
% BUILD_PMUS_MIQP_OPTIMIZER
% Builds a parametric MIQP optimizer object for estimating Pmus.
% R and E are treated as parameters.

%% Default parameters
if nargin < 2, initial_delay = false; end
if nargin < 3, l2_reg = false; end
if nargin < 4, delay_length = 20; end
if nargin < 5, epsilon = 1e-3; end

flow = waveforms.flow;
volume = waveforms.volume;
paw = waveforms.paw;
pmus_true = waveforms.pmus;
insexp = waveforms.insexp;

% Start of exhalation
k_soe = find(diff(insexp) <= -0.5, 1) + 1;
padding = zeros(delay_length,1);
if initial_delay
    flow = [padding; flow];
    volume = [padding; volume];
    paw = [padding; paw];
    pmus_true = [padding; pmus_true];
    insexp = [padding; insexp];
    k_soe = k_soe + delay_length;
end

% Convert units (L/min → mL/s)
flow = flow / 60 * 1000;
N = length(flow);

if initial_delay
    Ns = 3;
else
    Ns = 2;
end

%% Decision variables
tik = binvar(N, Ns);
pmus = sdpvar(N,1);

%% Parameters (R and E will be supplied later)
Rpar = sdpvar(1,1);
Epar = sdpvar(1,1);

%% Switching constraints
constraint_occur = [];
for i=2:Ns
    constraint_occur = [constraint_occur, (1:N)*tik(:,i-1) <= (1:N)*tik(:,i)];
end
constraint_unique = [];
for i=1:Ns
    constraint_unique = [constraint_unique, sum(tik(:,i)) == 1];
end

%% Regions and pmus dynamics
constraint_regions = [];
if initial_delay
    for i=1:N-1
        constraint_regions = [constraint_regions, ...
            implies(sum(tik(1:i,1)) <= 0.5, [pmus(i)==0, pmus(i+1)==0]), ...
            implies(0.5 <= (sum(tik(1:i,1))-sum(tik(1:i,2))), pmus(i+1)+epsilon <= pmus(i)), ...
            implies(0.5 <= (sum(tik(1:i,2))-sum(tik(1:i,3))), pmus(i)+epsilon <= pmus(i+1)), ...
            implies(0.5 <= sum(tik(1:i,3)), [pmus(i)==0, pmus(i+1)==0]) ];
    end
else
    for i=1:N-1
        constraint_regions = [constraint_regions, ...
            implies(sum(tik(1:i,1)) <= 0.5, pmus(i+1)+epsilon <= pmus(i)), ...
            implies(0.5 <= (sum(tik(1:i,1))-sum(tik(1:i,2))), pmus(i)+epsilon <= pmus(i+1)), ...
            implies(0.5 <= sum(tik(1:i,2)), pmus(i)==pmus(i+1)) ];
    end
end

%% Exhalation switching constraint
tau_soe = 50;
constraint_exhalation = (1:N)*tik(:,end) <= k_soe + tau_soe;

%% Cost function (R,E are parameters!)
paw_hat = pmus + Rpar * flow + Epar * volume;
residual = paw - paw_hat;
cost = residual' * residual;
if l2_reg
    cost = cost + 1.0e-3 * (pmus' * pmus);
end

%% pmus bounds
constraint_real = [ones(N,1)*(-20) <= pmus <= ones(N,1)*1.0];
constraints = [constraint_occur, constraint_unique, ...
                constraint_regions, constraint_real, ...
                constraint_exhalation];

%% Solver options
options = sdpsettings('solver','gurobi', ...
    'gurobi.OutputFlag', 0, ...
    'gurobi.TimeLimit', 30, ...
    'gurobi.Threads', 8, ...
    'verbose', 0);

%% Build optimizer
% Input parameters: {Rpar, Epar}
% Outputs: {pmus, paw_hat}
miqp_optimizer = optimizer(constraints, cost, options, {Rpar, Epar}, {pmus, paw_hat});

end
