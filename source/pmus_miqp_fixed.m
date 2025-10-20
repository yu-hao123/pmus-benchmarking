function [waveforms_true, waveforms_hat, params_true, params_fixed, solver_time, switching_times] = ...
    pmus_miqp_fixed(waveforms, initial_delay, l2_reg, delay_length, epsilon, R_fixed, E_fixed)
% PMUS_MIQP_FIXED Estimates the Pmus waveform for fixed R and E values
% using a mixed-integer quadratic programming formulation.
%
% Input:
%   waveforms   - struct with fields: flow, volume, paw, pmus, insexp
%   initial_delay
%   l2_reg
%   delay_length - integer, delay before respiratory effort begins
%   epsilon     - small number for inequality margin
%   R_fixed     - resistance [(cmH2O.s)/mL]
%   E_fixed     - elastance [cmH2O/mL]
%
% Output:
%   waveforms_true - ground truth signals
%   waveforms_hat  - estimated signals
%   params_true    - ground truth parameters (from LSE on true pmus)
%   params_fixed   - fixed R and E values used
%   solver_time    - solver runtime
%   switching_times - estimated switching instants

%% Default parameters
if nargin < 2, initial_delay = false; end
if nargin < 3, l2_reg = false; end
if nargin < 4, delay_length = 20; end
if nargin < 5, epsilon = 1e-3; end

%% Extract waveforms
flow = waveforms.flow;
volume = waveforms.volume;
paw = waveforms.paw;
pmus_true = waveforms.pmus;
insexp = waveforms.insexp;

k_soe = find(diff(insexp) <= -0.5, 1) + 1;

padding = zeros(delay_length,1);
if initial_delay
    flow   = [padding; flow];
    volume = [padding; volume];
    paw    = [padding; paw];
    pmus_true = [padding; pmus_true];
    insexp = [padding; insexp];
    k_soe = k_soe + delay_length;
end

flow = flow / 60 * 1000; % (L/min -> mL/s)

N = length(flow);
if initial_delay
    Ns = 3;
else
    Ns = 2;
end
tik = binvar(N, Ns);

constraint_occur = [];
for i=2:Ns
    constraint_occur = [constraint_occur, ...
        (1:N)*tik(:,i-1) <= (1:N)*tik(:,i)];
end

constraint_unique = [];
for i=1:Ns
    constraint_unique = [constraint_unique, sum(tik(:,i)) == 1];
end

%% regions
pmus = sdpvar(N,1);
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

%% cost function (R and E fixed!)
paw_hat = pmus + R_fixed * flow + E_fixed * volume;
residual = paw - paw_hat;
cost = residual' * residual;
if l2_reg
    cost = cost + 1.0e-3 * (pmus' * pmus);
end

constraint_real = [ones(N,1)*(-20) <= pmus <= ones(N,1)*1];
options = sdpsettings('solver','gurobi', ...
                      'gurobi.TimeLimit', 60, ...
                      'gurobi.TuneTimeLimit', 0, ...
                      'gurobi.OutputFlag',0, ...
                      'gurobi.threads',1,...
                      'verbose',0);

solution = optimize([constraint_occur, constraint_unique, ...
                     constraint_regions, constraint_real, ...
                     constraint_exhalation], cost, options);

if solution.problem ~= 0
    fprintf('Solver failed for status=%d, R = %.4f; C = %.4f\n', solution.problem, R_fixed * 1000, 1/E_fixed);
end

solver_time = solution.solvertime;

waveforms_true.flow   = flow;
waveforms_true.volume = volume;
waveforms_true.paw    = paw;
waveforms_true.pmus   = pmus_true;
waveforms_true.insexp = insexp;

waveforms_hat = table();
waveforms_hat.pmus = value(pmus);
waveforms_hat.paw  = waveforms_hat.pmus + R_fixed*flow + E_fixed*volume;

% true parameters from LSE on pmus_true
params_lse = (([flow volume]'*[flow volume]) \ ([flow volume]')) * (paw - pmus_true);
params_true.resistance = params_lse(1)*1000;
params_true.elastance  = params_lse(2)*1000;

params_fixed.resistance = R_fixed*1000;
params_fixed.elastance  = E_fixed*1000;

switching_times = zeros(Ns,1);
for i=1:Ns
    switching_times(i) = (1:N)*value(tik(:,i));
end

end
