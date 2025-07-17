function [waveforms_true, waveforms_hat, params_true, params_hat,...
    solver_time, switching_times] = pmus_miqp(waveforms, initial_delay,...
    l2_reg, delay_length, epsilon)
% PMUS_MIQP Estimates the Pmus waveform and the respiratory dynamics' 
% parameters (resistance and elastance) by employing a mixed-integer 
% quadratic programming formulation.
% Input:
%   filename - name of the binary file from the simulator.
%   initial_delay - if an initial delay should be considered before the 
%                   respiratory effort begins. Defaults to false.
%   l2_reg - if L2 regularization should be used. Defaluts to false.
%   delay_length - the length of the delay before the respiratory effort
%                  begins.
% Output:
%   waveforms_true - ground truth values of the waveforms.
%   waveforms_hat - estimated values of the waveforms.
%   params_true - ground truth values of the parameters.
%   params_hat - estimated values of the parameters.
%   solver_time - time to solve the optimization problem.
%   switching_times - switching points of regions.

%% Default parameters

if nargin < 2
    initial_delay = false;
end

if nargin < 3
    l2_reg = false;
end

if nargin < 4
    delay_length = 20;
end

if nargin < 5
    epsilon = 1e-3;
end

%% Uncomment if you want to load Gurobi and YALMIP every time (makes the script slower)
% set_gurobi
% set_yalmip

%% Loading waveforms
% filename = 'samplefile015.bin';
%[flow, volume, pao, pmus_true, insex] =  load_python_bin(filename);

flow = waveforms.flow;
volume = waveforms.volume;
pao = waveforms.pressure;
pmus_true = waveforms.pmus;
insex = waveforms.insex;

% k_soe stands for the index where the start of exhalation occurs
k_soe = find(diff(insex) <= -0.5);
k_soe = k_soe(1) + 1; % since we lose 1 index by using diff
padding = zeros(delay_length, 1);
if initial_delay
    flow = [padding; flow];
    volume = [padding; volume];
    pao = [padding; pao];
    pmus_true = [padding; pmus_true];
    insex = [padding; insex];
    k_soe = k_soe + delay_length; % to take into account the padding
end

% Converting units of airflow
flow = flow / 60 * 1000;

N = length(flow); % number of samples
if initial_delay
    Ns = 3; % number of time breaks
else
    Ns = 2; % number of time breaks
end
tik = binvar(N, Ns); % binary variables related to the switching instants

%% Constraining later switching times to occur after the earlier ones
constraint_occur = []; % t1 <= t2 <= ... <= tM
for i=2:1:Ns
	constraint_occur = constraint_occur +...
        [(1:N)* tik(1:end, i-1) <= (1:N)* tik(1:end, i)];
end

%% Constraining the switching instants to occur only once
constraint_unique = []; % only one '1' in t1, t2, ..., tM
for i=1:1:Ns
	constraint_unique = constraint_unique +...
        [sum(tik(1:end, i)) == 1];
end

%% Constraints to define the regions
pmus = sdpvar(N, 1);

constraint_regions = [];
if initial_delay
    for i=1:1:N-1
        constraint_regions = constraint_regions +...
            [implies(sum(tik(1:i, 1)) <= 0.5, [pmus(i) == 0, pmus(i+1) == 0]), ...
             implies(0.5 <= (sum(tik(1:i, 1)) - sum(tik(1:i, 2))), pmus(i+1) + epsilon <= pmus(i) ), ...
             implies(0.5 <= (sum(tik(1:i, 2)) - sum(tik(1:i, 3))), pmus(i) + epsilon <= pmus(i+1)), ...
             implies(0.5 <= sum(tik(1:i, 3)), [pmus(i) == 0, pmus(i+1) == 0])];
    end
else
    for i=1:1:N-1
        constraint_regions = constraint_regions +...
                 [implies(sum(tik(1:i, 1)) <= 0.5, pmus(i+1) + epsilon <= pmus(i)), ...
                  implies(0.5 <= (sum(tik(1:i, 1)) - sum(tik(1:i, 2))), pmus(i) + epsilon <= pmus(i+1)), ...
                  implies(0.5 <= sum(tik(1:i, 2)), [pmus(i) == pmus(i+1)])];
    end
end

%% Constraining the exhalation switching instant

tau_soe = 25;
constraint_exhalation = [(1:N) * tik(1:N, end) <= k_soe + tau_soe];

%% Cost function definition

resistance = sdpvar(1, 1);
elastance = sdpvar(1, 1);

cost = (pao - (pmus + resistance * flow + volume * elastance))' * ...
	   (pao - (pmus + resistance * flow + volume * elastance));
       
if l2_reg
    cost = cost + 1.0e-3 * (pmus' * pmus);
end

%% Constraining the real-valued variables

constraint_real = [ones(N,1)*(-20) <= pmus <= ones(N,1)*1, ...
    0 <= resistance, resistance <= 0.1, 0.005 <= elastance, elastance <= 1];       

%% Solving the optimization problem

options = sdpsettings;
solution = optimize([constraint_occur, constraint_unique, constraint_regions, constraint_real, constraint_exhalation],...
    cost, options);

%% Collecting outputs

solver_time = solution.solvertime;

waveforms_true.flow = flow;
waveforms_true.volume = volume;
waveforms_true.pao = pao;
waveforms_true.pmus = pmus_true;
waveforms_true.insex = insex;

waveforms_hat.pmus = value(pmus);
waveforms_hat.pao = waveforms_hat.pmus + value(resistance) * flow + volume * value(elastance);

params_lse = (([flow volume]' * [flow volume]) \ ([flow volume]')) * (pao - pmus_true);
params_true.resistance = params_lse(1) * 1000;
params_true.elastance = params_lse(2) * 1000;

params_hat.resistance = value(resistance) * 1000;
params_hat.elastance = value(elastance) * 1000;

switching_times = zeros(Ns, 1);
for i=1:Ns
    switching_times(i) = (1:N) * value(tik(1:N, i));
end

end