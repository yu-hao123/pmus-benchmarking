function [waveforms_true, waveforms_hat, params_true, params_hat, solinfo] = ...
         pmus_cubic(waveforms, k0, km, opts)
% PMUS_CUBIC  Jointly estimate cubic pmus and R,E given effort indices
%
% Inputs:
%   waveforms - struct with fields: flow, volume, pressure, pmus (true)
%   t         - time vector (length N)
%   k0        - index of effort start
%   km        - index of effort finish (relaxation)
%   opts      - optional struct with fields:
%               R_bounds = [R_min R_max]  (default [0 0.1])
%               E_bounds = [E_min E_max]  (default [0.005 1])
%               pmus_bounds = [pmin pmax] (default [-50 50])
%               solver = 'quadprog' or 'gurobi' (default 'quadprog')

if nargin < 5, opts = struct(); end
if ~isfield(opts,'R_bounds'),   opts.R_bounds = [0 0.1]; end
if ~isfield(opts,'E_bounds'),   opts.E_bounds = [0.005 0.1]; end
if ~isfield(opts,'pmus_bounds'),opts.pmus_bounds = [-50 1]; end
if ~isfield(opts,'solver'),     opts.solver = 'quadprog'; end

t = waveforms.time;
N = length(t);
if k0 < 1 || k0 >= km || km > N
    error('k0 and km must be valid indices with 1 <= k0 < km <= N');
end

flow = waveforms.flow(:) / 60 * 1000; % ensure column (from L/min to mL/s)
paw  = waveforms.pressure(:);
volume = waveforms.volume(:);

% normalized time s in [0,1]
s = zeros(N,1);
idx_effort = k0:km;
tm = t(km) - t(k0); % real duration
if tm <= 0
    error('Invalid indices: t(km) must be > t(k0)');
end

% normalized s on effort window
s(idx_effort) = ((0:(length(idx_effort)-1))') / (length(idx_effort)-1);  % 0..1
% outside effort window s stays zero and forcing pmus=0

a0 = sdpvar(1,1);
a1 = sdpvar(1,1);
a2 = sdpvar(1,1);
a3 = sdpvar(1,1);
R  = sdpvar(1,1);
E  = sdpvar(1,1);

% build pmus(s) = a0 + a1*s + a2*s^2 + a3*s^3 (vector for all samples)
phi = [ones(N,1), s, s.^2, s.^3];   % N x 4
pmus_all = phi * [a0; a1; a2; a3];

constraints = [];
% pmus outside effort window is zero
outside_idx = setdiff(1:N, idx_effort);
if ~isempty(outside_idx)
    constraints = [constraints, pmus_all(outside_idx) == 0];
end

% endpoint constraints in normalized domain:
% pmus(s=0) = a0 = 0
% pmus(s=1) = a0 + a1 + a2 + a3 = 0
constraints = [constraints, a0 == 0];
constraints = [constraints, a0 + a1 + a2 + a3 == 0];

% p'(t0) -> a1 < 0
% p''(t0) -> a2 > 0
% p'''(t0) -> a3 < 0
% small numeric margins to avoid strict inequalities
constraints = [constraints, a1 <= -1e-6];
constraints = [constraints, a2 >= 1e-6];
constraints = [constraints, a3 <= -1e-6];

constraints = [constraints, opts.R_bounds(1) <= R, R <= opts.R_bounds(2)];
constraints = [constraints, opts.E_bounds(1) <= E, E <= opts.E_bounds(2)];

pmin = opts.pmus_bounds(1);
pmax = opts.pmus_bounds(2);
constraints = [constraints, pmin <= pmus_all(idx_effort), pmus_all(idx_effort) <= pmax];

% model predicted paw
paw_model = pmus_all + R*flow + E*volume;
res = paw - paw_model;
objective = res' * res;

if strcmpi(opts.solver,'gurobi')
    options = sdpsettings('solver','gurobi', 'gurobi.OutputFlag', 0);
else
    options = sdpsettings('solver','quadprog','verbose', 0);
end
diagn = optimize(constraints, objective, options);

solinfo = struct();
solinfo.problem = diagn.problem;
solinfo.info = yalmiperror(diagn.problem);
solinfo.solvertime = diagn.solvertime;

if diagn.problem ~= 0
    warning('Solver returned problem: %d (%s)', diagn.problem, yalmiperror(diagn.problem));
end

% extract results
A = [value(a0); value(a1); value(a2); value(a3)];
R_hat = value(R);
E_hat = value(E);
pmus_hat = double(phi * A);

waveforms_true.flow = flow;
waveforms_true.volume = volume;
waveforms_true.paw = paw;
if isfield(waveforms,'pmus')
    waveforms_true.pmus = waveforms.pmus;
end

waveforms_hat = table();
waveforms_hat.pmus = pmus_hat;
waveforms_hat.paw  = pmus_hat + R_hat*flow + E_hat*volume;

% params_true (least squares baseline) computed using available true pmus
params_lse = (([flow volume]' * [flow volume]) \ ([flow volume]')) * (paw - waveforms.pmus);
params_true.resistance = params_lse(1) * 1000;
params_true.elastance = params_lse(2) * 1000;

params_hat.resistance = R_hat * 1000;
params_hat.elastance  = E_hat * 1000;
params_hat.A = A;

end