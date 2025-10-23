function [waveforms_true, waveforms_hat, params_true, params_fixed, solver_time, switching_times] = ...
    pmus_cubic_fixed(waveforms, k0, km, R_fixed, E_fixed, opts)
% PMUS_CUBIC_FIXED  Fit a cubic pmus on [k0,km] with R and E fixed.

if nargin < 6, opts = struct(); end
if ~isfield(opts,'pmus_bounds'), opts.pmus_bounds = [-50 1]; end
if ~isfield(opts,'solver'), opts.solver = 'quadprog'; end
if ~isfield(opts,'verbose'), opts.verbose = false; end

t = waveforms.time(:);
N = length(t);
if k0 < 1 || k0 >= km || km > N
    error('k0 and km must be valid indices with 1 <= k0 < km <= N');
end

flow = waveforms.flow(:) / 60 * 1000; % L/min -> mL/s
paw  = waveforms.pressure(:);
volume = waveforms.volume(:);

idx_effort = k0:km;

% normalized time s in [0,1] only inside effort window
s = zeros(N,1);
L = length(idx_effort);
s(idx_effort) = (0:(L-1))' ./ (L-1);  % 0..1

a0 = sdpvar(1,1);
a1 = sdpvar(1,1);
a2 = sdpvar(1,1);
a3 = sdpvar(1,1);

phi = [ones(N,1), s, s.^2, s.^3];
pmus_all = phi * [a0; a1; a2; a3];

constraints = [];

% outside effort window pmus == 0
outside_idx = setdiff(1:N, idx_effort);
if ~isempty(outside_idx)
    constraints = [constraints, pmus_all(outside_idx) == 0];
end

% endpoint constraints in normalized domain
% a0 == 0  (pmus at s=0)
% a0 + a1 + a2 + a3 == 0  (pmus at s=1)
constraints = [constraints, a0 == 0];
constraints = [constraints, a0 + a1 + a2 + a3 == 0];

% derivative sign heuristics with a small margin to avoid equality
constraints = [constraints, a1 <= -1e-6];
constraints = [constraints, a2 >= 1e-6];
constraints = [constraints, a3 <= -1e-6];

% pmus bounds inside effort
pmin = opts.pmus_bounds(1);
pmax = opts.pmus_bounds(2);
constraints = [constraints, pmin <= pmus_all(idx_effort), pmus_all(idx_effort) <= pmax];

R = R_fixed; E = E_fixed;
paw_model = pmus_all + R*flow + E*volume;
res = paw - paw_model;
objective = res' * res;  % least-squares

if strcmpi(opts.solver,'gurobi')
    options = sdpsettings('solver','gurobi', ...
                          'verbose', 0, ...
                          'gurobi.TimeLimit', 60, ...
                          'gurobi.OutputFlag', 0);
else
    options = sdpsettings('solver','quadprog', 'verbose', 0);
end

% Solve
diagn = optimize(constraints, objective, options);
solver_time = diagn.solvertime;

% prepare outputs in case of failure
waveforms_hat = [];
switching_times = [k0; km];

if diagn.problem ~= 0
    % solver failed or infeasible
    warning('pmus_cubic_fixed: solver returned problem %d (%s)', diagn.problem, yalmiperror(diagn.problem));

    % outputs: empty waveforms_hat, solver_time already set
    waveforms_true.flow = flow;
    waveforms_true.volume = volume;
    waveforms_true.paw = paw;
    if isfield(waveforms,'pmus'), waveforms_true.pmus = waveforms.pmus; end
    params_true = struct();
    params_fixed = struct('resistance', R*1000, 'elastance', E*1000);
    return
end

% results
A = [value(a0); value(a1); value(a2); value(a3)];
R_hat = R;
E_hat = E;
pmus_hat = double(phi * A);

% outputs
waveforms_true.flow = flow;
waveforms_true.volume = volume;
waveforms_true.paw = paw;
if isfield(waveforms,'pmus')
    waveforms_true.pmus = waveforms.pmus;
end

waveforms_hat = table();
waveforms_hat.pmus = pmus_hat;
waveforms_hat.paw  = pmus_hat + R_hat*flow + E_hat*volume;

% params_true computed by LSE on true pmus if available
params_lse = (([flow volume]' * [flow volume]) \ ([flow volume]')) * (paw - waveforms.pmus);
params_true.resistance = params_lse(1) * 1000;
params_true.elastance  = params_lse(2) * 1000;

params_fixed.resistance = R_hat * 1000;
params_fixed.elastance  = E_hat * 1000;
params_fixed.A = A;
params_fixed.k0 = k0;
params_fixed.km = km;
params_fixed.tm = t(km) - t(k0);

end
