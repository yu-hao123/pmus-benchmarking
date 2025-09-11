clc;
clear all;
close all;

addpath('source');

%run set_gurobi.m
%run set_yalmip.m

load("data/ASL_spont_01.mat");
set(0, 'DefaultLineLineWidth', 0.8);

%%
time = acq_table.time;
pressure = acq_table.pressure;
flow = acq_table.flow;
volume = acq_table.volume;
pmus = acq_table.pmus;
pmus_estimate = acq_table.pmus_estimate; % PMUS-MAG

[ins_marks, exp_marks] = retrieve_parity_marks(volume * 10);
fprintf("retrieved %d ins/exp marks from volume parity bit\n", length(ins_marks));

%%
offset = 50;
start_idx = 15171;
finish_idx = 15176;
peep = 5;
plot_dataset(acq_table(ins_marks(start_idx)-offset:ins_marks(finish_idx+1)-offset, :));

%% pmus miqp estimation
waveforms = table();
waveforms_hat = table();
params_true = table();
params_hat = table();
cost_hat = [];
% estimate multiple sequential cycles
for i=start_idx:finish_idx
    fprintf(" optimizing cycle of idx: %d \n", i);

    cycle = extract_cycle(acq_table, ins_marks(i), ins_marks(i + 1), exp_marks(i), peep);
    [cycle_true, cycle_hat, cycle_params_true, cycle_params_hat] = pmus_miqp(cycle, false, true, 0);
    cycle_cost_hat = cost(cycle_params_hat.resistance / 1000, cycle_params_hat.elastance / 1000, ...
        cycle.paw, cycle.flow, cycle.volume, cycle_hat.pmus);

    cost_hat = [cost_hat; cycle_cost_hat];

    resistance = cycle_params_true.resistance; % cmH2O / (L * s)
    compliance = 1000 / cycle_params_true.elastance; % cmH2O / mL
    pmus_recalculated = cycle.paw - resistance / 60 * cycle.flow - cycle.volume / (compliance);
    cycle.pmus_recalculated = pmus_recalculated;

    if (i == start_idx)
        waveforms = cycle;
        waveforms_hat = cycle_hat;
        params_true = cycle_params_true;
        params_hat = cycle_params_hat;
    else
        waveforms = [waveforms; cycle];
        waveforms_hat = [waveforms_hat; cycle_hat];
        params_true = [params_true; cycle_params_true];
        params_hat = [params_hat; cycle_params_hat];
    end

end


paw_est = waveforms_hat.paw;
pmus_optimized =  waveforms_hat.pmus;

%% pmus cubic estimation

%[waveforms_true, waveforms_hat, params_true, params_hat, solinfo] = pmus_cubic(waveforms, 43, 175);
%pmus_optimized = waveforms_hat.pmus;
%paw_est = waveforms_hat.paw;

%% plot and display selected R, C parameters
[f, t, linkplot] = plot_dataset(waveforms);

plot(linkplot(3), waveforms.time - waveforms.time(1), pmus_optimized);
plot(linkplot(1), waveforms.time - waveforms.time(1), paw_est);
legend(linkplot(1), 'paw ASL', 'paw MIQP est');

plot(linkplot(3), waveforms.time - waveforms.time(1), waveforms.pmus_recalculated);
plot(linkplot(3), waveforms.time - waveforms.time(1), waveforms.pmus_mag);
plot(linkplot(3), waveforms.time - waveforms.time(1), waveforms.insexp);
legend('pmus ASL', 'pmus MIQP', 'pmus LS recalculated', 'pmus MAG', 'expiration switch');

cidx = 3;
fprintf("true parameters - \n");
fprintf("    resistance: %.2f \n", params_true(cidx).resistance);
fprintf("    compliance: %.2f \n", 1000 / params_true(cidx).elastance);

fprintf("estimated parameters - \n");
fprintf("    resistance: %.2f \n", params_hat(cidx).resistance);
fprintf("    compliance: %.2f \n", 1000 / params_hat(cidx).elastance);


% considering PEEP = 0
% expects resistance in cmH2O / (mL * s) and elastance in cmH2O / mL
function J = cost(resistance, elastance, paw, flow, volume, pmus)
    flow = flow * 1000 / 60;
    paw_hat = pmus + resistance .* flow + (volume * elastance);
    resid = paw - paw_hat;
    J = sum(resid.^2);
end

function waveforms = extract_cycle(acq_table, ins_mark, next_ins_mark, exp_mark, peep, offset)
arguments
    acq_table
    ins_mark
    next_ins_mark
    exp_mark
    peep
    offset = 50
end

interval = (ins_mark-offset):(next_ins_mark-offset-1);
interval_table = acq_table(interval, :);
interval_table.pmus = interval_table.pmus;
interval_table.flow = fir_filter(8, 0.2, 100, interval_table.flow);
interval_table.pressure = fir_filter(8, 0.2, 100, interval_table.pressure);

exp_start = exp_mark - ins_mark + offset;

waveforms = table();
waveforms.time = interval_table.time;
waveforms.paw = interval_table.pressure - peep;
waveforms.pressure = waveforms.paw; % not ideal
waveforms.flow = interval_table.flow;
waveforms.volume = interval_table.volume;
waveforms.pmus = interval_table.pmus;
waveforms.pmus_mag = interval_table.pmus_estimate;
insexp = ones(length(waveforms.paw), 1);
for i=1:length(waveforms.paw)
    if i >= exp_start
        insexp(i) = 0;
    end
end
waveforms.insexp = insexp;

end
