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
offset = 60;
idx1 = 9000;
idx2 = 9001;

interval = (ins_marks(idx1)-offset):(ins_marks(idx2));
interval_table = acq_table(interval, :);
interval_table.pmus = interval_table.pmus;
interval_table.flow = fir_filter(8, 0.2, 100, interval_table.flow);
interval_table.pressure = fir_filter(8, 0.2, 100, interval_table.pressure);

plot_dataset(interval_table);

exp_start = exp_marks(idx1) - ins_marks(idx1) + offset;

peep = 5;
waveforms = table();
waveforms.time = interval_table.time;
waveforms.pressure = interval_table.pressure - peep;
waveforms.flow = interval_table.flow;
waveforms.volume = interval_table.volume;
waveforms.pmus = interval_table.pmus;
insexp = ones(length(waveforms.pressure), 1);
for i=1:length(waveforms.pressure)
    if i >= exp_start
        insexp(i) = 0;
    end
end
waveforms.insexp = insexp;

%% pmus miqp estimation
[waveforms_true, waveforms_hat, params_true, params_hat] = pmus_miqp(waveforms, false, true, 0);

pmus_optimized =  waveforms_hat.pmus;
pmus_true_miqp = waveforms_true.pmus;

%% pmus cubic estimation

%[waveforms_true, waveforms_hat, params_true, params_hat, solinfo] = pmus_cubic(waveforms, 43, 138);
%pmus_optimized = waveforms_hat.pmus;

%% plot
[f, t, linkplot] = plot_dataset(interval_table);

plot(linkplot(3), interval_table.time - interval_table.time(1), pmus_optimized);

resistance = params_true.resistance; % cmH2O / (L * s)
compliance = 1000 / params_true.elastance; % cmH2O / mL
pmus_recalculated = waveforms.pressure - resistance / 60 * waveforms.flow - waveforms.volume / (compliance);
plot(linkplot(3), interval_table.time - interval_table.time(1), pmus_recalculated);

legend('pmus ASL', 'pmus MIQP', 'pmus recalculated');

params_true
params_hat

cost_hat = cost(params_hat.resistance / 1000, params_hat.elastance / 1000, ...
    waveforms.pressure, waveforms.flow, waveforms.volume, waveforms_hat.pmus)

% considering PEEP = 0
% expects resistance in cmH2O / (mL * s) and elastance in cmH2O / mL
function J = cost(resistance, elastance, pressure, flow, volume, pmus)

    flow = flow * 1000 / 60;
    pressure_hat = pmus + resistance .* flow + (volume * elastance);
    resid = pressure - pressure_hat;
    J = sum(resid.^2);

end
