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

interval = (ins_marks(idx1)-offset):(ins_marks(idx2)-offset);
interval_table = acq_table(interval, :);
interval_table.pmus = interval_table.pmus;
interval_table.flow = fir_filter(8, 0.2, 100, interval_table.flow);
interval_table.pressure = fir_filter(8, 0.2, 100, interval_table.pressure);

plot_dataset(interval_table);

exp_start = exp_marks(idx1) - ins_marks(idx1) + offset;

peep = 5;
waveforms = table();
waveforms.time = interval_table.time;
waveforms.paw = interval_table.pressure - peep;
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

%% pmus miqp estimation
[waveforms_true, waveforms_hat, params_true, params_hat] = pmus_miqp(waveforms, false, true, 0);

paw_est = waveforms_hat.paw;
pmus_optimized =  waveforms_hat.pmus;
pmus_true_miqp = waveforms_true.pmus;

%% pmus cubic estimation

%[waveforms_true, waveforms_hat, params_true, params_hat, solinfo] = pmus_cubic(waveforms, 43, 175);
%pmus_optimized = waveforms_hat.pmus;
%paw_est = waveforms_hat.paw;

%% plot
[f, t, linkplot] = plot_dataset(interval_table);

plot(linkplot(3), interval_table.time - interval_table.time(1), pmus_optimized);

resistance = params_true.resistance; % cmH2O / (L * s)
compliance = 1000 / params_true.elastance; % cmH2O / mL
pmus_recalculated = waveforms.paw - resistance / 60 * waveforms.flow - waveforms.volume / (compliance);
plot(linkplot(1), interval_table.time - interval_table.time(1), paw_est + peep);
legend(linkplot(1), 'paw ASL', 'paw MIQP est');

plot(linkplot(3), interval_table.time - interval_table.time(1), pmus_recalculated);
plot(linkplot(3), interval_table.time - interval_table.time(1), waveforms.pmus_mag);
plot(linkplot(3), interval_table.time - interval_table.time(1), waveforms.insexp);
legend('pmus ASL', 'pmus MIQP', 'pmus LS recalculated', 'pmus MAG', 'expiration switch');

fprintf("true parameters - \n");
fprintf("    resistance: %.2f \n", params_true.resistance);
fprintf("    compliance: %.2f \n", 1000 / params_true.elastance);

fprintf("estimated parameters - \n");
fprintf("    resistance: %.2f \n", params_hat.resistance);
fprintf("    compliance: %.2f \n", 1000 / params_hat.elastance);

cost_hat = cost(params_hat.resistance / 1000, params_hat.elastance / 1000, ...
    waveforms.paw, waveforms.flow, waveforms.volume, waveforms_hat.pmus);

% considering PEEP = 0
% expects resistance in cmH2O / (mL * s) and elastance in cmH2O / mL
function J = cost(resistance, elastance, paw, flow, volume, pmus)
    flow = flow * 1000 / 60;
    paw_hat = pmus + resistance .* flow + (volume * elastance);
    resid = paw - paw_hat;
    J = sum(resid.^2);
end
