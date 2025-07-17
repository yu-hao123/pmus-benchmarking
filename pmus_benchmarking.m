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
idx1 = 8000;
idx2 = 8001;
interval = (ins_marks(idx1)-offset):(ins_marks(idx2));
interval_table = acq_table(interval, :);

plot_dataset(interval_table);

exp_start = exp_marks(idx1) - ins_marks(idx1) + offset;

waveforms = table();
waveforms.pressure = interval_table.pressure;
waveforms.flow = interval_table.flow;
waveforms.volume = interval_table.volume;
waveforms.pmus = interval_table.pmus;

insexp = ones(length(waveforms.pressure), 1);
for i=1:length(waveforms.pressure)
    if i >= exp_start
        insexp(i) = 0;
    end
end
%%
waveforms.insexp = insexp;
[waveforms_true, waveforms_hat, params_true, params_hat] = pmus_miqp(waveforms, true, false, 0);

pmus_estimate_miqp =  waveforms_hat.pmus;
pmus_true_miqp = waveforms_true.pmus;

%%
[f, t, linkplot] = plot_dataset(interval_table);

plot(linkplot(3), interval_table.time - interval_table.time(1), pmus_estimate_miqp);
plot(linkplot(3), interval_table.time - interval_table.time(1), waveforms.insexp);
legend('pmus ASL', 'pmus MIQP', 'insexp')

