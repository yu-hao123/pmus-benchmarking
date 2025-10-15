clc;
clear all;
close all;

addpath('source');

% run set_gurobi.m
% run set_yalmip.m

load("data/ASL_spont_01.mat");
set(0, 'DefaultLineLineWidth', 0.8);

%% define estimator
estimator = "cubic"; % miqp or cubic
cubic_k0 = 35;
cubic_km = 130;
heatmap_dim = 48;

%% Prepare waveforms
time = acq_table.time;
pressure = acq_table.pressure;
flow = acq_table.flow;
volume = acq_table.volume;
pmus_true = acq_table.pmus;

[ins_marks, exp_marks] = retrieve_parity_marks(volume * 10);
fprintf("retrieved %d ins/exp marks from volume parity bit\n", length(ins_marks));

offset = 50;
idx1 = 8099;
idx2 = 8100;
interval = (ins_marks(idx1)-offset):(ins_marks(idx2)-offset);
interval_table = acq_table(interval, :);

interval_table.flow = fir_filter(8, 0.2, 100, interval_table.flow);
interval_table.pressure = fir_filter(8, 0.2, 100, interval_table.pressure);

exp_start = exp_marks(idx1) - ins_marks(idx1) + offset;

waveforms = struct();
waveforms.time = interval_table.time;
waveforms.paw   = interval_table.pressure - 5;
waveforms.pressure = waveforms.paw; % not ideal
waveforms.flow  = interval_table.flow;
waveforms.volume = interval_table.volume;
waveforms.pmus  = interval_table.pmus;

insexp = ones(length(waveforms.paw), 1);
insexp(exp_start:end) = 0;
waveforms.insexp = insexp;

plot_dataset(interval_table);

%% Grid search for R and C
if isempty(gcp("nocreate"))
    parpool(8);
end

R_values = linspace(5, 50, heatmap_dim); % Resistance [(cmH20.s) / mL]
C_values = linspace(10, 80, heatmap_dim); % Compliance [ml / cmH20]
error_matrix = nan(length(C_values), length(R_values));

% Setup progress tracking
dq = parallel.pool.DataQueue;
N_total = length(R_values) * length(C_values);

afterEach(dq, @(~)updateProgress(N_total));
    function updateProgress(N_total)
        persistent progress
        if isempty(progress)
            progress = 0;
        end
        progress = progress + 1;
        fprintf("Progress: %d / %d (%.1f%%) \n", progress, N_total, 100 * progress/N_total);
    end

fprintf("initiating parallel estimation method: " + estimator + "\n");
parfor iR = 1:length(R_values)
    R = R_values(iR);
    temp_error_col = nan(length(C_values), 1);

    for iC = 1:length(C_values)
        C = C_values(iC);
        if (estimator == "miqp")
            [~, waveforms_hat, ~, ~, ~, ~] = pmus_miqp_fixed(waveforms, false, true, 0, 1e-3, R/1000, 1/C);
        else
            [~, waveforms_hat, ~, ~, ~, ~] = pmus_cubic_fixed(waveforms, cubic_k0, cubic_km, R/1000, 1/C);
        end

        if isempty(waveforms_hat)
            send(dq, 1);
            continue;
        end

        err = norm(waveforms.paw - waveforms_hat.paw);
        temp_error_col(iC) = err;

        % update progress
        send(dq, 1);
    end

    error_matrix(:, iR) = temp_error_col;
end
%%
[best_error, linear_idx] = min(error_matrix(:));
[best_iC, best_iR] = ind2sub(size(error_matrix), linear_idx);

%% Heatmap
best_R = R_values(best_iR);
best_C = C_values(best_iC);

fprintf("Best R: %.2f (cmH2O.s)/mL \n", best_R);
fprintf("Best C: %.2f mL/cmH2O\n", best_C);

figure;
imagesc(R_values, C_values(1:end), log10(error_matrix)); hold on;
set(gca, 'YDir', 'normal');
colorbar;
xlabel('Resistance R (cmH2O.s)/mL');
ylabel("Compliance C (mL/cmH2O)");
title(estimator + " Log10 (residual cost) surface");

if (estimator == "miqp")
    [~, waveforms_hat, params_true, params_hat, ~, ~] = pmus_miqp_fixed(waveforms, false, true, 0, 1e-3, best_R/1000, 1/best_C);
else
    [~, waveforms_hat, params_true, params_hat, ~, ~] = pmus_cubic_fixed(waveforms, cubic_k0, cubic_km, best_R/1000, 1/best_C);
end

fprintf("params true - R = %.2f, C = %.2f\n", params_true.resistance, 1000/params_true.elastance);
fprintf("cost (R = %.2f, C = %.2f) : %.2f \n", best_R, best_C, norm(waveforms_hat.paw - waveforms.paw));
plot(params_true.resistance, 1000/params_true.elastance, 'r*', 'MarkerSize', 8, 'LineWidth', 2);
plot(best_R, best_C, 'wo', 'MarkerSize', 8, 'LineWidth', 2);

pmus_optimized = waveforms_hat.pmus;
pmus_candidate = waveforms.pressure ...
            - best_R * waveforms.flow / 60 ...
            - waveforms.volume / best_C;
paw_est = waveforms_hat.paw;
[f, t, linkplot] = plot_dataset(waveforms);
plot(linkplot(3), waveforms.time - waveforms.time(1), pmus_optimized);
%plot(linkplot(3), waveforms.time - waveforms.time(1), pmus_candidate);
plot(linkplot(1), waveforms.time - waveforms.time(1), paw_est);
legend(linkplot(1), "paw ASL", "paw est");
legend(linkplot(3), "pmus ASL", "pmus est");
sgtitle("best estimation (" + estimator + ")");