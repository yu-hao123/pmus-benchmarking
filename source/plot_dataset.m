function [f, t, link_plot] = plot_dataset(acq_table, ins_marks, exp_marks, options)
%PLOT_DATASET Plot all the waveforms of an acquisition table (paw, flow, pmus/vol)
%   acq_table: contains
%           t: tiled layout
%           f: resulting plot object

arguments
    acq_table
    ins_marks = []
    exp_marks = []
    options.LastPlot {mustBeMember(options.LastPlot, ["vol", "pmus"])} = "pmus"
    options.Position = [70 170 800 600]
end

time = acq_table.time;
time = time - time(1);
pressure = acq_table.pressure;
volume = acq_table.volume;
flow = acq_table.flow;
pmus = acq_table.pmus;

ie_marks_present = true;

if isempty(ins_marks) || isempty(exp_marks)
    ie_marks_present = false;
end

f = figure('Renderer', 'painters', 'Position', options.Position);
t = tiledlayout(3,1);

link_plot(1) = nexttile;
plot(time, pressure, 'k'); grid on; hold on;
set(gca,'TickLabelInterpreter','latex','FontSize',12);
set(gca,'XTickLabel',[]);
ylabel('$P_{aw}$ $[cmH_{2}O]$','Interpreter','latex');

if ie_marks_present
    for i = 1:length(ins_marks)
        xline(time(ins_marks(i)),'--r', 'HandleVisibility', 'off');
    end
    for i = 1:length(exp_marks)
        xline(time(exp_marks(i)),'--b', 'HandleVisibility', 'off');
    end
end

link_plot(2) = nexttile;
plot(time, flow, 'k', 'HandleVisibility', 'off'); grid on; hold on;
set(gca,'TickLabelInterpreter','latex','FontSize',12);
set(gca,'XTickLabel',[]);
ylabel('$\dot{V}$ $[L/min]$','Interpreter','latex');
if ie_marks_present
    for i = 2:length(ins_marks)
        xline(time(ins_marks(i)),'--r', 'HandleVisibility', 'off');
    end
    for i = 2:length(exp_marks)
        xline(time(exp_marks(i)),'--b', 'HandleVisibility', 'off');
    end
    xline(time(ins_marks(1)),'--r', 'HandleVisibility', 'on');
    xline(time(exp_marks(1)),'--b', 'HandleVisibility', 'on');
end

link_plot(3) = nexttile;
hold on; grid on;

switch options.LastPlot
    case "pmus"
        plot(time, pmus, 'k', 'HandleVisibility', 'on');
        ylabel('$P_{mus}$ $[cmH_{2}O]$','Interpreter','latex');
    case "vol"
        plot(time, volume, 'k', 'HandleVisibility', 'on');
        ylabel('$V$ $[mL]$','Interpreter','latex');
end

set(gca,'TickLabelInterpreter','latex','FontSize',12);

if ie_marks_present
    for i = 1:length(ins_marks)
        xline(time(ins_marks(i)),'--r', 'HandleVisibility', 'off');
    end
    for i = 1:length(exp_marks)
        xline(time(exp_marks(i)),'--b', 'HandleVisibility', 'off');
    end
end

linkaxes(link_plot, 'x');

t.Padding = 'compact';
t.TileSpacing = 'compact';

%align_Ylabels(f);
end
