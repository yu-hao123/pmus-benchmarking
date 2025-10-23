function waveforms = extract_single_cycle(acq_table, ins_mark, next_ins_mark, exp_mark, peep, offset)
arguments
    acq_table
    ins_mark
    next_ins_mark
    exp_mark
    peep
    offset = 30
end

    interval = (ins_mark-offset):(next_ins_mark-offset-1);
    interval_table = acq_table(interval, :);
    interval_table.pmus = interval_table.pmus;
    interval_table.flow = fir_filter(8, 0.2, 100, interval_table.flow);
    interval_table.pressure = fir_filter(8, 0.2, 100, interval_table.pressure);
    interval_table.volume = fir_filter(8, 0.2, 100, interval_table.volume);

    exp_start = exp_mark - ins_mark + offset;

    waveforms = table();
    waveforms.time = interval_table.time;
    waveforms.paw = interval_table.pressure - peep;
    waveforms.pressure = waveforms.paw; % not ideal
    waveforms.flow = interval_table.flow;
    waveforms.volume = interval_table.volume - interval_table.volume(offset); % remove leaks, volume starts cycle in zero
    waveforms.pmus = interval_table.pmus;
    insexp = ones(length(waveforms.paw), 1);
    for i=1:length(waveforms.paw)
        if i >= exp_start
            insexp(i) = 0;
        end
    end
    waveforms.insexp = insexp;
end