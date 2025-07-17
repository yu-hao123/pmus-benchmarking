function [ins_marks, exp_marks] = retrieve_parity_marks(volume)

ins_marks = [];
exp_marks = [];

parity = mod(volume(1), 2);

for i=2:length(volume)
    % sanity verification given some datasets (MAG-ASL) corrupted volume data
    if volume(i) == 0.0
        continue
    end

    new_parity = mod(volume(i), 2);
    if (new_parity ~= parity)
        if (parity == 0)
            ins_marks = [ins_marks; i];
        else
            exp_marks = [exp_marks; i];
        end
    end
    parity = new_parity;
end

end