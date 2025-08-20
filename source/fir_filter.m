function y = fir_filter(order, cutoff, fs, x)

wn = cutoff / (fs/2);
filter_array = fir1(order, wn, 'low', hann(order+1));
y = filtfilt(filter_array, 1, x);

end