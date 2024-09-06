function [interferogram_res, interferogram_non] = gen_interferogram(df, dz, poly)
N = 640; % Number of frequency sampling points, 640 for 400kHz source

% This polynomial describes the frequency sweep of the actual source
if nargin <= 2
    poly = [-1.00660644521354e+43	3.43649507353042e+37	-4.29710160447565e+31...
    2.39878076665655e+25	1.40171350209608e+19	270665517408362];
end

% Uncomment below line to simulate a linear sweep instead;
% poly = [1.9403e+19 270665517408362];

% Sampling frequency should be high enough to minimise numerical error
fs = 4e9;
% 400kHz source
sweep_period = 1 / 400e3;

% Time axis
t = 0:1/fs:(sweep_period / 2)-1/fs;

c = 299792458;
% Factor of 2 required to account for double-pass
dt = dz * 2 / c;

freqs = polyval(poly, t);

% Linear frequency vector, N is number of frequency sampling points
freqs_lin = linspace(freqs(1), freqs(end), N);
diff_freqs = polyval(polyder(poly), t);

% This represents times at which the interferogram is sampled by the
% k-clock
t_new = interp1(freqs, t, freqs_lin);

lambs = c ./ freqs;

lambda = median(lambs);
bandwidth = lambs(1) - lambs(end);

% Maximum detectable distance, and distance sampling interval
z_max = lambda^2 / bandwidth * N / 4;
z_sampling = z_max / (0.5 * N);

plotlim = z_max * 1e3;

% poly1 describes frequency in reference arm, which is offset by df
poly1 = poly;
poly1(end) = poly(end) + df;

% phase of reference arm (we are integrating the frequency sweep to get phase)
p1 =  2 * pi * polyval(polyint(poly1), t);

% phase of reflectors in sample arm (same as above)
p2 = 2 * pi * polyval(polyint(poly), t - dt);
p3 = 2 * pi * polyval(polyint(poly), t - (dt + 30e-13));
p4 = 2 * pi * polyval(polyint(poly), t - (dt - 30e-13));

% interferogram is formed by taking the difference between reference arm
% phase and reflector phase, see Drexler Fujimoto Eqn 2.9 pp.74
interferogram = (cos(p1 - p2) + cos(p1 - p3) + cos(p1 - p4));

interferogram_non = interp1(linspace(0, 1, length(interferogram)), interferogram.', linspace(0, 1, N), "makima");

% Interferogram is resampled at fixed frequency intervals
interferogram_res = interp1(t, interferogram.', t_new, "makima");

ascan = abs(fftshift(ifft(interferogram_res .* hann(length(interferogram_res)).')));
ascan_non = abs(fftshift(ifft(interferogram_non .* hann(length(interferogram_non)).')));
end

