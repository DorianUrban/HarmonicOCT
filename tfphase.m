function [ridge_out] = tfphase(frame_in, I_in)
%tfphase Summary of this function goes here
%   Detailed explanation goes here

window = gausswin(65, 2);
window = hann(65, "symmetric");
% window = gausswin(32, 2);
rng = 1;
rngr = -rng:1:rng;

maxima = max(abs(frame_in), [], 1);
sorted = sort(maxima, "descend");

if nargin < 2
    I = find(maxima == sorted(3));
else
    I = find(maxima == sorted(I_in));
end

fringe = frame_in(:, floor(I));
pwv = wvd(fringe, "smoothedPseudo", window, window, "MinThreshold", 0);
% pwv = stft(fringe, [], window=window, OverlapLength=length(window) - 1, FFTLength=2 * length(window));

mean_corr = zeros(size(tfxcorr(pwv, 0)));

% for ii = rngr

fringe = frame_in(:, I);

[pwv, f_pwv, t_pwv] = wvd(hilbert(real((fringe))), "smoothedPseudo", window, window);

pwv = abs(pwv).^2;
mean_corr = mean_corr + tfxcorr(pwv, 0);

% end
% 
figure
pcolor(t_pwv, f_pwv, pwv)
ylabel("Frequency (rad/sample)")
xlabel("Sample Number")
shading flat

f_corr = linspace(0, 1, size(mean_corr, 1));

ridge = tfridge(mean_corr, f_corr, .1, "NumRidges", 3, "NumFrequencyBins", 30);
ridge = ridge(:, 1);
ridge_i = interp1(1:length(ridge), ridge, linspace(1, length(ridge), size(mean_corr, 2)));
ridge_out = cumsum(ridge_i - mean(ridge_i));

% figure
% plot(ridge_i - mean(ridge_i))

f_pwv = linspace(0, 1, size(mean_corr, 1));
t_pwv = interp1(linspace(0, 1, length(t_pwv)), t_pwv, linspace(0, 1, size(mean_corr, 2)));
% 
figure
pcolor(t_pwv, f_pwv, mean_corr)
shading flat
hold on
plot(t_pwv, ridge.', "r--")
ylabel("Lag (A.U.)")
xlabel("Sample Number")

end

