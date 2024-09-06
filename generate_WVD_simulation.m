dz = 2e-3;
df = -135e6;

poly = [-1.00660644521354e+43	3.43649507353042e+37	-4.29710160447565e+31...
    2.39878076665655e+25	1.40171350209608e+19	270665517408362];

fs = 4e9;
sweep_period = 1 / 400e3;
t = 0:1/fs:(sweep_period / 2)-1/fs;
freqs = polyval(poly, t);
c = 299792458;
N = 640;

freqs_lin = linspace(freqs(1), freqs(end), N);
t_new = interp1(freqs, t, freqs_lin);

w_m = 2 * pi * df;
alpha = (freqs(end) - freqs(1)) / (t(end) - t(1));
dalpha = polyval(polyder(poly), t_new);

pref = (1 ./ (dalpha + alpha)) - (1 / alpha);

[interferogram_res, interferogram_non] = gen_interferogram(df, dz, poly);

ascan_res = abs(fftshift(ifft(interferogram_res .* hann(length(interferogram_res)).')));
ascan_non = abs(fftshift(ifft(interferogram_non .* hann(length(interferogram_non)).')));

ascan_res = ascan_res(1:end/2);
ascan_non = ascan_non(1:end/2);

% ---

window = gausswin(129, 2);
[sst, f_sst] = wvd((interferogram_res), "smoothedPseudo", window, window, "MinThreshold", 0);
sst_non = wvd(hilbert(interferogram_non), "smoothedPseudo", window, window, "MinThreshold", 0);

err = (w_m * pref * c / 2) * N / sweep_period / 1e6;
err = err - mean(err);
err = interp1(linspace(0, 1, length(err)), err, linspace(0,1,size(sst, 2)), "makima");

k_min = 2 * pi / 1110e-9;
k_max = 2 * pi / 1010e-9;

k_axis = linspace(k_min, k_max, size(sst, 2)) / 1e6;

lambs = c ./ freqs;
l0 = median(lambs);
dl = lambs(1) - lambs(end);
z_max = l0^2 / dl * N / 4;

sst_ref = sst(:, end/2);
corr = normxcorr2(abs(sst_ref), abs(sst));
corr = corr(1:2:end, :);
corr_z_range = linspace(-z_max, z_max, size(corr, 1)) * 1e3;
corr_f_range = interp1(1:length(f_sst), f_sst, linspace(1, length(f_sst), size(corr, 1)));
corr_f_range = corr_f_range - mean(corr_f_range);

ridge = tfridge(corr, corr_f_range, 1e-3, "NumRidges", 1);
ridge = ridge - mean(ridge);

ridge_i = interp1(linspace(0, 1, length(ridge)), ridge, linspace(0, 1 ,length(interferogram_res)));
ridge_i = cumsum(ridge_i);
phi = exp(2i * pi * ridge_i);
interferogram = hilbert(interferogram_res) .* conj(phi);
ascan = abs(fftshift(ifft(interferogram .* hann(length(interferogram)).')));
ascan = ascan(1:end/2);
sst_cor = wvd(real(interferogram), "smoothedPseudo", window, window, "MinThreshold", 0);

maxx = max(ascan);
ascan = ascan / maxx;
ascan_res = ascan_res / maxx;
ascan_non = ascan_non / maxx;

z_range = linspace(0, z_max, size(sst, 1)) * 1e3;
ascan_z_range = linspace(0, z_max, length(ascan_res));
sample_range = linspace(1, N, size(sst, 2));

t_range = linspace(t(1), t(end), size(sst, 2)) * 1e6;
rad_range = linspace(0, 2*pi, size(sst, 1));

figure
hold on
f = tiledlayout(4, 3, "TileSpacing", "tight", "Padding", "tight");
ax = nexttile([1, 2]);
pcolor(t_range, z_range, abs(sst_non)/max(abs(sst_non(:))))
shading flat
xlabel("time (\mus)")
text(ax, t_range(10), z_range(end-70), "a", "FontSize", 12, "Color", [0.99 0.99 0.99])
cbar = colorbar;
cbar.Location = "south";
cbar.Position = cbar.Position + [0.55 0.12 -0.55 0];
cbar.Ticks = [0, 1];
cbar.Color = [1 1 1];

ax = nexttile;
plot(fliplr(ascan_non), ascan_z_range)
xlim([0, 1])
xticks([])
yticks([])

ax = nexttile([1, 2]);
pcolor(k_axis, z_range, abs(sst))
shading flat
text(ax, k_axis(10), z_range(end-70), "b", "FontSize", 12, "Color", [0.99 0.99 0.99])

ax = nexttile;
plot(fliplr(ascan_res), ascan_z_range)
xlim([0, 1])
xticks([])
yticks([])

ax = nexttile([1, 2]);
pcolor(k_axis, corr_z_range, corr)
shading flat
hold on
plot(k_axis, err * max(corr_z_range) / max(corr_f_range), "m--")
plot(k_axis, ridge * max(corr_z_range) / max(corr_f_range), "r")
l = legend("", "Theoretical Ridge", "Extracted Ridge");
l.Box = "off";
fontsize(l, 8, "points")
l.Position = l.Position + [0.28 0.0 0 0];

text(ax, k_axis(10), corr_z_range(end-70), "c", "FontSize", 12, "Color", [0.99 0.99 0.99])

nexttile
axis off

ax = nexttile([1, 2]);

pcolor(k_axis, z_range, sst_cor)
shading flat
xlabel("Wavenumber k (rad/\mum)")
text(ax, k_axis(10), z_range(end-70), "d", "FontSize", 12, "Color", [0.99 0.99 0.99])

nexttile
plot(fliplr(ascan), ascan_z_range)
xlim([0, 1])
yticks([])
xticks([0, 1])
% xticklabels({"0", "1"})
xlabel("Amplitude (a.u.)")

ylabel(f, "Distance z (mm)", "FontSize", 12)
% exportgraphics(f, "C:\Users\durban\OneDrive - Optos Inc\Documents\MATLAB\MatLab\Process B-Scans\Figures\SimulationWVD.png", "Resolution", 600)