close all

directory = "\\skat\research team\DUrban\MatLab\Process B-Scans\Flattening Grape\20240408-164353-00";

% f = Frame(directory, 3);
% f = f.disp_comp();

f_blurred = get_average_frame(directory, 10, 6, false);
f_sharp = get_average_frame(directory, 10, 6);

h = findall(groot, type="figure");

stftx = h(2).Children.Children.XData;
stfty = h(2).Children.Children.YData;
stftz = h(2).Children.Children.CData;

xcorrx = h(1).Children.Children(2).XData;
xcorry = h(1).Children.Children(2).YData;
xcorrz = h(1).Children.Children(2).CData;

ridgex = h(1).Children.Children(1).XData;
ridgey = h(1).Children.Children(1).YData;


img_sharp = logScaling(f_sharp(end/2:end, end/2:end));
img_blurred = logScaling(f_blurred(end/2:end, 1:end/2), f_sharp(end/2:end, end/2:end));

x_extent = 400:1700;
px_per_mm_x = 4000/12/1e3; % BscanLen / Imaging Range / mm
px_per_mm_y = 640/4/1e3; % AscanLen / AscanRange / mm
aspect_ratio = px_per_mm_x / px_per_mm_y;

z_max = 2;
k_min = 2 * pi / 1110e-9;
k_max = 2 * pi / 1010e-9;

k_axis = linspace(k_min, k_max, size(stftz, 2)) / 1e6;

figure
f = tiledlayout(2, 2, tilespacing="tight", padding="tight");

ax = nexttile;
pcolor(k_axis, stfty * z_max, flipud(stftz./vecnorm(stftz, "Inf", 1)))
shading flat
ylabel("z (mm)", "FontSize", 10)
yticks([0, 1, 2])
yticklabels({"0", "1", "2"})
xlabel("Central Wavenumber (rad/\mum)", "FontSize", 10)
% daspect([1, 5 1])

hold on
[~, I] = max(max(abs(stftz), [],  1));
xline(k_axis(I), "r--")

text(ax, 5.665, 1.85, "a", "FontSize", 12, "Color", [0.99 0.99 0.99])

ax = nexttile;
pcolor(k_axis, ((xcorry - 0.5) * 4), flipud(xcorrz))
shading flat
hold on
plot(k_axis, -1 * ((ridgey - 0.5) * 4), "r--")
ylabel("\Delta z (mm)", "FontSize", 10)
yticks([-2, 0, 2])
yticklabels({"-2", "0", "2"})
% daspect([1, 5, 1])
xlabel("Central Wavenumber (rad/\mum)", "FontSize", 10)

text(ax, 5.665, 1.75, "b", "FontSize", 12, "Color", [0.99 0.99 0.99])

cbar = colorbar(ax);
cbar.Label.String = "";
cbar.Limits = [0, 1];
cbar.Ticks = [0, 1];
% cbar.Position = cbar.Position + [-0 0 0 0.22];

ax = nexttile;
imshow(img_blurred(:, x_extent), colormap=turbo, DisplayRange=[min(img_sharp(:)), max(img_sharp(:))])
daspect([1, 1/aspect_ratio 1])

text(ax, 10, 31, "c", "FontSize", 12, "Color", [0.99 0.99 0.99])

ax = nexttile;
imshow(img_sharp(:, x_extent), colormap=turbo, DisplayRange=[])
daspect([1, 1/aspect_ratio 1])

text(ax, 10, 31, "d", "FontSize", 12, "Color", [0.99 0.99 0.99])

cbar2 = colorbar(ax);
cbar2.Label.String = "";
cbar2.Limits = [0, 1];
cbar2.Ticks = [0, 1];
% cbar2.Position = cbar2.Position + [-0.00 0 0 0.22];

bar_length = 1000;

hbar = scalebar(ax, "y");
hbar.Visible = "off";
hbar.Axis = "y";
hbar.Color = [0.99 1 1];
hbar.ScalebarLength = bar_length;
hbar.ConversionFactor = 640/4/1e3; % AscanLen / AscanRange / mm
hbar.UnitLabel = char("");

hbar2 = scalebar(ax, "x");
hbar2.Visible = "off";
hbar2.Axis = "x";
hbar2.Color = [0.99 1 1];
hbar2.ScalebarLength = bar_length;
hbar2.ConversionFactor = 4000/12/1e3; % BscanLen / Bscan Range / mm
hbar2.UnitLabel = char("");

