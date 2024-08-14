N = 10;

directory = "\\skat\research team\DUrban\MatLab\Process B-Scans\Flattening Grape\20240408-163327-00";
f_real = get_average_frame(directory, N, 3);
f_real = f_real(end/2:end, 1:end/2);

directory = "\\skat\research team\DUrban\MatLab\Process B-Scans\Flattening Grape\20240408-164353-00";
f_harmonic = get_average_frame(directory, N, 8);
f_harmonic = f_harmonic(end/2:end, 1:end/2);

directory = "\\skat\research team\DUrban\MatLab\Process B-Scans\Flattening Grape\20240408-165438-00";
f_flat = get_average_frame(directory, N, 8);
f_flat = f_flat(end/2:end, 1:end/2);

%% Generate Figure

range = [min(logScaling(f_real(:))) max(logScaling(f_real(:)))];

x_extent = 200:1900;
px_per_mm_x = 2000/12/1e3; % BscanLen / Imaging Range / mm
px_per_mm_y = 640/4/1e3; % AscanLen / AscanRange / mm

aspect_ratio = px_per_mm_x / px_per_mm_y;

figure;
f = tiledlayout("flow","TileSpacing","none", "padding", "tight");

ax = nexttile;
imshow(logScaling(f_real(:, x_extent)), displayRange=range)
colormap turbo
% clim([0, 1])
daspect([1, 1/aspect_ratio, 1])
% a1 = annotation("textbox", [0.25 0.8 0.1 0.1], "String", "A", "LineStyle","none");
text(ax, 10, 50, "a", "FontSize", 20, "Color", [0.99 0.99 0.99])

ax = nexttile;
imshow(logScaling(f_harmonic(:, x_extent), f_real), displayRange=range)
colormap turbo
% clim([0, 1])
daspect([1, 1/aspect_ratio, 1])
text(ax, 10, 50, "b", "FontSize", 20, "Color", [0.99 0.99 0.99])

ax = nexttile;
imshow(logScaling(f_flat(:, x_extent), f_real), displayRange=range)
colormap turbo
% clim([0, 1])
daspect([1, 1/aspect_ratio, 1])
% a1 = annotation("textbox", [0.25 0.2 0.1 0.1], "String", "C", "LineStyle","none");
text(ax, 10, 50, "c", "FontSize", 20, "Color", [0.99 0.99 0.99])

c = colorbar;
c.Label.String = "";
c.Limits = [0, 1];
c.Layout.Tile = "east";
c.Ticks = [0, 1];

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
hbar2.ConversionFactor = 2000/12/1e3; % BscanLen / Imaging Range / mm
hbar2.UnitLabel = char("");


% exportgraphics(f, "Figures\GrapeFlattening.png", "Resolution", 600)
