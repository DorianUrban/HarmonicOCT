folders = dir("\\skat\Research Team\DUrban\MatLab\Process B-Scans\PMLongRange\Grape Set 1");
folder = folders.folder;

folders = folders(3:end);
names = {folders(:).name};

% Get background - last file in folder
directory = strcat(folder, "\",  cell2mat(names(end)));
f = Frame(directory);

% f = f.select_roi();
% pos = f.roi;

pos = [922 31 198 51];
noise = f.noise;

names = names(1:end-1);

results = {};
N = 10;

for ii = 1:length(names)
    snrs = zeros(1, N);

        % mean_bscan_corr = zeros(f.siz);
        mean_bscan = zeros(f.siz);
        mean_bscan_unregistered = zeros(f.siz);
        mean_bscan_corr = zeros(f.siz);
    
    for jj = 1:N
        directory = strcat(folder, "\",  cell2mat(names(ii)));
        f = Frame(directory, jj);
    
        if jj == 1
            f_corr = f.disp_comp(4);
            phs = f_corr.phase;
        else
            f_corr = f.apply_phase(phs);

        end

        mean_bscan = mean_bscan + abs(f_corr.bscan_original);
        mean_bscan_corr = mean_bscan_corr + abs(f_corr.bscan_real);
            
        f.noise = noise;
    
        f = f.get_psnr(pos);
        snrs(jj) = f.snr;
        % f.show();
        % title(strcat(f.info.comment, ", SNR = ", num2str(f.snr), " j= ", num2str(jj)))
       
        % sprintf("%i, %i", ii, jj)
        % f.comment

    end

    results{ii, 1} = f.comment;
    results{ii, 2} = str2double((f.comment));
    results{ii, 3} = mean(snrs);
    results{ii, 4} = std(snrs);
    results{ii, 5} = (f.info.reference_position);
    results{ii, 6} = fftshift(mean_bscan, 2);
    results{ii, 7} = mean_bscan_corr;
    % results{ii, 8} = mean_bscan_unregistered;

end

%%

idx = circshift(linspace(1, size(results, 1), size(results, 1)), 1);

results = sortrows(results, 2);
results_rearranged = results(idx, :);
% results_rmv = results(1:end, :);

freqs = cell2mat(results_rearranged(:, 5));
opd = 1.6e-3 .* freqs - 6.55;

snrs = cell2mat(results_rearranged(:, 3));
p = polyfit(opd(2:end), snrs(2:end), 1);
snrs_fit = polyval(p, opd(2:end));


img_ref = results{1, 7};
img_ref = img_ref(end/2:end, end/2:end);
maxx = max(max(logScaling(img_ref)));
minn = min(min(logScaling(img_ref)));
range = [minn, maxx];

x_extent = 200:1900;
px_per_mm_x = 2000/12/1e3; % BscanLen / Imaging Range / mm
px_per_mm_y = 640/4/1e3; % AscanLen / AscanRange / mm
aspect_ratio = px_per_mm_x / px_per_mm_y;

for jj = [6, 7]
    
    f = figure;
    tl = tiledlayout(8, 2, TileSpacing="tight", Padding="tight");
    
    axs = [];
    
    for ii = 1:size(results, 1)

        axs = [axs, nexttile];

        if ii == 1
            img = results_rearranged{ii, jj};
            imshow(logScaling(img(end/2:end, end/2:end), img_ref), Colormap=turbo, DisplayRange=[]);
            % title(sprintf("\\Deltaz %0.1f mm (real image)", opd(ii)))
            text(axs(ii), 1, 50, sprintf("z_s - z_r = %0.1f mm (conv. image)", opd(ii)), "Color", [1 1 1], "FontSize", 10)
        else
            img = results_rearranged{ii, jj};
            imshow(logScaling(img(end/2:end, end/2:end), img_ref),Colormap=turbo, DisplayRange=range);
            % title(sprintf("\\Deltaz %0.1f mm", opd(ii)))
            text(axs(ii), 1, 50, sprintf("z_s - z_r = %0.1f mm", opd(ii)), "Color", [1 1 1], "FontSize", 10)
        end
    
        daspect([1 1/aspect_ratio 1])
    end

    cbar = colorbar;
    cbar.Label.String = "";
    cbar.Limits = [0, 1];
    cbar.Ticks = [0, 1];
    cbar.Layout.Tile = "east";

    linkaxes(axs, "xy")
    
    bar_length = 1000;
    
    hbar = scalebar(axs(end), "y");
    hbar.Visible = "off";
    hbar.Axis = "y";
    hbar.Color = [0.99 1 1];
    hbar.ScalebarLength = bar_length;
    hbar.ConversionFactor = 640/4/1e3; % AscanLen / AscanRange / mm
    hbar.UnitLabel = char("");
    
    hbar2 = scalebar(axs(end), "x");
    hbar2.Visible = "off";
    hbar2.Axis = "x";
    hbar2.Color = [0.99 1 1];
    hbar2.ScalebarLength = bar_length;
    hbar2.ConversionFactor = 2000/12/1e3; % BscanLen / Bscan Range / mm
    hbar2.UnitLabel = char("");

    ax = nexttile([2, 2]);
    plot(opd(2:end), snrs(2:end), ".-")
    hold on
    plot(opd(1), snrs(1), "*")
    plot(opd(1:2), snrs(1:2), "--", "Color", [0.8500 0.3250 0.0980])

    yline(max(snrs(2:end))-3, "--")
    legend("1^{st} Harmonic", "Conventional Image", "", "-3dB")
    % errorbar(opd, cell2mat(results_rearranged(:, 3)), cell2mat(results_rearranged(:, 4)))
    xlabel("z_s - z_r (mm)")
    ylabel("SNR (dB)")
    ax.FontSize = ax.FontSize + 2;
    

    % if jj == 6
    %     exportgraphics(f, "Figures\GrapeLongRangeDispersed.png", "Resolution", 600)
    % else
    %     exportgraphics(f, "Figures\GrapeLongRange.png", "Resolution", 600)
    % end
end






