classdef Frame
    properties
        frame_original
        frame
        phase
        ascan_no
        roi
        mod_freq
        subdir
        bscan_original
        noise
        siz
        bscan_real
        stft_spectum
        window
        stft_ref
        snr
        tukey_win
        info
        comment
        mean_bscan
        mean_constant
        focus_phase
        w0
        def_locs
        def_dists
        bscan_refocused
        z_component

    end

    methods

        function obj = Frame(directory, frame_no)

            if nargin < 2
                frame_no = 1;
            end

            [frame_in, obj.info] = get_frame(directory, frame_no, false);

            if size(frame_in, 1) == 768
                frame_in = frame_in(1:640, :);
            end

            obj.comment = obj.info.comment;

            obj.mean_bscan = mean(frame_in, 2);
            rolling_mean = movmean(frame_in - obj.mean_bscan, 300, 1);
            obj.mean_constant = mean(frame_in(:));
            obj.window = hanning(size(frame_in, 1));

            obj.tukey_win = tukeywin(size(frame_in, 1), 0.05) + eps;
            obj.frame = fft(ifftshift(fftshift(ifft(frame_in - rolling_mean  - obj.mean_bscan)).* obj.tukey_win));
            obj.frame_original = fft(ifftshift(fftshift(ifft(frame_in -rolling_mean - obj.mean_bscan)).* obj.tukey_win));
            obj.phase = ones(size(obj.frame, 1), 1);
            obj.bscan_original = fftshift(ifft(obj.frame .* obj.window), 1);
            obj.bscan_real = obj.bscan_original;
            [~, obj.ascan_no] = max(sum(abs(obj.bscan_original(:, 1:end/2)).^3, 1));
            obj.siz = [size(obj.bscan_original, 1), size(obj.bscan_original, 2)];

            obj.w0 = 6.2e-6;

%             obj.noise = [  -0.0000    6.0055    0.1916  143.6159
%                            -0.0000    2.6989    0.0030    3.1231];           
        
        end

        function obj = apply_phase(obj, phase)
            % [obj.frame, phase_out] = apply_coefficients(obj.frame, coeffs);

            obj.bscan_real = fftshift(ifft((obj.frame .* conj(phase) .* obj.window)));
            obj.phase = obj.phase .* phase;
        end

        function obj = show(obj, show_phase)

            % bscan = fftshift(ifft(obj.frame .* obj.window));
            bscan = obj.bscan_real;

            figure

            if isreal(obj.frame_original)
                imshowpair(logScaling(obj.bscan_original(:, 1:end/2)), logScaling(bscan(:, 1:end/2)), "method", "montage")
            else
                imshowpair(logScaling(obj.bscan_original(end/2:end, end/2:end)), logScaling(bscan(end/2:end, end/2:end)), "method", "montage")
            end

            if nargin >= 2
                if show_phase
                    figure
                    plot(unwrap(angle(obj.phase)));
                end
            end

        end

        function obj = optimise(obj)

            if isempty(obj.y_roi)
                obj = obj.select_roi;
            end

%             [frame_corrected, coeffs] = brute_NDC(obj.frame_original, obj.ascan_no, false, obj.y_roi);

            fringe = obj.frame(:, obj.ascan_no);

            orders = 2:8;
            coeffs_initial = [0 -300, 0, -500, -9.2e4, 0, 0];
%             coeffs_initial = zeros(length(orders), 1);
%             coeffs_initial = coeffs;
            coeffs = zeros(length(orders), 1);

            for i = 1:length(orders)
                [coeffs(i), ~] = finder(fringe, orders(i), coeffs_initial(orders(i) - 1), obj.y_roi);
                fringe = apply_coefficients(fringe, coeffs(i));
            end

            obj = obj.apply_phase(coeffs);
        end

        function obj = disp_comp(obj, I)


            if nargin > 1
                ridge = tfphase(obj.frame, I);
                
            else
                ridge = tfphase(obj.frame);
            end

            ridge_i = interp1(1:length(ridge), ridge, linspace(1, length(ridge), size(obj.frame, 1)), "pchip");

            % figure
            % plot(ridge_i)
            
            % obj.phase = obj.phase .* exp(-1i * pi .* ridge_i).';
            phs = exp(-1i * pi .* ridge_i).';
            obj = obj.apply_phase(phs);

            % frame_corr = obj.frame .* conj(obj.phase);
            % obj.frame = frame_corr;
            % 
            % obj.bscan_real = fftshift(ifft((frame_corr .* obj.window)), 1);

            % figure;
            % plot(unwrap(angle(obj.phase)))

        end
 
        function obj = select_roi(obj)
            f1 = figure;
            
            % if isreal(obj.frame_original)
                imshow(logScaling(obj.bscan_real(:, 1:end/2)));
            % else
                % imshow(logScaling(obj.bscan_real(1:end/2, 1:end/2)));
            % end
            
            % xline(obj.ascan_no, "r-.")
            title("Select Region of Interest")

            % if ~isempty(obj.y_roi)
            %     yline(obj.y_roi(1), "c-.")
            %     yline(obj.y_roi(2), "c-.")
            % end

            obj.roi = round(drawrectangle().Position);
            pos = obj.roi;
            
            bscan_cropped = obj.bscan_real(pos(2):pos(2)+pos(4) - 1, pos(1):pos(1)+pos(3) - 1);
            obj.noise = var(abs(bscan_cropped(:)));

            close(f1)
        end

        % function obj = estimate_noise(obj)
        % 
        %     f = figure;
        % 
        %     if isreal(obj.frame_original)
        %         imshow(logScaling(obj.bscan_original(:, :)));
        %     else
        %         imshow(logScaling(obj.bscan_original(1:end/2, :)));
        %     end
        % 
        %     title("Select Noise")
        % 
        %     roi = drawrectangle();
        %     pos = round(roi.Position);
        % 
        %     close(f)
        % 
        %     bscan_cropped = obj.bscan_original(pos(2):pos(2)+pos(4) - 1, pos(1):pos(1)+pos(3) - 1);
        % 
        %     m = mean(mean(real(bscan_cropped)));
        %     s = std(real(bscan_cropped), 0, "all");
        %     skew = skewness(real(bscan_cropped), 1, "all");
        %     kurt = kurtosis(real(bscan_cropped), 1, "all");
        % 
        %     obj.noise(1, :) = [m, s, skew, kurt];
        % 
        %     m = mean(mean(imag(bscan_cropped)));
        %     s = std(imag(bscan_cropped), 0, "all");
        %     skew = skewness(imag(bscan_cropped), 1, "all");
        %     kurt = kurtosis(imag(bscan_cropped), 1, "all");
        % 
        %     obj.noise(2, :) = [m, s, skew, kurt];
        % 
        % end

        function obj = get_snr(obj)

            f = figure;

            if isreal(obj.frame_original)
                imshow(logScaling(obj.bscan_real(:, :)));
            else
                imshow(logScaling(obj.bscan_real(1:end/2, :)));
            end

            title("Select Noise")
            
            roi = drawrectangle();
            pos = round(roi.Position);
            
            close(f)

            bscan_cropped = obj.bscan_real(pos(2):pos(2)+pos(4) - 1, pos(1):pos(1)+pos(3) - 1);

            s = std(abs(bscan_cropped(:)));
%             m = mean(imag(bscan_cropped(:)));

            %

            f = figure;

            if isreal(obj.frame_original)
                imshow(logScaling(obj.bscan_real(:, :)));
            else
                imshow(logScaling(obj.bscan_real(1:end/2, :)));
            end

            title("Select reflection")
            
            roi = drawrectangle();
            pos = round(roi.Position);
            
            close(f)

            bscan_cropped = obj.bscan_real(pos(2):pos(2)+pos(4) - 1, pos(1):pos(1)+pos(3) - 1);

            ma = median(max(abs(bscan_cropped)));
%             ma = mean(max(abs(obj.bscan_real)));

            ma = prctile(max(abs(obj.bscan_real).^2, [], 1), 99);

        end

        function obj = get_psnr(obj, roi)

            if nargin < 2
                f = figure;

                imshow(logScaling(obj.bscan_real(:, :)));

    
                title("Select Noise")
                
                roi = drawrectangle();
                pos = round(roi.Position);
                
                close(f)
    
                bscan_cropped = obj.bscan_real(pos(2):pos(2)+pos(4) - 1, pos(1):pos(1)+pos(3) - 1);
                s = var(abs(bscan_cropped(:)));
            
            else

                if isempty(obj.noise)
                    bscan_cropped = obj.bscan_real(roi(2):roi(2)+roi(4) - 1, roi(1):roi(1)+roi(3) - 1);
                    s = var(abs(bscan_cropped(:)));
                else
                    s = obj.noise;
                end

            end
%             m = mean(imag(bscan_cropped(:)));

            % ma = median(max(abs(obj.bscan_real).^2, [], 1));
            % ma = prctile(max(abs(obj.bscan_real).^2, [], 1), 95);
            ma = prctile(max(abs(obj.bscan_real(1:end/3, :)).^2, [], 1), 95);

            obj.snr = 10 * log10(ma / s);

        end


        function obj = downsample_laterally(obj, factor)

            obj.frame = downsample(obj.frame.', factor).';
    
            obj.mean_bscan = mean(obj.frame, 2);
            rolling_mean = movmean(obj.frame - obj.mean_bscan, 300, 1);


                
            % obj.mean_bscan = mean(obj.frame, 2);
            obj.mean_constant = mean(obj.frame(:));
            % rolling_mean = movmean(obj.frame, 10);
            obj.window = hanning(size(obj.frame, 1));

            obj.tukey_win = tukeywin(size(obj.frame, 1), 0.05) + eps;
            obj.frame = fft(ifftshift(fftshift(ifft(obj.frame - obj.mean_bscan - rolling_mean)).* obj.tukey_win));
            obj.frame_original = fft(ifftshift(fftshift(ifft(obj.frame - obj.mean_bscan)).* obj.tukey_win));
            obj.phase = ones(size(obj.frame, 1), 1);
            obj.bscan_original = fftshift(ifft(obj.frame .* obj.window));
            obj.bscan_real = obj.bscan_original;
            [~, obj.ascan_no] = max(sum(abs(obj.bscan_original(:, 1:end/2)).^3, 1));
            obj.siz = [size(obj.bscan_original, 1), size(obj.bscan_original, 2)];  

        end

        function obj = coherent_average_laterally(obj, factor)
            obj.frame = squeeze(mean(reshape(obj.frame, obj.siz(1), factor, []), 2));

            obj.mean_bscan = mean(obj.frame, 2);
            rolling_mean = movmean(obj.frame - obj.mean_bscan, 300, 1);


            % obj.mean_bscan = mean(obj.frame, 2);
            obj.mean_constant = mean(obj.frame(:));
            % rolling_mean = movmean(obj.frame, 10);
            obj.window = hanning(size(obj.frame, 1));

            obj.tukey_win = tukeywin(size(obj.frame, 1), 0.05) + eps;
            obj.frame = fft(ifftshift(fftshift(ifft(obj.frame - obj.mean_bscan - rolling_mean)).* obj.tukey_win));
            obj.frame_original = fft(ifftshift(fftshift(ifft(obj.frame - obj.mean_bscan)).* obj.tukey_win));
            obj.phase = ones(size(obj.frame, 1), 1);
            obj.bscan_original = fftshift(ifft(obj.frame .* obj.window));
            obj.bscan_real = obj.bscan_original;
            [~, obj.ascan_no] = max(sum(abs(obj.bscan_original(:, 1:end/2)).^3, 1));
            obj.siz = [size(obj.bscan_original, 1), size(obj.bscan_original, 2)];  
        end

        % function bscan_average = incoherent_average_laterally(obj, factor)
        %     bscan_average = squeeze(mean(reshape(abs(obj.bscan_real), obj.siz(1), factor, []), 2));
        % end

        function obj = incoherent_average_laterally(obj, factor)

            obj.bscan_real = squeeze(mean(reshape(abs(obj.bscan_real), obj.siz(1), factor, []), 2));
            obj.bscan_original = squeeze(mean(reshape(abs(obj.bscan_original), obj.siz(1), factor, []), 2));

            obj.frame = fft(ifftshift(obj.bscan_real));
            obj.frame_original = fft(ifftshift(obj.bscan_original));

            obj.mean_bscan = mean(obj.frame, 2);

            obj.window = hanning(size(obj.frame, 1));
            % 
            [~, obj.ascan_no] = max(sum(abs(obj.bscan_original(:, 1:end/2)).^3, 1));
            obj.siz = [size(obj.bscan_original, 1), size(obj.bscan_original, 2)];  
        end


        function obj = apply_noise_mask(obj)

            if isempty(obj.noise)
                obj = obj.estimate_noise();
            end

            f = figure;

            bscan = fftshift(ifft(obj.frame));
            bscan_show = fftshift(ifft(obj.frame .* obj.window));

            if isreal(obj.frame_original)
                imshow(logScaling(bscan_show(:, :)));
            else
                imshow(logScaling(bscan_show(1:end/2, :)));
            end

            title("Select Real Image")

            roi = drawrectangle();
            pos = round(roi.Position);

            close(f)

            bscan_cropped = bscan(pos(2):pos(2)+pos(4) - 1, pos(1):pos(1) + pos(3) - 1);

            gen_noise = pearsrnd(obj.noise(1, 1), obj.noise(1, 2), obj.noise(1, 3), obj.noise(1, 4), obj.siz(1), obj.siz(2));
            gen_noise_i = pearsrnd(obj.noise(2, 1), obj.noise(2, 2), obj.noise(2, 3), obj.noise(2, 4), obj.siz(1), obj.siz(2));

            gen_noise = gen_noise + gen_noise_i;
            
            bscan(pos(2):pos(2)+pos(4) - 1, pos(1):pos(1) + pos(3) - 1) = gen_noise(pos(2):pos(2)+pos(4) - 1, pos(1):pos(1) + pos(3) - 1);

            gen_noise(pos(2) : pos(2) + pos(4) - 1, pos(1):pos(1) + pos(3) - 1) = bscan_cropped;

            obj.bscan_real = gen_noise;
            
            obj.frame = fft(fftshift(bscan));
            
        end

        function obj = stft_correction(obj)

            bscan = fftshift(ifft(obj.frame));

            % Small offset is added to window to avoid 0s in the image
            tukey_win = tukeywin(size(obj.frame, 1), 0.05) + eps;
            
            bscan_filtered = bscan .* tukey_win;
            
            obj.frame = fft(fftshift(bscan_filtered));
        
            window = hann(128, "periodic");

            fringe = obj.frame(:, floor(end/4) + 1);

            [sst, ~] = stft(fringe, [], "window", window, "OverlapLength", length(window) - 1, "FFTLength", 2 * length(window));

            sst_half = sst(1:end/2, :);
            template = sst_half(:, floor(end/2) - 100);

            corr = normxcorr2(abs(template), abs(sst_half));
            [fridge, ~] = tfridge(corr, 1:size(corr, 1), .1, "NumRidges", 1);

%             figure
%             imagesc(abs(sst))
%             set(gca, "YDir", "normal")
            
            figure
            hold on
            imagesc(abs(corr))
            plot(fridge, "r")
            set(gca, "YDir", "normal")

            fridge_i = interp1(1:length(fridge), fridge, linspace(1, length(fridge), obj.siz(1)), "pchip");
            fridge_p = cumsum(fridge_i - mean(fridge_i));

            scale = 2 / obj.siz(1);
            disp_phase = conj(exp(2i * pi .* fridge_p * scale));

            obj.frame = obj.frame.' .* conj(disp_phase);obj.frame = obj.frame.';
            obj.phase = obj.phase .* conj(disp_phase);

            obj.bscan_real = fftshift(ifft(obj.frame .* obj.window));

            obj.stft_spectrum = abs(sst);
        end

        function obj = stft_manual(obj)

            [sst, f] = fsst(obj.frame(:, 220), 600e6, hann(256, "periodic"), "yaxis");
            obj.stft_spectrum = abs(sst);
            
            figure
            hold on
            x = 1:size(obj.stft_spectrum,2);
            y = linspace(-250e6, 250e6, size(obj.stft_spectrum, 1));

            pcolor(x, y, obj.stft_spectrum)
            shading flat
%             ylim([0, y(end)])

            roi = drawpolyline();

            p = interp1(roi.Position(:, 1), roi.Position(:, 2), x, "pchip");

            plot(x, p, "r.-");

            p_s = cumsum(p - mean(p));

            x = x.';

            x_new = interp1(p_s, x, linspace(p_s(1), p_s(end), length(obj.frame_original)), "pchip");

            x = 1:length(obj.frame_original);

%             obj.frame = interp1(x, obj.frame_original, x_new, "pchip");
            
            phase_corr = exp(-2j * pi * p_s / 600e6);
            obj.frame = obj.frame_original .* phase_corr.';

%             figure
%             pcolor(x_new, y, obj.stft_spectrum)
%             shading flat
%             ylim([0, y(end)])

            
            
        
        end

        function obj = stft_xcorr(obj)
            
            sizz = size(obj.stft_spectrum);
            flat_spectrum = obj.stft_spectrum(:);

            kernel = abs(obj.bscan_real(:, obj.ascan_no)); % This needs to be the same as fringe used to calculate stft
            kernel = obj.stft_ref;

            corr = normxcorr2(kernel, flat_spectrum);

            corr = corr(1:obj.siz(1)^2);
            
            [m, I] = max(corr);

            figure
            plot(corr(I:I+4096))

            [r, c] = ind2sub(size(obj.stft_spectrum), I);

            corr = reshape(corr, sizz(1), sizz(2));

            [~, locs] = max(corr(1:end/2, :), [], 1);
            size(locs)

            figure
            imshow(corr);
            colorbar
            hold on
            plot(locs, "r")

        end

        function obj = fsst(obj)

            win = hann(256, "periodic");
            method = "linear";

            fs = 600e6;

            fringe_p = obj.frame(:, 212);

            [sst, f] = fsst(fringe_p, fs, win, "yaxis");

            z = img_real_distance(4096).';

            z = interp1(z, z, linspace(min(z), max(z), size(sst, 1)));

            z = f;

            t = (-length(fringe_p) / 2) : (length(fringe_p) / 2) - 1;
            t = t.';
            
            sst(end, :) = min(sst(:));
            
            mx = max(abs(sst(:)))*ones(size(t));
            
%             sst_filt = imgaussfilt(abs(sst), 3);
%             sst_filt = sst_filt - mean(sst_filt(:));
%             sst_filt(end, :) = min(sst_filt(:));

            sst_filt = sst;
            [fridge, iridge] = tfridge(sst_filt, f, 3, NumRidges=2, NumFrequencyBins=16);      
            
            figure;
            subplot(2, 1, 1)
            pcolor(t, f, abs(sst))
            title(method)
            shading flat
            hold on
            plot(t, fridge, "--")
            ylim([min(z) max(z)])

%             roi = drawpolyline();
            
%               #############
            ridge_s = cumsum(fridge(:, 2));
            t_new = interp1(ridge_s, t, linspace(ridge_s(1), ridge_s(end), length(ridge_s)), method);

            fridge_new = interp1(t, fridge(:, 2), t_new, method);

            plot(t, fridge_new', "g--")
%               #############

            obj.frame = interp1(t, obj.frame, t_new, method);

            subplot(2,1,2)
            [sst_new, f] = fsst(obj.frame(:, 124), fs, win, "yaxis");
            pcolor(t, f, abs(sst_new))
            shading flat
            ylim([0 max(f)])


%             figure
%             hold on
%             plot(t, ridge_s)
%             plot(t, cumsum(fridge_new))
            
        end     

        function obj = refocus(obj, z_in)
            % 
            % if nargin <= 2
            %     w0 = 6.2e-6;
            % end
    
            % Fourier transform of image along x
            imgfft = fft(obj.bscan_real(:, 1:end/2), [], 2);
            
            % system properties
            N = obj.siz(1);
            M = obj.siz(2) / 2;
            lambda_0 = 1060e-9;
            dlambda = 95e-9;
            mag = 1.33;
            n0 = 1.45;
            Nz = 1:N;
            % w0 = 6.2e-6;            % Mode field diameter
        
            % 
            % Calculating system resolution
            axial_res = (2 * log2(2) / pi) * (lambda_0^2 / dlambda);    % m/px
            
            z_max = (lambda_0^2 * 2 * N / (4 * dlambda));
            
            lateral_res = sqrt(2 * log2(2)) * obj.w0;   % Need to add M here?
        
        %     lateral_res = 4.45e-6 * scale_factor;
            
            % System spatial frequencies
            
            axial_freq = 1 / axial_res;
            lateral_freq = 1 / lateral_res;
        
            % Find frequency increment
            dz = 1 / (axial_res * N);
            dx = 1 / (lateral_res * M);
            
            % Spatial frequency axes
            z_axis = -axial_freq / 2 : dz : (axial_freq / 2) - dz;
            x_axis = -lateral_freq / 2 : dx : (lateral_freq / 2) - dx;
            
            % Converting to angular wavenumber
            z_axis = z_axis * 2 * pi;
            x_axis = x_axis * 2 * pi;
            
            z = Nz .* lambda_0^2 / (2 * n0 * dlambda);z = z.';
            % z = defocus; % + z
            z = z_in;
            zmat = repmat(z, 1, M);
            
            % obj.focus_phase = exp(1i .* zmat .* lambda_0 * mag^2 .* x_axis.^2 / (4 * pi * n0));
            
            % obj.focus_phase = exp(2i * pi * zmat.' .* x_axis.^2 / 1e9);
            obj.focus_phase = exp(1i * zmat .* x_axis.^2 * lambda_0 * mag^2 / (4 * pi * n0));
            % obj.focus_phase = exp(1i * zmat.' .* x_axis.^2 / 1e8);
    
            % refocus = exp(1i * pi * lambda_0 .* zmat .* x_axis.^2 * scalefactor / 2);
            % 
            % figure
            % imagesc(angle(obj.focus_phase))
    
            imgfft2 = fft(obj.bscan_real(:, end/2+1:end), [], 2);
        %     
            img_ref_l = ifft(imgfft .* fftshift(obj.focus_phase), [], 2);
            img_ref_r = ifft(imgfft2 .* fftshift(obj.focus_phase), [], 2);
    
            obj.bscan_refocused = [img_ref_l img_ref_r];
            % obj.refocus = refocus;
    
        end

        function obj = find_focus(obj)

            if isempty(obj.roi)
                obj = obj.select_roi;
            end

            section = logScaling((obj.bscan_real(obj.roi(2):obj.roi(2) +obj.roi(4), obj.roi(1):1:obj.roi(1)+obj.roi(3))));  
    
            z_peaks = sum(section, 2);
            maxx = max(z_peaks);
            [pks, locs] = findpeaks(z_peaks, "MinPeakProminence", 0.1 * maxx, "NPeaks", 8);
            % locs = locs + floor(obj.roi(2));

            cmapp = parula(length(locs));
            
            figure;
            imshow(section)
            axis on
            colorbar
            
            hold on
            
            for ii = 1:length(locs)
                yline(locs(ii), "Color", cmapp(ii, :), "LineWidth", 1)
            end

            optim_l = 50;
            focus_distance = linspace(-5e-3, 5e-3, optim_l);

            data = zeros(1, length(locs));

            for jj = 1:length(locs)

                cubes = zeros(1, optim_l);
                maxs = zeros(1, optim_l);
                widths = zeros(1, optim_l);

                for ii = 1:optim_l
    
                    obj_temp = obj.refocus(focus_distance(ii));
                    section = abs((obj_temp.bscan_refocused(obj_temp.roi(2):obj_temp.roi(2) +obj_temp.roi(4), obj_temp.roi(1):1:obj_temp.roi(1)+obj_temp.roi(3)))).^3;
                    profile = section(locs(jj), :);
                    
                    cubes(ii) = sum(profile);
                    maxs(ii) = max(profile);
    
                    [~, ~, w] = findpeaks(profile, "MinPeakProminence", 0.1 * maxs(ii), "WidthReference","halfheight", "Annotate","extents");
                    widths(ii) = mean(w);
                     
                    % if mod(ii, 2) == 0
                    % figure
                    % findpeaks(profile, "MinPeakProminence", 0.1 * maxs(ii), "WidthReference","halfheight", "Annotate","extents");
                    % title(sprintf("Z pos: %f", z_range(ii)))
                    % end
                    
    
                end
    
                % figure
                % % yyaxis left
                % plot(1000 * focus_distance, normalize(cubes, "range", [0 1]));
                % hold on
                % plot(1000 * focus_distance, normalize(maxs, "range", [0 1]));
                % % hold on
                % % yyaxis right
                % plot(1000 * focus_distance, normalize(widths, "range", [0 1]));
                % legend("Cubes", "Max", "Width")

                % [~, I] = min(widths);
                [~, I] = max(cubes);
                data(jj) = focus_distance(I);


    
            end

            figure
            plot(locs + obj.roi(2), 1e3 * data, "o-")

            obj.def_locs = locs;
            obj.def_dists = data;
            N = obj.siz(1);
            M = obj.siz(2) / 2;
            Nz = 1:N;
            lateral_res = sqrt(2 * log2(2)) * obj.w0;   % Need to add M here?
            lateral_freq = 1 / lateral_res;
            dx = 1 / (lateral_res * M);
            x_axis = (-lateral_freq / 2 : dx : (lateral_freq / 2) - dx) * 2 * pi;

            lambda_0 = 1060e-9;
            dlambda = 95e-9;
            mag = 1.33;
            n0 = 1.4;

            z = Nz .* lambda_0^2 / (2 * n0 * dlambda);
            % z = z * lambda_0 * mag^2 / (4 * pi * n0);
            
            figure;
            plot(locs+obj.roi(2), 1000 * data, "o-")
            hold on
            plot(1:N, -1000*z)
            obj.z_component = polyfit(locs+obj.roi(2), data, 1);
            focus_distance = polyval(obj.z_component, 1:N);
            hold on
            plot(1:N, 1000 * focus_distance);
            focus_distance = repmat(focus_distance, length(x_axis), 1);

            % obj.focus_phase = exp(2i * pi .* x_axis.^2 .* focus_distance.' / 1e9);
            obj.focus_phase = exp(1i * focus_distance.' .* x_axis.^2 * lambda_0 * mag^2 / (4 * pi * n0));
            % obj.focus_phase = exp(1i * focus_distance.' .* x_axis.^2 / 1e8);


            imgfft = fft(obj.bscan_real(:, 1:end/2), [], 2);
            imgfft2 = fft(obj.bscan_real(:, end/2+1:end), [], 2);
        %     
            img_ref_l = ifft(imgfft .* fftshift(obj.focus_phase, 2), [], 2);
            img_ref_r = ifft(imgfft2 .* fftshift(obj.focus_phase, 2), [], 2);
    
            obj.bscan_refocused = [img_ref_l img_ref_r];


        end

        function obj = plus(obj1, obj2)

            obj = obj1;
            obj.bscan_real = obj1.bscan_real + obj2.bscan_real;
            obj.bscan_original = obj1.bscan_original + obj2.bscan_original;

            obj.frame = obj1.frame + obj2.frame;
            obj.frame_original = obj1.frame_original + obj2.frame_original;

        end

        function obj = incoh_plus(obj1, obj2)

            obj = obj1;
            obj.bscan_real = abs(obj1.bscan_real) + abs(obj2.bscan_real);
            obj.bscan_original = abs(obj1.bscan_original) + abs(obj2.bscan_original);

            obj.frame = fft(ifftshift(obj.bscan_real));
            obj.frame_original = fft(ifftshift(obj.bscan_original));

        end

        function obj = generate_stft_video(obj)

            downsample_factor = 5;
            window = hann(33, "symmetric");
            frame = obj.frame_original;

            fringe = frame(:, 1);
            wv = wvd(fringe, "smoothedPseudo", window, window, "MinThreshold", 0);

            M_wv = zeros([size(wv) size(frame, 2) / downsample_factor]);

            for ii = 1:size(frame, 2) / downsample_factor

                wv = wvd(hilbert(real(frame(:, ii * downsample_factor))), "smoothedPseudo", window, window, "MinThreshold", 0);
                M_wv(:, :, ii) = logScaling(abs(wv).^2+eps);

                (size(frame, 2) / downsample_factor) - ii
        
            end

            vid_wv = VideoWriter("videos\FlatGrape2.avi", "Uncompressed AVI");
            vid_wv.FrameRate = 15;
            
            open(vid_wv)
            
            for ii = 1:size(M_wv, 2)
            
                writeVideo(vid_wv, squeeze(M_wv(:, ii, :)));
            end
            
            close(vid_wv)

        end

    end    
end

