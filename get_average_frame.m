function f_mean = get_average_frame(directory, N, I, disp_comp)

if nargin < 3
    I = 1;
end

if nargin < 4
    disp_comp = true;
end

f = Frame(directory, 3);

if disp_comp
    f_phs = f.disp_comp(I);
end

i = 1;
n = 1;

while true

    f_n = Frame(directory, i);
    
    if disp_comp
        f_n = f_n.apply_phase(f_phs.phase);
    end

    if i == 1
        f_mean = abs(f_n.bscan_real).^2;
    else
        [reg, count] = register_bscans(f_n.bscan_real, f_mean);
        if count > 85
            f_mean = f_mean + abs(reg).^2;
            n = n + 1;
        end
    end

    if n >= N
        break
    end
    
    i = i + 1;
    


end
end