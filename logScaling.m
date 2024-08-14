function image = logScaling(imageIn, ref_img)
%     image = imageIn;
    image = 20 * log10(abs(imageIn));
    if nargin < 1
        ref_img = 20 * log10(abs(imageIn));
    end
%     image = image(30:end-30,:);

    if nargin < 2
        immax = max(image,[],'all');
        noiseFloor = mean(mean(image));
    else
        ref_img = 20 * log10(abs(ref_img));
        immax = max(ref_img(:));
        noiseFloor = mean(ref_img(:));
    end


    immin = 0;
    image = image - noiseFloor;
    image(image < 0) = 0;
    image = (image - immin) / (immax - immin);
end