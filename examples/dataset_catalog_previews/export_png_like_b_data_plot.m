function export_png_like_b_data_plot(b_data, filepath, dynamic_range_db)
%EXPORT_PNG_LIKE_B_DATA_PLOT  Grayscale PNG matching b_data.plot(..., dr, 'log') scaling.
if nargin < 3 || isempty(dynamic_range_db)
    dynamic_range_db = 60;
end

% Match uff.beamformed_data.get_image('log'); collapse to one 2-D image (first tx / frame).
envelope = b_data.get_image('log');
while ndims(envelope) > 2
    envelope = envelope(:, :, 1);
end
max_value = 0;
min_value = -dynamic_range_db;
g = (envelope - min_value) / (max_value - min_value + eps);
g = max(0, min(1, g));
img = g;
m = max(size(img));
if m > 800 && exist('imresize', 'file') == 2 %#ok<EXIST>
    img = imresize(img, 800 / m);
elseif m > 800
    step = ceil(m / 800);
    img = img(1:step:end, 1:step:end);
end

rgb = uint8(255 * repmat(img, [1, 1, 3]));
[p, ~] = fileparts(filepath);
if ~isfolder(p)
    mkdir(p);
end
imwrite(rgb, filepath);
end
