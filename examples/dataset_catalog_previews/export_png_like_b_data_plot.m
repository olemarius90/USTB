function export_png_like_b_data_plot(b_data, filepath, dynamic_range_db)
%EXPORT_PNG_LIKE_B_DATA_PLOT  Grayscale PNG matching b_data.plot(..., dr, 'log') scaling.
if nargin < 3 || isempty(dynamic_range_db)
    dynamic_range_db = 60;
end

data = b_data.data;
if size(data, 4) > 1
    data = data(:, :, :, 1);
end
envelope = abs(data);
envelope = 20 * log10(envelope ./ max(envelope(:)) + realmin('single'));
max_value = 0;
min_value = -dynamic_range_db;
g = (envelope - min_value) / (max_value - min_value + eps);
g = max(0, min(1, g));

switch class(b_data.scan)
    case 'uff.linear_scan'
        img = reshape(g, [b_data.scan.N_z_axis, b_data.scan.N_x_axis, size(g, 3)]);
    case 'uff.sector_scan'
        img = reshape(g, [b_data.scan.N_depth_axis, b_data.scan.N_azimuth_axis, size(g, 3)]);
    otherwise
        error('export_png: unsupported scan %s', class(b_data.scan));
end

img = squeeze(img(:, :, 1));
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
