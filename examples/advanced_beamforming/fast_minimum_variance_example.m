%% Fast Minimum Variance Beamformer — Validation Example
%
% This example compares the new vectorized fast_minimum_variance implementation
% against the reference capon_minimum_variance, using the PICMUS CPWC dataset
% (same as the Generalized Beamformer publication).
%
% _by Ole Marius Hoel Rindal and Cursor Agent_

clear all; close all;

%% Download and load channel data
url = tools.zenodo_dataset_files_base();
filename = 'PICMUS_numerical_calib_v2.uff';
tools.download(filename, url, data_path);

channel_data = uff.channel_data();
channel_data.read([data_path, filesep, filename], '/channel_data');

fprintf('Dataset: %s\n', filename);
fprintf('N_channels = %d, N_waves = %d\n', channel_data.N_channels, channel_data.N_waves);

%% Define scan (same as Generalized Beamformer)
scan = uff.linear_scan('x_axis', linspace(-20e-3, 20e-3, 128).', ...
                       'z_axis', linspace(10e-3, 50e-3, 128).');

%% DAS beamforming (dimension.none — keep channels for MV)
das = midprocess.das();
das.channel_data = channel_data;
das.scan = scan;
das.dimension = dimension.none;
das.transmit_apodization.window = uff.window.tukey50;
das.transmit_apodization.f_number = 1.75;
das.receive_apodization.window = uff.window.tukey25;
das.receive_apodization.f_number = 1.75;

fprintf('\nRunning DAS (dimension.none)...\n');
tic;
b_data_delayed = das.go();
t_das = toc;
fprintf('DAS completed in %.2f seconds.\n', t_das);

%% MV parameters (matching Generalized Beamformer publication)
L_elements = floor(channel_data.probe.N_elements / 2);
K_in_lambda = 1.5;
regCoef = 1/100;

%% Reference: capon_minimum_variance
fprintf('\n=== Reference: capon_minimum_variance ===\n');
mv_ref = postprocess.capon_minimum_variance();
mv_ref.dimension = dimension.receive;
mv_ref.channel_data = channel_data;
mv_ref.scan = scan;
mv_ref.K_in_lambda = K_in_lambda;
mv_ref.L_elements = L_elements;
mv_ref.regCoef = regCoef;
mv_ref.doForwardBackward = 1;
mv_ref.receive_apodization = das.receive_apodization;
mv_ref.transmit_apodization = das.transmit_apodization;
mv_ref.input = b_data_delayed;

tic;
b_data_ref = mv_ref.go();
t_ref = toc;
fprintf('Reference MV completed in %.2f seconds.\n', t_ref);

%% New: fast_minimum_variance
fprintf('\n=== New: fast_minimum_variance ===\n');
mv_fast = postprocess.fast_minimum_variance();
mv_fast.dimension = dimension.receive;
mv_fast.channel_data = channel_data;
mv_fast.scan = scan;
mv_fast.K_in_lambda = K_in_lambda;
mv_fast.L_elements = L_elements;
mv_fast.regCoef = regCoef;
mv_fast.doForwardBackward = 1;
mv_fast.receive_apodization = das.receive_apodization;
mv_fast.transmit_apodization = das.transmit_apodization;
mv_fast.input = b_data_delayed;

tic;
b_data_fast = mv_fast.go();
t_fast = toc;
fprintf('Fast MV completed in %.2f seconds.\n', t_fast);

%% Compare results
fprintf('\n=== Comparison ===\n');
fprintf('Reference time: %.2f s\n', t_ref);
fprintf('Fast time:      %.2f s\n', t_fast);
fprintf('Speedup:        %.1fx\n', t_ref / t_fast);

% Numerical comparison
ref_data = b_data_ref.data(:);
fast_data = b_data_fast.data(:);

% Relative error (excluding near-zero pixels)
mask = abs(ref_data) > max(abs(ref_data)) * 1e-6;
rel_err = norm(fast_data(mask) - ref_data(mask)) / norm(ref_data(mask));
fprintf('Relative error: %.2e\n', rel_err);

% Correlation
corr_val = abs(ref_data' * fast_data) / (norm(ref_data) * norm(fast_data));
fprintf('Correlation:    %.6f\n', corr_val);

if rel_err < 1e-10
    fprintf('\n*** PASS: Results are numerically identical. ***\n');
elseif rel_err < 1e-3
    fprintf('\n*** PASS: Results match within tolerance (rel_err < 1e-3). ***\n');
else
    fprintf('\n*** FAIL: Results differ significantly! ***\n');
end

%% Compound and display
cc_ref = postprocess.coherent_compounding();
cc_ref.input = b_data_ref;
img_ref = cc_ref.go();

cc_fast = postprocess.coherent_compounding();
cc_fast.input = b_data_fast;
img_fast = cc_fast.go();

figure;
subplot(1,3,1);
img_ref.plot(gca, 'Reference MV');
subplot(1,3,2);
img_fast.plot(gca, 'Fast MV');
subplot(1,3,3);
diff_img = uff.beamformed_data(img_ref);
diff_img.data = abs(img_fast.data) - abs(img_ref.data);
imagesc(scan.x_axis*1e3, scan.z_axis*1e3, ...
    reshape(db(abs(diff_img.data)), scan.N_z_axis, scan.N_x_axis));
colorbar; title('Difference (dB)');
xlabel('x [mm]'); ylabel('z [mm]');
