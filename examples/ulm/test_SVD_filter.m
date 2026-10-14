%% Set paths
local_path = fullfile(ustb_path(), 'data');

input_folder = fullfile(local_path, 'RF_channeldata');
beamformed_path_CF = fullfile(local_path, 'beamformed_slidingSVD_CF');
if ~exist(beamformed_path_CF, 'dir')
    mkdir(beamformed_path_CF);
end

% beamformed_path_DAS = 

% Creating file for generated figures
figure_path = fullfile(ustb_path(), 'examples', 'ulm', 'Figures');
if ~exist(figure_path, 'dir')
    mkdir(figure_path);
end

% Find all channeldata uff files
files = dir(fullfile(input_folder, '*.uff'));

uffFiles = fullfile({files.folder}, {files.name});

nFiles = numel(uffFiles);

fprintf('Found %d UFF-files.\n', nFiles);

for i = 1:nFiles
    fprintf('%2d: %s\n', i, files(i).name);
end


%% Scan
invivo_scan = uff.read_object(fullfile(local_path, 'InVivoRatBrain_scan.uff'), '/scan');

%% Find lambda
invivo_ch_data = uff.read_object(uffFiles{1}, '/channel_data');
lambda = invivo_ch_data.lambda;
clear invivo_ch_data

%% Step 2: Set up beamforming
% Upscale scan resolution
invivo_scan = tools.scan_integer_upscale(invivo_scan, 2);

% DAS
das = midprocess.das();
das.dimension = dimension.transmit;
das.scan = invivo_scan;
das.receive_apodization.f_number = 1.4; 
das.receive_apodization.window = uff.window.hamming;

% CF
cf = postprocess.coherence_factor();
cf.dimension = dimension.receive;

% Sliding SVD settings
windowFiles = 2; % Red two files at a time
cutoff = 5; % Removes ranks 1:5
batch_size = 40;
frame_rate = 1000;

for targetFileIndex = 1:nFiles
    svd_ch_data = tools.slidingTemporalSVD(uffFiles, targetFileIndex, windowFiles, cutoff);

    total_frames = size(svd_ch_data.data, 4);
    
    % Hold all beamformed frames from this input file
    b_CF_file = [];
    % b_DAS_file = [];

    for start_frame = 1:batch_size:total_frames
        end_frame = min(start_frame + batch_size - 1, total_frames);
        fprintf('Processing frames %d to %d of %d...\n', start_frame, end_frame, total_frames);
        
        % Slice current batch of frames
        %invivo_ch_data.data = raw_frames(:,:,:, start_frame:end_frame);
        frame_idx = start_frame:end_frame;

        % Extract batch from SVD filtered file
        svd_ch_data_batch = uff.channel_data(svd_ch_data);
        svd_ch_data_batch.data = svd_ch_data.data(:,:,:,frame_idx);
    
        % DAS
        das.channel_data = svd_ch_data_batch;
        b_tx_batch = das.go();
    
        % CF
        cf.input = b_tx_batch;
        b_CF_batch = cf.go();
        b_CF_batch.frame_rate = frame_rate; % WHY from param.mat?
    
        % Coherent compounding (ADD)
        % das_rx.input = b_data_tx;
        % invivo_b_das = das_rx.go();
    
        % Append batches from this file 
        if isempty(b_CF_file)
            b_CF_file = b_CF_batch;
            %b_DAS_file = b_tx_batch;
        else
            b_CF_file.data = cat(4, b_CF_file.data, b_CF_batch.data);
            %b_DAS_file.data = cat(4, b_DAS_file.data, b_tx_batch.data);
        end
    
        clear svd_ch_data_batch b_tx_batch b_CF_batch;
    
    end
   
    output_file_CF = fullfile(beamformed_path_CF, ...
       sprintf('beamformed_CF_slidingSVD_%03d.uff', targetFileIndex));

    uff.write_object(output_file_CF, b_CF_file,'b_data');

    clear svd_ch_data b_CF_file

end

%% Combine beamformed files
files_CF = dir(fullfile(beamformed_path_CF, '*uff'));

for i = 1:numel(files_CF)
    filename = sprintf('beamformed_CF_slidingSVD_%03d.uff', i);
    fprintf('Loading beamformed file %d/%d:\n%s\n', i, numel(files), filename);

    b_data_part = uff.read_object(fullfile(beamformed_path_CF, filename), '/b_data');
    if i==1
        b_data_all_CF = b_data_part;
    else
        b_data_all_CF.data = cat(4, b_data_all_CF.data, b_data_part.data);
    end
    clear b_data_part
end


%% Saving B-mode images
% Save SVD filtered Coherence factor beamformed B-mode image
fig_CF = figure('Visible','off');
b_data_all_CF.plot(fig_CF, 'In Vivo Rat Brain sliding SVD - CF (5 files)', 60);
exportgraphics(gca, fullfile(figure_path, 'Beamformed_CF_sliding_SVD_5.png'));

% % Save SVD filtered DAS beamformed B-mode image
% fig_DAS = figure('Visible','off');
% invivo_b_DAS_full.plot(fig_DAS, 'In Vivo Rat Brain sliding SVD - DAS (5 files)', 60);
% exportgraphics(gca, fullfile(figure_path, 'Beamformed_DAS_sliding_SVD_5.png'));
% 
% % Save Maximum Intensity Projection (MIP) of SVD-filtered CF over all 40 frames
% fig_mip = figure('Visible', 'off');
% mip_data = uff.beamformed_data(invivo_b_CF_full);
% mip_data.data = max(abs(invivo_b_CF_full.data), [], 4);
% mip_data.plot(fig_mip, 'InVivo Rat Brain SVD-Filtered CF (MIP over 40 frames)', 60);
% exportgraphics(gca, fullfile(figure_path, 'Beamformed_CF_SVD_filtered_MIP.png'));

%% ULM Process
u = ulm.ULM();
u.framerate = 500;
u.fwhm = [5 5];
u.numberOfParticles = 40;
u.NLocalMax = 2;
u.max_linking_distance = 3; 
u.min_length = 15;
u.algorithm = ulm.algorithm.radial;
u.tracking = ulm.tracking.tracks;
u.lambda = lambda;
u.input = b_data_all_CF;
u.scan = invivo_scan;

% Execute ULM process
tracks = u.go();


%% Creating ULM Images

ulm_img = u.create_image(tracks, "tracks");
fig1 = figure('Visible', 'off');
imagesc(invivo_scan.x_axis * 1e3, invivo_scan.z_axis * 1e3, ulm_img);
xlabel("X [mm]");ylabel("z [mm]");
title('ULM InVivo Rat Brain - Sliding SVD (5 files)')
colormap turbo;
exportgraphics(gca, fullfile(figure_path, 'ULM_InVivo_Rat_Brain_sliding_SVD_5.png'));

%% With interpolation
u.tracking = ulm.tracking.velocity_interpolation;
ulm_img = u.create_image(u.go(), "tracks");
fig2 = figure('Visible', 'off');
imagesc(invivo_scan.x_axis * 1e3, invivo_scan.z_axis * 1e3, ulm_img);
xlabel("X [mm]");ylabel("z [mm]");
title('ULM InVivo Rat Brain with interpolation - Sliding SVD (5 files)')
colormap turbo;
exportgraphics(gca, fullfile(figure_path, 'ULM_InVivo_Rat_Brain_interpolated_sliding_SVD_5.png'));


