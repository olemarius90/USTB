%% In silico correction of blocked array using RTB and REFoCUS
%
% Creating the simulated results in the paper:
%
% A. E. Vrålstad, S.-E. Måsøy, T. G. Bjåstad, A. R. Sørnes and O. M. H. Rindal,
% "Retrospective transmit correction of blocked arrays applied to cardiac
% ultrasound imaging," IEEE Trans. Ultrason. (T-USON) (to appear).
%
% This script creates Figs. 3-4 in the paper.
%
% This code uses the UltraSound ToolBox (USTB). Clone or add USTB to the
% MATLAB path before running. See <https://github.com/unioslo/USTB> and
% the project website.
%
% Authors: Anders E. Vrålstad, Ole Marius Høel Rindal
%% Clear environment
clear all; close all;
% Headless / publish / CI: no interactive ROI drawing; no demod figure
headless = ~usejava('desktop');
if exist('batchStartupOptionUsed', 'file') == 2 %#ok<*EXIST>
    try
        headless = headless | batchStartupOptionUsed;
    catch
    end
end
if headless
    set(groot, 'DefaultFigureVisible', 'off');
end

% Figures that must appear in MATLAB publish()+snapnow (-nodisplay) need Visible 'on' even in
% headless mode; invisible + colorbar reliably drops snapshots for publish (triptych, diff maps).
fv_publish = tools.headless_publish_figure_visible();
%% Load data
% Prefer Zenodo (CI); fall back to ustb.no mirrors — speckle_sim_FI_P4_* are not on Zenodo 19550715.
local_path = [ustb_path(), '/data/'];
addpath(local_path);

%% Speckle chamber ROI — three simulated blockage levels
% Red rectangle matches gCNR `center_rectangle` = [-0.006, 0.07, 0.007, 0.04] m on the
% Cartesian scan-converted grid (shown here in mm): [-6, 70, 7, 40].
% Three panel images in **one figure** (`subplot(1,3,:)`): one embedded snapshot for this cell.
% MATLAB publish() does not reliably embed snapshots for `beamformed_data.plot(ax, …)`
% when `ax` belongs to an ad-hoc invisible figure axis handle; passing a explicit figure
% (as in Save PNGs) works — see +tools/publish_beamformed_snap.
speckle_chamber_mm = [-6, 70, 7, 40];
dataset_roi = {
    'speckle_sim_FI_P4_probe_apod_1_speckle_long_many_angles.uff', 'Full aperture (no element blockage)'
    'speckle_sim_FI_P4_probe_apod_2_speckle_long_many_angles.uff', 'One-third of aperture blocked'
    'speckle_sim_FI_P4_probe_apod_3_speckle_long_many_angles.uff', 'Half of aperture blocked'
    };

figure('Visible', fv_publish);
set(gcf, 'Position', [100, 100, 1200, 400]);
for r = 1:size(dataset_roi, 1)
    fn_roi = dataset_roi{r, 1};
    roi_title = dataset_roi{r, 2};
    correction_download_uff_with_fallbacks(fn_roi, data_path());
    cd_roi = uff.read_object([data_path(), filesep, fn_roi], '/channel_data');
    cd_roi.data = cd_roi.data ./ max(cd_roi.data(:));

    depth_roi = linspace(0e-3, 110e-3, 512).';
    az_roi = zeros(cd_roi.N_waves, 1);
    for n = 1:cd_roi.N_waves
        az_roi(n) = cd_roi.sequence(n).source.azimuth;
    end
    scan_roi = uff.sector_scan('azimuth_axis', az_roi, 'depth_axis', depth_roi);
    Fnr = cd_roi.sequence(1).source.distance / (max(cd_roi.probe.x) * 2);

    mid_roi = midprocess.das();
    mid_roi.channel_data = cd_roi;
    mid_roi.dimension = dimension.both();
    mid_roi.scan = scan_roi;
    if isempty(which('das_c'))
        mid_roi.code = code.matlab;
    else
        mid_roi.code = code.mex;
    end
    mid_roi.receive_apodization.window = uff.window.boxcar;
    mid_roi.receive_apodization.f_number = 1.7;
    mid_roi.transmit_apodization.window = uff.window.hamming;
    mid_roi.transmit_apodization.f_number = Fnr;
    mid_roi.transmit_apodization.minimum_aperture = 3e-3;
    b_roi = mid_roi.go();

    [RTB_sc_roi, Xr, Zr] = tools.scan_convert(b_roi.get_image(), b_roi.scan.azimuth_axis, b_roi.scan.depth_axis, 1024, 1024);

    subplot(1, 3, r);
    imagesc(Xr * 1e3, Zr * 1e3, RTB_sc_roi, [-60 0]);
    colormap(gca, gray);
    colorbar;
    axis image;
    xlabel('x [mm]');
    ylabel('z [mm]');
    title(roi_title, 'Interpreter', 'none');
    xlim([-20 20]);
    ylim([60 110]);
    hold on;
    rectangle(gca, 'Position', speckle_chamber_mm, 'EdgeColor', 'r', 'LineWidth', 2);
    hold off;
end
set(findall(gcf, '-property', 'FontSize'), 'FontSize', 12);
drawnow;
snapnow;
% Do not close figures here — `publish` snaps at cell boundaries and the next `close all`
% would discard the ROI triptych before it is flushed to HTML.

%% Primary case — full pipeline (paper figures, gCNR, differences): half aperture blocked
filename = 'speckle_sim_FI_P4_probe_apod_3_speckle_long_many_angles.uff';
tag = 'half';
correction_download_uff_with_fallbacks(filename, data_path());

channel_data = uff.read_object([data_path(), filesep, filename], '/channel_data');
channel_data.data = channel_data.data ./ max(channel_data.data(:));

storefolder = ['./Figures/simulated_gCNR_',tag, '/'];
mkdir(storefolder);
%% Run REFoCUS preprocess
tic
REFoCUS = preprocess.refocus();
REFoCUS.input = channel_data;
REFoCUS.use_filter = 0;
REFoCUS.filter_N = 10;
REFoCUS.filter_Wn = [0.05,0.4];
REFoCUS.regularization = @Hinv_tikhonov;
REFoCUS.decode_parameter = 0.01;
REFoCUS.post_pad_samples = 0;
channel_data_REFoCUS = REFoCUS.go();
toc()

%% Create Sector scan
depth_axis=linspace(0e-3,110e-3,512).';
azimuth_axis=zeros(channel_data.N_waves,1);
for n=1:channel_data.N_waves
    azimuth_axis(n) = channel_data.sequence(n).source.azimuth;
end
scan=uff.sector_scan('azimuth_axis',azimuth_axis,'depth_axis',depth_axis);
%% RTB Beamforming
% Calculate F-number
Fnumber = channel_data.sequence(1).source.distance/(max(channel_data.probe.x)*2);

mid=midprocess.das();
mid.channel_data=channel_data;
mid.dimension = dimension.both();
mid.scan=scan;
if isempty(which('das_c'))
    mid.code = code.matlab;
else
    mid.code = code.mex;
end
mid.receive_apodization.window=uff.window.boxcar;
mid.receive_apodization.f_number=1.7;
mid.transmit_apodization.window=uff.window.hamming;
mid.transmit_apodization.f_number=Fnumber;
mid.transmit_apodization.minimum_aperture = 3e-3;
b_data_RTB = mid.go();
b_data_RTB.frame_rate = 20;
tools.publish_beamformed_snap(b_data_RTB, 'RTB');
%% Store the Original RTB weights for plotting later
b_data_tx_apod = uff.beamformed_data(b_data_RTB);
b_data_tx_apod.data = mid.transmit_apodization.data;
tools.publish_beamformed_snap(b_data_tx_apod, ['Tx Weights no shift'], [], 'none');
colormap default;
%% Change beam geometry for RTB processing
switch tag
    case 'full'
        origin_x = 0;
    case 'third'
        origin_x = max(channel_data.probe.x)*(2*2/3-1);
        Fnumber = Fnumber*3/2;
    case 'half'
        origin_x = max(channel_data.probe.x)*0.5;
        Fnumber = Fnumber*2;
end
channel_data_shifted = uff.channel_data(channel_data);
for seq = 1:channel_data_shifted.N_waves
    channel_data_shifted.sequence(seq).origin.x = origin_x;
end
mid.channel_data = channel_data_shifted;
mid.transmit_apodization.f_number=Fnumber;
b_data_RTB_comp = mid.go();
b_data_RTB_comp.frame_rate = 20;
tools.publish_beamformed_snap(b_data_RTB_comp, 'Proposed RTB');
%% Store the Compensated RTB weights for plotting later
b_data_tx_apod = uff.beamformed_data(b_data_RTB_comp);
b_data_tx_apod.data = mid.transmit_apodization.data;
tools.publish_beamformed_snap(b_data_tx_apod, ['Tx Weights with shift'], [], 'none');
colormap default;
%% Demodulate REFoCUS RF-data before DAS
demod = preprocess.fast_demodulation();
demod.modulation_frequency = 2.5*10^6;
demod.input = channel_data_REFoCUS;
demod.plot_on = ~headless;
channel_data_STAI_demod = demod.go();
%% REFoCUS Beamforming
mid_REFoCUS = midprocess.das();
mid_REFoCUS.channel_data=channel_data_STAI_demod;
mid_REFoCUS.dimension = dimension.receive();
mid_REFoCUS.scan=scan;
if isempty(which('das_c'))
    mid_REFoCUS.code = code.matlab;
else
    mid_REFoCUS.code = code.mex;
end
mid_REFoCUS.transmit_apodization.window=uff.window.boxcar;
mid_REFoCUS.receive_apodization.f_number=1.7;
mid_REFoCUS.receive_apodization.window=uff.window.boxcar;
mid_REFoCUS.transmit_apodization.f_number=1;
b_data_delayed_REFoCUS = mid_REFoCUS.go();
tools.publish_beamformed_snap(b_data_delayed_REFoCUS, 'REFoCUS');
cc = postprocess.coherent_compounding;
cc.input = b_data_delayed_REFoCUS;
b_data_REFoCUS = cc.go();
tools.publish_beamformed_snap(b_data_REFoCUS);

%% Save PNGs
f = figure('Visible', 'off');
b_data_RTB.plot(f,'RTB');
rectangle(gca,'Position',[-6 5 7 105],'EdgeColor','r','LineWidth',2)
clim([-60 0]);xlim([-20 20]);
savefig(f,[storefolder,'RTB_', tag,'.fig']);
saveas(f,[storefolder,'RTB_', tag,'.png']);
drawnow;
snapnow;

b_data_RTB_comp.plot(f,'RTB Compensated');
rectangle(gca,'Position',[-6 5 7 105],'EdgeColor','r','LineWidth',2)
clim([-60 0]);xlim([-20 20]);
savefig(f,[storefolder,'RTB_compensated_', tag,'.fig']);
saveas(f,[storefolder,'RTB_compensated_', tag,'.png']);
drawnow;
snapnow;

b_data_REFoCUS.plot(f,'REFoCUS');
rectangle(gca,'Position',[-6 5 7 105],'EdgeColor','r','LineWidth',2)
clim([-60 0]);xlim([-20 20])
savefig(f,[storefolder,'REFoCUS_', tag,'.fig']);
saveas(f,[storefolder,'REFoCUS_', tag,'.png']);
drawnow;
snapnow;

%% Comparison cube (difference images; optional GIF when not headless)
b_data_compare = uff.beamformed_data(b_data_RTB);
b_data_compare.data(:, 1, 1, 1) = b_data_RTB.data ./ max(b_data_RTB.data(:));
b_data_compare.data(:, 1, 1, 2) = b_data_RTB_comp.data ./ max(b_data_RTB_comp.data(:));
b_data_compare.data(:, 1, 1, 3) = b_data_REFoCUS.data ./ max(b_data_REFoCUS.data(:)) / 3;
all_images_cmp = squeeze(b_data_compare.get_image());
b_data_compare.data(:, 1, 1, 2) = b_data_compare.data(:, 1, 1, 2) .* median(all_images_cmp(:, :, 1) ./ all_images_cmp(:, :, 2), 'all', 'omitnan');
b_data_compare.data(:, 1, 1, 3) = b_data_compare.data(:, 1, 1, 3) .* median(all_images_cmp(:, :, 1) ./ all_images_cmp(:, :, 3), 'all', 'omitnan');

if ~headless
    b_data_compare.frame_rate = 1;
    b_data_compare.plot([]);
    rectangle(gca, 'Position', [-6 5 7 105], 'EdgeColor', 'r', 'LineWidth', 2)
    clim([-60 0]); xlim([-20 20]);
    b_data_compare.frame_rate = 1;
    b_data_compare.save_as_gif(['Figures/Comparison_', tag, '.gif']);
else
    fprintf('[Correction_of_simulated_blockage] Skipping animated GIF (headless publish).\n');
    % Still embed the three matched images in published HTML (publish ignores invisible figures without snapnow).
    cmp_labs = {'RTB (amp. matched)', 'Proposed RTB', 'REFoCUS'};
    for k = 1:3
        bk = uff.beamformed_data(b_data_RTB);
        bk.data = b_data_compare.data(:, 1, 1, k);
        fh_cmp = figure('Visible', fv_publish);
        bk.plot(fh_cmp, cmp_labs{k});
        rectangle(gca, 'Position', [-6 5 7 105], 'EdgeColor', 'r', 'LineWidth', 2);
        clim([-60 0]); xlim([-20 20]);
        drawnow;
        delete(findall(fh_cmp, 'Type', 'uicontrol'));
        snapnow;
        close(fh_cmp);
    end
end


%% Measure contrast (interactive drawrectangle; skipped in headless publish/CI)
if ~headless
    [RTB_sc, Xs,Zs] = tools.scan_convert(b_data_RTB.get_image(),b_data_compare.scan.azimuth_axis,b_data_compare.scan.depth_axis, 1024,1024);
    [RTB_comp_sc, Xs,Zs] = tools.scan_convert(b_data_RTB_comp.get_image(),b_data_compare.scan.azimuth_axis,b_data_compare.scan.depth_axis, 1024,1024);
    [REFoCUS_sc, Xs,Zs] = tools.scan_convert(b_data_REFoCUS.get_image()-12,b_data_compare.scan.azimuth_axis,b_data_compare.scan.depth_axis, 1024,1024);
    img_cell = {RTB_sc,RTB_comp_sc,REFoCUS_sc};
    name_cell = {'RTB','RTB Compensated', 'REFoCUS'};
    das_handle = figure(); imagesc(Xs,Zs,img_cell{3});


    center_rectangle = [-0.006,0.07,0.007,0.04];
    v1_area =drawrectangle('Position',center_rectangle-[0.008,0,0,0]);
    v2_area =drawrectangle('Position',center_rectangle+[0.008,0,0,0]);
    c_area =drawrectangle('Position',center_rectangle);

    [GCNR, v1_binary, v2_binary, c_binary] = contrast_calc_insilico(img_cell,name_cell, Xs, Zs, das_handle, 60,storefolder,c_area,v1_area,v2_area);
else
    fprintf('[Correction_of_simulated_blockage] Skipping interactive gCNR (headless run).\n');
end

%% Make Difference Images: of RTBs
all_images = b_data_compare.get_image();
diff = all_images(:,:,2)-all_images(:,:,1)-1.6;
[diff_sc, Xs,Zs] = tools.scan_convert(diff, scan.azimuth_axis, scan.depth_axis,1024,1024);
f = figure('Visible', fv_publish);
imagesc(Xs*1e3,Zs*1e3,diff_sc,[-30,30]);colormap(bluewhitered);colorbar
xlabel('x[mm]')
ylabel('z[mm]')
axis image
xlim([-20,20]);
ylim([60,110]);
set(findall(gcf,'-property','FontSize'),'FontSize',15)
savefig(f,[storefolder,'RTBminusRTBcomp_',tag,'.fig']);
saveas(f,[storefolder,'RTBminusRTBcomp_',tag,'.png']);
drawnow;
snapnow;

%% Make Difference Images: RTB comp and REFoCUS 
all_images = b_data_compare.get_image();
diff = all_images(:,:,3)-all_images(:,:,2)-7;
[diff_sc, Xs,Zs] = tools.scan_convert(diff, scan.azimuth_axis, scan.depth_axis,1024,1024);
f = figure('Visible', fv_publish);
imagesc(Xs*1e3,Zs*1e3,diff_sc,[-30,30]);colormap(bluewhitered);colorbar
xlabel('x[mm]')
ylabel('z[mm]')
axis image
xlim([-20,20]);
ylim([60,110]);

set(findall(gcf,'-property','FontSize'),'FontSize',15)
savefig(f,[storefolder,'REFoCUSminusRTBcomp_',tag,'.fig']);
saveas(f,[storefolder,'REFoCUSminusRTBcomp_',tag,'.png']);
drawnow;
snapnow;

%% Make Difference Images: RTB and REFoCUS 
all_images = b_data_compare.get_image();
diff = all_images(:,:,3)-all_images(:,:,1)-7;
[diff_sc, Xs,Zs] = tools.scan_convert(diff, scan.azimuth_axis, scan.depth_axis,1024,1024);
f = figure('Visible', fv_publish);
imagesc(Xs*1e3,Zs*1e3,diff_sc,[-30,30]);colormap(bluewhitered);colorbar
xlabel('x[mm]')
ylabel('z[mm]')
axis image
xlim([-20,20]);
ylim([60,110]);

set(findall(gcf,'-property','FontSize'),'FontSize',15)
savefig(f,[storefolder,'REFoCUSminusRTB_',tag,'.fig']);
saveas(f,[storefolder,'REFoCUSminusRTB_',tag,'.png']);
drawnow;
snapnow;

