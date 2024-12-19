%% Getting channel data
%
clear all;
close all;


%% Choose adaptive algoritm
MINIMUM_VARIANCE = true;
PHASE_COHERENCE_FACTOR = false;
SHORT_LAG_SPATIAL_COHERENCE = false;

% data location
url='http://ustb.no/datasets/';      % if not found data will be downloaded from here

filename='PICMUS_numerical_calib_v2.uff'; data_label = 'CPWC';
%filename='FieldII_STAI_dynamic_range.uff';
%filename='L7_FI_IUS2018.uff'; data_label = 'RTB';
%filename='P4_2v_DW134847.uff'; data_label = 'DW';

% checks if the data is in your data path, and downloads it otherwise.
% The defaults data path is under USTB's folder, but you can change this
% by setting an environment variable with setenv(DATA_PATH,'the_path_you_want_to_use');
tools.download(filename, url, data_path);   

%filename = [data_path() filesep '']; 
channel_data = uff.channel_data();
channel_data.read([data_path,filesep,filename],'/channel_data');


%%
if strcmp(data_label,'RTB')
    MLA = 4;
    scan=uff.linear_scan('x_axis',linspace(channel_data.sequence(1).source.x,channel_data.sequence(end).source.x,channel_data.N_waves.*MLA).', 'z_axis', linspace(10e-3,50e-3,256).');
elseif strcmp(data_label,'DW')
    depth_axis=linspace(20e-3,90e-3,512).';
    angles_axis = linspace(deg2rad(-40),deg2rad(40),128);
    scan=uff.sector_scan('azimuth_axis',angles_axis.','depth_axis',depth_axis);
else
    scan=uff.linear_scan('x_axis',linspace(-20e-3,20e-3,256).', 'z_axis', linspace(10e-3,50e-3,512).');
end

% pipeline
pipe = pipeline();

pipe.channel_data=channel_data;
pipe.scan=scan;

pipe.transmit_apodization.probe = channel_data.probe;
pipe.transmit_apodization.sequence = channel_data.sequence;
pipe.transmit_apodization.focus= scan;
if strcmp(data_label,'RTB') || strcmp(data_label,'DW')
    pipe.transmit_apodization.window = uff.window.hamming;
else
    pipe.transmit_apodization.window=uff.window.boxcar;
end
pipe.transmit_apodization.minimum_aperture = [3.07000e-03 3.07000e-03];
pipe.transmit_apodization.f_number=1.75;

pipe.receive_apodization.probe = channel_data.probe;
pipe.receive_apodization.focus= scan;
pipe.receive_apodization.window=uff.window.hamming;
pipe.receive_apodization.f_number=1.75;


% By calling these we do the calculation of the apodization here
% and thus not affecting the timing later on.
pipe.transmit_apodization.data;
pipe.receive_apodization.data;


das=midprocess.das();
%das.pw_margin = 1e-3;
%%
if MINIMUM_VARIANCE
    adapt_rx = postprocess.capon_minimum_variance();
    adapt_rx.dimension = dimension.receive;
    adapt_rx.scan = scan;
    adapt_rx.channel_data = channel_data;
    adapt_rx.K_in_lambda = 1.5;
    adapt_rx.L_elements = round(3*channel_data.probe.N/4);
    adapt_rx.regCoef = 1/100;
    
    adapt_tx = postprocess.capon_minimum_variance();
    adapt_tx.dimension = dimension.transmit;
    adapt_tx.scan = scan;
    adapt_tx.channel_data = channel_data;
    adapt_tx.K_in_lambda = 1.5;
    adapt_tx.L_elements = channel_data.probe.N/2;
    adapt_tx.regCoef = 1/100;
    
    compression = 'log';
    adapt_label = 'MV'
    dynamic_range = 60;
end

if PHASE_COHERENCE_FACTOR
    adapt_rx = postprocess.phase_coherence_factor();
    adapt_rx.dimension = dimension.receive;
    
    adapt_tx = postprocess.phase_coherence_factor();
    adapt_tx.dimension = dimension.transmit;
    
    compression = 'log';
    adapt_label = 'PCF'
    dynamic_range = 80;
end

if SHORT_LAG_SPATIAL_COHERENCE
    adapt_rx = postprocess.short_lag_spatial_coherence();
    adapt_rx.dimension = dimension.receive;
    adapt_rx.channel_data = channel_data;
    adapt_rx.maxM = 10;
    adapt_rx.scan = scan;
    adapt_rx.K_in_lambda = 1;
    
    adapt_tx = postprocess.short_lag_spatial_coherence();
    adapt_tx.dimension = dimension.transmit;
    adapt_tx.channel_data = channel_data;
    adapt_tx.maxM = 10;
    adapt_tx.scan = scan;
    adapt_tx.K_in_lambda = 1;
    
    %pipe.transmit_apodization.window=uff.window.none;
    pipe.receive_apodization.window=uff.window.none;
    % By calling these we do the calculation of the apodization here
    % and thus not affecting the timing later on.
    pipe.transmit_apodization.data;
    pipe.receive_apodization.data;
    
    compression = 'none';
    adapt_label = 'SLSC'
    dynamic_range = 60;
end

%% DAS on RX -> DAS on Tx : Conventional
das.dimension = dimension.both;

% Reset the hash so that everything is recalculated
das.reset_hash();

tic()
img{1} = pipe.go({das});
if strcmp(data_label,'RTB')
     img{1}.data = img{1}.data.*1./sum(pipe.transmit_apodization.data,2);
end
time{1}  = toc()
label{1} = ['DASonRX->DASonTX'];


%% Adapt on Rx -> DAS on Tx
das.dimension = dimension.none;

if MINIMUM_VARIANCE
    pipe.transmit_apodization.window=uff.window.none;
end

das_tx=postprocess.coherent_compounding();
das_tx.dimension = dimension.transmit;

tic()
img{2} = pipe.go({das adapt_rx das_tx});
if strcmp(data_label,'RTB')
     img{2}.data = img{2}.data.*1./sum(pipe.transmit_apodization.data,2);
end
time{2}  = toc();
label{2} = [adapt_label,'onRX->DASonTX'];

%% DAS on Rx -> Adapt on Tx
das.dimension = dimension.receive;

if MINIMUM_VARIANCE
    adapt_tx.L_elements = floor(channel_data.N_waves/3);
    adapt_tx.regCoef = 1/100;
    pipe.transmit_apodization.window=uff.window.none;
end

% Reset the hash so that everything is recalculated
das.reset_hash();

tic()
img{3} = pipe.go({das adapt_tx});
if strcmp(data_label,'RTB')
     %img{3}.data = img{3}.data.*1./sum(pipe.transmit_apodization.data,2);
end
time{3}  = toc()
label{3} = ['DASonRX->',adapt_label,'onTX'];

%% Adapt on Rx -> Adapt on Tx : Double adaptive
das.dimension = dimension.none;

if MINIMUM_VARIANCE
    adapt_tx.L_elements = floor(channel_data.N_waves/3);
    adapt_tx.regCoef = 1/100;
    pipe.transmit_apodization.window=uff.window.none;
end

% Reset the hash so that everything is recalculated
das.reset_hash();
adapt_rx.reset_hash();
adapt_tx.reset_hash();

tic();
img{4} = pipe.go({das adapt_rx adapt_tx});
time{4}  = toc()

label{4} = [adapt_label,'onRX->',adapt_label,'onTX'];

%% DAS on TX -> Adapt on RX
das.dimension = dimension.transmit;

if MINIMUM_VARIANCE
    adapt_rx.L_elements = floor(channel_data.probe.N_elements/2);
    adapt_rx.regCoef = 1/100;
end
    

tic()
img{5} = pipe.go({das adapt_rx});
if strcmp(data_label,'RTB')
     img{5}.data = img{5}.data.*1./sum(pipe.transmit_apodization.data,2);
end
time{5} = toc()
label{5} = ['DASonTX->',adapt_label, 'onRX'];

%%



img{5}.plot(figure(123),label{5},dynamic_range,compression);



%% Plot Figure
figure(3);
img{1}.plot(subplot(321),label{1},dynamic_range,compression)
ax(1) = gca;
img{2}.plot(subplot(322),label{2},dynamic_range,compression)
ax(2) = gca;
img{3}.plot(subplot(323),label{3},dynamic_range,compression);
ax(3) = gca;
img{4}.plot(subplot(324),label{4},dynamic_range,compression);
ax(4) = gca;
img{5}.plot(subplot(325),label{5},dynamic_range,compression);
ax(5) = gca;
linkaxes(ax);
subplot(3,2,[6]);
bar([time{:}]/60)
xticklabels(label)
ylabel('Computation time (m)');
set(gca,'FontSize',14);
set(gcf,'Position',[29 33 1109 715]);

%%
folder_path = ['figures'];
mkdir(folder_path)

b_data_compare = uff.beamformed_data(img{1})
b_data_compare.data(:,1,1,1) = img{1}.data;
close all;
for i = 1:length(img)
    b_data_compare.data(:,1,1,i) = img{i}.data./max(img{i}.data);
    f = figure(i);clf;
    imagesc(scan.x_axis*1000,scan.z_axis*1000,img{i}.get_image(compression));
    colormap gray; caxis([-dynamic_range 0]);axis image;
    xlabel('x [mm]');ylabel('y [mm]');
    set(gca,'FontSize',14)
    img{i}.plot(f,[],dynamic_range,compression)
    saveas(f,[folder_path,filesep,strrep(label{i}, '->', '_')],'eps2c')
    axis([0 10 20 28]);
    saveas(f,[folder_path,filesep,strrep(label{i}, '->', '_'),'_zoomed'],'eps2c')
end

%%
b_data_compare.plot()
%%
running_time_in_m = [time{:}]/60;
f99 = figure(100);clf;
%subplot(1,2,1);
bar([time{:}]/60)

labels_latex{1} = ['$b^{\overline{T_{x DAS}}~\overline{R_{x DAS}}}$'];
labels_latex{2} = ['$b^{\overline{T_{x DAS}}~\overline{R_{x MV}}}$'];
labels_latex{3} = ['$b^{\overline{T_{x MV}}~\overline{R_{x DAS}}}$'];
labels_latex{4} = ['$b^{\overline{T_{x MV}}~\overline{R_{x MV}}}$'];
labels_latex{5} = ['$b^{\overline{R_{x MV}}~\overline{T_{x DAS}}}$'];

ylim([0 max(running_time_in_m)+1]);
ylabel('Computation time (m)');
set(gca,'FontSize',18);
x_pos = [1 2 3 4 5];
text(x_pos,running_time_in_m,num2str(running_time_in_m(:),'%.2f'),'HorizontalAlignment','center','VerticalAlignment','bottom','FontSize',18)
set(gca, 'XTickLabel', labels_latex, 'TickLabelInterpreter', 'latex','FontSize',24);
set(gcf,'Position',[250 344 781 447]);grid on
saveas(f99,[folder_path,filesep,'adaptive_timing'],'eps2c')
