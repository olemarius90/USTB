%% Computation of a FI dataset with Field II using a 2D array and beamforming with USTB
%
% authors:  Ole Marius Hoel Rindal <olemarius@olemarius.net>
%
% Last updated: 12.11.2017
clear all;
close all;

%sequence_type = 'PW';
sequence_type = 'STA';
save_folder = '2D_array_STAI_simulation'
mkdir(save_folder)
N_transmit = 31*31;
%% basic constants
c0=1540;     % Speed of sound [m/s]
fs=100e6;    % Sampling frequency [Hz]
dt=1/fs;     % Sampling step [s]

%% field II initialisation
field_init(0);
set_field('c',c0);              % Speed of sound [m/s]
set_field('fs',fs);             % Sampling frequency [Hz]
set_field('use_rectangles',1);  % use rectangular elements

f0=2.56e6;          % Transducer center frequency [Hz]
bw=0.67;            % probe bandwidth [1]
lambda=c0/f0;       % Wavelength [m]
pulse_duration=2.5; % pulse duration [cycles]
%% Create idealized Matrix Array
x_size = 30*lambda/2;
y_size = 30*lambda/2;

element_center_x = -x_size/2 : lambda/2 : x_size/2
element_center_y = -y_size/2 : lambda/2 : y_size/2

N_x = length(element_center_x);
N_y = length(element_center_y);
full_geometry(:,1) = repmat(element_center_x',N_y,1);
for y = 1:N_y
    idx = (y-1)*N_x+1:(y)*N_x;
    full_geometry(idx,2) = ones(N_x,1)*element_center_y(y);
end
full_geometry(:,3) = zeros(N_x*N_y,1);

kerf=lambda/5;
probe = uff.matrix_array();
probe.pitch_x = lambda/2;
probe.pitch_y = lambda/2;
probe.element_width = probe.pitch_x-kerf;
probe.element_height = probe.pitch_x-kerf;
probe.N_x = N_x;
probe.N_y = N_y;
probe.geometry = full_geometry;
center_point = uff.point();
center_point.xyz = [mean(probe.x) mean(probe.y) mean(probe.z)];
probe.origin = center_point;
probe.plot()

%% pulse definition
pulse = uff.pulse();
pulse.center_frequency = f0;
pulse.fractional_bandwidth = 0.65;        % probe bandwidth [1]
t0 = (-1/pulse.fractional_bandwidth/f0): dt : (1/pulse.fractional_bandwidth/f0);
impulse_response = gauspuls(t0, f0, pulse.fractional_bandwidth);
impulse_response = impulse_response-mean(impulse_response); % To get rid of DC

te = (-pulse_duration/2/f0): dt : (pulse_duration/2/f0);
excitation = square(2*pi*f0*te+pi/2);
one_way_ir = conv(impulse_response,excitation);
two_way_ir = conv(one_way_ir,impulse_response);
lag = length(two_way_ir)/2+1;

% show the pulse to check that the lag estimation is on place (and that the pulse is symmetric)
figure;
plot((0:(length(two_way_ir)-1))*dt -lag*dt,two_way_ir); hold on; grid on; axis tight
plot((0:(length(two_way_ir)-1))*dt -lag*dt,abs(hilbert(two_way_ir)),'r')
plot([0 0],[min(two_way_ir) max(two_way_ir)],'g');
legend('2-ways pulse','Envelope','Estimated lag');
title('2-ways impulse response Field II');

%% Phantom
if 0
    %% Create a phantom to image
    phantom_positions(1,:) = [0 0 60/1000];
    phantom_positions(2,:) = [0 0/1000 20/1000];
    phantom_positions(3,:) = [0 0/1000 40/1000];
    phantom_positions(4,:) = [0 0 80/1000];
    phantom_positions(5,:) = [0 0/1000 100/1000];
    phantom_amplitudes(1:size(phantom_positions,1)) = 1;
else
  %  addpath = [ustb_path,'\publications\DynamicRange\make_simulation\'];
    %% speckle
    [sca,amp] = simulatedPhantomGradientBlock3D(5,-20/1000,20/1000,45/1000,75/1000,-20/1000,20/1000,0,0)

    % Define geometry of nonechoic cyst
    r=10e-3;           % Radius of cyst [m]
    xc=0e-3;           % Position of cyst in x [m]
    zc=60e-3;           % Position of cyst in z [m]
    yc=0;
    %Find the indexes inside cyst
    inside_nonechoic = (((sca(:,1)-xc).^2 + (sca(:,3)-zc).^2) < r^2);
    %clear sca_hc;
    sca_temp = sca;
    clear sca;
    
     sca(:,1) = sca_temp(inside_nonechoic==0,1);
     sca(:,2) = sca_temp(inside_nonechoic==0,2);
     sca(:,3) = sca_temp(inside_nonechoic==0,3);
     amp = amp(inside_nonechoic==0);

    

    phantom_positions = sca;
    phantom_amplitudes = amp';
    figure;
    plot3(sca(:,1)*1e3,sca(:,2)*1e3,sca(:,3)*1e3,'b.'); hold on; axis equal; grid on;
    %plot3(sca(~mask,1)*1e3,sca(~mask,3)*1e3,20*log10(amp(~mask)),'r.');
    xlabel('x[mm]');
    ylabel('y[mm]');
    zlabel('z[mm]');
    
end

%% N_transmit hack

F_number = 1.7;
alpha_max = deg2rad(-3);%atan(1/2/F_number);
Na=N_transmit;                                      % number of plane waves 
F=1;                                        % number of frames
alpha=linspace(-alpha_max,alpha_max,Na);    % vector of angles [rad]
 
% %% output data
% cropat=round(1.1*2*sqrt((max(phantom_positions(:,1))-min(probe.x))^2+max(phantom_positions(:,3))^2)/c0/dt);   % maximum time sample, samples after this will be dumped
% t_out=0:dt:((cropat-1)*dt);                 % output time vector
% %STA=zeros(cropat,probe.N_x*probe.N_y,probe.N_x*probe.N_y);    % impulse response channel data
% STA=zeros(cropat,probe.N_x*probe.N_y,N_transmit);
%% Compute STA signals
disp('Field II: Computing FI dataset');
cropat=round(1.1*2*sqrt((max(phantom_positions(:,1))-min(probe.x))^2+max(phantom_positions(:,3))^2)/c0/dt);   % maximum time sample, samples after this will be dumped
%STA=zeros(cropat,probe.N_x*probe.N_y,N_transmit);
parfor n=1:N_transmit%probe.N_elements
    %% output data
       t_out=0:dt:((cropat-1)*dt);                 % output time vector
    STA=zeros(cropat,probe.N_x*probe.N_y);    % impulse response channel data
    

    fprintf('Simulating transmit %d / %d\n',n,probe.N_elements);
    % Define Th and Rh in loop to be able to do parfor
    % definition of the mesh geometry
    noSubAz=round(probe.element_width/(lambda/8));        % number of subelements in the azimuth direction
    noSubEl=round(probe.element_height/(lambda/8));       % number of subelements in the elevation directio
    enabled = ones(probe.N_x,probe.N_y);
    field_init(0);
    Th = xdc_2d_array (probe.N_x, probe.N_y, probe.element_width, probe.element_height, kerf, kerf, enabled, noSubAz, noSubEl, [0 0 Inf]);
    Rh = xdc_2d_array (probe.N_x, probe.N_y, probe.element_width, probe.element_height, kerf, kerf, enabled, noSubAz, noSubEl, [0 0 Inf]);

    % setting excitation, impulse response and baffle
    xdc_excitation (Th, excitation);
    xdc_impulse (Th, impulse_response);
    xdc_baffle(Th, 0);
    xdc_impulse (Rh, impulse_response);
    xdc_baffle(Rh, 0);
    xdc_center_focus(Rh,[0 0 0]);

    %% Set transmit element
    %xdc_apodization(Th, 0, [zeros(1,n-1) 1 zeros(1,probe.N_x*probe.N_y-n)]);
    
    xdc_apodization(Th, 0, ones(1,probe.N_x*probe.N_y));
    xdc_times_focus(Th, 0, probe.geometry(:,1)'.*sin(alpha(n))/c0);
    %xdc_focus_times(Th, 0, zeros(1,probe.N_x*probe.N_y));
    
    % receive aperture
    xdc_apodization(Rh, 0, ones(1,probe.N_x*probe.N_y));
    xdc_focus_times(Rh, 0, zeros(1,probe.N_x*probe.N_y));

    % do calculation
    [v,t]=calc_scat_multi(Th, Rh, phantom_positions, phantom_amplitudes');

    % save data -> with parloop we need to pad the data
%     if size(v,1)<cropat
%         STA(:,:,n)=padarray(v,[cropat-size(v,1) 0],0,'post');
%     else
%         STA(:,:,n)=v(1:cropat,:);
%     end
    STA(1:size(v,1),:)=v;

    % Sequence generation
    seq(n)=uff.wave();
    seq(n).probe=probe;
    if contains(sequence_type,'STA')
        seq(n).source.xyz=[probe.x(n) probe.y(n) probe.z(n)];
        seq(n).origin.xyz=seq(n).source.xyz;
        seq(n).delay = probe.r(n)/c0-lag*dt+t; % t0 and center of pulse compensation
    elseif contains(sequence_type,'PW')
        seq(n).source.azimuth=alpha(n);
        seq(n).source.distance=Inf;
        seq(n).delay = -lag*dt+t;
    end
    seq(n).sound_speed=c0;

    
    %% CHANNEL DATA
    channel_data_temp = uff.channel_data();
    channel_data_temp.sampling_frequency = fs;
    channel_data_temp.sound_speed = c0;
    channel_data_temp.initial_time = 0;
    channel_data_temp.pulse = pulse;
    channel_data_temp.probe = probe;
    channel_data_temp.sequence = seq(n);
    channel_data_temp.data = STA;
    channel_data_temp.write([save_folder,'/channel_data_tx_',num2str(n)],'/channel_data')
end

%%
channel_data = uff.channel_data();
channel_data.read([save_folder,'/channel_data_tx_',num2str(1)],'/channel_data')

for n=1:N_transmit
    channel_data_temp = uff.channel_data();
    channel_data_temp.read([save_folder,'/channel_data_tx_',num2str(n)],'/channel_data')
    channel_data.sequence(n) = channel_data_temp.sequence;
    channel_data.data(:,:,n) = channel_data_temp.data;
end

%% Create Sector scan
scan = uff.linear_scan();
scan.x_axis = linspace(-35/1000,35/1000,256)';
scan.z_axis = linspace(40/1000,90/1000,256)';
scan.y =  0.*ones(size(scan.y))*channel_data.probe.y(1);%-5/1000;

%% BEAMFORMER
mid=midprocess.das();
mid.channel_data=channel_data;
mid.dimension = dimension.both()
mid.scan=scan;
mid.spherical_transmit_delay_model = spherical_transmit_delay_model.spherical;
mid.receive_apodization.window=uff.window.hamming;
mid.receive_apodization.f_number=1.5;
mid.transmit_apodization.window=uff.window.none;
mid.transmit_apodization.f_number=3;

% Delay the data
b_data_3D_RTB = mid.go();
b_data_3D_RTB.plot()

%% Visualize apodization
b_data_tx_apod = uff.beamformed_data(b_data_3D_RTB)
b_data_tx_apod.data = mid.transmit_apodization.data;

%% Test with 3D scan
scan_3D_linear = uff.linear_scan_3D;
scan_3D_linear.x_axis = linspace(-10/1000,10/1000,25)';
scan_3D_linear.y_axis = linspace(-10/1000,10/1000,25)';
scan_3D_linear.z_axis = linspace(0/1000,110/1000,256)';
scan_3D_linear.plot()

%%
scan_3D_sector = uff.sector_scan_3D;
scan_3D_sector.azimuth_axis = deg2rad(linspace(-30,30,32))';
scan_3D_sector.elevation_axis = deg2rad(linspace(-30,30,32))';
scan_3D_sector.depth_axis = linspace(0/1000,110/1000,256)';
scan_3D_sector.plot()

%%
mid.scan = scan_3D_linear;
mid.dimension = dimension.both();
b_data_3D_linear = mid.go()
img_3D_linear = reshape(b_data_3D_linear.data,scan_3D_linear.N_z_axis,scan_3D_linear.N_x_axis,scan_3D_linear.N_y_axis);
img_3D_linear = img_3D_linear./max(img_3D_linear(:));

%%
mid.scan = scan_3D_sector;
mid.dimension = dimension.both();
b_data_3D_sector = mid.go()

%%
img_3D_sector= reshape(b_data_3D_sector.data,scan_3D_sector.N_depth_axis,scan_3D_sector.N_azimuth_axis,scan_3D_sector.N_elevation_axis);
img_3D_sector = img_3D_sector./max(img_3D_sector(:));

%%

x_cube_spher = reshape(b_data_3D_sector.scan.x,b_data_3D_sector.scan.N_depth_axis,b_data_3D_sector.scan.N_azimuth_axis,b_data_3D_sector.scan.N_elevation_axis);
y_cube_spher = reshape(b_data_3D_sector.scan.y,b_data_3D_sector.scan.N_depth_axis,b_data_3D_sector.scan.N_azimuth_axis,b_data_3D_sector.scan.N_elevation_axis);
z_cube_spher = reshape(b_data_3D_sector.scan.z,b_data_3D_sector.scan.N_depth_axis,b_data_3D_sector.scan.N_azimuth_axis,b_data_3D_sector.scan.N_elevation_axis);
% 
x_cube_cart = reshape(b_data_3D_linear.scan.x,b_data_3D_linear.scan.N_z_axis,b_data_3D_linear.scan.N_x_axis,b_data_3D_linear.scan.N_y_axis);
y_cube_cart = reshape(b_data_3D_linear.scan.y,b_data_3D_linear.scan.N_z_axis,b_data_3D_linear.scan.N_x_axis,b_data_3D_linear.scan.N_y_axis);
z_cube_cart = reshape(b_data_3D_linear.scan.z,b_data_3D_linear.scan.N_z_axis,b_data_3D_linear.scan.N_x_axis,b_data_3D_linear.scan.N_y_axis);

%% Plot 3D plot
f = figure(90);clf;hold all;
idx_x = 20;
idx_y = 20;
idx_z = 200;
surface(squeeze(x_cube_cart(:,:,idx_x)),squeeze(y_cube_cart(:,:,idx_x)),squeeze(z_cube_cart (:,:,idx_x)),db(abs(squeeze(img_3D_linear(:,:,idx_x)./max(max(img_3D_linear(:)))))));
surface(squeeze(x_cube_cart(:,idx_y,:)),squeeze(y_cube_cart(:,idx_y,:)),squeeze(z_cube_cart (:,idx_y,:)),db(abs(squeeze(img_3D_linear(:,idx_y,:)./max(max(img_3D_linear(:)))))));
%alpha 0.5
%surface(squeeze(z_cube_spher(idx_z,:,:)),squeeze(z_cube_spher(idx_z,:,:)),squeeze(z_cube_spher (idx_z,:,:)),db(abs(squeeze(img_3D_linear(idx_z,:,:)./max(max(img_3D_linear(:)))))));
surface(squeeze(x_cube_spher(:,idx_y,:)),squeeze(y_cube_spher(:,idx_y,:)),squeeze(z_cube_spher (:,idx_y,:)),db(abs(squeeze(img_3D_sector(:,idx_y,:)./max(max(img_3D_sector(:)))))));
axis image; title('3D','Color','white');
xlabel('x [mm]');ylabel('y [mm]');zlabel('z [mm]');
shading('flat');colormap gray; clim([-60 0]);view(40,0)
set(gca,'YColor','white','XColor','white','ZColor','white','Color','k');
set(gcf,'color','black'); set(gca,'Color','k')
f_handle = gcf;set(gca,'ZDir','reverse')
f_handle.InvertHardcopy = 'off';

%% Plot 3D plot
f = figure(91);clf;hold all;
idx_x = size(y_cube_spher,3)/2;
idx_y = size(y_cube_spher,2)/2;
idx_z = 200;
surface(squeeze(x_cube_spher(:,idx_y,:)),squeeze(y_cube_spher(:,idx_y,:)),squeeze(z_cube_spher (:,idx_y,:)),db(abs(squeeze(img_3D_sector(:,idx_y,:)./max(max(img_3D_sector(:)))))));
surface(squeeze(x_cube_spher(:,:,idx_x)),squeeze(y_cube_spher(:,:,idx_x)),squeeze(z_cube_spher (:,:,idx_x)),db(abs(squeeze(img_3D_sector(:,:,idx_x)./max(max(img_3D_sector(:)))))));
axis image; title('3D','Color','white');
xlabel('x [mm]');ylabel('y [mm]');zlabel('z [mm]');
shading('flat');colormap gray; clim([-60 0]);view(40,0)
set(gca,'YColor','white','XColor','white','ZColor','white','Color','k');
set(gcf,'color','black'); set(gca,'Color','k')
f_handle = gcf;set(gca,'ZDir','reverse')
f_handle.InvertHardcopy = 'off';