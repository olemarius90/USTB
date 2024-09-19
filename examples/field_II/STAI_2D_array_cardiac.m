%% Computation of a STAI dataset with Field II using a 2D array and beamforming with USTB
%
% authors:  Ole Marius Hoel Rindal <olemarius@olemarius.net>
%           Stefano Fiorentini <stefano.fiorentini@ntnu.no>
%
% Last updated: 04.12.2023

clear all
close all
clc

%% basic constants
c0=1540;        % Speed of sound [m/s]
fs=100e6;       % Sampling frequency [Hz]
dt=1/fs;        % Sampling step [s]
downfact = 20;  % channel data downsampling factor;

f0 = 2.7e6;               % Transducer center frequency [Hz]
bw = 0.9;                 % Transducer bandwidth [1]
lambda = c0/f0;           % Wavelength [m]
Nc = 1.5;                 % Number of cycles
scat_density = 10;

%% Create idealized Matrix Array
kerf = lambda/5;

probe = uff.matrix_array();
probe.pitch_x = lambda/2;
probe.pitch_y = lambda/2;
probe.element_width = probe.pitch_x-kerf;
probe.element_height = probe.pitch_x-kerf;
probe.N_x = 40;
probe.N_y = 32;
probe.plot();

%% pulse definition
ti = -2/f0 : 1/fs : 2/f0;
impulse_response = 1e9 * gauspuls(ti, f0, bw);

te = 0:1/fs:(Nc/f0-1/fs);
excitation = square(2*pi*f0*te);
ir = conv(conv(impulse_response,excitation), impulse_response);
[maxVal , lag] = max(abs(hilbert(ir)));

% show the pulse to check that the lag estimation is on place (and that the pulse is symmetric)
figure()
hold on
plot((0:(length(ir)-1))/fs*1e6, ir, 'k')
plot((0:(length(ir)-1))/fs*1e6, abs(hilbert(ir)),'--r')
plot(lag/fs*1e6, maxVal,'g*')
hold off
box on
grid on
axis tight
xlabel('Fast time [\mus]')
ylabel('two Impulse response')
legend('RF','Envelope','Estimated lag');
title('2-way impulse response Field II');

%% Phantom
% === Define region filled with scatterers ===
zr = [3e-2, 8e-2];
xr = [-2.5e-2, 2.5e-2];
yr = [-2.5e-2, 2.5e-2];

% === estimate PSF size to calculate scatterer number
f_n = [3, 3.5];
xPSF = 1.2*f_n(1)*lambda;
yPSF = 1.2*f_n(2)*lambda;
zPSF = Nc*lambda/2;

% === Define geometry of nonechoic cyst ===
r = 1e-2;      % Radius of cyst [m]
xc = 0;        % Position of cyst in x [m]
yc = 0;
zc = 5.5e-2;     % Position of cyst in z [m]

N_scatterers = round(scat_density * diff(zr)*diff(xr)*diff(yr) / (pi/6*xPSF*yPSF*zPSF))
positions = (rand([N_scatterers, 3]) .* ...
        [diff(xr), diff(yr), diff(zr)]) + [xr(1), yr(1), zr(1)];

%Find the indexes inside cyst
positions(sqrt(sum((positions - [xc, yc, zc]).^2, 2)) < r, :) = [];
amplitudes = randn([size(positions, 1), 1]);

figure()
plot3(positions(:,1)*1e3,positions(:,2)*1e3,positions(:,3)*1e3,'b.')
hold on
axis equal
grid on
xlabel('x[mm]')
ylabel('y[mm]')
zlabel('z[mm]')

%% Compute STAI channel data
to = 2*zr(1)/c0 : 1/fs : 2*zr(2)/c0;

ch = complex(zeros([length(to(1:downfact:end)), probe.N_elements, probe.N_elements, 1], 'single'));   

% === Start Field II ===
field_init(0)
set_field('c',c0);              % Speed of sound [m/s]
set_field('fs',fs);             % Sampling frequency [Hz]
set_field('threads', 12)
set_field('show_times', 0)

noSubAz = round(probe.element_width/(lambda/8));        % number of subelements in the azimuth direction
noSubEl = round(probe.element_height/(lambda/8));       % number of subelements in the elevation direction
enabled = ones(probe.N_x, probe.N_y);

% === Define transmit aperture ===
Th = xdc_2d_array(probe.N_x, probe.N_y, probe.element_width, probe.element_height, ...
                    kerf, kerf, enabled, noSubAz, noSubEl, [0, 0, Inf]);

xdc_excitation (Th, excitation);
xdc_impulse(Th, impulse_response);

xdc_focus_times(Th, 0, zeros([1, probe.N_elements]));

% === Define receive aperture ===
Rh = xdc_2d_array(probe.N_x, probe.N_y, probe.element_width, probe.element_height, kerf, kerf, ...
                    enabled, noSubAz, noSubEl, [0, 0, Inf]);
xdc_apodization(Rh, 0, ones([1, probe.N_elements]));
xdc_focus_times(Rh, 0, zeros([1, probe.N_elements]));

h = waitbar(0, 'Simulating channel data...');

for n = 1:probe.N_elements
    waitbar(n/probe.N_elements, h, sprintf('Simulating channel data...Tx %d of %d', n, probe.N_elements))
    
    % === Set active element ===
    xdc_apodization(Th, 0, [zeros([1, n-1]), 1, zeros([1, probe.N_elements-n])]);

    % Generate channel data
    [scat, start_time] = calc_scat_multi(Th, Rh, positions, amplitudes);

    ti = start_time + (0:length(scat)-1)/fs - lag/fs;

    % Interpolation
    tmp = interp1(ti, scat, to, 'linear', 0);
        
    % === Hilbert transform-based demodulation and downsampling ===
    tmp = hilbert(tmp) .* exp(-1j*2*pi*f0*to(:));
    ch(:,:,n) = tmp(1:downfact:end,:);
end
    
% === Close Field II ===
close(h)
xdc_free(Th)
xdc_free(Rh)
field_end()


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
mid.dimension = dimension.both();
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
b_data_3D_linear = mid.go();
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