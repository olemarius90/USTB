
%% Define the parameters

c0=1540;     % Speed of sound [m/s]
fs = 20.867e6;    % DSB L7 probe parameters
dt=1/fs;     % Sampling step [s]

%% field II initialisation

field_init(0);
set_field('c',c0);              % Speed of sound [m/s]
set_field('fs',fs);             % Sampling frequency [Hz]
set_field('use_rectangles',1);  % use rectangular elements

%% Transducer definition L11-4v, 128-element linear array transducer
probe = uff.linear_array();
f0                      = 5e+06;     
lambda                  = c0/f0;           % Wavelength [m]
probe.pitch             = 0.1*1e-3 ;

probe.element_width     = lambda/4 ;  kerf = probe.pitch - probe.element_width ;
probe.N                 = 128 ;             % Number of elements
pulse_duration          = 2.5 ;             % pulse duration [cycles]

%% Pulse definition

pulse = uff.pulse();
pulse.center_frequency = f0;
pulse.fractional_bandwidth = 0.65;        % probe bandwidth [1]
t0 = (-1/pulse.fractional_bandwidth/f0): dt : (1/pulse.fractional_bandwidth/f0);
impulse_response = 1e6*gauspuls(t0, f0, pulse.fractional_bandwidth);
impulse_response = impulse_response-mean(impulse_response); % To get rid of DC

te = (-pulse_duration/2/f0): dt : (pulse_duration/2/f0);
excitation = square(2*pi*f0*te+pi/2)*1e20;
one_way_ir = conv(impulse_response,excitation);
two_way_ir = conv(one_way_ir,impulse_response);
lag = length(two_way_ir)/2+1;

%% Aperture Objects

noSubAz=round(probe.element_width/(lambda/8));        % number of subelements in the azimuth direction
noSubEl=round(probe.element_height/(lambda/8));       % number of subelements in the elevation direction
Th = xdc_linear_array (probe.N, probe.element_width, probe.element_height, kerf, noSubAz, noSubEl, [0 0 Inf]);
Rh = xdc_linear_array (probe.N, probe.element_width, probe.element_height, kerf, noSubAz, noSubEl, [0 0 Inf]);

% We also set the excitation, impulse response and baffle as below:
xdc_excitation (Th, excitation);
xdc_impulse (Th, impulse_response);
% xdc_baffle(Th, 0);
xdc_center_focus(Th,[0 0 0]);
xdc_impulse (Rh, impulse_response);
% xdc_baffle(Rh, 0);
xdc_center_focus(Rh,[0 0 0]);

%% Phantom to check angular coverage

angle = (-60:20:60)'*pi/180 ;
dist = 15*1e-3 ;
point_position = dist.*[sin(angle), zeros(size(angle)), cos(angle)] ;

angle = (-40:20:40)'*pi/180 ;
dist = 25*1e-3 ;

point_position = [point_position; dist.*[sin(angle), zeros(size(angle)), cos(angle)]] ;

point_position = [point_position; [zeros(4, 1), zeros(4, 1), (30:10:60).']*1e-3] ;

point_amplitude = ones(size(point_position,1),1);


%% Output data
% We define the variables to store our output data
cropat=round(1.1*2*sqrt((max(point_position(:,1))-min(probe.x))^2+max(point_position(:,3))^2)/c0/dt);   % maximum time sample, samples after this will be dumped

%% Compute STA signals

disp('Field II: Computing STA dataset');
wb = waitbar(0, 'Field II: Computing STA dataset');
STA=zeros(cropat,probe.N,probe.N);    % impulse response channel data
clear seq
for n=1:probe.N
    waitbar(n/probe.N, wb);

    % transmit aperture
    xdc_apodization(Th, 0, [zeros(1,n-1) 1 zeros(1,probe.N-n)]);
    xdc_focus_times(Th, 0, zeros(1,probe.N));

    % receive aperture
    xdc_apodization(Rh, 0, ones(1,probe.N));
    xdc_focus_times(Rh, 0, zeros(1,probe.N));

    % do calculation
    [v,t]=calc_scat_multi(Th, Rh, point_position, point_amplitude);

    % build the dataset
    STA(1:size(v,1),:,n)=v;

    % Sequence generation
    seq(n)=uff.wave();
    seq(n).probe=probe;
    seq(n).source.xyz=[probe.x(n) probe.y(n) probe.z(n)];
    seq(n).sound_speed=c0;
    seq(n).delay = probe.r(n)/c0-lag*dt+t; % t0 and center of pulse compensation

    fprintf("Transmit no %2.0f is finished \n", n)
end
close(wb);

%% Channel Data

channel_data = uff.channel_data();
channel_data.sampling_frequency = fs;
channel_data.sound_speed = c0;
channel_data.initial_time = 0;
channel_data.modulation_frequency = 0 ;
channel_data.pulse = pulse;
channel_data.probe = probe;
channel_data.sequence = seq;
channel_data.data = STA ;

%% Reconstruction scan

scan=uff.linear_scan('x_axis', linspace(-30e-3, 30e-3, 512).', 'z_axis', linspace(0e-3, 65e-3, 512).');

%% Demodulation
demod = preprocess.demodulation();
demod.input = channel_data;
channel_data_demod = demod.go();

%% Wavenumber

mid = midprocess.Fourier_beamforming ;

mid.channel_data  = channel_data_demod ;
mid.scan          = scan ;
mid.refocus       = false ;
mid.spatial_padding = 3 ;
mid.temporal_padding = 3 ;
mid.DAS_consistent = false ;
mid.USTB_scan = false ;
mid.temp_origin = 40e-3 ;
mid.angle_apodization = 60 ;

w_data = mid.go() ;

w_data.plot([], w_data.name, []) ;
wavenumber_image = w_data.get_image('abs') ;
wavenumber_image = wavenumber_image/max(wavenumber_image(:)) ;

% % % % % % % % % % % %
% When reconstructing images using the Wavenumber algorithm and DAS
% wavenumber algorithm, grating lobes in channel data were also
% reconstructed.
% Till now (01/12/2025), midprocess files Wavenumber_algorithm.m and
% DAS_wavenumber_algorithm.m do not reconstruct grating lobes in final
% images. To reconstruct grating lobes, rotational interpolation is
% required, this can be done by simply repeating the frequency domain
% channel data in ku and kv domain.

%% DAS consistent wavenumber algorithm

mid = midprocess.Fourier_beamforming ;

mid.channel_data  = channel_data_demod ;
mid.scan          = scan ;
mid.refocus       = false ;
mid.spatial_padding = 3 ;
mid.temporal_padding = 3 ;
mid.DAS_consistent = true ;
mid.USTB_scan = false ;
mid.temp_origin = 40e-3 ;
mid.angle_apodization = 60 ;
DCWA_data = mid.go() ;

DCWA_data.plot([], DCWA_data.name, []) ;
DCWA_image = DCWA_data.get_image('abs') ;
DCWA_image = DCWA_image/max(DCWA_image(:)) ;

%% DAS method

mid = midprocess.das ;
mid.channel_data=channel_data_demod;
mid.scan = w_data.scan;
apod_angle = 60 ;
f_number = cot(deg2rad(apod_angle))/2 ;
mid.transmit_apodization.f_number = f_number ;
mid.transmit_apodization.window = uff.window.tukey25 ;
mid.transmit_apodization.minimum_aperture = [0 0] ;
mid.transmit_apodization.maximum_aperture = [1e5 1e5] ;
mid.receive_apodization.f_number = f_number ;
mid.receive_apodization.window = uff.window.tukey25 ;
mid.receive_apodization.minimum_aperture = [0 0] ;
mid.receive_apodization.maximum_aperture = [1e5 1e5] ;

b_data=mid.go() ;
b_data.plot([],'DAS',[] );
das_image = b_data.get_image('abs') ;
das_image = das_image/max(das_image(:)) ;

%% Final image

scan = w_data.scan ;
Y_lim = [0 55] ;
X_lim = [-20 20] ;

figure;
t1 = tiledlayout(1,3,'TileSpacing','Compact','Padding','Compact');

nexttile
imagesc(scan.x_axis*1e3, scan.z_axis*1e3, db(wavenumber_image), [-60, 0])
% colorbar
colormap gray
axis equal tight
title("WA",'fontsize', 30)
ylim(Y_lim)
xticks([-15, 0, 15])
xlim(X_lim)

nexttile
imagesc(scan.x_axis*1e3, scan.z_axis*1e3, db(das_image), [-60, 0])
% colorbar
colormap gray
axis equal tight
title("DAS",'fontsize', 30)
yticklabels([])
xticks([-15, 0, 15])
xtickformat() ;
ylim(Y_lim)
xlim(X_lim)
xlabel(gca, 'x [mm]')


nexttile
imagesc(scan.x_axis*1e3, scan.z_axis*1e3, db(DCWA_image), [-60, 0])
% colorbar
colormap gray
axis equal tight
title("DCWA",'fontsize', 30 )
yticklabels([])
xticks([-15, 0, 15])
clr_bar = colorbar ;
ylim(Y_lim)
xlim(X_lim)

fontsize(20,'points')
% title(ax(1), 'Wavenumber algorithm','fontsize', 30)
% set(gca,'FontWeight','bold')

yl = ylabel(t1, 'y [mm]') ;

clr_bar.Layout.Tile = "south" ;
clr_bar.TickDirection = "out" ;

%% Legends
legends_(1) = "WA" ;
legends_(2) = "DAS" ;
legends_(3) = "DCWA" ;

%% Vertical drop

Vertical_points = find(point_position(:, 1)==0) ;
N_points = length(Vertical_points) ;

Lat_dist = 0 ;
lat_index = get_index(Lat_dist, scan.x_axis);
Thickness_lat = 0.4e-3 ; %% in mm
Thickness_idx_lat = round(Thickness_lat/scan.x_step) ;

Axi_dist = point_position(:, 3) ;
Axi_index = get_index(Axi_dist, scan.z_axis);
Thickness_axi = 0.4e-3 ; %% in mm
Thickness_idx_axi = round(Thickness_axi/scan.z_step) ;

Full_image = cat(3, wavenumber_image, das_image, DCWA_image);

point_val = zeros(N_points, size(Full_image, 3)) ;
for ii = 1:N_points
    index_z = Vertical_points(ii) ;
    temp_mat = Full_image(Axi_index(index_z)+(-Thickness_idx_axi:Thickness_idx_axi), lat_index+(-Thickness_idx_lat:Thickness_idx_lat), :) ;
    point_val(ii, :) = max(temp_mat, [], [1, 2]) ;
end

z_ax = 1e3*point_position(Vertical_points, 3) ;

figure
plot(z_ax, point_val(:, 1), LineWidth=2, DisplayName=legends_(1), Color=[255 128 0]/255, Marker="o", MarkerSize=10)
hold on
plot(z_ax, point_val(:, 2), LineWidth=5, DisplayName=legends_(2), Color=0.8*[1 1 1], Marker="*", MarkerSize=15)
plot(z_ax, point_val(:, 3), LineWidth=2, DisplayName=legends_(3), Color='k', Marker="diamond", MarkerSize=10, LineStyle="--")

xlabel('Axial Distance [mm]')
ylabel('Normalized Amplitude')

fontsize(gcf, 20, 'points')
legend("Location","northeast", Box="off")
xlim([15 50])
% set(gcf, "Position", [589,469,715,481]) ;
set(gcf, "Position", [5.826e+02,5.098e+02,7.214e+02,4.402e+02]) ;


%% Angular drop
dist = 15e-3 ;
angle = (-60:20:60)'*pi/180 ;

% angle = (-40:20:40)'*pi/180 ;
% dist = 25*1e-3 ;
new_point_position = (dist.*[sin(angle), zeros(size(angle)), cos(angle)]) ;

N_points = size(new_point_position, 1) ;

Ax_dist = new_point_position(:, 1) ; 
Ax_index = get_index(Ax_dist, scan.x_axis);
Thickness_Ax = 0.4e-3 ; %% in mm
Thickness_idx_Ax = round(Thickness_Ax/scan.x_step) ;

Lat_dist = new_point_position(:, 3) ; 
Lat_index = get_index(Lat_dist, scan.z_axis);
Thickness_Lat = 0.4e-3 ; %% in mm
Thickness_idx_Lat = round(Thickness_Lat/scan.z_step) ;

Full_image = cat(3, wavenumber_image, das_image, DCWA_image);

point_val = zeros(N_points, size(Full_image, 3)) ;
for ii = 1:N_points
    temp_mat = Full_image(Lat_index(ii)+(-Thickness_idx_Lat:Thickness_idx_Lat), Ax_index(ii)+(-Thickness_idx_Ax:Thickness_idx_Ax), :) ;
    point_val(ii, :) = max(temp_mat, [], [1, 2]) ;
end

angle = rad2deg(angle) ;

figure
% subplot(2, 1, 2)
plot(angle, point_val(:, 1), LineWidth=2, DisplayName=legends_(1), Color=[255 128 0]/255, Marker="o", MarkerSize=10)

hold on
plot(angle, point_val(:, 2), LineWidth=5, DisplayName=legends_(2), Color=0.8*[1 1 1], Marker="*", MarkerSize=15)
plot(angle, point_val(:, 3), LineWidth=2, DisplayName=legends_(3), Color='k', Marker="diamond", MarkerSize=10, LineStyle="--")

xlabel('Angle [Degrees]')
ylabel('Normalized Amplitude')
% title('Amplitude comparison over angles at depth=15mm','FontSize', 30)
fontsize(gcf, 20, 'points')
legend("Location","south", Box="off")
% xlim([10 55])
set(gcf, "Position", [5.826e+02,5.098e+02,7.214e+02,4.402e+02]) ;


%% This function finds the index of closest number to 'dist' in 'in_axis'
function out_index = get_index(dist, in_axis, out_range)
if nargin<3
    out_range = 0;
end
if isscalar(dist)
    if dist>max(in_axis)||dist<min(in_axis)
        if out_range==0
            error("Querry value is not in the range of sample points")
        end
    end
    [~, out_index] = min(abs(dist - in_axis(:))) ;
else
    size_dist = size(dist) ;
    dist = reshape(dist, 1, []) ;
    in_axis = in_axis(:) ;
    if max(dist)>max(in_axis)||min(dist)<min(in_axis)
        if out_range==0
            error("Some of the querry values are not in the range of sample points. " + ...
                "If you want to proceed anyway use 1 as a third argument while calling the function")
        end
    end
    [~, out_index] = min(abs(dist - in_axis)) ;
    out_index =  out_index(:) ;
    out_index = reshape(out_index, size_dist) ;
end
end



