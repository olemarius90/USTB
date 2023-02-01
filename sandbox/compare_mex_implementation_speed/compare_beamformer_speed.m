%% CPWC simulation to compare speeds of the various USTB beamformers.
%
% In this example, we conduct a simple simulation to compare the speeds
% achieved with USTB's:
%
% # MEX CUDA beamformer
% # MEX FAST CPU beamformer
%
% This tutorial assumes familiarity with the contents of the
% <./CPWC_linear_array.html 'CPWC simulation with the USTB built-in Fresnel
% simulator'> tutorial. Please feel free to refer back to that for more
% details.
%
% Alfonso Rodriguez-Molares <alfonso.r.molares@ntnu.no>
% Stefano Fiorentini <stefano.fiorentini@ntnu.no>
% Last edited: 24.11.2022

clear all
close all
clc

testCase = 1; % 0 = linear/planewave, 1 = sector/focused
do_demodulation = true;
nFrames = 100;

%% Phantom
switch testCase
    case 0
        x_sca=[zeros(1,7), -15e-3:5e-3:15e-3];
        z_sca=[5e-3:5e-3:35e-3, 20e-3*ones(1,7)];
    case 1
        x_sca=[zeros(1,7), -45e-3:15e-3:45e-3];
        z_sca=[10e-3:15e-3:100e-3, 65e-3*ones(1,7)];
end

N_sca=length(x_sca);
pha=uff.phantom();
pha.sound_speed=1540;            % speed of sound [m/s]
pha.points=[x_sca.', zeros([N_sca,1]), z_sca.', ones([N_sca,1])];    % point scatterer position [m]

%% Probe
switch testCase
    case 0
        prb=uff.linear_array();
        prb.N=128;                  % number of elements
        prb.pitch=300e-6;           % probe pitch in azimuth [m]
        prb.element_width=270e-6;   % element width [m]
        prb.element_height=5e-3;    % element height [m]
    case 1
        prb=uff.linear_array();
        prb.N=90;                   % number of elements
        prb.pitch=220e-6;           % probe pitch in azimuth [m]
        prb.element_width=200e-6;   % element width [m]
        prb.element_height=8e-3;    % element height [m]

end

%% Pulse

switch testCase

    case 0
        pul=uff.pulse();
        pul.center_frequency=7e6;       % transducer frequency [MHz]
        pul.fractional_bandwidth=0.8;   % fractional bandwidth [unitless]

    case 1
        pul=uff.pulse();
        pul.center_frequency=4e6;       % transducer frequency [MHz]
        pul.fractional_bandwidth=0.8;   % fractional bandwidth [unitless]
end
%% Sequence generation

switch testCase
    case 0
        nTx=15;
        angles=linspace(-10, 10, nTx)/180*pi;
        seq=uff.wave();

        for n=1:nTx
            seq(n)=uff.wave();
            seq(n).probe=prb;
            seq(n).source.azimuth=angles(n);
            seq(n).wavefront = uff.wavefront.plane;
            seq(n).sound_speed=pha.sound_speed;
        end
    case 1
        nTx=192; % n transmits
        F = 6.5e-2; % focus speed
        angles=linspace(-40, 40, nTx)/180*pi;
        seq=uff.wave();

        for n=1:nTx
            seq(n)=uff.wave();
            seq(n).probe=prb;
            seq(n).source.azimuth=angles(n);
            seq(n).source.distance = F;
            seq(n).wavefront = uff.wavefront.spherical;
            seq(n).sound_speed=pha.sound_speed;
        end

    otherwise
        error("Case not supported")
end

%% Fresnel simulator
sim=fresnel();
sim.phantom=pha;                % phantom
sim.pulse=pul;                  % transmitted pulse
sim.probe=prb;                  % probe
sim.sequence=seq;               % beam sequence
sim.sampling_frequency=50e6;    % sampling frequency [Hz]

tic()
channel_data=sim.go();
toc()

if do_demodulation
    demod = preprocess.fast_demodulation();
    demod.input = channel_data;
    demod.modulation_frequency = pul.center_frequency;
    demod.downsample_frequency = 2*pul.center_frequency;

    channel_data = demod.go();
end

%% Scan
switch testCase
    case 0
        scan = uff.linear_scan('x_axis',linspace(-2e-2,2e-2,256).', 'z_axis', linspace(0, 4e-2, 512).');
    case 1
        scan = uff.sector_scan('azimuth_axis',linspace(angles(1),angles(end),256).', 'depth_axis', linspace(0, 12e-2, 768).');
end

%% Pipeline

pipe=pipeline();
pipe.channel_data=channel_data;
pipe.scan=scan;

switch testCase
    case 0
        pipe.receive_apodization.window=uff.window.hamming;
        pipe.receive_apodization.f_number=2;

        pipe.transmit_apodization.window=uff.window.hamming;
        pipe.transmit_apodization.f_number=2;
    case 1
        pipe.receive_apodization.window=uff.window.none;
        pipe.receive_apodization.f_number=3.5;
        pipe.receive_apodization.maximum_aperture = 2e-2 * pipe.receive_apodization.f_number(1)^2;

        pipe.transmit_apodization.window=uff.window.hamming;
        pipe.transmit_apodization.f_number=3.5;
        pipe.transmit_apodization.minimum_aperture = 5e-3 * pipe.transmit_apodization.f_number(1)^2;
end

proc            = midprocess.das();
proc.code       = code.mex();
proc.dimension  = dimension.both;
fprintf(1, 'Precalculating apodization\n')
pipe.go({proc});

%% Test beamforming speed

dOp_per_frame = 2*scan.N_pixels*channel_data.N_channels*channel_data.N_waves;

das_mexFast_time = zeros([length(nFrames), 1]);
das_mex_gpu_time = zeros([length(nFrames), 1]);

if isscalar(nFrames)
    profile on
end

for n=1:length(nFrames)
    % replicate frames

    channel_data.data=repmat(channel_data.data(:,:,:,1),[1, 1, 1, nFrames(n)]);

    % Time USTB's MEX GPU tex 2D implementation
    proc            = midprocess.das();
    proc.code       = code.mex_gpu_tex2D;
    proc.gpu_device = 0;
    proc.dimension  = dimension.both;
    fprintf(1, 'Processing %d frames: MEX CUDA tex 2D\n', nFrames(n))
    tic()
    bf_data_mex_gpu = pipe.go({proc});
    das_mex_gpu_time(n) = toc();

    % Time USTB's MEX FAST implementation
    proc            = midprocess.das();
    proc.code       = code.mexFast;
    proc.dimension  = dimension.both;
    fprintf(1, 'Processing %d frames: MEX FAST\n', nFrames(n))
    tic()
    bf_data_mexFast_cpu = pipe.go({proc});
    das_mexFast_time(n) = toc();
end

if isscalar(nFrames)
    profile off
    profile viewer
end
%% Plot the images for visual inspection of the results
switch class(scan)
    case "uff.linear_scan"
        dim = [scan.N_z_axis, scan.N_x_axis];

    case "uff.sector_scan"
        dim = [scan.N_depth_axis, scan.N_azimuth_axis];

end

X = reshape(scan.x, dim);
Y = reshape(scan.y, dim);
Z = reshape(scan.z, dim);


figure('Color', 'white')
tiledlayout(1, 2, "TileSpacing", "compact", "Padding", "compact")
hAx(1) = nexttile();
surface(X*1e2, Y*1e2, Z*1e2, ...
    20*log10(abs(reshape(bf_data_mex_gpu.data(:,end), dim)) / ...
    max(abs(bf_data_mex_gpu.data(:,end)))), "LineStyle", "none")
clim([-60, 0])
view([0, 0])
set(gca, "ZDir", "reverse")
grid on
box on
axis equal tight
xlabel("x [cm]")
ylabel("z [cm]")
title("mex CUDA tex 2D")

hAx(2) = nexttile();
surface(X*1e2, Y*1e2, Z*1e2, ...
    20*log10(abs(reshape(bf_data_mexFast_cpu.data(:,end), dim)) / ...
    max(abs(bf_data_mexFast_cpu.data(:,end)))), "LineStyle", "none")
clim([-60, 0])
view([0, 0])
set(gca, "ZDir", "reverse")
grid on
box on
axis equal tight
xlabel("x [cm]")
ylabel("z [cm]")
title("mex FAST")

linkaxes(hAx)
linkprop(hAx, {'CameraPosition','CameraUpVector'});

diff = mean(abs(bf_data_mex_gpu.data(:,1) - bf_data_mexFast_cpu.data(:,1)), 1)

%% Plot the runtimes
figure('Color', 'white')

cMap = lines(2);

hold on
plot(nFrames(1:n)*dOp_per_frame/1e9,das_mex_gpu_time(1:n),'s-','linewidth',1.5,'color',cMap(1,:));
plot(nFrames(1:n)*dOp_per_frame/1e9,das_mexFast_time(1:n),'o-','linewidth',1.5,'color',cMap(2,:));
hold off

for nn=1:length(nFrames)
    text(nFrames(nn)*dOp_per_frame/1e9,das_mex_gpu_time(nn),sprintf('%0.2f s', das_mex_gpu_time(nn)), ...
        'horizontalalignment', 'left', 'verticalalignment', 'top','color',cMap(1,:),'fontweight','bold');

    text(nFrames(nn)*dOp_per_frame/1e9,das_mexFast_time(nn),sprintf('%0.2f s', das_mexFast_time(nn)), ...
        'horizontalalignment', 'right', 'verticalalignment', 'bottom','color',cMap(2,:),'fontweight','bold');
end

grid on
box on

legend('MEX CUDA tex 2D', 'MEX FAST', 'Location','Best');
xlabel('Delay operations [1e9]');
ylabel('Elapsed time [s]');



