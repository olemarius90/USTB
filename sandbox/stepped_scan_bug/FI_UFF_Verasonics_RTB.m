%% Reading FI data from an UFF file recorded from a Verasonics Scanner
%
% In this example we show how to read channel data from a
% UFF (Ultrasound File Format) file recorded with a Verasonics scanner.
% You will need an internet connection to download data.
%
% _by Ole Marius Hoel Rindal <olemarius@olemarius.net>
%   and Alfonso Rodriguez-Molares <alfonso.r.molares@ntnu.no>
%
%   $Last updated: 2017/10/06$

%% Checking the file is in the path
%
% To read data from a UFF file the first we need is, you guessed it, a UFF
% file. We check if it is on the current path and download it from the USTB
% websever.

clear all; close all;

% data location
url='http://ustb.no/datasets/';      % if not found downloaded from here
filename='L7_FI_IUS2018.uff';
tools.download(filename, url, data_path);   

channel_data=uff.read_object([data_path filesep filename],'/channel_data');
channel_data_corrected=uff.read_object([data_path filesep filename],'/channel_data');
%% channel_data_corrected = copy(channel_data);
%% Set the origin in the channel data sequence to the focal point. Should
%this be default when storing the data in the Verasonics class?
for w = 1:channel_data.N_waves
    channel_data_corrected.sequence(w).origin.x = channel_data.sequence(w).source.x;
end

%% Define Scan
MLA = 4;
x_axis=zeros(channel_data.N_waves,1);
for n=1:channel_data.N_waves
    x_axis(n)=channel_data.sequence(n).source.x;
end
z_axis=linspace(5e-3,50e-3,512*2).';
scan_RTB = uff.linear_scan('x_axis',linspace(x_axis(1),x_axis(end),length(x_axis)*MLA)','z_axis',z_axis);

%% Retrospective beamforming
mid_RTB=midprocess.das();
mid_RTB.dimension = dimension.receive();
mid_RTB.channel_data=channel_data;
mid_RTB.scan=scan_RTB;
mid_RTB.spherical_transmit_delay_model = spherical_transmit_delay_model.hybrid;
mid_RTB.transmit_apodization.window=uff.window.tukey25;
mid_RTB.transmit_apodization.f_number = 2;
mid_RTB.transmit_apodization.MLA = MLA;
mid_RTB.transmit_apodization.MLA_overlap = MLA;
mid_RTB.transmit_apodization.minimum_aperture = [3.0000e-03 3.0000e-03];
mid_RTB.receive_apodization.window=uff.window.boxcar;
mid_RTB.receive_apodization.f_number=1.7;

b_data_RTB=mid_RTB.go();
b_data_RTB.plot([])
b_data_RTB.frame_rate = 20;
b_data_RTB.save_as_gif("example_transmit_apod_error.gif")

%%
mid_RTB.channel_data=channel_data_corrected;
b_data_RTB_corrected=mid_RTB.go();
b_data_RTB_corrected.plot([])
b_data_RTB_corrected.frame_rate = 20;
b_data_RTB_corrected.save_as_gif("example_transmit_apod_corrected.gif")