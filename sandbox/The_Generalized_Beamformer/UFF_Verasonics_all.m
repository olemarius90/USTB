%% Reading data from an UFF file recorded with the Verasonics CPWC_L7 example
clear all;
close all;

%% Checking the file is in the path
% data location
url='http://ustb.no/datasets/';      % if not found data will be downloaded from here

% all_filenames{1}='L7_FI_TheGB.uff'; selected_tx(1) = 20; tag{1} = 'FI'; tag_title{1} = 'Focused';
% all_filenames{2}='L7_DW_TheGB.uff'; selected_tx(2) = 1; tag{2} = 'DW'; tag_title{2} = 'Diverging';
% all_filenames{3}='L7_CPWC_TheGB.uff'; selected_tx(3) = 1; tag{3} = 'PW'; tag_title{3} = 'Plane';
% all_filenames{4}='L7_STA_TheGB.uff'; selected_tx(4) = 32; tag{4} = 'STA'; tag_title{4} = 'Single Element Diverging';
%%
all_filenames{1}='L7_DW_TheGB.uff'; selected_tx(1) = 1; tag{1} = 'DW'; tag_title{1} = 'Diverging';
all_filenames{2}='L7_CPWC_TheGB.uff'; selected_tx(2) = 1; tag{2} = 'PW'; tag_title{2} = 'Plane';

b_data_compare = uff.beamformed_data();

for f = 1:length(all_filenames)
    filename = all_filenames{f};
    tools.download(filename, url, data_path);
    channel_data=uff.read_object([data_path filesep filename],'/channel_data');
    channel_data.N_frames = 1; %Only reconstruct one frame
    channel_data.sound_speed = 1460;
   
    if contains(filename,'DW') || contains(filename,'FI')
        offset = eps;
        for seq = 1:channel_data.N_waves
            channel_data.sequence(seq).sound_speed = channel_data.sound_speed;
            %channel_data.sequence(seq).delay = channel_data.sequence(seq).delay + channel_data.sequence.t0_origin;
            channel_data.sequence(seq).delay
        end
    end

    scan=uff.linear_scan();
    scan.x_axis = linspace(channel_data.probe.x(1),channel_data.probe.x(end),512).';
    scan.z_axis = linspace(3e-3,50e-3,512).';
    scan.plot();
    title('Linear Scan');
    %saveas(gcf,'Figures/linear_scan.png')
    %%
    %
    % and beamform
    mid=midprocess.das();
    mid.dimension = dimension.receive;
    mid.channel_data=channel_data;
    mid.scan=scan;
    if contains(filename,'FI')
        MLA = scan.N_x_axis/channel_data.N_waves;
        mid.spherical_transmit_delay_model = spherical_transmit_delay_model.hybrid;
        mid.pw_margin = 2/1000;
        mid.transmit_apodization.window=uff.window.tukey25;
        mid.transmit_apodization.f_number = 2.5;
        mid.transmit_apodization.MLA = MLA;
        mid.transmit_apodization.MLA_overlap = MLA;
        mid.transmit_apodization.minimum_aperture = [2.5e-03 2.5e-03];
    else
        mid.transmit_apodization.window=uff.window.none;
        mid.transmit_apodization.f_number=1.7;
    end

    mid.receive_apodization.window=uff.window.hamming;
    mid.receive_apodization.f_number=1.7;
    b_data_single_tx=mid.go();

    b_data = uff.beamformed_data(b_data_single_tx);
    b_data.data = reshape(b_data.data,size(b_data.data,1),1,1,size(b_data.data,3)); %Hack to just plot one transmit, make the transmit Frames...

    %%
%     fig_handle = figure();clf;
%     channel_data.probe.plot(fig_handle);hold on;
%     channel_data.sequence(selected_tx(f)).source.plot(fig_handle);
%     view(30,30)
    %% Illustrate t0 compensation
    single_tx_images = b_data.get_image();

    fig = figure;clf;
    subplot(5,1,1);hold all;
    if contains(filename,'STA')
        active_elements = zeros(1,128);
        active_elements(selected_tx(f)) = 1;
        plot(active_elements,'*')
        xlim([1 128]);xlabel('Elements');
        xlabel('Elements');title('Transmitting Element');%ylabel(["Active Elements"])
    elseif contains(filename,'FI')
        %% Illustrate t0 compensation
        transmit_delays = channel_data.sequence(selected_tx(f)).delay_values;
        transmit_delays(50:end) = transmit_delays(49);

        plot(transmit_delays.*10^6,'HandleVisibility','off','LineWidth',2);hold on;
        plot(channel_data.sequence(selected_tx(f)).delay_values.*10^6,'HandleVisibility','off','LineWidth',2)
        plot(49,transmit_delays(49).*10^6,'o','DisplayName','Conventional t_0','LineWidth',2);hold on;
        plot(64,channel_data.sequence(selected_tx(f)).delay_values(end/2)*10^6,'o','DisplayName','Generalized t_0','LineWidth',2)
        xlim([1 128]);xlabel('Elements');
        xlabel('Elements');ylabel(["Delay [ms]"]);title('Transmit Wavefront Delay');
    else
        plot(channel_data.sequence(selected_tx(f)).delay_values*10^6,'HandleVisibility','off','LineWidth',2)
        plot([channel_data.sequence(selected_tx(f)).delay_values+abs(min(channel_data.sequence(selected_tx(f)).delay_values))]*10^6,'HandleVisibility','off','LineWidth',2)
        plot(channel_data.N_elements,[channel_data.sequence(selected_tx(f)).delay_values(end)+abs(min(channel_data.sequence(selected_tx(f)).delay_values))]*10^6,'o','DisplayName','Conventional t_0','LineWidth',2)
        plot(channel_data.N_elements/2,channel_data.sequence(selected_tx(f)).delay_values(end/2)*10^6,'ro','DisplayName','Generalized t_0','LineWidth',2)
        xlabel('Elements');ylabel(["Delay [ms]"]);title('Transmit Wavefront Delay');
    end
    set(gca,'FontSize',13)
    axis tight
    legend show

    b_data.plot(subplot(5,1,[2:5]),['Single ',tag_title{f},' Transmit Image'],[],[],[],[selected_tx(f)]);
    colorbar off
    set(gcf,'Position',[659 87 581 891])
    saveas(fig,['Figures/single_tx_illustrations/',tag{f},'.eps'],'eps2c')

    %%
    b_data_illustrate_tx_delay = uff.beamformed_data(b_data);
    b_data_illustrate_tx_delay.data = mid.transmit_delay(:,selected_tx(f))*10^6;
    b_data_illustrate_tx_delay.plot([],['Transmit Distance'],[],['none']);
    title(['Transmit Distance']);
    colormap jet;
    cb = colorbar(); 
    ylabel(cb,'Distance [mm]','FontSize',12)
    %%
    fig = figure()
    b_data_illustrate_rx_delay = uff.beamformed_data(b_data);
    b_data_illustrate_rx_delay.data = mid.receive_delay(:,64)*10^6;
    b_data_illustrate_rx_delay.plot(fig,['Receive  Distance'],[],['none']);
    title(['Receive  Distance']);
    colormap jet;
    cb = colorbar
    ylabel(cb,'Distance [mm]','FontSize',12)
    saveas(fig,['Figures/DelayCalculation/rx_64.png'])

    %% Coherent Compounding
    cc = postprocess.coherent_compounding();
    cc.input = b_data_single_tx;
    b_data_CC_postprocess = cc.go();

    % Compensate for different number of TX in area for FI
    if contains(filename,'FI')
        tx_comp = sum(mid.transmit_apodization.data,2);
        b_data_CC_postprocess.data = b_data_CC_postprocess.data.*(1./tx_comp);
    end
    fig = figure();
    b_data_CC_postprocess.plot(fig,['Coherent Compounded ',tag_title{f},' Transmit Image']);
    set(gcf,'Position',[638 434 602 544])
    saveas(fig,['Figures/compounded_images/',tag{f},'_compounded.eps'],'eps2c')

   % ['/beamformed_data_',tag{f},'_SOS_',num2str(channel_data.sound_speed)]
   % b_data_CC_postprocess.write(['./compare_transmit_types','_SOS_',num2str(channel_data.sound_speed),'.uff'],['/beamformed_data_',tag{f}]);

    b_data_compare.scan =  b_data_CC_postprocess.scan;
    b_data_compare.data(:,1,1,f) = b_data_CC_postprocess.data./max(b_data_CC_postprocess.data);
end
%%
%b_data_compare.data(:,1,1,2) = [];

fig = figure();
b_data_compare.plot(fig,['1=FI,2=DW,3=PW,4=STA']);
b_data_compare.save_as_gif('no_compensation.gif');