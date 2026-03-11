filename = 'P4_2v_DW134847.uff';
DW_channel_data = uff.read_object([data_path filesep filename],'/channel_data');

%%
depth_axis=linspace(0,90e-3,512).';
angles_axis = linspace(deg2rad(-50),deg2rad(50),256);
scan=uff.sector_scan('azimuth_axis',angles_axis.','depth_axis',depth_axis);
    
figure;plot(DW_channel_data.sequence(1).delay_values)

das = midprocess.das();
das.dimension = dimension.receive();
das.channel_data = DW_channel_data;
das.scan = scan;
das.receive_apodization.window = uff.window.tukey25;
das.receive_apodization.f_number = 1.75;
das.transmit_apodization.window = uff.window.none;

b_data_dw = das.go()

%%
conv_transmit = DW_channel_data.sequence(1).delay_values*10^6-min(DW_channel_data.sequence(1).delay_values)*10^6;
[~,min_conv_t_0] = min(conv_transmit)
figure;
subplot(5,3,2);hold all;
plot(DW_channel_data.sequence(1).delay_values*10^6,'HandleVisibility','off')
plot(conv_transmit,'HandleVisibility','off')
plot(min_conv_t_0,conv_transmit(min_conv_t_0),'o','DisplayName','Conventional t_0');hold on;
plot(32,channel_data.sequence(1).delay_values(end/2)*10^6,'o','DisplayName','Generalized t_0')
xlabel('Elements');ylabel(["Delay [ms]"]);title('Tx Wavefront');
legend
set(gca,'FontSize',12)
ylim([-2 5])
b_data_dw.plot(subplot(5,3,[4:15]),['Single DW Image'])
cbar = colorbar
set(gcf,'Position',[680 300 750 678])
%set(cbar,'Visible','off')
%%
f = gcf;
saveas(f,'Figures/DW.png');