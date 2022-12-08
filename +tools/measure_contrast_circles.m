function [CR CNR GCNR CR_LC] = measure_contrast_circles(b_data, xc_nonecho, zc_nonecho, xc_echo, zc_echo, r, plot_flag, title_text, file_tag)
if nargin < 7
    plot_flag = 0;
end

if nargin < 8
    title_text = 'ROIs indicated';
end
%% Non echo Cyst contrast

xc_speckle = xc_nonecho;
zc_speckle = zc_nonecho;
if isa(b_data.scan,'uff.sector_scan')
    positions = reshape(b_data.scan.xyz,b_data.scan.N_depth_axis,b_data.scan.N_azimuth_axis,3);
else
    positions = reshape(b_data.scan.xyz,b_data.scan.N_z_axis,b_data.scan.N_x_axis,3);
end
points = ((positions(:,:,1)-xc_nonecho).^2) + (positions(:,:,3)-zc_nonecho).^2;
idx_cyst = (points < (r)^2);                     %ROI inside cyst
points = ((positions(:,:,1)-xc_echo).^2) + (positions(:,:,3)-zc_echo).^2;
idx_speckle = (points < (r)^2);       


%%
b_data.plot([],[title_text],[],[],[],[],'m','dark')
caxis([-50 0]);
axi = gca;
viscircles(axi,[xc_nonecho,zc_nonecho],r,'EdgeColor','r');
viscircles(axi,[xc_echo,zc_echo],r,'EdgeColor','b');
b_data.frame_rate = 1;
b_data.save_as_gif(['Figures/compare_ROI_',file_tag,'.gif'])

%%
img_signal = b_data.get_image('none');
%%
for f = 1:b_data.N_frames
    img_signal_current = img_signal(:,:,f);

    % Estimate the mean and the background of all images and calculate the CR
    mean_background = mean(abs(img_signal_current(idx_speckle(:))).^2)
    mean_ROI = mean(abs(img_signal_current(idx_cyst(:))).^2)
    
    CR(f) = 10*log10(mean_ROI/mean_background)

    sigma_background = std( abs(img_signal_current(idx_speckle(:))).^2 );
    sigma_ROI = std( abs(img_signal_current(idx_cyst(:))).^2 );
    
    CNR(f) = abs(mean_ROI - mean_background) / sqrt(sigma_ROI^2 + sigma_background^2)

    x=linspace(min(db(abs(img_signal_current(:)))),max(db(abs(img_signal_current(:)))),100);

    [pdf_i]=hist(db(abs(img_signal_current(idx_cyst(:)))),x);
    [pdf_o]=hist(db(abs(img_signal_current(idx_speckle(:)))),x);

%     mask_img_i = reshape(mask_i,size(img,1),size(img,2));
%     figure();
%     subplot(121);
%     imagesc(img.*mask_i);
%     subplot(122);
%     imagesc(img.*mask_o);

%% Plot probability density function

    fi = figure()
    plot(x,pdf_i./sum(pdf_i),'r-', 'linewidth',2); hold on; grid on;
    plot(x,pdf_o./sum(pdf_o),'b-', 'linewidth',2);
    hh=area(x,min([pdf_i./sum(pdf_i); pdf_o./sum(pdf_o)]), 'LineStyle','none');
    hh.FaceColor = [0.6 0.6 0.6];
    xlabel('dB(|s|)');
    ylabel('Probability');
    legend('p_i','p_o','OVL','Location','best');
    set(gca,'FontSize', 14);
    xlim([-150 0]);
    saveas(fi,['Figures/',file_tag,'_',num2str(f),'.png']);
    
    idx = find(min([pdf_i./sum(pdf_i); pdf_o./sum(pdf_o)])>0);
    xlim([-200 0])


    OVL =sum(min([pdf_i./sum(pdf_i); pdf_o./sum(pdf_o)]));
    MSR  = 1 - OVL/2;
    GCNR(f) = 1 - OVL;

    
%     if plot_flag
%         %%
%         figure;
%         subplot(211);
%         imagesc(db(abs(img_signal_current)).*idx_cyst)
%         caxis([-60 0]);
%         colorbar
%         subplot(212);
%         imagesc(db(abs(img_signal_current)).*idx_speckle)
%         caxis([-60 0]);
%         colorbar;
%     else
end
end

