% -------------------------------------------------------------------------
% Title: fieldmap_compare.m
% Author: Bernhard Gruber
%
% Purpose:
%
% Versionhistory:
%   12/2020 -
%   06/2021 - Adapt to B0 field map analyization
%
% Notes:
%   cite: Fa-Hsuan Lin, "Magnetic field by Biot-Savart's Law"
%         http://maki.bme.ntu.edu.tw/?page_id=333
% -------------------------------------------------------------------------
% Units of T, V/m, and S/m
%
% SNR = s'*s/sqrt(s'*N*s);
% -------------------------------------------------------------------------
clear all, close all, clc

curr_path = pwd;
outdir=[curr_path,'\Output\'];
%outdir = [uigetdir(),'\'];

addpath('DATA\');
addpath('C:\Users\BGruber\Documents\!Work\02.2 - MUW\02_Projects\Measurements\3TPrisma_20042021\s4_dats');
addpath('mapVBVD\');

%% Load Data
os = 2; %oversampling factor

%SNR_noSHIM = mapVBVD('DATA\meas_MID01258_FID03672_sn_noshim.dat');
%SNR_SHIM = mapVBVD('DATA\meas_MID01250_FID03664_sn_shim_0A.dat');
%FM_0A = mapVBVD('DATA\meas_MID01254_FID03668_fieldmap_shim_0A.dat');
%FM_IC_1A = mapVBVD('DATA\meas_MID01256_FID03670_fieldmap_shim_IC_1A_cableaway.dat');
%FM_OC_1A = mapVBVD('DATA\meas_MID01257_FID03671_fieldmap_shim_OC_1A_cableaway.dat');

for setup=1:1%3
    switch setup
        case 1
            twix = mapVBVD('s4_dats\meas_MID01113_FID22284_ke_gre_aspireTest_B0map_3D_fullFOV_inhaled.dat');
            twix=twix{2};
%            twix.image.dataSize(1:11)
        case 2
            twix = mapVBVD('DATA\meas_MID01256_FID03670_fieldmap_shim_IC_1A_cableaway.dat');
            twix=twix{2};
%            twix.image.dataSize(1:11)
        case 3
            twix = mapVBVD('DATA\meas_MID01257_FID03671_fieldmap_shim_OC_1A_cableaway.dat');
            twix=twix{2};
%            twix.image.dataSize(1:11)
    end
    twix.image.flagDoAverage = true; % gets rid of the average dimension by automatically averaging them during the read procedure
    twix.image.flagRemoveOS = true;


    %% Data reconstruction (fft)
    %img = zeros([twix.image.NCol/os, twix.image.NLin, twix.image.NSli, twix.image.NEco], 'single');    
    pha_sc = zeros([twix.image.NCol/os, twix.image.NCha/2, twix.image.NLin, 1+0*twix.image.NSli, twix.image.NEco], 'single'); %single channel phases
    img_sc = pha_sc;
    TE = single(twix.hdr.Meas.alTE(1:twix.image.NEco)/1000000);
    disp(['setup = ',num2str(setup)]);    
    for echo=1:2%twix.image.NEco
        disp(['echo = ',num2str(echo)]);
        for sli = 1:2%1:twix.image.NSli
            disp(['slice = ',num2str(sli)]);
            tic
            data = twix.image(:,2:2:4,:,1,sli,1,1,echo); % read in the data for slice 'sli'
            data(1:5,:,:)=0;
            data(124:128,:,:)=0;
            fft_dims = [1 3]; % fft in col and lin dimension
            for f = fft_dims
                data= ifftshift(ifft(fftshift(data,f),[],f),f);
            end
            img_sc(:,:,:,sli,echo)=abs(data);
            pha_sc(:,:,:,sli,echo)=sunwrap(data);
            toc
        end
    end
    
    img_sc=permute(img_sc,[1 3 2 4 5]);
    pha_sc=permute(pha_sc,[1 3 2 4 5]);

    %obtain b0 map from combined phase image
    tmp=0;
    for ch=1:2
        tmp=tmp+img_sc(:,:,ch,:,2).*img_sc(:,:,ch,:,1).*exp(1i*(pha_sc(:,:,ch,:,2)-pha_sc(:,:,ch,:,1))); 
    end
    disp('unwrapping combined b0 map');
    tic
    b0(:,:,:,setup)=sunwrap(squeeze(tmp))/(2*pi*(TE(2)-TE(1)));
    toc
    %obtain combined magnitude image   
    img_combined(:,:,:,:,setup)=squeeze(sum(img_sc.^2,3));
end
%%
disp('subtracting baseline b0 map');
b0_IC=squeeze(b0(:,:,:,2)-b0(:,:,:,1));
b0_OC=squeeze(b0(:,:,:,3)-b0(:,:,:,1));
%%
disp('applying magnitude treshold');
magnitude_threshold=0.02;
mask=ones(size(squeeze(img_combined(:,:,:,1,1))));
mask(img_combined(:,:,:,1)<magnitude_threshold*max(img_combined(:)))=NaN;

figure(4);
slice=5;
for setup=1:3
    subplot(1,3,setup);
    imagesc(mask(:,:,slice).*squeeze(img_combined(:,:,slice,setup)));axis equal tight square;colormap gray;
end

figure(5);
subplot(2,3,5); imagesc(mask.*b0_IC,[-300 300]);axis equal tight square;
subplot(2,3,6); imagesc(mask.*b0_OC,[-300 300]);axis equal tight square;

for setup=1:3
    subplot(2,3,setup);
    imagesc(mask.*squeeze(b0(:,:,:,setup)),[-300 300]);axis equal tight square;
end

%%
%     figure(setup);
%     sli=1;
%     for ch=1:2
%         subplot(2,5,5*(ch-1)+1);imagesc(squeeze(img_sc(:,:,ch,sli,1)));axis equal square tight off;colormap gray;
%         for echo=1:3
%             subplot(2,5,5*(ch-1)+1+echo);imagesc(squeeze(pha_sc(:,:,ch,sli,echo)));axis equal square tight off;colormap gray;
%         end
%         subplot(2,5,5*(ch-1)+5);imagesc(squeeze(b0(:,:,ch,sli,setup)));axis equal square tight off;colormap gray;
%     end

% figure(300)
% imagesc(rot90(img(:,:)), [0, 0.7*max(img(:))]);
% colormap jet, h=colorbar(),set(h,'FontSize',14), axis image off;
% title('field map without applied SHIM');

return
%% Process Field Maps

nx=256; ny=256; nz=53;    % number of field maps points along x, y, and z
r         = 38;           % set coil radius
fm_dim    = [256 256 90]; % size of field map FOV along x, y, and z

b1_z_scaled = -b1_z.*42.57e6; % scale the field map

%%%%%%%%%%%%%%%%%%%%%%%% Create MASK for phantom %%%%%%%%%%%%%%%%%%%%%%
imageSizeX = 128;   imageSizeY = 128;   imageSizeZ = 45;
[columnsInImage rowsInImage pagesInImage] = meshgrid(1:imageSizeX, 1:imageSizeY,1:imageSizeZ);

centerX = 70; centerY = 78.5;   centerZ = 45;   radius = 33.5; % EDIT

cylinderVoxels = (rowsInImage - centerY).^2 + (columnsInImage - centerX).^2 <= radius.^2;
mask_double = double(cylinderVoxels);

Mask_new = imerode(mask_double, [0 1 0; 1 1 1; 0 1 0]);
Mask_new = imerode(Mask_new, [0 1 0; 1 1 1; 0 1 0]);
Mask_new = imerode(Mask_new, [0 1 0; 1 1 1; 0 1 0]);
Mask_nan = Mask_new; %white background
Mask_nan(Mask_new==0)=nan; %white background

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear temp_row indices
temp_plot = []; temp_row = [];
plot_range = [-80 +80]; % define the plot range of the field map figures

slice_to_plot = [ 7 28 35];    % ORIG 15, 28, 45
indices_1 = [35:110]'; indices_2 = [35:110]'; % Image crop to fit images in array

for ii=1:numel(slice_to_plot)
    % Field Maps acquired in experiment
    experimental_masked(:,:,slice_to_plot(ii))=field_map(:,:,slice_to_plot(ii)).*Mask_nan(:,:,slice_to_plot(ii));
    % Field Maps simulated using the loop's estimated position
    simulated_masked(:,:,slice_to_plot(ii))=b1_z_scaled(:,:,slice_to_plot(ii)).*Mask_nan(:,:,slice_to_plot(ii));
    % the result of the subtraction from both above shall ideally be 0
    difference(:,:,slice_to_plot(ii))=experimental_masked(:,:,slice_to_plot(ii))-simulated_masked(:,:,slice_to_plot(ii)).*Mask_nan(:,:,slice_to_plot(ii));
    
    %%mean_1(ii)= abs(mean(mean(mean(difference(:,:,slice_to_plot(ii))))))
    
    % chain plots/data together
    temp_row = [experimental_masked(indices_1,indices_2,slice_to_plot(ii)) simulated_masked(indices_1,indices_2,slice_to_plot(ii)) difference(indices_1,indices_2,slice_to_plot(ii))];
    temp_plot = [temp_plot;temp_row];
end

figure(10)
imagesc2(temp_plot,plot_range) %white background using imagesc2 function!! (External function)
axis image, axis off
colormap(jet), h=colorbar(),set(h,'FontSize',20)

%% SNR Images
disp('...generate SNR images from RAW data.')
SNR_noSHIM.image.flagDoAverage = true; SNR_noSHIM.image.flagRemoveOS = true;
SNR_SHIM.image.flagDoAverage = true; SNR_SHIM.image.flagRemoveOS = true;

% View rawData
%SNR_noSHIM.image.dataSize(1:11)
%                                   Col Cha Lin Par Sli Ave Phs Eco Rep Set Seg
temp_SNR_noSHIM = SNR_noSHIM.image( :,  :,  :,  1,  1,  1,  1,  1,  1,  1,  :);
figure(99), imagesc(abs(squeeze(temp_SNR_noSHIM(:,1,:, 1))).^0.2);

% Data reconstruction (fft)
img_SNR_noShim = zeros([SNR_noSHIM.image.NCol/os, SNR_noSHIM.image.NLin, SNR_noSHIM.image.NSli], 'single');
for sli = 1:SNR_noSHIM.image.NSli
    temp_SNR_noSHIM = SNR_noSHIM.image(:,:,:,1,sli); % read in the data for slice 'sli'
    fft_dims = [1 3]; % fft in col and lin dimension
    for f = fft_dims
        temp_SNR_noSHIM = ifftshift(ifft(fftshift(temp_SNR_noSHIM,f),[],f),f);
    end
    img_SNR_noShim(:,:,sli) = squeeze(sqrt(sum(abs(temp_SNR_noSHIM).^2,2))); % sum-of-square coil combination:
end

img_SNR_Shim = zeros([SNR_SHIM.image.NCol/os, SNR_SHIM.image.NLin, SNR_SHIM.image.NSli], 'single');
for sli = 1:SNR_SHIM.image.NSli
    temp_SNR_SHIM = SNR_SHIM.image(:,:,:,1,sli); % read in the data for slice 'sli'
    fft_dims = [1 3]; % fft in col and lin dimension
    for f = fft_dims
        temp_SNR_SHIM = ifftshift(ifft(fftshift(temp_SNR_SHIM,f),[],f),f);
    end
    img_SNR_Shim(:,:,sli) = squeeze(sqrt(sum(abs(temp_SNR_SHIM).^2,2))); % sum-of-square coil combination:
end

figure(300)
subplot(2,1,1), imagesc(rot90(img_SNR_noShim(:,:)), [0, 0.7*max(img_SNR_noShim(:))]);
colormap jet, h=colorbar(),set(h,'FontSize',14), axis image off;
title('SNR without SHIM');
subplot(2,1,2), imagesc(rot90(img_SNR_Shim(:,:)), [0, 0.7*max(img_SNR_Shim(:))])
h=colorbar(),set(h,'FontSize',14), colormap jet, axis image off;
title('SNR with SHIM');

%saveas(figure(300),[outdir,['SNRmaps_1.fig']]);
%saveas(figure(300),[outdir,['SNRmaps_1.png']]);
