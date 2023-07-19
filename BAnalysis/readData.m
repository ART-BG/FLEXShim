% -------------------------------------------------------------------------
% Title: readData.m
% Author: Bernhard Gruber
%
% Purpose: Read in MRI data
%
% Versionhistory:
%   03/2021 - Read in MRI DICOMs and display them using montage()
%   06/2021 - 
%
%
% -------------------------------------------------------------------------
clear all;close all;
curr_path = pwd;
outdir=[curr_path,[filesep 'Output' filesep]];

Archive = ('C:\Users\BGruber\Documents\!Work\02.2 - MUW\02_Projects\Measurements\');
Datafolder = ([Archive,'3TPrisma_20042021\']);
Subject = ([Datafolder,'s4_dicoms']);

show_plots = 1;
% save_data=1;
% save_fig=1;

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%curr_path = [uigetdir(),filesep];
%p = genpath(curr_path); addpath(p); % add all subfolders within path
%dicom_temp = dir([curr_path, '*.IMA']);

%% Read B0 mapping data
%prompt="Enter delta TE [e.g. 3T = 0.00246s, 7T = 0.00102s] = ";
%delta_TE = input(prompt);
%delta_TE = 0.00246; % for 3T - for 7T: 1.02 ms
%delta_TE = 0;

disp('Read in dicom data...');

%fm_path = [uigetdir(),filesep];
fm_path = ([Subject, 'B0_exhaled\P4_P4\FLEXSHIM_TEST_SEQS_20210420_115632_636000\ASPIRE_P_KE_GRE_ASPIRETEST_B0MAP_3D_FULLFOV_EXHALED_0011\']);
[dicom_temp,field_map_matrix,field_map_res,TE] = dicom_fm_import2(strcat(fm_path));
unshimmed_phase = dicom_temp;%/(2*pi*(TE(2)-TE(1)));

figure(1)
if show_plots == 1;
    montage(unshimmed_phase, [], 'DisplayRange', [0 100]);
end
%%
fm_path = ('C:\Users\BGruber\Documents\!Work\02.2 - MUW\02_Projects\Measurements\3TPrisma_20042021\s4_dicoms\B0_exhaled\P4_P4\FLEXSHIM_TEST_SEQS_20210420_115632_636000\KE_GRE_ASPIRETEST_B0MAP_3D_FULLFOV_EXHALED_0012\');
[dicom_temp,field_map_matrix,field_map_res,TE] = dicom_fm_import2(strcat(fm_path),1,0);
magnitude = dicom_temp;

figure(2)
if show_plots == 1;
    montage(magnitude, [], 'DisplayRange', [0 100]);
end

field_map_FoV = field_map_res.*double(field_map_matrix); % field map Field of View

%delta_TE = (TE(2)-TE(1));

disp('=============================================================');
%disp(['delta TE                   : ', num2str(delta_TE)]);
disp(['Field map matrix (ex)      : ', num2str(field_map_matrix)]);
disp(['Field map FoV (ex)         : ', num2str(field_map_FoV)]);
disp(['Field map Resolution (ex)  : ', num2str(field_map_res(1:2))]);
disp(['Field map Slice Thick. (ex): ', num2str(field_map_res(3))]);
disp('=============================================================');

%delta_TE = (TE(2)-TE(1));

% %% Magnitude Images
% [dicom_temp] = dicom_fm_import2(strcat(fm_path),delta_TE,0);
% magimg = dicom_temp;
% figure(2)
% if show_plots == 1;
%     montage(magimg, [], 'DisplayRange', [0 100]);
% end