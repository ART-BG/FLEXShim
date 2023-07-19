% function out = dicom_fm_import(pathname,delta_TE,scale_switch)
% Select all Siemens dicoms in selected folder (extension *.ima)
% Cris LaPierre 01/20/2014
% modified by Jason Stockmann in 2014
% modified by Bernhard Gruber in 2021

% set 'scale_switch = 1' for scaling field maps into the range [-pi, pi]
% set zero otherwise and use '0' for delta_TE

function [dicom_temp,field_map_matrix,field_map_res,TE] = dicom_fm_import2(pathname,varargin)
nvar = (length(varargin));

if nvar == 1
    delta_TE = varargin{1};
elseif nvar == 2
    delta_TE = varargin{1};     % delta_TE as input
    scale_switch = varargin{2}; % scaling of field maps into the range [-pi, pi]
else
    delta_TE = 1;
    scale_switch = 0; % no scaling of field maps into the range [-pi, pi]
end

% file option: standard [0]*.ima , [1]*.dcm, [2]*.nii
file_opt = 0;
files = dir(fullfile(pathname,'*.IMA'));

if isempty(files)
    files = dir(fullfile(pathname,'*.dcm'));
    file_opt = 1
end
if isempty(files)
    files = dir(fullfile(pathname,'*.nii'));
    file_opt = 2
end
% if no file type is detected
 if isempty(files)
%     files = dir(fullfile(pathname));
%     file_opt = 0
    file_opt = 3;
 end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Loop through all file names to load every sli
if file_opt == 0   % .IMA case
    %disp(['Number of images = ',num2str(length(files))]);
    for file_num = 1:length(files), 
        % Uses a regular expression to extract the slice number from the file
        % name (slice is the 5th spot if periods are used as delimiters.
        splitStr = regexp(files(file_num).name,'\.','split');

        % Use dicomread to load the data. Eval lets me create a dynamic name
        % for each variable that includes the current slice number
        eval(['temp1(:,:,' splitStr{5} ')=dicomread(fullfile(pathname,files(file_num).name));']);
    %     % Simple plot command to view each slice
    %     eval(['imagesc(abs(' 'slice' splitStr{5} '))'])
%         eval(['temp2(:,:,' splitStr{5} ')=dicominfo(fullfile(pathname,files(file_num).name));']);
    end
    % read echo times
    info_tmp=dicominfo(fullfile(pathname,files(file_num-1).name));
    TE(1)=info_tmp.EchoTime*1e-3; %read first echo time
    eval(['temp3(:,:,' splitStr{5} ')=dicominfo(fullfile(pathname,files(file_num).name));']);
    info_tmp2=dicominfo(fullfile(pathname,files(file_num).name));
    TE(2)=info_tmp2.EchoTime*1e-3; %read second echo time
    
    % read field map matrix
    eval(['temp3(:,:,' splitStr{5} ')=dicominfo(fullfile(pathname,files(file_num).name));']);
    field_map_matrix(1)= info_tmp.Rows;
    field_map_matrix(2)= info_tmp.Columns;
    field_map_matrix(3) = file_num;
    
    % read field map resolution
    field_map_res(1:2)=info_tmp.PixelSpacing/1000;
    field_map_res(3)=info_tmp.SliceThickness/1000;

elseif file_opt == 1   % .dcm case
    
    for file_num = 1:length(files),
    % Uses a regular expression to extract the slice number from the file
    % name (slice is the 5th spot if periods are used as delimiters.
    splitStr_1 = regexp(files(file_num).name,'\-','split');
    tempStr = splitStr_1{3};
    splitStr_2 = regexp(tempStr,'\.','split');
    
    % Use dicomread to load the data. Eval lets me create a dynamic name
    % for each variable that includes the current slice number
    eval(['temp1(:,:,' splitStr_2{1} ')=dicomread(fullfile(pathname,files(file_num).name));']);
%     % Simple plot command to view each slice
%     eval(['imagesc(abs(' 'slice' splitStr{5} '))'])
    end

elseif file_opt == 2   % .nii case
      
    hdr = load_nifti([files(1).name]);
    temp1 = hdr.vol;
    
    FOV_actual_fm(1) = hdr.pixdim(4)*hdr.dim(4);
    FOV_actual_fm(2) = hdr.pixdim(2)*hdr.dim(2);
    FOV_actual_fm(3) = hdr.pixdim(3)*hdr.dim(3);

%         for aa = 1:hdr.dim(2)
%             for bb = 1:hdr.dim(3)
%                 for cc = 1:hdr.dim(4)
%                     pos(aa,bb,cc) = hdr.vox2ras*[aa-1; bb-1; cc-1; 1];
%                 end
%             end
%         end
    temp1 = permute(temp1,[1 3 2]);
    %temp1 = permute(temp1,[1 2 3]); % BG
    scale_switch = 0;
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% set 'scale_switch = 1' for scaling field maps into the range [-pi, pi]
if file_opt ~= 3
    if scale_switch == 1
        dicom_temp = pi*((double(temp1)-2048)/2048)/delta_TE/2/pi;
    else
        dicom_temp = double(temp1);
    end
else
   disp('Please check to make sure DICOM files are in appropriate folder and that the folder is in the Matlab path.') 
   dicom_temp = temp1;
end


