clear; close all;

%  Compute Total Counts in a Digital Image 
function TIFF_Image_Quantification()

clear;
close all;
dirs.images = '/Users/brentfisher/Documents/MATLAB/UVQ/matDat/2024_04_17--19_LabImages/blue450nm_calibration_of_Chameleon(sn 19415615)'

% Select Files (MANY FILES)
%     cd(dirs.images);
%     [IMGfile_list,pathname,filterindex] =uigetfile({'*.tif','*.tiff'},'Multiselect','on');    
%     if isstr(IMGfile_list)
%         Nfiles = 1; tmp = IMGfile_list;  clear IMGfile_list;
%         IMGfile_list{1} = tmp;
%     else
%         Nfiles = length(IMGfile_list);
%     end

    cd(dirs.images);
    IMGfile_list = dir('*.tiff' )
    Nfiles = length(IMGfile_list);
    
        
    %Loop over all files 
    for kf=1:Nfiles
        %filename = [pathname,IMGfile_list{kf}];
        filename =  [IMGfile_list(kf).folder,'/',IMGfile_list(kf).name]; 
        I0 = imread(filename);
        
        %ID saturated pixels
        satPixList = find(I0>=255);
        
        
        %statistics
        TotalCounts(kf,1)       = sum(I0(:));
        NumSatPix(kf,1)         = length(satPixList);
        BGLevelEstimate(kf,1)   = mean(mean(I0(1:10,1:10))); %uupper corner)

        disp([IMGfile_list(kf).name,',',num2str(TotalCounts(kf,1)),',', num2str(NumSatPix(kf,1)),',', num2str(BGLevelEstimate(kf,1)) ]);
        
    end

figure; imagesc(Isum);
