% 1. Open TIFF containing a single peak 
 
%
function TIFF_Image_Stack_Analysis()
    close all;
    clear;
    program_name = 'image analysis';
    
    start_path = 'G:\Shared drives\Corp Main\Engineering\LAB EXPERIMENTS\2024_09_23_LabData\2024_09_24_UvqA-sd7_Lwgwh520id19.1__SHG_Search_UVImages';
    
%% select IMAGE STACK
   
    %Select by Hand
    cd(start_path);
    [filenames,path] = uigetfile('*.tiff','Multiselect','on');
    Nf = length(filenames);
    if isstr(filenames)
        filename = filenames; clear filenames;
        filenames{1} = filename;
        Nf=1;
    end
    
%% Read All Files
    for kf = 1:Nf        
       [DAT.I0{kf},DAT.peaks{kf}] = get_image([path,filenames{kf}]) ;                  
    end
    
%% select BG IMAGE
    %Select by Hand
    cd(start_path);
    [BG_filename,BG_path] = uigetfile('*.tiff','Multiselect','off');
    if BG_filename==0   BGflag = 0;
    else                BFflag = 1;
    end
    
    if isstr(BG_filename)   BG_filename = BG_filename;  
    end
    
    
%% Read BG Image
    [BG.I0,BG.peak] = get_image([BG_path,BG_filename]) ;   
 
%% Define ROI
ROI.r1 = BG.peak.r-50; 
ROI.r2 = BG.peak.r+50; 
ROI.c1 = BG.peak.c-50; 
ROI.c2 = BG.peak.c+50; 

%% Generate Animation
hf = figure;
for kf = 1:Nf
    figure(hf);
    I1 = DAT.I0{kf}(ROI.r1:ROI.r2,ROI.c1:ROI.c2) - BG.I0(ROI.r1:ROI.r2,ROI.c1:ROI.c2);
    imagesc(I1); title(num2str(kf));
    drawnow;
    F(kf) = getframe;
end
    
figure;    
movie(F,2,1)


end  

%% Load TIFF Image + Find Center Peak

function [I0,peak] = get_image(full_filename)
    I0  = imread(full_filename);
    if size(I0,3)==4
        I0  = I0 (:,:,1:3);
    end
    I0 = double(I0);

    %convolve and find peak of iamge    
    A = conv2(I0,ones(10))/(100);
    [peak.maxval,imax] = max(A(:));
    [peak.r,peak.c] = ind2sub(size(A),imax);
    
    
    s = split(full_filename,'\');
    filename = s(end);
    disp(['Got Image: ',filename]);
end

