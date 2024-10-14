
%  Compute Total Counts in a Digital Image 
function Demosaic_Images()

clear;
close all;
% dirs.images = '/Users/brentfisher/Documents/MATLAB/UVQ/matDat/2024_04_17--19_LabImages/InputImages_with_COLOR_camera/'
% dirs.images = '/Users/brentfisher/Documents/MATLAB/UVQ/matDat/2024_04_17--19_LabImages/reverse_illumination_using_50umFiber(E)/';
dirs.images = '/Users/brentfisher/Documents/MATLAB/UVQ/matDat/2024_04_17--19_LabImages/reverse_Illumination_using_LensedFiber/'

    cd(dirs.images);
    IMGfile_list_tiff = dir('*.tiff' );
    IMGfile_list_bmp = dir('*.bmp' );
    IMGfile_list_png = dir('*.png' );
    IMGfile_list_jpeg = dir('*.jpeg' );
    IMGfile_list = cat(1,IMGfile_list_tiff,IMGfile_list_bmp,IMGfile_list_png, IMGfile_list_jpeg);
    Nfiles = length(IMGfile_list);
    
        
    %Loop over all files 
    for kf=1:Nfiles
        %filename = [pathname,IMGfile_list{kf}];
        filename =  [IMGfile_list(kf).folder,'/',IMGfile_list(kf).name]; 
        I0 = imread(filename);
        
        BayerPattern = 'rggb'; %'gbrg' | 'grbg' | 'bggr' | 'rggb'
        Irgb = demosaic(I0,BayerPattern);
        
        imshow(Irgb)
        pause(0.5);
        
        %save
        if ~isdir([dirs.images,'/demosaiced_rgb'])
            mkdir([dirs.images,'/demosaiced_rgb'])
        end
        s = split(IMGfile_list(kf).name,'.');
        output_filename = [dirs.images,'/demosaiced_rgb/',s{1},'_rgb.tiff'];
        imwrite(Irgb,output_filename)
        
        %progrss statment
        disp([IMGfile_list(kf).name,'   ------>  ',s{1},'_rgb.tiff' ]);
            
    end
end