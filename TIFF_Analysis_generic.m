% 1. Open TIFF containing a single peak

%
function TIFF_Analysis_generic()
    close all;
    clear;
    program_name = 'image analysis'
%% get data
##    dirs.f  = 'C:\Users\brent\MATLAB\UVQ';
##    dirs.f = 'G:\My Drive\Analyses(BF)\Matlab\UVQ';
##    start_path = 'C:\Users\brent\MATLAB\matDat\';
##    start_path = 'C:\Users\brent\MATLAB\matDat\2024_07_26_AlN1C_WG6.4_Throughput_JaiGoImages\';
##    start_path = 'C:\Users\brent\Pictures\'
##    start_path = 'G:\Shared drives\Corp Main\Engineering\LAB EXPERIMENTS\';
##    start_path = 'G:\Shared drives\Corp Main\Engineering\LAB EXPERIMENTS\2024_08_19_LabData\UV_Throughput_Pics';
##    start_path = 'G:\Shared drives\Corp Main\Engineering\LAB EXPERIMENTS\2024_09_10_LabData\20240911 Data';
    start_path = 'G:\Shared drives\Corp Main\Engineering\LAB EXPERIMENTS\2024_09_10_LabData\20240911 Data';

    %Select by Hand
    cd(start_path);
    [filenames,path] = uigetfile('*.tiff','Multiselect','on');
    Nf = length(filenames);
    if isstr(filenames)
        filename = filenames; clear filenames;
        filenames{1} = filename;
        Nf=1;
    end

%% Analyze Image
    I0 = imread([path,filenames{1}]);

    if size(I0,3)==4
        I0 = I0(:,:,1:3);
    end
    disp('Got the File')

    disp(['max(image) = ',num2str(max(I0(:)))]);

%     imtool(I0);

    figure; imagesc(I0); title(['max(image) = ',num2str(max(I0(:)))]);



%% Select ROI
    OK=0;
%     figure; imshow(I0)
    figure; imagesc(log10(double(I0)));

    while OK==0

        h = msgbox('Click two corners to define ROI',program_name);
        uiwait(h);
        [Ps_x(1),Ps_y(1)] = ginput(1);
        hold all; plot(Ps_x(1),Ps_y(1),'o','Color','r');
        [Ps_x(2),Ps_y(2)] = ginput(1);
        hold all; plot(Ps_x(2),Ps_y(2),'o','Color','r');
        hold all; plot([Ps_x(1),Ps_x(2),Ps_x(2),Ps_x(1),Ps_x(1)],...
                       [Ps_y(1),Ps_y(1),Ps_y(2),Ps_y(2),Ps_y(1)],'Color','r');

        ans = questdlg('Is the ROI OK?','verify','OK','No','OK')
        if strcmp(ans,'OK')    OK=1;  end
    end
    Iroi = I0(floor(min(Ps_y)):floor(max(Ps_y)),floor(min(Ps_x)):floor(max(Ps_x)),:);
    outfilename = [filenames{1}(1:end-4),'_ROI.tiff'];
    imwrite(Iroi,[path,outfilename]);

end
