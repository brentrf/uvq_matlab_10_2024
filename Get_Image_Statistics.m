% Open Image, Select ROI,  Deliver MONOCHROME Statistics
% 18 July 2024

%%
clear;
close all
program_name = 'Get Image Statistics';


%select file
dirs.d = 'C:\Users\brent\MATLAB\matDat';
cd(dirs.d); [FileName,PathName] = uigetfile({'*.tiff'},program_name);
filenametitle = replace(FileName,'_',' '); 
 
%open file 
if strcmp(FileName(end-3:end),'tiff')  %for TIFF file input
    I0 = imread([PathName,FileName]);
    [Nr,Nc,Nlayers] = size(I0);    
    if size(I0,3)>=3
        I1 = rgb2gray(I0); %convert to monochrom
    else
        I1 = I0;        %take as monochrome
    end
elseif strcmp(FileName(end-3:end),'.mat')   %for HDR (*.mat) file input
    INDAT = load([PathName,FileName]);
    I0 = uint32( INDAT.I_HDR );
    [Nr,Nc,Nlayers] = size(I0);   
    I1 = I0(:,:,1);            
end

%show monochrome image and define points for ROI
hf1 = figure; 
himg  = imagesc(log10(double(I1))); colorbar;
h = msgbox('Click Corners of ROI',program_name);
uiwait(h);
[Ps_x,Ps_y] = ginput(1);
figure(hf1); hold all;  plot(Ps_x,Ps_y,'o','Color','y');
[Ps_x(2),Ps_y(2)] = ginput(1);
figure(hf1); hold all;  plot(Ps_x(2),Ps_y(2),'o','Color','y');
figure(hf1); hold all;  plot([Ps_x(1),Ps_x(2),Ps_x(2),Ps_x(1),Ps_x(1)],... 
                             [Ps_y(1),Ps_y(1),Ps_y(2),Ps_y(2),Ps_y(1)],'Color','y');

%Get ROI and compute Statistics
Ps_x = sort(round(Ps_x,0));
Ps_y = sort(round(Ps_y,0));
I_roi  = I1( [Ps_y(1):Ps_y(2)], [Ps_x(1):Ps_x(2)] );

I_roi_mean = mean(I_roi(:));
I_roi_stdev = std(double(I_roi(:)));
disp(FileName);
disp(I_roi_mean);
disp(I_roi_stdev);




