
function Plot_LIV_Files()

clear;
close all;

dirs.h = '/Users/brentfisher/Documents/MATLAB';

% dirs.d = '/Volumes/GoogleDrive/Team Drives/Corporate Main/Lidar Lab Data (Share)/L-I-V Data/2018_01_30_Emitter03';
dirs.d = '/Volumes/GoogleDrive/Team Drives/Corporate Main/Lidar Lab Data (Share)/L-I-V Data/2018_01_30_VixarStandards';

%% get files
cd(dirs.d);
filelist = dir('*.txt');

% isel = find(filelist(:).name 
% Assumed Format:  [Voltage (V),	Current (A),	Optical Power (Watts)]

% [1] Examine all files & Determine common CURRENT axis
Imin=[]; Imax = []; dI=[];
for kf=1:length(filelist)
    filename =  filelist(kf).name;
    A = dlmread([dirs.d,'/',filename],'\t',1,0);
    Imin = min([Imin;A(:,2)]);
    Imax = max([Imax;A(:,2)]);
    dI = min([dI;median(diff(A(:,2)))]);    
end

Iaxis = [Imin:dI:Imax]';

% [2] Read All Data and Interpolate to *same* CURRENT axis
for kf=1:length(filelist)
    filename =  filelist(kf).name;
    A = dlmread([dirs.d,'/',filename],'\t',1,0);
    Vmat(:,kf) = interp1(A(:,2),A(:,1),Iaxis);  %Voltage v Current
    Lmat(:,kf) = interp1(A(:,2),A(:,3),Iaxis);  %Laser Power v Current
    legnames{kf} = replace(filename,'_',' ');
    
    %plot (single file)
    figure; 
    [hAx,hLine1,hLine2] = plotyy(A(:,2)*1000,A(:,1),A(:,2)*1000,A(:,3)*1000);
    title(replace(filename,'_',' ')); xlabel('Current [mA]'); grid on;
    ylabel(hAx(1),'Voltage [V]') % left y-axis 
    ylabel(hAx(2),'Power [mW]') % right y-axis
end

% [3] Plot Data
figure; plotyy(Iaxis,Vmat,Iaxis,Lmat); legend(legnames)
hold on; plot(Iaxis,Lmat);

% [4] Save to CSV
csvwrite([dirs.d,'/temp.csv'],[Iaxis,Vmat,Lmat]);



% % determine common VOLTAGE axis
% Vmin=[]; Vmax = []; dV=[];
% for kf=1:length(filelist)
%     filename =  filelist(kf).name;
%     A = dlmread([dirs.d,'/',filename],'\t',1,0);
%     Vmin = min([Vmin;A(:,1)]);
%     Vmax = max([Vmax;A(:,1)]);
%     dV = min([dV;median(diff(A(:,1)))]);    
% end














