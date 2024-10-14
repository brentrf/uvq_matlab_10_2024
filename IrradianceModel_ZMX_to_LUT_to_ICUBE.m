% IRRADIANCE ZEMAX DATA --> LOOK UP TABLE --> DATA CUBE
%
%  *[1]* MAIN():   
%        - reads in LUT.mat OR  Generates LUT  -->  ICUBE.mat
%
%  *[2]* readAll_ZemaxDataTXT_Files(zdat,dirs):   
%        - loop that fetches ZMX data from TXT Files --> A{} 
%
%  *[3]* Generate_Irradiance_LUT(zdat,dirs):      
%        - ZMX detector data stack A{}   -->  LUT.mat
%  
%  **** copy_of_read_Zemax_TXT(filename)
%        - this is just a valuable function that resides also in another file
%  **** Scale_Irradiance_LUT(LUT)
%        - Needs work.  Meant to scale LUT....

function MAIN()
    close all;
    clear;
    dirs.f  = '/Users/brentfisher/Documents/MATLAB/UVQ';
    start_path = '/Users/brentfisher/Documents/MATLAB/UVQ/irradiance_LUTS/';
    zdat.path = '/Users/brentfisher/Documents/MATLAB/UVQ/irradiance_LUTS/uniformCone_30deg/';
    zdat.path = [uigetdir(start_path),'/'];
%% Inputs
    %datasource = 'matfile'; 
    datasource = 'ZemaxTXT';
    
    
    
    %% GET LUT
    %get Data: from mat file or Build LUT
    switch datasource
        case 'matfile'
            cd(zdat.path)
            load('LUT_cm.mat');
        case 'ZemaxTXT'
            LUT = Generate_Irradiance_LUT(zdat,dirs);
            savename = 'LUT_cm.mat';
            save([zdat.path,savename],'LUT');
                
        otherwise
            error('No data source found');
    end
    
        
    %%  BUILD DATA CUBE:   Generate Data Cube Based on LUT
    % LUT --->  Data Cube
    
    %Define CUBE Dimensions -------------------------
    ICUBE.dxyz = 1.0;      %[cm]
    ICUBE.Lx = max(LUT.XX_cm(:)) - min(LUT.XX_cm(:));     %[cm]
    ICUBE.Ly = max(LUT.YY_cm(:)) - min(LUT.YY_cm(:));     %[cm]
    ICUBE.Lz = max(LUT.Z_cm(:));     %[cm]
    ActualSourceSize.L_cm = 0;     %actual Length of source (tube) in question
    ActualSourceSize.W_cm = 0;     %actual WIDTH of source (tube)in question
    ActualSourceSize.Ptot_W = LUT.Ptot;    %Ptot    
    prompt = {'enter Voxel dimension in [cm]:','Enter CUBE WIDTH(X) in [cm]:  ','Enter CUBE LENGTH(Y) in [cm]:','Enter CUBE HEIGHT(ZY) in [cm]:',...
            'source size LENGTH in [cm]:','source size WIDTH in [cm]:','source POWER in [W]:'};
    dlg_title = 'Source size';
    num_lines = 1;
    defaultans = {num2str(ICUBE.dxyz), num2str(ICUBE.Lx),num2str(ICUBE.Ly),num2str(ICUBE.Lz),num2str(ActualSourceSize.L_cm),num2str(ActualSourceSize.W_cm),num2str(ActualSourceSize.Ptot_W)};
    answer = inputdlg(prompt,dlg_title,num_lines,defaultans);
    ICUBE.dxyz = str2num(answer{1});
    ICUBE.Lx = str2num(answer{2});
    ICUBE.Ly = str2num(answer{3});
    ICUBE.Lz = str2num(answer{4});
    ActualSourceSize.L_cm = str2num(answer{5});
    ActualSourceSize.W_cm = str2num(answer{6});
    ActualSourceSize.Ptot_W = str2num(answer{7});
         
    %BUILD DATA CUBE
    Nz = size(LUT.Irradiance,3);
    zv = [1:ICUBE.dxyz:ICUBE.Lz]';
    xv = [-ICUBE.Lx/2:ICUBE.dxyz:ICUBE.Lx/2];
    yv = [-ICUBE.Ly/2:ICUBE.dxyz:ICUBE.Ly/2]';   
    [ICUBE.XXX,ICUBE.YYY,ICUBE.ZZZ] = meshgrid(xv,yv,zv);
    xxx0 = repmat(LUT.XX_cm,1,1,Nz);
    yyy0 = repmat(LUT.YY_cm,1,1,Nz);
    zzz0 = zeros(size(yyy0));
    for kz=1:Nz zzz0(:,:,kz) = ones*LUT.Z_cm(kz); end
    %execute interpolation to Cube in *** mW/cm2 ***
    ICUBE.III_mWcm2 = interp3(xxx0, yyy0, zzz0, LUT.Irradiance*(ActualSourceSize.Ptot_W/LUT.Ptot)*1000,  ICUBE.XXX,ICUBE.YYY,ICUBE.ZZZ);

    
    %%SAVE DATA CUBE
    ICUBE.Ptot_W = ActualSourceSize.Ptot_W;
    savename = 'ICUBE.mat';
    ICUBE.Pinput_allLamps_Wopt  = 0.3;
    save([zdat.path,savename],'ICUBE','-v7.3');
    
    %% ANALYZE CUBE:  PLOT & Average Value of Cube
%     %3D Plot --- numbers need to be adjusted to give a good result
%     figure; 
%     h = slice(ICUBE.XXX,ICUBE.YYY,ICUBE.ZZZ, ICUBE.III_mWcm2,[0],[0],[-1]);
%     alpha('color')
%     set(h,'EdgeColor','none','FaceColor','interp','FaceAlpha','interp')
%     alphamap('increase',.2)
%      alphamap('rampdown'); 
%     colorbar();
%     xlim([-50,50]);    
%     ylim([-50,50]);    
%     zlim([0,50]);
%         
%     
%     %TEST VOLUME:  volumetric average irradiance
%     Zmin_cm = 0.0;
%     Zmax_cm = 304.5;
%     LplotXY = 152.4;  %full width 
%     LplotXY = 304.5;  %full width 
%     
%     %     ikeep_in_average = find(~isnan(ICUBE.III_mWcm2(:)));
%     ikeep_in_average = find(~isnan(ICUBE.III_mWcm2(:)) & abs(ICUBE.XXX(:))<LplotXY/2 & abs(ICUBE.YYY(:))<LplotXY/2  &  ICUBE.ZZZ(:)< Zmax_cm  &  ICUBE.ZZZ(:) > Zmin_cm);
%     Irr_vol_average_mWcm2 = mean( ICUBE.III_mWcm2(ikeep_in_average) );  %[mW/cm2]
%     
%     %Report Result in Title
%     title({['volume avg Irr: ',num2str(Irr_vol_average_mWcm2),'mW/cm2'];...
%         ['Optical Power = ',num2str(ActualSourceSize.Ptot_W*1000),'mW ...ZemaxFiles=','','deg'];...
%         ['Test Volume(Lx, Ly, Hz): (',num2str(LplotXY),', ',num2str(LplotXY),', ',num2str(Zmax_cm-Zmin_cm),') cm'  ];...
%         ['dX, dZ = (',num2str(ICUBE.dxyz),', ',num2str(ICUBE.dxyz),') cm'  ]});
%     
%    
%   %2D Plot ---     
%     Zsel_cm = 304.5;
%     [dummy,isel] = min( abs(ICUBE.ZZZ(1,1,:) -  Zsel_cm) );
%     XX = ICUBE.XXX(:,:,isel); YY=ICUBE.YYY(:,:,isel); II = ICUBE.III_mWcm2(:,:,isel);
%     iselctr7cm = find( sqrt(XX(:).^2+YY(:).^2)<3.5);
%     IIavg_ctr_mWcm2 = mean(II(iselctr7cm));
%     figure; imagesc(ICUBE.XXX(1,:,isel),ICUBE.YYY(:,1,isel),ICUBE.III_mWcm2(:,:,isel)); colorbar;  axis equal
%     title({'Irradiance [mW/cm2] in this plane';' ** !not accounting for cosine loss @ plane! **';...
%         ['Optical Power = ',num2str(ActualSourceSize.Ptot_W*1000),'mW ...FWHM=','','deg... plane: Z=',num2str(Zsel_cm),'cm'];...
%         ['Avg Irradiance @ Ctr: ',num2str(IIavg_ctr_mWcm2*1e3),'uW/cm2']});    
%     Power_plotted = sum(II(:))*(ICUBE.dxyz)^2    
%     
%     
end






%% Read ALL Zemax TXT files from a directory

function [Adat,meta] = readAll_ZemaxDataTXT_Files(zdat,dirs)
    fileslist = dir([zdat.path,'/*.TXT']);
        Nf = length(fileslist);
        for kf=1:Nf %loop over files
            filename = [zdat.path,fileslist(kf).name];
            [Adat{kf},meta0] =  read_Zemax_TXT(filename);  %read this file.
            meta.Z_cm(kf,1) = meta0.Z_cm;
            meta.Lx_cm(kf,1) = meta0.Lx_cm;
            meta.Ly_cm(kf,1) = meta0.Ly_cm;
            meta.Ptot_W(kf,1) = meta0.Ptot_W;            
        end %done with this file

end


%% ================

function out = Generate_Irradiance_LUT(zdat,dirs)

    %Read Detector Data from Zemax
    csv_filelist    = dir([zdat.path,'/*.csv']);
    ZMXTXT_filelist = dir([zdat.path,'/*.TXT']);
    
    if length(ZMXTXT_filelist)>0  %AS TEXT FILES (automated output)
        [A,meta] = readAll_ZemaxDataTXT_Files(zdat,dirs);     %DEFAULT !!!
        Z_cm = meta.Z_cm;
        Lx_cm = meta.Lx_cm;
        Ly_cm = meta.Ly_cm;
        LUT.Ptot = max(meta.Ptot_W);  %record total power
        
    elseif length(csv_filelist)>0 %AS CSV FILES (special case: manually built with Z_cm in filename)
        Nf = length(fileslist)
        for kf=1:Nf
            filename = fileslist(kf).name;
            tmp = split(filename,{'_','cm.'});
            Z_cm(kf,1) = str2num(tmp{2});
            Lx_cm = Z_cm*6;                             %check this!!
            Ly_cm = Lx_cm;
            A{kf} = csvread([zdat.path,'/',filename],24,1);
            ans = inputdlg('Enter Ptot from Zemax raytrace in [W]:','Ptot',1,{'1'});
            LUT.Ptot = str2num(ans{1})
        end
        
    else %neither TXT or CSV found in this path
       warning(['Neither TXT nor CSV files found at this path:  ',zdat.path])
       return
    end
    
    disp(' '); disp('****************');
    disp(['Data Collected from:  ',zdat.path]);
    disp('****************'); disp(' ');
    
    
    
    %record SOURCE SIZE of LUT data
    LUT.scale_sourceLength_cm = 0;
    LUT.scale_sourceWidth_cm = 0;
    prompt = {'Enter source WIDTH(X) in [cm] (point source = 0):', 'Enter source LENGTH(Y) in [cm] (point source = 0):  '};
    dlg_title = 'Enter Source size (used in ZMX)';
    num_lines = 1;
    defaultans = {num2str(LUT.scale_sourceWidth_cm),num2str(LUT.scale_sourceLength_cm)};
    answer = inputdlg(prompt,dlg_title,num_lines,defaultans);
    LUT.scale_sourceWidth_cm = str2num(answer{1});
    LUT.scale_sourceLength_cm = str2num(answer{2});
 
    
    
    
    % **** Process into Data Cube *****
    
    %sort Z 
    Nf = length(Z_cm);
    [Z_cm_srt, Indc_srt] = sort(Z_cm);
    for kf=1:Nf
        B{kf} = A{Indc_srt(kf)};
    end 
    Z_cm = Z_cm_srt;
    A = B; %  clear B  Z_cm_srt;


    
    
    %Interpolate Inputs to common XY Samples
    dxy = 0.5; %[cm]
    X = [-max(Lx_cm)/2: dxy :max(Lx_cm)/2];
    Y = [-max(Ly_cm)/2: dxy :max(Ly_cm)/2]';
    [LUT.XX_cm,LUT.YY_cm] = meshgrid(X,Y);
    LUT.Z_cm = Z_cm;
    
    [nr,nc]  = size(LUT.XX_cm);
    LUT.Irradiance = zeros(nr,nc,Nf);
    hf1 = figure;
    for kf=1:Nf
        [nr,nc] = size(A{kf});
        xA = linspace(-Lx_cm(kf),Lx_cm(kf),nc);
        yA = linspace(-Ly_cm(kf),Ly_cm(kf),nr)'; 
%         C = interp2(xA,yA,A{kf},X,Y);
        LUT.Irradiance(:,:,kf)= interp2(xA,yA,A{kf},X,Y);
        figure(hf1); imagesc(LUT.Irradiance(:,:,kf));
        title(['Z=',num2str(Z_cm(kf)),'cm']); pause(0.4); 
    end
    
    out=LUT;
    
end




%%  *** NOT SURE THIS WORKS AS INTENDED ****
function  Ivals = Scale_Irradiance_LUT(x,y,z,ActualTubeSize,LUT)
    warning(' THIS FUNCTION NEEDS REVIEW');
    warning(' THIS FUNCTION NEEDS REVIEW');
    warning(' THIS FUNCTION NEEDS REVIEW');
    warning(' THIS FUNCTION NEEDS REVIEW');
    %scale factor := (size of tube in question) / (size of Tube used for LUT)
    
    %determine scale factor
    scalefactorL = ActualTubeSize.L_cm / LUT.scale_tubeLength_cm;
    scalefactorR = ActualTubeSize.R_cm / LUT.scale_tubeRadius_cm;
    if scalefactorL >= scalefactorR
        scalefactor = scalefactorL;
        warning(['Aspect Ratio of desired Tube not matched to LUT.  Changed: (LxR = ',num2str(ActualTubeSize.L_cm),'x',num2str(ActualTubeSize.R_cm)...
                , ')   to (LxR = ',num2str(ActualTubeSize.L_cm),'x',num2str(ActualTubeSize.R_cm*scalefactor/scalefactorR),')']);
    elseif scalefactorL < scalefactorR
        scalefactor = scalefactorR;
        warning(['Aspect Ratio of desired Tube not matched to LUT.  Changed: (LxR = ',num2str(ActualTubeSize.L_cm),'x',num2str(ActualTubeSize.R_cm)...
                , ')   to (LxR = ',num2str(ActualTubeSize.L_cm*scalefactor/scalefactorL),'x',num2str(ActualTubeSize.R_cm),')']);
    end
    
    %interp look up
    x_scaled2lut = x / scalefactor;
    y_scaled2lut = y / scalefactor;
    z_scaled2lut = z / scalefactor;
    [xxx,yyy,zzz] = meshgrid(LUT.XX_cm(1,:),LUT.YY_cm(:,1),LUT.Z_cm);
    
    %execute interpolation
    Ivals = interp3(xxx,yyy,zzz, LUT.Irradiance,  x_scaled2lut, y_scaled2lut, z_scaled2lut  );
    
end


%% ================ 
% This function resides in another file, but a copy is kept here......
%  this function reads 2D irradiance data from Zemax Text File
% input:  string:  "path/filename"
% output1:  DAT = 2D matrix of numbers
% output2:  meta = struct containing meta data taken from header

function [DAT,meta] = copy_of_read_Zemax_TXT(filename)
            
    fID = fopen(filename);
    if fID~=-1
            %--------------------------------------------------------------

            %** Z *** ---get Z value for this file
            Nskip = 33;
            frewind(fID);
            for kk=1:Nskip
                tline1 = fgetl(fID);
            end
              %disp(tline1);
            tline2 = tline1(2:2:end); %remove alternating spaces from strong
            sline = split(tline2,':');
            meta.Z_cm = str2num(sline{2});

            %**L_x,L_y*** --- get L_x and L_y value for this file
            Nskip = 17;
            frewind(fID);
            for kk=1:Nskip
                tline1 = fgetl(fID);
            end
                    %disp(tline1)
            tline2 = tline1(2:2:end); %remove alternating spaces from strong
            sline = split(tline2,' ');
            meta.Lx_cm  = str2num(sline{2}); 
            meta.Ly_cm  = str2num(sline{5}); 
            
            
            %**TOTAL POWER*** --- get Total Power for this file
            Nskip = 23;
            frewind(fID);
            for kk=1:Nskip
                tline1 = fgetl(fID);
            end
                        %disp(tline1)
            tline2 = tline1(2:2:end); %remove alternating spaces from string
            sline = split(tline2,{':',' '});
            meta.Ptot_W  = str2num(sline{9});   
                       
            %--------------------------------------------------------------
            
            
            %*** READ 2D Data ****  
                % % %             Nskip = 49;  %skip to this line
                % % %             frewind(fID);
                % % %             for kk=1:Nskip
                % % %                 tline1 = fgetl(fID);
                % % %             end
                % % %             disp(tline1);
                % %             
            frewind(fID);
            tline1=fgetl(fID);
            tline2=tline1(2:2:end);
            n=0;
            while n<1000 & ~(length(strfind(tline2,'Units'))~=0 & length(strfind(tline2,'Watts')) ~= 0)  
             	tline1=fgetl(fID);
                tline2=tline1(2:2:end);
                n=n+1;
            end
            Nskip=6;    %SKIP 6 LINES
            for kk=1:Nskip	 tline1=fgetl(fID); end
            
            %disp(tline1);
            
            DAT = [];  
            while feof(fID)==0  %while NOT end of file
                %remove spaces from this line
                tline2 = tline1(2:2:end); %remove alternating spaces from string
                %convert text to numeric & store
                if length(tline2)>0
                    dat_thisline = textscan(tline2,'%f');
                    if length(DAT)==0
                        DAT = dat_thisline{1}';
                    else
                        DAT = [DAT;dat_thisline{1}'];
                    end
                end
                %read next line 
                tline1 = fgetl(fID);   
            end
            DAT = DAT(:,2:end);  %trim first column off  (not data, just indices)

            %--------------------------------------------------------------
            
            fclose(fID); 
    else
        warning('File not found')
    end%(if) 
                    
end





