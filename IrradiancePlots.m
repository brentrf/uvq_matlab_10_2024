% NORMAL INCIDENCE IRRADIANCE IN A PLANE BENEATH EXTENDED SOURCE
%
%   This is the *Normal Incidence* Irradiation at each point in the plane
%
%  These values to NOT correct for cosine loss over the flat area of a
%  surface. 
%  Instead this represents the intensity that is passing through a 
%  normally oriented "voxel" in this plane
%
%

%%  TBD
%   [1.3] - Analytic Tile
%    1.6 - check
%   [1.9] - Tube - Extended Source, Combined Near & Far
%   Add Tilt  
%   Add Superposition (i.e. many lamps)



%  2D:
%    -eACH = f(II)
%    - TLV = f(II)
%    - Prob_Sick = f(II) wells riley
%  make 3D
%
%
%%

function IrradiancePlots()
    clear;
     close all;
        %first UVC.... 
        IN.Power_tot        = 200;    %[mW]
        IN.L_tile           = 11;     %[cm]
        IN.W_tile           = 6;     %[cm]
        IN.lambertian_FWHM  = 60;    %full angle (1/2-max), for LAMBERTIAN (TILE)
        IN.alpha_tube       = 60;    %full angle (uniform), for UNIFORM    (TUBE)
        
        %Small Tile Lambertian ... 
        IN.Power_tot        = 12;    %[mW]
        IN.L_tile           = 1;     %[cm]
        IN.W_tile           = 2;     %[cm]
        IN.lambertian_FWHM  = 60;    %full angle (1/2-max), for LAMBERTIAN (TILE)
        
       %USHIO Lambertian ... 
        IN.Power_tot        = 150;    %[mW]
        IN.L_tile           = 6.0;     %[cm]
        IN.W_tile           = 4.5;     %[cm]
        IN.lambertian_FWHM  = 60;    %full angle (1/2-max), for LAMBERTIAN (TILE)
        
    IN.L_tube           = 100 ;   %[cm] ... only used for Tube_AnalyticInt & Ideal Near Field
    IN.Z                = 198.69;   %[cm]
    IN.A_tile           = IN.W_tile * IN.L_tile ;  

    %Calculation Models...
%      IN.Type = 'PointSrc_2D_1_over_COS3';                % 0.xx
%    IN.Type = 'PointSrc_2DLamb_Truncated';         %working (analytic)
%    IN.Type = 'PointSrc_2DLamb_alpha';                    %working (analytic)
%      IN.Type = 'PointSrc_Uniform';                    %working (analytic)
%     IN.Type = 'PointSrc_1DLambertian_1DUniform';      %working (analytic)   
%     IN.Type = 'Tile_Convolve_PointSrc_2DLamb';        %working      
%    IN.Type = 'Tile_Convolve_PointSrc_1DLambertian_1DUniform';   
%      IN.Type = 'Tile_Analytic_2DLambertian';        %  ** not ** working yet (not  urgent)
%     IN.Type = 'Tube_AnalyticInt_2DLamb';      %  Working, but assumes ALL of tube light goes downward
%   IN.Type = 'Tube_Ideal_NearField';            %working;  only valid for very very near field

    IN.Type = 'TEST_VOLUMETRIC_IRRADIANCE_TEST';  %working (analytic)


    %% (cos^alpha): SETTABLE FWHM Pattern Width LUT
    IN.lambertian_a     = get_cosine_power(IN.lambertian_FWHM);   %a = 65022.7000431579*(1./FWHM).^2.3243817615;
    
                                % show quality of model to generalize FWHM LUT ....
                                %     FWHMS = [5:5:180];
                                %     alphas_calcd = get_cosine_power(FWHMS);  %a = 65022.7000431579*(1./FWHM).^2.3243817615;
                                %     solid_angles = 2*pi ./ (alphas_calcd + 1);    %definite_integral_[0 to pi/2]_of_( cos^alpha *sin ) = 1 / (alpha+1)  
                                %     figure; plot(2*LUT.HWHM,LUT.alphas,'.');
                                %     hold all; plot(FWHMS,alphas_calcd)

                                % show quality of Model: FWHM ~  cos^alpha  
                                % thetas = [-90:90]'; FWHMsel=20; a = get_cosine_power(FWHMsel);
                                % figure; plot(thetas,(cosd(thetas)).^a);

   %% VOLUMETRIC_IRRADIANCE_TEST....
    if strcmp(IN.Type,'TEST_VOLUMETRIC_IRRADIANCE_TEST')
%         IN.testpattern = 'PointSrc_Uniform'; 
        IN.testpattern = 'PointSrc_2DLamb_alpha';
        TESTRES = TEST_VOLUMETRIC_IRRADIANCE(IN);                
        return;
    end
                     
    
    %% HORIZONTAL CROSS SECTION @ Z = ______
    Lplot_X	 = 1.5*tand(IN.lambertian_FWHM)*IN.Z;% + IN.W_tile;
    Lplot_Y	 = 1.5*tand(IN.lambertian_FWHM)*IN.Z;% + IN.L_tile;
    Nx = 100;
    [IN.XX,IN.YY] = getXXYY(Lplot_X,Lplot_Y,Nx);
      % * * * * * * * * * * * * * * * * * * * * * * *
    RES = get_II_2D(IN);   %% GET 2D IRRADIANCE PATTERN
      % * * * * * * * * * * * * * * * * * * * * * * *
  
    %gather results  %% AND NORMALIZE %%
    XX = RES.XX;        YY = RES.YY; 
    dstep = RES.dstep;  dA = dstep^2;           %[cm2]
    IIraw = RES.II;
    Praw = sum(IIraw(:))*dA
%     II = IIraw ./ (sum(IIraw(:))*dA) * IN.Power_tot;    %renormalize  (assuming all power in calculated field!!)
    II = IIraw;
    
    %PLOT
    figure; imagesc(XX(1,:),YY(:,1),II); colorbar;  axis equal
    title({'Irradiance [mW/cm2] in this plane';' ** !not accounting for cosine loss @ plane! **';...
        ['Optical Power = ',num2str(IN.Power_tot),'mW ...FWHM=',num2str(IN.lambertian_FWHM),'deg... plane: Z=',num2str(IN.Z),'cm'];...
        ['Model: ',replace(IN.Type,'_',' ')]});    
    Power_plotted = sum(II(:))*dA;
    
    
    %% GENERATE Z-Curve (Center) (uses FULL CUBE)
    Nzsteps  = 20;   Zmin_cm = 10;      	%[cm]
    LplotXY = 1;                            %[cm] size of XY meas.plane
    Zmax_cm = IN.Z;

    Nx = 3;
    [IN.XX,IN.YY] = getXXYY(Lplot_X,Lplot_Y,Nx);  

    %function: generate III_Cube
    [ICUBE,Zvals] = get_II_2DI_Cube(IN, Nzsteps, Zmin_cm, Zmax_cm, LplotXY);
    
        
    % **** plot v. Z ****
    Zvals = flipud(Zvals);
    sel.x=0;
    sel.y=0;
    
    [ dummy ,isel]  = min( sqrt( (ICUBE.XX(:)-sel.x).^2 + (ICUBE.YY(:)-sel.y).^2  )   );
    [iy,ix] = ind2sub(size(ICUBE.XX),isel);
    I_v_Z_list = ICUBE.III(iy,ix,:);    
    figure; semilogy(Zvals,I_v_Z_list(:),'.-'); xlabel('Z [cm]'); ylabel('Irradiance [mW/cm2]');
    title({'Irradiance [mW/cm2] v Z at Peak of Plane';...
        ['Optical Power = ',num2str(IN.Power_tot),'mW ...FWHM=',num2str(IN.lambertian_FWHM),'deg']});  
    
     
    
    %% 3D DATA CUBE +  Pattern Rotation
    Nzsteps     = 20;   
    Zmax_cm     = sqrt( IN.Z^2 +  max(IN.XX(:))^2 + max(IN.YY(:))^2 );         %[cm]
    Zmin_cm     = min([10,Zmax_cm-1]);                          %[cm]
    LplotXY     = max( [ abs(max(IN.XX(:))-min(IN.XX(:))) , abs(max(IN.YY(:))-min(IN.YY(:))) ] );                                            %[cm] size of XY meas.plane
    
    dX_cm = 0.5; 
    Nx = ceil(LplotXY / dX_cm);  
    
    if Nx>200 warning(['Nx = ',num2str(Nx)]); end
    [IN.XX,IN.YY] = getXXYY(Lplot_X,Lplot_Y,Nx);  
    
    
    %call  function: generate III_Cube for single pattern
    [ICUBE,Zvals] = get_II_2DI_Cube(IN, Nzsteps, Zmin_cm, Zmax_cm, LplotXY);

    
    
%     %APPLY ROTATION   --- WORK IN PROGRESS ...
%     INP.XX = RES.XX;
%     INP.YY = RES.YY;
%     INP.II = RES.II;
%     INP.Z  = IN.Z ; 
%     
%      RES = rotate_Pattern_Cube(INP,RotAngles)
    
    
    
    
    %3D Plot --- numbers need to be adjusted to give a good result
    figure; 
    %     h = slice(ICUBE.XXX,ICUBE.YYY,ICUBE.ZZZ, ICUBE.III,[0],[0],[0,200])    
    h = slice(ICUBE.XXX,ICUBE.YYY,ICUBE.ZZZ, ICUBE.III,[0],[0],[0]);
%     xlim([-20,20]); ylim([-20,20]);  zlim([0,50]);
    alpha('color')
    set(h,'EdgeColor','none','FaceColor','interp',...
    'FaceAlpha','interp')
    alphamap('increase',.5)
     alphamap('rampdown');
    
    colorbar();
    
    %volumetric average irradiance
    Rcyl_cm=100; Zcyl_cm = 200;
    ikeep_in_average = find( sqrt(ICUBE.XXX(:).^2 + ICUBE.YYY(:).^2) < Rcyl_cm  &  ICUBE.ZZZ(:)<Zcyl_cm);
    Irr_vol_average_mWcm2 = mean( ICUBE.III(ikeep_in_average) );  %[mW/cm2]
    title({['volume avg Irr: ',num2str(Irr_vol_average_mWcm2),'mW/cm2'];...
        ['Optical Power = ',num2str(IN.Power_tot),'mW ...FWHM=',num2str(IN.lambertian_FWHM),'deg']});

    
    
    
end  %end MAIN

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%
function TESTRES = TEST_VOLUMETRIC_IRRADIANCE(IN); 
    IN.Type = IN.testpattern;   %
   
    %Test Volume Definition -------------------------
        %defaults
        Zmin_cm     = 0.0;       	%[cm]
        Zmax_cm     = 304.5       %[cm]
        LplotXY     = 152.4;      %[cm] size of XY meas.plane    
        dX_cm       = 1;           %step size 
        dZ_cm       = 1;
    %confirm entries
    prompt = {'Height (Hz, in cm):  ', 'X=Y Dimension (in cm): ','dZ (in cm): ','dX=dY (in cm): '};
    dlg_title = 'Describe Volume over which to Compute Average Irradiance';
    num_lines = 1;
    defaultans = {num2str(Zmax_cm),num2str(LplotXY),num2str(dX_cm),num2str(dZ_cm)};
    answer = inputdlg(prompt,dlg_title,num_lines,defaultans);
        Zmax_cm     = str2num(answer{1});     %[cm]
        LplotXY     = str2num(answer{2});    %[cm] size of XY meas.plane    
        dX_cm       = str2num(answer{3});
        dZ_cm       = str2num(answer{4});
    
    
    %-----------------------------------------------------
    Nzsteps     = round((Zmax_cm-Zmin_cm)/dZ_cm,0);   
    Nx = ceil(LplotXY / dX_cm);  
    
    if Nzsteps>100 warning(['Nzsteps = ',num2str(Nzsteps)]); end
    if Nx>200 warning(['Nx = ',num2str(Nx)]); end
    [IN.XX,IN.YY] = getXXYY(LplotXY,LplotXY,Nx);  
    
    
    %call  function: generate III_Cube for single pattern
    [ICUBE,Zvals] = get_II_2DI_Cube(IN, Nzsteps, Zmin_cm, Zmax_cm, LplotXY);         
    
    %3D Plot --- numbers need to be adjusted to give a good result
    figure; 
    h = slice(ICUBE.XXX,ICUBE.YYY,ICUBE.ZZZ, ICUBE.III,[0],[0],[-1]);
    alpha('color')
    set(h,'EdgeColor','none','FaceColor','interp','FaceAlpha','interp')
    alphamap('increase',.2)
     alphamap('rampdown'); 
    
    colorbar();
    xlim([-50,50]);    
    ylim([-50,50]);    
    zlim([0,50]);
        
    %volumetric average irradiance
    Irr_vol_average_mWcm2 = mean( ICUBE.III(:) );  %[mW/cm2]
        %     Rcyl_cm=100; Zcyl_cm = 200;
        %     ikeep_in_average = find( sqrt(ICUBE.XXX(:).^2 + ICUBE.YYY(:).^2) < Rcyl_cm  &  ICUBE.ZZZ(:)<Zcyl_cm);
    
    %Report Result in Title
    title({['volume avg Irr: ',num2str(Irr_vol_average_mWcm2),'mW/cm2'];...
        ['Optical Power = ',num2str(IN.Power_tot),'mW ...FWHM=',num2str(IN.lambertian_FWHM),'deg'];...
        ['Test Volume(Lx, Ly, Hz): (',num2str(LplotXY),', ',num2str(LplotXY),', ',num2str(Zmax_cm-Zmin_cm),') cm'  ];...
        ['dX, dZ = (',num2str(dX_cm),', ',num2str(dZ_cm),') cm'  ]});
    

    
    %2D: Generate II_Surface for single pattern
    IN.Z = Zmax_cm;
    
      % * * * * * * * * * * * * * * * * * * * * * * *
    RES = get_II_2D(IN);   %% GET 2D IRRADIANCE PATTERN
      % * * * * * * * * * * * * * * * * * * * * * * *
  
    %gather results  %% AND NORMALIZE %%
    XX = RES.XX;        YY = RES.YY; 
    dstep = RES.dstep;  dA = dstep^2;           %[cm2]
    IIraw = RES.II;
    Praw = sum(IIraw(:))*dA
%     II = IIraw ./ (sum(IIraw(:))*dA) * IN.Power_tot;    %renormalize  (assuming all power in calculated field!!)
    II = IIraw;
    
    %PLOT2D
    figure; imagesc(XX(1,:),YY(:,1),II); colorbar;  axis equal
    title({'Irradiance [mW/cm2] in this plane';' ** !not accounting for cosine loss @ plane! **';...
        ['Optical Power = ',num2str(IN.Power_tot),'mW ...FWHM=',num2str(IN.lambertian_FWHM),'deg... plane: Z=',num2str(IN.Z),'cm'];...
        ['Model: ',replace(IN.Type,'_',' ')]});    
    Power_plotted = sum(II(:))*dA;    
    
end






%%
function  [ICUBE,Zvals]  = get_II_2DI_Cube(IN0, Nzsteps, Zmin_cm, Zmax_cm, LplotXY)
    
%     [IN0.XX,IN0.YY] = getXXYY(LplotXY,LplotXY,3);       %XY samples
    Zvals = linspace(Zmax_cm, Zmin_cm ,Nzsteps)';       %Z in DESCENDING ORDER
    for kz = 1:Nzsteps
        IN0.Z = Zvals(kz);           
          % * * * * * * * * * * * * * * * * * * * * * * *
        RES = get_II_2D(IN0);   %% GET IRRADIANCE PATTERN  (mW/cm2)
          % * * * * * * * * * * * * * * * * * * * * * * *    
          
        if kz==1    % *** SET THE X,Y SAMPLING  ***
            [Ny,Nx] = size(RES.XX);
            ICUBE.XX = RES.XX; ICUBE.YY = RES.YY;
            ICUBE.III = zeros(Ny,Nx,Nzsteps);
            ICUBE.XXX = zeros(Ny,Nx,Nzsteps);
            ICUBE.YYY = zeros(Ny,Nx,Nzsteps);
            ICUBE.ZZZ = zeros(Ny,Nx,Nzsteps);
        end
        dstep = RES.dstep;  dA = dstep^2;           %[cm2]
        Praw = sum(RES.II(:))*dA;        
        ICUBE.III(:,:,Nzsteps-kz+1) = RES.II;  %record intensity for this Z slice
        ICUBE.XXX(:,:,Nzsteps-kz+1) = RES.XX;  %record intensity for this Z slice
        ICUBE.YYY(:,:,Nzsteps-kz+1) = RES.YY;  %record intensity for this Z slice
        ICUBE.ZZZ(:,:,Nzsteps-kz+1) = IN0.Z*ones(size(RES.II));  %record intensity for this Z slice
         
    end

end


%choose the right pattern and call selected function
function RES = get_II_2D(IN)

    switch IN.Type
        case 'PointSrc_2D_1_over_COS3'
            RES = get_PointSrc_2D_1_over_COS3(IN);
        case 'PointSrc_Uniform'
            RES = get_PointSrc_Uniform(IN);                        	%working
        case 'PointSrc_2DLamb_Truncated'
            RES = get_PointSrc_2DLamb_Truncated(IN);                %working 
        case 'PointSrc_2DLamb_alpha'
            RES = get_PointSrc_2DLamb_alpha(IN);                 	%working
        case 'Tile_Convolve_PointSrc_2DLamb'
            RES = get_Tile_Convolve_PointSrc_2DLamb(IN);            %working
        case 'Tile_Analytic_2DLambertian'
            RES = get_Tile_Analytic_2DLambertian(IN);              	%not working
        case 'PointSrc_1DLambertian_1DUniform'
            RES = get_PointSrc_1DLambertian_1DUniform(IN);          %working
        case 'Tile_Convolve_PointSrc_1DLambertian_1DUniform'
            RES = get_Tile_Convolve_PointSrc_1DLambertian_1DUniform(IN); %working
        case 'Tube_AnalyticInt_2DLamb' % Assumes ALL light from "tube" is (a) Lambertian + (b) directed downward    
            RES = get_Tube_AnalyticInt_2DLamb(IN);               	%working 
        case 'Tube_Ideal_NearField'
            RES = get_Ideal_Simple_Tube(IN);
   
        otherwise
            disp(' ');disp(' ');disp(' ');
            warning('NO MODEL INPUT FOUND'); disp('NO MODEL INPUT FOUND');disp(' ');disp(' '); disp(' ');
            return;
    end

end



% FWHM = f(alpha_power)
function a = get_cosine_power(FWHM)
    a = 65022.7000431579*(1./FWHM).^2.3243817615; %generalized model
    %reference data:
    LUT.HWHM = [75, 60, 45, 38, 33, 27, 21];      %HWHM matching...
    LUT.alphas = [0.5, 1,  2,  3,  4,  6, 10];    %raising cosine to this power gives desired HWHM
end



function [XX,YY] = getXXYY(Lplot_X,Lplot_Y,Nx)
        Xv = linspace(-Lplot_X/2,Lplot_X/2,Nx);
        dstep = median(diff(Xv));
        Yv =  [-Lplot_Y/2:dstep:Lplot_Y/2]';
        Ny = length(Yv);
        [XX,YY] = meshgrid(Xv,Yv);  
end

function RES = rotate_Pattern_Cube(INP,RotAngles)
    %make 3d Data Cubes
    Nz = 10;
    iii = repmat(INP.II,[1,1,Nz]);  %build data cube of intensity
    [xxx,yyy,zzz] = meshgrid(INP.XX(1,:),INP.YY(:,1), linspace(0,INP.Z,Nz) );
    
    %apply Z rotation to input points (X-Y)
    xxx_n  = xxx*cosd(RotAngles.Z) - yyy*sind(RotAngles.Z);
    yyy_n  = xxx*sind(RotAngles.Z) + yyy*cosd(RotAngles.Z);
    xxx = xxx_n;   
    yyy = yyy_n;   
  
    %apply Y rotation to latest points  (X-Z)
    xxx_n  = xxx*cosd(RotAngles.Z) - zzz*sind(RotAngles.Z);
    zzz_n  = xxx*sind(RotAngles.Z) + zzz*cosd(RotAngles.Z);
    xxx = xxx_n;   
    zzz = zzz_n; 
    
    %apply X rotation to latest positions  (Y-Z)
    yyy = yyy*cosd(RotAngles.Z) - zzz*sind(RotAngles.Z);
    zzz = yyy*sind(RotAngles.Z) + zzz*cosd(RotAngles.Z);
    yyy = yyy_n; 
    zzz = zzz_n; 
    clear xxx_n yyy_n zzz_n;
    
    %interpolate resulting x,y,z back to original grid
    [xxx0,yy0y,zzz0] = meshgrid(INP.XX(1,:),INP.YY(:,1), linspace(0,INP.Z,Nz) );   
    Ivals = interp3(xxx,yyy,zzz, LUT.Irradiance,  xxx0, yyy0, zzz0  );

    
    

    RES.XX=XX; RES.YY=YY; RES.II=II; RES.dstep = dstep;
end



%Working
function RES = get_PointSrc_2DLamb_Truncated(IN)            %[1.1]
 
    XX = IN.XX; YY=IN.YY; 
    dstep =  median(diff(XX(1,:)));
    
    %calculate II=Irradiance
  	solidangle_sr = (2*pi) * (0.5*sind( IN.lambertian_FWHM/2 )^2);       %[sr] equal in all azimuth --> 2pi 
	thetas = atan2d((XX.^2 + YY.^2).^0.5,IN.Z);  %lamp angle = f(x,y)
 	II = IN.Power_tot .* (cosd(thetas)) ./ (solidangle_sr * ( IN.Z^2 + YY.^2 + XX.^2) ) ...
                        .* (thetas<=IN.lambertian_FWHM/2 );    % THIS APPLIES TRUNCATION at CONE ANGLE
                    
    %fix (special case if no samples are less than cone angle()):
    if sum(thetas(:)<=IN.lambertian_FWHM/2 ) == 0
        isel = find(thetas(:) == min(thetas(:))) ;
        II(isel) = IN.Power_tot .* (cosd(thetas(isel)))  ./ (solidangle_sr * ( IN.Z^2 + (YY(isel)/2).^2 + (XX(isel)/2).^2) )  / length(isel);
    end

    RES.XX=XX; RES.YY=YY; RES.II=II; RES.dstep = dstep;
end

%Working
function RES = get_PointSrc_2D_1_over_COS3(IN)                      %[0.XX]
 
    XX = IN.XX; YY=IN.YY; 
    dstep =  median(diff(XX(1,:)));
    
    %calculate II=Irradiance
    a = -3;                        %power to raise cosine to in order to get desired FWHM
     % ******** NEED TO DOUBLE CHECK SOLID ANGLE OF 1/COS3 *************
  	solidangle_sr = (2*pi) * (1 / (abs(a)+1));       % [sr] equal in all azimuth --> 2pi 
     % ******** NEED TO DOUBLE CHECK SOLID ANGLE OF 1/COS3 *************
    
	thetas = atan2d((XX.^2 + YY.^2).^0.5,IN.Z);  %lamp angle = f(x,y)
 	II = IN.Power_tot .* ((cosd(thetas)).^a) ./ (solidangle_sr * ( IN.Z^2 + YY.^2 + XX.^2) );            
    RES.XX=XX; RES.YY=YY; RES.II=II; RES.dstep = dstep;
end

%Working
function RES = get_PointSrc_Uniform(IN)                      %[1.2]
 
    XX = IN.XX; YY=IN.YY; 
    dstep =  median(diff(XX(1,:)));
    
    %calculate II=Irradiance
  	solidangle_sr = (2*pi) * (1-cosd(IN.lambertian_FWHM/2) );       %[sr] equal in all azimuth --> 2pi 
	thetas = atan2d((XX.^2 + YY.^2).^0.5,IN.Z);  %lamp angle = f(x,y)
 	II = IN.Power_tot .* (1) ./ (solidangle_sr * ( IN.Z^2 + YY.^2 + XX.^2) )...
                         .* (thetas<=IN.lambertian_FWHM/2 );    % THIS APPLIES TRUNCATION at CONE ANGLE
    RES.XX=XX; RES.YY=YY; RES.II=II; RES.dstep = dstep;
end


%Working
function RES = get_PointSrc_2DLamb_alpha(IN)                      %[1.2]
 
    XX = IN.XX; YY=IN.YY; 
    dstep =  median(diff(XX(1,:)));
    
    %calculate II=Irradiance
    a = IN.lambertian_a;                        %power to raise cosine to in order to get desired FWHM
  	solidangle_sr = (2*pi) * (1 / (a+1));       %[sr] equal in all azimuth --> 2pi 
	thetas = atan2d((XX.^2 + YY.^2).^0.5,IN.Z);  %lamp angle = f(x,y)
 	II = IN.Power_tot .* ((cosd(thetas)).^a) ./ (solidangle_sr * ( IN.Z^2 + YY.^2 + XX.^2) );            
    RES.XX=XX; RES.YY=YY; RES.II=II; RES.dstep = dstep;
end

%Working
function RES = get_Tile_Convolve_PointSrc_2DLamb(IN)        %[1.2]convolved
 
    XX = IN.XX; YY=IN.YY; 
    dstep =  median(diff(XX(1,:)));
    
    %calculate II=Irradiance
      %1. get point source irradiance first
    IN.XX = XX; IN.YY=YY;
    RES = get_PointSrc_2DLamb_alpha(IN);
    XX = RES.XX;  YY = RES.YY; 
    II_ptsrc = RES.II; 
     %2.apply convolution    
    Tile_z0 = zeros(size(XX));
    isel = find( abs(XX(:)) <  IN.W_tile  &  abs(YY(:)) < IN.L_tile );
    Tile_z0(isel) = 1;
    II = conv2(II_ptsrc,Tile_z0,'same');  %2D numeric convolution with square tile source
    
    %return data
    RES.XX=XX; RES.YY=YY; RES.II=II; RES.dstep = dstep;    
end

%NOT Working
function RES = get_Tile_Analytic_2DLambertian(IN)       %[1.3]
 
    XX = IN.XX; YY=IN.YY; 
    dstep =  median(diff(XX(1,:)));
    
    %calculate II=Irradiance : Analytic point source *PLUS* ADD AREA OF SOURCE
 	a = IN.lambertian_a;  %power to raise cosine to in order to get desired FWHM
	solidangle_sr = (2*pi) * (1 / (a+1));       %[sr] equal in all azimuth --> 2pi 
	Ltile = IN.L_tile; Wtile=IN.W_tile;
	thetas = atan2d(((abs(XX)-Wtile).^2 + (abs(YY)-Ltile).^2).^0.5,IN.Z);  % lamp angle = f(x,y)
    
	isel = find( abs(XX(:))<Wtile/2 & abs(YY(:))<Ltile/2 );      % points that are directuly beneath tile
	thetas(isel) = 0; 
 	II = IN.Power_tot .* ((cosd(thetas)).^a) ./ (solidangle_sr * ( IN.Z^2 + YY.^2 + XX.^2) + IN.A_tile);
    %return data
    RES.XX=XX; RES.YY=YY; RES.II=II; RES.dstep = dstep;    
            
end

%Working
function  RES = get_PointSrc_1DLambertian_1DUniform(IN) %[1.5]  
 
    XX = IN.XX; YY=IN.YY; 
    dstep =  median(diff(XX(1,:)));

    a = IN.lambertian_a;  %power to raise cosine to in order to get desired FWHM
	solidangle_sr = (2*pi)*(IN.alpha_tube/360) * (1 / (a+1));       %[sr] equal in all azimuth --> 2pi 
	thetas_Y = atan2d(YY,IN.Z);  %lamp angle = f(x,y)
	thetas_X = atan2d(XX,IN.Z);  %lamp angle = f(x,y)
	II = IN.Power_tot .* ((cosd(thetas_Y)).^a).*(abs(thetas_X)<=IN.alpha_tube/2)...
                      ./ (solidangle_sr * ( IN.Z^2 + YY.^2 + XX.^2) );  
    %return data
    RES.XX=XX; RES.YY=YY; RES.II=II; RES.dstep = dstep;                      
end

%Working
function RES = get_Tile_Convolve_PointSrc_1DLambertian_1DUniform(IN)
                                                        %[1.5]convolved
 
    XX = IN.XX; YY=IN.YY; 
    dstep =  median(diff(XX(1,:)));
    
    %calculate II=Irradiance
      %1. get point source irradiance first
    IN.XX = XX; IN.YY=YY;
       % ****** ONLY DIFFERENCE from get_Tile_Convolve_PointSrc_2DLamb()***
       RES = get_PointSrc_1DLambertian_1DUniform(IN);  
       %*********************************************************************
    XX = RES.XX;  YY = RES.YY; 
    II_ptsrc = RES.II; 
     %2.apply convolution    
    Tile_z0 = zeros(size(XX));
    isel = find( abs(XX(:)) <  IN.W_tile  &  abs(YY(:)) < IN.L_tile );
    Tile_z0(isel) = 1;
    II = conv2(II_ptsrc,Tile_z0,'same');  %2D numeric convolution with square tile source
    
    %return data
    RES.XX=XX; RES.YY=YY; RES.II=II; RES.dstep = dstep;
            
end

%NEEDS ATTENTION -- Assumes ALL light from "tube" is (a) Lambertian + (b) directed downward
function RES = get_Tube_AnalyticInt_2DLamb(IN)          %[1.6]
 
    XX = IN.XX; YY=IN.YY; 
    dstep =  median(diff(XX(1,:)));

    %ANALYTIC CALCULATION FOR: 1D-extended Source + dP = Lambertian in X,Y ....  
	cc = (XX.^2 + IN.Z.^2);
	L = IN.L_tube;
	II = IN.Z/1 * ( (L/2-YY)./(cc.*(cc+(L/2-YY).^2).^0.5 ) ...
                    -(-L/2-YY)./(cc.*(cc+(-L/2-YY).^2).^0.5 )   );
    dPdy = 1;
    solidangle = pi*sind(IN.FWHM);
    
	II = dPdy*IN.Z/1 * ( (L/2-YY)./(cc.*(cc+(L/2-YY).^2).^0.5 ) ...
                    -(-L/2-YY)./(cc.*(cc+(-L/2-YY).^2).^0.5 )   );                
    %return data
    RES.XX=XX; RES.YY=YY; RES.II=II; RES.dstep = dstep;
end

%working - *near field only*
function RES = get_Ideal_Simple_Tube(IN)                %[1.7]
 
    XX = IN.XX; YY=IN.YY; 
    dstep =  median(diff(XX(1,:)));
    
    %calculate II=Irradiance
    if IN.Z < IN.L_tube/2
        thetas_X = atan2d(XX,IN.Z);  %lamp angle = f(x,y) 
        II = zeros(size(XX));
        isel = find( abs(thetas_X(:)) < IN.alpha_tube  &  abs(YY(:)) < IN.L_tube/2 );
        Area_Cylinder_Rs = IN.L_tube * 2*pi*((IN.Z^2 + XX.^2).^0.5) * (IN.alpha_tube/360);
        Pdens = IN.Power_tot./Area_Cylinder_Rs;
        II(isel) = Pdens(isel);
    else
        error('This Model not valid in far field (Z > Ltube) ')
    end
    %return
    RES.XX=XX; RES.YY=YY; RES.II=II; RES.dstep = dstep;
end



    
%% **************************** NOT USED *****************************
  
    % gaussian for angular pattern .... NOT WORKING ....
        %         FWHMs = [0:10:180]';
        %         sigmas = FWHMs * pi/180 / ( 2*sqrt(2*log(2)) );
        %         solidangles_0 = 2.8*pi*( 0.5*sqrt(pi/2)*exp(-0.5*sigmas.^2).*sigmas); 
        %         figure; plot(FWHMs,solidangles_0);
        % 
        %         FWHM = 120 * pi/180;
        %         sigma = FWHM / ( 2*sqrt(2*log(2)) ); 
        %         solidangle0 = (1/ (sigma*sqrt(2*pi)) ) ;   
        %         solidangle1 = solidangle0 * ( -0.5*sqrt(pi/2)*exp(-0.5*sigma^2).*sigma );
        %         solidangle2 = solidangle1 * ( 2*(-i)*erf(i*sigma/sqrt(2)) ...
        %                                         - (-i)*erf(i* (sigma^2-i*pi/2)/(sqrt(2)*sigma)) ...
        %                                         - (-i)*erf(i* (sigma^2+i*pi/2)/(sqrt(2)*sigma)) ...


function G = GAUS(x,fwhm,x0) 
    sigma =  fwhm / ( 2*sqrt(2*log(2)) );   %log()=ln()
    G = (1/ (sigma*sqrt(2*pi)) ) * exp( -0.5*(x - x0).^2 / sigma^2 );
end
      
function RES = get_IrradianceMap(IN)
    %set up X Y samples

% %     Lplot_X	 = 5*tand(IN.alpha_tube/2)*IN.Z + sqrt(IN.W_tile);
% %     Lplot_Y	 = 2*tand(IN.alpha_tube/2)*IN.Z + sqrt(IN.L_tile) + 2*IN.L_tube*(strcmp(IN.Type,'Tube'));
%     Lplot_X	 = 1.5*tand(IN.lambertian_FWHM)*IN.Z;% + IN.W_tile;
%     Lplot_Y	 =  1.5*tand(IN.lambertian_FWHM)*IN.Z;% + IN.L_tile;
%      Nx = 100;
%     Xv = linspace(-Lplot_X/2,Lplot_X/2,Nx);
%     dstep = median(diff(Xv));
%     Yv =  [-Lplot_Y/2:dstep:Lplot_Y/2]';
%     Ny = length(Yv);
%     [XX,YY] = meshgrid(Xv,Yv);
%     II=1;
    
%     switch IN.Type
%         
%             
%         %Working
%         case 'PointSrc_2DLamb_alpha'
%             a = IN.lambertian_a;  %power to raise cosine to in order to get desired FWHM
%             solidangle_sr = (2*pi) * (1 / (a+1));       %[sr] equal in all azimuth --> 2pi 
%             thetas = atan2d((XX.^2 + YY.^2).^0.5,IN.Z);  %lamp angle = f(x,y)
%             II = IN.Power_tot .* ((cosd(thetas)).^a) ./ (solidangle_sr * ( IN.Z^2 + YY.^2 + XX.^2) );            
%             
%         % NOT WORKING  ...   
%         case 'Tile_Analytic_2DLambertian'   % NOT WORKING
%             % Analytic point source *PLUS* ADD AREA OF SOURCE
%             a = IN.lambertian_a;  %power to raise cosine to in order to get desired FWHM
%             solidangle_sr = (2*pi) * (1 / (a+1));       %[sr] equal in all azimuth --> 2pi 
%             Ltile = IN.L_tile;
%             thetas = atan2d(((abs(XX)-Ltile).^2 + (abs(YY)-Ltile).^2).^0.5,IN.Z);  % lamp angle = f(x,y)
%             isel = find( abs(XX(:))<sqrt(Ltile) & abs(YY(:))<sqrt(Ltile) );      % points that are directuly beneath tile
%             thetas(isel) = 0; 
%             II = IN.Power_tot .* ((cosd(thetas)).^a) ./ (solidangle_sr * ( IN.Z^2 + YY.^2 + XX.^2) + Ltile);
%             
%         %Working ....
%         case  'Tile_Convolve_PointSrc_2DLamb' 
%             a = IN.lambertian_a;  %power to raise cosine to in order to get desired FWHM
%             solidangle_sr = (2*pi) * (1 / (a+1));       %[sr] equal in all azimuth --> 2pi 
%             thetas = atan2d((XX.^2 + YY.^2).^0.5,IN.Z);  %lamp angle = f(x,y)
%             II_ptsrc = IN.Power_tot .* ((cosd(thetas)).^a) ./ (solidangle_sr * ( IN.Z^2 + YY.^2 + XX.^2) );
%             %apply convolution    
%             Tile_z0 = zeros(size(XX));
%             isel = find( abs(XX(:)) <  IN.W_tile  &  abs(YY(:)) < IN.L_tile );
%             Tile_z0(isel) = 1;
%             II = conv2(II_ptsrc,Tile_z0,'same');  %2D numeric convolution with square tile source
% 
%         %Working
%         case 'PointSrc_1DLambertian_1DUniform'     
%             a = IN.lambertian_a;  %power to raise cosine to in order to get desired FWHM
%             solidangle_sr = (2*pi)*(IN.alpha_tube/360) * (1 / (a+1));       %[sr] equal in all azimuth --> 2pi 
%             thetas_Y = atan2d(YY,IN.Z);  %lamp angle = f(x,y)
%             thetas_X = atan2d(XX,IN.Z);  %lamp angle = f(x,y)
%             II = IN.Power_tot .* ((cosd(thetas_Y)).^a).*(abs(thetas_X)<=IN.alpha_tube/2)...
%                             ./ (solidangle_sr * ( IN.Z^2 + YY.^2 + XX.^2) );            
%             
%        %WIP ....
%         case  'Tile_Convolve_PointSrc_1DLambertian_1DUniform' 
%             a = IN.lambertian_a;  %power to raise cosine to in order to get desired FWHM
%             solidangle_sr = (2*pi)*(IN.alpha_tube/360) * (1 / (a+1));       %[sr] equal in all azimuth --> 2pi 
%             thetas_Y = atan2d(YY,IN.Z);  %lamp angle = f(x,y)
%             thetas_X = atan2d(XX,IN.Z);  %lamp angle = f(x,y)
%             II_ptsrc =  IN.Power_tot .* ((cosd(thetas_Y)).^a).*(abs(thetas_X)<=IN.alpha_tube/2)...
%                             ./ (solidangle_sr * ( IN.Z^2 + YY.^2 + XX.^2) );  
%             %apply convolution    
%             Tile_z0 = zeros(size(XX));
%             isel = find( abs(XX(:)) <  IN.W_tile &  abs(YY(:)) <  IN.L_tile );
%             Tile_z0(isel) = 1;
%             II = conv2(II_ptsrc,Tile_z0,'same');  %2D numeric convolution with square tile source
%             
%             
%         case 'Tube_AnalyticInt_2DLamb'
%              %Tube =>> 1D-extended Source; dP = Lambertian in X,Y ....  
%             cc = (XX.^2 + IN.Z.^2);
%             L = IN.L_tube;
%             II = IN.Z/1 * ( (L/2-YY)./(cc.*(cc+(L/2-YY).^2).^0.5 ) ...
%                            -(-L/2-YY)./(cc.*(cc+(-L/2-YY).^2).^0.5 )   );
%             II = II./sum(II) * IN.Power_tot;
%             
%             
%         case 'Tube_NumericInt_1DLamb'    
%                %Tube  =>> 1D-extended Source; dP = Lambertian in Y,  Uniform in X ....  
%             dL = 1;
%             Nsteps = round(IN.L_tube / dL,0);
%             dP = IN.Power_tot / Nsteps;
% 
%             II = zeros(size(XX));
%             for k=1:Nsteps
%                 yk = -IN.L_tube/2 + k*dL
%                 %lambertian
%                 II = II + ( dP * cosd(atand((YY-yk)./IN.Z)) ./(2*deg2rad(IN.alpha_tube)*(IN.Z^2 + (YY-yk).^2+ XX.^2))  );        
%             end
% 
%             IIa = IN.Power_tot ./(deg2rad(IN.alpha_tube)*IN.L_tube*(IN.Z^2 + XX.^2).^0.5).*(abs(YY)<IN.L_tube/2 .* abs(atand(XX./IN.Z))<IN.alpha_tube/2); 
%     %         IIb = IN.Power_tot*cosd(atand(YY./IN.Z)) ./(2*deg2rad(IN.alpha_tube)*(IN.Z^2 + YY.^2+ XX.^2));
%     %         II = min(cat(3,IIa,IIb),[],3);
%             figure; imagesc(XX(1,:),YY(:,1), II); axis equal
%             
%     end %switch        

    RES.XX = XX; RES.YY = YY;  RES.II = II;
    RES.dstep = dstep;
    
end %calc_II

            