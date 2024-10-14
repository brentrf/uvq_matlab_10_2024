
%% Blah
function Analyze_ICUBE()
    close all;
    clear;
    dirs.f  = '/Users/brentfisher/Documents/MATLAB/UVQ';
    start_path = '/Users/brentfisher/Documents/MATLAB/UVQ/irradiance_LUTS/';
   	[filename,path] = uigetfile(start_path);
      
%     path  =  '/Users/brentfisher/Documents/MATLAB/UVQ/irradiance_LUTS/UniformCone_30deg _Ptot1W_N=1_Tilt=0_Lx304_Ly304cm/';
%     filename = 'ICUBE.mat';
%     
    %load ICUBE
    load([path,'/',filename]);
    pathsplit = split(path,'/');
    ICUBE
    
    %eACH  & Wells Riley Inputs
    phi_pathogen_mJcm2perlog = 1/3;     %D90 =  dose required for 1 log kill (pathogen dependent)
    phi_pathogen_mJcm2perlog = 0.4;     %D90 = dose required for 1 log kill (pathogen dependent)
    F_AC = 0.63;                        %definition of 1 Air Change: fraction of volume exchanged for 
    
    
    %calculate Average Irradiance & eACH FOR DESIRED VOL & PTOT
    res = ICUBE_Iavg(ICUBE);
    res.eACH_avg =  res.Irr_vol_average_mWcm2 * 3600 / (phi_pathogen_mJcm2perlog * (-log10(1-F_AC)) )  ;
    res.D90   = phi_pathogen_mJcm2perlog; 
  
    disp(pathsplit{end});
    disp(['I_avg_vol =',num2str(res.Irr_vol_average_mWcm2),'mW/cm2']);
    disp(['eACH_full_room = ',num2str(res.eACH_avg)]);
    disp(['D90[assumed] = ',num2str(res.D90),'mJ/cm2 per log10_kill ']); 
    
    disp(res.titletext);
    
    %3D Plot
    dummy= ICUBE_3D_Plot(ICUBE,res.titletext);
    
    %2D Plot
    Z_plane_to_plot_cm = 225-170;
    dummy = ICUBE_2DPlot(ICUBE, Z_plane_to_plot_cm, res );
    Z_plane_to_plot_cm = 225-100;
    dummy = ICUBE_2DPlot(ICUBE, Z_plane_to_plot_cm, res );

    
    % ------------------------------------------------------------
    %Wells Riley - assumptions from "COVID-19_Aerosol_Transmission_Estimator"
    Ptot_W_thiscase = 0.15;            %[W] for this particular case
    rate_new_qph = 10;        %[quanta/h] rate of new pathogen added (1 quanta = 1 infective dose)
    rate_kill_nonUV = 1.62;     %[1/h]  sum of ( air[0.62] + gravity[0.3] + ventilation[0.7] )
    BreathRate_mLph = 4500*60;         %[mL/h] www.mdapp.co/minute-ventilation-equation-calculator-416/
    ICUBE.rate_kill_UV_xyz ...        %[1/h] 
            = (Ptot_W_thiscase / ICUBE.Ptot_W)*(ICUBE.III_mWcm2 / phi_pathogen_mJcm2perlog) * 3600 * log(10);  %log=natural log
    rate_kill_UV_avg ...        %[1/h] 
            = (res.Irr_vol_average_mWcm2 / phi_pathogen_mJcm2perlog) * 3600 * log(10);  %log=natural log
        
    Tresidence_h = 8;             %[h] time that person is present
    Vratio = (ICUBE.dxyz^3) / (ICUBE.Lx*ICUBE.Ly*ICUBE.Lz);
    
    Prob_ref_noUV  = 1 - exp(-1*(rate_new_qph*Vratio*(rate_kill_nonUV).^-1 ...
                            * BreathRate_mLph *Tresidence_h )); 
    %high Mix Limit:  Probability 
    Prob_avg_high_mix  = 1 - exp(-1*(rate_new_qph*Vratio*(rate_kill_nonUV+ rate_kill_UV_avg).^-1 ...
                            * BreathRate_mLph *Tresidence_h )) 
    %CONSERVATIVE Limit:  calc probability of each voxel, THEN volume average
    ICUBE.Prob  = 1 - exp(-1*(rate_new_qph*Vratio*(rate_kill_nonUV+ICUBE.rate_kill_UV_xyz).^-1 ...
                            * BreathRate_mLph *Tresidence_h )); 
    isel = find( ~isnan(ICUBE.Prob(:)) );                 
    Prob_xyz_avg = mean(ICUBE.Prob(isel));  %CONSERVATIVE
 	
    disp(['Prob_WR (ref,noUV):    ',num2str(Prob_ref_noUV)]);
    disp(['Prob_WR (@High Mix):   ',num2str(Prob_avg_high_mix)]);
    disp(['Prob_WR (conservative):',num2str(Prob_xyz_avg)]);
    
    %% 
    %calculate - ICUBE.eACH cube
    
    ICUBE.eACH  =  ICUBE.III_mWcm2 * 3600 / (phi_pathogen_mJcm2perlog * (-log10(1-F_AC)) )  ;
    
    %calculate 
    
    
    
end


    

   
   
    
%% 3D Plot
function hf = ICUBE_3D_Plot(ICUBE,titletext) 
    hf = figure; 
    h = slice(ICUBE.XXX,ICUBE.YYY,ICUBE.ZZZ, ICUBE.III_mWcm2,[0],[0],[-1]);
    alpha('color')
    set(h,'EdgeColor','none','FaceColor','interp','FaceAlpha','interp')
    alphamap('increase',.2)
     alphamap('rampdown'); 
    colorbar();
    xlim([-50,50]);    
    ylim([-50,50]);    
    zlim([0,50]);
    xlabel('mW/cm2');
    title(titletext);
    

end

%% 2D Plot --- 
function dummy = ICUBE_2DPlot(ICUBE,Z_plane_to_plot_cm,result)
    
    %Z_plane_to_plot_cm = 304.5;
    
    [dummy,isel] = min( abs(ICUBE.ZZZ(1,1,:) -  Z_plane_to_plot_cm) );
    XX = ICUBE.XXX(:,:,isel); YY=ICUBE.YYY(:,:,isel); 
    II = ICUBE.III_mWcm2(:,:,isel) *  (result.newVals.Ptot_W / ICUBE.Ptot_W);
    iselctr7cm = find( sqrt(XX(:).^2+YY(:).^2)<3.5);
    IIavg_ctr_mWcm2 = mean(II(iselctr7cm));
    
    figure; imagesc(ICUBE.XXX(1,:,isel),ICUBE.YYY(:,1,isel),ICUBE.III_mWcm2(:,:,isel)); colorbar;  axis equal
    title({['Irradiance [mW/cm2] in this plane: (z=',num2str(Z_plane_to_plot_cm),')'] ;' ** !not accounting for cosine loss @ plane! **';...
        ['Optical Power = ',num2str(result.newVals.Ptot_W*1000),'mW ...FWHM=','','deg... plane: Z=',num2str(Z_plane_to_plot_cm),'cm'];...
        ['Avg Irradiance @ Ctr: ',num2str(IIavg_ctr_mWcm2*1e3),'uW/cm2']});    
    Power_mW_this_nplane = sum(II(:))*(ICUBE.dxyz)^2 
end


function result = ICUBE_Iavg(ICUBE)
     %Define VOLUME DIMENSIONS for AVERAGING -------------------------
    Zmin_cm = 0.0;
    Zmax_cm = ICUBE.Lz;
    LplotX = ICUBE.Lx;  %full width 
    LplotY = ICUBE.Ly;  %full width  
    Ptot_W  = ICUBE.Ptot_W;  %Ptot    
    prompt = {'Enter *minimum* (Z) in [cm]:  ','Enter *maximum* (Z) in [cm]:  ','Enter CUBE WIDTH(X) in [cm]:  ',...
              'Enter CUBE LENGTH(Y) in [cm]:','source POWER in [W]:'};
    dlg_title = 'VOL_to_AVG';
    num_lines = 1;
    defaultans = {num2str(Zmin_cm), num2str(Zmax_cm),num2str(LplotX),num2str(LplotY),num2str(Ptot_W)};
    answer = inputdlg(prompt,dlg_title,num_lines,defaultans);
    Zmin_cm = str2num(answer{1});
    Zmax_cm = str2num(answer{2});
    LplotX = str2num(answer{3});
    LplotY = str2num(answer{4});
    Ptot_W = str2num(answer{5}); 
    result.newVals.LplotX = LplotX;
    result.newVals.LplotY = LplotY;
    result.newVals.Ptot_W = Ptot_W;    
    
    % ikeep_in_average = find(~isnan(ICUBE.III_mWcm2(:)));
    ikeep_in_average = find(~isnan(ICUBE.III_mWcm2(:)) & abs(ICUBE.XXX(:))<LplotX/2 & abs(ICUBE.YYY(:))<LplotY/2  &  ICUBE.ZZZ(:)< Zmax_cm  &  ICUBE.ZZZ(:) > Zmin_cm);
    result.Irr_vol_average_mWcm2 = mean( ICUBE.III_mWcm2(ikeep_in_average) ) *  Ptot_W / ICUBE.Ptot_W;  %[mW/cm2]
    
    %Title Text
    result.titletext = {['volume avg Irr: ',num2str(result.Irr_vol_average_mWcm2),'mW/cm2'];...
        ['Optical Power = ',num2str(Ptot_W*1000),'mW ...ZemaxFiles=','','deg'];...
        ['Test Volume(Lx, Ly, Hz): (',num2str(LplotX),', ',num2str(LplotY),', ',num2str(Zmax_cm-Zmin_cm),') cm'  ];...
        ['dX, dZ = (',num2str(ICUBE.dxyz),', ',num2str(ICUBE.dxyz),') cm'  ]};
    
end