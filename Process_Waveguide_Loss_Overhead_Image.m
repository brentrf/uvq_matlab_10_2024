% TOP DOWN WAVEGUIDE Image Analysis
% 6 May 2024


function A = Process_Waveguide_Loss_Overhead_Image( )
    %clear & get files
    clear; close all;
    program_name = 'Process_Waveguide_Loss_Image'
    dirs.d = 'G:\Shared drives\Corp Main\Engineering\LAB EXPERIMENTS\';
%     cd(dirs.d); [FileName,PathName] = uigetfile({'*.mat';'*.tiff'},program_name);
    cd(dirs.d); [FileName,PathName] = uigetfile({'*.tiff'},program_name);
    filenametitle = replace(FileName,'_',' '); 
 
    %open file 
    if strcmp(FileName(end-3:end),'tiff')  %for TIFF file input
        I0 = imread([PathName,FileName]);
        [Nr,Nc,Nlayers] = size(I0);
        if (Nlayers>=3)
            I1 = I0(:,:,3);                  %take Blue layer
        else
            I1 = I0;
        end
    elseif strcmp(FileName(end-3:end),'.mat')   %for HDR (*.mat) file input
        INDAT = load([PathName,FileName]);
        I0 = uint32( INDAT.I_HDR );
        [Nr,Nc,Nlayers] = size(I0);   
        I1 = I0(:,:,1);            
    end
    
    %show image and define points of (linear)Waveguide
    hf1 = figure; 
    %himg  =  imshow(I0); 
    himg  = imagesc(log10(double(I1))); colorbar;
    h = msgbox('Click origin(1) and end(2) of waveguide',program_name);
    uiwait(h);
    [Ps_x,Ps_y] = ginput(2);
    
    %How long is the segment defined?
    WG_Selection_Length_cm = 0.2; %[cm]
    answer = inputdlg('What length is selected WG segment [cm]',program_name,1,{num2str(WG_Selection_Length_cm)})
    WG_Selection_Length_cm = str2num(answer{1});
    cm_per_pixel =  WG_Selection_Length_cm / sqrt ( (Ps_x(2)-Ps_x(1))^2 +  (Ps_y(2)-Ps_y(1))^2  );
    
    
    %how close must pixel be to qualify?
    Dmax = 5;  %pixels
    answer = inputdlg('Dmax: how many pixels away from line should define signal?',program_name,1,{num2str(Dmax)})
    Dmax = str2num(answer{1});
    
    
    
    %calculate distance (DD)  of each pixel from waveguide
    [XX,YY] = meshgrid([1:Nc],[1:Nr]);
    
    DD = abs( (Ps_x(2) - Ps_x(1)) * (YY - Ps_y(1))  - (XX - Ps_x(1))*(Ps_y(2) - Ps_y(1))  ) ...
        /  sqrt( (Ps_x(2) - Ps_x(1))^2 + (Ps_y(2) - Ps_y(1))^2 );
    
    iix = find (  XX(:) <= min(Ps_x)-Dmax |  max(Ps_x)+Dmax <= XX(:) ) ;
        DD(iix) = Dmax*10;
    iiy = find (   YY(:) <= min(Ps_y)-Dmax |  max(Ps_y)+Dmax <= YY(:) ) ;
        DD(iiy) = Dmax*10;
        
    hf2 = figure;
    imagesc(DD); 
    
    %select specific pixels  -- this code *SHOULD* work .. but it doesnt
    %     isel = find(DD(:) <= Dmax);
    %     [rsel,csel] = ind2sub(size(DD),isel);
    %     length(isel);
    %     I1 = I0; I1(rsel,csel,[1,2])=100;
    %     h_img  =  imshow(I1); 
    
    %select specific pixels 
    I2 = cat(3, zeros(Nr,Nc),zeros(Nr,Nc),I1);
    I1 = I0; figure; ii=1;
    for kr=1:Nr
        for kc=1:Nc
            if(DD(kr,kc) <=  Dmax ) % -- WAVEGUIDE
                I2(kr,kc,[1,2])=100;
                rsel(ii,1)=kr;
                csel(ii,1)=kc;
                R_pix_fromstart(ii,1) = sqrt((kc-Ps_x(1))^2 + (kr-Ps_y(1))^2); %[pix]
                R_cm_fromstart(ii,1) = R_pix_fromstart(ii,1) * cm_per_pixel;  %[cm]
                
                Val_blue(ii,1) = double(I2(kr,kc,3));                %signal (within Dmax)
                Baseline_blue(ii,1) = mean(I2(kr+2*Dmax:kr+3*Dmax,kc,3));  %baseline (at Dmax x 2
                ii=ii+1;
            end
        end
    end
    h_img  =  imshow(I2);   %show's which pixels are selected to include
    isel = sub2ind(size(DD),rsel,csel);
    
    %Sort by position along waveguide ( all Pixels used for analysis )
    [Rs,I] = sort(R_pix_fromstart);
    [Rs,I] = sort(R_cm_fromstart)
    Val_blue_sorted = Val_blue(I);
    Baseline_blue_sorted = Baseline_blue(I);

    figure; plot(Rs, Val_blue_sorted); 
    hold on; plot(Rs, Baseline_blue_sorted); legend('Value','Baseline')
    
    figure; plot(Rs, Val_blue_sorted - Baseline_blue_sorted); 
    title('Net = Value - Baseline');
    
    % Smooth the data
%     yy = smooth(   Val_blue_sorted - Baseline_blue_sorted, 100);
    yy  = window_smooth(Val_blue_sorted - Baseline_blue_sorted, 100);
%     yy = Val_blue_sorted - Baseline_blue_sorted;
    figure; plot(Rs,yy); title('Smoothed (Net = Value - Baseline)');
     
    %Fit to  10* Log10(Intensity)
    fit.xmin_cm = 0.000;
    fit.xmax_cm = 0.10;
    
    ifit = find( fit.xmin_cm <= Rs & Rs <= fit.xmax_cm & yy >= 0);
    
    [fit.p,fit.S] = polyfit( Rs(ifit), 10*log10(yy(ifit)), 1);
    fit.log10yydB = polyval(fit.p,Rs);
    fit.yy = 10.^(fit.log10yydB / 10);
    
    figure; plot(Rs,[10*log10(yy),fit.log10yydB] );     
        hold all; plot(Rs(ifit),[10*log10(yy(ifit)),fit.log10yydB(ifit)] );
        xlabel('X, [cm]'); ylabel('10*log10(y) [dB]')    
         title({[filenametitle]; ['fit value  = ',num2str(fit.p(1)),'dB/cm'] })
         
    figure; plot(Rs, [yy,fit.yy]); 
        hold all; plot(Rs(ifit), [yy(ifit),fit.yy(ifit)]); 
        xlabel('X, [cm]'); ylabel(' (y) [a.u.]')   
         title({[filenametitle]; ['fit value  = ',num2str(fit.p(1)),'dB/cm'] })

    
    
end

function  y = window_smooth(y0,N)
    y = y0;  %initialize
    Nelements = size(y,1);
    for k=1:Nelements
        y(k,1) = mean( y0( [ max([1,k-round(N/2,0)]): min([k+round(N/2,0),Nelements]) ],  :  )  ,1 );
    end


end

