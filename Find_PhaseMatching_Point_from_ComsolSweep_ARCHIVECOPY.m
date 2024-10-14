% 1. Open  Mode Indice Data from XLS
%  Expected Input Format:
%    cols = [WLs(Pparam2=1),  WLs(@Param2), .... ]
%    row1 = WLs values
%    row2 = Param2 values
%    rows = Modes   
%
%    ** Rows must be SORTED so that INDEX is DECREASING top to bottom ***
%
%    
% Key Parameters to Update (in hard code):
%  a.)   Nfund = 2 ... how many of top rows to assign as "fundamental"
%                  ... The fundamental modes will be assessed for crossings 
%                  at their given wavelength.  ( Keep this below 2 or 3, else # plots will explode)
%                  ... All other modes (rows below Nfund) will have
%                  wavelength doubled (2x) in order to look for crossing
%                  with fundamentals.
%  b.) dWL = ____   ... precision of WL in the identifying crossing
%  c.) methodstr = 'linear';  ... interpolation method 
%  d.) 'extrap'    ... whether or not to exp
%  
% 
% USES "intersections() function
%    ....  [x0,y0,iout,jout] = intersections(x1,y1,x2,y2,robust)  
%


function Find_PhaseMatching_Point_from_ComsolSweep()
    close all;
    clear;
    addpath 'G:\My Drive\Analyses(BF)\Matlab\UVQ';
    program_name = 'find phase matching point'
    
%% KEY INPUT PARAMETERS 
    Nfund = 2; %designate these two rows as "Fundamental"
    methodstr = 'spline'; %'linear'  'makima' 'pchip'
    dWL = 0.1; 
    P2_sel_vals_extra = [ ];  %leave this empty if you don't want to run additional values it....
    
%% get data
    dirs.f = 'G:\My Drive\Analyses(BF)\Matlab\UVQ';
    dirs.d = 'G:\My Drive\Analyses(BF)\Matlab\matDat';
    
    cd(dirs.d); [FileName,PathName] = uigetfile({'*.xlsx'},program_name);
    filenametitle = replace(FileName,'_',' '); 
    
%% read data    
    %rows = different mode indices 
       %assume HIGHEST to LOWEST mode orders
    %columns = different parameters (wavelength, geometry, etc)
    A = xlsread([PathName,FileName],'Sheet1','a2:ff50');   %no more than 50 modes
%     A = xlsread([PathName,FileName],'real_part_only','a2:ff50');   %no more than 50 modes
    % complex values should not be read (appear as "NaN"
    % first two rows are parameter values (for each column)

%% pull out the parameter values
    WLs = A(1,:);
    P2s = A(2,:);
    WLs_u = unique(WLs(~isnan(WLs)));    Nvals1 = length(WLs_u);    %Wavelength values
    P2s_u = unique(P2s(~isnan(P2s)));    Nvals2 = length(P2s_u);    %number of param1 values
    Nm = size(A,1)-2;                                               %number of modes
    
%% Stack to 3D Cube on Parameter2
    B = zeros(size(A,1)-2, Nvals1, Nvals2);
    for k = 1:Nvals2
        p2 = P2s_u(k);
        icols_sel = find(P2s==p2)
        Nsel = length(icols_sel);
        B(:,1:Nsel,k) = A(3:end,icols_sel);
    end
        
%% Find WL Crossing of Fundamental Modes with Other Modes
    Nfund;  %(input above) ... designate these two rows as "Fundamental"
    C = B(Nfund+1:end,:,:);   %these modes will have WL doubled for comparison to fundamental dispersion curve

    

%% Do Checks -- Check Size of Output
    P2_sel_vals = [ P2s_u, P2_sel_vals_extra  ];  %These are the Parameter 2 Values where crossings (on param1 will be reported

    Nplots_expected = Nvals2 * Nfund;
    N_shg_modes     = Nm - Nfund;
    qstring = {['Num plots expected = ',num2str(Nplots_expected),'  (<15 is best) '  ];...
                ['Num lines on plots = ',num2str(N_shg_modes) ,'  (<20 is best)' ];...
                'Proceed ?' };
    button = questdlg(qstring,'Output Check','Yes','No','Yes')
    
    if min(P2_sel_vals) < min(P2s_u) | max(P2_sel_vals) > max(P2s_u)
        qstring = {[' Desired Parameter2 Values ( ',num2str(P2_sel_vals),'  ) '  ];...
                   [' are outside the range provided by comsol results: ( ',num2str(P2s_u) ,'  )' ];...
                    'Proceed ?' };
         button = questdlg(qstring,'Output Check','Yes','No','Yes')   ;             
    end
        
%% INTERPOLATION (w/ EXTRAPOLATION) - to find crossings

    close all; 
    WLs = [min(WLs_u):0.1:max(WLs_u)]';
    WLs = [200:dWL:max([WLs_u,1100])]';
    
%     WLs = [200:dWL:1500]';
    hf0=figure;

%%   *** Find Crossing along WL(param1) at each DESIRED Param2 Value ***
    
    % ** Find Crossing along  WL(Param1) at each USER DEFINED Param2  **
    for kf=1:Nfund                              %LOOP:  fundamental modes
        for kp=1:length(P2_sel_vals)            %LOOP:  over Param2
            legnames = [];
            
            P2value = P2_sel_vals(kp);
            
            if ismember(P2value,P2s_u)
                %No P2 interpolation ... 
                n_fund_points_tmp = B(kf,:,find(P2value==P2s_u));
            else
                %Interp over P2:  get Index v WL points at desired P2 Value by 2D Interpolation --- Vq = interp2(X,Y,V,Xq,Yq);    
                tmp = reshape(  B(kf,:,:) , Nvals1, Nvals2);  %rows = WLs, cols=Param2
                n_fund_points_tmp =  interp2(P2s_u,   WLs_u,   tmp,  P2value, WLs_u); % OR: griddata(P2s_u,   WLs_u,   tmp,  P2value, WLs_u);            
            end
            
            %Interp over P1 (WLs) with Extrapolation 
            isel_isan = find(~isnan(n_fund_points_tmp ));
            n_fund = interp1(WLs_u(isel_isan), n_fund_points_tmp(isel_isan), WLs, methodstr,'extrap') ;

            
            %find WL crossing at each Param2 Value (using 1D interpolation)
            for kc=1:size(C,1)                  %LOOP: non-fundamental modes
                
                if ismember(P2value,P2s_u)
                    %No P2 interpolation 
                    n_shg_points_tmp = C(kc,:,find(P2value==P2s_u));     
                else
                    %Interp over P2:
                    tmp = reshape(  C(kc,:,:) , Nvals1, Nvals2);  %rows = WLs, cols=Param2
                    n_shg_points_tmp = interp2(P2s_u,   WLs_u,   tmp,  P2value, WLs_u);
                    %                     n_shg_points_tmp = griddata(P2s_u,   WLs_u,   tmp,  P2value, WLs_u); % OR:
                
                end
                
                % ** INTERSECTION FINDER **
                %Interp over P1 (WLs) with Extrapolation  & Find Crossing ...
%                 [found_crossing,cross_interp_WL,cross_interp_n, n_shg] = find_crossing(WLs,n_fund, WLs_u, n_shg_points_tmp, methodstr);
                
%                 x1 = WLs;
%                 y1 = n_fund;
                x1 = WLs_u;
                y1 = n_fund_points_tmp;
                x2 = 2*WLs_u;
                y2 = n_shg_points_tmp;
                
%                 hf0=figure; 
%                 plot(x1,y1,'o-',x2,y2,'o-');
%                 title(['CHECK INPUT LINES kf = ',num2str(kf),'    kp = ',num2str(kp),...
%                         ' (',num2str(P2_sel_vals(kp)),') ']);
                robust = 0;
                [x0,y0,iout,jout] = intersections(x1,y1,x2,y2,robust);
                cross_interp_WL = x0;
                cross_interp_n = y0;
                found_crossing = length(x0)>0;
                
%                 close(hf0);
                
                                                 
                %record crossing result
                cross_lam_interp(kc,kp,kf) =  max([-1,cross_interp_WL]);   % [ row=SHGmode;  ool=Param2;  layer=FundMode ]
                cross_n_interp(kc,kp,kf)   = max([-1,cross_interp_n]);
              
                
                %FIGURE - for checking results
                if kc==1 
                    hff(kf,kp)  = figure(kp + Nvals2*(kf-1));  
%                     plot(WLs,n_fund,'.'); legnames{1}='fund';
                    plot(WLs_u,n_fund_points_tmp,'o-'); legnames{1}='fund';
                end
                if found_crossing==1
%                     figure(hff(kf,kp));  hold all; plot(WLs,n_shg,'-');  hold on;  plot(2*WLs_u,n_shg_points_tmp,'o');
%                     Nlegnames = length(legnames); legnames{Nlegnames+1}=num2str(kc);   legnames{Nlegnames+2}=num2str(kc);                       
                    figure(hff(kf,kp));  hold all; plot(2*WLs_u,n_shg_points_tmp,'o-'); % hold on;  plot(2*WLs_u,n_shg_points_tmp,'o');
                    Nlegnames = length(legnames); legnames{Nlegnames+1}=num2str(kc);   legnames{Nlegnames+2}=num2str(kc);                       
                    hold on;  plot(cross_interp_WL,cross_interp_n,'^','Color','b'); 
                    Nlegnames = length(legnames); legnames{Nlegnames+1}=[num2str(kc),'-X'];
                end
                legend(legnames);
                title([' kf = ',num2str(kf),'    kp = ',num2str(kp),...
                        ' (',num2str(P2_sel_vals(kp)),') ']);
 
                
                %Progress indicator 
                disp(['kf=',num2str(kf),'  kp=',num2str(kp),'  kc=',num2str(kc)]);
                
            end %loop: shg modes
            
        end     %loop: parameter 2
        
        
        disp(cross_lam_interp(:,:,kf));
        disp(cross_n_interp(:,:,kf));
        disp('%(-1): only one point, cant interpolate  n v WL')
        disp('%(-2): No Intersection within domain')
        disp('%(-3): No Sign Change with Crossing')
        
    end         %loop fundamental modes
    
    
    
%%   *** Find Crossing along WL(param1) at each DESIRED Param2 Value ***

%     for kf=1:Nfund
%         %Fundamental Curves 
%         tmp = reshape(  B(kf,:,:) , Nvals1, Nvals2);  %rows = WLs, cols=Param2
%             %2D Interpolation --- Vq = interp2(X,Y,V,Xq,Yq);    
%         n_fund_2d = interp2(P2s_u,   WLs_u,   tmp,  P2_sel_vals,  WLs);
%         
%         for kp=1:length(P2_sel_vals)                         %loop:  over Param2
%             legnames = [];
%             isel_isan = find(~isnan(B(kf,:,kp)));
%             n_fund = interp1(WLs_u(isel_isan), B(kf,isel_isan,kp), WLs, methodstr,'extrap') ;
%             
%             %find crossing at each Param2 Value (using 1D interpolation)
%             for kc=1:size(C,1)                  %loop: non-fundamental modes
%                 C_indices = C(kc,:,kp);
%                 [found_crossing,cross_interp_WL,cross_interp_n, n_shg] = find_crossing(WLs,n_fund, WLs_u, C_indices, methodstr);
% 
%         [found_crossing,cross_interp_WL,cross_interp_n, n_shg] = find_crossing(WLs,n_fund, WLs_u, C_indices, methodstr);
% 
%         
% %         isel_isan = find(~isnan(tmp(:) ) );
% %         [isel_rows,isel_cols] = ind2sub(size(tmp), isel_isan);
% %         disp([isel_rows,isel_cols]);
%         
%         
%     end
    
    
    
    

end %END --  MAIN FUNCTION
    
function [found_crossing,cross_interp_WL,cross_interp_n, n_shg] = find_crossing(WLs,n_fund, WLs_u,C_indices, methodstr)
                n_shg=[];
                found_crossing=0;  %reset      
                isel_isanC = find(~isnan(C_indices));
                
                %INTERPOLATION INTERSECTION
                if length(isel_isanC)<=1
                    cross_interp_WL = -1;                      %(-1) only one point, can't interpolate  n v WL
                    cross_interp_n  = -1;
                else
                    n_shg  = interp1(2*WLs_u(isel_isanC), C_indices(isel_isanC), WLs, methodstr, 'extrap') ;
                    [mindiff, imin]  = min(abs(n_shg - n_fund));
                    if ~(1 < imin & imin < length(n_shg)) %imin==1 | imin==length(n_shg)        
                        cross_interp_WL = -2;                  %(-2): No Intersection within domain
                        cross_interp_n  = -2;
                    else
                        if sign(n_shg(imin-1)-n_fund(imin-1)) ~= sign(n_shg(imin+1)-n_fund(imin+1))
                            cross_interp_WL = WLs(imin);
                            cross_interp_n = n_shg(imin);
                            found_crossing = 1;
                        else
                            cross_interp_WL = -3;              %(-3): no Crossing Detected
                            cross_interp_n = -3;
                        end
                    end
                end  
                
end
    
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%  CODE BELOW IS NOT WORKING .... 
%%  POLY FIT -- to find crossings

% 
%     for kf=1:Nfund
%         for kp=1:Nvals2 %loop over Param2
%             isel_isan = find(~isnan(B(kf,:,kp)));
%             [fund.a, fund.S] = polyfit(WLs_u(isel_isan), B(kf,isel_isan,kp), 2);
%             
%             for kc=1:size(C,1)  %loop over non-fundamental modes
%                  isel_isanC = find(~isnan(C(kc,:,kp)));
%                 
%                 %POLYFIT: intersection wavelength (solve quadratic)
%                 if length(isel_isanC)>=2  
%                     if length(isel_isanC)>=3   %(quadratic)
%                         [shg.b, shg.S] = polyfit(2*WLs_u(isel_isanC), C(kc,isel_isanC,kp), 2);
%                     else  %linear
%                         shg.b(3) = shg.b(2);
%                         shg.b(2) = shg.b(1);
%                         shg.b(1) = 0;   
%                     end
%                     %(solve quadratic)
%                     % !! FIND AN ONLINE CODE FOR HIGHER ORDER POLYNOMIAL
%                     % INTERSECTION POINTS !!
%                     crossfit.lam1 = ( -(fund.a(2)-shg.b(2))  + sqrt( (fund.a(2)-shg.b(2))^2 - 4*(fund.a(1)-shg.b(1))*(fund.a(3)-shg.b(3)) ) ) / ...
%                                     ( 2*(fund.a(1)-shg.b(1))  );
%                     crossfit.lam2 = ( -(fund.a(2)-shg.b(2)) - sqrt( (fund.a(2)-shg.b(2))^2 - 4*(fund.a(1)-shg.b(1))*(fund.a(3)-shg.b(3)) ) )/ ...
%                                     ( 2*(fund.a(1)-shg.b(1))  );
%                     %                 cross.lam = -(a1-b1) + sqrt((a1-b2)^2 - 4*(a0-b0)*(a2-b2) ) / ...
%                     %                                 ( 2*(a0-b0)  );
%                     
%                 else
%                     crossfit.lam1 = 0;
%                     crossfit.lam2 = 0;
%                 end
% 
%                 % figure;   %for debugging & checking result
%                 plot(WLs,polyval(fund.a,WLs)); hold all; plot(WLs,polyval(shg.b,WLs));
%                 title(['kf=',num2str(kf),'  kp=',num2str(kp),'  kc=',num2str(kc)]);
%                 hold off;
%                 %pause(0.1); 
% 
%                 %record crossing result
%                 if isreal(crossfit.lam1) & WLs(1) <= crossfit.lam1 & crossfit.lam1 <= WLs(end)  
%                     crosslam_fit1(kf,kp,kc) =  crossfit.lam1;
%                 else crosslam_fit1(kf,kp,kc) = -3;
%                 end
%                 if isreal(crossfit.lam2) & WLs(1) <= crossfit.lam2 & crossfit.lam2 <= WLs(end)  
%                     crosslam_fit2(kf,kp,kc) =  crossfit.lam2;
%                 else crosslam_fit2(kf,kp,kc) = -3;
%                 end
%                 
%                 disp(['kf=',num2str(kf),'  kp=',num2str(kp),'  kc=',num2str(kc)]);
%                 
%             end %loop: shg modes
%         end     %loop: parameter 2
%     end         %loop fundamental modes


%% incomplete

    
   
    
