% 1. Open  Mode Indice Data from XLS
%  Expected Input Format:
%    cols = [P1s(Pparam2=1),  P1s(@Param2), .... ]
%    row1 = Param1 values  (e.g. wavelength)
%    row2 = Param2 values  (e.g. width)
%    rows = Modes   
%
%    ** Rows SHOULD be SORTED so that INDEX is DECREASING top to bottom ***
%         however, an option is given to enforce this sorting 
%         the purpose of this option is to catch small things
%         that the use might have missed from the input data.
%
%
%    
 


function Sort_Output_of_ComsolSweep_by_Overlap()
    close all;
    clear;
    addpath 'G:\My Drive\Analyses(BF)\Matlab\UVQ';
    program_name = 'sort_comsol_output'
    
%% KEY INPUT PARAMETERS 
    Nfund               = 2;           %divider: designate these first (Nfund) rows as "Line1" (e.g. "Fundamental")
                                        %(remaining  rows will be "Line2" and checked against each of these )
    P1_limits           = [220,750];   % if empty just use the domain given by P1 values that are input
                                        % if not empty, extrapolate to these values for Line1
  
%% get data
    dirs.f = 'G:\My Drive\Analyses(BF)\Matlab\UVQ';
    dirs.d = 'G:\My Drive\Analyses(BF)\Matlab\matDat\Comsol_Output__Mode_Indices_from_ParamSweeps (INPUT DATA)\';
    
    cd(dirs.d); [FileName,PathName] = uigetfile({'*.xlsx'},program_name);
    filenametitle = replace(FileName,'_',' '); 
    
%% read data    
    %rows = different mode indices 
       %assume HIGHEST to LOWEST mode orders
    %columns = different parameters (wavelength, geometry, etc)
    A.neff = xlsread([PathName,FileName],'n_eff','c2:ff50');   %no more than 50 modes
    A.Ex = xlsread([PathName,FileName],'Ex-Ey','c2:ff50');   %no more than 50 modes
    A.ovlp = xlsread([PathName,FileName],'overlap','c2:ff50');   %no more than 50 modes
%     A = xlsread([PathName,FileName],'real_part_only','a2:ff50');   %no more than 50 modes
    % complex values should not be read (appear as "NaN"
    % first two rows are parameter values (for each column)

%% pull out the parameter values
    P1s = A.neff(1,:);
    P2s = A.neff(2,:);
    P1s_u = unique(P1s(~isnan(P1s)));    Nvals1 = length(P1s_u);    %Wavelength values
    P2s_u = unique(P2s(~isnan(P2s)));    Nvals2 = length(P2s_u);    %number of param1 values
    Nm = size(A,1)-2;                                               %number of modes    
    
    
 
%% Stack to 3D Cube on Parameter2   + SORT EACH ROW & COL: HIGH TO LOW
    B.neff = zeros(size(A.neff,1)-2, Nvals1, Nvals2);
    B.Ex = zeros(size(A.Ex,1)-2, Nvals1, Nvals2);
    B.ovlp = zeros(size(A.ovlp,1)-2, Nvals1, Nvals2);
    for k = 1:Nvals2
        p2 = P2s_u(k);
        icols_sel = find(P2s==p2);
        Nsel = length(icols_sel);
        B.neff(:,:,k) =  A.neff(3:end,icols_sel);
        B.Ex(:,:,k)   =  A.Ex(3:end,icols_sel);
        B.ovlp(:,:,k) =  A.ovlp(3:end,icols_sel);
    end
        
%% Divide Data into Fundamental ("Line 1") and SHG ("Line2") 
    Nfund;  %(input above) ... designate these two rows as "Fundamental"
    C.neff = B.neff(Nfund+1:end,:,:);   %these  rows comprise the "SHG Lines" and will assessed separately from fundamental
    B.neff = B.neff(1:Nfund,:,:);
    C.Ex = B.Ex(Nfund+1:end,:,:);   %these  rows comprise the "SHG Lines"  and will assessed separately
    B.Ex = B.Ex(1:Nfund,:,:);
    C.ovlp = B.ovlp(Nfund+1:end,:,:);   %these  rows comprise the "SHG Lines"and will sorted  
    B.ovlp = B.ovlp(1:Nfund,:,:);    
    
    
 %% SORT 
 % loop over columns
 % for each column, shuffle the (vertical) order of neighbor column to match that of the current column
 
 D.neff = C.neff(:,:,1);
 D.Ex = C.Ex(:,:,1);
 D.ovlp = C.ovlp(:,:,1);
    
 %work from Left to Right  or from Left to Right
[val,imax] = max(sum(~isnan(D.neff),1));  %imin = which column (first or last) has most real values
if imax~=1
    flag_reverse = 1;
    neff = fliplr(D.neff);
    TBL1 = fliplr(D.Ex);
    TBL2 = fliplr(D.ovlp);
else
    flag_reverse = 0;
    neff = D.neff;
    TBL1 = D.Ex;
    TBL2 = D.ovlp;
end


%Loop:  work Left to Right through table
TBL1_new = zeros(size(TBL1));  TBL1_new(:,1) = TBL1(:,1);
TBL2_new = zeros(size(TBL2));  TBL2_new(:,1) = TBL2(:,1);
indcs_jj = zeros(size(TBL1));  indcs_jj(:,1) = [1:size(TBL1,1)]';
neff_new = zeros(size(TBL1));  neff_new(:,1) = neff(:,1);

for kcol = 1:size(neff,2)-1
    isel_v = find( ~isnan(neff(:,kcol)) );      % identify all elements in target vector that are not NAN
    isel_u = find( ~isnan(neff(:,kcol+1)) );    % identify all elements in target vector that are not NAN
    v1 = TBL1_new(isel_v,kcol);                   %col vector that is target ('master') -- this column was shuffled on last pass of loop     
    u1 = TBL1(isel_u,kcol+1);                     %col vector that will be modified   
    v2 = TBL2_new(isel_v,kcol);                   %col vector that is target ('master') -- this column was shuffled on last pass of loop     
    u2 = TBL2(isel_u,kcol+1);                     %col vector that will be modified   
    
%     [u_s, j] = sort_U_to_match_V(v,u);          %execute sorting on vector next to this one....
%         %sorting concept:  align each element of next vector to "master" vector
        
    [u1_s, u2_s, j] = sort_U1u2_to_match_V1V2(v1,u1, v2,u2)  ; %execute sorting on vector next to this one....
        %sorting concept:  align each element of "next vector" to pair of "master" vectors
    
    
    TBL1_new(isel_u,kcol+1) = u1_s;
    TBL2_new(isel_u,kcol+1) = u2_s;
    indcs_jj(isel_u,kcol+1) = j;
    neff_new(isel_u,kcol+1) = neff(j, kcol+1);  %apply shuffling to each column of neff     
end

if flag_reverse==1
    TBL_new = fliplr(TBL_new);
    indcs_kk = fliplr(indcs_jj);
    neff_new = fliplr(neff_new);
end

 
disp(D.neff)
disp(neff_new)

%check families of modes -- n_effective v. parameter
figure; plot(P1s_u,D.neff');
figure; plot(P1s_u,neff_new');


 

        
%% Show complete family of input data curves ( for each P2 )
figctr=1;
for kp=1:length(P2_sel_vals)
    hf0 = figure;    figctr=figctr+1; 
    plot(P1s_u, B(1:Nfund,:,kp),'.-','Linewidth',3);            %fundamental Lines 
    hold on;   plot(P1s_FACTOR*P1s_u, C(:,:,kp),'o-');  xlabel('Param1');   %SHG Lines
    title(['Input: Family of Lines for P2=', num2str(P2_sel_vals(kp))] );                 
end

    
    
end %END --  MAIN FUNCTION



function [u_s, j] = sort_U_to_match_V(v,u)
    j = [];
    for k=1:length(v)                                               %loop over every element of v(k)
        i_remain_u = setxor(j,[1:length(u)]);                       %indices of u that have not yet been used
        if length(i_remain_u)>0
            [val,ii_sel] = min( abs( (u(i_remain_u)-v(k)) /v(k))  );    %find element of u that matches v(k) most closely
            imin = i_remain_u(ii_sel);                                  %record which element of u
            j(k,1) = imi n;                                              %index to u that matches v(k)
            u_s(k,1) = u(imin);                                         %record value of u that is the best match to v(k)
        end                                                                    %this is the sorted vector
    end
    
% figure; 
% subplot(2,1,1); semilogy([1:length(v)],v,'.-'); title('sort elements of vector "u" to match vector "v" ')
% hold all;       semilogy([1:length(u)],u,'.-'); legend('v = targ','u = to be sorted');
% subplot(2,1,2); semilogy([1:length(v)],v,'.-');  xlabel('element #')
% hold all;       semilogy([1:length(u_s)],u_s,'.-'); legend('v = targ',' **sorted** ');    
    
end

function [u1_s, u2_s, j] = sort_U1u2_to_match_V1V2(v1,u1, v2,u2)
    if length(v1)~=length(v2)  warning('input vector lengths do not match! (v1  v2) - master vectors '); end
    if length(u1)~=length(u2)  warning('input vector lengths do not match! (u1  u2) - to be shuffled '); end
    
    j = [];
    for k=1:length(v1)                                                %loop over every element of v(k)
        i_remain_u = setxor(j,[1:length(u1)]);                         %indices of u that have not yet been used
        if length(i_remain_u)>0
                  %find element(row) of [u1,u2] that matches [v1(k),v2(k)] most closely
            [val,ii_sel] = min( abs( (u1(i_remain_u)-v1(k)) / v1(k) )  ...
                             +  abs( (u2(i_remain_u)-v2(k)) / v2(k) )  );  %match u1->v1  and   u2 --> v2  (weighted equally)
            
            
            imin = i_remain_u(ii_sel);                                 %record which element of u
            j(k,1) = imin;                                             %index to [u1,u2] that matches [v1(k),v2(k)]
            u1_s(k,1) = u1(imin);                                      %record value of u1 that best matches
            u2_s(k,1) = u2(imin);                                      %record value of u2 that best matches
            
            
        end                                                            %this is the sorted vector
    end
end












% 
% 
% function MM = sort_v(MM,sort_direction)
%     for ll=1:size(MM,3)
%         M = MM(:,:,ll); 
%         for cc=1:size(M,2)
%             isel=~isnan(M(:,cc,:));
%             M(isel,cc) = sort(M(isel,cc,:),1,sort_direction); % !! 1 !!
%         end
%         MM(:,:,ll) = M;
%     end
% end
% function MM = sort_h(MM,sort_direction)
%     for ll=1:size(MM,3)
%         M = MM(:,:,ll); 
%         for rr=1:size(M,1)
%             isel=~isnan(M(rr,:,:));
%             M(rr,isel) = sort(M(rr,isel,:),2, sort_direction); % !! 2 !!
%         end
%         MM(:,:,ll) = M;
%     end
% end

  