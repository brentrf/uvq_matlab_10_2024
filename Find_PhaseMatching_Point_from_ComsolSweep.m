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
% Key Parameters to Update (in hard code):
%  a.)   Nfund = 2 ... how many of top rows to assign as "fundamental"
%                  ... The fundamental modes will be assessed for crossings 
%                  at their given wavelength.  ( Keep this below 2 or 3, else # plots will explode)
%                  ... All other modes (rows below Nfund) will have
%                  wavelength doubled (2x) in order to look for crossing
%                  with fundamentals.
%  b.) dosort_v = 1 ... This will enforce each column to be in *descending* order
%      dosort_h = 1 ... This will enforce each ROW to be in *descending* order (performed second)
%  c.) P2_sel_vals_extra   = [ ] 
%                    .... this lets you find intersections of lines for 
%                         values of P2 that are between those that are given  
%                    ....  leave this empty if you don't want to check line intersection for additional P2 values

%  c2.) P1_limits   =  []  if empty just use the domain given by P1 values that are input
%                         if not empty, extrapolate to these values for Line1
                                    
%  d.) P1s_FACTOR = (1,2) ... if this is 2, then the x-value of second line
%                          will be doubled (multipled by FACTOR) 
%                        (if P1 = Wavelength, then set P1s_FACTOR=2,  else 1 )
%
%  e.) dWL = ____   ... precision of WL in the identifying crossing
%  f.) methodstr = 'linear';  ... interpolation method 
%  g.) 'extrap'    ... whether or not to exp
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
    Nfund               = 2;           %divider: designate these first (Nfund) rows as "Line1" (e.g. "Fundamental")
                                        %(remaining  rows will be "Line2" and checked against each of these )
    Nfund_to_test       = 2;            %How many of the "Line1" (e.g. "Fundamental") to check against the remaining lines
                                        %(remaining  rows will be "Line2" and checked against each of these )
    P1_limits           = [220,750];   % if empty just use the domain given by P1 values that are input
                                        % if not empty, extrapolate to these values for Line1
    dosort_v              = 0;        %this enforces sorting of effective index from HIGH to LOW along columns & rows of input
    dosort_h              = 0;        %this enforces sorting of effective index from HIGH to LOW along columns & rows of input
    sort_direction        = 'descend';
    P1s_FACTOR            = 2;        %if this is 2, then the x-value of second line will be doubled (multipled by FACTOR)
                                    % (if P1 = Wavelength,  then set P1s_FACTOR=2,  else set = 1 )
                                    
    P2_sel_vals_extra   = [];           %leave this empty if you don't want to run additional values it....
    
    methodstr           = 'pchip'; %'spline'; %'linear'  'makima' 'pchip'
   
    
    
%     Nfund               = 15;           %divider: designate these first (Nfund) rows as "Line1" (e.g. "Fundamental")
%                                         %(remaining  rows will be "Line2" and checked against each of these )
%     Nfund_to_test       = 2;            %How many of the "Line1" (e.g. "Fundamental") to check against the remaining lines
%                                         %(remaining  rows will be "Line2" and checked against each of these )
%     P1_limits           = []; %[200,1100];   % if empty just use the domain given by P1 values that are input
%                                         % if not empty, extrapolate to these values for Line1
%     dosort_v              = 0;        %this enforces sorting of effective index from HIGH to LOW along columns & rows of input
%     dosort_h              = 0;        %this enforces sorting of effective index from HIGH to LOW along columns & rows of input
%     sort_direction        = 'ascend';  %'descend';
%     P1s_FACTOR            = 1;        %if this is 2, then the x-value of second line will be doubled (multipled by FACTOR)
%                                     % (if P1 = Wavelength,  then set P1s_FACTOR=2,  else set = 1 )
%                                     
%     P2_sel_vals_extra   = [];           %leave this empty if you don't want to run additional values it....
% 


%% get data
    dirs.f = 'G:\My Drive\Analyses(BF)\Matlab\UVQ';
    dirs.d = 'G:\My Drive\Analyses(BF)\Matlab\matDat\Comsol_Output__Mode_Indices_from_ParamSweeps (INPUT DATA)\';
    
    cd(dirs.d); [FileName,PathName] = uigetfile({'*.xlsx'},program_name);
    filenametitle = replace(FileName,'_',' '); 
    
%% read data    
    %rows = different mode indices 
       %assume HIGHEST to LOWEST mode orders
    %columns = different parameters (wavelength, geometry, etc)
    A = xlsread([PathName,FileName],'n_eff','c2:ff50');   %no more than 50 modes
%     A = xlsread([PathName,FileName],'real_part_only','a2:ff50');   %no more than 50 modes
    % complex values should not be read (appear as "NaN"
    % first two rows are parameter values (for each column)

%% pull out the parameter values
    P1s = A(1,:);
    P2s = A(2,:);
    P1s_u = unique(P1s(~isnan(P1s)));    Nvals1 = length(P1s_u);    %Wavelength values
    P2s_u = unique(P2s(~isnan(P2s)));    Nvals2 = length(P2s_u);    %number of param1 values
    Nm = size(A,1)-2;                                               %number of modes    
    
%% Stack to 3D Cube on Parameter2   + SORT EACH ROW & COL: HIGH TO LOW
    B = zeros(size(A,1)-2, Nvals1, Nvals2);
    for k = 1:Nvals2
        p2 = P2s_u(k);
        icols_sel = find(P2s==p2);
        Nsel = length(icols_sel);
        B(:,:,k) =  A(3:end,icols_sel);
    end
        
%% Divide Data into Fundamental ("Line 1") and SHG ("Line2") 
    Nfund;  %(input above) ... designate these two rows as "Fundamental"
    C = B(Nfund+1:end,:,:);   %these  rows comprise the "SHG Lines" and will be compared to the "Fundamental lines" for intersections
    B = B(1:Nfund,:,:);
    
    % show input
    figure; plot(P1s_u,C(:,:,1)');  xlabel('P1s_u'); title('Layer = 1');
    figure; plot(P1s_u,C(:,:,2)');  xlabel('P1s_u'); title('Layer = 2');   
    
%% Sort If Desired
    if dosort_v==1              %sort_VERTICAL
        B = sort_v(B, sort_direction);
        C = sort_v(C, sort_direction);
    end
    if dosort_h==1              %sort_HORIZONTAL
        B = sort_h(B, sort_direction);
        C = sort_h(C, sort_direction);        
    end

    
    
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
        
%% Show complete family of input data curves ( for each P2 )
figctr=1;
for kp=1:length(P2_sel_vals)
    hf0 = figure;    figctr=figctr+1; 
    plot(P1s_u, B(1:Nfund,:,kp),'.-','Linewidth',3);            %fundamental Lines 
    hold on;   plot(P1s_FACTOR*P1s_u, C(:,:,kp),'o-');  xlabel('Param1');   %SHG Lines
    title(['Input: Family of Lines for P2=', num2str(P2_sel_vals(kp))] );                 
end

    
%%   *** Find Crossing along WL(param1) at each DESIRED Param2 Value ***
    kf=1;
    kp=1;    
    % ** Find Crossing along  WL(Param1) at each USER DEFINED Param2  **
    for kf=1:Nfund_to_test                      %LOOP:  fundamental modes
        for kp=1:length(P2_sel_vals)            %LOOP:  over Param2
    
            legnames = [];
            
            P2value = P2_sel_vals(kp);
            
            
            %USE 2D-Interpolation to find Line between two given P2 points
            if ismember(P2value,P2s_u)
                %No P2 interpolation ... 
                n_fund_points_tmp = B(kf,:,find(P2value==P2s_u));
            else
                %Interp over P2:  get Index v WL points at desired P2 Value by 2D Interpolation --- Vq = interp2(X,Y,V,Xq,Yq);    
                tmp = reshape(  B(kf,:,:) , Nvals1, Nvals2);  %rows = P1s, cols=Param2
                n_fund_points_tmp =  interp2(P2s_u,   P1s_u,   tmp,  P2value, P1s_u); % OR: griddata(P2s_u,   P1s_u,   tmp,  P2value, P1s_u);            
            end
            
             %Interp over P1 (P1s) with Extrapolation 
            if length(P1_limits)==2
                P1s  = [min(P1_limits), P1s_u, max(P1_limits) ]
                isel_isan = find(~isnan(n_fund_points_tmp ));
                n_fund = interp1(P1s_u(isel_isan), n_fund_points_tmp(isel_isan), P1s, methodstr,'extrap') ;
            else disp('NO Extrapolation')
                P1s = P1s_u; n_fund = n_fund_points_tmp;
            end
            

            %at this Param2 Value: find P1 crossing with each NON-FUND. Mode (using 1D interpolation)
            for kc=1:size(C,1)                  %LOOP: non-fundamental modes
                
            
                %USE 2D-Interpolation to find Line between two given P2 points
                if ismember(P2value,P2s_u)
                    %No P2 interpolation 
                    n_shg_points_tmp = C(kc,:,find(P2value==P2s_u));     
                else
                    %Interp over P2:
                    tmp = reshape(  C(kc,:,:) , Nvals1, Nvals2);  %rows = P1s, cols=Param2
                    n_shg_points_tmp = interp2(P2s_u,   P1s_u,   tmp,  P2value, P1s_u);
                    %  n_shg_points_tmp = griddata(P2s_u,   P1s_u,   tmp,  P2value, P1s_u); % (alternative to interp2)
                end
                
                
                % ** INTERSECTION FINDER **
                %Interp over P1 (P1s) with Extrapolation  & Find Crossing ...
                %                 [found_crossing,cross_interp_WL,cross_interp_n, n_shg] = find_crossing(P1s,n_fund, P1s_u, n_shg_points_tmp, methodstr);
                
                x1 = P1s;
                y1 = n_fund;
                x2 = P1s_u * P1s_FACTOR;
                y2 = n_shg_points_tmp;
                robust = 0;
                [x0,y0,iout,jout] = intersections(x1,y1,x2,y2,robust);      %CALL KEY FUNCTION 
                %calculate slope diff at intersection
                slopes1 = ( y1(ceil(iout)) - y1(floor(iout)) ) / ...
                          ( x1(ceil(iout)) - x1(floor(iout)) );
                slopes2 = ( y2(ceil(jout)) - y2(floor(jout)) ) / ...
                          ( x2(ceil(jout)) - x2(floor(jout)) );
                if length(slopes1)>0   Dslopes = slopes2 - slopes1; else Dslopes = -1; end
                
                %record crossing results
                found_crossing = length(x0)>0;
                if length(x0)>1  warning('more than 1 intersection found...'); 
                    cross_interp_WL = x0(1);
                    cross_interp_n = y0(1);
                else
                    cross_interp_WL = x0;
                    cross_interp_n = y0;
                end                                                                                 
                cross_lam_interp(kc,kp,kf) =  max([-1,cross_interp_WL]);   % [ row=SHGmode;  ool=Param2;  layer=FundMode ]
                cross_n_interp(kc,kp,kf)   = max([-1,cross_interp_n]);
                cross_Dslopes(kc,kp,kf)    = Dslopes ;
                
                %FIGURE - for checking results
                if kc==1                     
                    hff(kf,kp) = figure(figctr); figctr=figctr+1;
                    plot(x1,y1,'o-'); legnames{1}='fund';
                end
                if found_crossing==1
                    figure(hff(kf,kp));  hold all; plot(x2,y2,'o-'); % hold on;  plot(x2,y2,'o');
                    Nlegnames = length(legnames); legnames{Nlegnames+1}=num2str(kc);  % legnames{Nlegnames+2}=num2str(kc);                       
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
        
        
        disp(['Wavelengths of Crossings of kf=',num2str(kf),'  (1 column for each P2 values):'])
        disp(cross_lam_interp(:,:,kf));
        disp(['Index of Crossings of kf=',num2str(kf),'  (1 column for each P2 values):'])
        disp(cross_n_interp(:,:,kf));
        disp(['Relative Slope Difference of Crossings of kf=',num2str(kf),'  (1 column for each P2 values):'])        
        cross_Dslopes(:,:,kf)
        disp('%(-1): only one point, cant interpolate  n v WL')
        disp('%(-2): No Intersection within domain')
        disp('%(-3): No Sign Change with Crossing')
        
    end         %loop fundamental modes
    
    
end %END --  MAIN FUNCTION

function MM = sort_v(MM,sort_direction)
    for ll=1:size(MM,3)
        M = MM(:,:,ll); 
        for cc=1:size(M,2)
            isel=~isnan(M(:,cc,:));
            M(isel,cc) = sort(M(isel,cc,:),1,sort_direction); % !! 1 !!
        end
        MM(:,:,ll) = M;
    end
end
function MM = sort_h(MM,sort_direction)
    for ll=1:size(MM,3)
        M = MM(:,:,ll); 
        for rr=1:size(M,1)
            isel=~isnan(M(rr,:,:));
            M(rr,isel) = sort(M(rr,isel,:),2, sort_direction); % !! 2 !!
        end
        MM(:,:,ll) = M;
    end
end



function [found_crossing,cross_interp_WL,cross_interp_n, n_shg] = find_crossing(P1s,n_fund, P1s_u,C_indices, methodstr)
                n_shg=[];
                found_crossing=0;  %reset      
                isel_isanC = find(~isnan(C_indices));
                
                %INTERPOLATION INTERSECTION
                if length(isel_isanC)<=1
                    cross_interp_WL = -1;                      %(-1) only one point, can't interpolate  n v WL
                    cross_interp_n  = -1;
                else
                    n_shg  = interp1(2*P1s_u(isel_isanC), C_indices(isel_isanC), P1s, methodstr, 'extrap') ;
                    [mindiff, imin]  = min(abs(n_shg - n_fund));
                    if ~(1 < imin & imin < length(n_shg)) %imin==1 | imin==length(n_shg)        
                        cross_interp_WL = -2;                  %(-2): No Intersection within domain
                        cross_interp_n  = -2;
                    else
                        if sign(n_shg(imin-1)-n_fund(imin-1)) ~= sign(n_shg(imin+1)-n_fund(imin+1))
                            cross_interp_WL = P1s(imin);
                            cross_interp_n = n_shg(imin);
                            found_crossing = 1;
                        else
                            cross_interp_WL = -3;              %(-3): no Crossing Detected
                            cross_interp_n = -3;
                        end
                    end
                end  
                
end
    
    

%% INTERSECTIONS


function [x0,y0,iout,jout] = intersections(x1,y1,x2,y2,robust)
%INTERSECTIONS Intersections of curves.
%   Computes the (x,y) locations where two curves intersect.  The curves
%   can be broken with NaNs or have vertical segments.
%
% Example:
%   [X0,Y0] = intersections(X1,Y1,X2,Y2,ROBUST);
%
% where X1 and Y1 are equal-length vectors of at least two points and represent
% curve 1.  Similarly, X2 and Y2 represent curve 2.  X0 and Y0 are column
% vectors containing the points at which the two curves intersect.
%
% ROBUST (optional) set to 1 or true means to use a slight variation of the
% algorithm that might return duplicates of some intersection points, and then
% remove those duplicates.  The default is true, but since the algorithm is
% slightly slower you can set it to false if you know that your curves don't
% intersect at any segment boundaries.  Also, the robust version properly
% handles parallel and overlapping segments.
%
% The algorithm can return two additional vectors that indicate which segment
% pairs contain intersections and where they are:
%
%   [X0,Y0,I,J] = intersections(X1,Y1,X2,Y2,ROBUST);
%
% For each element of the vector I, I(k) = (segment number of (X1,Y1)) + (how
% far along this segment the intersection is).  For example, if I(k) = 45.25
% then the intersection lies a quarter of the way between the line segment
% connecting (X1(45),Y1(45)) and (X1(46),Y1(46)).  Similarly for the vector J
% and the segments in (X2,Y2).
%
% You can also get intersections of a curve with itself.  Simply pass in only
% one curve, i.e.,
%
%   [X0,Y0] = intersections(X1,Y1,ROBUST);
%
% where, as before, ROBUST is optional.

% Version: 3.0, 30 May 2024
% Author:  Douglas M. Schwarz
% Email:   dmschwarz=ieee*org, dmschwarz=urgrad*rochester*edu
% Real_email = regexprep(Email,{'=','*'},{'@','.'})


% Theory of operation:
%
% Given two line segments, L1 and L2,
%
%   L1 endpoints:  (x1(1),y1(1)) and (x1(2),y1(2))
%   L2 endpoints:  (x2(1),y2(1)) and (x2(2),y2(2))
%
% we can write four equations with four unknowns and then solve them.  The four
% unknowns are t1, t2, x0 and y0, where (x0,y0) is the intersection of L1 and
% L2, t1 is the distance from the starting point of L1 to the intersection
% relative to the length of L1 and t2 is the distance from the starting point of
% L2 to the intersection relative to the length of L2.
%
% So, the four equations are
%
%    (x1(2) - x1(1))*t1 = x0 - x1(1)
%    (x2(2) - x2(1))*t2 = x0 - x2(1)
%    (y1(2) - y1(1))*t1 = y0 - y1(1)
%    (y2(2) - y2(1))*t2 = y0 - y2(1)
%
% Rearranging and writing in matrix form,
%
%  [x1(2)-x1(1)       0       -1   0;      [t1;      [-x1(1);
%        0       x2(2)-x2(1)  -1   0;   *   t2;   =   -x2(1);
%   y1(2)-y1(1)       0        0  -1;       x0;       -y1(1);
%        0       y2(2)-y2(1)   0  -1]       y0]       -y2(1)]
%
% Let's call that A*T = B.  We can solve for T with T = A\B.
%
% Once we have our solution we just have to look at t1 and t2 to determine
% whether L1 and L2 intersect.  If 0 <= t1 < 1 and 0 <= t2 < 1 then the two line
% segments cross and we can include (x0,y0) in the output.
%
% In principle, we have to perform this computation on every pair of line
% segments in the input data.  This can be quite a large number of pairs so we
% will reduce it by doing a simple preliminary check to eliminate line segment
% pairs that could not possibly cross.  The check is to look at the smallest
% enclosing rectangles (with sides parallel to the axes) for each line segment
% pair and see if they overlap.  If they do then we have to compute t1 and t2
% (via the A\B computation) to see if the line segments cross, but if they don't
% then the line segments cannot cross.  In a typical application, this technique
% will eliminate most of the potential line segment pairs.


% Input checks.
if verLessThan('matlab','7.13')
	error(nargchk(2,5,nargin)) %#ok<NCHKN>
else
	narginchk(2,5)
end

% Adjustments based on number of arguments.
switch nargin
	case 2
		robust = true;
		x2 = x1;
		y2 = y1;
		self_intersect = true;
	case 3
		robust = x2;
		x2 = x1;
		y2 = y1;
		self_intersect = true;
	case 4
		robust = true;
		self_intersect = false;
	case 5
		self_intersect = false;
end

% x1 and y1 must be vectors with same number of points (at least 2).
if sum(size(x1) > 1) ~= 1 || sum(size(y1) > 1) ~= 1 || ...
		length(x1) ~= length(y1)
	error('X1 and Y1 must be equal-length vectors of at least 2 points.')
end
% x2 and y2 must be vectors with same number of points (at least 2).
if sum(size(x2) > 1) ~= 1 || sum(size(y2) > 1) ~= 1 || ...
		length(x2) ~= length(y2)
	error('X2 and Y2 must be equal-length vectors of at least 2 points.')
end


% Force all inputs to be column vectors.
x1 = x1(:);
y1 = y1(:);
x2 = x2(:);
y2 = y2(:);

% Compute number of line segments in each curve and some differences we'll need
% later.
n1 = length(x1) - 1;
n2 = length(x2) - 1;
xy1 = [x1 y1];
xy2 = [x2 y2];
dxy1 = diff(xy1);
dxy2 = diff(xy2);


% Determine the combinations of i and j where the rectangle enclosing the i'th
% line segment of curve 1 overlaps with the rectangle enclosing the j'th line
% segment of curve 2.

% Original method that works in old MATLAB versions, but is slower than using
% binary singleton expansion (explicit or implicit).
% [i,j] = find( ...
% 	repmat(mvmin(x1),1,n2) <= repmat(mvmax(x2).',n1,1) & ...
% 	repmat(mvmax(x1),1,n2) >= repmat(mvmin(x2).',n1,1) & ...
% 	repmat(mvmin(y1),1,n2) <= repmat(mvmax(y2).',n1,1) & ...
% 	repmat(mvmax(y1),1,n2) >= repmat(mvmin(y2).',n1,1));

% Select an algorithm based on MATLAB version and number of line segments in
% each curve.  We want to avoid forming large matrices for large numbers of line
% segments.  If the matrices are not too large, choose the best method available
% for the MATLAB version.
if n1 > 1000 || n2 > 1000 || verLessThan('matlab','7.4')
	% Determine which curve has the most line segments.
	if n1 >= n2
		% Curve 1 has more segments, loop over segments of curve 2.
		ijc = cell(1,n2);
		min_x1 = mvmin(x1);
		max_x1 = mvmax(x1);
		min_y1 = mvmin(y1);
		max_y1 = mvmax(y1);
		for k = 1:n2
			k1 = k + 1;
			ijc{k} = find( ...
				min_x1 <= max(x2(k),x2(k1)) & max_x1 >= min(x2(k),x2(k1)) & ...
				min_y1 <= max(y2(k),y2(k1)) & max_y1 >= min(y2(k),y2(k1)));
			ijc{k}(:,2) = k;
		end
		ij = vertcat(ijc{:});
		i = ij(:,1);
		j = ij(:,2);
	else
		% Curve 2 has more segments, loop over segments of curve 1.
		ijc = cell(1,n1);
		min_x2 = mvmin(x2);
		max_x2 = mvmax(x2);
		min_y2 = mvmin(y2);
		max_y2 = mvmax(y2);
		for k = 1:n1
			k1 = k + 1;
			ijc{k}(:,2) = find( ...
				min_x2 <= max(x1(k),x1(k1)) & max_x2 >= min(x1(k),x1(k1)) & ...
				min_y2 <= max(y1(k),y1(k1)) & max_y2 >= min(y1(k),y1(k1)));
			ijc{k}(:,1) = k;
		end
		ij = vertcat(ijc{:});
		i = ij(:,1);
		j = ij(:,2);
	end
	
elseif verLessThan('matlab','9.1')
	% Use bsxfun.
	[i,j] = find( ...
		bsxfun(@le,mvmin(x1),mvmax(x2).') & ...
		bsxfun(@ge,mvmax(x1),mvmin(x2).') & ...
		bsxfun(@le,mvmin(y1),mvmax(y2).') & ...
		bsxfun(@ge,mvmax(y1),mvmin(y2).'));
	
else
	% Use implicit expansion.
	[i,j] = find( ...
		mvmin(x1) <= mvmax(x2).' & mvmax(x1) >= mvmin(x2).' & ...
		mvmin(y1) <= mvmax(y2).' & mvmax(y1) >= mvmin(y2).');
	
end


% Find segments pairs which have at least one vertex = NaN and remove them.
% This line is a fast way of finding such segment pairs.  We take advantage of
% the fact that NaNs propagate through calculations, in particular subtraction
% (in the calculation of dxy1 and dxy2, which we need anyway) and addition.  At
% the same time we can remove redundant combinations of i and j in the case of
% finding intersections of a line with itself.
if self_intersect
	remove = isnan(sum(dxy1(i,:) + dxy2(j,:),2)) | j <= i + 1;
else
	remove = isnan(sum(dxy1(i,:) + dxy2(j,:),2));
end
i(remove) = [];
j(remove) = [];

% Initialize matrices.  We'll put the T's and B's in matrices and use them one
% column at a time.  AA is a 3-D extension of A where we'll use one plane at a
% time.
n = length(i);
T = zeros(4,n);
AA = zeros(4,4,n);
AA([1 2],3,:) = -1;
AA([3 4],4,:) = -1;
AA([1 3],1,:) = dxy1(i,:).';
AA([2 4],2,:) = dxy2(j,:).';
B = -[x1(i) x2(j) y1(i) y2(j)].';

% Loop through possibilities.  If robust is true (any number > 0), trap
% singularity warning and then use lastwarn to see if that plane of AA is near
% singular.  Process any such segment pairs to determine if they are colinear
% (overlap) or merely parallel.  That test consists of checking to see if one of
% the endpoints of the curve 2 segment lies on the curve 1 segment.  This is
% done by checking the cross product
%
%   (x1(2),y1(2)) - (x1(1),y1(1)) x (x2(2),y2(2)) - (x1(1),y1(1)).
%
% If this is close to zero then the segments overlap.
%
% If robust is < 0, do the above using the older, slower algorithm.
%
% If the robust option is false (0) then we assume no two segment pairs are
% parallel and just go ahead and do the computation.  If A is ever singular, a
% warning will appear.  This is faster and obviously you should use it only when
% you know you will never have overlapping or parallel segment pairs.

if robust > 0
	overlap = false(n,1);
	warning_state = warning('off','MATLAB:singularMatrix');
	% Use try-catch to guarantee original warning state is restored.
	try
		lastwarn('')
		% Attempt to use pagemldivide (new in R2022a) with a fallback to the old
		% algorithm if it doesn't exist.
		try
			T = permute(pagemldivide(AA,permute(B,[1 3 2])),[1 3 2]);
		catch
			for k = 1:n
				[L,U] = lu(AA(:,:,k));
				T(:,k) = U\(L\B(:,k));
			end
		end
		% Look more carefully at intersections where any element of each column
		% of T is NaN.
		suspect_list = find(any(isnan(T),1));
		for k = suspect_list
			T(:,k) = AA(:,:,k)\B(:,k);
			[unused,last_warn] = lastwarn; %#ok<ASGLU>
			lastwarn('')
			if strcmp(last_warn,'MATLAB:singularMatrix')
				% Force in_range(k) to be false.
				T(1,k) = NaN;
				% Determine if these segments overlap or are just parallel.
				overlap(k) = rcond([dxy1(i(k),:);xy2(j(k),:) - xy1(i(k),:)]) < eps;
			end
		end
		warning(warning_state)
	catch err
		warning(warning_state)
		rethrow(err)
	end
	% Find where t1 and t2 are between 0 and 1 and return the corresponding
	% x0 and y0 values.
	in_range = (T(1,:) >= 0 & T(2,:) >= 0 & T(1,:) <= 1 & T(2,:) <= 1).';
	% For overlapping segment pairs the algorithm will return an intersection
	% point that is at the center of the overlapping region.
	if any(overlap)
		ia = i(overlap);
		ja = j(overlap);
		% set x0 and y0 to middle of overlapping region.
		T(3,overlap) = (max(min(x1(ia),x1(ia+1)),min(x2(ja),x2(ja+1))) + ...
			min(max(x1(ia),x1(ia+1)),max(x2(ja),x2(ja+1)))).'/2;
		T(4,overlap) = (max(min(y1(ia),y1(ia+1)),min(y2(ja),y2(ja+1))) + ...
			min(max(y1(ia),y1(ia+1)),max(y2(ja),y2(ja+1)))).'/2;
		selected = in_range | overlap;
	else
		selected = in_range;
	end
	xy0 = T(3:4,selected).';
	
	% Remove duplicate intersection points.
	[xy0,index] = unique(xy0,'rows');
	x0 = xy0(:,1);
	y0 = xy0(:,2);
	
	% Compute how far along each line segment the intersections are.
	if nargout > 2
		sel_index = find(selected);
		sel = sel_index(index);
		iout = i(sel) + T(1,sel).';
		jout = j(sel) + T(2,sel).';
	end

elseif robust < 0 % use the algorithm from version 2.0 of this m-file
	overlap = false(n,1);
	warning_state = warning('off','MATLAB:singularMatrix');
	% Use try-catch to guarantee original warning state is restored.
	try
		lastwarn('')
		for k = 1:n
			T(:,k) = AA(:,:,k)\B(:,k);
			[unused,last_warn] = lastwarn; %#ok<ASGLU>
			lastwarn('')
			if strcmp(last_warn,'MATLAB:singularMatrix')
				% Force in_range(k) to be false.
				T(1,k) = NaN;
				% Determine if these segments overlap or are just parallel.
				overlap(k) = rcond([dxy1(i(k),:);xy2(j(k),:) - xy1(i(k),:)]) < eps;
			end
		end
		warning(warning_state)
	catch err
		warning(warning_state)
		rethrow(err)
	end
	% Find where t1 and t2 are between 0 and 1 and return the corresponding
	% x0 and y0 values.
	in_range = (T(1,:) >= 0 & T(2,:) >= 0 & T(1,:) <= 1 & T(2,:) <= 1).';
	% For overlapping segment pairs the algorithm will return an intersection
	% point that is at the center of the overlapping region.
	if any(overlap)
		ia = i(overlap);
		ja = j(overlap);
		% set x0 and y0 to middle of overlapping region.
		T(3,overlap) = (max(min(x1(ia),x1(ia+1)),min(x2(ja),x2(ja+1))) + ...
			min(max(x1(ia),x1(ia+1)),max(x2(ja),x2(ja+1)))).'/2;
		T(4,overlap) = (max(min(y1(ia),y1(ia+1)),min(y2(ja),y2(ja+1))) + ...
			min(max(y1(ia),y1(ia+1)),max(y2(ja),y2(ja+1)))).'/2;
		selected = in_range | overlap;
	else
		selected = in_range;
	end
	xy0 = T(3:4,selected).';
	
	% Remove duplicate intersection points.
	[xy0,index] = unique(xy0,'rows');
	x0 = xy0(:,1);
	y0 = xy0(:,2);
	
	% Compute how far along each line segment the intersections are.
	if nargout > 2
		sel_index = find(selected);
		sel = sel_index(index);
		iout = i(sel) + T(1,sel).';
		jout = j(sel) + T(2,sel).';
	end

else % non-robust option
	% Attempt to use pagemldivide (new in R2022a) with a fallback to the old
	% algorithm if it doesn't exist.
	try
		T = permute(pagemldivide(AA,permute(B,[1 3 2])),[1 3 2]);
	catch
		for k = 1:n
			[L,U] = lu(AA(:,:,k));
			T(:,k) = U\(L\B(:,k));
		end
	end
	
	% Find where t1 and t2 are between 0 and 1 and return the corresponding
	% x0 and y0 values.
	in_range = (T(1,:) >= 0 & T(2,:) >= 0 & T(1,:) < 1 & T(2,:) < 1).';
	x0 = T(3,in_range).';
	y0 = T(4,in_range).';
	
	% Compute how far along each line segment the intersections are.
	if nargout > 2
		iout = i(in_range) + T(1,in_range).';
		jout = j(in_range) + T(2,in_range).';
	end
end

end %intersections function

% Plot the results (useful for debugging).
% plot(x1,y1,x2,y2,x0,y0,'ok');

function y = mvmin(x)
% Faster implementation of movmin(x,k) when k = 1.
y = min(x(1:end-1),x(2:end));
end

function y = mvmax(x)
% Faster implementation of movmax(x,k) when k = 1.
y = max(x(1:end-1),x(2:end));
end



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%  CODE BELOW IS NOT WORKING .... 
%%  POLY FIT -- to find crossings

% 
%     for kf=1:Nfund
%         for kp=1:Nvals2 %loop over Param2
%             isel_isan = find(~isnan(B(kf,:,kp)));
%             [fund.a, fund.S] = polyfit(P1s_u(isel_isan), B(kf,isel_isan,kp), 2);
%             
%             for kc=1:size(C,1)  %loop over non-fundamental modes
%                  isel_isanC = find(~isnan(C(kc,:,kp)));
%                 
%                 %POLYFIT: intersection wavelength (solve quadratic)
%                 if length(isel_isanC)>=2  
%                     if length(isel_isanC)>=3   %(quadratic)
%                         [shg.b, shg.S] = polyfit(2*P1s_u(isel_isanC), C(kc,isel_isanC,kp), 2);
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
%                 plot(P1s,polyval(fund.a,P1s)); hold all; plot(P1s,polyval(shg.b,P1s));
%                 title(['kf=',num2str(kf),'  kp=',num2str(kp),'  kc=',num2str(kc)]);
%                 hold off;
%                 %pause(0.1); 
% 
%                 %record crossing result
%                 if isreal(crossfit.lam1) & P1s(1) <= crossfit.lam1 & crossfit.lam1 <= P1s(end)  
%                     crosslam_fit1(kf,kp,kc) =  crossfit.lam1;
%                 else crosslam_fit1(kf,kp,kc) = -3;
%                 end
%                 if isreal(crossfit.lam2) & P1s(1) <= crossfit.lam2 & crossfit.lam2 <= P1s(end)  
%                     crosslam_fit2(kf,kp,kc) =  crossfit.lam2;
%                 else crosslam_fit2(kf,kp,kc) = -3;
%                 end
%                 
%                 disp(['kf=',num2str(kf),'  kp=',num2str(kp),'  kc=',num2str(kc)]);
%                 
%             end %loop: shg modes
%         end     %loop: parameter 2
%     end         %loop fundamental modes
