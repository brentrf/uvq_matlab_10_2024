%  SHG: Phase Matching Calculations
%
%
% *Uniaxial* 
% crystal axes define coordinate system 
%   axis3 = optical axis (Z)
%   axis2 = perpendicular plane (Y)
%   axis1 = perpendicular plane (X)
%   
%

% maxwell
% curl(H) = dD/dt   [D = eps*E]
% curl(E) = dB/dt   [B = mu*H]
% div(D) = 0    (source free) 
% div(B) = 0    (source free)
%
%

function Analyze_NLO()

    %% INPUTS

    Xtal_Tens = ...              % BBO
    [   0,	0,	7.10E-24;
        0,	1.40E-23,	0;
        7.10E-24,	0,	0;
        0,      0,      0;
        0,      0,      0;
        0,      0,      0];
  
    Xtal_Dispersion = ...       % DEFAULT = 0
    [200,   250,    300,    350,    400,    450,    500;...
     2.5,   2.5,    2.5,    2.5,    2.5,    2.5,    2.5;...  %ordinary
     2.5,   2.5,    2.5,    2.5,    2.5,    2.5,    2.5];    %extraordinary
        
    S0_in = [0,1,0];  %Poynting vector of incident ray (wavefront)
    E0_in = [1,0,0];  %polarization of input beam
    
    
   
    
    
    
    
    %[1] Vary Polarization Angle at **** Normal Incidence ****
    
    S0 = S0_in; %Poynting vector of incident ray (wavefront)
    E0 = E0_in; %polarization of input beam
    H0 = cross(S0,E0);
    theta_pol_list = linspace(0,360,100)';
    
    for k_pol = 1:length(theta_pol_list)
        theta_pol = theta_pol_list(k_pol);
        
        E0 = [cosd(theta_pol),sind(theta_pol),0];     %update polarization angle
        
        %  S0 = cross(H0,E0)                          %ensure poynting is perpenticular to polarization plane
        [d,d_mag_v_pol(k_pol)] = calc_NLO(S0,E0, Xtal_Tens);
        
    end
    %PLOT v. POLARIZATION
    figure; subplot(1,2,1); plot(theta_pol_list,d_mag_v_pol); 
    xlabel('polarization angle ( 0deg || axis1 )');    ylabel('d = nonlinear response')
    subplot(1,2,2);  polarplot(deg2rad(theta_pol_list),d_mag_v_pol);  grid on;
    title('Nonlinear Response Coefficient v. Polarisation at NORMAL INCIDENCE'); grid on;
    
    
    
    
    
    %[2] Vary Angle of Incidence at **** Constant Polarization ****
    
    S0ref = S0_in;   %Poynting vector of incident ray (wavefront)
    E0ref = E0_in;   %polarization of input beam
    H0ref = cross(S0ref,E0ref); 
    alpha_list = linspace(0,80,80)'; %ANGLES OF INCIDENCE
    phi_list = linspace(0,90,9)';   %phi of AOI
    phi_deg = 0;        %phi = alignment of rotation axis around which Angle of Incidence is formed  (phi=0 || axis1)    

    
    for kp = 1:length(phi_list)    
        phi_deg = phi_list(kp);      
        rotPHI  = [[1,0,0 ];[0, cosd(phi_deg), -sind(phi_deg)]; [0, sind(phi_deg), cosd(phi_deg)]];
        legnames{kp} = (['phi=',num2str(phi_deg),'deg']);
        for ka = 1:length(alpha_list)    
            alpha_deg = alpha_list(ka);
            rotALPH = [[cosd(alpha_deg), -sind(alpha_deg),0];[sind(alpha_deg), cosd(alpha_deg),0 ];[0,0,1]];
            %apply rotations (first: rotPHI, second: rotALPH)
            S1 = rotALPH*rotPHI*S0ref';
            E1 = rotALPH*rotPHI*E0ref';
            %[1] calculate Nonlinear Response
            [d,d_mag_v_AOI(ka,kp)] = calc_NLO(S1',E1', Xtal_Tens);
            
            %[2] calculate phase mismatch
            
            %[3] calculate conversion Efficiency
            

            
        end
    end    
    
    %PLOT d v. ANGLE OF INCIDENCE
    figure; plot(alpha_list,d_mag_v_AOI); legend(legnames);
    xlabel('angle of incidence (deg)'); ylabel('d = nonlinear response');
    title('Nonlinear Response Coefficient v. ANGLE OF INCIDENCE'); grid on;
    
    figure; pcolor(d_mag_v_AOI); shading flat;

    
    
end

    
    
%% function: Calc NonLinear Response
function [d,d_mag] = calc_NLO(S0,E0, Xtal_Tens)
    %normalize inputs (unit vectors)
    S0 = S0 / sqrt(dot(S0,S0));
    E0 = E0 / sqrt(dot(E0,E0));

    %define wave1 and wave2
    E1 = E0;    %fundamental1
    S1 = S0;
    E2 = E1;    %fundamental2 = fundamental1  (SHG = 3 wave mixing where 1&2 are same)
    S2 = S1;
    %compute magnetic field unit vector
    H0 = cross(S0,E0);
    H1 = cross(S1,E1);
    H2 = cross(S2,E2);

    %convert Crystal Tensor to coordinate system
    Xtal_M = convert_Tensor(Xtal_Tens);

    % compute d = NonLinear Response
    d(1) = E2*Xtal_M(:,:,1)*E1';
    d(2) = E2*Xtal_M(:,:,2)*E1';
    d(3) = E2*Xtal_M(:,:,3)*E1';
    d_mag = sqrt(dot(d,d));
    
    S3 = S1;
    


end %(calc Nonlinear Response)


%%
function M = convert_Tensor(T)
    M(:,:,1) = [T(1,1),T(6,1),T(5,1);...
                T(6,1),T(2,1),T(4,1);...
                T(5,1),T(4,1),T(3,1)];
    M(:,:,2) = [T(1,2),T(6,2),T(5,2);...
                T(6,2),T(2,2),T(4,2);...
                T(5,2),T(4,2),T(3,2)];
    M(:,:,3) = [T(1,3),T(6,3),T(5,3);...
                T(6,3),T(2,3),T(4,3);...
                T(5,3),T(4,3),T(3,3)];
end





