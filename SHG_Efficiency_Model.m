%  SHG Efficiency Model(s)
%  based on Breunig & Sturman
%  Parts
%    [a] Efficiency Model for Ring + Bus
%    [b] Compare Guo Model (non-depletion)
%    [c] Mode Index = f(Geometry) .... 
%    [d] Overlap = 
%

close all;
clear;
c.c0                = 3e8; %[m/s]
c.eps0              = 8.85E-12;  %[C/(V.m)]


dirs.h = 'C:\Users\brent\MATLAB\matDat';
inputfile =  'SHGEfficModel(matlab)_INPUT.xlsx';
outputfile = 'temp.csv';


%%
INs = get_INPUT(dirs,inputfile);



%% set up / test
IN = INs{1};
mod = IN;       %copy input values to this model ('mod')
mod = calcAll(c,IN,mod);



%% SWEEP over Selected Parameters
sweep1.fieldname  =  'Pp_avg';
sweep1.values      = logspace(-3,1,20)';  %[W]

sweep2.fieldname  =  'detune1_nm';
sweep2.values      = [0: 0.2 :2]';  %[nm]

MR = []; 
for kc=1:length(INs)
    IN = INs{kc}; 
    legnames{kc} = INs{kc}.name;
    
    %Run this Case
    mod = IN; 
    mod = calcAll(c,IN,mod);
    
    %Extract Results of Interest
    MR = ResultsMatrix_AddColumn(MR, mod);   
    
    %Parameter Sweep1 ...........................
    for k = 1:length(sweep1.values);
      %change Parameter
        IN.(sweep1.fieldname)    =  sweep1.values(k,1)';          
      %run Model
        sweep1.mod(k)            =  calcAll(c,IN,mod);           
      %Record key results
        sweep1.Effic(k,kc)          = sweep1.mod(k).Efficiency_avg;        
        sweep1.Effic_detuned(k,kc)  = sweep1.mod(k).det.Efficiency_avg;    
    end
  
    %Parameter Sweep2 ..........................
    for k = 1:length(sweep2.values);
      %change Parameter
        IN.(sweep2.fieldname)    =  sweep2.values(k,1)';          
      %run Model
        sweep2.mod(k)            =  calcAll(c,IN,mod);           
      %Record key results
        sweep2.Effic(k,kc)          = sweep2.mod(k).Efficiency_avg;        
        sweep2.Effic_detuned(k,kc)  = sweep2.mod(k).det.Efficiency_avg;      
    end
    
    
end
csvwrite(outputfile, MR);
disp(['Wrote output to:', outputfile ]);

figure; semilogx(sweep1.values, sweep1.Effic_detuned,'Linewidth',2); grid on;
xlabel(sweep1.fieldname); ylabel('Efficiency (avg) [0...1]');
legend(legnames);

figure; plot(sweep2.values, sweep2.Effic_detuned,'Linewidth',2); grid on;
xlabel(sweep2.fieldname); ylabel('Efficiency (avg) [0...1]');
legend(legnames);



 



%%  Calc All
function mod = calcAll(c,IN,mod);
    mod = calc_P0_prefactor(c,IN,mod);
    mod = ModeVolume(c,IN,mod);
    mod = Q_calcs(c,IN,mod);
    mod = DeTuning(c,IN,mod);
    mod = Efficiency_OnRes(c,IN,mod);
    mod = Efficiency_DeTune(c,IN,mod);
end


%% Calc Coupling Coefficients

% mod.b.gamma2 = mod.w_2 * mod.d_xx * ( 2/(c.eps0 *c.c0^3 * mod.n_1^2 * mod.n_2 ) )^0.5 ;  %[1/m] eq.35 in stirman,breunig

% mod.g.g = 



%%  Prefactor Calc
function mod = calc_P0_prefactor(c,IN,mod) 
    mod.w_1             = 2*pi* (c.c0/IN.n_1) / (  IN.lam1    * 1e-9); %[Hz]
    mod.w_2             = 2*pi* (c.c0/IN.n_2) / ( (IN.lam1/2) * 1e-9); %[Hz]
    mod.P0_prefactor    = (0.5*mod.w_1)*c.eps0* IN.n_1^4*IN.n_2^2/(8*(IN.d_xx*1e-12)^2); %[W/m3]
    mod.lam2            = mod.lam1 / 2;
end

%% Mode Volume & Overlap
function mod = ModeVolume(c,IN,mod)
    mod.L_um            = 2*pi*IN.R_um;     %ring circumference
    mod.m1              = mod.L_um*1000 / (IN.lam1/IN.n_1);  % number of optical cycles per L  @ lam1	#	224.98
    mod.m2              = 2 * mod.m1;   % enforce phase matching:  must be 2 * m1
    % mod.m2              = mod.L_um*1000 / ((IN.lam1/2)/IN.n_2);  % number of optical cycles per L  @ lam1	#	224.98
    
    mod.V_1_um3         = IN.x_waveguide_um * IN.y_waveguide_um * (2*pi*IN.R_um);
    mod.V_2_um3         = mod.V_1_um3*0.8;  %?????????????          %
    mod.V_ratio         = ( mod.V_1_um3^2 * mod.V_2_um3 ) / ( IN.V_overlap^2 * mod.V_1_um3 * mod.V_2_um3 ) * (1e-6)^3;
                         % Volume Ratio  = (Vs * Vp^2) / (Vps)   [m3]
end

%% Q Estimates      
function mod = Q_calcs(c,IN,mod)
% Intrinsic Q  
    mod.alpha1_cm1      = IN.alpha1_dbcm / 10*log(10);
    mod.alpha2_cm1      = IN.alpha2_dbcm / 10*log(10);
    mod.kappa_abs_p     = sqrt( mod.alpha1_cm1 * (mod.L_um*0.0001) ); %material abs loss rate (unitless) at lam1
    mod.kappa_abs_s     = sqrt( mod.alpha2_cm1 * (mod.L_um*0.0001) ); %material abs loss rate (unitless) at lam2
    mod.kappa_bend_p    = IN.kappa_bend_p; %bending loss rate (unitless) at lam1
    mod.kappa_bend_s    = IN.kappa_bend_s; %bending  loss rate (unitless) at lam2
    mod.kappa_scatt_p   = IN.kappa_scatt_p; %scattering loss rate (unitless) at lam1
    mod.kappa_scatt_s   = IN.kappa_scatt_s; %scattering loss rate (unitless) at lam2
    mod.kappa0_p        = (mod.kappa_abs_p^2 + mod.kappa_bend_p^2 + mod.kappa_scatt_p^2)^0.5;  %[total] intrinsic loss rate (unitless) at lam1
    mod.kappa0_s        = (mod.kappa_abs_s^2 + mod.kappa_bend_s^2 + mod.kappa_scatt_s^2)^0.5;  %[total] intrinsic loss rate (unitless) at lam2
    mod.delta_W0_p      = (2*pi)*(c.c0/mod.n_1)/(mod.L_um*0.000001) * mod.kappa0_p^2;          % intrinsic linewidth in Hz
    mod.delta_W0_s      = (2*pi)*(c.c0/mod.n_2)/(mod.L_um*0.000001) * mod.kappa0_s^2;          % intrinsic linewidth in Hz
    
    mod.Q0_p            = mod.w_1 / mod.delta_W0_p; 
    mod.Q0_s            = mod.w_2 / mod.delta_W0_s; 

  % Coupled Q  
    mod.kappa_couple_p  = sqrt(IN.rp)*mod.kappa0_p;    %[unitless]  coupling rate at lam1
    mod.kappa_couple_s  = sqrt(IN.rs)*mod.kappa0_s;    % [unitless]coupling rate at lam2
    mod.delta_Wc_p      = mod.delta_W0_p*(1+IN.rp);      % intrinsic linewidth in Hz
    mod.delta_Wc_s      = mod.delta_W0_s*(1+IN.rs);   %intrinsic linewidth in Hz
    %mod.delta_NU0_p     = mod.delta_W0_p / (2*pi) * 0.000001;   %intrinsic linewidth in Hz   (= W0 / 2pi)
    %mod.delta_NU0_s     = mod.delta_W0_s / (2*pi) * 0.000001; %intrinsic linewidth in Hz   (= W0 / 2pi)
    mod.Qc_p            = mod.Q0_p / (1+IN.rp); 
    mod.Qc_s            = mod.Q0_s / (1+IN.rs);    
    
end

%%  Detuning
function mod = DeTuning(c,IN,mod)
    IN.detune2_nm = IN.detune1_nm/2;
    mod.detune1_Hz  = (2*pi)*c.c0/(mod.n_1*mod.m1)*(1/mod.lam1 - 1/(mod.lam1 +  IN.detune1_nm))*1e9;  %??? divide by m1  ?    detuning of pump in Hz
    mod.detune2_Hz  = (2*pi)*c.c0/(mod.n_2*mod.m2)*(2/mod.lam1 - 2/(mod.lam1 +  IN.detune2_nm))*1e9;  %??? divide by m1  ?    detuning of SHG in Hz

    mod.detune_ratio_p =  mod.detune1_Hz / mod.delta_Wc_p; %detuning ratio (@ pump)
    mod.detune_ratio_s =  mod.detune2_Hz / mod.delta_Wc_s; %detuning ratio (@ pump)
   
    mod.Qc_detun_p  = mod.Q0_p / (1+IN.rp + mod.detune_ratio_p^2);   %????
    mod.Qc_detun_s  = mod.Q0_s / (1+IN.rs + mod.detune_ratio_s^2);   %????    
end

%% Efficiency v Pump Power -- On Resonance
function mod = Efficiency_OnRes(c,IN,mod)
    mod.Pp              = IN.Pp_avg / IN.DF;  % external pump (into waveguide), instantaneous
    mod.P0              = mod.P0_prefactor * mod.V_ratio / (mod.Qc_p^2 * mod.Qc_s);   % prefactor * Vratio / (Qp^2 * Qs)
    mod.Max_Efficiency  = (IN.rp/(1+IN.rp))*(IN.rs/(1+IN.rs))*IN.DeRate;

    mod.Pp_max          = 4*(1+IN.rp)/IN.rp * mod.P0;    %Pump Power Req. to reach Max Efficiency

    mod.GAMMA           = 4*4*(IN.rp/(1+IN.rp))*(IN.rs/(1+IN.rs))*(1/mod.P0);  % := (Slope Efficiency / Pp at low pump)   =  1/Pp_max
    mod.A               = mod.Pp / mod.P0 * IN.rp/(1 + IN.rp); %= Pp / Po * (rp/(1+rp))
    
    mod.Xt              = (3*sqrt(81*mod.A^2+12*mod.A)+27*mod.A+2)^(1/3)  /  (2^(1/3)); %= X_term1 = f(A)
    mod.X               = (mod.Xt + 1/mod.Xt - 2 )/3;  % =  f(A)    :=  Ppump_intracavity_normalized
    mod.X_dblchk        = mod.X*(1+mod.X)^2  -  mod.A ; % X*(1+X)^2 = A 
                    
    if (mod.X_dblchk < 1e-12) disp('OK'); end

    mod.Psh             = 4*IN.rs / (1+IN.rs) * mod.P0 * mod.X^2 * IN.DeRate; %output SHG power, instantaneous
    mod.Efficiency      = mod.Psh / mod.Pp;   %[%] abs. efficiency
    mod.Efficiency_dblchk = 4*(IN.rs/(1+IN.rs)) * (IN.rp/(1+IN.rp)) * mod.X / (1+mod.X)^2 * IN.DeRate;
                    
    if (abs(mod.Efficiency_dblchk - mod.Efficiency) < 1e-12) disp('OK'); end
    
    mod.Psh_avg         = mod.Psh*IN.DF;                %[W] average output SHG power
    mod.Efficiency_avg  = mod.Psh_avg / IN.Pp_avg;

end

%% Efficiency v Pump Power -- with Detuning
function mod = Efficiency_DeTune(c,IN,mod)
    detun1 = mod.detune_ratio_p;% temp variables
    detun2 = mod.detune_ratio_s;

    mod.det.Pp              = mod.Pp;      % external pump (into waveguide), instantaneous
    mod.det.P0              = mod.P0;     % prefactor * Vratio / (Qp^2 * Qs)
    
    mod.det.Max_Effic       = (IN.rp/(1 + IN.rp + detun1^2))*(IN.rs/(1 + IN.rs + detun2^2))*IN.DeRate;
    mod.det.Pp_max          = 4*(1+ IN.rp + detun1^2)/IN.rp * mod.P0;    %Pump Power Req. to reach Max Efficiency
    mod.det.GAMMA           = 4*4*(IN.rp/(1+IN.rp + detun1^2))*(IN.rs/(1+IN.rs + detun2^2))*(1/mod.P0);  % := (Slope Efficiency / Pp at low pump)   =  1/Pp_max
    mod.det.A               = mod.det.Pp / mod.det.P0 * IN.rp/(1 + IN.rp + detun1^2); %= Pp / Po * (rp/(1+rp+detun1^2))
    
    mod.det.Xt              = (3*sqrt(81*mod.det.A^2+12*mod.det.A)+27*mod.det.A+2)^(1/3)  /  (2^(1/3)); %= X_term1 = f(A)
    mod.det.X               = (mod.det.Xt + 1/mod.det.Xt - 2 )/3;  % =  f(A)    :=  Ppump_intracavity_normalized
    X_dblchk        = mod.det.X*(1+mod.det.X)^2  -  mod.det.A ; % X*(1+X)^2 = A 
                   
    if (X_dblchk < 1e-12) disp('OK'); end

    mod.det.Psh             = 4*IN.rs / (1+IN.rs + detun2^2) * mod.det.P0 * mod.det.X^2 * IN.DeRate; %output SHG power, instantaneous
    mod.det.Efficiency      = mod.det.Psh / mod.det.Pp;   %[%] abs. efficiency
    Efficiency_dblchk = 4*(IN.rs/(1+IN.rs+ detun2^2)) * (IN.rp/(1+IN.rp + detun1^2)) * mod.det.X / (1+mod.det.X)^2 * IN.DeRate;
                    
    if (abs(Efficiency_dblchk - mod.det.Efficiency) < 1e-12) disp('OK'); end
    
    mod.det.Psh_avg         = mod.det.Psh*IN.DF;                %[W] average output SHG power
    mod.det.Efficiency_avg  = mod.det.Psh_avg / IN.Pp_avg;
  
           
end


%% get INPUT 
function INs = get_INPUT(dirs,inputfile)
    sheetname = 'Input_EfficModel';
    %get input
    [tmp.num,tmp.txt,tmp.raw]       =  xlsread([dirs.h,'/',inputfile],sheetname, 'd3:z100' );
    [tmp2.num,tmp2.txt,tmp2.raw]    =  xlsread([dirs.h,'/',inputfile],sheetname, 'd2:z2' );
    
    
    if size(tmp.num,2)>=1 && size(tmp.num,1)>=20

        for kc=1:size(tmp.num,2)
            INs{kc}.name             = replace(tmp2.txt{kc},'_',' ');
            INs{kc}.DeRate           = tmp.num(1,kc);   % (to match measured result)
            INs{kc}.Pp_avg           = tmp.num(2,kc);   %[W] external pump (into waveguide), average power
            INs{kc}.DF               = tmp.num(3,kc);   % (1  = CW)
            INs{kc}.rp               = tmp.num(4,kc);   
            INs{kc}.rs               = tmp.num(5,kc);   

            INs{kc}.d_xx             = tmp.num(6,kc);   %pm/V
            INs{kc}.lam1             = tmp.num(7,kc);   %[nm]
            INs{kc}.n_1              = tmp.num(8,kc);   %unitless
            INs{kc}.n_2              = tmp.num(9,kc);  %unitless

            INs{kc}.x_waveguide_um   = tmp.num(10,kc);  
            INs{kc}.y_waveguide_um   = tmp.num(11,kc); 
            INs{kc}.R_um             = tmp.num(12,kc);           %ring radius
            INs{kc}.V_overlap        = tmp.num(13,kc); %percentage

            INs{kc}.alpha1_dbcm      = tmp.num(14,kc);    %intrinsic loss @  450nm
            INs{kc}.alpha2_dbcm      = tmp.num(15,kc);    %intrinsic loss @  225nm
            INs{kc}.kappa_bend_p     = tmp.num(16,kc);   %bending loss rate (unitless) at lam1
            INs{kc}.kappa_bend_s     = tmp.num(17,kc);   %bending  loss rate (unitless) at lam2
            INs{kc}.kappa_scatt_p    = tmp.num(18,kc);   %scattering loss rate (unitless) at lam1
            INs{kc}.kappa_scatt_s    = tmp.num(19,kc);   %scattering loss rate (unitless) at lam2

            INs{kc}.detune1_nm       = tmp.num(20,kc);         %detuning Ring rel. to  pump in lamba[nm]
        end
        Ncases = length(INs);
        
    else
        warning('Excel Input was not obtained... Returning Default INPUTs ')

        IN.Pp_avg           = 0.1;           %[W] external pump (into waveguide), average power
        IN.DF               = 1;             % (1  = CW)
        IN.rp               = 1;
        IN.rs               = 1;
        IN.DeRate           = 0.5;          % (to match measured result)

        IN.d_xx             = 6.2;          %pm/V
        IN.lam1             = 1550;         %[nm]
        IN.n_1              = 1.85          %unitless
        IN.n_2              = 2.1;          %unitless

        IN.x_waveguide_um   = 1.3;   
        IN.y_waveguide_um   = 1.0054;
        IN.R_um             = 30;           %ring radius
        IN.V_overlap         = 0.009; 

        IN.alpha1_dbcm       = 0.025;       %intrinsic loss
        IN.alpha2_dbcm       = 0.11;        %
        IN.kappa_bend_p      = 0; %bending loss rate (unitless) at lam1
        IN.kappa_bend_s      = 0; %bending  loss rate (unitless) at lam2
        IN.kappa_scatt_p     = 0; %scattering loss rate (unitless) at lam1
        IN.kappa_scatt_s     = 0; %scattering loss rate (unitless) at lam2

        IN.detune1_nm        = 0.1;         %detuning Ring rel. to  pump in lamba[nm]
        IN.detune2_nm        = 0.05;        %detuning Ring rel. to SHG in lamba[nm]
        
        INs{1} = IN;
    end
    
end

%% Generate Output Matrix
function M_out = ResultsMatrix_AddColumn(MR, mod);
    newres = zeros(13,1);
    newres(1,1) = mod.Q0_p; 
    newres(2,1) = mod.Q0_s;
    newres(3,1) = mod.Qc_p;
    newres(4,1) = mod.Qc_s;
    newres(5,1) = mod.V_ratio;
    newres(6,1) = mod.P0_prefactor;

    newres(7,1) = mod.Max_Efficiency;
    newres(8,1) = mod.Pp_max;
    newres(9,1) = mod.GAMMA;
    newres(10,1) = mod.A;
    newres(11,1) = mod.X;
    newres(12,1) = mod.Psh;
    newres(13,1) = mod.Efficiency;  

    newres(14,1) = mod.det.Max_Effic;
    newres(15,1) = mod.det.Pp_max;
    newres(16,1) = mod.det.GAMMA;
    newres(17,1) = mod.det.A;
    newres(18,1) = mod.det.X;
    newres(19,1) = mod.det.Psh;
    newres(20,1) = mod.det.Efficiency; 
    
   [nr,nc] = size(MR);
   if nr==0
       M_out = newres;
   else
       if nr == length(newres)
           M_out = [MR,newres];
       else
           error('Results List does not match size of Matrix, MR');
       end
   end   
end





