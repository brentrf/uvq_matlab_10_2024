%  SHG Efficiency Model(s)
%  based : Boyd 2.2.19 (with Area = Area_waveguide)
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

mod = calcAll(c,IN,mod,1);

%% SWEEP over Selected Parameters

% sweep1.fieldname  =  'Pp';
% sweep1.values      = logspace(-3,1,20)';  %[W]

sweep1.fieldname  =  'L_waveguide_um';
sweep1.values      = [20:20:2000]';  %[W]

sweep2.fieldname  =  'linewidth_nm';
sweep2.values      = [0: 0.1 :2]';  %[nm]

sweep3.fieldname  =  'Pump_Detune_nm';
sweep3.values     =  linspace(0,mod.linewidth_nm,10)';  %[nm]

sweep4.fieldname  =  'dispersion_slope_diff';
sweep4.values      = linspace(0,2*mod.dispersion_slope_diff,10)';  %[1/nm]

sweep5.fieldname  =  'Loss_intrinsic_dBcm';
sweep5.values      = [0:10:200]';  %[dB/cm]


 %LOOP OVER CASES
MR = []; do_plots=0;
for kc=1:length(INs)
    IN = INs{kc}; 
    legnames{kc} = INs{kc}.name;
    
    %Run this Case
    mod = IN; 
    mod = calcAll(c,IN,mod,do_plots);
    
    %Extract Results of Interest
    MR = ResultsMatrix_AddColumn(MR, mod);   
    
    %Parameter Sweep 1 ...........................
    for k = 1:length(sweep1.values)
      %change Parameter
        IN.(sweep1.fieldname)           =  sweep1.values(k,1)';
        mod = IN;
      %run Model
        sweep1.mod(k)                   =  calcAll(c,IN,mod,do_plots);           
      %Record key results
        sweep1.Effic(k,kc)          	= sweep1.mod(k).Efficiency;        
        sweep1.Effic_ZeroLinewidth(k,kc)= sweep1.mod(k).Efficiency_ZeroLinewidth;
        sweep1.DeRate(k,kc)             = sweep1.mod(k).DeRate_NonZeroLinewidth;
        sweep1.Lcoh_um(k,kc)            = sweep1.mod(k).Lcoh_nm * 1e-3;
    end
  
    %Parameter Sweep 2 ..........................
    for k = 1:length(sweep2.values);
      %change Parameter
        IN.(sweep2.fieldname)    =  sweep2.values(k,1)';          
        mod = IN;
      %run Model
        sweep2.mod(k)            =  calcAll(c,IN,mod,do_plots);             
      %Record key results
        sweep2.Effic(k,kc)                  = sweep2.mod(k).Efficiency;        
        sweep2.Effic_ZeroLinewidth(k,kc)    = sweep2.mod(k).Efficiency_ZeroLinewidth;
        sweep2.Lcoh_um(k,kc)                = sweep2.mod(k).Lcoh_nm * 1e-3;
    end
    
     %Parameter Sweep 3 ..........................
    for k = 1:length(sweep3.values);
      %change Parameter
        IN.(sweep3.fieldname)    =  sweep3.values(k,1)';          
        mod = IN;
      %run Model
        sweep3.mod(k)            =  calcAll(c,IN,mod,do_plots);             
      %Record key results
        sweep3.Effic(k,kc)                  = sweep3.mod(k).Efficiency;        
        sweep3.Effic_ZeroLinewidth(k,kc)    = sweep3.mod(k).Efficiency_ZeroLinewidth;
        sweep3.Lcoh_um(k,kc)                = sweep3.mod(k).Lcoh_nm * 1e-3;
    end   
    
     %Parameter Sweep 4 ..........................
    for k = 1:length(sweep4.values);
      %change Parameter
        IN.(sweep4.fieldname)    =  sweep4.values(k,1)';          
        mod = IN;
      %run Model
        sweep4.mod(k)            =  calcAll(c,IN,mod,do_plots);             
      %Record key results
        sweep4.Effic(k,kc)                  = sweep4.mod(k).Efficiency;        
        sweep4.Effic_ZeroLinewidth(k,kc)    = sweep4.mod(k).Efficiency_ZeroLinewidth;
        sweep4.Lcoh_um(k,kc)                = sweep4.mod(k).Lcoh_nm * 1e-3;
    end   
    
     %Parameter Sweep 5..........................
    for k = 1:length(sweep5.values);
      %change Parameter
        IN.(sweep5.fieldname)    =  sweep5.values(k,1)';          
        mod = IN;
      %run Model
        sweep5.mod(k)            =  calcAll(c,IN,mod,do_plots);             
      %Record key results
        sweep5.Effic(k,kc)                  = sweep5.mod(k).Efficiency;        
        sweep5.Effic_ZeroLinewidth(k,kc)    = sweep5.mod(k).Efficiency_ZeroLinewidth;
        sweep5.Lcoh_um(k,kc)                = sweep5.mod(k).Lcoh_nm * 1e-3;
    end     
    
end

%Write Results for All Cases
csvwrite(outputfile, MR);
disp(['Wrote output to:', outputfile ]);


%--------------- Plot: EFFICIENCY (for Each Sweep)

figure; semilogx(sweep1.values, sweep1.Effic,'Linewidth',2); grid on;
xlabel(replace(sweep1.fieldname,'_',' ')); ylabel('Efficiency (avg) [0...1]');
legend(legnames);  title(replace(sweep1.fieldname,'_',' '));

figure; plot(sweep2.values, sweep2.Effic,'Linewidth',2); grid on;
xlabel(replace(sweep2.fieldname,'_',' ')); ylabel('Efficiency (avg) [0...1]');
legend(legnames);  title(replace(sweep2.fieldname,'_',' '));

figure; plot(sweep3.values, sweep3.Effic,'Linewidth',2); grid on;
xlabel(replace(sweep3.fieldname,'_',' ')); ylabel('Efficiency (avg) [0...1]');
legend(legnames);  title(replace(sweep3.fieldname,'_',' '));

figure; plot(sweep4.values, sweep4.Effic,'Linewidth',2); grid on;
xlabel(replace(sweep4.fieldname,'_',' ')); ylabel('Efficiency (avg) [0...1]');
legend(legnames);  title(replace(sweep4.fieldname,'_',' '));

figure; plot(sweep5.values, sweep5.Effic,'Linewidth',2); grid on;
xlabel(replace(sweep5.fieldname,'_',' ')); ylabel('Efficiency (avg) [0...1]');
legend(legnames);  title(replace(sweep5.fieldname,'_',' '));

%--------------- Plot: Coherence Length

figure; plot(sweep2.values, sweep2.Lcoh_um,'Linewidth',2); grid on;
xlabel(replace(sweep2.fieldname,'_',' ')); ylabel('Coherence Length [um]');
legend(legnames);  title(replace(sweep2.fieldname,'_',' '));
figure; plot(sweep3.values, sweep3.Lcoh_um,'Linewidth',2); grid on;
xlabel(replace(sweep3.fieldname,'_',' ')); ylabel('Coherence Length [um]');
legend(legnames);  title(replace(sweep3.fieldname,'_',' '));
figure; plot(sweep4.values, sweep4.Lcoh_um,'Linewidth',2); grid on;
xlabel(replace(sweep4.fieldname,'_',' ')); ylabel('Coherence Length [um]');
legend(legnames);  title(replace(sweep4.fieldname,'_',' '));




%%  Calc All
function mod = calcAll(c,IN,mod,do_plots)
    mod = IN;
    mod = calc_prefactor(c,mod); 
    mod = Efficiency_ZeroLinewidth(c,mod);
    mod = DeRate_NonZeroLinewidth(c,mod,do_plots);
    
end

%%
function mod = Efficiency_ZeroLinewidth(c,mod)
    mod.A_waveguide_um2  = mod.x_waveguide_um * mod.y_waveguide_um;  %[um^2]
    mod.GAMMA_pctPerWatt = mod.prefactor * mod.overlap_frac^2 ...
                            * (mod.L_waveguide_um^2 / mod.A_waveguide_um2);  
    mod.Efficiency_ZeroLinewidth = mod.GAMMA_pctPerWatt * mod.Pp;
end

%%
function mod = DeRate_NonZeroLinewidth(c,mod,do_plots)
    
    %Coherence Length (of PumpLinewidth)
    mod.Lcoh_nm = (2 * 1.39156) /(2*pi) * (mod.lam1/mod.linewidth_nm) / mod.dispersion_slope_diff;
        
    %Bandwidth of Phase Matching  
    HWHM = mod.linewidth_nm/2;
    deltas_lam_nm = [-2*HWHM:0.001:2*HWHM]'; %[nm]
    for k = 1:length(deltas_lam_nm)
        mod = PhaseMismatch_Derate(c,mod,deltas_lam_nm(k));
        DeRate(k,1) = mod.DeRate;
        DeRateAppx(k,1) = mod.DeRateAppx;
    end
    
    %Overlap:  Linewidth(pump) x Bandwidth(Phasematch)
    Ypump = exp(-0.5*((deltas_lam_nm - mod.Pump_Detune_nm).^2)/(HWHM)^2 );
    mod.DeRate_NonZeroLinewidth = sum(Ypump.*DeRate)/sum(Ypump);
    mod.Efficiency = mod.Efficiency_ZeroLinewidth * mod.DeRate_NonZeroLinewidth;    
    
    
    %plot: Linewidth(Pump) v. Bandwidth(PhaseMatching)
    if do_plots==1
        figure; plot(deltas_lam_nm,[DeRate,DeRateAppx]);
        hold on; plot(deltas_lam_nm,Ypump);
        legend('DeRate: PhaseMatching','Pump_Linewidth')
        legend('DeRate: PhaseMatching (exact)', 'DeRate: PhaseMatching (appx)','Pump_Linewidth')
        title({'PhaseMatch Bandwidth v. Pump Linewidth';...
            ['derate value = ',num2str(mod.DeRate_NonZeroLinewidth)]});
    end
    

    
end

%%
function mod = PhaseMismatch_Derate(c,mod, delta_lam_nm)
    delta_lam   = delta_lam_nm; %(mod.linewidth_nm/2) * 2;
    delta_n     = delta_lam * mod.dispersion_slope_diff;    %unitless
    delta_k_cm  = 4*pi/(mod.lam1*1e-7) * delta_n;           %[cm^-1]
    mod.alpha_cm = mod.Loss_intrinsic_dBcm / log(10);       % alpha[cm-1]   =  [dB/cm] /ln(10)
    
    dkL  = delta_k_cm   * (mod.L_waveguide_um * 1e-4);
    aL   = mod.alpha_cm * (mod.L_waveguide_um * 1e-4); 
 
    if dkL>0==0 && aL==0
        mod.DeRate = 1;
        mod.DeRateAppx=  1;
    else
        mod.DeRate =  2 * ( cosh(aL) - cos(dkL) ) * exp(-aL) / ( (dkL)^2 + (aL)^2 ) ;
        mod.DeRateAppx = 2 *( sinc(dkL) )^2 * (exp(-aL));
    end
    
%     if dkL>0
%         mod.DeRate =  2 * ( cosh(aL) - cos(dkL) ) * exp(-aL) / ( (dkL)^2 + (aL)^2 ) ;
%         mod.DeRateAppx = 2 *( sinc(dkL) )^2 * (exp(-aL));
%     else
%         if aL>0
%             mod.DeRate =  exp(-aL);
%             mod.DeRateAppx =  exp(-aL);
%         else
%             mod.DeRate = 1;
%             mod.DeRateAppx = 1;
%         end
%     end
    
    
end


%%  Prefactor Calc
function mod = calc_prefactor(c,mod) 
    mod.lam2            = mod.lam1 / 2;
    mod.w_1             = 2*pi* (c.c0/mod.n1) / (  mod.lam1 * 1e-9); %[Hz]
    mod.w_2             = 2*pi* (c.c0/mod.n2) / (  mod.lam2 * 1e-9); %[Hz]
    mod.prefactor       = 2*mod.w_2^2*(mod.d_xx*1e-12)^2 / (c.eps0*mod.n1^2*mod.n2*c.c0^3)
end


%% get INPUT 
function INs = get_INPUT(dirs,inputfile)
    sheetname = 'Linear';
    %get input
    [tmp.num,tmp.txt,tmp.raw]       =  xlsread([dirs.h,'/',inputfile],sheetname, 'd5:z20' );
    [tmp2.num,tmp2.txt,tmp2.raw]    =  xlsread([dirs.h,'/',inputfile],sheetname, 'd4:z4' );

    if size(tmp.num,2)>=1 && size(tmp.num,1)>=10

        for kc=1:size(tmp.num,2)
            INs{kc}.name             = replace(tmp2.txt{kc},'_',' ');
            INs{kc}.Pp               = tmp.num(1,kc);   %[W] external pump (into waveguide), average power
            INs{kc}.linewidth_nm     = tmp.num(2,kc);   % [nm]
            INs{kc}.d_xx             = tmp.num(3,kc);   %pm/V
            INs{kc}.lam1             = tmp.num(4,kc);   
            INs{kc}.n1               = tmp.num(5,kc);   
            INs{kc}.n2               = tmp.num(6,kc);   
            INs{kc}.overlap_frac     = tmp.num(7,kc);  %fraction (of 1)
            INs{kc}.x_waveguide_um   = tmp.num(8,kc);  
            INs{kc}.y_waveguide_um   = tmp.num(9,kc); 
            INs{kc}.L_waveguide_um   = tmp.num(10,kc); 
            INs{kc}.dispersion_slope_diff   = tmp.num(11,kc); 
            INs{kc}.Loss_intrinsic_dBcm  = tmp.num(12,kc); 
            INs{kc}.Pump_Detune_nm   = tmp.num(13,kc); 
        end
        Ncases = length(INs);
        
    else
        error('Excel Input was not obtained... Returning Default INPUTs ')
    end
    
end


%% Generate Output Matrix
function M_out = ResultsMatrix_AddColumn(MR, mod)
   newres      = zeros(5,1);
   newres(1,1) = mod.Lcoh_nm*1e-3;
   newres(2,1) = mod.prefactor;
   newres(3,1) = mod.GAMMA_pctPerWatt; 
   newres(4,1) = mod.Efficiency_ZeroLinewidth;
   newres(5,1) = mod.DeRate_NonZeroLinewidth;
   newres(6,1) = mod.Efficiency;
    
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





