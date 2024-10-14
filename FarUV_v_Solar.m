% SPECTRUM:  FAR-UV v. SOLAR 

function FarUV_v_Solar()

%INPUTs
dirs.d = '/Users/brentfisher/Documents/MATLAB/UVQ/matDat/';
LampFile = '/UniLam_EX150S07_222nm_SPECTRUM.csv';  %spectrum: KrCl
% LampFile = '/KrBr_spec1.csv';        %spectrum: KrBr
LampFile = '/KrCl_spec0.csv';          %spectrum: KrCl

LampIrrad_fullspectrum_mWcm2 = 0.005;   %[mW/cm2] assume this is IRRADIANCE of the *FULL SPECTRUM*
SunIrrad_fullspectrum_Wm2    = 100;   %[W/m2] assumed sun intensity
lam_upper_lim                = 400;      

%%
% READ LAMP Data .....
tmp = csvread([dirs.d,LampFile],1,0);
lamp0.nm         = tmp(:,1);
lamp0.relPow     = tmp(:,2);

%TLV Data ....
tmp = csvread([dirs.d,'/TLV.csv'],1,0);
TLV0.nm         = tmp(:,1);
TLV0.mWcm2.old  = tmp(:,2);
TLV0.mWcm2.new  = tmp(:,3);
indx270         = find(TLV0.nm == 270);
TLV0.S.old      = TLV0.mWcm2.old(indx270)./TLV0.mWcm2.old;
TLV0.S.new      = TLV0.mWcm2.new(indx270)./TLV0.mWcm2.new;
% figure; plot(TLV.nm,TLV.S.old); hold all; plot(TLV.nm,TLV.S.new);
% xlim([200,300]); title('Raw TLV, as read');

solar = load([dirs.d,'/AM15.mat']);



% INTERPOLATE TLV & Lamp
lambdas         = [180:1:lam_upper_lim]';
TLV.nm          = lambdas;
TLV.mWcm2.old   = interp1(TLV0.nm,TLV0.mWcm2.old,lambdas);
    TLV.mWcm2.old(find(isnan(TLV.mWcm2.old))) = 0;
TLV.mWcm2.new   = interp1(TLV0.nm,TLV0.mWcm2.new,lambdas);
    TLV.mWcm2.new(find(isnan(TLV.mWcm2.new))) = 0;
indx270         = find(TLV.nm == 270);
TLV.S.old       = TLV.mWcm2.old(indx270)./TLV.mWcm2.old;
TLV.S.new       = TLV.mWcm2.new(indx270)./TLV.mWcm2.new;        % S[mW/cm2] = inverse sensitivity  
lamp.nm         = lambdas;
lamp.relPow     = interp1(lamp0.nm,lamp0.relPow,lambdas);  
    lamp.relPow(find(isnan(lamp.relPow))) = 0;
TLV.lamp.relPow  = lamp.relPow;

TLV.solar.absPow_Wm2nm = interp1(solar.AM15.lam.nm, solar.AM15.G.Wm2_per_nm, lambdas) * (SunIrrad_fullspectrum_Wm2/1000)  ;  % @ 1000W/m2  
dlams = diff(lambdas); dlams = [dlams;dlams(end)];

 
%% Compute TLV-Weighted Spectra  (from: 200-300nm)










% %interpolate & return result
% disp(TLV.lamp.relPow);

%apply ideal filter (10nm wide around peak)
FilterFullWidth_nm      = 15;
[maxval,imax]           = max(lamp.relPow); lamp.peak_nm = lamp.nm(imax);
TLV.lampwfilt.relPow    = TLV.lamp.relPow ...
                         .* ( abs(TLV.nm - lamp.peak_nm) <= FilterFullWidth_nm/2 );  %apply filter window




%% Compute Effective Irradiance (from: 200-300nm)

%[1] % Normalize Lamp Spectra  (to Total Power, i.e. Area Under Curve)
TLV.lamp.relPow         = TLV.lamp.relPow       / sum(TLV.lamp.relPow);
TLV.lampwfilt.relPow    = TLV.lampwfilt.relPow  / sum(TLV.lampwfilt.relPow);

%[2a] % get EFFECTIVE Lamp Spectra (to Total LAMP Power, i.e. Area Under Curve)
 
ikeep = find( (TLV.nm <=lam_upper_lim) & ~isnan(TLV.lamp.relPow + TLV.lampwfilt.relPow +  TLV.mWcm2.old + TLV.mWcm2.new ));
peaksum = sum(TLV.lampwfilt.relPow(ikeep));  
EffectiveIrradiance_ofLAMP_mWcm2__old    = sum(TLV.S.old(ikeep) .* TLV.lampwfilt.relPow(ikeep) )/peaksum * LampIrrad_fullspectrum_mWcm2;
EffectiveIrradiance_ofLAMP_mWcm2__new    = sum(TLV.S.new(ikeep) .* TLV.lampwfilt.relPow(ikeep) )/peaksum * LampIrrad_fullspectrum_mWcm2;

%[2b] % get EFFECTIVE *Solar*  Spectra  

ikeep = find( (TLV.nm <=lam_upper_lim) & ~isnan(TLV.solar.absPow_Wm2nm  +  TLV.mWcm2.old + TLV.mWcm2.new ));
EffectiveIrradiance_ofSOLAR_mWcm2__old    = sum(TLV.S.old(ikeep) .* TLV.solar.absPow_Wm2nm(ikeep) .* dlams(ikeep) ) * (1000/10000);
EffectiveIrradiance_ofSOLAR_mWcm2__new    = sum(TLV.S.new(ikeep) .* TLV.solar.absPow_Wm2nm(ikeep) .* dlams(ikeep) ) * (1000/10000);
 
%% Generate Comparison Spectrum:  Lamp v. Solar

figure; subplot(3,1,1);
    semilogy(lambdas,  TLV.S.new);   ylabel('Relative Risk Level (1/TLV)');
    yyaxis right;
    plot(lambdas ,               TLV.solar.absPow_Wm2nm ); hold on;
    plot(lambdas,  TLV.lampwfilt.relPow  * (LampIrrad_fullspectrum_mWcm2/peaksum) );
    ylabel('mW/cm2/nm Irradiance');
    title({['Total(spec.integrated) LAMP Irradiance = ',num2str(LampIrrad_fullspectrum_mWcm2),'mW/cm2'];...
           ['Total(spec.integrated) SOLAR Irradiance = ',num2str(SunIrrad_fullspectrum_Wm2),'W/m2']})  
    xlim([200,lam_upper_lim]);
    
subplot(3,1,2);
    plot(lambdas , TLV.S.new  .* TLV.solar.absPow_Wm2nm ); hold on;
    plot(lambdas,  TLV.S.new  .* TLV.lampwfilt.relPow  * (LampIrrad_fullspectrum_mWcm2/peaksum) );
    ylabel('relative intensity of photon risk')
    xlim([200,lam_upper_lim]);
    legend('Solar','Lamp');
    
subplot(3,1,3);
    c_labels = categorical({'FarUV Lamp','Solar UV'});
    EIrr = [EffectiveIrradiance_ofLAMP_mWcm2__new,EffectiveIrradiance_ofSOLAR_mWcm2__new];
    bar(c_labels,EIrr);
    ylabel('Effective Irradiance [mW_2_7_0_n_m / cm^2] ')

% % % % % % % % 
% % % % % % % % %show interplated values
% % % % % % % % figure; plot(TLV.nm,TLV.lamp.relPow,''); ylabel('Lamp Spectral Density [relative]');
% % % % % % % % hold all; plot(TLV.nm,TLV.lampwfilt.relPow);  %ylim([0,1]);
% % % % % % % % hold all; plot(1,1);
% % % % % % % % yyaxis right;
% % % % % % % % hold all; plot(TLV.nm, TLV.S.old ); ylabel('Sensitivity rel.270nm');
% % % % % % % % hold all; plot(TLV.nm, TLV.S.new );  %ylim([0,500]);
% % % % % % % % title({ ['Input: Full Spectrum Dose = ',num2str(LampIrrad_fullspectrum_mWcm2),' mJ/cm^2'];[''];
% % % % % % % %         ['RESULTS - "Effective Irradiance" '];
% % % % % % % %         ['     *** Requirement:  Eff.Dose <= 3 mJ_2_7_0/cm^2  ***']; 
% % % % % % % %         ['Eff.Dose(old,Raw):   ',num2str(round(E_old_per_RawSpectrum,2)),' mJ_2_7_0/cm^2'];...
% % % % % % % %         ['Eff.Dose(old, Filt):   ',num2str(round(E_old_per_FiltSpectrum,2)),' mJ_2_7_0/cm^2 (i.e. ',...
% % % % % % % %                                 num2str(round(E_old_per_RawSpectrum/E_old_per_FiltSpectrum,1)),'x lower)'];...
% % % % % % % %         ['Eff.Dose(new,Raw):   ',num2str(round(E_new_per_RawSpectrum,2)),' mJ_2_7_0/cm^2'];...
% % % % % % % %         ['Eff.Dose(new, Filt):   ',num2str(round(E_new_per_FiltSpectrum,2)),' mJ_2_7_0/cm^2 (i.e. ',...
% % % % % % % %                                 num2str(round(E_new_per_RawSpectrum/E_new_per_FiltSpectrum,1)),'x lower)']},'HorizontalAlignment','left');
% % % % % % % % xlim([200,lam_upper_lim]); xlabel('wavelength [nm]');
% % % % % % % % legend('Raw Lamp','Filt Lamp','','TLV(old)','TLV(new)')
% % % % % % % % yyaxis left;  yaxis_limits = ylim();
% % % % % % % % 
% % % % % % % % 
% % % % % % % % %% compute Weighted Value of TLV (from: 200-300nm)
% % % % % % % % E_eff_Limit = 3;
% % % % % % % % peaksum = sum(TLV.lamp.relPow(ikeep));  
% % % % % % % % TLV.lamp.OldLimit_mJcm2         =  E_eff_Limit / ( sum(TLV.S.old(ikeep) .* TLV.lamp.relPow(ikeep) )/peaksum );
% % % % % % % % TLV.lamp.NewLimit_mJcm2         =  E_eff_Limit / ( sum(TLV.S.new(ikeep) .* TLV.lamp.relPow(ikeep) )/peaksum );
% % % % % % % % peaksum = sum(TLV.lampwfilt.relPow(ikeep));  
% % % % % % % % TLV.lampwfilt.OldLimit_mJcm2    =  E_eff_Limit / ( sum(TLV.S.old(ikeep) .* TLV.lampwfilt.relPow(ikeep) )/peaksum );
% % % % % % % % TLV.lampwfilt.NewLimit_mJcm2    =  E_eff_Limit / ( sum(TLV.S.new(ikeep) .* TLV.lampwfilt.relPow(ikeep) )/peaksum );
% % % % % % % % 
% % % % % % % % yyaxis left; 
% % % % % % % % text(260,0.4*yaxis_limits(2),{'Effective TLVs: '; ...
% % % % % % % %     ['[Old, RawLamp ]  = ',num2str(round(TLV.lamp.OldLimit_mJcm2,1)),'mJ/cm2'];...
% % % % % % % %     ['[Old, Filtered]       = ',num2str(round(TLV.lampwfilt.OldLimit_mJcm2,1)),'mJ/cm2'];...
% % % % % % % %     ['[New, RawLamp ] = ',num2str(round(TLV.lamp.NewLimit_mJcm2,1)),'mJ/cm2'];...
% % % % % % % %     ['[New, Filtered]     = ',num2str(round(TLV.lampwfilt.NewLimit_mJcm2,1)),'mJ/cm2'];...
% % % % % % % %     });

end
 