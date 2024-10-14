% TLV_Spectral_Calculations

function TLV_Spectral_Calculations_v1()
clear; close all;

%INPUTs
dirs.d = '/Users/brentfisher/Documents/MATLAB/UVQ/matDat';
LampFile = '/KrCl_spec0.csv';       %spectrum: KrCl
% LampFile = '/KrBr_spec1.csv';     %spectrum: KrBr
LampDose_mJcm2_at_Peak = 218;        %assume this is dose of the (filtered) *Peak Only*

%%
%TLV Data ....
tmp = csvread([dirs.d,'/TLV.csv'],1,0);
TLV.nm          = tmp(:,1);
TLV.mWcm2.old   = tmp(:,2);
TLV.mWcm2.new   = tmp(:,3);
indx270            = find(TLV.nm == 270);
TLV.S.old       = tmp(indx270,2)./tmp(:,2);
TLV.S.new       = tmp(indx270,3)./tmp(:,3);
% figure; plot(TLV.nm,TLV.S.old); hold all; plot(TLV.nm,TLV.S.new); xlim([200,300]);

%LAMP Data .....
tmp = csvread([dirs.d,LampFile],1,0);
lamp.nm         = tmp(:,1);
lamp.relPow     = tmp(:,2);

%interpolate & return result
TLV.lamp.relPow  = interp1(tmp(:,1),tmp(:,2),TLV.nm);
disp(TLV.lamp.relPow);

%apply ideal filter (10nm wide around peak)
FilterFullWidth_nm = 15;
[maxval,imax] = max(lamp.relPow); lamp.peak_nm = lamp.nm(imax);
TLV.lampwfilt.relPow  = TLV.lamp.relPow ...
                     .* ( abs(TLV.nm - lamp.peak_nm) <= FilterFullWidth_nm/2 );  %apply filter window


%normalize
lamp.relPow = lamp.relPow / max(lamp.relPow);
TLV.lamp.relPow      = TLV.lamp.relPow / max(TLV.lamp.relPow);
TLV.lampwfilt.relPow = TLV.lampwfilt.relPow / max(TLV.lampwfilt.relPow);


%% compute Weighted Value of TLV (from: 200-300nm)
ikeep = find( (TLV.nm <=300) & ~isnan(TLV.lamp.relPow + TLV.lampwfilt.relPow +  TLV.mWcm2.old + TLV.mWcm2.new ));
% wrong...
% TLV_old_per_RawSpectrum     = sum(TLV.mWcm2.old(ikeep) .* TLV.lamp.relPow(ikeep)      ) / sum(TLV.lamp.relPow(ikeep));
% TLV_old_per_FiltSpectrum    = sum(TLV.mWcm2.old(ikeep) .* TLV.lampwfilt.relPow(ikeep) ) / sum(TLV.lampwfilt.relPow(ikeep));
% TLV_new_per_RawSpectrum     = sum(TLV.mWcm2.new(ikeep) .* TLV.lamp.relPow(ikeep)      ) / sum(TLV.lamp.relPow(ikeep));
% TLV_new_per_FiltSpectrum    = sum(TLV.mWcm2.new(ikeep) .* TLV.lampwfilt.relPow(ikeep) ) / sum(TLV.lampwfilt.relPow(ikeep));
% right...
peaksum = sum(TLV.lampwfilt.relPow(ikeep))  
E_old_per_RawSpectrum     = sum(TLV.S.old(ikeep) .* TLV.lamp.relPow(ikeep)      )*LampDose_mJcm2_at_Peak/peaksum ;
E_old_per_FiltSpectrum    = sum(TLV.S.old(ikeep) .* TLV.lampwfilt.relPow(ikeep) )*LampDose_mJcm2_at_Peak/peaksum ;
E_new_per_RawSpectrum     = sum(TLV.S.new(ikeep) .* TLV.lamp.relPow(ikeep)      )*LampDose_mJcm2_at_Peak/peaksum ;
E_new_per_FiltSpectrum    = sum(TLV.S.new(ikeep) .* TLV.lampwfilt.relPow(ikeep) )*LampDose_mJcm2_at_Peak/peaksum ;



%show interplated values
figure; plot(lamp.nm,lamp.relPow); 
hold all; plot(TLV.nm,TLV.lamp.relPow,'o');
hold all; plot(TLV.nm,TLV.lampwfilt.relPow,'x');  ylim([0,1])
yyaxis right;
hold all; plot(TLV.nm, TLV.mWcm2.old );
hold all; plot(TLV.nm, TLV.mWcm2.new );  ylim([0,500]);
title({ ['Dose of Peak (only): ',num2str(LampDose_mJcm2_at_Peak),' mJ/cm^2'];
        ['Eff.Dose(old,Raw):  ',num2str(E_old_per_RawSpectrum),' mJ_2_7_0/cm^2'];...
        ['Eff.Dose(old,Filt): ',num2str(E_old_per_FiltSpectrum),' mJ_2_7_0/cm^2'];...
        ['Eff.Dose(new,Raw):  ',num2str(E_new_per_RawSpectrum),' mJ_2_7_0/cm^2'];...
        ['Eff.Dose(new,Filt): ',num2str(E_new_per_FiltSpectrum),' mJ_2_7_0/cm^2']});
xlim([200,300]);
legend('Raw Lamp','Raw Lamp','Filt Lamp','TLV(old)','TLV(new)')


end
