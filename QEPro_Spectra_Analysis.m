% Open QEPro Data Files
%
%
function MAIN()
%% get data
    close all;
    clear;
    dirs.f  = 'G:\My Drive\Analyses(BF)\Matlab\UVQ';
    start_path = 'G:\My Drive\Analyses(BF)\Matlab\matDat\';
    filepath = 'G:\My Drive\Analyses(BF)\Matlab\matDat\';
    filepath = 'G:\Shared drives\Corp Main\Engineering\LAB EXPERIMENTS\2024_09_23_LabData\test\'
%     filepath = [uigetdir(start_path),'/'];

    cd(start_path)
    [filenames,path] = uigetfile('*.txt','Multiselect','on');
     Nf = length(filenames);
     if isstr(filenames)
         filename = filenames; clear filenames;
         filenames{1} = filename;
         Nf=1;
     end


% Loop Over Files
    for kf=1:Nf

        DAT  =  readQEPro_SpectrumFile([path,filenames{kf}]);
        M(:,kf) = DAT.A(:,2); %copy spectrum to matrix
        BG(kf,1) = mean( DAT.A(end-10:end,2) );
        BGsigma = std(DAT.A(end-10:end,2) );
        BGrange = [DAT.A(end-10,1),DAT.A(end,1)];

        thresh = 200;   %BG(kf,1) + 50*BGsigma
        isel = find(DAT.A(:,2) > thresh);
        signal(kf,1) = sum(DAT.A(isel,2));
        numbins(kf,1) = length(isel);
        Tint(kf,1) = DAT.Tint;
        if Nf <= 20
            figure; plot(DAT.A(:,1),DAT.A(:,2)); hold all; plot(DAT.A(isel,1),DAT.A(isel,2));
            title({replace(filenames{kf},'_',' '); ['sig = ',num2str(signal(kf,1))]});
        end
    end
    tmp = split(filenames{kf},'QEP064681');
    titlestring = replace(tmp{1},'_',' ');
    figure; imagesc(M);  colorbar;   title({titlestring;'linear scale'});
    colorbar;    ylabel('wavelength[nm]');  xlabel('spectrum count');
    figure;  imagesc([1:size(M,2)],DAT.A(:,1),abs(log10(M))); title({titlestring;'log10 scale'});
    colorbar;    ylabel('wavelength[nm]');  xlabel('spectrum count');

    Mavg  = mean(M,2);
    figure; plot(DAT.A(:,1),Mavg);  xlabel('wavelength [nm]');  title({titlestring;'linear scale'});
    grid on;  ylim([50,150]*15); xlim([250,400]);



    csvwrite([path,'temp.csv'],[signal,numbins,BG,Tint]);
    csvwrite([path,'allspectra.csv'],[DAT.A(:,1),M]);
    disp(filenames');

end


 %% Load QE Pro File

function DAT  =  readQEPro_SpectrumFile(filepath_and_name)
    fID = fopen(filepath_and_name);
    clear DAT;
    if fID~=-1

        %Skip Header
        Nskip = 14;
        frewind(fID);     	for kk=1:Nskip
            tline1 = fgetl(fID);
            [newstrings,matches] = split(tline1,'(sec):');
            if length(matches)>0
                DAT.Tint = str2num(newstrings{2});
            end
        end

        %read spectral data
         kk=1;
        while feof(fID)==0  %while NOT end of file
            tline1 = fgetl(fID);

            %[1] split values by delimiter
            %               [vals,matches] = split(tline1,'\t');  %DOESNT WORK
            %[2] split values by 8 characters for each %DOESNT WORK
            %             val1 = num2str(tline1(1:7));
            %             val2 = num2str(tline1(9:end));
            %[3] split using "sscanf"
            vals = sscanf(tline1,'%f') ; %works
            wl =    vals(1);
            counts =  vals(2);
            DAT.A(kk,1) =   wl;
            DAT.A(kk,2) =  counts;
            kk=kk+1;
        end

%         figure; plot(DAT(:,1),DAT(:,2));




    end %if fid~=1
end
