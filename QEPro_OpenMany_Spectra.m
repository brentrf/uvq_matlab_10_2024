% Open QEPro Data Files
% export single matrix CSV
%
function QEPro_OpenMany_Spectra()
%% get data
    close all;
    clear;
    dirs.f  = 'G:\My Drive\Analyses(BF)\Matlab\UVQ';
    start_path = 'G:\My Drive\Analyses(BF)\Matlab\matDat\';
    start_path = 'G:\Shared drives\Corp Main\Engineering\LAB EXPERIMENTS\';
%     filepath = [uigetdir(start_path),'/'];

    cd(start_path)
    [filenames,path] = uigetfile('*.txt','Multiselect','on');
     Nf = length(filenames)
     if isstr(filenames)
         filename = filenames; clear filenames;
         filenames{1} = filename;
         Nf=1;
     end


% Loop Over Files
    for kf=1:Nf
        DAT  =  readQEPro_SpectrumFile([path,filenames{kf}]);
        M(:,kf) = DAT.A(:,2); %copy spectrum to matrix

        Tint(kf,1) = DAT.Tint;
%         if Nf <= 20
%             figure; plot(DAT.A(:,1),DAT.A(:,2)); hold all; plot(DAT.A(isel,1),DAT.A(isel,2));
%             title({replace(filenames{kf},'_',' '); ['sig = ',num2str(signal(kf,1))]});
%         end
    end

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
