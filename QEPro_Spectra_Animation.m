% 1. Open Spectra Files 
 
%
function QEPro_Spectra_Animation()
    close all;
    clear;
    program_name = 'image analysis';
    
    start_path = 'G:\Shared drives\Corp Main\Engineering\LAB EXPERIMENTS\2024_09_23_LabData\';
    
%% select SPECTRA STACK
   
    %Select by Hand
    cd(start_path);
    [filenames,path] = uigetfile('*.txt','Multiselect','on');
    Nf = length(filenames);
    if isstr(filenames)
        filename = filenames; clear filenames;
        filenames{1} = filename;
        Nf=1;
    end
%% Read All SPECTRA Files
    for kf = 1:Nf        
       DAT{kf} =  readQEPro_SpectrumFile([path,filenames{kf}])
    end
    
%% select BG SPECTRUM
    %Select by Hand
    cd(start_path);
    [BG_filename,BG_path] = uigetfile('*.txt','Multiselect','off');
    if BG_filename==0   BGflag = 0;
    else                BFflag = 1;
    end
    
    if isstr(BG_filename)   BG_filename = BG_filename;  
    end
%% Read BG SPECTRUM
    BG =  readQEPro_SpectrumFile([BG_path,BG_filename]);
    
    
s = split(filenames{1},'set') 
s2 = split(s{2},'__') 
WL = str2num(s2{1});

%% Generate Animation
%set up video
v = VideoWriter('test3.avi');
v.FrameRate = 10;
open(v);
hf = figure;
for kf = 1:Nf
    s = split(filenames{kf},'set');
    s2 = split(s{2},'__') ;
    pumpWL = str2num(s2{1});

    %figure update
    figure(hf);
    WL = DAT{kf}.A(:,1);
    Y = DAT{kf}.A(:,2) - BG.A(:,2);
    plot(WL,Y); xlim([280,340]); 
    title(['pump WL=',num2str(pumpWL),'nm   frame=',num2str(kf)]);
    ylim([0,40]);
    
    %generate video
    movie_F(kf) = getframe(gcf);
    writeVideo(v,movie_F(kf));
    pause(0.25);    
end

close(v);

movie(gcf,movie_F,10,1);


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