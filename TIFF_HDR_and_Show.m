% Read in a stack of Images and generate HDR
%

function TIFF_HDR_and_Show()

    clear;
    close all;
    dirs.f = 'G:\My Drive\Analyses(BF)\Matlab\UVQ';
    
    %........ DIRECTORIES where relevant data is 
%      dirs.input = '/Users/brentfisher/Documents/MATLAB/UVQ/matdat';
     dirs.input = 'G:\Shared drives\Corp Main\Engineering\LAB EXPERIMENTS\2024_07_30_lab\2024_07_31_ALN2B(round1)_TopDown_BlueStreaks';
     
     dirs.input = [dirs.input, ''];
     cd(dirs.input);
     
    SatLimit_DN = 255;
    
     
    %select files to plot
    [IMGfile_list,IMGpath,dummy] = uigetfile({'*.tiff'},'Select DATA File to ','MultiSelect','on');
    if isstr(IMGfile_list)
        Nfiles = 1; tmp = IMGfile_list;
        IMGfile_list{1} = tmp;
    else
        Nfiles = length(IMGfile_list);
    end
    
    %SORT: lowest to highest integration time & get integration times
    for kf=1:Nfiles
        IMGfile = IMGfile_list{kf};
        I0 = imread([IMGpath,'/',IMGfile]);
        I1 = rgb2gray(I0(:,:,1:3));
        isel = find(1<I1 & I1<SatLimit_DN);
        Isum(kf,1) = sum(I1(isel));
        %get exposure time of each file
        ss1 = split(IMGfile,{'Te=','Texp=','Tint=','Tint'});
        ss2 = split(ss1{2},{'ms','us'});
        Texp_list(kf,1) = str2num(ss2{1});
    end
    filename_base = ss1{1};
    
    %sort order of intensity
    [sortedsums,file_order] = sort(Isum);
    [Texp_sortedlist,file_order_Texp] = sort(Texp_list);
    
    fhdr = figure;
    %LOOP & Build HDR
    for kf=1:Nfiles
        ksel = file_order_Texp(kf);  %which file to pull
        IMGfile = IMGfile_list{ksel}; 
        I0 = imread([IMGpath,'/',IMGfile]); 
        tmp = replace(IMGfile(1:end-5),'temp_profiles',' ');
        file_names{kf} = replace(tmp,'_',' '); 
        
        
        %get scale_factor ( for this file )
        Texp = Texp_list(ksel);
        imgsum =  sortedsums(kf);
        if Texp == max(Texp_list) 
            min_imgsum = min( Isum(find(Texp_list==max(Texp_list))) );  %find all images taken with maximum Texp
            scale_factors(kf,1) = min_imgsum / imgsum;                  %scale by total sum intensity of images that have "maximum Texp"
        else
            scale_factors(kf,1) = max(Texp_list) / Texp;
        end
        
        % process images into HDR
        if kf==1
            [Nr,Nc,Ncolors] = size(I0);
            I_HDR = zeros(Nr,Nc);
        end
        I1 = double( rgb2gray(I0(:,:,1:3)) );
        isel = find( I1 < SatLimit_DN );        %which pixels are less 
        I_HDR(isel) = I1(isel)*scale_factors(kf,1);
        
        %show progress
        figure(fhdr); imagesc(log10(I_HDR)); colorbar;
        dispstring = sprintf(['%8f   -- ',file_names{kf}], scale_factors(kf,1));
        disp(dispstring)
        pause(0.2);
    end 
    
    %PLOTS
    figure; surf(log10(I_HDR)); shading flat;
    
    [maxval,ipk] = max(I_HDR(:));
    [rpk,cpk] = ind2sub(size(I_HDR),ipk);
    figure; semilogy([1:Nr]-rpk,mean(I_HDR(:,cpk-5:cpk+5),2));
     hold all; semilogy([1:Nc]-cpk,mean(I_HDR(rpk-5:rpk+5,:),1));    
     grid on;  xlabel('pixels');
    
    % SAVE FILES
    save([filename_base,'_HDR.mat'],'I_HDR','scale_factors','file_names',...
                'Texp_list','file_order_Texp','Texp_sortedlist', 'Isum', 'file_order', 'sortedsums');
    csvwrite([filename_base,'_HDR.csv'],I_HDR);
    
    imwrite(uint16(I_HDR),[filename_base,'_HDR.tiff'],'tiff')

    
end