%  this function reads 2D irradiance data from Zemax Text File
% input:  string:  "path/filename"
% output1:  DAT = 2D matrix of numbers
% output2:  meta = struct containing meta data taken from header

function [DAT,meta] = read_Zemax_TXT(filename)
            
    fID = fopen(filename);
    if fID~=-1
            %--------------------------------------------------------------

            %** Z *** ---get Z value for this file
            Nskip = 33;
            frewind(fID);
            for kk=1:Nskip
                tline1 = fgetl(fID);
            end
              %disp(tline1);
            tline2 = tline1(2:2:end); %remove alternating spaces from strong
            sline = split(tline2,':');
            meta.Z_cm = str2num(sline{2});

            %**L_x,L_y*** --- get L_x and L_y value for this file
            Nskip = 17;
            frewind(fID);
            for kk=1:Nskip
                tline1 = fgetl(fID);
            end
                    %disp(tline1)
            tline2 = tline1(2:2:end); %remove alternating spaces from strong
            sline = split(tline2,' ');
            meta.Lx_cm  = str2num(sline{2}); 
            meta.Ly_cm  = str2num(sline{5}); 
            
            
            %**TOTAL POWER*** --- get Total Power for this file
            Nskip = 23;
            frewind(fID);
            for kk=1:Nskip
                tline1 = fgetl(fID);
            end
                        %disp(tline1)
            tline2 = tline1(2:2:end); %remove alternating spaces from string
            sline = split(tline2,{':',' '});
            meta.Ptot_W  = str2num(sline{9});   
                       
            %--------------------------------------------------------------
            
            
            %*** READ 2D Data ****  
                % % %             Nskip = 49;  %skip to this line
                % % %             frewind(fID);
                % % %             for kk=1:Nskip
                % % %                 tline1 = fgetl(fID);
                % % %             end
                % % %             disp(tline1);
                % %             
            frewind(fID);
            tline1=fgetl(fID);
            tline2=tline1(2:2:end);
            n=0;
            while n<1000 & ~(length(strfind(tline2,'Units'))~=0 & length(strfind(tline2,'Watts')) ~= 0)  
             	tline1=fgetl(fID);
                tline2=tline1(2:2:end);
                n=n+1;
            end
            Nskip=6;    %SKIP 6 LINES
            for kk=1:Nskip	 tline1=fgetl(fID); end
            
            %disp(tline1);
            
            DAT = [];  
            while feof(fID)==0  %while NOT end of file
                %remove spaces from this line
                tline2 = tline1(2:2:end); %remove alternating spaces from string
                %convert text to numeric & store
                if length(tline2)>0
                    dat_thisline = textscan(tline2,'%f');
                    if length(DAT)==0
                        DAT = dat_thisline{1}';
                    else
                        DAT = [DAT;dat_thisline{1}'];
                    end
                end
                %read next line 
                tline1 = fgetl(fID);   
            end
            DAT = DAT(:,2:end);  %trim first column off  (not data, just indices)

            %--------------------------------------------------------------
            
            fclose(fID); 
    else
        warning(['File not found:  (filename= ',filename])
    end%(if) 
                    
end