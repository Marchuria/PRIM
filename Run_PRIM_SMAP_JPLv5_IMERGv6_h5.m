%-------------------------------------------------------------------------%
% This script was developed by Maria Jacob, Kyla Drushka, Bill Asher,     %
% Linwood Jones, Andrea Santos-Garcia                                     %
% Additional contact information: march.jacob@gmail.com &                 %
% kdrushka@apl.uw.edu                                                     %
%-------------------------------------------------------------------------%
% PRIM for SMAP                                                           %
% Inputs:                                                                 %
% IMERG v6B Precipitation                                                 %
% ERA5 wind speed                                                         %
% SMAP JPL v5: lat, lon, time, sss and                                    %
% RG Argo Salinity                                                        %
% It calculates:                                                          %
% PRIM_S0m: PRIM salinity estimate at 0m depth                            %
% PRIM_S1m: PRIM salinity estimate at 1m depth                            %
% PRIM_S5m: PRIM salinity estimate at 5m depth                            %
% Kz: vertical diffusivity coefficient                                    %
% PSS: Probability of salinity stratification between 0 and 10 m          %
%-------------------------------------------------------------------------%

function Run_PRIM_SMAP

% yyyy = integer input for the year in 4 digits 
% mm = vector input for the months% dd = vector input for the total of day per month
% flag = 0 means running in linux, 1 in windows

    yyyy = 2021; %year to run PRIM
    mm = 1:12;%months to run PRIM
    dd = 1:31; %days to run PRIM

    PRIM_SMAP_setting(yyyy,mm,dd,flag); %call to function to read the paths to data
    
end

function PRIM_SMAP_setting(yyyy,mm,dd,flag)

    if flag==0
        diskname = '/disk/'; %main disk name in linux
    elseif flag==1
        diskname = 'C:\'; %windows computer disk name
        addpath('Z:/OceanSalinity/RIM/kz_fromAPLUW'); 
    end
 
    path_smap = [diskname 'Data/path_to_gridded_smap_data']; %path to SMAP 
    %gridded data, generated with Reading_Gridding_SMAP_JPLV5.m
    %Note: we'll need to change this if we change the code name

    path_imerg = [diskname 'Data/path_to_gridded_daily_imerg_data']; 
    %path to imerg data, gridded data and concatenated daily, out of 
    %Concatenating_subdaily_IMERG_V06B.m
    %Note: we'll need to change this if we change the code name

    path_argo = [diskname 'Data/path_to_argo_gridded_data']; %path to argo 
    %data, generated with Reading_Gridding_Argo_until2018.m if data is 
    %previous to 2018 or Concatenating_Argo_after2018.m if data is after 2018
    %Note: we'll need to change this if we create one Argo file
        
    out_path = [diskname 'Data/path_to_output']; %location of the output 
    %path    

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Reading Argo data
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    argofile = [path_argo '/Argo_RG_Salinity_2021.mat'];
%     argofile = [path_argo '/Argo_RG_Salinity_2015-2018.mat'];
    load(argofile);
    
%% For loop design for running for multiple days and months
    for m = 1:size(mm,2)

        argo_month = time(:,2) == mm(1,m) & time(:,3) == yyyy;
    
        argo = salinity(:,:,argo_month);
        
        for n = 1:size(dd,2)  
            PRIM_SMAP_JPLv5_IMERGv6B(yyyy,mm(1,m),dd(1,n),path_imerg,path_smap,out_path,argo,flag,test_name);
%             PRIM_SMAP_JPLv5_IMERGv6B(yyyy,mm(1,m),dd(1,n),path_IMERG_Mat,pathSMAP,outputPath,flag,test_name);
        end
    end

end

function PRIM_SMAP_JPLv5_IMERGv6B(year,month,day,path_IMERG_Mat,pathSMAP,outputPath,argo,flag,test_name)
% function PRIM_SMAP_JPLv5_IMERGv6B(year,month,day,path_IMERG_Mat,pathSMAP,outputPath,flag,test_name)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Reading IMERG data
% The algorith needs to have available rain data from the previous day and
% the whole current day (48 hours total) to cover all the posible
% calculations
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Reading IMERG from Previous day
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    [year_1,month_1,day_1] = datevec(datenum(year,month,day)-1);

    fileRain1 = ['IMERG_V06B_30min_quarter_'  num2str(year_1) num2str(month_1,'%0.2d') num2str(day_1,'%0.2d')];

    try
        load([path_IMERG_Mat '/' fileRain1])
    catch
        fprintf(1,'%s: no file, skipping this date\n',[path_IMERG_Mat '/' fileRain1])
        return
    end

    prevIMERG = imerg;

    clear imerg

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Reading IMERG from Current day
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    fileRain2 = ['IMERG_V06B_30min_quarter_'  num2str(year) num2str(month,'%0.2d') num2str(day,'%0.2d') ];

    try
        load([path_IMERG_Mat '/' fileRain2])
    catch
        fprintf(1,'%s: no file, skipping this date\n',[path_IMERG_Mat '/' fileRain2])
        return
    end

    IMERG = imerg;

    precipitation_new = cat(3,prevIMERG,IMERG);

    clear IMERG imerg prevIMERG
	
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Calling function to generate WS from ERA5
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

wsera5 = era5ForRim(year_1,month_1,day_1,year,month,day,flag);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Depth Parameters
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    zsurf = 0.005; %This is the depth used to calculate RIM
    z1 = 1;
    z5 = 5;
    z10 = 10;   %This is depth used to calculate a reference at 10 meters to later calculate the probability of stratification per orbit

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Antenna parameters
% Weighting parameters to calculate the average over the SMAP foot print
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    w1 = 1;
    w2 = 1;
    wVec = [w2;w2;w1;w2;w2];
    ww = w2*4+w1;
    wMat = repmat(wVec,1,49);
    wMat2 = repmat(wVec,1,25);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Time vector for IMERG
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    time = repmat((24*3600:-1800:1800),5,1); %This variable is in seconds

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Creating list of files to process
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    s = dir([pathSMAP,'/SMAP_L2B_SSS_JPL_v5_',num2str(year,'%0.4d'),num2str(month,'%0.2d'),num2str(day,'%0.2d'),'_*.mat']);
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Running for the list of AQ files of the day
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

for rev = 2%x1:size(s,1)
    tic
    if s(rev).isdir==0
        
        %%%% Reading SMAP data %%%%
        gridfile=[pathSMAP '/' s(rev).name];
        load(gridfile); % SMAP file with grid indexes
        
        %%%% Creating sample time from IMERG %%%%
        sample_time = (0:0.5:23.5); % 0-24 hr every 30 min
				
		%%%% Creating sample time from ERA5 %%%%
		sample_timeWs = (0:1:23.5); % 0-24 hr every 30 min
		
        %%%% Creating and initializing new variables %%%%
        
        [rows,cols] = size(data.sss);
        
        kk = [1];% 25 37 43]; %43 (3h prior) - 37 (6 h prior) - 25 (12 h prior) - 13 (18 h  prior)
        
        data.PRIM_S0 = nan(rows,cols,length(kk),'single');
%         data.RIM_ins = nan(rows,cols,length(kk),'single');
%         data.RIM_prev = nan(rows,cols,length(kk),'single');
        data.PRIM_S1m = nan(rows,cols,length(kk),'single');
        data.PRIM_S5m = nan(rows,cols,length(kk),'single');
        data.PRIM_S10m = nan(rows,cols,length(kk),'single');
        
        data.S_ref = nan(rows,cols,length(kk),'single');
                
        data.rain_rate = nan(rows,cols,48,'single'); 
        data.rain_rate_inst = nan(rows,cols,'single'); 
        
        data.PSS = nan(rows,cols,'single');
%         data.bf_IRR = nan(rows,cols,1,'single');
%         data.bf_RA = nan(rows,cols,48,'single');
        
%         data.RIF = nan(rows,cols,48,'single');
%         data.RIFins = nan(rows,cols,'single');
		
% 		data.d0_w = nan(rows,cols,49,'single');
        data.Kz = nan(rows,cols,49,'single');
        data.wind_speed = nan(rows,cols,48,'single');
        
        flag2 = zeros(rows,cols,'int8');
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %% Main process
		%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        x = data.lat_ind; %extracting lat indexes from 0.25 degree resolution
        y = data.lon_ind; %extracting lon indexes from 0.25 degree resolution
        
        [rowInd,colInd] = find(x>=(1+1) & x<=(720-1) & y>=(1+1) & y<=(1440-1));% Finds indexes to be able to do the 5 pixels average
        
        flag1 = ~(x>=(1+1) & x<=(720-1) & y>=(1+1) & y<=(1440-1)); %quality flag that show points where IMERG is not available > we should take this, IMERG is always
            
        for m = 1:size(rowInd,1)
            
            a = (rowInd(m)); %lat ind
            b = (colInd(m)); %lon ind
            aa = (x(a,b)); %lat
            bb = (y(a,b)); %lon

            smap_time = datevec(data.time(a,b)); %SMAP time associated to the sample
            smap_t = smap_time(4) + (smap_time(5)/60)+ (smap_time(6)/3600); %Time in hours
%             data.secs = smap_time(4)*3600 + (smap_time(5))*60+ smap_time(6);

            t_dif = (sample_time - smap_t); % Difference between sample time and SMAP time (hours)
            [pick,id] = (max(t_dif(t_dif <= 0))); %To find the closest sample to SMAP time
%             [pick,id] = (min(abs(t_dif))); %test 2
			
			t_difWs = (sample_timeWs - smap_t); % Difference between sample time and SMAP time (hours) 
            [pickWs,idWs] = (max(t_difWs(t_difWs <= 0)));%To find the closest sample to SMAP time (hours)
%             [pickWs,idWs] = (min(abs(t_difWs))); %test 2
            
            cc = 1800*ones(size(sample_time)); %time vector for calculating RA (seconds)
            cc = repmat(cc,5,1)/(3600); %coefficient matrix in hours
            
%             hyc = data.hyc(a,b);
% 
%             if hyc < 0
%                 % no hycom data, so skip this point
%                 continue
%             end

%             argos = nanmean([argo(aa-1,bb,1);argo(aa,bb-1,1);...
%                 argo(aa,bb,1);argo(aa,bb+1,1); ...
%                 argo(aa+1,bb,1)]);
            argos = argo(aa,bb);
            argos(isnan(argos)==1) = -999;
            
            if argos < 0
                % no argo data, so skip this point
                continue
            end
            
            % CALCULATION OF RAIN HISTORY 
            % This process is done using the 5 boxes weighted averaged
            
%             m
            
            precip = reshape([precipitation_new(aa-1,bb,id:id+48);precipitation_new(aa,bb-1,id:id+48);...
                precipitation_new(aa,bb,id:id+48);precipitation_new(aa,bb+1,id:id+48); ...
                precipitation_new(aa+1,bb,id:id+48)],5,49);
            
%             if pick < 0
%                 timeIns = -pick*3600; %Time since IRR calculation, in secs
%             else
%                 timeIns = pick*3600; %Time since IRR calculation, in secs
%             end            
            
            timeIns = -pick*3600; %Time since IRR calculation, in secs
            indNan = ~isnan(precip); %finding non-nan positions of rain
            rain_w = nansum(precip.*wMat,1)./(sum(wMat.*indNan,1)); % weighted rain (not for RIM calculation, just for saving)
            
            %%%% calculating RA for kz %%%% 
            accum = precip(:,1:48).*cc./1000; %rain accumulation for 0.5 - 24 hs prior to satellite obs. in meters (mm/hr * hr * m/mm)
            accum(accum < 0) = 0;
			pick = 0.5;
            
            % pick is by construction negative, so we take abs
            accumins = abs(precip(:,49).*pick./1000); % instant accumulation in meters (mm/hr * hr * m/mm) 
            accum49 = [accum accumins]; % previous and instantaneous ccumulation in m
           
            indNan2=logical(sum(indNan,2)); %nan mask for rain
			
			% weighted instant accum, for storing
            accumins_w = nansum(accumins .* wMat(:,49))./sum(wMat(:,49).*indNan2);
            
            % --- if all rain data are nan, skip the remaining computation
            % for this pixel-
            if ~any(indNan2)
                continue
            end
            
            % --- if all rain data are nan *or zero*, fill the zero values
            % and skip the remaining computation
            if all(isnan(rain_w) | rain_w==0)
                data.rain_rate(a,b,:) = rain_w(1,1:48);
                data.rain_rate_inst(a,b) = rain_w(1,49);
                if isnan(data.rain_rate_inst(a,b))==1
                    flag2(a,b) = 1;
                end
                continue
            end
			
% 			% ***** if no rain in the most recent hour, skip this point ***
%             % we can comment this part out later, but it's useful for now
%             % to speed up the calculation since we mainly focus the
%             % statistics on times with instantaneous rain
%             if nansum(rain_w(end-2:end))==0
%                 1;
%                 continue
%             end
            
			%Calculation of the interpolated WS from GDAS
            winds = reshape([wsera5(aa-1,bb,idWs:idWs+24);wsera5(aa,bb-1,idWs:idWs+24);...
                wsera5(aa,bb,idWs:idWs+24);wsera5(aa,bb+1,idWs:idWs+24); ...
                wsera5(aa+1,bb,idWs:idWs+24)],5,25); % m/s
            
            indNanws = (winds<25 | winds>0);%finding bad positions
            ws_wtemp = sum(winds.*wMat2.*indNanws,1)./(sum(wMat2.*indNanws,1));
            %ws_wtemp(5)=data.scat_ws(a,b);% too unreliable to use!!!
            ws_wtemp(25) = data.NCEP_ws(a,b); % use ancillary wind speed as most recent
            ws_era5 = interp1([sample_timeWs 24],ws_wtemp,sample_time);%should have a 48 dimension
            % KD add: interpolate unweighted winds (t_ws) to 0:0.5:24
            % [sample_time 24] to get winds on the 13x49 grid of accum49
            wind_temp = winds;      
            wind_temp(:,25) = data.SMAP_ws(a,b); % instantaneous wind (radiometer)
            wind_for_kz = interp2([sample_timeWs 24],(1:5)',wind_temp,[sample_time 24],(1:5)');
			
			%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %% Vertical Difusivity
            %  ** Kz from kz_function_v2
			%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            
            kz = zeros(size(accum49));
            d0 = zeros(size(accum49));
			
			t_acum = (24:-0.5:0)*3600; %1800; %30 min of RA in sec
             
            for hp = 1:size(accum49,1)
                for hr = 1:size(accum49,2)
                    if isnan(wind_for_kz(hp,hr)+accum49(hp,hr)) || accum49(hp,hr)==0
                        % skip the kz calculation if we know it will be nan
                        kz(hp,hr)=nan;
                        d0(hp,hr)=nan;
                    else
                        % only one of kz_lookup or kz_function_v2 are used
                        % - leave one of them commented out for future versions.
                        
                        %[kz_ref, d0_ref] = kz_lookup(wind_for_kz(hp,hr), accum49(hp,hr));
                        [kz_ref, d0_ref] = kz_function_v2(wind_for_kz(hp,hr), accum49(hp,hr),t_acum(hr));   
                        kz(hp,hr) = kz_ref;
                        d0(hp,hr) = d0_ref;
                    end
                end                
            end

			% *** this is where we scale kz for testing, if needed ***
			% e.g., kz=kz*0.01;
            
            % *** Rev15: set d0 = 1 ***
            d0(~isnan(d0))=1; 

            % average kz and d0 over a grid box (not for the RIM
            % calculation, just for storing)
            % - Rev15: instead of using "indNan" (which is only the index
            % of nans for rain), use the index of nan values for kz and d0:
            indNan_kz=~isnan(kz);% finding non-nan positions of kz            
            kz_w = nansum(kz.*wMat,1)./(sum(wMat.*indNan_kz,1)); 
            indNan_d0=~isnan(d0);% finding non-nan positions of d0            
            d0_w = nansum(d0.*wMat,1)./(sum(wMat.*indNan_d0,1)); 
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %Calculating the Rain Impulse Functions for RIM
			%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            
            RIF49 = accum49.*d0; % multiplied by (integrated over) d0
            indNan_RIF = ~isnan(RIF49);% finding non-nan positions of RIF49
            RIF_w = nansum(RIF49.*wMat,1)./(sum(wMat.*indNan_RIF,1)); % weighted RIF
			
            
            for num = 1:length(kk) %remember adding 4 values to 3 dim in RIM values
                k = kk(num);
                
                %RIM @ the surface
                [rim,rimins,rimprev] = RIM_Calc(k,RIF49,kz,d0,zsurf,time,timeIns,wVec,argos);
%                 [rim,rimins,rimprev] = RIM_Calc(k,RIF49,kz,d0,zsurf,time,timeIns,wVec,hyc);
                
                if ~isnan(rim) && rim>-100 && nanmean(precip(:,49))>10
                    1;
                end
                
                %RIM @ 1m

                rim1 = RIM_Calc(k,RIF49,kz,d0,z1,time,timeIns,wVec,argos);
%                 rim1 = RIM_Calc(k,RIF49,kz,d0,z1,time,timeIns,wVec,hyc);
                
                %RIM @ 15m

                rim5 = RIM_Calc(k,RIF49,kz,d0,z5,time,timeIns,wVec,argos);
%                 rim5 = RIM_Calc(k,RIF49,kz,d0,z5,time,timeIns,wVec,hyc);

                %RIM @ 10m

                rim10 = RIM_Calc(k,RIF49,kz,d0,z10,time,timeIns,wVec,argos);
%                 rim10 = RIM_Calc(k,RIF49,kz,d0,z10,time,timeIns,wVec,hyc);

                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%                              
                %%% saving variables
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

                rain_w = single(rain_w);

                data.PRIM_S0(a,b) = rim;
%                 data.RIM_ins(a,b) = rimins;
%                 data.RIM_prev(a,b) = rimprev;
                data.PRIM_S1m(a,b) = rim1;
                data.PRIM_S5m(a,b) = rim5;
                data.PRIM_S10m(a,b) = rim10;
                
                data.S_ref(a,b) = argos;

                data.rain_rate(a,b,:)= rain_w(1,1:48) ;
                data.rain_rate_inst(a,b) = rain_w(1,49);
%                 data.RIF(a,b,:) = RIF_w(1,1:48);
%                 data.RIFins(a,b) = RIF_w(1,49);
% 				data.accumins(a,b) = accumins_w;
                data.wind_speed(a,b,:) = ws_era5;
                data.wind_speed_inst = data.SMAP_ws;
%                 d0_w = single(d0_w);
%                 data.d0_w(a,b,:) = d0_w(1,1:49);
                kz_w = single(kz_w);
                data.Kz(a,b,:) = kz_w(1,1:49);

                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                %calculating beam fraction BF
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%                 bf = nansum((precip>0.25).*wMat,1)./ww;
%                 data.bf_IRR(a,b) = bf(1,49);
%                 data.bf_RA(a,b,:) = bf(1,1:48);

                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                % Calculating flags
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

                %Idetifiying IRR bad values
                if isnan(data.rain_rate(a,b)) == 1
                flag2(a,b) = 1;
                end

                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                % Calculating the probability of stratification PS
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                deltaRim = data.PRIM_S10m - data.PRIM_S0;
                mindelta = min(min(deltaRim));
                maxdelta = max(max(deltaRim));
                data.PSS = (deltaRim - mindelta)*((1 - 0)/(maxdelta - mindelta))+0;
            end %k loop end
        end
        
        flagB1 = bi2de([flag1(:,:);flag2(:,:)]');
        flag = uint8(flagB1');
        fl = flag>0;
                     
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %Organizing and populating flags
		%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%        

        data.rain_rate_inst(fl) = -999;
        data.rain_rate_inst(isnan(data.rain_rate_inst)) = -999;
        data.rain_rate_inst(fl) = -999;
        data.rain_rate(isnan(data.rain_rate_inst)) = -999;
        
        data.PRIM_S0(fl) = nan;
        data.PRIM_S0(isnan(data.PRIM_S0)) = -999;
%         data.RIM_ins(fl) = nan;
%         data.RIM_ins(isnan(data.RIM_SSS)) = nan;
%         data.RIM_prev(fl) = nan;
%         data.RIM_prev(isnan(data.RIM_SSS)) = nan;
        data.PRIM_S1m(fl) = nan;
        data.PRIM_S1m(isnan(data.PRIM_S1m)) = -999;
        data.PRIM_S5m(fl) = nan;
        data.PRIM_S5m(isnan(data.PRIM_S5m)) = -999;
        
%         data.S_ref(fl) = nan;
%         data.S_ref(isnan(data.S_ref)) = -999;
        
        data.PSS(fl) = nan;
        data.PSS(isnan(data.PSS)) = -999;
%         data.bf_IRR(fl) = nan;
%         data.bf_IRR(isnan(data.bf_IRR)) = nan;
%         data.bf_RA(fl) = nan;
%         data.bf_RA(isnan(data.bf_RA)) = nan;
        
%         data.RIFins(fl) = nan;
%         data.RIFins(isnan(data.RIFins)) = nan;
%         data.RIF(fl) = nan;
%         data.RIF(isnan(data.RIF)) = nan;
		
% 		data.d0_w(fl) = nan;
% 		data.d0_w(isnan(data.d0_w)) = nan;
        data.Kz(fl) = nan;
		data.Kz(isnan(data.Kz)) = -999;

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %% Saving output per orbit
		%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        name_file = ['PRIM_' ,num2str(year), num2str(month,'%0.2d'),num2str(day,'%0.2d'),'_', num2str(rev,'%0.2d'),'_JPLv5_IMERGv6B'];
        outputfile = [outputPath,'/' name_file ,'.h5'];
        outputfilem = [outputPath,'/' name_file ,'.mat'];
        save(outputfilem,'data');

        hdf5write(outputfile,'SMAP_variables/time',data.time,...
            'SMAP_variables/lat',data.lat,'SMAP_variables/lon',data.lon,...
            'PRIM_Inputs/S_ref',data.S_ref,...
            'PRIM_Inputs/rain_rate_inst',data.rain_rate_inst,'PRIM_Inputs/rain_rate',data.rain_rate,...
            'PRIM_Inputs/wind_speed_inst',data.wind_speed_inst,'PRIM_Inputs/wind_speed',data.wind_speed,...            
            'PRIM_Outputs/PRIM_S0',data.PRIM_S0,'PRIM_Outputs/PRIM_S1m',data.PRIM_S1m,...
            'PRIM_Outputs/PRIM_S5m',data.PRIM_S5m, ...
            'PRIM_Outputs/PSS',data.PSS,'PRIM_Outputs/Kz',data.Kz);
        
        today = datevec(now);
        doy = juliandate(today(1:3)) - juliandate(today(1),0,0);
        proctime = [num2str(today(1),'%0.4d'),num2str(doy,'%0.3d'),num2str(today(4),'%0.2d'),num2str(today(5),'%0.2d'),num2str(floor(today(6)),'%0.2d'),num2str(round((today(6)-floor(today(6)))*1000),'%0.3d')];
        
        %% General Atributes
%         h5writeatt(outputfile, '/Ancillary_Info/','description','ancillary information');
        h5writeatt(outputfile, '/','Start Year',AllAttributes.attributes(1,2).Value);
        h5writeatt(outputfile, '/','Start Day of the Year',AllAttributes.attributes(1,3).Value);
        h5writeatt(outputfile, '/','Start Time',AllAttributes.attributes(1,6).Value.Data);
        h5writeatt(outputfile, '/','End Time',AllAttributes.attributes(1,7).Value.Data);
        h5writeatt(outputfile, '/','Orbit Number',AllAttributes.attributes(1,1).Value.Data);
        h5writeatt(outputfile, '/','Product Name',[name_file,'.h5']);
        h5writeatt(outputfile, '/','Processing Time',proctime);
        productInfo = 'TBD';
        h5writeatt(outputfile, '/','Product_Info',productInfo);
        h5writeatt(outputfile, '/','Precipitation Source','IMERG_V60B');
        
        %% SMAP_data Attributes
        h5writeatt(outputfile, '/SMAP_variables/','description','parameters from SMAP JPL dataset level 2 v5');        
        %lat
        h5writeatt(outputfile, '/SMAP_variables/lat','units','degrees');
        h5writeatt(outputfile, '/SMAP_variables/lat','long_name','latitude');
        %lon
        h5writeatt(outputfile, '/SMAP_variables/lon','units','degrees');
        h5writeatt(outputfile, '/SMAP_variables/lon','long_name','longitude');
        %sec
        h5writeatt(outputfile, '/SMAP_variables/time','units','seconds');
        h5writeatt(outputfile, '/SMAP_variables/time','long_name','second of day');        
%         %asc_index
%         h5writeatt(outputfile, '/SMAP_data/asc_index','units','no units');
%         h5writeatt(outputfile, '/SMAP_data/asc_index','long_name','ascending/descending index');
%         %sst
%         h5writeatt(outputfile, '/SMAP_data/sst','units','kelvin');
%         h5writeatt(outputfile, '/SMAP_data/sst','long_name','sea surface temperature');
        
        %% Ancillary_info attributes
        h5writeatt(outputfile, '/PRIM_Inputs/','description','ancillary datasets)');
        %anc_SSS
        h5writeatt(outputfile, '/PRIM_Inputs/S_ref','units','psu');
        h5writeatt(outputfile, '/PRIM_Inputs/S_ref','long_name','salinity from Argo');
        %% Rain_Product
        % IRR
        h5writeatt(outputfile, '/PRIM_Inputs/rain_rate_inst','units','millimeter/hour');
        h5writeatt(outputfile, '/PRIM_Inputs/rain_rate_inst','long_name','closest IMERG rain rate to SMAP time & collocated to the center of each SMAP IFOV');
        h5writeatt(outputfile, '/PRIM_Inputs/rain_rate_inst','Not calculated or invalid','-999');

        % RR
        h5writeatt(outputfile, '/PRIM_Inputs/rain_rate','units','millimeter/hour');
        h5writeatt(outputfile, '/PRIM_Inputs/rain_rate','long_name','IMERG Rain Rate for the previous 24 hs to SMAP time');
        h5writeatt(outputfile, '/PRIM_Inputs/rain_rate','Not calculated or invalid','-999');
        %SMAP_ws
        h5writeatt(outputfile, '/PRIM_Inputs/wind_speed_inst','units','meter/second');
        h5writeatt(outputfile, '/PRIM_Inputs/wind_speed_inst','long_name','SMAP radiometer wind speed');
        %anc_ws
        h5writeatt(outputfile, '/PRIM_Inputs/wind_speed','units','meter/second');
        h5writeatt(outputfile, '/PRIM_Inputs/wind_speed','long_name','ERA5 wind speed for the previous 24 hs to SMAP time');
%         % status_flag
%         h5writeatt(outputfile, '/Ancillary_Info/status_flag','long_name','status flag');
%         h5writeatt(outputfile, '/Ancillary_Info/status_flag','value','2 bits flag, first bit in 1  indicates CMORPH is not available, second bit in 1 indicates CMORPH for IRR has been remove for low quality');
        h5writeatt(outputfile, '/PRIM_Outputs/','description','rain impact parameters based on IMERG, ERA5 and PRIM (Parametrized Rain Impact Model)');
%         % BF_IRR
%         h5writeatt(outputfile, '/Rain_Impact/BF_IRR','units','%');
%         h5writeatt(outputfile, '/Rain_Impact/BF_IRR','long_name','Rain Beam Fill Fraction for the Instantaneous Rain Rate');
%         h5writeatt(outputfile, '/Rain_Impact/BF_IRR','Not calculated or invalid','-999');
% 
%         % BF_RA
%         h5writeatt(outputfile, '/Rain_Impact/BF_RA','units','%');
%         h5writeatt(outputfile, '/Rain_Impact/BF_RA','long_name','Rain Beam Fill Fraction for the precipitation history');
%         h5writeatt(outputfile, '/Rain_Impact/BF_RA','Not calculated or invalid','-999');

        % PSS
        h5writeatt(outputfile, '/PRIM_Outputs/PSS','valid_max','1');
        h5writeatt(outputfile, '/PRIM_Outputs/PSS','valid_min','0');
        h5writeatt(outputfile, '/PRIM_Outputs/PSS','long_name','Probability of Salinity Stratification');
        h5writeatt(outputfile, '/PRIM_Outputs/PSS','Not calculated or invalid','-999');

        % Kz
        h5writeatt(outputfile, '/PRIM_Outputs/Kz','units','m2/s');
        h5writeatt(outputfile, '/PRIM_Outputs/Kz','long_name','Vertical Diffussivity Coefficient');
        h5writeatt(outputfile, '/PRIM_Outputs/Kz','Not calculated or invalid','-999');

        % RIM
        h5writeatt(outputfile, '/PRIM_Outputs/PRIM_S0','units','psu');
        h5writeatt(outputfile, '/PRIM_Outputs/PRIM_S0','long_name','Sea Surface Salinity (SSS) estimated using PRIM');
        h5writeatt(outputfile, '/PRIM_Outputs/PRIM_S0','Not calculated or invalid','-999');
        
        h5writeatt(outputfile, '/PRIM_Outputs/PRIM_S1m','units','psu');
        h5writeatt(outputfile, '/PRIM_Outputs/PRIM_S1m','long_name','Salinity at 1 meter depth estimated using RIM');
        h5writeatt(outputfile, '/PRIM_Outputs/PRIM_S1m','Not calculated or invalid','-999');
        
%         h5writeatt(outputfile, '/Rain_Impact/RIM_3m','units','psu');
%         h5writeatt(outputfile, '/Rain_Impact/RIM_3m','long_name','Salinity at 3 meters depth estimated using RIM');
%         h5writeatt(outputfile, '/Rain_Impact/RIM_3m','Not calculated or invalid','-999');
%         
        h5writeatt(outputfile, '/PRIM_Outputs/PRIM_S5m','units','psu');
        h5writeatt(outputfile, '/PRIM_Outputs/PRIM_S5m','long_name','Salinity at 5 meters depth estimated using RIM');
        h5writeatt(outputfile, '/PRIM_Outputs/PRIM_S5m','Not calculated or invalid','-999');
        
%         h5writeatt(outputfile, '/Rain_Impact/RIM_anom','units','psu');
%         h5writeatt(outputfile, '/Rain_Impact/RIM_anom','long_name','Sea Surface Salinity (SSS) anomaly estimated using RIM');
%         h5writeatt(outputfile, '/Rain_Impact/RIM_anom','Not calculated or invalid','-999');
%         
%         h5writeatt(outputfile, '/Rain_Impact/RIM_anom1m','units','psu');
%         h5writeatt(outputfile, '/Rain_Impact/RIM_anom1m','long_name','Salinity anomaly at 1 meter depth estimated using RIM');
%         h5writeatt(outputfile, '/Rain_Impact/RIM_anom1m','Not calculated or invalid','-999');
%         
%         h5writeatt(outputfile, '/Rain_Impact/RIM_anom3m','units','psu');
%         h5writeatt(outputfile, '/Rain_Impact/RIM_anom3m','long_name','Salinity anomaly at 3 meters depth estimated using RIM');
%         h5writeatt(outputfile, '/Rain_Impact/RIM_anom3m','Not calculated or invalid','-999');
%         
%         h5writeatt(outputfile, '/Rain_Impact/RIM_anom5m','units','psu');
%         h5writeatt(outputfile, '/Rain_Impact/RIM_anom5m','long_name','Salinity anomaly at 5 meters depth estimated using RIM');
%         h5writeatt(outputfile, '/Rain_Impact/RIM_anom5m','Not calculated or invalid','-999');
             
        clearvars -except year month day path_IMERG_Mat pathSMAP path_argo outputPath rev s precipitation_new argo w1 w2 w3 w4 wVec ww wMat wMat2 time kz1 d01 zsurf z1 z5 z10 wsera5 test_name;
    end
    toc
end

end

function [rim,rimins,rimprev] = RIM_Calc(k,RIF49,kz,d0,z,time,timeIns,wVec,argos)

% k is the time index between 1 and 48 over which to compute RIM
%   k=1:48 corresponds to times of 24:-0.5:0.5 hours from present
%   (use k=1 for the past 24 hours)
% RIF49 is the rain impact function (size 5x49 : 5 corresponds to
%    elements within the footprint; 49 corresponds to time, with the 49th
%    value being instantaneous RIF)
% timeIns is the time in seconds between the last rain observation & SMAP
% wVec is the weighting vector for the satellite footpring (size 5)
% argo is hycom reference salinity (size 1)
% indNan2 is the index of non-nan rain values (size 5)
%
% The RIM calculation is 
% S = So * do * ( do + R/sqrt(Kz*t)*exp(-z^2/*4*Kz*t) )^-1
%   where So is the initial bulk salinity (hycom)
%   do is the mixing depth and is now 1
%   R (m^2) is the rain surface impulse function (RIF) = rain at time t0 integrated over do (?)
%   z = depth at with RIM is calculated, t is time in sec
%
% The way RIM works is that So is successively multiplied by the 
% ( do + R/sqrt(Kz*t)*exp(-z^2/*4*Kz*t) )^-1 term , which is computed at
% each time.

    % dilM is the dilution matrix for all previous times (i.e., everything but
    % instananeous)
    % i.e., ( do + R/sqrt(Kz*t)*exp(-z^2/*4*Kz*t) )

    dilM = (d0(:,k:48)+(((RIF49(:,k:48))./(sqrt(time(:,k:48).*kz(:,k:48)))).*exp((-z.^2)./(4*kz(:,k:48).*time(:,k:48)))));
    
    % next we compute the terms of instantaneous dilution:
    dilIns = (d0(:,49)+(((RIF49(:,49))./(sqrt(timeIns.*kz(:,49)))).*exp((-z.^2)./(4*kz(:,49).*timeIns))));
   
    % merge the previous (dilM) and instantanous (dilIns) terms into one matrix
    % and divide by d0 because later we take the inverse, so we are ultimately multiplying by d0
    dil_Hyst = [(dilM./d0(:,k:48)) (dilIns./d0(:,49))];
    % just the instaneous part:
    dil_Hystins = (dilIns./d0(:,49));
    % just the previous part
    dil_Hystprev = (dilM./d0(:,k:48));
    % set nan values to 1 to avoid problems taking the inverse, since 1 will
    % have no impact in the product (?? is this a good strategy?)
    dil_Hyst(isnan(dil_Hyst)) = 1;
    dil_Hystins(isnan(dil_Hystins)) = 1;
    dil_Hystprev(isnan(dil_Hystprev)) = 1;
    % take the inverse, to compute do * ( do + R/sqrt(Kz*t)*exp(-z^2/*4*Kz*t) )^-1
    % and apply the product operator (over multiple times)
    %   previous and instantaneous times:
    dil = (prod(dil_Hyst,2)).^-1; 
    %   just the instantaneous part (??? does the product do anything ???)
    dil2 = (prod(dil_Hystins,2)).^-1;
    %   just the previous part:
    dil1 = (prod(dil_Hystprev,2)).^-1;
    
    % use the nan indices for each term
    %   previous and instantaneous times:
    indNan_dil=logical(sum(~isnan(dil),2)); % finding non-nan positions of dil   
    dil = nansum(dil.*(wVec.*indNan_dil),1)./(sum(wVec.*indNan_dil,1));
    %   just the instantaneous part  
    indNan_dil2=logical(sum(~isnan(dil2),2)); % finding non-nan positions of dil    
    dil2 = nansum(dil2.*(wVec.*indNan_dil2),1)./(sum(wVec.*indNan_dil2,1));
    %   just the previous part:
    indNan_dil1=logical(sum(~isnan(dil1),2)); % finding non-nan positions of dil   
    dil1 = nansum(dil1.*(wVec.*indNan_dil1),1)./(sum(wVec.*indNan_dil1,1));
    % now, compute RIM SSS by multiplying by initial salinity
    % S = So * do * ( do + R/sqrt(Kz*t)*exp(-z^2/*4*Kz*t) )^-1
    %   previous and instantaneous times (this is RIM SSS)
    rim = argos.*dil;
    %   just the instantaneous part 
    rimins = argos.*dil2;
    %   just the previous part:
    rimprev = argos.*dil1;
    % convert to singles and return 
    rim = single(rim);
    rimins = single(rimins);
    rimprev = single(rimprev);    
    
end

function [wsera5] = era5ForRim(year_1,month_1,day_1,year,month,day,flag)
%Developed by Andrea Santos-Garcia at Central Florida Remote Sensing Lab
%updated by March Jacob, on 10/15/2020 to read ERA5 data
%This function

temp = num2str(year_1);
if flag==0
    dirname = ['/data4/OceanSalinity/RIM/ERA5/Data/GriddedData'];
elseif flag==1
    dirname = ['Z:/OceanSalinity/RIM/ERA5/Data/GriddedData'];
end


% - load first find file:
fname_1 = [dirname '/ERA5_' temp num2str(month_1,'%0.2d') num2str(day_1,'%0.2d') '.mat'];
% fname_1 = [dirname '/' temp num2str(month_1,'%0.2d') num2str(day_1,'%0.2d') '_surf' '.mat']; %Need to change name of saving ERA5
% wwarning if this wind file doesn't exist
if ~exist(fname_1,'file')
    fprintf(1,'** wind file %s doesnt exist **\n',fname_1);
end
load(fname_1,'wind');
wsera5(:,:,1:24) = flipud(wind);
temp = num2str(year);

% l- oad second wind file:
fname = [dirname '/ERA5_' temp num2str(month,'%0.2d') num2str(day,'%0.2d') '.mat'];
% fname = [dirname '/' temp num2str(month,'%0.2d') num2str(day,'%0.2d') '_surf' '.mat'];
% warning if this wind file doesn't exist
if ~exist(fname,'file')
    fprintf(1,'** wind file %s doest exist **\n',fname);
end
load(fname,'wind');
wsera5(:,:,25:48) = flipud(wind);

end % end function