% This code reads the Argo Scripp files, grids and saves in Matlab format
% developed by March Jacob, Rev1: 11/16/2021
% The function takes as input 2004-2018 RG Argo Salinity Climatology
% update: 12/15/2021
% changing lon to 0-360 (originally is 20-380)
% and lat is flipped to match rim data (in RIM the index 1 is 90, and the
% 720 is -90. In Argo, originally, index 1 is -65 and 720 is 80)
% update: 12/20/2021
% making lat goes to -90 to 90, instead of original -65 to 80
% update: 03/12/2022
% to run on server or pc
% also, corrected a typo on latitude

function Reading_Gridding_Argo_Rev3(flag)

if flag==0
    inDir = '/data4/OceanSalinity/RIM/ARGO'; %input folder
    outDir = '/data4/OceanSalinity/RIM/ARGO/Data/GriddedData'; % output folder
elseif flag==1
    inDir = 'Z:\OceanSalinity\RIM\ARGO'; %input folder
    outDir = 'Z:\OceanSalinity\RIM\ARGO\Data\GriddedData'; % output folder
end

    res = 0.25; %desired resolution

    r = round(180/res); %rows in gridding matrix
    c = round(360/res); %columns in gridding matrix

    dirin = [inDir '/Data/RawData/']; % input folder
    files = dir([dirin,'RG_ArgoClim_S*']); %files to process

    name = [dirin,files.name]; %name of the file that is being processed
    disp(name)

    ncdisp(name)           

%     lati = double(ncread(name,'/LATITUDE')); %latitude vector
    loni = double(ncread(name,'/LONGITUDE')); %longitude vector
    time =  double(ncread(name,'/TIME')); %time vector, months since 2004-01-01 00:00:00
    sal = double(ncread(name,'/ARGO_SALINITY_MEAN')); %Argo salinity, pss-78 (-999 missing value), Size: 360x145x58 (lon,lat,pressure)
    anom = double(ncread(name,'/ARGO_SALINITY_ANOMALY')); %Argo salinity anomaly, pss-78 (-999 missing value), Size: 360x145x58x180 (lon,lat,pressure,time)
    sal = squeeze(sal(:,:,1));
    anom = squeeze(anom(:,:,1,:));    

    lati = -89.5:1:89.5; %original lat is from -65 to 80 (keeping the same format as original data)
    lati = lati';
    lati = flipud(repmat(lati,1,360)); %copy lat column for all lon and flip so it goes like IMERG
    loni = loni'; %convert long column to row
    loni = repmat(loni,180,1); %copy long rows to all lat
    anomi = permute(anom,[2 1 3]); %permute lat with long
    sali = permute(sal,[2 1 3]);
    
    anmly = [nan(25,360,180);anomi;nan(10,360,180)]; %expanding matrices to lat -90 to 90
    salty = [nan(25,360);sali;nan(10,360)]; %expanding matrices to lat -90 to 90
       
    lat = imresize(lati,[r c],'bicubic'); %resizing to matrix that is rows x cols
    lon = imresize(loni,[r c],'bicubic');
    lon(lon > 360) = lon(lon > 360) - 360; %from 20-380 to 0-360
    lon = [lon(:,1361:end) lon(:,1:1360)]; %so the first col is 0 lon
    anom = imresize(anmly,[r c],'bicubic');
    sal = imresize(salty,[r c],'bicubic');
       
    mydates = datetime(2004,1,1):calmonths(1):datetime(2018,12,1); %creating a vector with months and years for time
    m = month(mydates); %reading the month
    y = year(mydates);  %reading the yeat
    mm = m'; %fliping to col
    yy = y'; %fliping to col
    times = [time mm yy]; %adding month # and year # to time
    time = times(135:end,:); %saving just march/2015 till the end
    
    salinity = flipud(sal + anom); %calculating Argo salinity
    salinity = [salinity(:,1361:end,:) salinity(:,1:1360,:)]; %so it matches the lon
    salinity = salinity(:,:,135:end); %saving just march/2015 till the end

    name_file = ['Argo_RG_Salinity_2015-2018']; %setting name of file
    outputfile = [outDir,'/' name_file ,'.mat']; %setting path for saving
    save(outputfile,'lat','lon','salinity','time'); %saving output file
    disp('done')

end
