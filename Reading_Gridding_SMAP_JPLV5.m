%Reading JPL SMAP Data
% Developed by March at CFRSL - January 18 2016
% Updated on Nov 10th 2016 to read SMAP CAP v2.2_beta
% Updated on Nov 14th to read SMAP CAP v3
% Updated on Aug 28th 2018 to correct some typos >> time
% Updated on 11/28/2018 to save all flags needed
% Updated on 12/06/2018 to filter bad values (see comments below)
% Updated on 04/02/2021 to read JPL v4.3
% Updated on 07/20/2021 to read JPL v5
% and to use XCAl gridding formula
% Updated on 02/08/2022 for 10 km res
% Updated on 04/01/2022, so it deletes the processed file, to save space

%Reading variables in SMAP file

clearvars
clc

yyyy = 2018;

dirInput = '/disk4/OceanSalinity/RIM/SMAP/Data/RawData/JPL/v5';
% dirInput = 'Z:\OceanSalinity\RIM\SMAP\Data\RawData\JPL\v5';
dirOutput = '/disk4/OceanSalinity/RIM/SMAP/Data/GriddedData/JPL/v5/10km';
% dirOutput = 'Z:\OceanSalinity\RIM\SMAP\Data\GriddedData\JPL\v5\10km';
s = dir([dirInput '/SMAP_L2B_SSS*','_',num2str(yyyy,'%0.4d'),'*']);
files = ls([dirInput,'/SMAP*']);

for k = 1:size(s,1)
    
    name = [dirInput '/' s(k).name];
    hinfo = hdf5info(name);
    
    attributes = hinfo.GroupHierarchy.Attributes;      
    
    year = str2num(s(k).name(20:23)); % year from the file name
    month = str2num(s(k).name(24:25)); % month from the file name
    day = str2num(s(k).name(26:27)); % day from the file name
    revno = str2num(s(k).name(14:18)); % revno from the file name
    
    lat = double(hdf5read(hinfo.GroupHierarchy.Datasets(15)));
    lon = double(hdf5read(hinfo.GroupHierarchy.Datasets(16)));
    
    sss = double(hdf5read(hinfo.GroupHierarchy.Datasets(34)));
    hyc = double(hdf5read(hinfo.GroupHierarchy.Datasets(3)));
    sst = double(hdf5read(hinfo.GroupHierarchy.Datasets(4)));
    
    NCEP_ws = double(hdf5read(hinfo.GroupHierarchy.Datasets(2)));
    SMAP_ws = double(hdf5read(hinfo.GroupHierarchy.Datasets(33)));
    
    SMAP_year = double(hdf5read(hinfo.GroupHierarchy.Attributes(2)));
    SMAP_doy = double(hdf5read(hinfo.GroupHierarchy.Attributes(3)));
    secs_RowTime = double(hdf5read(hinfo.GroupHierarchy.Datasets(27))); %observation time for each row as UTC secs since 2015-1-1 00:00:00 UTC
    t_ref = datenum(2015,1,1); %reference time 
    row_time = datevec((secs_RowTime./86400) + t_ref);
    
    SMAP_dateVector = datenum(row_time);
%     SMAP_dateVector=datenum(ones(size(SMAP_RowTime))*SMAP_year,zeros(size(SMAP_RowTime)),ones(size(SMAP_RowTime))*SMAP_doy,zeros(size(SMAP_RowTime)),zeros(size(SMAP_RowTime)),SMAP_RowTime);
    time = repmat(SMAP_dateVector,1,76);
    
    % Designate ascending scans as 0 and descending scans as 1
    ascIndex = ones(size(lat));
    ascIndex(lat(1:end-1)<lat(2:end))=0;
    if(ascIndex(end-1)==0)
        ascIndex(end)=0;
    end
    
    %%Copying other variables for Eric Hackert HDF5 data
    anc_dir = double(hdf5read(hinfo.GroupHierarchy.Datasets(1)));
    anc_swh = double(hdf5read(hinfo.GroupHierarchy.Datasets(5)));
    azi_aft = double(hdf5read(hinfo.GroupHierarchy.Datasets(8)));
    azi_fore = double(hdf5read(hinfo.GroupHierarchy.Datasets(9)));
    inc_aft = double(hdf5read(hinfo.GroupHierarchy.Datasets(11)));
    inc_fore = double(hdf5read(hinfo.GroupHierarchy.Datasets(12)));
    land_fraction_aft = double(hdf5read(hinfo.GroupHierarchy.Datasets(13)));
    land_fraction_fore = double(hdf5read(hinfo.GroupHierarchy.Datasets(14)));
    n_h_aft = double(hdf5read(hinfo.GroupHierarchy.Datasets(17)));
    n_h_fore = double(hdf5read(hinfo.GroupHierarchy.Datasets(18)));
    n_v_aft = double(hdf5read(hinfo.GroupHierarchy.Datasets(19)));
    n_v_fore = double(hdf5read(hinfo.GroupHierarchy.Datasets(20)));
    nedt_h_aft = double(hdf5read(hinfo.GroupHierarchy.Datasets(21)));
    nedt_h_fore = double(hdf5read(hinfo.GroupHierarchy.Datasets(22)));
    nedt_v_aft = double(hdf5read(hinfo.GroupHierarchy.Datasets(23)));
    nedt_v_fore = double(hdf5read(hinfo.GroupHierarchy.Datasets(24)));
    num_ambiguities = double(hdf5read(hinfo.GroupHierarchy.Datasets(25)));
    smap_ambiguity_dir = double(hdf5read(hinfo.GroupHierarchy.Datasets(28)));
    smap_ambiguity_spd = double(hdf5read(hinfo.GroupHierarchy.Datasets(29)));
    smap_high_dir = double(hdf5read(hinfo.GroupHierarchy.Datasets(30)));
    smap_high_dir_smooth = double(hdf5read(hinfo.GroupHierarchy.Datasets(31)));
    smap_high_spd = double(hdf5read(hinfo.GroupHierarchy.Datasets(32)));
    smap_sss_uncertainty = double(hdf5read(hinfo.GroupHierarchy.Datasets(35)));
    tb_h_aft = double(hdf5read(hinfo.GroupHierarchy.Datasets(36)));
    tb_h_bias_adj = double(hdf5read(hinfo.GroupHierarchy.Datasets(37)));
    tb_h_fore = double(hdf5read(hinfo.GroupHierarchy.Datasets(38)));
    tb_v_aft = double(hdf5read(hinfo.GroupHierarchy.Datasets(39)));
    tb_v_bias_adj = double(hdf5read(hinfo.GroupHierarchy.Datasets(40)));
    tb_v_fore = double(hdf5read(hinfo.GroupHierarchy.Datasets(41)));
              
    %%
    %Flags
    Q_Flag = double(hdf5read(hinfo.GroupHierarchy.Datasets(26)));
    
    % Radiometer flags on is 16 bits, each bit represents a condition
    % Convert the flags from Unsigned integers to binary
    Flag = dec2bin(Q_Flag,16);
    
    %Usable data flag (bit#0) sss qf, 0 is good
    loc = 0 + 1; % Index of bit
    bit_num = 16 - loc + 1; % Matlab counts from left to right
    F0 = double(Flag(:,bit_num)) - 48; %1char = 49 double & 0char = 48 double
    F0 = reshape(F0,size(lat));
    
    %Usable data flag (bit#2) eia qf, 0 is Nominal incidence angles
    loc = 2 + 1; % Index of bit
    bit_num = 16 - loc + 1; % Matlab counts from left to right
    F2 = double(Flag(:,bit_num)) - 48; %1char = 49 double & 0char = 48 double
    F2 = reshape(F2,size(lat));
    
    %Usable data flag (bit#4) galaxy corr qf, 0 is All galaxy corrections < 5 K
    loc = 4 + 1; % Index of bit
    bit_num = 16 - loc + 1; % Matlab counts from left to right
    F4 = double(Flag(:,bit_num)) - 48; %1char = 49 double & 0char = 48 double
    F4 = reshape(F4,size(lat));
    
    %Usable data flag (bit#5) roughness qf, 0 is Ancillary wind speed < 20 m/s
    loc = 5 + 1; % Index of bit
    bit_num = 16 - loc + 1; % Matlab counts from left to right
    F5 = double(Flag(:,bit_num)) - 48; %1char = 49 double & 0char = 48 double
    F5 = reshape(F5,size(lat));
    
    %SST data flag (bit#6) low sst qf, 0 is SST > 5C
    loc = 6 + 1; % Index of bit
    bit_num = 16 - loc + 1; % Matlab counts from left to right
    F6 = double(Flag(:,bit_num)) - 48; %1char = 49 double & 0char = 48 double
    F6 = reshape(F6,size(lat));
    
    %Usable data flag (bit#7) land qf, 0 is No land detected in SWC
    loc = 7 + 1; % Index of bit
    bit_num = 16 - loc + 1; % Matlab counts from left to right
    F7 = double(Flag(:,bit_num)) - 48; %1char = 49 double & 0char = 48 double
    F7 = reshape(F7,size(lat));
    
    %Usable data flag (bit#8) ice qf, 0 is No ice detected in SWC
    loc = 8 + 1; % Index of bit
    bit_num = 16 - loc + 1; % Matlab counts from left to right
    F8 = double(Flag(:,bit_num)) - 48; %1char = 49 double & 0char = 48 double
    F8 = reshape(F8,size(lat));
        
    in_flag = find(F0~=0 | F2~=0 | F4~=0 | F5~=0 | F6~=0 | F7~=0 | F8~=0 | abs(lat)>90 | abs(lon)>180 | SMAP_ws==-9999 | NCEP_ws==-9999 | sst==-9999 | hyc==-9999 | sss==-9999);
    
    lon(in_flag) = nan; 
    lat(in_flag) = nan;
    SMAP_ws(in_flag) = nan;
    NCEP_ws(in_flag) = nan;
    sst(in_flag) = nan;
    hyc(in_flag) = nan;
    sss(in_flag) = nan;
    time(in_flag) = nan;
    ascIndex(in_flag) = nan;

    %% Finding Grid Indexing
    res = 0.1; %define resolution
    latt = lat;
    lont = lon;
%     lat_ind = round((-latt + 90.125)/0.25);
%     lat_ind(lat_ind==721) = 720;
%     lont(lont < 0) = lont(lont < 0) + 360;
%     lon_ind = round((lont + 0.125)/0.25);  
    
    lat_ind = round((90 - latt)/res) + 1;
%     lat_ind(lat_ind == 721) = 720; %for 25 km res
    lat_ind(lat_ind == 1801) = 1800; %for 10 km res
    lont(lont < 0) = lont(lont < 0) + 360;
    lon_ind = round(lont/res);
    
    col = round(360/res);
    
    if lon_ind == col
       lon_ind = 0;
    end
    lon_ind = lon_ind + 1; %so lon_ind starts at 1
    
    lat_ind(lat == -9999) = nan; %since lat_ind shouldn't be defined when
    %lat is not, using nan so it's not used in RIM (otherwise, if using 
    %-9999, it will have a value and I'm not sure what it does in RIM)
    lon_ind(lon == -9999) = nan; %same for lon_ind    
    
    data = struct('lon',lon,'lat',lat,'SMAP_ws',SMAP_ws,'NCEP_ws',NCEP_ws,'sst',sst,...
                'hyc',hyc,'sss',sss,...
                'lat_ind',lat_ind,'lon_ind',lon_ind,'time',time,'ascIndex',ascIndex);
    theRest = struct('anc_dir',anc_dir,'anc_swh',anc_swh,'azi_aft',azi_aft,'azi_fore',azi_fore,'inc_aft',inc_aft,...
        'inc_fore',inc_fore,...
        'num_ambiguities',num_ambiguities,'Q_Flag',Q_Flag,'row_time',row_time,'smap_ambiguity_dir',smap_ambiguity_dir,...
        'smap_ambiguity_spd',smap_ambiguity_spd,'smap_high_dir',smap_high_dir,'smap_high_dir_smooth',smap_high_dir_smooth,...
        'smap_high_spd',smap_high_spd,'smap_sss_uncertainty',smap_sss_uncertainty);
    AllAttributes = struct('attributes',attributes); %for Eric H.
            
    name_file=['SMAP_L2B_SSS_JPL_v5_',num2str(year), num2str(month,'%0.2d'),num2str(day,'%0.2d'),'_',num2str(revno,'%0.5d'),'_10kmgrd'];
    outputfile=[dirOutput,'/' name_file ,'.mat'];
    save(outputfile,'data','theRest','AllAttributes'); %for Eric H.
%     save(outputfile,'data');
    
	delete([name]);
	
    disp(['Processed file' num2str(k) ':' s(k).name])
end

