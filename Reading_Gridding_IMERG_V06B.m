    %% This code reads the V06B IMERG files and saves in Matlab format
% generated by March Jacob, Rev1: 04/03/2019
%updated on 08/18/2020 to plot before and after gridding, since the
%processed files doesn't seem to be the same as the originals
%changed <0 = nan to after gridding, in case nan are propagating
%the code was change to do the regular gridding and averaging

function Reading_Gridding_IMERG_V06B_Rev2(yyyy)

% yyyy = year of the data to process
% dd1 = initial day of the processing data interval
% dd2 = final day of the processing data interval
clc
% dir1=pwd;
%dir1='/disk4/OceanSalinity/RIM/IMERG/Data/GriddedData/PerOrbit';
%dir2='/disk4/OceanSalinity/RIM/IMERG/Data/RawData/V06B';
dir1 = 'Z:\OceanSalinity\RIM\IMERG\Data\GriddedData\PerOrbit\V06B';
dir2 = 'Z:\OceanSalinity\RIM\IMERG\Data\RawData\V06B';
dir3 = 'Z:\OceanSalinity\RIM\IMERG\Analysis\beforeafter';

outdir = [dir1 '/']; % output folder
dir_y = [dir2 '/'];%input folder
figdir = [dir3 '/'];%to save the figures

res = 0.25; %grid size resolution

%% Define grids
col = round(360/res);
row = round(180/res);

imerg = zeros(row,col);
lati = zeros(row,col);
longi = zeros(row,col);
count = zeros(row,col);

% for d = mm1:mm2 

    files = dir([dir_y,'3B-HHR.MS.MRG.3IMERG.',num2str(yyyy),'*']);
    
    for i = 1:size(files,1) % specify what orbit of the day
		if files(i).bytes<100000
            % bad/small file - skip
            continue
        end
        
        name = [dir_y,files(i).name];
        disp(name)
        
         %%Reading variables
         lon = double(hdf5read(name,'/Grid/lon')); % Longitude 0.1 deg res
         lat = double(hdf5read(name,'/Grid/lat')); % Latitude 0.1 deg res
        
        %precipitation
        prec = double(hdf5read(name,'/Grid/precipitationCal')); %imerg precip in mm/hr for 0.1 deg res
        precip = reshape(prec,[1,6480000]); 
        
        mm = str2num(files(i).name(26:27)); %month from the file name
        dd = str2num(files(i).name(28:29)); %day from the file name
        rev = str2num(files(i).name(47:50)); %rev from the file name (each rev every 30 min)
        
		%precip(precip<0)=nan; % kd add
        
        %% Finding Grid Indexing
        latt = repmat(lat,1,3600);
        lont = lon';
        lont = repmat(lont,1800,1);
        lat_ind = round((-latt + 90.125)/res);
        lat_ind(lat_ind==721) = 720;
        lont(lont < 0) = lont(lont < 0) + 360;
        lon_ind = round((lont + 0.125)/res);
        
        imerg = zeros(row,col);
        lati = zeros(row,col);
        longi = zeros(row,col);
        count = zeros(row,col);
        
        %%Gridding variables to 0.25 deg res
        for j = 1:length(precip)
            r = lat_ind(j);
            c = lon_ind(j);
            
            imerg(r,c) = imerg(r,c) + precip(j);
            lati(r,c) = lati(r,c) + latt(j);
            longi(r,c) = longi(r,c) + lont(j);
            
            count(r,c) = count(r,c) + 1;            
        end
        
        imerg(count > 0) = imerg(count > 0)./count(count > 0);
        lati(count > 0) = lati(count > 0)./count(count > 0);
        longi(count > 0) = longi(count > 0)./count(count > 0);
       
        %imerg = imresize(precip,[720 1440]);
        
        imerg(imerg<0)=nan; % kd add > change it in case it's propagating
        
        figure,clf
        set(gcf, 'Position', get(0, 'Screensize'));
        load coastlines
        coastlon_shift=coastlon;
        coastlon_shift(coastlon_shift<0)=coastlon_shift(coastlon_shift<0)+360;
        dc=diff(coastlon_shift);
        coastlon_shift(abs(dc)>350)=nan;
         prec(prec<0)=nan;
        subplot(211)
        cax = [0 2];
        turbo_cscat(lont,latt,prec,'.');grid on
        shading flat
        caxis(cax)
        title('imerg original')
        hold on
        plot(coastlon_shift,coastlat,'w','LineWidth',2)
        ylim([-90 90]);
        xlim([0 360]);
        
        subplot(212)
        turbo_cscat(longi,lati,imerg,'.');grid on
        shading flat
        caxis(cax)
        title('imerg resampled to quarter deg grid')
        hold on
        plot(coastlon_shift,coastlat,'w','LineWidth',2)
        ylim([-90 90]);
        xlim([0 360]);
               
        %% Saving Output
        nameFile = ['IMERG_V06B_30min_quarter_' num2str(yyyy) num2str(mm,'%0.2d') num2str(dd,'%0.2d') '_' num2str(rev,'%0.4d')];
        save([outdir nameFile '.mat'], 'imerg')
        
%         %nameFig1 = [figdir 'IMERG_V06B_' num2str(yyyy) num2str(mm,'%0.2d') num2str(dd,'%0.2d') '_' num2str(rev,'%0.4d') '.fig'];
%         nameFig2 = ['IMERG_V06B_' num2str(yyyy) num2str(mm,'%0.2d') num2str(dd,'%0.2d') '_' num2str(rev,'%0.4d')];
%         %savefig(nameFig1);
%         saveas(gcf,strcat(figdir, nameFig2), 'png');
%         %saveas(gcf,nameFig2);
        
        close all;
        
        disp(['done ',num2str(yyyy),num2str(mm,'%0.2d'),num2str(dd,'%0.2d'),'-',num2str(rev,'%0.4d')])
    end
% end
disp('done')
end