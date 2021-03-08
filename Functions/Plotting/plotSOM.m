function [h,Fig,cmap2] = plotSOM( new_data,month,nr_classes)%,p)
%Function to plot a map

%{
Parameters:
    new_data (matrix): Map data with dimension n x 180 x 360
    month (int): section of the map new_data that should be plotted.
        Usually is the month. Thus month is between 1 and n
    nr_classes (int): Number of unique classes/labels new_data. If 
        nr_classes is given then discrete colors are shown. If NaN is given
        the colorbar will be continuous
   
 Output:
    h (Figure object): Figure handle
    Fig (matrix): 2D version of the matrix that is plotted
    cmap2 (matrix): Colormap used to plot the map

%}

new_classes


lons = -179.5:179.5;
lats = -89.5:89.5;
%latlim = [-89.5 89.5]; % set map limits for plots
latlim = [-78 80];
%lonlim = [-179.5 179.5];
lonlim = [20.5 379.5];
h = figure('color','white');
%subplot(1,2,1,'Parent',p)
hold on

% axesm('eqdcylin','MapLatLimit',latlim,'MapLonLimit',lonlim, ...
%    'Frame','on',...
%    'FEdgeColor',[0.5 0.5 0.5],'FLineWidth',0.5,...
%    'MeridianLabel','on','ParallelLabel','on', 'FEdgeColor',[0.5 0.5 0.5],'FLineWidth',0.5,'FontSize',12);
%'eckert4'

axesm('robinson','MapLatLimit',latlim,'MapLonLimit',lonlim, ...
   'Frame','on','Grid','on','MeridianLabel','on','ParallelLabel','on')%...
%    'FEdgeColor',[0.5 0.5 0.5],'FLineWidth',0.5,...
%    'MeridianLabel','on','ParallelLabel','on', 'FEdgeColor',[0.5 0.5 0.5],'FLineWidth',0.5,'FontSize',12);
axis off
setm(gca,'MLabelLocation',60,'LabelFormat','none','MLabelParallel','South');
setm(gca,'PLabelLocation',45,'LabelFormat','none');
%geoshow('landareas.shp','EdgeColor',[0.5 0.5 0.5], 'FaceColor', [0.5 0.5 0.5]); % grey land mass
load coastlines
plotm(coastlat,coastlon,'k')

levels=[0.5 1 1.5 2.0];

nlat = length(lats);
nlon = length(lons);

Fig = NaN(nlat,nlon); %add one for the definition of the point 0.5 or 180.5
if(month < 13)
    for i=1:nlat
        for j=1:nlon
            Fig(i,j) = new_data(month,i,j);
        end
    end
    %Fig(:,361) = new_data(month,:,1);
else
    for i=1:nlat
        for j=1:nlon
            Fig(i,j) = new_data(1,i,j);
        end
    end
    %Fig(:,361) = new_data(1,:,1);
end
sum(isnan(Fig(:,1)));
sum(isnan(Fig(:,end)));

[lat,lon] = meshgrat(lats,round(lons));

t = geoshow(double(lat),double(lon),Fig,'DisplayType','texturemap') 
set(t,'FaceAlpha','texturemap','AlphaData',double(~isnan(Fig)));
% geoshow(double(lat),double(lon),Fig,'DisplayType','contour','LevelList',0.5,'LineColor','red')
% geoshow(double(lat),double(lon),Fig,'DisplayType','contour','LevelList',1,'LineColor','magenta')
% land = shaperead('landareas', 'UseGeoCoords', true)

% geoshow(land, 'FaceColor', [1 1 1])

aa = unique(Fig(~isnan(Fig)));
size(Fig);
% [lat,lon] = meshgrat(latchl,[lonchl;180.5]);
% surfm(lat,lon,Fig);
%surfm(double(latchl),double([lonchl]), Fig);
%surf(double(latchl),double([lonchl]),Fig);
%jj = input('dfd')
% contour(lon,lat, Fig,[0.5,1,1.5,2])
%colorbar
%
if(isnan(nr_classes))
    colormap(parula);
    %colormap(parula(nr_classes))
    %caxis([1 nr_classes])
elseif(strcmp(nr_classes,'black'))
    colormap(gray(1));
else
    cmap = parula(nr_classes);
    [cmap_tmp] = shuffle_colormap(cmap);
    [cmap_tmp] = shuffle_colormap(cmap_tmp);
    %colormap(cmap)
    cmap2 = colormap(cmap_tmp);
%colorbar%('Ticks', [aa])
caxis([1 nr_classes]);
end
%caxis([aa(1) aa(end)])
%cb = colorbar;
%aa = unique(Fig(~isnan(Fig)));
%set(cb, 'ytick', 1:aa(end), 'yticklabel', cellstr(num2str(aa))); 
%colorbar('Ticks',[-5,-2,1,4,7]
%}


