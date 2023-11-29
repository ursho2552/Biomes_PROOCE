function [h,Fig,cmap2] = plotSOM( new_data,month,nr_classes)
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

lons = -179.5:179.5;
lats = -89.5:89.5;
latlim = [-78 80];
lonlim = [20.5 379.5];

h = figure('color','white');
hold on

axesm('robinson','MapLatLimit',latlim,'MapLonLimit',lonlim, ...
   'Frame','on','Grid','on','MeridianLabel','on','ParallelLabel','on')%...
   
axis off
setm(gca,'MLabelLocation',60,'LabelFormat','none','MLabelParallel','South');
setm(gca,'PLabelLocation',45,'LabelFormat','none');

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
    
else
    for i=1:nlat
        for j=1:nlon
            Fig(i,j) = new_data(1,i,j);
        end
    end
end

[lat,lon] = meshgrat(lats,round(lons));

t = geoshow(double(lat),double(lon),Fig,'DisplayType','texturemap') 
set(t,'FaceAlpha','texturemap','AlphaData',double(~isnan(Fig)));

if(isnan(nr_classes))
    colormap(parula);
    
elseif(strcmp(nr_classes,'black'))
    colormap(gray(1));
    
else
    cmap = parula(nr_classes);
    [cmap_tmp] = shuffle_colormap(cmap);
    [cmap_tmp] = shuffle_colormap(cmap_tmp);
    %colormap(cmap)
    cmap2 = colormap(cmap_tmp);

caxis([1 nr_classes]);
end
