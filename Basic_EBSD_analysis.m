%% EBSD Analysis

%% Clear and close all

clear all; close all; clc

%% Define parameters

% Define x and y directions to match the SEM images
plotx2east; plotzIntoPlane

% Phase to be analyzed
Phase_chosen = 'Magnesium';

% Degree threshold of grain boundary
GrainAngle = 10;

% Grain too small to be excluded
grain_exclude = 10;

% When (length of twin boundary)/(length of grain boundary) > threshold, the grain is treated as a twin
threshold = 0.6;

%% Load EBSD data

[fname, pname] = uigetfile('*', 'Choose .crc EBSD data');
cd(pname)

ebsd = EBSD.load([pname '\' fname],'convertEuler2SpatialReferenceFrame','setting 2');
ebsd_raw = ebsd;
CS = ebsd.CS;

% plot raw map
figure
[~,mP] = plot(ebsd,ebsd.orientations);
print(gcf,[fname(1:end-4) '_RawMap'],'-dpng','-r400');
% saveas(gcf,[fname(1:end-4) '_RawMap'],'fig');

mP.micronBar.visible = 'off';
print(gcf,[fname(1:end-4) '_RawMap_NoScaleBar'],'-dpng','-r400');

% plot ipf key
ipfKey = ipfColorKey(ebsd(Phase_chosen));
figure;
plot(ipfKey)
print(gcf,[fname(1:end-4) '_ipfKey'],'-dpng','-r400');

%% Denose and reconstruct grains
% Here fill, no grains, indexed pixels (all will be filled) are choosen for reconstruction

% make sure the data is the original data imported
ebsd=ebsd_raw;

% smooth and fill missing
prompt = ['Level to denoise EBSD data. Input 1:High;  Input 2:Middle;  Input 3:Low'];
titlename = ['EBSD analysis'];
option_denoise = inputdlg(prompt,titlename,[1 100]);
option_denoise = str2num(option_denoise{1});
F = splineFilter;
% F = infimalConvolutionFilter;
if option_denoise == 1
    [grains,ebsd('indexed').grainId] = calcGrains(ebsd('indexed'),'angle',GrainAngle*degree);
    % remove very small grains
    ebsd(grains(grains.grainSize<=10)) = [];
    % redo grains reconstruction
    [grains,ebsd('indexed').grainId] = calcGrains(ebsd('indexed'),'angle',GrainAngle*degree);
    % smooth and fill missing
    F = splineFilter;
    ebsd = smooth(ebsd('indexed'),F,'fill');
    [grains,ebsd('indexed').grainId] = calcGrains(ebsd('indexed'),'angle',GrainAngle*degree);
elseif option_denoise == 2
    % reconstruct the grains structure
    [grains,ebsd('indexed').grainId,ebsd.mis2mean] = calcGrains(ebsd('indexed'));
    % remove very small grains
    ebsd(grains(grains.grainSize<3)) = [];
    % smooth and fill missing
    F = splineFilter;
    ebsd = smooth(ebsd('indexed'),F,'fill',grains);
    % redo the reconstruction of grains(why do not change the mis2mean?)
    [grains,ebsd('indexed').grainId] = calcGrains(ebsd('indexed'));
else
    [grains,ebsd('indexed').grainId] = calcGrains(ebsd('indexed'));
    ebsd(grains(grains.grainSize<=3)) = [];
    [grains,ebsd('indexed').grainId] = calcGrains(ebsd('indexed'));
    F = splineFilter;
    ebsd = smooth(ebsd('indexed'),F);
    [grains,ebsd('indexed').grainId] = calcGrains(ebsd('indexed'));
end

grains = grains.smooth(5);

%% Plot orientation map

%%%% plot X direction orientation map %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% defines an ipf color key for the Magnesium phase
ipfKey = ipfColorKey(ebsd(Phase_chosen));

% set the referece direction to X
ipfKey.inversePoleFigureDirection = vector3d.X;

% compute the colors
colors = ipfKey.orientation2color(ebsd(Phase_chosen).orientations);

% plot orientation map
figure;
[~,mP] = plot(ebsd(Phase_chosen),colors)
hold on
plot(grains.boundary,'linewidth',1.5)
hold off
print(gcf,[fname(1:end-4) '_IPFX'],'-dpng','-r400');

mP.micronBar.visible = 'off';
print(gcf,[fname(1:end-4) '_IPFX_NoScaleBar'],'-dpng','-r400');

%%%% plot Y direction orientation map %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% set the referece direction to Y
ipfKey.inversePoleFigureDirection = vector3d.Y;

% compute the colors
colors = ipfKey.orientation2color(ebsd(Phase_chosen).orientations);

% plot orientation map
figure;
[~,mP] = plot(ebsd(Phase_chosen),colors)
hold on
plot(grains.boundary,'linewidth',1.5)
hold off
print(gcf,[fname(1:end-4) '_IPFY'],'-dpng','-r400');

mP.micronBar.visible = 'off';
print(gcf,[fname(1:end-4) '_IPFY_NoScaleBar'],'-dpng','-r400');

%%%% plot Z direction orientation map %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% set the referece direction to Z
ipfKey.inversePoleFigureDirection = vector3d.Z;

% compute the colors
colors = ipfKey.orientation2color(ebsd(Phase_chosen).orientations);

% plot orientation map
figure;
[~,mP] = plot(ebsd(Phase_chosen),colors)
hold on
plot(grains.boundary,'linewidth',1.5)
hold off
print(gcf,[fname(1:end-4) '_IPFZ'],'-dpng','-r400');

mP.micronBar.visible = 'off';
print(gcf,[fname(1:end-4) '_IPFZ_NoScaleBar'],'-dpng','-r400');

%% Plot pole/inverse pole figure

%%%% make sure ED direction in the center of pole figure %%%%%%%%%%%%%%%%
% x,y,x are direction of sample; a,b,c are direction of crystal
% rotate the direction to make ED at the center
%rot = rotation.byAxisAngle(yvector,90*degree);
%ebsd = rotate(ebsd,rot,'keepEular');

%%%% plot pole figure %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% set annotation in pole figure
pfAnnotations = @(varargin) text([vector3d.X,vector3d.Y],{'ED','TD'},...
  'BackgroundColor','w','tag','axesLabels',varargin{:});
setMTEXpref('pfAnnotations',pfAnnotations);

rot = rotation.byAxisAngle(yvector,0*degree);
ebsd_rotated = rotate(ebsd_raw,rot);

% choose the plane for mapping
h = [Miller(0,0,0,1,ebsd(Phase_chosen).CS),...
    Miller(1,0,-1,0,ebsd(Phase_chosen).CS),...
    Miller(1,0,-1,1,ebsd(Phase_chosen).CS),...
    Miller(1,1,-2,2,ebsd(Phase_chosen).CS)];%,Miller(1,0,-1,2,ebsd(Phase_chosen).CS)];

% plot pole figure
figure;
plotPDF(ebsd_rotated(Phase_chosen).orientations,h,'contourf','minmax','antipodal');
% add colorbar
CLim(gcm,'equal');
mtexColorbar
% mtexColorMap parula
% mtexColorMap summer
% mtexColorMap hsv
% mtexColorMap cool
print(gcf,[fname(1:end-4) '_Pole'],'-dpng','-r400');

figure;
plotPDF(ebsd_rotated(Phase_chosen).orientations,h(1),'contourf','minmax','antipodal');
mtexColorbar
print(gcf,[fname(1:end-4) '_Pole1'],'-dpng','-r400');

figure;
plotPDF(ebsd_rotated(Phase_chosen).orientations,h(2),'contourf','minmax','antipodal');
mtexColorbar
print(gcf,[fname(1:end-4) '_Pole2'],'-dpng','-r400');

figure;
plotPDF(ebsd_rotated(Phase_chosen).orientations,h(3),'contourf','minmax','antipodal');
mtexColorbar
print(gcf,[fname(1:end-4) '_Pole3'],'-dpng','-r400');

figure;
plotPDF(ebsd_rotated(Phase_chosen).orientations,h(4),'contourf','minmax','antipodal');
mtexColorbar
print(gcf,[fname(1:end-4) '_Pole4'],'-dpng','-r400');

%%%% plot inverse pole figure %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure;
plotIPDF(ebsd(Phase_chosen).orientations,[xvector,yvector,zvector],'contourf');
% add colorbar
% CLim(gcm,'equal');
mtexColorbar
print(gcf,[fname(1:end-4) '_InvPole'],'-dpng','-r400');

%% Plot KAM map
% % specifying a threshold 2.5 degree to remove subgrain boundaries
%	after specifying the threshold, the remaining KAM is very sensitive to measurement errors and often noisy
% % specifying higher order of neighbors to reduce the meansurement errors
%   higher order smoothes away local dislocation structures
% % it is better to denoise first and compute the KAM from first order neighbors
figure;
[~,mP] = plot(ebsd,ebsd.KAM('threshold',2.5*degree,'order',1) ./ degree,'off')
caxis([0,2])
mtexColorbar
mtexColorMap LaboTeX
hold on
plot(grains.boundary,'lineWidth',1.5)
hold off
print(gcf,[fname(1:end-4) '_KAM'],'-dpng','-r400');

mP.micronBar.visible = 'off';
print(gcf,[fname(1:end-4) '_KAM_NoScaleBar'],'-dpng','-r400');


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%  Grain boundary analysis
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% extract new grain boundaries
gB = grains.boundary;
gB_MgMg = gB(Phase_chosen,Phase_chosen);

%% plot grain boundaries with misorientation angles

figure
[~,mP] = plot(gB_MgMg,gB_MgMg.misorientation.angle./degree,'linewidth',2);
mtexColorbar
saveas(gcf,[fname(1:end-4) '_grain boundary'],'fig');
print(gcf,[fname(1:end-4) '_grain boundary'],'-dpng','-r400');

mP.micronBar.visible = 'off';
print(gcf,[fname(1:end-4) '_grain boundary_NoScaleBar'],'-dpng','-r400');

%% distribution of grain angles

figure
histogram(gB_MgMg.misorientation.angle./degree,40,'Normalization','probability');
xlabel('Misorientation angle (degree)')
ylabel('Frequency (%)')
set(gca,'FontSize',15,'FontName','Times New Roman')
boxaxes
saveas(gcf,[fname(1:end-4) '_angles distribution'],'fig');
print(gcf,[fname(1:end-4) '_angles distribution'],'-dpng','-r400');

%% Grain size distribution

r_all = equivalentRadius(grains);
ind_grain = find(r_all > grain_exclude);
r_chosen = r_all(ind_grain);
grain_size0 = mean(r_chosen);

area_grains = sum(area(grains))
relative_area = area(grains)./area_grains;
grain_size = sum(r_all.*relative_area);

figure
a = histogram(r_chosen,15,'Normalization','probability');
xlabel('Equivelent diameter (\mum)')
ylabel('Frequency (%)')
boxaxes
set(gca,'FontSize',15,'FontName','Times New Roman')
print(gcf,[fname(1:end-4) '_GrainSize'],'-dpng','-r400');

%% Extract the tensile twinning

% define the tensile twin angle of Mg
tensile_twinning = gB.misorientation.angle>81.3*degree & gB.misorientation.angle<91.3*degree;
twinBoundary_tensile = gB(tensile_twinning);

% extract twinning grains
twin_id0 = unique(twinBoundary_tensile.grainId);
twin_tensile_id = [];
twin_tensile_row = [];
twin_tensile_row_all = [];
for i = 1:length(twin_id0)
    grain_id = twin_id0(i);
    [twin_row,~] = find(twinBoundary_tensile.grainId == grain_id);
    twin_row = unique(twin_row);
    twin_length_percent = sum(twinBoundary_tensile(twin_row).segLength)./...
        sum(grains(grain_id).boundary.segLength);
    if twin_length_percent > threshold
        twin_tensile_id = [twin_tensile_id grain_id];
        twin_tensile_row = [twin_tensile_row;twin_row];
% % Unfinished part to improve the detect of twins
    elseif twin_length_percent > 0.35
        possible_id = unique(twinBoundary_tensile(twin_row).grainId);
        if ~max(ismember(twin_tensile_id,possible_id)) & length(possible_id) < 4 & min(area(grains(possible_id)))/max(area(grains(possible_id))) < 0.25
            [~,b] = min(area(grains(possible_id)));
            twin_tensile_id = [twin_tensile_id possible_id(b)];
            twin_tensile_row = [twin_tensile_row;twin_row];
        end
    end
    twin_tensile_row_all = [twin_tensile_row_all;twin_row];
end

% The percentage of twins
area_grains = sum(area(grains));
area_tensile_twins = sum(area(grains(twin_tensile_id)));
twin_tensile_percent = area_tensile_twins / area_grains *100;
disp(['The percentage of tensile twins is ' num2str(twin_tensile_percent) '%'])

%% Extract compress twinning

% define the compress twin angle of Mg
compress_twinning = gB.misorientation.angle>51*degree & gB.misorientation.angle<61*degree;
twinBoundary_compress = gB(compress_twinning);

% extract twinning grains
twin_id0 = unique(twinBoundary_compress.grainId);
twin_compress_id = [];
twin_compress_row = [];
twin_compress_row_all = [];
for i = 1:length(twin_id0)
    grain_id = twin_id0(i);
    [twin_row,~] = find(twinBoundary_compress.grainId == grain_id);
    twin_row = unique(twin_row);
    twin_length_percent = sum(twinBoundary_compress(twin_row).segLength)./...
        sum(grains(grain_id).boundary.segLength);
    if twin_length_percent > threshold
        twin_compress_id = [twin_compress_id grain_id];
        twin_compress_row = [twin_compress_row;twin_row];
    end
    twin_compress_row_all = [twin_compress_row_all;twin_row];
end

% The percentage of twins
area_grains = sum(area(grains));
area_compress_twins = sum(area(grains(twin_compress_id)));
twin_compress_percent = area_compress_twins / area_grains *100;
disp(['The percentage of compress twins is ' num2str(twin_compress_percent) '%'])

%% Plot twins

figure;
plot(ebsd,ebsd.orientations);
hold on
plot(grains.boundary,'linewidth',0.5);
plot(twinBoundary_tensile(twin_tensile_row_all),'linecolor','w','linewidth',3);
if ~isempty(twin_compress_row_all)
    plot(twinBoundary_compress(twin_compress_row_all),'linecolor','y','linewidth',3);
end
hold off
% h1 = get(gca);
% set(h1.Legend,'color','none','edgecolor','none')
saveas(gcf,[fname(1:end-4) '_TwinBoundaries'],'fig');
print(gcf,[fname(1:end-4) '_TwinBoundaries'],'-dpng','-r400');

figure;
plot(ebsd,ebsd.orientations);
hold on
plot(grains.boundary,'linewidth',0.5);
plot(grains(twin_tensile_id).boundary,'linecolor','w','linewidth',3);
if ~isempty(twin_compress_row)
    plot(grains(twin_compress_id).boundary,'linecolor','y','linewidth',3);
end
hold off
saveas(gcf,[fname(1:end-4) '_TwinGrains'],'fig');
print(gcf,[fname(1:end-4) '_TwinGrains'],'-dpng','-r400');

%% Save results

% filename = [fname(1:end-4) '_Twin_result.mat'];
% save(filename)

% Write results into text file
fid=fopen([fname(1:end-4) '_Grain_statistics.txt'],'w');
fprintf(fid,' ---  Results of tension twins --- \n');
fprintf(fid,'%s %f\n','Tension area:',area_tensile_twins);
fprintf(fid,'%s %f \n','Area percentage:',twin_tensile_percent);
fprintf(fid,'\n');
fprintf(fid,' ---  Results of compression twins --- \n');
fprintf(fid,'%s %f\n','Compression area:',area_compress_twins);
fprintf(fid,'%s %f \n','Area percentage:',twin_compress_percent);
fprintf(fid,'\n');
fprintf(fid,'%s %f\n','Area of whole sample:',area_grains);
fprintf(fid,'%s %d \n','Number of grains:', length(ind_grain));
fprintf(fid,'%s %f \n','Threshold degree of grain boundary:', threshold);
fprintf(fid,'\n');
fprintf(fid,'%s %f\n','Grain size based on radius:',grain_size0);
fprintf(fid,'%s %f \n','Grain size weighted by area:', grain_size);
fprintf(fid,'%s %f \n','Exclude grain size <', grain_exclude);

fclose(fid);

%% Local function

function [a,b]=boxaxes
% get handle to current axes
% a = gca;
% % % get the original xlim, ylim
% % xlim_ori = get(a,'XLim');
% % ylim_ori = get(a,'YLim');
% % set box property to off and remove background color
% set(a,'box','off','color','none')
% % create new, empty axes with box but without ticks
% b = axes('Position',get(a,'Position'),'box','on','xtick',[],'ytick',[],'xlim',get(a,'XLim'),'ylim',get(a,'YLim'));
% % set original axes as active
% axes(a)
% % link axes in case of zooming
% linkaxes([a b])

box off;
a = gca;
x = get(a,'XLim');
y = get(a,'YLim');
line([x(1) x(2)],[y(2) y(2)],'Color','k','linewidth',0.54)
line([x(2) x(2)],[y(1) y(2)],'Color','k','linewidth',0.54)
xlim(x)
ylim(y)
end
