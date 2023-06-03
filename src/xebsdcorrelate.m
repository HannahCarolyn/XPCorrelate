close all;
clear;
clc
%load Mtex
Mtex_Loc='Z:\2022 Pt IIs\Hannah Cole\mtex-5.8.0\startup_mtex';
run(Mtex_Loc);

%Before indents xebsd output file path and .mat filename
folder_before='Z:\2022 Pt IIs\Hannah Cole\Merlin\sample62701\'; % folder wth input file - output folder to here too
file_before='HRGrid1sample6_XEBSD.mat';%input file with xEBSD(v3) data
%After indents xebsd output file path and .mat filename
folder_after='Z:\2022 Pt IIs\Hannah Cole\Merlin\Sample61702\';
file_after='HRGrid1sample6withindents_XEBSD.mat';


input_file_before=[folder_before file_before]; %file path for inputting before data
%load xebsd output from before indents
load(input_file_before, 'Map_stress_sample', 'Map_strain_sample' , 'Map_rotation_sample', 'Map_RefID',...
       'Map_EBSD_MTEX', 'Grains_MTEX', 'Map_MAE2', 'Map_PH2', 'Data_InputMap', 'MicroscopeData', 'iSet', 'GrainID_Setup', 'GND'); 
%rename each variable loaded to make it distinguish the data from being
%before the indents xebsd.
Map_stress_sample_before = Map_stress_sample;
clear('Map_stress_sample');
Map_strain_sample_before = Map_strain_sample;
clear('Map_strain_sample');
Map_rotation_sample_before = Map_rotation_sample;
clear('Map_rotation_sample');
Map_RefID_before = Map_RefID;
clear('Map_RefID');
Map_EBSD_MTEX_before = Map_EBSD_MTEX;
clear('Map_EBSD_MTEX');
Grains_MTEX_before = Grains_MTEX;
clear('Grains_MTEX');
Map_MAE2_before = Map_MAE2;
clear('Map_MAE2');
Map_PH2_before = Map_PH2;
clear('Map_PH2');
Data_InputMap_before = Data_InputMap;
clear('Data_InputMap');
MicroscopeData_before = MicroscopeData;
clear('MicroscopeData');
iSet_before = iSet;
clear('iSet');
GrainID_Setup_before = GrainID_Setup;
clear('GrainID_Setup');
GND_before = GND;
clear('GND');





input_file_after=[folder_after file_after]; %file path for inputting before data
%load xebsd output from after indents
load(input_file_after, 'Map_stress_sample', 'Map_strain_sample' , 'Map_rotation_sample', 'Map_RefID',...
       'Map_EBSD_MTEX', 'Grains_MTEX', 'Map_MAE2', 'Map_PH2', 'Data_InputMap', 'MicroscopeData', 'iSet', 'GrainID_Setup', 'GND'); 
%rename each variable loaded to make it distinguish the data from being
%after the indents xebsd.
Map_stress_sample_after = Map_stress_sample;
clear('Map_stress_sample');
Map_strain_sample_after = Map_strain_sample;
clear('Map_strain_sample');
Map_rotation_sample_after = Map_rotation_sample;
clear('Map_rotation_sample');
Map_RefID_after = Map_RefID;
clear('Map_RefID');
Map_EBSD_MTEX_after = Map_EBSD_MTEX;
clear('Map_EBSD_MTEX');
Grains_MTEX_after = Grains_MTEX;
clear('Grains_MTEX');
Map_MAE2_after = Map_MAE2;
clear('Map_MAE2');
Map_PH2_after = Map_PH2;
clear('Map_PH2');
Data_InputMap_after = Data_InputMap;
clear('Data_InputMap');
MicroscopeData_after = MicroscopeData;
clear('MicroscopeData');
iSet_after = iSet;
clear('iSet');
GrainID_Setup_after = GrainID_Setup;
clear('GrainID_Setup');
GND_after = GND;
clear('GND');
%% Plotting pattern quality before and after indents
figafterselectionplot=figure;
figure(figafterselectionplot)


image(Data_InputMap_after.RadonQuality,'CDataMapping','scaled');
colormap("gray")
title('Reference Picks After Indents Pattern Quality >3')
axis image
axis ij 
[x_after,y_after] = getpts;
axis off; hold off



figafterselected=figure;
figure(figafterselected)

image(Data_InputMap_after.RadonQuality,'CDataMapping','scaled');
colormap("gray")
title('Reference Picks After Indents Pattern Quality >3')
axis image
axis ij 
axis off; hold on
scatter(x_after,y_after,'rx')
for i=1:size(x_after,1)
    text(x_after(i), y_after(i),num2str(i));   %put up the text at the given locations
end
title('References Selected On After Pattern Map')
xlabel('\mum')
ylabel('\mum')
axis image
c=colorbar;
c.Label.String = 'Pattern Quality';
% figname=['refselect_after_input figure'];
% saveas(gcf,fullfile(resultsdir, figname),'png')



figafterelectionplot=figure;
figure(figafterelectionplot)

image(Data_InputMap_before.RadonQuality,'CDataMapping','scaled');
colormap("gray")
title('Before Indents Pattern Quality')
axis image
axis ij 
[x_before,y_before] = getpts;
axis off; hold off

figafterselected=figure;
figure(figafterselected)
image(Data_InputMap_before.RadonQuality,'CDataMapping','scaled');
colormap("gray")
title('Reference Picks Before Indents Pattern Quality >4')
axis image
axis ij 
axis off; hold on
scatter(x_before,y_before,'rx')
for i=1:size(x_before,1)
    text(x_before(i), y_before(i),num2str(i));   %put up the text at the given locations
end
title('References Selected On Before Pattern Map')
xlabel('\mum')
ylabel('\mum')
axis image
c=colorbar;
c.Label.String = 'Pattern Quality';
% figname=['refselect_after_input figure'];
% saveas(gcf,fullfile(resultsdir, figname),'png')
%% transformation calculation

movingPoints  = [x_before y_before];
fixedPoints = [x_after    y_after];


transtype='similarity'; %select what type of transformation you want 
tform = fitgeotrans(movingPoints,fixedPoints,transtype); 
% transform the data
A = tform.T;
PQ_before=Data_InputMap_before.RadonQuality;
PQ_after=Data_InputMap_after.RadonQuality;

xbeforevector=Data_InputMap_before.XSample(:);
ybeforevector=Data_InputMap_before.YSample(:);

[xbeforecorr,ybeforecorr] = transformPointsForward(tform,xbeforevector,ybeforevector);
 
xbeforecorrmatrix=reshape(xbeforecorr,[342,456]);%unhard cide this
ybeforecorrmatrix=reshape(xbeforecorr,[342,456]);


%random samples of points
randomsamplexindex=randi([1 155952],1,5);
xbeforesample=xbeforevector(randomsampleindex);
ybeforesample=ybeforevector(randomsampleindex);


%transform the x and y coordinates of the nanoindentation data
%using the affine transformation
[Xsampleaftercorr,Ysampleaftercorr] = transformPointsForward(tform,xbeforesample,ybeforesample);


figure;
image(PQ_after,'CDataMapping','scaled');
axis image
axis ij 
axis off; hold on
scatter(Xsampleaftercorr,Ysampleaftercorr,'rx')
for i=1:size(Xsampleaftercorr,1)
    text(Xsampleaftercorr(i), Ysampleaftercorr(i),num2str(i));   %put up the text at the given locations
end
title('Before point sample Transformed on After Pattern Quality Map (Checking Transformation)')
xlabel('\mum')
ylabel('\mum')
axis image
c=colorbar;
c.Label.String = 'Pattern Quality';
% figname=['refselect_after_input figure'];
% saveas(gcf,fullfile(resultsdir, figname),'png')

figure;
image(PQ_before,'CDataMapping','scaled');
axis image
axis ij 
axis off; hold on
scatter(xbeforesample,ybeforesample,'rx')
for i=1:size(xbeforesample,1)
    text(xbeforesample(i), ybeforesample(i),num2str(i));   %put up the text at the given locations
end
title('Before point sample on Before Pattern Quality Map (Checking Transformation)')
xlabel('\mum')
ylabel('\mum')
axis image
c=colorbar;
c.Label.String = 'Pattern Quality';
%% Choosing references for the nanoindentation corrlation
figafterselectionplot=figure;
figure(figafterselectionplot)


image(Data_InputMap_after.RadonQuality,'CDataMapping','scaled');
colormap("gray")
title('Reference Picks for XP correlate Pattern Quality >4')
axis image
axis ij 
[x_after_indent,y_after_indent] = getpts;
axis off; hold off



figafterselected=figure;
figure(figafterselected)

image(Data_InputMap_after.RadonQuality,'CDataMapping','scaled');
colormap("gray")
title('Reference Picked for XP correlate Pattern Quality >4')
axis image
axis ij 
axis off; hold on
scatter(x_after_indent,y_after_indent,'rx')
for i=1:size(x_after_indent,1)
    text(x_after_indent(i), y_after_indent(i),num2str(i));   %put up the text at the given locations
end
title('References Selected On After Pattern Map')
xlabel('\mum')
ylabel('\mum')
axis image
c=colorbar;
c.Label.String = 'Pattern Quality';
% figname=['refselect_after_input figure'];
% saveas(gcf,fullfile(resultsdir, figname),'png')

%% Revert references back onto original map


[xbeforecorr_indent,ybeforecorr_indent] = transformPointsInverse(tform,x_after_indent,y_after_indent);
figure;
image(Data_InputMap_before.RadonQuality,'CDataMapping','scaled');
colormap("gray")
axis image
axis ij 
axis off; hold on
scatter(xbeforecorr_indent,ybeforecorr_indent,'rx')
for i=1:size(xbeforecorr_indent,1)
    text(xbeforecorr_indent(i), ybeforecorr_indent(i),num2str(i));   %put up the text at the given locations
end
title('References Selected Indent On Before Pattern Map')
xlabel('\mum')
ylabel('\mum')
axis image
c=colorbar;
c.Label.String = 'Pattern Quality';
% figname=['refselect_after_input figure'];
% saveas(gcf,fullfile(resultsdir, figname),'png')