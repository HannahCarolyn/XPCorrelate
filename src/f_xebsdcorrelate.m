%%you need the Data_InputMap before and after results directory

%%you need out the transform

function [tform]=f_xebsdcorrelate(Data_InputMap_before,Data_InputMap_after,resultsdir);
%% Section for ensuring the data is orientated such that it will go through the code
% The y and x grids need to be ndgrids which go in ascending order inorder
% to acheive this the following rotations are needed. 

x_ebsd_after=flipud(((Data_InputMap_after.XSample)')); %Reoreintating the x coordinates after data 
y_ebsd_after=(Data_InputMap_after.YSample)'; % Reorientating the y coordinates after data
PQ_ebsd_after=(fliplr(Data_InputMap_after.RadonQuality))'; % Reorientating the Pattern Quality data after


x_ebsd_before=flipud(((Data_InputMap_before.XSample)'));%Reoreintating the x coordinates before data 
y_ebsd_before=(Data_InputMap_before.YSample)';%Reoreintating the y coordinates before data 
PQ_ebsd_before=(fliplr(Data_InputMap_before.RadonQuality))'; %Reoreintating the Pattern Quality before data 

%% Getting Points for Correlating the Before and After indents Pattern Quality Maps

%Plotting the After Pattern Quality Map for choosing points
figafterselectionplot=figure;
figure(figafterselectionplot)


hplot=contourf(x_ebsd_after,y_ebsd_after,PQ_ebsd_after,100,'LineColor','None'); %This function plots the PQ data note this function sometimes takes time to run.
colormap("gray") 
title('Reference Picks After Indents Pattern Quality >4') %For an affine transformation at least 4 points to carry out
axis image
[x_after,y_after] = getpts; %This is function that allows for the clicking in order to get points that can be recorded and used later on in the code.
axis off; hold off

% This plots the After Pattern Quality Map with the choosen points labelled
% such that it helps with choosing the same points on the after map
figafterselected=figure;
figure(figafterselected)

hplot=contourf(x_ebsd_after,y_ebsd_after,PQ_ebsd_after,100,'LineColor','None'); %This function plots the PQ data note this function takes some time
colormap("gray")
title('References Selected On After Pattern Map') 
xlabel('\mum')
ylabel('\mum')
axis image
hold on
scatter(x_after,y_after,'rx') %This takes gets each of the points and plots them as red crosses.
for i=1:size(x_after,1)
    text(x_after(i), y_after(i),num2str(i));   %put up the text at the given locations
end
c=colorbar;
c.Label.String = 'Pattern Quality';
figname=['refselect_after figure'];
saveas(gcf,fullfile(resultsdir, figname),'png') %This saves the figure as a png to your results directory


%Plotting the Before Pattern Quality Map for choosing points
figbeforeelectionplot=figure;
figure(figbeforeelectionplot)

hplot=contourf(x_ebsd_before,y_ebsd_before,PQ_ebsd_before,100,'LineColor','None');  %This function plots the PQ data note this function takes some time
colormap("gray")
title('Before Indents Pattern Quality')
axis image
[x_before,y_before] = getpts; %This gets the reference points for the before points.
axis off; hold off


% This plots the Before Pattern Quality Map with the choosen points
% labelled.

figafterselected=figure;
figure(figafterselected)
hplot=contourf(x_ebsd_before,y_ebsd_before,PQ_ebsd_before,100,'LineColor','None');  %This function plots the PQ data note this function takes some time
colormap("gray")
hold on
scatter(x_before,y_before,'rx') %This plots the selected points on the before PQ indents
for i=1:size(x_before,1)
    text(x_before(i), y_before(i),num2str(i));   %put up the text at the given locations
end
title('References Selected On Before Pattern Map')
xlabel('\mum')
ylabel('\mum')
axis image
c=colorbar;
c.Label.String = 'Pattern Quality';
figname=['refselect_before figure'];
saveas(gcf,fullfile(resultsdir, figname),'png') %This saves the figure as png in the results directory.
%% transformation calculation

movingPoints  = [x_before y_before]; %These are the moving points which have been choosen as the selected points from the before map.
fixedPoints = [x_after    y_after]; %These are you fixed points which have been choosen as the selected points from the after map.



transtype='affine'; %select what type of transformation you want %you likely want affine but could change to simiarity (look at the MatLab page for more information)
tform = fitgeotrans(movingPoints,fixedPoints,transtype); %This is the function that does creates a tform which is matrix which translate the points from one reference frame to the other.
% transform the data
A = tform.T; %This doesn't really do anything but was in Chris' code.


%transform the x and y coordinates of the nanoindentation data
%using the affine transformation
[Xbeforecorr,Ybeforecorr] = transformPointsForward(tform,x_before,y_before); %This transforms the before coordinates in the reference frame of the before image to the before coordinates in the reference frame of the after image.

%% A plot to check that the transformation has worked

% This plots the after pattern quality with the transformed before
% coordinates that are now in the after reference frame in order to ensure
% that transfomation has worked.
figure;
hplot=contourf(x_ebsd_after,y_ebsd_after,PQ_ebsd_after,100,'LineColor','None'); %Plots the after data. Note that not colour scaled to grey to help make it clear it's a different graph
axis off; hold on
scatter(Xbeforecorr,Ybeforecorr,'rx') %Plot the transformed after points
for i=1:size(Xbeforecorr,1)
    text(Xbeforecorr(i), Ybeforecorr(i),num2str(i));   %put up the text at the given locations
end
title('Before Coordinates transformed to After Pattern Map')
xlabel('\mum')
ylabel('\mum')
axis image
c=colorbar;
c.Label.String = 'Pattern Quality';
figname=['backcorrelatedpositions figure'];
saveas(gcf,fullfile(resultsdir, figname),'png') %save the figures 
close all %close all the figures before moving on
end