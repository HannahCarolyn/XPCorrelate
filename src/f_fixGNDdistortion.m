function [datastack,tform,tformxebsd] = f_fixGNDdistortion(ebsd,X,Y,fullres,Data_InputMap_before,Data_InputMap_after,Map_PH2_before,Map_PH2_after,GND_before,GND_after,microscope,primphase, EBSDREFq,resultsdir,fractionofspacing,stepsize);
% Script for correcting the distortion in the EBSD map based on what we
% presume to be a reliable hardness map. Selection of points are made 
% (at least 4) in both maps, and an affine transformation that links the
% hardness points onto the ebsd map is found. All coordinates on the 
% hardness map are then transformed to form "query coordinates"
% on the ebsd map. The result is then gridded, and outputted to datastack
% which contains the X, Y, Express, and EBSD data in the same size arrays.


%CMM and IH and AJW 2020.
% HC updated 2023

% if want to add extra data streams also add it too the function list to
% draw it in
%Datastack structure: 
% datastack.X               %X position
% datastack.Y               %Y position
% datastack.S               %Surface displacement
% datastack.D               %Depth
% datastack.L               %Load
% datastack.M               %Modulus
% datastack.St              %Stiffness^2/Load
% datastack.H               %Hardness
% datastack.phi1            %phi1  (ebsd)
% datastack.Phi             %Phi   (ebsd)
% datastack.phi2            %phi2  (ebsd)
% datastack.phase           %phase (ebsd)
% datastack.GNDtotalafter   %GNDtotal average after (xebsd) 
% datastack.GNDtotalbefore  %GNDtotal average before (xebsd) 

%% Correlate the before and after ebsd images extracting a transform
[tformxebsd]=f_xebsdcorrelate(Data_InputMap_before,Data_InputMap_after,resultsdir);


%% hardness input
%First make sure the data is appropriately structured. FROM HERE ON OUT,
%data should be structured in grids which start at 0,0, and X(2,1) is
%larger than X(1,1), and Y(1,2) is larger than Y(1,1). This goes against
%Matlab's conventional Row Column notation, but follows the cartesian form
%to relate position in matrices to positions on the cartesian grid
%described by X and Y. 

%quick check:
i=1;
while isnan(X(i,1)) || isnan(X(i+1,1))
    i=i+1;
end
if X(i+1,1)<X(i,1)
    X=flipud(X);
    Y=flipud(Y);
    fullres=flipud(fullres);
end
i=1;
while isnan(Y(1,i)) || isnan(Y(1,i+1))
    i=i+1;
end
if Y(1,i+1)<Y(1,i)
    X=fliplr(X);
    Y=fliplr(Y);
    fullres=fliplr(fullres);
end

if X(1,1)==0 && Y(1,1)==0 
    disp('H zerod correctly')
else
    X= X-min(min(X));
    Y= Y-min(min(Y));
end
Xspacing=X(2,1)-X(1,1); %Gives the spacing in x
Yspacing=Y(1,2)-Y(1,1); %Gives the spacing Y
%The following two lines extract the vectors that are used to describe all
%the possible points. See ndgrids matlab documentation for more help.
Xvector=X(:,1); 
Yvector=Y(1,:);


% Himage
H=fullres(:,:,6);%Extract hardness
H(H>1000)=0; %This cleans the data
figure 
hplot=contourf(X,Y,H,45,'LineColor','None'); %This plots the hardness data
xticks=(0:Xspacing:max(Xvector)); %This was aimed to help wit
yticks=(0:Yspacing:max(Yvector));
axis image
grid on
caxis([nanmean(H(:))-1*nanstd(H(:)) nanmean(H(:))+1*nanstd(H(:))]) %This gives a scale of the axis
title('Reference Selection in Hardness(>4)')
[x_og,y_og] = getpts; %obtain the "og" coordinates: the original, absolute reference frame for all future points.
close all

%plot it again with the points on for helpful guide, and save.
figure
hplot=contourf(X,Y,H,45,'LineColor','None'); %This plots the hardness data
hold on
scatter(x_og,y_og,'rx') %This plots the selected points on the hardness image.
for i=1:size(x_og,1)
    text(x_og(i)+0.01*max(X(:)), y_og(i)+0.01*max(Y(:)),num2str(i));   %put up the text at the given locations
end
title('References Selected On Nanoindentation Map')
xlabel('\mum')
ylabel('\mum')
axis image
c=colorbar;
c.Label.String = 'Hardness /GPa';
figname=['refselect_H figure'];
saveas(gcf,fullfile(resultsdir, figname),'png') %This saves the figure as a png

%% EBSD manipulation % HC need to work through the orientating for different microscopes and stuff at the momoment this section does nothing.
%We deal with EBSD data in the same way as hardness: based on 0,0 and going
%up in the correct way.

%check: are the x and y values changing in the same way as the
%hardness? i.e. does x(2,1) correspond to the pixel on the right of x(1,1)
% if x_ebsd(2,1)> x_ebsd(1,1) && y_ebsd(1,2)>y_ebsd(1,1) %if everythings fine
%     %do nothing
% elseif x_ebsd(2,1)<= x_ebsd(1,1) && y_ebsd(1,2)>y_ebsd(1,1)
%   
%     GNDtotalafter_ebsd=flipud(GNDtotalafter_ebsd);
%     
% 
% elseif x_ebsd(2,1)> x_ebsd(1,1) && y_ebsd(1,2)<=y_ebsd(1,1)
%     x_ebsd=fliplr(x_ebsd);
%     y_ebsd=fliplr(y_ebsd);
%     GNDtotalafter_ebsd=fliplr(GNDtotalafter_ebsd);
% 
% elseif x_ebsd(2,1)<= x_ebsd(1,1) && y_ebsd(1,2)<=y_ebsd(1,1)
%     x_ebsd=transpose(x_ebsd);
%     y_ebsd=transpose(y_ebsd);
%     GNDtotalafter_ebsd=fliplr(GNDtotalafter_ebsd);
% end
% 
% if x_ebsd(1,1)==0 && y_ebsd(1,1)==0 
%     %do nothing
% else
%     x_ebsd= x_ebsd-min(min(x_ebsd));
%     y_ebsd= y_ebsd-min(min(y_ebsd));
% end

%HOWEVER, each ebsd system will decide which was this is. I've tried to get
%the right info for each microscope, but IF THIS ISN'T RIGHT JUST CLOSE THE
%REF FIGURE AND IT WILL MANUALLY LET YOU PICK THE DIRECTION
% put the first point on (0, 0)
if strcmp(microscope,'merlin') 
    %now deprecated by using mtex conversion
elseif strcmp(microscope,'evo') 
    %if it's a evo:
    GNDtotal_ebsd_after=fliplr(GNDtotal_ebsd_after);
elseif strcmp(microscope,'tescan') 
    %if it's oxford instruments
    GNDtotal_ebsd_after=fliplr(GNDtotal_ebsd_after);
elseif strcmp(microscope,'xbeam') 
    %if it's the crossbeam instruments 
    %the xbeam is frustrating %I (HC) am masters student who never uses xbeam
    %so hasn't done the intergration on this very sorry :(
    [x_ebsd,y_ebsd,phi1, phi2, Phi, phase,BCebsd]=f_fixXbeam(x_ebsd,y_ebsd,phi1, phi2, Phi, phase,BCebsd);
end

%quick check:
% if x_ebsd(1,1)==0 && y_ebsd(1,1)==0 && x_ebsd(1,2)>x_ebsd(1,1) && y_ebsd(2,1)>y_ebsd(1,1)
%     disp('ebsd zerod correctly')
% else
%     error('ebsd data not correctly structured')
% end
%
%% K means clustering the PH2 before map to act as a phase map

% Map_PH2_before_KMC=imsegkmeans(uint16(Map_PH2_before),10,"NumAttempts",10);
% Map_PH2_before_KMC_double=(im2double(Map_PH2_before_KMC)*255);

%% Entropy filter of the map

% EntropyPH2=entropyfilt(Map_PH2_before);
% EntropyPH2im=rescale(EntropyPH2);
% figure;
% Eimage=imagesc(EntropyPH2);
% title("entropy filtered PH2")
% figure;
% PQ=imagesc(Data_InputMap_before.RadonQuality);
% title("PQ")
% Pqbefore=Data_InputMap_before.RadonQuality;
% figure;
% I=imfuse(EntropyPH2,Pqbefore);
% imshow(I);
% 
% 
% 
% 
% 
% Map_PH2_before_filtered = imbinarize(EntropyPH2im,0.8); % unhard code
% Map_PH2_before_filtered_dbl=double(Map_PH2_before_filtered);
% Map_PH2_before_filtered_dbl=Map_PH2_before_filtered_dbl+1;

%% select points on the manipulated ebsd. HC at the moment you can't reorientate
% Set up here the after data. Firstly you may question this oreintation the
% reason for this orientation is because the later gridded Interplot
% function doesn't work unless the ndgrid which is going in ascending
% order. Feel free to try yourself but I couldn't figure a better way of
% doing it.

%These are the x and y coordinates on the after image which is orietnated.
x_ebsd_after=flipud(((Data_InputMap_after.XSample)'));
y_ebsd_after=(Data_InputMap_after.YSample)';

%These are the reorietnated data for each of after data sets.
GNDtotal_ebsd_after = (fliplr(GND_after.total))';
PQ_ebsd_after=(fliplr(Data_InputMap_after.RadonQuality))';
PH2_ebsd_after=(fliplr(Map_PH2_after))';
%Feel free to add other data streams here for example strain. The format
%needed is:
%New_Data_ebsd_after=(fliplr(New_Data))';
%add in here the PH2 and other data

%These are the x and y coordinates on the before image which is orietnated.
x_ebsd_before=flipud(((Data_InputMap_before.XSample)'));
y_ebsd_before=(Data_InputMap_before.YSample)';


%These are the reorietnated data for each of before data sets.
GNDtotal_ebsd_before = (fliplr(GND_before.total))';
PQ_ebsd_before=(fliplr(Data_InputMap_before.RadonQuality))';
PH2_ebsd_before=(fliplr(Map_PH2_before))';
% Pearlite_ebsd_before=(fliplr(Pearlite_before))';
%Feel free to add other data streams here for example strain. The format
%needed is:
%New_Data_ebsd_before=(fliplr(New_Data))';
%add in here the PH2 and other data

%% plotting k means cluster map
% hplot=contourf(x_ebsd_before,y_ebsd_before,PH2_ebsd_before,45,'LineColor','None');
% title('Entropy Filtered PH2')
% xlabel('\mum')
% ylabel('\mum')
% axis xy
% colorbar
% figname=['Entropy filter map'];
% saveas(gcf,fullfile(resultsdir, figname),'png') %This saves the figure as a png in the results directory.
%% Plotting the after PQ in order to correlate the data to the nanoindentation

 figure
    hplot=contourf(x_ebsd_after,y_ebsd_after,PQ_ebsd_after,45,'LineColor','None'); %This plots the after ebsd data bare in mind that this takes time to run.
    colormap("gray")
    c=colorbar;
    c.Label.String = 'Band Contrast /arb units';
    hold on
    title('Reference Selection in EBSD (>4)')
    [x_transREF,y_transREF] = getpts; %This gets the points that are the same points on the nanoindentaiton data
    
    figname=['refselecting_EBSD figure'];
    

%plot it again with the points on for helpful guide, and save.
figure
     hplot=contourf(x_ebsd_after,y_ebsd_after,PQ_ebsd_after,45,'LineColor','None'); %This plots the after ebsd data bare in mind that this takes time to run.
    colormap("gray")
    c=colorbar;
    c.Label.String = 'Band Contrast /arb units';
hold on
scatter(x_transREF,y_transREF) %This plots the selected points on the after ebsd for reference.
for i=1:size(x_transREF,1)
    text(x_transREF(i)+0.01*max(X(:)), y_transREF(i)+0.01*max(Y(:)),num2str(i));   %put up the text at the given locations
end
title('References Selected On EBSD Map')
xlabel('\mum')
ylabel('\mum')
axis xy
figname=['refselect_EBSD figure'];
saveas(gcf,fullfile(resultsdir, figname),'png') %This saves the figure as a png in the results directory.
close all %This closes all the figures.

%% transformation calculation

movingPoints  = [x_og    y_og   ]; %The moving points are the reference choosen on the nanoindentation points
fixedPoints = [x_transREF y_transREF]; %The fixed points are the reference choosen on the after ebsd points

transtype='affine'; %select what type of transformation you want 
tform = fitgeotrans(movingPoints,fixedPoints,transtype); % This creates a tform which is matrix that correlates the after ebsd points to the nanoindentation points.
% transform the data
A = tform.T; %THis was in Chris' code I don't think it does anything.

%transform the x and y coordinates of the nanoindentation data
%using the affine transformation
[locebsdcorr(:,1),locebsdcorr(:,2)] = transformPointsForward(tform,X(:),Y(:)); %This transforms the nanoindentation points in the nanoindentation reference frame to the nanoindentation points in the after ebsd reference frame.
locebsdcorrafter=locebsdcorr; %This just helps me remember that these are the nanoindentation points in the after ebsd reference frame.

%% Plot transformed nanoindentation points in the frame of reference of the after ebsd on the after PQ map.
%plot it again with the points on after for helpful guide, and save.

figure
 hplot=contourf(x_ebsd_after,y_ebsd_after,PQ_ebsd_after,45,'LineColor','None'); %This plots the after PQ map. Please bare in mind that this function takes time to run.
colormap("gray")
c=colorbar;
c.Label.String = 'Band Contrast /arb units';
axis xy %this is the default 
hold on
scatter(locebsdcorrafter(:,1),locebsdcorrafter(:,2),'rx') %plot the transformed nanoindentation points in the reference frame of after ebsd on the after ebsd image.
title('Transformed nanoindentation points on after PQ Map')
xlabel('\mum')
ylabel('\mum')
figname=['Transformed_nanoindentation_after  figure'];
saveas(gcf,fullfile(resultsdir, figname),'png')

%% Transform the transformed nanoindentation points in the reference frame of the after ebsd to the reference frame of the before ebsd
%This uses the tform which is created using the f_xebsdcorrelate which
%gives a matrix which transforms between the before and after image. Okay
%so this function takes the nanoindenation points in the after reference
%frame and uses the transfrom points inverse to get these point in the
%reference frame of the before ebsd.

[locebsdcorrbefore(:,1),locebsdcorrbefore(:,2)] = transformPointsInverse(tformxebsd,locebsdcorrafter(:,1),locebsdcorrafter(:,2)); 

%This plots the transformed nanoindentation points in the reference frame
%of the before ebsd and plots it onto the before PQ image.
figure
hplot=contourf(x_ebsd_before,y_ebsd_before,PQ_ebsd_before,45,'LineColor','None'); %This plots the before PQ image. Bare in mind that this function takes time to run.
colormap("gray")
axis xy %this is the default 
hold on
scatter(locebsdcorrbefore(:,1),locebsdcorrbefore(:,2),'rx') %This plots the transformed nanoindenation points
title('Transformed nanoindentation points on before PQ Map')
xlabel('\mum')
ylabel('\mum')
c=colorbar;
c.Label.String = 'Band Contrast /arb units';
figname=['Transformed_nanoindentation_before  figure'];
saveas(gcf,fullfile(resultsdir, figname),'png') %This saves the png of the figure in the results directory.

%% This interplots the ebsd so there is a value for every coordinate within the map
%This function requires ndgrids in ascending order
GNDtotalbefore_ebsdinterp = griddedInterpolant(x_ebsd_before,y_ebsd_before,GNDtotal_ebsd_before,'nearest'); %this is the before ebsd data
GNDtotalafter_ebsdinterp = griddedInterpolant(x_ebsd_after,y_ebsd_after,GNDtotal_ebsd_after,'nearest'); %this is the after ebsd data
PH2_before_ebsdinterp = griddedInterpolant(x_ebsd_before,y_ebsd_before,PH2_ebsd_before,'nearest');
PH2_after_ebsdinterp = griddedInterpolant(x_ebsd_after,y_ebsd_after,PH2_ebsd_after,'nearest');
% Pearlite_before_ebsdinterp = griddedInterpolant(x_ebsd_before,y_ebsd_before,Pearlite_ebsd_before,'nearest');

% Feel free to add other data streams. The format you need is:
%Data_New_before_ebsdinterp = griddedInterpolant(x_ebsd_before,y_ebsd_before,Data_New_ebsd_before,'nearest');
%Data_New_after_ebsdinterp = griddedInterpolant(x_ebsd_after,y_ebsd_after,Data_New_ebsd_after,'nearest');
%% This sections finds a single sample point on the ebsd for each nanoindentation point.
%I made this section since for the first popin the only the very tip of the
%nanoindentation is in and this is around the step size of the ebsd. This
%is what is done in the f_fixEBSDdistortion section.
GNDbeforesinglepoint = GNDtotalbefore_ebsdinterp(locebsdcorrbefore(:,1), locebsdcorrbefore(:,2));
GNDbeforesinglepointG = f_gridify_vector(GNDbeforesinglepoint,size(X,1),size(Y,2))';

GNDaftersinglepoint = GNDtotalafter_ebsdinterp(locebsdcorrafter(:,1), locebsdcorrafter(:,2));
GNDaftersinglepointG = f_gridify_vector(GNDbeforesinglepoint,size(X,1),size(Y,2))';

%Feel free to add other data streams. The format you need is:
% Data_New_beforesinglepoint = Date_New_before_ebsdinterp(locebsdcorrbefore(:,1), locebsdcorrbefore(:,2));
% Data_New_beforesinglepointG = f_gridify_vector(Data_New_beforesinglepoint,size(X,1),size(Y,2))';
% 
% Data_New_aftersinglepoint = Data_New_after_ebsdinterp(locebsdcorrafter(:,1), locebsdcorrafter(:,2));
% Data_New_aftersinglepointG = f_gridify_vector(Data_New_aftersinglepoint,size(X,1),size(Y,2))';

%% This sections sets up the more sampling points in the reference frame of the nanoindentation.
%The reason for this section is that it allows for the fact more than just
%one EBSD point impacts the nanoindentation area in fact it will be the
%an area that impacts it.

Xvector=X(:,1); %extract ndgrid vectors
Yvector=Y(1,:);

Xspacing=X(2,1)-X(1,1); %get the spacing
Yspacing=Y(1,2)-Y(1,1);
% fractionofspacing=1.0;
% stepsize=0.25;% At the moment this is hard coded since all my ebsd data had a step size of 0.25um will change at some point
Noofsamplesxhalf=round(((Xspacing*fractionofspacing)/2)/stepsize); %Uses half the spacing between nanoindentation indents and divides it by the step size to get sampling points
Noofsamplesyhalf=round(((Yspacing*fractionofspacing)/2)/stepsize); %This does the same for y.


for noofxpoints=1:numel(Xvector) %For each of the x data points in nanoindentation ndgrid vector
    Xvectorlower=[];
for sampling=0:Noofsamplesxhalf %For each sampling point in x below x values
    indexproblem=sampling+1; %fixs zero index point
    Xvectorlower(indexproblem)=Xvector(noofxpoints)-(stepsize*sampling); %this makes a vector with lower half values of x sampling wanted
end
    Xvectorlower=fliplr(Xvectorlower); %This slips the vector so going from lowest to highest
    Xvectorlowersize=numel(Xvectorlower); %Gives the size of this vector
    
for sampling=1:Noofsamplesxhalf %for each sampling point in x above x value.
    
    Xvectorlower(Xvectorlowersize+sampling)=Xvector(noofxpoints)+(stepsize*sampling); %This adds the sampling points above the value to the vector
end
XNewcoordinate(noofxpoints,:)=Xvectorlower; %This puts the expanded x coordinates into a matrix.
end

%This is just a repeat of the above section but for the y coordinates.
for noofypoints=1:numel(Yvector)
    Yvectorlower=[];
for sampling=0:Noofsamplesyhalf
    indexproblem=sampling+1;
    Yvectorlower(indexproblem)=Yvector(noofypoints)-(stepsize*sampling);
end
    Yvectorlower=fliplr(Yvectorlower);
    Yvectorlowersize=numel(Yvectorlower);
    
for sampling=1:Noofsamplesyhalf

    Yvectorlower(Yvectorlowersize+sampling)=Yvector(noofypoints)+(stepsize*sampling);
end
YNewcoordinate(noofypoints,:)=Yvectorlower;
end

for noofxpoints=1:numel(Xvector) %for each of the x coordinates in the og nanoindentation coordinate system
    Xboxvector=XNewcoordinate(noofxpoints,:); %grab the new expanded coordinates
    for noofypoints=1:numel(Yvector) %for each of the y coordinates in the og nanoindentation coordinate system
    Yboxvector=YNewcoordinate(noofypoints,:); %grab the new expanded coordinates


[xboxmatrix,yboxmatrix]=ndgrid(Xboxvector,Yboxvector); %for each x and y coordinate combination create a ndgrid of the points from the vectors.

%take the box ndgrid matrix and convert to vector and transform the points
%from the reference frame of the nanoindenter to the reference frame of the
%after ebsd.
[locebsdcorrboxafter(:,1),locebsdcorrboxafter(:,2)] = transformPointsForward(tform,xboxmatrix(:),yboxmatrix(:));

 %This next line converts the expanded sampling nanoindentation data that
 %has been transformed to the reference frame of the after ebsd and
 %transfers it to the before ebsd reference frame.
 [locebsdcorrboxbefore(:,1),locebsdcorrboxbefore(:,2)] = transformPointsInverse(tformxebsd,locebsdcorrboxafter(:,1),locebsdcorrboxafter(:,2));
 
 % This section is for the GND data
 GNDtotalafter_ebsdnew_box = GNDtotalafter_ebsdinterp(locebsdcorrboxafter(:,1), locebsdcorrboxafter(:,2)); %This finds the point that is nearest the transformed box nanoindentation points in the reference frame of the after ebsd which has been interplated.
 GNDtotalafter_ebsdnew_box(GNDtotalafter_ebsdnew_box == 0)=NaN; %This takes any values that are 0 and converts them to NAN.
 GNDtotalafter_ebsdnew_boxaverage(noofxpoints,noofypoints)=nanmean(GNDtotalafter_ebsdnew_box); %This takes all the ebsd values at the transformed nanoindentation points that are not nans and averages them. Then puts the values into a matrix that is in the coordinate system of the og nanoindentation.
 

 %This section does the same as the after ebsd refernce frame
 %nanoindenation expanded sample points but in the before reference frame.
 GNDtotalbefore_ebsdnew_box = GNDtotalbefore_ebsdinterp(locebsdcorrboxbefore(:,1), locebsdcorrboxbefore(:,2));
 GNDtotalbefore_ebsdnew_box(GNDtotalbefore_ebsdnew_box == 0)=NaN;
 GNDtotalbefore_ebsdnew_boxaverage(noofxpoints,noofypoints)=nanmean(GNDtotalbefore_ebsdnew_box);

 %This section is for the peak height 2 data

 PH2_after_ebsdnew_box = PH2_after_ebsdinterp(locebsdcorrboxafter(:,1), locebsdcorrboxafter(:,2)); 
 PH2_after_ebsdnew_box(PH2_after_ebsdnew_box == 0)=NaN; 
 PH2_after_ebsdnew_boxaverage(noofxpoints,noofypoints)=nanmean(PH2_after_ebsdnew_box);

 PH2_before_ebsdnew_box = PH2_before_ebsdinterp(locebsdcorrboxbefore(:,1), locebsdcorrboxbefore(:,2)); 
 PH2_before_ebsdnew_box(PH2_before_ebsdnew_box == 0)=NaN; 
 PH2_before_ebsdnew_boxaverage(noofxpoints,noofypoints)=nanmean(PH2_before_ebsdnew_box);
 
%This section is for Pearlite map
% Pearlite_before_ebsdnew_box = Pearlite_before_ebsdinterp(locebsdcorrboxbefore(:,1), locebsdcorrboxbefore(:,2)); 
%  Pearlite_before_ebsdnew_box(Pearlite_before_ebsdnew_box == 0)=NaN; 
%  Pearlite_before_ebsdnew_boxaverage(noofxpoints,noofypoints)=nanmean(Pearlite_before_ebsdnew_box);
%  
 %If you would like to add another data stream the format that is needed
 %is:
%  Data_New_after_ebsdnew_box = Data_New_after_ebsdinterp(locebsdcorrboxafter(:,1), locebsdcorrboxafter(:,2)); 
%  Data_New_after_ebsdnew_box(Data_New_after_ebsdnew_box == 0)=NaN; 
%  Data_New_after_ebsdnew_boxaverage(noofxpoints,noofypoints)=nanmean(Data_New_after_ebsdnew_box);

%  Data_New_before_ebsdnew_box = Data_New_after_ebsdinterp(locebsdcorrboxbefore(:,1), locebsdcorrboxbefore(:,2)); 
%  Data_New_before_ebsdnew_box(Data_New_before_ebsdnew_box == 0)=NaN; 
%  Data_New_before_ebsdnew_boxaverage(noofxpoints,noofypoints)=nanmean(Data_New_before_ebsdnew_box);
    end
end
%% Rotating final matrix's from microscope rotation to be the same rotation as the nanoindentation data matrixs.



%output into a class datastack with reasonable names
datastack.X     = X;                 %X position
datastack.Y     = Y;                 %Y position
datastack.S     = fullres(:,:,1);   %Surface displacement
datastack.D     = fullres(:,:,2);   %Depth
datastack.L     = fullres(:,:,3);   %Load
datastack.M     = fullres(:,:,4);   %Modulus
datastack.St    = fullres(:,:,5);   %Stiffness^2/Load
datastack.H     = fullres(:,:,6);   %Hardness
datastack.GNDtotalafterpoint = GNDaftersinglepointG;    %GND after single point (xebsd)
datastack.GNDtotalbeforepoint = GNDbeforesinglepointG;  %GND before single point (xebsd)
datastack.GNDaverageafter = GNDtotalafter_ebsdnew_boxaverage;       %GNDtotal average after (xebsd) 
datastack.GNDaveragebefore = GNDtotalbefore_ebsdnew_boxaverage;       %GNDtotal average before (xebsd) 
datastack.PH2averageafter = PH2_after_ebsdnew_boxaverage;       %PH2 average after
datastack.PH2averagebefore = PH2_before_ebsdnew_boxaverage;     %PH2 average before
% datastack.Pearliteaveragebefore=Pearlite_before_ebsdnew_boxaverage; %Pearlite average before.
%If you would like to add extra data streams feel free the format that is
%needed is: %don't need to have single and box average if only added one
%type of stream.
% datastack.Data_Newafterpoint = Data_New_aftersinglepointG;    
% datastack.Data_Newbeforepoint = Data_New_beforesinglepointG;  
% datastack.Data_Newaverageafter = Data_New_after_ebsdnew_boxaverage;       
% datastack.Data_Newaveragebefore = Data_New_before_ebsdnew_boxaverage;

end