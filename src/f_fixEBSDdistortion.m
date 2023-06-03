function [datastack] = f_fixEBSDdistortion(ebsd,X,Y,fullres,Data_InputMap_before,Data_InputMap_after,GND_before,GND_after,microscope,primphase, EBSDREFq,resultsdir);
% Script for correcting the distortion in the EBSD map based on what we
% presume to be a reliable hardness map. Selection of points are made 
% (at least 4) in both maps, and an affine transformation that links the
% hardness points onto the ebsd map is found. All coordinates on the 
% hardness map are then transformed to form "query coordinates"
% on the ebsd map. The result is then gridded, and outputted to datastack
% which contains the X, Y, Express, and EBSD data in the same size arrays.

%CMM and IH and AJW 2020.

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
Xspacing=X(2,1)-X(1,1);
Yspacing=Y(1,2)-Y(1,1);
Xvector=X(:,1);
Yvector=Y(1,:);

% Himage
H=fullres(:,:,6);%Extract hardness
H(H>1000)=0;
figure
hplot=contourf(X,Y,H,45,'LineColor','None');
xticks=(0:Xspacing:max(Xvector));
yticks=(0:Yspacing:max(Yvector));
axis image
grid on
caxis([nanmean(H(:))-1*nanstd(H(:)) nanmean(H(:))+1*nanstd(H(:))])
title('Reference Selection in Hardness(>4)')
[x_og,y_og] = getpts; %obtain the "og" coordinates: the original, absolute reference frame for all future points.
close all

%plot it again with the points on for helpful guide, and save.
figure
hplot=contourf(X,Y,H,45,'LineColor','None');
hold on
scatter(x_og,y_og,'rx')
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
saveas(gcf,fullfile(resultsdir, figname),'png')
%% ebsd input 
% EBSD settings and cleaning:
%

if strcmp(primphase,'default')==1
    primphase=char(ebsd.mineralList(2));
end
setMTEXpref('xAxisDirection','west');

if EBSDREFq==0 || EBSDREFq==2
    ebsd=ebsd(primphase); 
end
%clean things up
[Grains,ebsd.grainId] = calcGrains(ebsd,'boundary','convexhull','angle',10*degree); %calc grains - this is useful for cleaning up
large_grains = Grains(Grains.grainSize >= 3); % HC may not want this?
ebsd_clean = ebsd(large_grains);
%fill EBSD
% delete nonindexed and fill them with nn existing phase 
ebsd_clean('n')=[] ;
ebsd_clean = fill(ebsd_clean) ;
ebsd=ebsd_clean; %HC have a look into the cleaning
%% EBSD manipulation
%We deal with EBSD data in the same way as hardness: based on 0,0 and going
%up in the correct way.
%[ebsdGrid] = gridify(ebsd);
% now set the x and y coordinates as what they are in the ebsd map
% (cropped and cleaned)
x_ebsd = Data_InputMap_after.XSample;
y_ebsd = Data_InputMap_after.YSample;

        
%extract the things we care about: angles % HC extract other things from
%here.
% if EBSDREFq==0 || EBSDREFq==2 || EBSDREFq==3
%     phi1 = ebsdGrid.orientations.phi1;
%     Phi  = ebsdGrid.orientations.Phi;
%     phi2 = ebsdGrid.orientations.phi2;
% elseif EBSDREFq==1
%     allphi1s=zeros(size(x_ebsd,1),size(y_ebsd,2),length(ebsdGrid.mineralList)-1);
%     allPhis=zeros(size(x_ebsd,1),size(y_ebsd,2),length(ebsdGrid.mineralList)-1);
%     allphi2s=zeros(size(x_ebsd,1),size(y_ebsd,2),length(ebsdGrid.mineralList)-1);
% 
%     for i = 1:length(ebsdGrid.mineralList)-1 %go through all the phases and get their orientations
%         try
%             phasenebsd=gridify(ebsdGrid(ebsdGrid.mineralList(i+1)));
%         catch
%             i=i+1;
%             if i<=length(ebsdGrid.mineralList)-1
%                 phasenebsd=gridify(ebsdGrid(ebsdGrid.mineralList(i+1)));
%             end
%         end
%         allphi1s(:,:,i) = phasenebsd.orientations.phi1;
%         allPhis(:,:,i)  = phasenebsd.orientations.Phi;
%         allphi2s(:,:,i) = phasenebsd.orientations.phi2;
%     end
%     phi1=mean(allphi1s,3,'omitnan'); delete allphi1s %get the phi of the relevant phase
%     Phi=mean(allPhis,3,'omitnan'); delete allPhis %don't worry about mean
%     phi2=mean(allphi2s,3,'omitnan'); delete allphi2s %if it's not that phase it's nan and this ignores nan
% end
% phase= ebsdGrid.phase;
 % loading extra details from HREBSD code variables
GNDtotalbefore_ebsd = GND_before.total;
GNDtotalafter_ebsd = GND_after.total;
%check: are the x and y values changing in the same way as the
%hardness? i.e. does x(2,1) correspond to the pixel on the right of x(1,1)
if x_ebsd(2,1)> x_ebsd(1,1) && y_ebsd(1,2)>y_ebsd(1,1) %if everythings fine
    %do nothing
elseif x_ebsd(2,1)<= x_ebsd(1,1) && y_ebsd(1,2)>y_ebsd(1,1)
    x_ebsd=flipud(x_ebsd);
    y_ebsd=flipud(y_ebsd);
    GNDtotalafter_ebsd=flipud(GNDtotalafter_ebsd);

elseif x_ebsd(2,1)> x_ebsd(1,1) && y_ebsd(1,2)<=y_ebsd(1,1)
    x_ebsd=fliplr(x_ebsd);
    y_ebsd=fliplr(y_ebsd);
    GNDtotalafter_ebsd=fliplr(GNDtotalafter_ebsd);

elseif x_ebsd(2,1)<= x_ebsd(1,1) && y_ebsd(1,2)<=y_ebsd(1,1)
    x_ebsd=transpose(x_ebsd);
    y_ebsd=transpose(y_ebsd);
    GNDtotalafter_ebsd=transpose(GNDtotalafter_ebsd);

end

if x_ebsd(1,1)==0 && y_ebsd(1,1)==0 
    %do nothing
else
    x_ebsd= x_ebsd-min(min(x_ebsd));
    y_ebsd= y_ebsd-min(min(y_ebsd));
end

%HOWEVER, each ebsd system will decide which was this is. I've tried to get
%the right info for each microscope, but IF THIS ISN'T RIGHT JUST CLOSE THE
%REF FIGURE AND IT WILL MANUALLY LET YOU PICK THE DIRECTION
% put the first point on (0, 0)
if strcmp(microscope,'merlin') 
    %now deprecated by using mtex conversion
    %need to rezero since not using mtex anymore 
elseif strcmp(microscope,'evo') 
    %if it's a evo:
    GNDtotalafter_ebsd=fliplr(GNDtotalafter_ebsd);
elseif strcmp(microscope,'tescan') 
    %if it's oxford instruments
    GNDtotalafter_ebsd=fliplr(GNDtotalafter_ebsd);
elseif strcmp(microscope,'xbeam') 
    %if it's the crossbeam instruments 
    %the xbeam is frustrating %I (HC) am masters student who never uses xbeam
    %so hasn't done the intergration on this very sorry :(
    [x_ebsd,y_ebsd,phi1, phi2, Phi, phase,BCebsd]=f_fixXbeam(x_ebsd,y_ebsd,phi1, phi2, Phi, phase,BCebsd);
end

%quick check:
if x_ebsd(1,1)==0 && y_ebsd(1,1)==0 && x_ebsd(2,1)>x_ebsd(1,1) && y_ebsd(1,2)>y_ebsd(1,1)
    disp('ebsd zerod correctly')
else
    error('ebsd data not correctly structured')
end
%
%% select points on the manipulated ebsd
% IF THIS IS THE WRONG WAY, just close the figure
try
    disp('If the axes are the wrong way, close figure to select the right ones')
    ebsdfig=figure;
    hplot=image(Data_InputMap_after.RadonQuality,'CDataMapping','scaled');
    colormap("gray")
    axis ij
    hold on
    axis image
    hold off
    title('Reference Selection in EBSD (>4)')
    [x_transREF,y_transREF] = getpts;
    close(ebsdfig)
catch
    %[x_ebsd,y_ebsd,phi1,Phi,phi2,phase,BCebsd,x_transREF,y_transREF]=f_ebsdwrongway(x_ebsd,y_ebsd,phi1,Phi,phi2,phase,BCebsd,EBSDREFq);
end

%plot it again with the points on for helpful guide, and save.
figure
% if EBSDREFq==0 
%     hplot=contourf(x_ebsd,y_ebsd,Phi,45,'LineColor','None');
%     c=colorbar;
%     c.Label.String = 'Declination angle /^{o}';
hplot=contourf(x_ebsd,y_ebsd,Data_InputMap_after.RadonQuality,45,'LineColor','None');
colormap("gray")
c=colorbar;
c.Label.String = 'Band Contrast /arb units';
hold on
scatter(x_transREF,y_transREF)
for i=1:size(x_transREF,1)
    text(x_transREF(i)+0.01*max(X(:)), y_transREF(i)+0.01*max(Y(:)),num2str(i));   %put up the text at the given locations
end
title('References Selected On EBSD Map')
xlabel('\mum')
ylabel('\mum')
axis image
figname=['refselect_EBSD figure'];
saveas(gcf,fullfile(resultsdir, figname),'png')
close all
%% transformation calculation

movingPoints  = [x_og    y_og   ];
fixedPoints = [x_transREF y_transREF];

transtype='affine'; %select what type of transformation you want % HC may want this to be different from corresponding the EBSD patterns
tform = fitgeotrans(movingPoints,fixedPoints,transtype); % HC needs updating
% transform the data
A = tform.T;

%transform the x and y coordinates of the nanoindentation data
%using the affine transformation
[locebsdcorr(:,1),locebsdcorr(:,2)] = transformPointsForward(tform,X(:),Y(:));
locebsdcorrafter=locebsdcorr;

%% Plot transformed points back for checking
%plot it again with the points on after for helpful guide, and save.

figure
hplot=contourf(x_ebsd,y_ebsd,Data_InputMap_after.RadonQuality,45,'LineColor','None');
colormap("gray")
c=colorbar;
c.Label.String = 'Band Contrast /arb units';
hold on
scatter(locebsdcorrafter(:,1),locebsdcorrafter(:,2),'rx')
title('Transformed nanoindentation points on after PQ Map')
figname=['Transformed_nanoindentation_after  figure'];
saveas(gcf,fullfile(resultsdir, figname),'png')

%plot it again with the points on before for helpful guide, and save.
[locebsdcorrbefore(:,1),locebsdcorrbefore(:,2)] = transformPointsForward(tformxebsd,locebsdcorrafter(:),locebsdcorrafter(:));
figure
hplot=contourf(x_ebsd,y_ebsd,Data_InputMap_after.RadonQuality,45,'LineColor','None');
colormap("gray")
c=colorbar;
c.Label.String = 'Band Contrast /arb units';
hold on
scatter(locebsdcorrbefore(:,1),locebsdcorrbefore(:,2),'rx')
title('Transformed nanoindentation points on before PQ Map')
axis image
figname=['Transformed_nanoindentation_before  figure'];
saveas(gcf,fullfile(resultsdir, figname),'png')

%% More sampling points
Xvector=X(:,1);
Yvector=Y(1,:);

Xspacing=X(2,1)-X(1,1);
Yspacing=Y(1,2)-Y(1,1);
stepsize=0.25;% get this from array
Noofsamplesxhalf=((Xspacing/2)/stepsize);
Noofsamplesyhalf=((Yspacing/2)/stepsize);

x_ebsd_raw_after_vector=Data_InputMap_after.X_axis;
x_ebsd_flipped_after_vector=flipud(x_ebsd_raw_after_vector);
y_ebsd_raw_after_vector=Data_InputMap_after.Y_axis;
GNDtotalafter_ebsd=GND_after.total';

[x_ebsd_after,y_ebsd_after]=ndgrid(x_ebsd_flipped_after_vector,y_ebsd_raw_after_vector);

x_ebsd_raw_before_vector=Data_InputMap_before.X_axis;
x_ebsd_flipped_before_vector=flipud(x_ebsd_raw_before_vector);
y_ebsd_raw_before_vector=Data_InputMap_before.Y_axis;
GNDtotalbefore_ebsd=GND_before.total';

[x_ebsd_before,y_ebsd_before]=ndgrid(x_ebsd_flipped_before_vector,y_ebsd_raw_before_vector);

GNDtotalbefore_ebsdinterp = griddedInterpolant(x_ebsd_before,y_ebsd_before,flipud(GNDtotalbefore_ebsd),'nearest');
GNDtotalafter_ebsdinterp = griddedInterpolant(x_ebsd_after,y_ebsd_after,flipud(GNDtotalafter_ebsd),'nearest');

GNDtotalbefore_ebsdinterp=flipud(GNDtotalbefore_ebsdinterp);
GNDtotalbefore_ebsdinterp=GNDtotalbefore_ebsdinterp';

GNDtotalafter_ebsdinterp=flipud(GNDtotalafter_ebsdinterp);
GNDtotalafter_ebsdinterp=GNDtotalafter_ebsdinterp';

for noofxpoints=1:numel(Xvector)
    Xvectorlower=[];
for sampling=0:Noofsamplesxhalf
    indexproblem=sampling+1;
    Xvectorlower(indexproblem)=Xvector(noofxpoints)-(stepsize*sampling);
end
    Xvectorlower=fliplr(Xvectorlower);
    Xvectorlowersize=numel(Xvectorlower);
    
for sampling=1:Noofsamplesxhalf
    Xvectorlower(Xvectorlowersize+sampling)=Xvector(noofxpoints)+(stepsize*sampling);
end
XNewcoordinate(noofxpoints,:)=Xvectorlower;
end


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

for noofxpoints=1:numel(Xvector)
    Xboxvector=XNewcoordinate(noofxpoints,:);
    for noofypoints=1:numel(Yvector)
Yboxvector=YNewcoordinate(noofypoints,:);


[xboxmatrix,yboxmatrix]=ndgrid(Xboxvector,Yboxvector);

[locebsdcorrboxafter(:,1),locebsdcorrboxafter(:,2)] = transformPointsForward(tform,xboxmatrix(:),yboxmatrix(:));


 GNDtotalafter_ebsdnew_box = GNDtotalafter_ebsdinterp(locebsdcorrboxafter(:,1), locebsdcorrboxafter(:,2));
 GNDtotalafter_ebsdnew_box(GNDtotalafter_ebsdnew_box == 0)=NaN;
 GNDtotalafter_ebsdnew_boxaverage(noofxpoints,noofypoints)=nanmean(GNDtotalafter_ebsdnew_box);
 
 [locebsdcorrboxbefore(:,1),locebsdcorrboxbefore(:,2)] = transformPointsInverse(tformxebsd,locebsdcorrboxafter(:,1),locebsdcorrboxafter(:,2));
 
 GNDtotalbefore_ebsdnew_box = GNDtotalbefore_ebsdinterp(locebsdcorrboxbefore(:,1), locebsdcorrboxbefore(:,2));
 GNDtotalbefore_ebsdnew_box(GNDtotalbefore_ebsdnew_box == 0)=NaN;
 GNDtotalbefore_ebsdnew_boxaverage(noofxpoints,noofypoints)=nanmean(GNDtotalbefore_ebsdnew_box);
    end
end
%% interpolation section
%using griddedInterpolant, setup a variable that contains all the
%information, and then probe at the locations of the distorted x and y
%locations from the nanoindentation data
% 'nearest' is used as we don't want to invent orientations

 phi1interp = griddedInterpolant(x_ebsd,y_ebsd,phi1,'nearest');
 phi1new = phi1interp(locebsdcorr(:,1), locebsdcorr(:,2));
 
 Phiinterp = griddedInterpolant(x_ebsd,y_ebsd,Phi,'nearest');
 Phinew = Phiinterp(locebsdcorr(:,1), locebsdcorr(:,2));
  
 phi2interp = griddedInterpolant(x_ebsd,y_ebsd,phi2,'nearest');
 phi2new = phi2interp(locebsdcorr(:,1), locebsdcorr(:,2));
 
 phaseinterp = griddedInterpolant(x_ebsd,y_ebsd,phase,'nearest');
 phasenew = phaseinterp(locebsdcorr(:,1), locebsdcorr(:,2));
 
 BCebsdinterp = griddedInterpolant(x_ebsd,y_ebsd,BCebsd,'nearest');
 BCebsdnew = BCebsdinterp(locebsdcorr(:,1), locebsdcorr(:,2));
 

 
%Gridify into a matrix 
phi1newG = f_gridify_vector(phi1new,size(X,1),size(Y,2))';
PhinewG = f_gridify_vector(Phinew,size(X,1),size(Y,2))';
phi2newG = f_gridify_vector(phi2new,size(X,1),size(Y,2))';
phasenewG = f_gridify_vector(phasenew,size(X,1),size(Y,2))';
BCebsdnewG = f_gridify_vector(BCebsdnew,size(X,1),size(Y,2))';


%output into a class datastack with reasonable names
datastack.X     = X;                 %X position
datastack.Y     = Y;                 %Y position
datastack.S     = fullres(:,:,1);   %Surface displacement
datastack.D     = fullres(:,:,2);   %Depth
datastack.L     = fullres(:,:,3);   %Load
datastack.M     = fullres(:,:,4);   %Modulus
datastack.St    = fullres(:,:,5);   %Stiffness^2/Load
datastack.H     = fullres(:,:,6);   %Hardness
datastack.phi1  = phi1newG;          %phi1  (ebsd)
datastack.Phi   = PhinewG;           %Phi   (ebsd)
datastack.phi2  = phi2newG;          %phi2  (ebsd)
datastack.phase = phasenewG;         %phase (ebsd)
datastack.BCebsd= BCebsdnewG;        %Band Contrast (ebsd)
datastack.GNDtotalbefore = GNDtotalafter_ebsdnew_boxaverage;       %GNDtotal average after (xebsd) 
datastack.GNDtotalafter = GNDtotalbefore_ebsdnew_boxaverage; 

end
