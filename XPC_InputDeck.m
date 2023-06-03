%Input deck for multi-dimensional spatially resolved analysis of EBSD,
%Express, and other maps. CMM 2019. %HC version that has been editted this
%is XPCorrelate Electric Boogaloo.
%% 1 Sorting dependencies
clear all
close all
clc
home
%adapted from MTEX - load mtex
%add the mtex path - this will need to be changed for the pc you run on!
%Mtex_Loc='Z:\2022 Pt IIs\Hannah Cole\mtex-5.8.0\startup_mtex';
Mtex_Loc='Z:\2022 Pt IIs\Hannah Cole\mtex-5.8.0\startup_mtex';
run(Mtex_Loc);

addpath(genpath('src'))
addpath(genpath('external'))

%% 2 User Inputs
%This script runs entirely within the location of the results - it does not
%save any file in the scripts location. 

filepath    = 'Z:\2022 Pt IIs\Hannah Cole\nanoindentation\testingcode\x65chargedlocation\';
%filepath    = %'Z:\2022 Pt IIs\Hannah Cole\nanoindentation\testingcode\XPcorrelatelocation\';  %Location of files
Premieroutputname = 'Workspace_Output_2_postcal.mat';
ebsdname    = 'sample61702withindents.ctf'; %name of the ebsd with extension after indents
xebsdnamebefore ='HRsample2charged_XEBSD.mat';
xebsdnameafter ='HRsample2charged_XEBSD.mat'; %name of xebsd .mat output


epmaq       =  0; %epma analysis? 1=yes, 0=no
epmaname    = 'Map 1_ds1_O  Ka_O  Ka (Sp 1)_item2.tiff'; %name of epma file if needed 
epmabsename = 'Map 1_ds1Vs1 BSE Z_BSE Z_item1'; %name of epma BSE file if needed 


gdcalcq     = 1; %grain boundary analysis? 1=yes, 0=no
microscope  ='merlin';%what microscope was the ebsd map taken on? 'evo', 'merlin', 'xbeam', or 'tescan'
primphase   ='default';%what's the primary phase in the ebsd map? default uses 1st material phase
%ebsd analysis questions
EBSDREFq    = 3;%0 if just one phase: uses IPF for EBSD reference selection
                %1 if multiphase: uses phase map
                %2 if you want to use the pattern quality map
                %3 if you want to use the GND plot
EPMArefq    = 1;%0 uses the chemical map for referencing
                %1 if present, uses the BSE map for referencing

hexmat      = 0; %if it's hexagonal, make reflect Phi values above 90 about 90deg
saveasebsdq = 1; %activate the writedata to a file that can be read by mtex
saveasfigq  = 1; %save images as .fig as well as .png?
resolution  = ['-r' num2str(600)]; %resolution for .pngs
%% 3 Loading and some manipulation

% idx = strfind(ebsdname,'.'); %get the extension and load ebsd
% if strcmp(ebsdname(idx+1:end),'h5') %if it is an h5 file
%     ebsd=loadEBSD_h5(fullfile(filepath,ebsdname));
% elseif strcmp(ebsdname(idx+1:end),'ctf') %if it is a ctf file
%     ebsd=EBSD.load(fullfile(filepath,ebsdname),ebsdname(idx+1:end),'convertSpatial2EulerReferenceFrame');
% end

load(fullfile(filepath,Premieroutputname));%load the nanoindentation data

% loading the HR-EBSD data
input_file_before=[filepath xebsdnamebefore];
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

input_file_after=[filepath xebsdnameafter]; %file path for inputting before data
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


currdate=datestr(datetime); %use the date to distinguish between analysis runs
currdate=currdate(1:11);

resultsdir=[filepath,'Results_',currdate]; % resultsdirold(idxres(1)+1:idxres(2))

%results directory creation 
if ~exist(resultsdir, 'dir')
   mkdir(resultsdir)
end

%Cropping of nanoindentation data if need be:
cropytop=0;
cropybot=0;
cropxright=0; 
cropxleft=0;
X=X(cropxleft+1:end-cropxright,cropytop+1:end-cropybot);
Y=Y(cropxleft+1:end-cropxright,cropytop+1:end-cropybot);
fullres=fullres(cropxleft+1:end-cropxright,cropytop+1:end-cropybot,:);


%% 4 Running the EBSD registration

%% GND function
ebsd=1; % hard code to fix the EBSD not being loaded
fractionofspacing=0.6; % The size of the sampling box as a fraction of the spacing between nanoindentation indents. (It will round the number of sampling points to an integar)
stepsize=0.25;%step size of HR-EBSD data this is coded to be the same for both before and after
[datastack,tform,tformxebsd] = f_fixGNDdistortion(ebsd,X,Y,fullres,Data_InputMap_before,Data_InputMap_after,Map_PH2_before,Map_PH2_after,GND_before,GND_after,microscope,primphase, EBSDREFq,resultsdir,fractionofspacing,stepsize);

%%

[datastack] = f_fixEBSDdistortion(ebsd,X,Y,fullres,Data_InputMap_before,Data_InputMap_after,GND_before,GND_after,microscope,primphase, EBSDREFq,resultsdir);
%Datastack structure: 
% datastack.X               %X position
% datastack.Y               %Y position
% datastack.S               %Surface displacement %HC removed
% datastack.D               %Depth
% datastack.L               %Load
% datastack.M               %Modulus
% datastack.St              %Stiffness^2/Load
% datastack.H               %Hardness
% datastack.phi1            %phi1  (ebsd)
% datastack.Phi             %Phi   (ebsd)
% datastack.phi2            %phi2  (ebsd)
% datastack.phase           %phase (ebsd)
% datastack.BCebsd          %Band Contrast (ebsd)
% datastack.GNDtotal        %GNDtotal (xebsd)

%CLEANING:
datastack.S(datastack.S>1e100)=0;     
datastack.D(datastack.D>1e100)=0;  
datastack.L(datastack.L>1e100)=0;    
datastack.M(datastack.M>1e100)=0;     
datastack.St(datastack.St>1e100)=0; 
datastack.H(datastack.H>1e100)=0;
%housekeeping:
X=datastack.X;
Y=datastack.Y;

if hexmat==1
    %first fix the phi to be between 0 and 90 FOR HCP
    datastack.Phirefl=datastack.Phi;
    datastack.Phirefl(datastack.Phirefl>(pi/2))=pi-datastack.Phirefl(datastack.Phirefl>(pi/2));
else
    datastack.Phirefl=datastack.Phi; %if it isn't hexmat, still make a variable for the plotters
end
%% 5 Saving figures
v_plotfig(datastack,resultsdir,ebsdname,xebsdnamebefore,xebsdnameafter,saveasfigq);

if EBSDREFq==1
    v_multiphaseplotter(datastack,resultsdir,ebsdname)
end


%% 6 EPMA input
if epmaq==1
    cropq=0; %DO WE NEED TO CROP THE DATA? 1=yes 0=no
    %import the epma data (currently an image)
    [C, XC, YC] = f_loadEPMA(epmaname, filepath, epmabsename);
    
    %alignment
    [datastack]=f_fixEPMAdistortion(datastack,C,XC,YC,cropq, EPMArefq);
    
    %some smoothing
    datastack.EPMAOS = smoothdata(datastack.EPMAO,2,'gaussian',5.5);
    datastack.EPMAOS = smoothdata(datastack.EPMAOS,1,'gaussian',5.5);
    datastack.EPMAO = datastack.EPMAOS;
    
    v_plotfigepma(datastack,resultsdir, ebsdname)
end
%% 7 save as ebsd and xebsd
if saveasebsdq==1
    disp('Saving registered EBSD file')
    f_writeEBSDdata(fullfile(resultsdir, [ebsdname(1:(max(size(ebsdname)-4))) 'corrected' currdate]),datastack)
end 


%% 8 grain distance analysis
if gdcalcq==1
    disp('Running grain boundary distance analysis')
    datastack=f_dist2grainb(resultsdir, ebsdname,datastack,currdate,ebsd);
end

%% 9 Save things
close all
save([fullfile(resultsdir,[filename(1:length(filename)-4) '_XPCorrelate_results' currdate]) '.mat']);
disp('Analysis complete')
 