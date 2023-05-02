close all
clear
clc
%% 
addpath 'G:\D\EXP2 matlab function'
tic
% Given: edges, Threshold, lamda3 and Force

edges = 0:10:90;

Threshold = 300;

%%%%%%%%%%%%%%%%%%%%%%%%%%%% Healthy Variable %%%%%%%%%%%%%%%%%%%%%%%%%%%%
lamda3 = [1.1111,1.1222,1.1334,1.1448 ,1.1562,1.1678,1.1795,1.1913,1.2032];
Force = [8.4,14.2, 20.6, 32.0, 45.0, 73.5, 116.9, 156.0, 234.0];

%%%%%%%%%%%%%%%%%%%%%%%%%%  pathological Variable %%%%%%%%%%%%%%%%%%%%%%%%%
% lamda3 = [1.4206 1.4348 1.4491 1.4636  1.4782 1.4930 1.5080 1.5230 1.5383];
% Force = [2.3 2.5 3.0 3.4 4.2 5.1 6.3 7.8 9.2];

C = length(lamda3);
S = zeros(C,1);    % SIGMA template at each stress level
B = zeros(C,1) ;   % BITA template at each stress level
Result = zeros(C,5);
TV = zeros(C,1);   % Total Volume of objects template at each stress level
ElasticModulus= nan(C, 1); % to fill the unused spaces with NAN

% myDir = 'D:\T027\New folder'; 
myDir ='G:\Fatehia\T027\New folder Healthy\Hessian Original';% healthy
% myDir ='G:\Fatehia\T024\New folder'; % the direction of the volume images
ext_img = '*.nii';              %file extension of interest
dircontent = dir(fullfile(myDir, ext_img));
assert(numel(dircontent) > 0, 'No file was found. Check that the path is correct');
my_img = struct('img', cell(size(dircontent)));  %preallocation of the structure
 path = 'G:\Fatehia\ACL results 11.1.2023\Healthy';%to save the Healty/pathological figuers
% path = 'D:\T027\New folder\Healthy Figures' ; %to save the Healthy figuers
% path = 'D:\T024\New folder\Pathological Figures' ; %to save the pathological figuers

for fileidx = 1:numel(dircontent)
% we get the image from the folder at Path determine already (myDir)line 27  
my_img(fileidx).img = niftiread(fullfile(myDir, dircontent(fileidx).name));

    [Diameter,Volume,Theta]= ACL_FeatEx(my_img(fileidx).img, Threshold);
    Data = [Diameter,Volume,Theta];
    Data= abs(Data);
    Data = sortrows(Data,3); % to order the orientation
    Diameter = Data(:,1);  % in mm
    Volume  = Data(:,2);  % in mm cube
    Theta  = 90-Data(:,3);   % in Degree
    
    TV(fileidx) = sum(Volume);
    
%%%%%%%%%%%%%%%%%%%%  PlotDiameter Distribution (mm)  %%%%%%%%%%%%%%%%%%%%%
f = figure;
h= histogram(Diameter,20 ,'Normalization','probability' );
h.BinEdges = 0:.05:.6;
xlabel('Diameter (mm)')
ylabel('Relative Frequency')
% title(['Healthy T',num2str(fileidx-1), '- Diameter Distribution (mm)'])
title(['Pathological at T',num2str(fileidx-1), ' - Diameter Distribution (mm)'])

saveas (f, fullfile (path, [ 'Diameter T', num2str(fileidx-1),'.jpg']));
 
%%%%%%%%%%%%%%%%%%%%%%% Plot Orientation in degree  %%%%%%%%%%%%%%%%%%%%%%%
f1 = figure;
h1 = histogram(Theta,10 ,'Normalization','probability' );
h1.BinEdges = 0:10:90;
xlabel('Orientation ')
ylabel('Relative Frequency')
% title(['Healthy T',num2str(fileidx-1),'- Orientation Distribution (Degree)'])
title(['Pathological T',num2str(fileidx-1),' - Orientation Distribution (Degree)'])

saveas (f1, fullfile (path, [ 'Orientation T', num2str(fileidx-1),'.jpg']));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[Vf,Sigma,Beta,Results] = ACL_Strain(Diameter,Volume,Theta,lamda3(fileidx),Force(fileidx),edges);
S(fileidx) = Sigma;
B(fileidx) = Beta;
Result(fileidx,:) = Results;

%%%%%%%%%%%%%%%%%%%%% Plot Volume fraction in degree %%%%%%%%%%%%%%%%%%%%%
f2 = figure;
xy1 = Vf.* Theta; 
h11 = histogram(xy1,10 ,'Normalization','probability' );
h11.BinEdges = 0:10:90;
% edges = 0:10:90;  % 0:0.1745:1.6; = deg2rad (BinEdges) ,BinWidth = 0.1745
% h11.BinEdges = deg2rad(edges);
xlabel('Volume fraction')
ylabel('Relative Frequency')
% title(['Healthy T',num2str(fileidx-1),'- Volume fraction'])    
title(['Pathological T',num2str(fileidx-1),' - Volume fraction'])

saveas (f1, fullfile (path, [ 'Weighted Orientation T', num2str(fileidx-1),'.jpg']));


end

%% Calculate the Elastic Module
% S(lamda3)= Ef* B(lamda3)
% Perform simple linear regression using the \ operator.

Ef= B\S; 
sCalc = Ef*B;
%%
%%%%%%%%%%%%%%%%%%%%%  Plot and save Stress & Beta  %%%%%%%%%%%%%%%%%%%%%
f3 = figure;
scatter(B,S)
hold on
plot(B,sCalc)
xlabel('Beta')
ylabel('Strees')
title('Linear Regression Relation Between Stress & Beta')
grid on
saveas (f3, fullfile (path, 'Elastic Modulus .jpg'));

% Find the better fit of the two fits by comparing values of R2. 
Rsq1 = 1 - sum((S - sCalc).^2)/sum((S - mean(S)).^2);

%%%%%%%%%%%%%%%%%%%%%  Plot and save Stress & Lamda  %%%%%%%%%%%%%%%%%%%%%
f4 = figure;
scatter(lamda3,S)
hold on
plot(lamda3,sCalc)
xlabel(' Fibers axial stretch-ratio,')
ylabel('Tissue Stress')
title('Effect of Microstructural Organization on Tissue Stres')
grid on
saveas (f4, fullfile (path, 'Stress & Lamda .jpg'));
%%
%%%%%%%%%%%%%%%%%%%%%%% Save Results in Excele File %%%%%%%%%%%%%%%%%%%%%%%

StressLevel = {'T0';'T1';'T2';'T3';'T4';'T5';'T6';'T7';'T8'};
TotalVolume = TV;
Sigma = S;
Beta =  B;
Area = Result(:,1);
Entropy = Result(:,2);
Energy = Result(:,3);
MeanDiameter = Result(:,4);
VarDiameter= Result(:,5);
ElasticModulus(1) = Ef;
T = table(StressLevel,TotalVolume,Sigma,Beta,Area,Entropy,Energy,MeanDiameter,VarDiameter,ElasticModulus);
% filename = fullfile(path,'Healty Data .xlsx'); 
filename = fullfile(path,'Pathological Data .xlsx'); 
     writetable( T, filename)

toc