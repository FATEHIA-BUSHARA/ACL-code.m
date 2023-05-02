function [Diameter,Volume, Theta]= ACL_FeatEx(V, Threshold)
%% FEATURES EXTRACTION: Volume, Orientation and Diameter

% Convert the intensity volume into a 3-D binary volume.
% Sensitivity factor for adaptive thresholding,
% specified as a value in the range [0, 1]. 
% A high sensitivity value leads to thresholding more pixels as foreground, 
% at the risk of including some background pixels.
%% Eliminate small objects volume with different threshold values dx
J = imbinarize(V,'adaptive','Sensitivity',0.1);
cc = bwconncomp(J,18); 
stats = regionprops3(cc, 'Volume'); 
idx = find([stats.Volume] > Threshold); 
J = ismember(labelmatrix(cc), idx);

% Extract volume, diameter and orientation and saved in excel sheet.
% s = regionprops3(J,"Volume",'EquivDiameter',"PrincipalAxisLength",'Orientation','SurfaceArea');
s = regionprops3(J,"Volume","PrincipalAxisLength",'Orientation');
D=s.PrincipalAxisLength;
Diameter = mean(D(:,2:3),2);
% Diameter = mean(s.PrincipalAxisLength,2);
Diameter =  Diameter* 9e-3; 
Volume = (s.Volume) * 0.729e-6;  % 1 voxel = 0.729e-6 mm^3
Theta_all = s.Orientation;
Theta = Theta_all(:,2);
% Theta = 90 - Theta;

% My_Meas = [Diameter, volume,Theta];

% filename = fullfile('D:\mesearment t27\Diffusion' ,...
%      sprintf('Diffusion_t0_27 .xlsx',n));
%  writetable( s, filename,'Sheet',num2str )
%  
% FileName = fullfile('D:\mesearment t27\Diffusion',...
%         sprintf('Diffusion_Meas_t0_27 ' , num2str (dx(n)), '.xlsx'));
% 
% xlswrite(FileName, My_Meas, n)

% % % % save the binary volume into slices
% for sliceIndex = 1: size (J, 3)
%     thisSlice = J (:,:, sliceIndex);
%     FileName = fullfile('D:\T027\TK027_t0_1750-1292 in VOI\New Diffusion\BW',...
%         sprintf('TK027_t0_Diffusion BW # %d.jpg', sliceIndex));
%     imwrite (thisSlice, FileName);
% end
% clear thisSlice


% 'Volume' Count of the actual number of 'on' voxels in the region,(scalar)
% Volume represents the measure of the number of voxels 
% in the regions within the volumetric binary image.

% 'EquivDiameter'	Diameter of a sphere with the same volume as the region, 
% as a scalar. Computed as (6*Volume/pi)^(1/3).

% 'PrincipalAxisLength'	Length (in voxels) of the major axes of the
% ellipsoid that have the same normalized second central moments
% as the region, returned as 1-by-3 vector.

% 'Orientation' Euler angles, returned as a 1-by-3 vector.
% The angles are based on the right-hand rule.
% regionprops3 interprets the angles by looking at the origin along the 
% x-, y-, and z-axis representing roll, pitch, and yaw respectively.
% A positive angle represents a rotation in the counterclockwise direction

% 'SurfaceArea'	Distance around the boundary of the region, (scalar).