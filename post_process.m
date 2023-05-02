function [dilateE]= post_process(V)

% Post Processing of the hessian output
%As the opening removes small and isolated objects and preserves the background

I0 = imopen(V,strel3d(10));
% volumeViewer(I0)  

%-------- Binarize the processed volume by intensity thresholding---------
BW = imbinarize(V,'adaptive','Sensitivity',0.1);
BW=(I0 > 0.05);
% volumeViewer(BW)

%% morphological 

% To separate touching fibers, an  erosion method was used followed by
% dilation to repair fascicles which may have been split into different
% connected components during the binarization process

%----------------------morphological erosion-------------------------%

%  Remove structures having a radius less than 10 voxels by erosion 
%  Line3D (creates a 3D Lined-shaped structuring element).

se = Line3D(20,[ 0 0   1 ] ,1,1,1); %3D line structuring element
erodeBW = imerode(BW, se);
% volumeViewer(erodeBW)

%----------------------morphological dilation-------------------------%

% % A morphological dilation operation was used to connect nearby ridges. 
% % The radius of the sphere is 10 voxels. 
%  STREL3D (creates a 3D sphere-shaped structuring element).

dilateE = imdilate(erodeBW, strel3d(10)); 
% volumeViewer(dilateE)

end
