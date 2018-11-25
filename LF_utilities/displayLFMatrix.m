% Displays a LightField Matrix as a video
%Input :
% - LFMatrix : Light field given as a matrix where each column is a vectorized view.
% - imgSize : (vector of 2 elements) size of the views ([YResolution, XResolution]).
% - nU, nV : number of views in each angular dimension.
% - Scan : (optional) integer value defining the scanning order for the video :
%               0 : Zscan (views are displayed as they are ordered in the matrix).
%               1 : ZigZag scan (default).
%               2 : Scan vertically, horizontally and diagonally along the lightfield (does not show all the views).
% - NormBounds : (optional) parameter defining the minimum and maximum pixel values to display.
%                   - default (parameter not given) : display the range [0 1] for logicals / [0 255] otherwise.
%                   - [mini, maxi] : use mini and maxi as the minimum and maximum values.
%                   - false or empty array : Adjust to minimum and maximum values in the matrix.
%                   - true : Adjust to the minimum and maximum values for each frame.
% - colorMap : (optional, only for grayscale data) colormap (nx3 matrix) defining color display at each intensity level.
%              Set to empty to ignore.
%
%Output :
% - Video : the video given as an array concatenating the frames (suitable for matlab's implay function).
% - path : indices of the frames in the ordered they are displayed.

function [Video,path] = displayLFMatrix (LFMatrix, imgSize, nU, nV, Scan, NormBounds, colormap)

numImages=nU*nV;



path=1:numImages;
if(exist('Scan','var') && Scan==1)
    for idx=2:2:nU
        path(1+nV*(idx-1):nV*idx)=nV*idx:-1:1+nV*(idx-1);
    end
elseif(Scan==2)
    path=[[1:nV] ...
     [nV:nV:numImages]...
     [numImages:-1:numImages-nV+1]...
     [numImages-nV+1:-nV:1]...
     [1:nV+1:numImages]];
end

frame_id=0;

nChan = size(LFMatrix,3);
UseColorMap = nChan==1 && exist('colormap','var') && size(colormap,1)>0 && size(colormap,2)==3;

if(UseColorMap)
    maxDisp = 1;
    Video = zeros(imgSize(1),imgSize(2), nChan, length(path));
else
    maxDisp = 255;
    Video = zeros(imgSize(1),imgSize(2), nChan, length(path),'uint8');
end

if(islogical(LFMatrix)), minVal=0; maxVal=1;
else minVal=0; maxVal=255;
end
AdjustPerFrame = false;

if(exist('NormBounds','var'))
    if(length(NormBounds)==2)
        minVal=NormBounds(1); maxVal=NormBounds(2);
    elseif(isempty(NormBounds) || (isscalar(NormBounds) && NormBounds==false) )
        minVal = min(LFMatrix(:)); maxVal = max(LFMatrix(:));
    elseif(isscalar(NormBounds) && NormBounds==true)%Adjust bounds per frame
        AdjustPerFrame = true;
    end
end

if(AdjustPerFrame)
    for idx=path
        frame_id = frame_id+1;
        minVal = min(min(LFMatrix(:,idx,:)));
        maxVal = max(max(LFMatrix(:,idx,:)));
        mult = maxDisp/(maxVal-minVal);
        Video(:,:,:,frame_id) = reshape((double(LFMatrix(:,idx,:))-minVal)*mult,[imgSize,nChan]);
    end
else
mult = maxDisp/(maxVal-minVal);
    for idx=path
        frame_id = frame_id+1;
        Video(:,:,:,frame_id) = reshape((double(LFMatrix(:,idx,:))-minVal)*mult,[imgSize,nChan]);
    end
end

if(UseColorMap)
    R = colormap(:,1);
    G = colormap(:,2);
    B = colormap(:,3);
    Video = round(1 + Video * (size(colormap,1)-1));
    Video = cat(3,R(Video),G(Video),B(Video));
end

%{
if(~exist('luma_only','var'))
    luma_only = true;
end

if(luma_only)
    Video = zeros(imgSize(1),imgSize(2), length(path));
    for idx=path
        frame_id = frame_id+1;
        Video(:,:,frame_id) = reshape(LFMatrix(:,idx),imgSize);
    end
else
    Video = zeros(imgSize(1),imgSize(2), nChan, length(path));
    for idx=path
        frame_id = frame_id+1;
        Video(:,:,:,frame_id) = reshape(LFMatrix(:,3*idx-2:3*idx),[imgSize,3]);
    end
end
%}

implay(Video,16);