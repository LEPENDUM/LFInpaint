rootDir = fileparts(mfilename('fullpath'));
DataRoot=[rootDir '/LFData/PNG/'];
LabelsDir='GT/';
MaskDir='Mask/';

clear baselineUV

%%%%%%%%%%%%%%%%%%%%%%    LF parameters     %%%%%%%%%%%%%%%%%%%%%
%{
LFName='stilllife';%'BeeCrop';%
uRange=[0:8];%uRange=[0:0];%
vRange=[0:8];
objLblColor = [255 255 0];

dilateMask=1;%dilateMask=0;%
ROI = [37,45,484,386];%[];%[80,80,400,380];%    %[x0,y0,x1,y1]  / set to [] for full image.
rescale = 1;
luma_only = false;%luma_only=false;%
MaskFromLabel = false;
baselineUV = [-1,1];
%}
%{
LFName='butterfly';%'BeeCrop';%
uRange=[0:8];%uRange=[0:0];%
vRange=[0:8];
objLblColor = [0 0 255];
dilateMask = 2;
ROI = [262,192,654,708];%[262,390,654,708];%       %set to [] for full image.
rescale = 1;
luma_only = false;
MaskFromLabel = false;
baselineUV = [-1,1];
%}
%{
LFName='budha';
uRange=[0:8];%uRange=[0:0];%
vRange=[0:8];
objLblColor = [0 0 255];
dilateMask = 2;
ROI = [];%[51,1,558,768];       %set to [] for full image.
rescale = 1;
luma_only = false;
MaskFromLabel = false;
baselineUV = [-1,1];
%}
%%{
LFName='F01_TotoroWaterfall';%'F01_TotoroAlley'/'F01_TotoroWaterfall'/'F01_TapeMeasureAWB'/'F01_BeersAWB'/'F01_Guitar'/
uRange=[2:8];
vRange=[2:8];
dilateMask=0;
ROI = [2,2,378,378];%[146,131,333,377];%[49,145,293,377];%[];%[1,212,379,379];%    %[x0,y0,x1,y1]  / set to [] for full image.
rescale = 1;%2/(3*sqrt(3));%
luma_only = false;
MaskFromLabel = false;
%}
%{
LFName='F01_TotoroAlleyKAIST';%'F01_BeersKAIST'/'F01_TotoroAlleyKAIST'/'F01_TotoroWaterfallKAIST'/'F01_GuitarKAIST'/
uRange=[0:6];
vRange=[0:6];
dilateMask=0;
ROI = [];%[1,212,379,379];%    %[x0,y0,x1,y1]  / set to [] for full image.
rescale = 1;
luma_only = false;
MaskFromLabel = false;
%}
%{
LFName='Illum_Figurines';%'Illum_Field'/'Illum_TinyMoon'/'Illum_Bee2'/'Illum_Figurines'/'Illum_Bumblebee'/'Illum_Fruits'/'EPFL_PalaisDuLuxembourg'/'EPFL_Chain-link_Fence_2'/
uRange=[2:12];
vRange=[2:12];
dilateMask=0;
%ROI : Illum_Figurines=[96,2,527,433] / Illum_Field=[98,3,528,432] / Illum_TinyMoon=[5,5,434,612] / Illum_Bee2=[3,3,623,432] / Illum_Bumblebee=[192,128,421,357]
ROI = [96,2,527,433];%[192,128,421,357];%[98,3,528,432];%    %[x0,y0,x1,y1]  / set to [] for full image.
rescale = 1;%2/(3*sqrt(3));
luma_only = false;
MaskFromLabel = false;
%baselineUV = [-1,1];
%}
%{
LFName='Illum_Cactus';
uRange=[3:11];
vRange=[3:11];
dilateMask=0;
ROI = [];    %[x0,y0,x1,y1]  / set to [] for full image.
rescale = 1;%2/(3*sqrt(3));
luma_only = false;
MaskFromLabel = false;
baselineUV = [-1,1];
%}

%{
LFName='lego';
uRange=[0:4];
vRange=[0:4];
dilateMask=0;
ROI = [1,384,430,960];    %[x0,y0,x1,y1]  / set to [] for full image.
rescale = 1;%2/(3*sqrt(3));
luma_only = false;
MaskFromLabel = false;
%}

SingleFormat = false;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if(~exist('baselineUV','var'))
    baselineUV = [];
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
nV = length(vRange);
nU = length(uRange);
numImages = nU*nV;
centralView=ceil(numImages/2);
uvImgRefDefault = [min(uRange)+floor((max(uRange)-min(uRange))/2), min(vRange)+floor((max(vRange)-min(vRange))/2)];
idImgRefDefault = sum(uRange<uvImgRefDefault(1))*nV + sum(vRange<=uvImgRefDefault(2));

%Search available mask or inpainted view to determine its index (imgRefId)
fInp = dir([DataRoot LFName '/']);
if(~MaskFromLabel)
    objLblColor=[255 255 255];
    fMask = dir([DataRoot LFName '/' MaskDir]);
    fMask = regexp({fMask.name},'^\d+_\d+[.]png$','match');
    fMask = [fMask{:}];
    if(isempty(fMask))
        uvImgRef = uvImgRefDefault;
        idImgRef = idImgRefDefault;
        HasInpainedView = false;
        HasMask = false;
    else
        HasMask = true;
        uvImgRef=[NaN NaN];
        i=1;
        while(~(ismember(uvImgRef(1), uRange) && ismember(uvImgRef(2), vRange)) && i<=length(fMask))
            uv = regexp(fMask{i},'\d+','match');
            uvImgRef = [str2num(uv{2}) str2num(uv{1})];
            i=i+1;
        end
        if(~(ismember(uvImgRef(1), uRange) && ismember(uvImgRef(2), vRange)))
            error('No Mask found for any view of the light field in the specified range.');
        end
        idImgRef = sum(uRange<uvImgRef(1))*nV + sum(vRange<=uvImgRef(2));
        fInp = regexp({fInp.name},['^0*' num2str(uvImgRef(2)) '_0*' num2str(uvImgRef(1)) '_inp[.]png$'],'match');
        fInp = [fInp{:}];
        HasInpainedView = ~isempty(fInp) ;
        %ismember([num2str(uvImgRef(2)) '_' num2str(uvImgRef(1)) '_inp.png'],{fInp.name});
    end
else
    fInp = regexp({fInp.name},'^\d+_\d+_inp[.]png$','match');
    fInp = [fInp{:}];
    if(isempty(fInp))
        HasInpainedView = false;
        uvImgRef = uvImgRefDefault;
        idImgRef = idImgRefDefault;
    else
        uvImgRef = [NaN NaN];
        i=1;
        while(~(ismember(uvImgRef(1), uRange) && ismember(uvImgRef(2), vRange)) && i<=length(fInp))
            uv = regexp(fInp{i},'\d+','match');
            uvImgRef = [str2num(uv{2}) str2num(uv{1})];
            i=i+1;
        end
        HasInpainedView = ismember(uvImgRef(1), uRange) && ismember(uvImgRef(2), vRange);
        if(HasInpainedView)
            idImgRef = sum(uRange<uvImgRef(1))*nV + sum(vRange<=uvImgRef(2));
        else
            uvImgRef = uvImgRefDefault;
            idImgRef = idImgRefDefault;
        end
    end
    
end


fMask = dir([DataRoot LFName '/' MaskDir]);
if(HasMask) %Determine name of mask image.
    fMaskName = regexp({fMask.name},['^0*' num2str(uvImgRef(2)) '_0*' num2str(uvImgRef(1)) '[.]png$'],'match');
    fMaskName=[fMaskName{:}];fMaskName=fMaskName{1};
end
if(HasInpainedView) 
    fInpName = fInp{1}; %Get filename name of inpainted view.
    %Search for segmentation image of inpainted view.
    fInpSegName = regexp({fMask.name},['^0*' num2str(uvImgRef(2)) '_0*'  num2str(uvImgRef(1)) '_inp_label[.]png$'],'match');
    fInpSegName = [fInpSegName{:}];
    if(~isempty(fInpSegName))
        fInpSegName=fInpSegName{1};
        HasInpSeg=true;
    else
        HasInpSeg=false;
    end
end

if(length(idImgRef)>1), warning('More than one view is already inpainted'); end

if(~HasMask)
    warning('No Mask found for any view of the light field!');
end
if(~HasInpainedView)
    warning('No Inpainted image found for any view of the light field!');
end

clear fInp fMask exprNumU exprNumV uvImgRefDefault idImgRefDefault
clear M_org M_org_inp Mask
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Load Light Field images
fViews = dir([DataRoot LFName '/']);
idx=1;
for u=uRange
    for v=vRange
        %Load Data
        fName = regexp({fViews.name},['^0*' num2str(v) '_0*' num2str(u) '[.]png$'],'match');
        fName=[fName{:}];
        if(isempty(fName)), error(['Cant find image file of light field view : u=' num2str(u) ', v=' num2str(v) ', in folder ' pwd '/' DataRoot LFName '/']);end
        fName=fName{1};
        TmpImg = double(imread([DataRoot LFName '/' fName]));
        if(MaskFromLabel)
        	TmpMask = imread([DataRoot LFName '/' LabelsDir fName]);
        elseif(HasMask)
            TmpMask = imread([DataRoot LFName '/' MaskDir  fMaskName]);
        else
            TmpMask = false(size(TmpImg));
        end
        
        %Conversion steps
        TmpMask = TmpMask(:,:,1)==objLblColor(1) & TmpMask(:,:,2)==objLblColor(2) & TmpMask(:,:,3)==objLblColor(3);
        
        if(dilateMask>0)
            TmpMask = enlarge_GTmask(TmpMask,dilateMask);
        end
        
        if(luma_only)
            TmpImg = 0.2126 * TmpImg(:,:,1) + 0.7152 * TmpImg(:,:,2) + 0.0722 * TmpImg(:,:,3);
            %TmpImgInp = 0.2126 * TmpImgInp(:,:,1) + 0.7152 * TmpImgInp(:,:,2) + 0.0722 * TmpImgInp(:,:,3);
            nChan = 1;
        else
            nChan = 3;
        end
        
        %Size reduction / cropping
        
        if(~isempty(ROI))
            TmpImg = TmpImg(ROI(2):ROI(4),ROI(1):ROI(3),:);
            %TmpImgInp = TmpImgInp(ROI(2):ROI(4),ROI(1):ROI(3),:);
            TmpMask = TmpMask(ROI(2):ROI(4),ROI(1):ROI(3),:);
        end
        
        if(rescale~=1)
            TmpImg=imresize(TmpImg, rescale, 'nearest');
            %TmpImgInp=imresize(TmpImgInp, rescale);
            TmpMask=imresize(TmpMask, rescale,'nearest');
        end
        
        imgSize = [size(TmpImg,1) size(TmpImg,2)];
        
        % Rearrange Data into Matrices (Tensors for color) suitable for low rank matrix completion.
        if(SingleFormat)
            M_org(:,idx,:) = single(reshape(TmpImg, imgSize(1)*imgSize(2),nChan));
        else
            M_org(:,idx,:) = reshape(TmpImg, imgSize(1)*imgSize(2),nChan);
        end
        %M_org_inp(:,idx,:) = reshape(TmpImgInp, imgSize(1)*imgSize(2),nChan);
        Mask(:,idx) = reshape(TmpMask, imgSize(1)*imgSize(2),1);
        
        idx=idx+1;
    end
end
clear TmpImg TmpMask fViews

%%Load Area Mask if it exists
try
    AreaMask = imread([DataRoot LFName '/' MaskDir 'AreaMask.png']);
    AreaMask = logical(AreaMask(:,:,1));
    if(~isempty(ROI)), AreaMask = AreaMask(ROI(2):ROI(4),ROI(1):ROI(3),:);end
    if(rescale~=1), AreaMask=imresize(AreaMask, rescale,'nearest');end
    AreaMask = reshape(AreaMask, imgSize(1)*imgSize(2),1);
catch err
    AreaMask = vec( enlarge_GTmask(reshape(Mask(:,idImgRef),imgSize), round(30*rescale)) );
end
numPxArea = sum(AreaMask(:));

if(HasInpainedView)
    %%Load Segmentation of inpainted view if it exists
    if(HasInpSeg)
        InpSeg = imread([DataRoot LFName '/' MaskDir fInpSegName]);
        InpSeg = InpSeg(:,:,1);%InpSeg = logical(InpSeg(:,:,1));%
        if(~isempty(ROI)), InpSeg = InpSeg(ROI(2):ROI(4),ROI(1):ROI(3),:);end
        if(rescale~=1), InpSeg=imresize(InpSeg, rescale,'nearest');end
        if(~islogical(InpSeg))
            InpSeg = reshape(InpSeg, imgSize(1)*imgSize(2),1);
            InpSegLabels=find(histc(InpSeg(:),[0:255]))-1;
        end
    else
        %InpSeg = vec( enlarge_GTmask(reshape(Mask(:,idImgRef),imgSize), round(30*rescale)) );
    end
    
    if(SingleFormat)
        M_org_inp = single(imread([DataRoot LFName '/' fInpName]));
    else
        M_org_inp = double(imread([DataRoot LFName '/' fInpName]));
    end
    if(luma_only), M_org_inp = 0.2126 * M_org_inp(:,:,1) + 0.7152 * M_org_inp(:,:,2) + 0.0722 * M_org_inp(:,:,3);end
    if(~isempty(ROI)), M_org_inp = M_org_inp(ROI(2):ROI(4),ROI(1):ROI(3),:);end
    if(rescale~=1), M_org_inp=imresize(M_org_inp, rescale);end
    M_org_inp = reshape(M_org_inp, imgSize(1)*imgSize(2),nChan);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%Define random mask instead :
%% Fully random mask
%{
rng('default');rng(0);
Mask = rand(imgSize(1)*imgSize(2),numImages)>.5;
%}
%% Fixed block on 1 or several views
%{
rng('default');rng(0);
Mask=false(imgSize);
Mask(75:115,95:135)=1;%Mask(2:end-1,2:end-1)=1;%
Mask=repmat(Mask,1,1,numImages);
idImgRef = centralView;%SelectViews = rand(numImages,1)>.5;%
Mask = reshape(Mask,imgSize(1)*imgSize(2),numImages);
M_org_inp=M_org;
%}
