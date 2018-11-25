% Create a set of homography warpings of the input image with different 
% possible parmetrizations of the homographies. The resulting warped images 
% are vectorized as column vectors and arranged in the output matrix D.
%
% Input :
%   - img : the input image to warp;
%   - numCols : the number of homography warpings to process (=number of
%   column of the output matrix D).
%   - HomogCenter : string definig how the center C is determined for
%   'RectDisplace' and '3DRot' options (parameter ignored for '8rand' method):
%                   - 'random' : randomly chosen center for each homography.
%                   - 'grid'   : several possible centers defined on a regular grid.
%                   - 'fixTL'  : always top-left point of the image.
%                   - 'fixTC'  : always top-center point of the image.
%                   - 'fixTR'  : always top-right point of the image.
%                   - 'fixML'  : always middle-left point of the image.
%                   - 'fixMC'  : always middle-center point of the image.
%                   - 'fixMR'  : always middle-right point of the image.
%                   - 'fixBL'  : always bottom-left point of the image.
%                   - 'fixBC'  : always bottom-center point of the image.
%                   - 'fixBR'  : always bottom-right point of the image.
%   - mask : (optional), binary image to warp along with the input image. Can be set to empty to ignore.
%   - ParamMethod : parametrization method :
%                   - 'RectDisplace' : random parametrization based on the
%                     displacement of the 4 corners of a rectangle centered on
%                     a given point C (see HomogCenter parameter).
%                   - '3DRot' : Homorgaphies parameterized with 2 angles corresponding to a
%                     3D rotation (centered on a given point C) of the image plane and re-projection.
%                   - '8rand' : Generate random homography matrix from precomputed multivariate 
%                     distributions (gaussian or generalized gaussian).
%
% Output :
%   - D : Matrix containing all the warped images vectorized in columns.
%   - DMask : Matrix containing the warped masks vectorized in columns.
%   - H : Matrix of homography parameters (8 parameters per columns).
%   - InWarp: Matrix containing the masks of unknown values at the border of the image resulting from
%   the homography warping. 0 = unknown / 1 = known.
%
function [D,DMask,H,InWarp] = CreateHomogDico(img,numCols,HomogCenter,mask,ParamMethod)

if(~exist('ParamMethod','var') || strcmp(ParamMethod,'RectDisplace'))
    ParamMethod=0;
elseif(strcmp(ParamMethod,'3DRot'))
    ParamMethod=1;
elseif(strcmp(ParamMethod,'8rand'))
    ParamMethod=2;
else
    warning('Unknwon param method : Using ''RectDisplace'' by default.');
    ParamMethod=0;
end

[imgSize(1), imgSize(2), nChan] = size(img);

SingleHomogCenter = -1;
if(ParamMethod~=2) %SingleHomogCenter not valid for method '8rand'.
    switch HomogCenter
        case 'random',
            SingleHomogCenter = 1;
        case 'grid',
            SingleHomogCenter = 0;
        case 'fixTL',
            xc = 1;
            yc = 1;
        case 'fixTC',
            xc = round(imgSize(2)/2);
            yc = 1;
        case 'fixTR',
            xc = imgSize(2);
            yc = 1;
        case 'fixML',
            xc = 1;
            yc = round(imgSize(1)/2);
        case 'fixMC',
            xc = round(imgSize(2)/2);
            yc = round(imgSize(1)/2);
        case 'fixMR',
            xc = imgSize(2);
            yc = round(imgSize(1)/2);
        case 'fixBL',
            xc = 1;
            yc = imgSize(1);
        case 'fixBC',
            xc = round(imgSize(2)/2);
            yc = imgSize(1);
        case 'fixBR',
            xc = imgSize(2);
            yc = imgSize(1);
        otherwise,
            warning( [HomogCenter ' : is not a valid string for that parameter. Using ''random'' by default!']);
            SingleHomogCenter = 1;
    end
end

%if a mask is given (with same dimensions as the image), it will be warped with the same homographies.
MoveMask = exist('mask','var') && size(mask,1)==imgSize(1) && size(mask,2)==imgSize(2);

ExtOffsetX = imgSize(2)/2;
ExtOffsetY = imgSize(1)/2;

D = zeros(imgSize(1)*imgSize(2),numCols,nChan,class(img));
H = zeros(9,numCols,class(img));

if(MoveMask && ~islogical(mask))
    DMask = zeros(imgSize(1)*imgSize(2),numCols);
else
    DMask = false(imgSize(1)*imgSize(2),numCols);
end

alphaBounds = [0,pi];
betaBounds  = [0, pi/10];
%numAlpha = 10;
%numBeta = ceil(numCols/numAlpha);
%betaList = linspace(pi/48, pi/12, numBeta);


if(SingleHomogCenter)
    numC = numCols;
    numParams = 1;
else
    distHx = imgSize(2)*2/3;
    distHy = imgSize(1)*2/3;
    %distHx=100; distHy=100;
    
    xcRange = -ExtOffsetX : distHx : imgSize(2)+ExtOffsetX ;
    ycRange = -ExtOffsetY : distHy : imgSize(1)+ExtOffsetY ;
    numC = length(xcRange)*length(ycRange);
    numParams = ceil(numCols/numC);
end

i=1;
for C_id = 1:numC
    
    if(SingleHomogCenter==1)
        xc = -ExtOffsetX + rand * (imgSize(2)+2*ExtOffsetX);
        yc = -ExtOffsetY + rand * (imgSize(1)+2*ExtOffsetY);
    elseif(SingleHomogCenter==0)
        xc = xcRange(floor((C_id-1)/length(ycRange))+1);
        yc = ycRange(mod(C_id-1,length(ycRange))+1);
    end
    
    if(ParamMethod==0) %'RectDisplace' parametrization
        distAnchorPointx =  imgSize(2)/2;
        distAnchorPointy =  imgSize(1)/2;
        %distAnchorPointx=100; distAnchorPointy=100;
        
        x0 = xc-distAnchorPointx/2;
        y0 = yc-distAnchorPointy/2;
        startpoints = [x0 y0; x0+distAnchorPointx y0; x0+distAnchorPointx y0+distAnchorPointy; x0 y0+distAnchorPointy];
    end
    
    for Param_idx = 1:numParams

        if(ParamMethod==1) %'3DRot' parametrization.
            %idA = mod(idx-1,numAlpha)+1;
            %idB = floor((idx-1)/numAlpha)+1;
            %tform.T = inv(genHomography3DRot((idA-1)*pi/numAlpha,betaList(idB),[xc,yc]))';

            alpha = alphaBounds(1) + rand * (alphaBounds(2)-alphaBounds(1));
            beta = betaBounds(1) + rand * (betaBounds(2)-betaBounds(1));
            tform.T = inv(genHomography3DRot(alpha,beta,[xc,yc]))';

        elseif(ParamMethod==0) %'RectDisplace' parametrization.
            endpoints(:,1) = startpoints(:,1) + (rand(4,1)-.5) * distAnchorPointx/100;
            endpoints(:,2) = startpoints(:,2) + (rand(4,1)-.5) * distAnchorPointy/100;
            tform = fitgeotrans(startpoints,endpoints,'projective');

        elseif(ParamMethod==2) %'8rand' generation method.
            
            curDir=fileparts(mfilename('fullpath'));
            load([curDir '/Distrib/MGGD_Homog_F0.mat']);
            
            %Generate random homography matrix from multivariate
            %distributions (gaussian or generalized gaussian).
            %Homography distributions must be calibrated for images of size 1x1.
            if(exist('b','var'))
                %multivariate generalized gaussian distrib with shape parameter b.
                tform.T = MGGD_generation(1, 8, Cov, b);
                tform.T = tform.T + Means';
            else
                %multivariate gaussian distrib.
                tform.T = mvnrnd(Means,Cov);
            end
            
            tform.T(9) = 1;
            tform.T = reshape(tform.T,3,3);
            
            %Extend the generated homography matrix for the image size considered.
            szX = imgSize(2);
            szY = imgSize(1);
            tform.T(2,1) = tform.T(2,1) * szX / szY;
            tform.T(1,2) = tform.T(1,2) * szY / szX;
            tform.T(3,1) = tform.T(3,1) * szX;
            tform.T(3,2) = tform.T(3,2) * szY;
            tform.T(1,3) = tform.T(1,3) / szX;
            tform.T(2,3) = tform.T(2,3) / szY;

        end
        
        [imgWarp,warpMask] = homographyWarp(tform.T, img);
        D(:,i,:)   = reshape( imgWarp, [imgSize(1)*imgSize(2),nChan] );
        InWarp(:,i) = vec(warpMask);
        H(:,i) = vec(tform.T);
        
        if(MoveMask)
            if(islogical(mask))
                DMask(:,i) = ~InWarp(:,i) | vec( round(homographyWarp(tform.T, double(mask))) );
            else
                DMask(:,i) = vec( homographyWarp(tform.T, mask, 'nearest') );
                DMask(~InWarp(:,i),i)=nan;
            end
        else
            DMask(:,i) = ~InWarp(:,i);
        end
        i=i+1;
    end

end
