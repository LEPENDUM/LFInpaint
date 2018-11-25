SimuRoot = fileparts(mfilename('fullpath'));%'LowRankInpaint/';
InpResDir = [SimuRoot '/InpaintRes/' LFName '/'];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
rng('default');rng(0);


MoveAreaMask = false;	% Wrap the AreaMask for each warped view.
UseAllMasks  = false;	% Use the masks of all the views (if false, use only the mask of the inpainted view).
UseInpSeg    = 2;		% Use binary segmentation of inpainted view if it is given: 1=Warp only 1 segment and leave the other part empty (to fill with LRMC) / 2=Warp both segments (with different homographies).
UnknownExt   = floor(10*rescale); % Number of pixels to enlarge the area of unknown data for the Matrix Completion.

YUVConvert = 1; %0-false=none (RGB) / 1=Y'u'v' / 2=YCbCr / 3=HSV

skipDictGen = false;


%TgtViews = [1];
%TgtViews = setdiff(1:numImages,idImgRef);
TgtViews = 1:numImages;
numPxArea = sum(AreaMask(:));

                     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%    Generate Homographies    %%%%%%%%%%%%%%%%%%%%%%%%
                     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if(~skipDictGen)

    Dsize = 400 ;
    HomogCenter = 'random';
    ParamMethod = '8rand';%'RectDisplace' / '3DRot' / '8rand'
    
    ImgRef = M_org_inp;
    
    if(YUVConvert && nChan==3)
        if(YUVConvert==1)%Yu'v'
            ImgRef = squeeze( XYZ2Luv(sRGB_lin2XYZ(permute((ImgRef/255).^2.2,[1 3 2]))) );
            ImgRef(:,1) = ImgRef(:,1).^(1/2.2);
        elseif(YUVConvert==2) %YCbCr
            ImgRef = squeeze(RGB2YCbCr(permute(ImgRef,[1 3 2]),256));
        else%HSV
            ImgRef = squeeze(rgb2hsv(permute(ImgRef,[1 3 2])));
        end
    end
    
    if(UseInpSeg && HasInpSeg)
        if(UseInpSeg==1)
            Dsize1=200; Dsize2=200; Dsize = Dsize1+Dsize2;
        elseif(UseInpSeg==2)
            Dsize1=Dsize; Dsize2=Dsize;
        else
            error('UseInpSeg should be either 0 (=false),1 or 2.');
        end
        
        if(islogical(InpSeg))
            [D, DMask, Hd, InWarp] = CreateHomogDico(reshape(ImgRef,[imgSize nChan]),Dsize1,HomogCenter,reshape(InpSeg,imgSize),ParamMethod);
        	[D2, DMask2, Hd2, InWarp2] = CreateHomogDico(reshape(ImgRef,[imgSize nChan]),Dsize2,HomogCenter,~reshape(InpSeg,imgSize),ParamMethod);
        else
            [D0, DMask0, Hd, InWarp] = CreateHomogDico(reshape(ImgRef,[imgSize nChan]),Dsize1,HomogCenter,reshape(InpSeg,imgSize),ParamMethod);
        end

        if(UseInpSeg==1)
            if(islogical(InpSeg))
                D(:,Dsize1+1:Dsize1+Dsize2,:) = D2;
                DMask(:,Dsize1+1:Dsize1+Dsize2) = DMask2;

                D(:,end+1,:) = ImgRef;
                DMask(:,end+1) = 0;
                Dsize = Dsize+1;
            else
                error('Not implemented! Segmentation map of the inpainted view must be binary to use this mode.')
            end
        
        elseif(UseInpSeg==2)
            if(islogical(InpSeg))
                D=D.*repmat(DMask,1,1,3);
                D(repmat(InWarp2&DMask2,1,1,3))=D2(repmat(InWarp2&DMask2,1,1,3));
                DMask = ~( DMask.*InWarp|DMask2.*InWarp2 );
            else
                D=zeros(size(D0));
                DMask = true(size(DMask0));
                for i=1:Dsize
                    %hi = reshape(Hd(:,i),[3,3]);
                    %lbls = vec(homographyWarp(hi, reshape(InpSeg,imgSize),'nearest'));
                    for lbl_id=length(InpSegLabels):-1:1%randperm(length(InpSegLabels))%
                        Warp_id = mod(i+lbl_id-2,Dsize)+1;
                        lbl = DMask0(:,Warp_id) == InpSegLabels(lbl_id);
                        D(lbl,i,:) = D0(lbl,Warp_id,:);
                        DMask(lbl,i,:) = false;
                    end
                end
            end
        end
        
    else
        if(MoveAreaMask)
            MaskRef = reshape(~AreaMask,imgSize);
        else
            MaskRef = 0;
        end
        [D, DMask, Hd] = CreateHomogDico(reshape(ImgRef,[imgSize nChan]),Dsize,HomogCenter,MaskRef,ParamMethod);
    end
    clear MaskRef ImgRef InWarp



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Add Views to inpaint (+ mask) to the matrix
	if(YUVConvert && nChan==3)
		if(YUVConvert==1)%Yu'v'
			D(:,Dsize+1:Dsize+length(TgtViews),:) = XYZ2Luv(sRGB_lin2XYZ((M_org(:,TgtViews,:)/255).^2.2));
			D(:,Dsize+1:Dsize+length(TgtViews),1) = D(:,Dsize+1:Dsize+length(TgtViews),1).^(1/2.2);
		elseif(YUVConvert==2)%YCbCr
			D(:,Dsize+1:Dsize+length(TgtViews),:) = RGB2YCbCr(M_org(:,TgtViews,:),256);
		else%HSV
			D(:,Dsize+1:Dsize+length(TgtViews),:) = rgb2hsv(M_org(:,TgtViews,:));
		end
	else
		D(:,Dsize+1:Dsize+length(TgtViews),:) = M_org(:,TgtViews,:);
	end
	if(UseAllMasks)
		DMask(:,Dsize+1:Dsize+length(TgtViews)) = Mask(:,TgtViews);
	else
		DMask(:,Dsize+1:Dsize+length(TgtViews)) = repmat(Mask(:,idImgRef),1,length(TgtViews));
	end
	% Enlarge the Masks before matrix completion
	for idx=Dsize+1:size(D,2)
		DMask(:,idx) = vec(enlarge_GTmask(reshape(DMask(:,idx),imgSize),UnknownExt));
	end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

else
    warning('''skipDictGen'' is active : The dictionnary must be correctly initialized');
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear X_rec0 X_rec
%% Solve low rank matrix completion

%Low rank completion parameters
alpha = [1, 0, 0];
betaMultLum = 1.44;
betaMultColor = 5;
maxIter = 100;
epsilonLum = [1e-2, 3e-3];%[1e-2, 3e-3];%
epsilonColor = [3e-2, 3e-3];%epsilon = [2e-2, 1e-2];%epsilon = [1e-1, 3e-1];
TraceNormToRankParam = inf; %0->minimize trace norm / inf->minimize rank

tic

clear X
X_rec0 = zeros(numPxArea,size(D,2),nChan);
X = X_rec0;
startId=1;
for chan=1:nChan
    if(chan>1 && YUVConvert>0)
        betaMult = betaMultColor;
        epsilon = epsilonColor;
        startId=round(Dsize*3/4);%Use less homographies for chroma components (faster)
		NoisyCols = [Dsize+1:size(D,2)]+1-startId;
    else
        betaMult = betaMultLum;
        epsilon = epsilonLum;
		NoisyCols = [Dsize+1:size(D,2)];
    end
    [X_rec0(:,startId:end,chan), X(:,startId:end,chan), ~, rank, errList] = LRTCADMM(...
    single(D(AreaMask,startId:end,chan).*~DMask(AreaMask,startId:end)), ... % a tensor whose elements in Omega are used for estimating missing value
    ~DMask(AreaMask,startId:end,:),...  % the index set indicating the observed elements
    single(alpha),...             	    % the coefficient of the objective function
    single(betaMult),...                % Parameter (>1) to control the algorithm speed (close to 1 : slow but accurate / high values gives faster but more blurry results).
    maxIter,...                         % the maximum iteration
    single(epsilon),...            	    % the tolerance
    TraceNormToRankParam,...            % parameter (set to 0 for a trace norm minimization / higher to preserve the high frequencies of the matrix)
    uint32(NoisyCols));
end

toc
        
if(YUVConvert && nChan==3)
    if(YUVConvert==1)%Yu'v'
        X_rec0(:,:,1) = realpow(max(X_rec0(:,:,1),0), 2.2);
        X_rec0 = 255 * realpow(max(XYZ2sRGB_lin(Luv2XYZ(X_rec0)),0), 1/2.2);
    elseif(YUVConvert==2)%YCbCr
        X_rec0 = YCbCr2RGB(X_rec0,256);
    else%HSV
        X_rec0 = hsv2rgb(X_rec0);
    end
end
        



X_rec = M_org;
X_rec(AreaMask,TgtViews,:) = X_rec0(:,Dsize+1:Dsize+length(TgtViews),:);
%X_test = zeros(imgSize(1)*imgSize(2),Dsize,nChan);
%X_test(AreaMask,:,:)=X_rec0(:,1:Dsize,:);

if(~any(idImgRef==TgtViews))
    X_rec(:,idImgRef,:) = M_org_inp;
end

Video = displayLFMatrix(X_rec,imgSize,nU,nV,1);

%Video = displayLFMatrix(D.*repmat(AreaMask,1,size(D,2))*255, imgSize, nU, nV, 1);
%saveLFMatrix(uint8(X_rec),imgSize,uRange,vRange,[InpResDir 'SubFolder/inp'],1,2);
%R=refocus(X_rec/255,imgSize,uRange,vRange,uvImgRef,[-1:.1:1],baselineUV);
%saveLFMatrix(permute(reshape(255*R,[imgSize(1)*imgSize(2),nChan,size(R,4)]),[1,3,2]),imgSize,1,1:size(R,4),[InpResDir 'SubFolder/fstack/fstack'],0,0);