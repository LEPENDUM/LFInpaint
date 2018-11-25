% Refocus Light field data.
% Inputs :
%  - LFMatrix :  Matrix of views (each column is a vectorized view).
% The views are assumed to be in the order :
%     (u0,v0), ... (u0,vn), ..., (un,v0), ... (un,vn)
%
%  - imgSize : vector of size 2 : (y,x) spatial resolution of the views.
%  - uRange  : Range of values for angluar dimension u (low values = top views).
%  - vRange  : Range of values for angluar dimension v (low values = left views).
%  - uvRef    : vector of size 2 : (u,v) indices of the reference view (e.g. central view).
%  - dRange  : angle value for refocusing (if a vector is given, the
%  function generates a video containing the refocused images with all the
%  angles).
%  - baselineUV : - (Optional : default=[1 1]).
%                 - vector of size 2 : factor between u (resp. v) index in
%  uRange (resp. vRange) and physical u (resp. v).
%                 - For example, it should be set to [-1 1] if the order is
%                 inverted in the u dimension (i.e. if high values of uRange correspond to left views).
%  - invertUV : Invert (Optional : default=false). Must be set to true if
%  the indexing of u and v is inverted in the input LFMatrix.
%  - DisableDisplay : (Optional : default=false), does not display the video if set to true (only return the array).

function [I] = refocus(LFMatrix,imgSize,uRange,vRange,uvRef,dRange,baselineUV,invertUV,DisableDisplay,gamma)

if(~exist('invertUV','var') || isempty(invertUV))
    invertUV=false;
end

if(~exist('baselineUV','var') || isempty(baselineUV))
    baselineUV=[1 1];
end
if(~exist('DisableDisplay','var') || isempty(DisableDisplay))
    DisableDisplay=false;
end
if(~exist('gamma','var') || isempty(gamma))
    gamma=1;
end

nChan = size(LFMatrix,3);
stackSize = length(dRange);
I = zeros([imgSize nChan stackSize]);
White = ones(imgSize);

for id_d=1:stackSize
    d=dRange(id_d);
    Count = zeros(imgSize);
    idView=1;
    
    for u=uRange
        for v=vRange
            View  = reshape(LFMatrix(:,idView,:),[imgSize nChan]);
            
            if(invertUV)
                tr_vec = [(v-uvRef(2))*d*baselineUV(2), (u-uvRef(1))*d*baselineUV(1)];
            else
                tr_vec = [(u-uvRef(1))*d*baselineUV(1), (v-uvRef(2))*d*baselineUV(2)];
            end
            %attenuation : max(0, 1-((x.^2)/(r.^2-.1)).^7) :
%            apertureFact = max(0, 1 - (((u-uvRef(1))^2 + (v-uvRef(2))^2 ) / (.5*min(length(uRange),length(vRange))).^2 - .1).^3);
            %attenuation : exp(-x^2n) :
%            apertureFact = exp(-2*(((u-uvRef(1))^2 + (v-uvRef(2))^2) / (.5*min(length(uRange),length(vRange))).^2) .^5 );
            apertureFact=1;
            
            I(:,:,:,id_d) = I(:,:,:,id_d) + apertureFact*imtranslate(View,tr_vec,'cubic').^gamma;
            Count = Count + apertureFact*imtranslate(White,tr_vec,'linear');
            idView = idView+1;
        end
    end
    
    for ch=1:nChan
        I(:,:,ch,id_d) = (I(:,:,ch,id_d)./Count).^(1/gamma);
    end
    

end

if(~DisableDisplay)
    if(stackSize==1)
        figure,imshow(I);
    else
        implay(I);
    end
end