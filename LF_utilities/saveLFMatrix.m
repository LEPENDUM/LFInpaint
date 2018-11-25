% Save Ligh Field views in a folder
%Input :
% - LFMatrix : Matrix of views (each column is a vectorized view)
% - imgSize  : (y,x) spatial resolution of the views.
% - uRange   : Range of values for angluar dimension u.
% - vRange   : Range of values for angluar dimension v.
% - filename : filename to save to.
% - video    : 0 -> Save png image for each view ('filename_u_v.png').
%              1 (default) -> Save as a video file with specified codec;
% - ScanType : (only if video>0) : defines scanning order for the video.
%              - 0 = Zscan.
%              - 1 = ZigZag scan.
%              - 2 = Scan vertically, horizontally and diagonally along the lightfield (does not show all the views).
% - codec : (string) Codec to use for video encoding :
%                       - 'lossless' (=> use Uncompressed AVI for color data and Grayscale AVI otherwise).
%                       - 'MPEG-4' (default)    | MPEG-4 file with H.264 encoding
%                       - 'Motion JPEG AVI'     | AVI file using Motion JPEG encoding
%                       - 'Motion JPEG 2000'    | Motion JPEG 2000 file
%                       - 'Archival'            | Motion JPEG 2000 file with lossless compression
%                       - 'Indexed AVI'         | Uncompressed AVI file with indexed video
% - FrameRate : frame rate used if the LF is saved as a video.

function saveLFMatrix (LFMatrix, imgSize, uRange, vRange, filename, video, ScanType, Codec, FrameRate )

nU = length(uRange);
nV = length(vRange);

numImages=nU*nV;
nChan = size(LFMatrix,3);
if(numImages ~= size(LFMatrix,2)),
    error('Wrong matrix size : number of columns should be equal to number of views.')
end

%% Save in image format
if(exist('video','var') && ~video)
    idx=1;
    for u=uRange
        for v=vRange
            imwrite(reshape(uint8(LFMatrix(:,idx,:)),[imgSize,nChan]), [filename '_' num2str(v) '_' num2str(u) '.png']);
%            imwrite(reshape(uint8(LFMatrix(:,idx,:)),[imgSize,nChan]), [filename '_' num2str(idx,'%.3d') '.png']);
            idx=idx+1;
        end
    end
else
%% Save in video format
    if(~exist('Codec','var'))
        Codec = 'MPEG-4';
    end
    
    %% Define path along the views    
    path=1:numImages;
    if(exist('ScanType','var') && ScanType==1)
        for idx=2:2:nU
            path(1+nV*(idx-1):nV*idx)=nV*idx:-1:1+nV*(idx-1);
        end
    elseif(ScanType==2)
        path=[[1:nV] ...
         [nV:nV:numImages]...
         [numImages:-1:numImages-nV+1]...
         [numImages-nV+1:-nV:1]...
         [1:nV+1:numImages]];
    end

    frame_id=0;
    
    VidSize = ceil(imgSize/2)*2;
    
    Video = zeros(VidSize(1),VidSize(2), nChan, length(path));
    for idx=path
        frame_id = frame_id+1;
        Video(1:imgSize(1),1:imgSize(2),:,frame_id) = reshape(LFMatrix(:,idx,:),[imgSize,nChan]);
    end
    
    %% Save in uncompressed video format
    if(strcmp(Codec,'lossless'))
        if(size(LFMatrix,3)==1), format = 'Grayscale AVI';
        else                     format = 'Uncompressed AVI';end
    else
        format = Codec;
    end
    
    v = VideoWriter(filename,format);
    if(exist('FrameRate','var')), v.FrameRate=FrameRate; end
    open(v);
    for frame_id=1:length(path)
        writeVideo(v,uint8(Video(:,:,:,frame_id)));
    end
    close(v);
end
