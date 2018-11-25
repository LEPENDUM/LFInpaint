function [HList,ImgSizeList] = loadHomogData(folder)

HList=[];
ImgSizeList=[];

f = dir(folder);
f = regexp({f.name},'.*_k1_p1.*.mat$','match');%only load mat files with parameters k=1, p=1
f = [f{:}];

if(exist([folder '/_defaultImgSize.txt'],'file'))
    eval(fileread([folder '/_defaultImgSize.txt']));
    defaultImgSize = imgSize;
    givenDefault=true;
else
    givenDefault=false;
end

numHAll=0;
for fi=1:length(f)
    clear imgSize
    Vars = whos('-file',[folder '/' f{fi}]);
    Vars = {Vars.name};
    
    if(any(ismember(Vars,'xi')))
        if(any(ismember(Vars,'imgSize')))
            load([folder '/' f{fi}],'imgSize');
        elseif(givenDefault)
            imgSize = defaultImgSize;
        else
            warning(['File "_defaultImgSize.txt" not found and data variable ''imgSize'' not found in file :' f{fi} '. -> File ignored!']);
        end
        load([folder '/' f{fi}],'xi');
        
        numH = length(xi);        
        HList(:,end+1:end+numH) = zeros(8,numH);
        ImgSizeList(:,end+1:end+numH) = zeros(2,numH);
        
        
        err=zeros(1,numH);
        for i=[1:numH]
            HList(:,numHAll+i) = xi{i}(1:8);
            ImgSizeList(1,numHAll+i)=imgSize(1);
            ImgSizeList(2,numHAll+i)=imgSize(2);
            err(i) = sum((xi{i} - [1 0 0 0 1 0 0 0]').^2);
        end
        
        %Remove from the list the homography that is the closest to identity (corresponds to central reference view).
        [~, idMin] = min(err);
        HList(:,numHAll+idMin)=[];
        ImgSizeList(:,numHAll+idMin)=[];
        
        numHAll = size(HList,2);
        
    else
        warning(['Data variable ''xi'' not found in file :' f{fi} '. -> File ignored!']);
    end
end

%{
BouquetFlower1
BouquetFlower2
BouquetFlower3
%}