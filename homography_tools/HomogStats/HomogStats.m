% Computes statistics on homographhy matrix data.
%
% Requires Hd, imgSize.
%
% -Assumes homography data is stored in the matrix Hd
%(8xn matrix with n = number of homographies in the data.)
%
% -WARNING : Assumes that all of the n homographies where found for images
%of the same size given in imgSize (imgSize(1)=vertical definition /
%imgSize(1)=horizontal definition).


%Homography matrix elements chosen for display of bi-variate distrib :
Var1 = [1,1];% element 1
Var2 = [3,1];% element 2
plotBiVarFit = false;
plotBiVarData = false;
plot8Independant = true;
saveMGDStats = false;
saveMGGDStats = false;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%           Load Hd and image sizes from files                %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if(exist('Hd','var') && exist('imgSize','var') )
    numH = size(Hd,2);
    imgSizes = repmat(imgSize,1,size(Hd,2));
    numDataBases=1;databaseHRange=[1 numH];
end

%%{
curDir = fileparts(mfilename('fullpath'));
databaseFolders = {...
    'F0_homography_k1',...
    %'Illum_homography_k1',...
    };
numDataBases = length(databaseFolders);
HdCell = cell(numDataBases,1);
ImgSizesCell = cell(numDataBases,1);
Hd=[];
ImgSizes=[];
databaseHRange=ones(numDataBases+1,1);
for i=1:numDataBases
    [HdCell{i},ImgSizesCell{i}] = loadHomogData([curDir '/' databaseFolders{i}]);
    Hd = cat(2,Hd,HdCell{i});
    ImgSizes = cat(2,ImgSizes,ImgSizesCell{i});
    databaseHRange(i+1) = size(Hd,2);
end

numH = size(Hd,2);
%}
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


meanH = eye(3,3);
maxH=zeros(3,3);
minH=zeros(3,3);
bH=zeros(3,3);

%Convert homography matrix data to unit size image homography.
HdUnit = Hd;
for i=1:numH
    szX = ImgSizes(2,i);
    szY = ImgSizes(1,i);
    HdUnit(2,i) = HdUnit(2,i) * szY / szX;
    HdUnit(4,i) = HdUnit(4,i) * szX / szY;
    HdUnit(3,i) = HdUnit(3,i) / szX;
    HdUnit(6,i) = HdUnit(6,i) / szY;
    HdUnit(7,i) = HdUnit(7,i) * szX;
    HdUnit(8,i) = HdUnit(8,i) * szY;
end
%HdUnit=Hd;

HdCent = bsxfun(@minus,HdUnit,mean(HdUnit,2));%centered homography data

if(plot8Independant),figure,set(gcf,'color','white');end

for j=1:3
for i=1:3
    
    idx = i + (j-1)*3;
    idxT = j + (i-1)*3;

    if(idx==9)
        h_legend=legend('histogram',['generalized' 10 'gaussian fit'],'gaussian fit');
        set(h_legend, 'FontSize', 14.5, 'FontName', 'Times New Roman');
        continue;
    end

    hvec = HdUnit(idx,:);

    meanVal = mean(hvec);
    minVal=min(hvec);
    maxVal=max(hvec);

    maxH(i,j) = maxVal;
    minH(i,j) = minVal;
    meanH(i,j) = meanVal;
    
    %Plot univariate distribution with gaussian fitting for each element of the
    %homography matrix :
    if(plot8Independant)
        %Actual data histograms
        bins = [minVal:(maxVal-minVal)/200:maxVal];
        histo=histc(hvec,bins);
        histo = histo ./ sum(histo);
        
        %Gaussian fit
        sig=sqrt(var(hvec));
        gaussFit = normpdf(bins,meanVal,sig);
        gaussFit = gaussFit/sum(gaussFit);
        
        %Generalized gaussian fit
        [V, b, m, ~] = RA_FP(HdCent(idx,:),0);
        V=V*m;
        bH(i,j)=b;
        GDFit = exp(-.5*((bins-meanVal).^2/V).^b);
        GDFit = GDFit ./ sum(GDFit);
        
        h_plot=subplot(3,3,idxT);
        plot(bins,histo,'linewidth',1.4); hold on;
        plot(bins,GDFit,':r','linewidth',1.5);
        plot(bins,gaussFit,'-g','linewidth',1.1);
        
        xlabel(['h_{' num2str(i) ',' num2str(j) '}'],'FontWeight','bold','FontSize',14,'FontName','Times New Roman');
        xlabh = get(gca,'XLabel');
        set(xlabh,'Units','Normalized','Position',[.84 .98 -1]);
        set(h_plot,'FontName','Times New Roman');
        axis([meanVal-3.2*sig, meanVal+3.2*sig, 0, max(GDFit)*1.2]);
    end

end
end

%Plot bivariate gaussian distribution for the two specified variables in Var1, Var2.
Var1idx = Var1(1) + (Var1(2)-1)*3;
Var2idx = Var2(1) + (Var2(2)-1)*3;

if(plotBiVarFit)
    mu = [meanH(Var1(1),Var1(2)) meanH(Var2(1),Var2(2))];
    
    x1 = minH(Var1(1),Var1(2)):(maxH(Var1(1),Var1(2))-minH(Var1(1),Var1(2)))/100:maxH(Var1(1),Var1(2));
    x2 = minH(Var2(1),Var2(2)):(maxH(Var2(1),Var2(2))-minH(Var2(1),Var2(2)))/100:maxH(Var2(1),Var2(2));
    [X1,X2] = meshgrid(x1,x2);
    Xc = [X1(:)-mu(1) X2(:)-mu(2)];
    
    %bivariate Gaussian fit
    Cov = cov(HdCent([Var1idx Var2idx],:)');
    Gauss2 = mvnpdf([X1(:) X2(:)],mu,Cov);
    Gauss2 = reshape(Gauss2,length(x2),length(x1));
    figure,surf(x1,x2,Gauss2,'LineStyle', 'none');
    
    %bivariate generalized gaussian
    [V, b, m, ~] = RA_FP(HdCent([Var1idx Var2idx],:),0);
    V=V*m;
    const = (gamma(1)/(pi*gamma(1/b)*2^(1/b)))*b/sqrt(det(V));
    Gengauss2 = const*exp(-0.5*diag(Xc*inv(V)*Xc').^b);
    Gengauss2 = reshape(Gengauss2,length(x2),length(x1));
    figure,surf(x1,x2,Gengauss2,'LineStyle', 'none');
    
end

%Plot the point cloud for the two specified variables in Var1, Var2.
if(plotBiVarData)
    figure,
    for i=1:numDataBases
        plot(HdUnit(Var1idx,databaseHRange(i):databaseHRange(i+1)),HdUnit(Var2idx,databaseHRange(i):databaseHRange(i+1)),'o');hold on;
    end
end

curDir = fileparts(mfilename('fullpath'));
%8-var Gaussian distrib : Covariance matrix and means for the 8 homography variables.
Cov = cov(HdCent(1:8,:)');
Means = meanH(1:8);
if(saveMGDStats)
	save([curDir '/Distrib/MGD_Homog.mat'],'Cov','Means');
end
Cov0=Cov;

%8-var GGD : Covariance matrix, means for the 8 homography variables + beta parameters.
[Cov, b, m, ~] = RA_FP(HdCent(1:8,:),0);
Cov=Cov*m;
if(saveMGGDStats)
    save([curDir '/Distrib/MGGD_Homog.mat'],'Cov','Means','b');
end