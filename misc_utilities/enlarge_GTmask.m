%% enlarge the ground truth binary mask zone (segmentation mask)
function Im_enlarged = enlarge_GTmask(Im, iter)
%Im : ground truth binary mask
%iter : number of iterations. Each iteration enlarges the mask by 1 pixel

%{
[m, n] = size(Im);
Im_enlarged = Im;
for i = 1:iter
    for j = 1:m
        for k = 1:n
            if Im(j,k)==1
                % enlarge the mask zone
                if j-1>0
                   Im_enlarged(j-1,k)=1;
                end
                if j+1<=m
                   Im_enlarged(j+1,k)=1; 
                end
                if k-1>0
                   Im_enlarged(j,k-1)=1;
                end 
                if k+1<=n
                   Im_enlarged(j,k+1)=1; 
                end
            end
        end
    end
    Im = Im_enlarged;
end
%}

%Or alternatively (a bit slower with R2013b / much faster with R2015b):
Im_enlarged = imdilate(Im,strel('diamond',iter));

end

