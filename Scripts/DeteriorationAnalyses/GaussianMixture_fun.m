function [Exsmooth,Vsmooth,y_Data,Re]=GaussianMixture_fun(LambdaModel,LambdaObs,ExElm,VarElm,YObs,Rvals)

ExsmoothMat=nan(size(ExElm));
VsmoothMat=nan(size(VarElm));
y_DataMat=nan(size(YObs));
ReMat=nan(size(Rvals));
for i=1:size(LambdaModel,1)
    ExsmoothMat(:,:,i)=LambdaModel(i,:).*ExElm(:,:,i);
    VsmoothMat(:,:,:,i)=reshape(LambdaModel(i,:),1,1,length(LambdaModel(i,:))).*VarElm(:,:,:,i);
    y_DataMat(i,:)=LambdaObs(i,:).*YObs(i,:);
    ReMat(i,:)=LambdaObs(i,:).*Rvals(i,:);
end
Exsmooth=nansum(cat(3,ExsmoothMat),3);
Vsmooth=nansum(cat(4,VsmoothMat),4);
if size(y_DataMat,1)>1
    y_Data=nansum(y_DataMat);
    Re=nansum(ReMat);
else
    y_Data=y_DataMat;
    Re=ReMat;
end
VarMixMat=nan(size(VarElm));
VarMix_RMat=nan(size(Rvals));

for i=1:size(LambdaModel,1)
    for j=1:size(LambdaModel,2)
        VarMixMat(:,:,j,i)=LambdaModel(i,j).*(ExElm(:,j,i)-Exsmooth(:,j))*(ExElm(:,j,i)-Exsmooth(:,j))';
    end
    IndVal=find(~isnan(YObs(i,:))&(YObs(i,:)~=0));
    VarMix_RMat(i,IndVal)=LambdaObs(i,IndVal).*(YObs(i,IndVal)-y_Data(IndVal)).^2;
end
VarMix=nansum(cat(4,VarMixMat),4);
VarMix_R=nansum(VarMix_RMat);
Vsmooth=Vsmooth+VarMix;
Re=Re+VarMix_R;



end

% Exsmooth=cellfun(@times, ExElm, num2cell(Lambda), 'UniformOutput', false);
% Exsmooth=sum(cat(3,Exsmooth{:}),3);
% Vsmooth=cellfun(@times, VarElm, num2cell(Lambda), 'UniformOutput', false);
% Vsmooth=sum(cat(4,Vsmooth{:}),4);
% y_Data=cellfun(@times, YObs, num2cell(Lambda), 'UniformOutput', false);
% y_Data=sum(cat(3,y_Data{:}),3);
% Re=cellfun(@times, Rvals, num2cell(Lambda), 'UniformOutput', false);
% Re=sum(cat(3,Re{:}),3);