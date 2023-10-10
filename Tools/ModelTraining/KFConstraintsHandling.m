% [x,a]=KFConstraintsHandling(1,1,[1;1;1],[0;-3],1)
function [xTrunc,PTrunc,Success]=KFConstraintsHandling(xTrunc,PTrunc,D,d,param)
xTrunc=gather(xTrunc);
PTrunc=gather(PTrunc);
xTruncOut=zeros(size(xTrunc));
PTruncOut=zeros(size(PTrunc));
Success=1;
try
parfor ki=1:size(PTrunc,3)
    [xTruncOut(:,ki),PTruncOut(:,:,ki)]=KFConstraintsHandlingSeq(xTrunc(:,ki),PTrunc(:,:,ki),D,d,1);
end
catch
    sprintf('Imposing constraints have failed');
end
xTrunc=xTruncOut;
PTrunc=PTruncOut;
end