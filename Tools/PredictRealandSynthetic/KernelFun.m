function K = KernelFun(x,y,Param,KernelType)
N = size(x,1);
M = size(y,1);

switch KernelType
    case 'RBF' % RBF kernel including (ARD: Automatic Relevance Determination)
        K=zeros(N,M);
        for i=1:length(Param)
            K=K+CalculateDistance(x(:,i)/Param(i),y(:,i)/Param(i));
        end
        K=exp(-0.5*(K));
    case 'AAK' % RBF kernel including (ARD: Automatic Relevance Determination)
        % Aitchison and Aitken (1976) unordered discrete kernel
        h=Param;
        num_levels = length(unique(y(:,1)));
        K = ones(N,M).* h ./(num_levels - 1);
        for i=1:length(y)
            Ind=(y(i)==x(:,1));
            K(Ind,i) = (1 - h);
        end   
    case 'Matern12' % Matern32 kernel
        K=zeros(N,M);
        for i=1:length(Param)
            K=K+CalculateDistance(x(:,i)/Param(i),y(:,i)/Param(i));
        end
        K=sqrt(K);
        K=exp(-K);
    case 'Matern32' % Matern32 kernel
        K=zeros(N,M);
        for i=1:length(Param)
            K=K+CalculateDistance(x(:,i)/Param(i),y(:,i)/Param(i));
        end
        K=sqrt(3)*sqrt(K);
        K=(1 + K).*exp(-K);
        K=K./sum(K,2);
    case 'Matern52' % Matern52 kernel
        K=zeros(N,M);
        for i=1:length(Param)
            K=K+CalculateDistance(x(:,i)/Param(i),y(:,i)/Param(i));
        end
        K=sqrt(5)*sqrt(K);
        K=(1 + K.*(1 + K/3)).*exp(-K);
    case 'Anisotropic-RBF' % Anisotropic RBF kernel
        D = Param; % Anisotropic scaling factor
        x = x*D;
        y = y*D;
        K = NormKernel(x,y,1,'RBF');
    case 'Laplace' % Laplace kernel
        K=zeros(N,M);
        for i=1:length(Param)
            K=K+sqrt(CalculateDistance(x(:,i),y(:,i)))./(2*Param(i)^2);
        end
        K=exp(-K);
    case 'Polynomial' % Polynomial kernel
        a = Param(1);   %  order
        b = Param(2);   %  constant
        K = ((x*y' + b).^a);
    case 'Linear' % Linear kernel
        K = x*y';
end
end
function Distance=CalculateDistance(XN,XM)
Positive=0; % Force distance>=0
Distance = bsxfun(@plus,bsxfun(@plus,sum(XN.^2,2),-2*XN*XM'),sum(XM.^2,2)');
if Positive
    Distance = max(0,Distance);
end
end
