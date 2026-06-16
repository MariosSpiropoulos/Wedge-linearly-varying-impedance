function nvals = theta_equation(alpha,Z1,Z2,Nmodes,theta)
N = length(theta);
[x,A]=chebdif(N,2);
D2 = squeeze(A(:,:,2));
D = squeeze(A(:,:,1));
%size(x)
fac = pi-alpha;

x = x*fac;
D = D/fac;
D2 = D2/fac^2;

L = D2;
F = -eye(size(L));

%Z = 1*1i; % = Z/k

% BCs
if(isfinite(Z2))
    L(1,:) = 0; L(1,:) = 1i*Z2*D(1,:); L(1,1) = L(1,1) + 1; F(1,:) = 0;
else
    L(1,:) = 0; L(1,:) = D(1,:); F(1,:) = 0;
end
if(isfinite(Z1))
    L(N,:) = 0; L(N,:) = -1i*Z1*D(N,:); L(N,N) = L(N,N) + 1; F(N,:) = 0;
else
    L(N,:) = 0; L(N,:) = -D(N,:); F(N,:) = 0;
end
% solve
[V,lambda]=eig(L,F);
lambda = diag(lambda);
pos = find(isfinite(lambda));
lambda = lambda(pos); V = V(:,pos);

[~,ord]=sort(real(lambda));
lambda = lambda(ord);
V = V(:,ord);

kt = sqrt(lambda);
nvals = kt(1:Nmodes);
end


