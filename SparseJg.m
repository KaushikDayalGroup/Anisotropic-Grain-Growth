function [Jg]=SparseJg(dt,dims)

n1 = dims(1); % Size of computational grid.
n2 = dims(2); % Size of computational grid.
n3 = dims(3); % Size of computational grid.

[x,y,z]=meshgrid([1:n1],[1:n2],[1:n3]);
x=permute(x,[2 1 3]);
y=permute(y,[2 1 3]);
z=permute(z,[2 1 3]);
x=((x-n1/2-1)/n1/dt^0.5); y=((y-n2/2-1)/n2/dt^0.5); z=((z-n3/2-1)/n3/dt^0.5);

N=40; %number of point for Guass Quadrature integration

atheta=0; btheta=pi;
aphi=0; bphi=pi;
[pointt,wt]=lgwt(N,atheta,btheta);
[pointp,wp]=lgwt(N,aphi,bphi);
[phi1,theta1]=meshgrid(pointp,pointt);
eps=0.25;
h=n1*n2*n3*N*N;
%Approximation for Dawson function
dd=@(x) x.*(1 + 0.1048210880*x.^2 + 0.0423601431*x.^4 + 0.0072584274*x.^6+ 0.0005045886*x.^8 + 0.0001777233*x.^10)...
    ./(1 + 0.7713501425*x.^2 + 0.2908548286*x.^4 + 0.0693596200*x.^6+ 0.0139927487*x.^8 + 0.0008303974*x.^10 + 2*0.0001777233*x.^12);

%%
parfor m=1:n1
    l{m}=sparse(h,1);
    m
    for n=1:n2
        for p=1:n3
            aa=x(m,n,p)^2+y(m,n,p)^2+z(m,n,p)^2;
            
            dotxv=x(m,n,p)*sin(phi1).*cos(theta1)+y(m,n,p)*sin(phi1).*sin(theta1)+z(m,n,p)*cos(phi1);
            b=((1-eps^2)^0.5)*dotxv/2/eps;
            a=exp(-aa./4)./((eps^2).*((4*pi)^1.5)).*(1-(sqrt(1-eps^2))/eps*dotxv.*dd(b));
            a=a.*(abs(a)>1e-10);
            d=N*N*(m-1)+N*N*n1*(n-1)+N*N*n1*n2*(p-1);
            l{m}=l{m}+sparse(d+1:d+N*N,1,a(:),h,1);
        end 
    end 
end 

Jg=spalloc(n1*n2*n3*N*N,1,0);
for i=1:n1
    Jg=Jg+l{i};
end
end