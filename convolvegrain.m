function cval = convolvegrain(ind,val,dims,common,current_grain,fftK)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Convolves the characteristic function of a grains whose
% level set representation is stored in "ind" and "val" input vars
% with the appropriate anisotropic kernel based on its orientation "fftK".
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
ZP = -ones(dims);
n1 = dims(1); % Size of computational grid.
n2 = dims(2); % Size of computational grid.
n3 = dims(3); % Size of computational grid.

%% Convert level set data to volume of fluid representation:
[x,y,z] = ind2sub(dims,ind);
vf = ls2vf3D(int32(x),int32(y),int32(z),val,ZP,dims(1),dims(2),dims(3));
%% Carry out the convolution:
Ku = zeros(dims);
Ku(ind) = vf; % Characteristic function of the union of grains.
M=fftn(Ku);
cval=cell([length(common),1]);
%%
K2 = fftK{current_grain};

for ij=1:length(common)
        K1 = fftK{common(ij)};
        KERNEL = (K1 + K2)./2;
        cval{ij} = real(fftn(conj(M .* KERNEL)))/n1/n2/n3;
end