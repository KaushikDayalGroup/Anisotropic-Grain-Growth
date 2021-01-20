% load('Ni_An0_simu_CVZ.mat', 'rmlist', 'dims')
% i=3; % iteration step

grains=timestep{1,1};
id=timestep{1,3};
%% decrease the size of the domain to  the original old one
grain0=cell(1,3);
grain0{1,1}=rmlist;
grain0{1,2}=ones(size(rmlist,1),1);
grain0{1,3}=zeros(size(rmlist,1),1);
grains=[grain0;grains];
id=[0;id];
dimsold = [395 404 84];
P = -ones(dimsold);
N = size(grains,1); % Number of grains.
for k=2:N % Loop over grains.
    indbig = grains{k,1}; % Pixels within a nhd. of grain.
    [d1, d2, d3] = ind2sub(dims,indbig);
    ind = sub2ind(dimsold, d1-2, d2, d3);
    grains{k,1} = ind;
    val = grains{k,2}; % Lev. set. vals. at those pixels.
    posind = ind(val>0); % Pixels in the interior of grain.
    P(posind) = id(k);        
end
k=1;
ind = rmlistold; % Pixels within a nhd. of empty grain in original dimension
posind = ind; % Pixels in the interior of grain.
P(posind) = id(k);
dataid=zeros([1,dimsold],'int32');


dataid(1,:,:,:)=P;

h5write('C:\Users\k_nag\GrainGrowSimu\Runs\output\aniso\Aniso300.dream3d',...
    '/DataContainers/ImageDataContainer/CellData/FeatureIds',dataid);
