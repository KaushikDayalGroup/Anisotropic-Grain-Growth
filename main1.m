%% input
no=20;
do=1;
grain_number =size(grains,1);

ori = 2*pi*rand(grain_number,1); %orientation generate
id=(1:grain_number)';
evol=cell(no,3);
evol{1,1}=grains;
evol{1,2}=ori;
evol{1,3}=id;
avegrainsize = zeros(no,1); 
avegrainsize(1) = (dims(1)*dims(2)*dims(3) - size(rmlist,1))/(size(grains,1));
avegrainsize_goal= 4.8470e+03; %This is from An5_ou2 data 
%%
fftK=cell([30,1]);

MyDir='C:\Users\k_nag\GrainGrowSimu\Runs\Ni2D\KERNEL_An4_out2\';

for i=1:30
    i;
    fname1 = sprintf('fftK%d.mat', i);
    MyFile1=strcat(MyDir,fname1);
    K1 = struct2cell(load(MyFile1)); fftK{i} = (K1{1})/3.880620e+01;
    mmax=max(max(max(fftK{i})));
    mmin=min(min(min(fftK{i})));
    fprintf('For grain %d, maximum of fftkernel is %d and the minimum is %d. \n',i,mmax,mmin)
end
%%
fname1 = sprintf('OriFamily');
MyFile1=strcat(MyDir,fname1);
K1 = struct2cell(load(MyFile1)); OriFamily = (K1{1});

%% Orientatoin group list
d=zeros(30,200);
for i=1:30
    l=length(OriFamily{i,1});
    d(i,1:l) = OriFamily{i,1};
end

%% iteration
tic
MyDir2='C:\Users\k_nag\GrainGrowSimu\Runs\Ni2D\KERNEL_An4_out2\';
for i=2:no
    i
    [evol{i,1},evol{i,2},evol{i,3}] = gbm3diso(do,dt,grains,dims,ori,id,rmlist);    
    grains=evol{i,1};
    id=evol{i,3};
    avegrainsize(i) = (dims(1)*dims(2)*dims(3) - size(rmlist,1))/(size(grains,1));
    m=avegrainsize(i);
    (size(grains,1))
    if avegrainsize(i)>=avegrainsize_goal
        display('Done!');
%         pause
    end
    fname = sprintf('timestep%d.mat', i);
    MyFile=strcat(MyDir2,fname);
    timestep=evol(i,:);
    save(MyFile,'timestep', '-v7.3');
    evol{i,1}=[]; evol{i,2}=[]; evol{i,3}=[]; 
end
time=toc
