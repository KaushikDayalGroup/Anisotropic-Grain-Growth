function [ngrains,orio,id] = gbm3d(nt,grains,dims,ori,id,rmlist,d,fftK,MyDir2)
%% Auxiliary vars.
global WORKSPACE Z LONGARRAY;

n1 = dims(1); % Size of computational grid.
n2 = dims(2); % Size of computational grid.
n3 = dims(3); % Size of computational grid.
%n1n2n3 = n1*n2*n3;

WORKSPACE = int32(-ones(dims)); % Needed for dilation.c.
LONGARRAY = int32(zeros(120000000,1)); % Needed for dilation.c.
Z = -ones(dims); % Another Workspace var.
                % Needed in ls2vf3D.c.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% MAIN TIME LOOP STARTS:
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for t=1:nt % Main time iteration starts.

  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  %% CONVOLUTION STEP:
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  N = size(grains,1); % Number of grains.
  % Find grains present in a nhd. of each grid pt.:
  presence = get_nhd_grains(grains,dims(1)*dims(2)*dims(3));
  lp= length(presence);

  shared = cell([N, 1]);
  sharedgroup = cell([N, 1]);
  indd = cell([N, 1]);
  gcval = cell([N, 1]);
  LN = zeros(N,1);
  LNgroup = zeros(N,1);
  k=1;
  %%
  a=cell(N,1);
    for k=1:N % Loop over grains.
        indd{k} = grains{k,1};
        shared{k}=presence(indd{k},1); %local id
        shared{k}=[shared{k}{:}];
        shared{k} = unique(shared{k});      
        LN(k)=length(shared{k}); 
        a{k}=zeros(LN(k),1);
        b=zeros(LN(k),1);
        for i=1:LN(k)
                [a{k}(i),b(i)] = find(d==id(shared{k}(i)));
            %d is the matrix of groups, with each group in one row
        end
        sharedgroup{k} = unique(a{k});
        LNgroup(k) = length(sharedgroup{k});
        
    end
    %%
    clear presence %As it uses a lot of memory and slows down the next computation
    %%
    
    for k=1:N
        k
        a1 = grains(k,:);
        ii=indd{k};
        [gnum,~] = find(d==id(k));
        sharedgrouppar = sharedgroup{k};
        apar=a{k};
        h = convolvegrain(a1{1},a1{2},dims, sharedgrouppar, gnum,fftK); % Calulate convolution.
        b=zeros(length(ii),LN(k));
        for i=1:LN(k)
                no=find(apar(i)==sharedgrouppar);
                b(:,i) = h{no}(ii);
        end
        gcval{k} = b;
    end 
    %%        
    for j=1:N
        grains{j,3} = gcval{j};
    end
    fname = sprintf('gval%d.mat', j);
    MyFile=strcat(MyDir2,fname);
    save(MyFile,'gcval', '-v7.3'); %Because computing gcval is the most expensive part of the code so store it!
%%    
    sharedtemp=cell([N+1,1]);
    sharedtemp{1} = 1;
    for i=2:N+1
        sharedtemp{i} = id(shared{i-1})+1;
    end
    % sharedtemp contains global temp ids and not locals
    
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  %% REDISTRIBUTION STEP:
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  grain0=cell(1,3);
  grain0{1,1}=rmlist;
  grain0{1,2}=ones(size(rmlist,1),1);
  grain0{1,3}=zeros(size(rmlist,1),1);
  tempgrains=[grain0;grains];
  tempid=[1;id+1];
  tempori=[0;ori];
  presence = get_nhd_grains(tempgrains,dims(1)*dims(2)*dims(3));
%% 
  updatelevelsetdata(presence,tempgrains,tempid,tempori,sharedtemp);
  grains=tempgrains(2:end,:);

  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  %% REMOVE EMPTY GRAINS:
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  
  % Not necessary, but speeds up the algorithm.
  [grains,id] = removeemptygrains(grains,dims,id);

  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  %% REFRESH GRAIN BUFFERS:
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  % Level set data for each grain extends beyond the interface,
  % to a tubular neighborhood. That neighborhood has to be
  % updated once the interface has (potentially) moved.
  N = size(grains,1); % Number of grains.
  for k=1:N % Loop over grains.
    ind = grains{k,1}; % Pixels within a nhd. of grain.
    val = grains{k,2}; % Lev. set. vals. at those pixels.
    cval = grains{k,3}; % Convolution vals. at those pixels.
    Z(ind) = val;      % Lev. set. representation on grid.
    posind = ind(val>0); % Pixels in the interior of grain.
    [x,y,z] = ind2sub(dims,posind);
    [x3,y3,z3] = dilation_fixedbd(int32(x),int32(y),int32(z),5,WORKSPACE,LONGARRAY); % Dilation.
    ind3 = sub2ind(dims,x3,y3,z3);   
    ind2 = rmdilation1(ind3,rmlist );
%     ind2=ind3;
    val2 = Z(ind2);    % Level set vals.
    Z(ind2) = -1;      % Reset Z back to all -1's.
%     Z(ind) = cval() - 1; % Convolution values - 1.
    cval2 = zeros(length(ind2),size(cval,2));
    for j=1:size(cval,2)
        Z(ind) = cval(:,j) - 1;
        cval2(:,j) = Z(ind2);   % Convolution vals - 1.
        Z(ind2) = -1;
    end   % Convolution vals - 1.
%     Z(ind2) = -1;      % Reset Z back to all -1's.
    grains{k,1} = ind2;   % Refresh grain's data structure.
    grains{k,2} = val2;   % Ditto.
    grains{k,3} = cval2 + 1; % Ditto.
  end % (for k). Loop over grains ends.

end % (for t). Main time iteration ends.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% MAIN TIME LOOP ENDS.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%
ngrains=grains;
% Order of grains may have changed; return new orientation array.
N = size(grains,1); % Number of grains.
orio=zeros(N,1);
for k=1:N % Loop over grains.
%   dataorio(k) = dataori(id(k));
  orio(k) = ori(id(k));
end

end