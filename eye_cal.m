%% This is the code for eye movement correction in the paper "Subspace reverse-correlation estimation of receptive fields during free viewing"
%% Loading the data
clc;
clear;
experiment='np01_003_000'; % Experiment name
data=load(strcat('./data/',experiment,'.db'),'-mat'); % loading the log data

% dimensions of the visual field (the monitor). Here we downsampled our monitor size by a factor of 8.
RF_size=[480,135];

%% Computing the frequency domain receptive fields (fRFs)

[fRFs,norms]=krnls(data.log); % Computing the frequency domain receptive fields (fRFs)
fRFs=mean(fRFs,3);

%% Selecting cells to pass to the optimization algorithm

% FIRST METHOD
% If the RFs are strongly distorted by eye movements such that they are not detectable we may select the cells based on the fRFs and their frequency.
% [poris,freqs,frfm]=freq_extract(fRFs); % Extracting preferred oreintations (poris), preferred spatial frequencies (freqs)
% index=find(frfm>0.01 & (freqs<16)); % finding the cells with large enough mean response and preferred frequency below a threshold
% select=index; % The select indices will be used for optimization.
% mask=[]; % empty mask allows optimization to be applied to all the visual field domain.

% SECOND METHOD
% If the RFs are visible we can select cells based on the Riesz envlope sharpness of the RFs.
% If needed a mask can be used to limit the optimizer to particular areas of the receptive fields.
% The following code gives us mask and selected neurons in mouse data.
m=0.5;
[mask,index,env,Narea,Nmax]=getmasks(data.log,m,RF_size); % Riesz envelopes are returned as env. Masks are seleced by thresholding the envelopes above 50 percent of their max value.
select=index(1:40);
disp('selected=')
disp(select)
mask=[]; % empty mask allows optimization to be applied to all the visual field domain.

% THIRD METHOD
% You can select the neurons and masks by hand.
% select must be a one dimensional array with indices of selected neurons
% mask must be a 3D array of zeros and ones with first index for neurons and other two
% equivalent to dimensions of the screen. The optimization will be applied to mask==1.
%% running the optimization algorithm
A0=[4.5,0 ; 0,4.5]; % initial guess of the A matrix.
[A,fval,exitflag,output,grad]=EYE_CORRECT(data.log,select,mask,A0,RF_size);

%% Evaluating all RFs before and after correction
tic
rfsb=RECEPTIVE_FIELDS(data.log,[],[],0,RF_size); % before correction
rfsa=RECEPTIVE_FIELDS(data.log,[],A,0,RF_size); % after correction
toc
 

% rearranging the data
% rfsa= permute(rfsa,[1 2 4 3]);
% rfsb= permute(rfsb,[1 2 4 3]);



%% visualizing the RFs
ftsz=18;
N=length(select);
% disp(select);

for i=1:N

    ii=select(i);
    % rfb=squeeze(mean(squeeze(rfsb(:,ii,:,:)),1));
    % rfa=squeeze(mean(squeeze(rfsa(:,ii,:,:)),1));

    rfb=squeeze(squeeze(rfsb(ii,:,:)))';
    rfa=squeeze(squeeze(rfsa(ii,:,:)))';
    bmax=max(rfb(:));
    amax=max(rfa(:));
    cmax=max(bmax,amax);

    
    bmin=min(rfb(:));
    amin=min(rfa(:));
    cmin=min(bmin,amin);
    figure(1)
    ax1=subplot(2,2,2);
    imagesc(rfb);
    axis equal
    axis tight
    % axis off
    set(gca,'xticklabel',[])
    set(gca,'yticklabel',[])
    set(gca,'XTick',[])
    set(gca,'YTick',[])
    set(gca,'YDir','normal')

    colormap(ax1,parula)
    clim([cmin cmax])
    title('Naive','fontsize',ftsz)

    colorbar
    ax2=subplot(2,2,4);
    imagesc(rfa);
    axis equal
    axis tight
    % axis off
    set(gca,'xticklabel',[])
    set(gca,'yticklabel',[])
    set(gca,'XTick',[])
    set(gca,'YTick',[])
    set(gca,'YDir','normal')
    colormap(ax1,parula)
    clim([cmin cmax])
    title('Corrected','fontsize',ftsz)
    colorbar
     subplot(1,2,1)
    imagesc(squeeze(fRFs(:,:,ii)))
    axis equal
    axis tight
    xlabel('$k$','interpreter','latex','fontsize',ftsz)
    ylabel('$l$','interpreter','latex','fontsize',ftsz)
    title('fRF','fontsize',ftsz)
    pause;

   
end





%%


function [A,fval,exitflag,output,grad,mask]=EYE_CORRECT(DATA,select,mask,A0,RF_size)

% The cost function is based on the norm^2 of the receptive
% fields. 


% Inputs are:

% DATA=log

% "select" is the neurons selected to do the optimization based on their
% receptive fields. If select = [] all the neurons will be selected 

% mask allows focusing on the receptive fields is a binary array of zeros
% and ones. zero outside of the RF domain and one in the RF domain of each
% neuron. If mask=[] the code uses all the visual field to calibrate eye
% movements

% A0 is the initial guess for 2by2 matrix A. if A=[] then A0=[4,0 ; 0,4]
% will be used


% Outputs are:
% A, fval , exitflag , output , grad at the optimum point.

% Putting data to the format to pass to optimization algorithm
disp('Putting data to the format to pass to optimization algorithm');
ncells = size(DATA.signal{1},2);
nstim = size(DATA,1);
x_E=DATA.xpos-mean(DATA.xpos,'omitnan');
y_E=DATA.ypos-mean(DATA.ypos,'omitnan');

if(sum(isnan(x_E))+sum(isnan(y_E))>0)
    disp('--------------------------------')
    warning('There are trials with NAN values for eye positions that are set to zero')
    disp('--------------------------------')
end

idxnanx=isnan(x_E);
idxnany=isnan(y_E);
x_E(idxnanx)=0;
y_E(idxnany)=0;



S_tr=zeros(nstim,3);

S_tr(:,1)=DATA.kx;
S_tr(:,2)=DATA.ky;
S_tr(:,3)=DATA.sign;

R=cat(3,DATA.signal{:});
Lambda_tr=squeeze(mean(R(2:7,:,:),1));


if(isempty( select ))
    select=1:ncells;
end

if(isempty( A0 ))
    A0=[4,0 ; 0,4];
end

r=Lambda_tr(select,:);

if(isempty(mask))
    mask1=[];
else
    mask1=mask(select,:,:);
end



% running the minimization algorithm
X0=[ A0(1,1) , A0(1,2) , A0(2,1) , A0(2,2) ];
disp('starting optimization')
tic
opt = optimoptions('fminunc','SpecifyObjectiveGradient',true,'Display','iter','PlotFcns',@optimplotfval,'MaxIter',150,'TolFun',1e-5,'TolX',1e-5,'OutputFcn', @outfun);
[X,fval,exitflag,output,grad] = fminunc( @(X)EYE_CORRECT_COSTSD(X,x_E,y_E,S_tr,r,RF_size,mask1) , X0 , opt );
toc

A=[X(1),X(2);X(3),X(4)];

end




%% Functions

function U=KXLY(mon_size,k,l)
% meshgrid phase of the cas function in the visual field
Xmon=mon_size(1);
Ymon=mon_size(2);
ix=(1:Xmon);
iy=(1:Ymon);
[X,Y]=meshgrid(ix,iy);
U=(2*pi/Xmon)*(k*X+l*Y);
U=U';
end



%% This is a function that decouples the space and time components of the cas


function [Rcos,Rsin,DaRcos,DaRsin,DbRcos,DbRsin,DcRcos,DcRsin,DdRcos,DdRsin]=reorganize(Lambda_tr,S_tr,x_E,y_E,A,RF_size)
a=A(1,1);
b=A(1,2);
c=A(2,1);
d=A(2,2);

[N,~]=size(Lambda_tr);
k=unique(S_tr(:,1));
l=unique(S_tr(:,2));
p=unique(S_tr(:,3));


Nk=length(k);
Nl=length(l);
Np=length(p);

XRF=RF_size(1);

Rsin=zeros(Nk,Nl,Np,N);
Rcos=zeros(Nk,Nl,Np,N);
DaRcos=zeros(Nk,Nl,Np,N);
DaRsin=zeros(Nk,Nl,Np,N);
DbRcos=zeros(Nk,Nl,Np,N);
DbRsin=zeros(Nk,Nl,Np,N);
DcRcos=zeros(Nk,Nl,Np,N);
DcRsin=zeros(Nk,Nl,Np,N);
DdRcos=zeros(Nk,Nl,Np,N);
DdRsin=zeros(Nk,Nl,Np,N);

for ik=1:Nk
    kkk=k(ik);

    for il=1:Nl
        lll=l(il);
        for ip=1:Np
            pp=p(ip);
            idx=find(S_tr(:,1)==kkk & S_tr(:,2)==lll & S_tr(:,3)==pp);
            kk=kkk*2*pi/XRF;
            ll=2*pi*lll/XRF;

            xx=x_E(idx);
            yy=y_E(idx);
            
            Rcos(ik,il,ip,:)=pp*tensorprod( Lambda_tr(:,idx),cos(   kk*(a*xx+b*yy)+ll*(c*xx+d*yy)   ) , 2,1 );
            Rsin(ik,il,ip,:)=pp*tensorprod( Lambda_tr(:,idx),sin(   kk*(a*xx+b*yy)+ll*(c*xx+d*yy)   ) ,  2,1);

            DaRcos(ik,il,ip,:)=-pp*tensorprod( Lambda_tr(:,idx), (kk*xx) .* sin(   kk*(a*xx+b*yy)+ll*(c*xx+d*yy)   ), 2,1);
            DaRsin(ik,il,ip,:)= pp*tensorprod( Lambda_tr(:,idx), (kk*xx) .* cos(   kk*(a*xx+b*yy)+ll*(c*xx+d*yy)   ), 2,1);

            DbRcos(ik,il,ip,:)=-pp*tensorprod( Lambda_tr(:,idx), (kk*yy) .* sin(   kk*(a*xx+b*yy)+ll*(c*xx+d*yy)   ), 2,1);
            DbRsin(ik,il,ip,:)= pp*tensorprod( Lambda_tr(:,idx), (kk*yy) .* cos(   kk*(a*xx+b*yy)+ll*(c*xx+d*yy)   ), 2,1);

            DcRcos(ik,il,ip,:)=-pp*tensorprod( Lambda_tr(:,idx), (ll*xx) .* sin(   kk*(a*xx+b*yy)+ll*(c*xx+d*yy)   ), 2,1);
            DcRsin(ik,il,ip,:)= pp*tensorprod( Lambda_tr(:,idx), (ll*xx) .* cos(   kk*(a*xx+b*yy)+ll*(c*xx+d*yy)   ), 2,1);

            DdRcos(ik,il,ip,:)=-pp*tensorprod( Lambda_tr(:,idx), (ll*yy) .* sin(   kk*(a*xx+b*yy)+ll*(c*xx+d*yy)   ), 2,1);
            DdRsin(ik,il,ip,:)= pp*tensorprod( Lambda_tr(:,idx), (ll*yy) .* cos(   kk*(a*xx+b*yy)+ll*(c*xx+d*yy)   ), 2,1);



      
            
        end
    end
end


end







%% 
function [RFS,DRFS11,DRFS12,DRFS21,DRFS22]=RF_ESTSD(Lambda_tr,S_tr,x_E,y_E,RF_size,A)
% This function computes the variables that are used to evaluate the objective function and its gradient 
[Rcos,Rsin,DaRcos,DaRsin,DbRcos,DbRsin,DcRcos,DcRsin,DdRcos,DdRsin]=reorganize(Lambda_tr,S_tr,x_E,y_E,A,RF_size);
% [N,~]=size(Lambda_tr);
Xmon=RF_size(1);
Ymon=RF_size(2);
k=unique(S_tr(:,1));
l=unique(S_tr(:,2));
p=unique(S_tr(:,3));


Nk=length(k);
Nl=length(l);
% Np=length(p);



casU=zeros(Nk,Nl,Xmon,Ymon);
cas1U=zeros(Nk,Nl,Xmon,Ymon);

Rcos1=squeeze(sum(Rcos,3));
Rsin1=squeeze(sum(Rsin,3));

DaRcos1=squeeze(sum(DaRcos,3));
DaRsin1=squeeze(sum(DaRsin,3));

DbRcos1=squeeze(sum(DbRcos,3));
DbRsin1=squeeze(sum(DbRsin,3));


DcRcos1=squeeze(sum(DcRcos,3));
DcRsin1=squeeze(sum(DcRsin,3));

DdRcos1=squeeze(sum(DdRcos,3));
DdRsin1=squeeze(sum(DdRsin,3));

for ik=1:Nk
    kk=k(ik);
    for il=1:Nl
        ll=l(il);
        U=KXLY(RF_size,kk,ll);
        casU(ik,il,:,:)=cas(U);
        cas1U(ik,il,:,:)=cas1(U);

    end

end

RFS=tensorprod(Rcos1,casU,[1,2],[1,2]) + tensorprod(Rsin1,cas1U,[1,2],[1,2]);
DRFS11=tensorprod(DaRcos1,casU,[1,2],[1,2]) + tensorprod(DaRsin1,cas1U,[1,2],[1,2]);
DRFS12=tensorprod(DbRcos1,casU,[1,2],[1,2]) + tensorprod(DbRsin1,cas1U,[1,2],[1,2]);
DRFS21=tensorprod(DcRcos1,casU,[1,2],[1,2]) + tensorprod(DcRsin1,cas1U,[1,2],[1,2]);
DRFS22=tensorprod(DdRcos1,casU,[1,2],[1,2]) + tensorprod(DdRsin1,cas1U,[1,2],[1,2]);




end

%%

function RFS=RF_ESTS(Lambda_tr,S_tr,x_E,y_E,RF_size,A)

[Rcos,Rsin,~,~,~,~,~,~,~,~]=reorganize(Lambda_tr,S_tr,x_E,y_E,A,RF_size);

Xmon=RF_size(1);
Ymon=RF_size(2);
k=unique(S_tr(:,1));
l=unique(S_tr(:,2));
p=unique(S_tr(:,3));


Nk=length(k);
Nl=length(l);

casU=zeros(Nk,Nl,Xmon,Ymon);
cas1U=zeros(Nk,Nl,Xmon,Ymon);
Rcos1=squeeze(sum(Rcos,3));
Rsin1=squeeze(sum(Rsin,3));
for ik=1:Nk
    kk=k(ik);
    for il=1:Nl
        ll=l(il);
        U=KXLY(RF_size,kk,ll);
        casU(ik,il,:,:)=cas(U);
        cas1U(ik,il,:,:)=cas1(U);

    end

end



RFS=tensorprod(Rcos1,casU,[1,2],[1,2]) + tensorprod(Rsin1,cas1U,[1,2],[1,2]);


end



%%


function stop = outfun(x, optimValues, state)
stop = false;
A=[x(1), x(2) ; x(3) , x(4)];
disp('A=')
disp(A)
drawnow
end

%%
function [ELL,DELL]=EYE_CORRECT_COSTSD(X,x_E,y_E,S_tr,Lambda_tr,RF_size,mask)
% This function returns 
% ELL = -E(A) where E(A) is our opjective function
% DELL = the gradient 
 
A(1,1)=X(1);
A(1,2)=X(2);
A(2,1)=X(3);
A(2,2)=X(4);

[RFS,DRFS11,DRFS12,DRFS21,DRFS22]=RF_ESTSD(Lambda_tr,S_tr,x_E,y_E,RF_size,A);

if(~isempty(mask))
    RFS=RFS.*mask;
end
DELL = -2 * [sum(RFS.*DRFS11,'all') , sum(RFS.*DRFS12,'all') , sum(RFS.*DRFS21,'all') , sum(RFS.*DRFS22,'all')]/(RF_size(1)*RF_size(2));
ELL = -norm(RFS(:))^2 /(RF_size(1)*RF_size(2));

end



%%
function y=cas1(x)

y=sin(x)-cos(x);

end

%%

function y=cas(x)

y=sin(x)+cos(x);

end

%%
function [mask,index,envelopes,Narea,N1]=getmasks(log,m,RF_size)
% This function gets the data (log table), the threshold for masking 0<m<1
% and RF size [length , height].

% Returns: 
% masks, 
% index = cell indices sorted according to the naive RF sharpness, 
% envelopes = The Riesz Transforms
% Narea = number of connected components after thresholding
% N1 = Number of RFs with only one connected component
rfs=RECEPTIVE_FIELDS(log,[],[],0,RF_size);

[N,~,~]=size(rfs);
mask=zeros(size(rfs));
Cost=zeros(N,1);
envelopes=zeros(size(rfs));
Narea=zeros(N,1);
for i=1:N
    E=rieszEnvelope(squeeze(rfs(i,:,:)));
    Eth=m*max(E,[],'all');
    % Eth=mean(E,'all')+1*std(E,[],'all');
    idx=E>Eth;
 
    
    CC = bwconncomp(idx);
    p = regionprops(CC,"Area");
    [maxArea,maxIdx] = max([p.Area]);
    [areas,IDX]=sort([p.Area],'descend');
    Narea(i)=CC.NumObjects;
    idx2= cc2bw(CC,ObjectsToKeep=maxIdx); 
    idx2 = imdilate(idx2,strel("disk",8));
    ratio=maxArea/sum(areas);



    envelopes(i,:,:)=E;
    mask(i,:,:)=idx2;


    
    Cost(i)=ratio*sum(idx2.*squeeze(rfs(i,:,:).^2),'all');
    Cost(i)=Cost(i)/(100000^(Narea(i)-1));
end


ind=Narea==1;
N1=sum(ind);

[~,index]=sort(Cost,'descend');



end



%%



function RFST=RECEPTIVE_FIELDS(DATA,select,A,flag,RF_size)

% Inputs are:

% DATA=log

% "select" is an array indices of the neurons you want their receptive fields.
% If select = [] all the neurons will be selected


% Outputs are:
% rfsb and rfsa
% Receptive fields for neurons with and without corrections based on A.

% flag=0 RFs are calculated for all time bins. falg=1 RFs are calculated for
% the average 2:7 time bins.

nstim = size(DATA,1);
ncells = size(DATA.signal{1},2);
if(isempty( select ))
    select=1:ncells;
end

if(isempty( A ))
    A=[0,0 ; 0,0];
end

x_E=DATA.xpos-mean(DATA.xpos,'omitnan');
y_E=DATA.ypos-mean(DATA.ypos,'omitnan');

idxnanx=isnan(x_E);
idxnany=isnan(y_E);
x_E(idxnanx)=0;
y_E(idxnany)=0;

S_tr=zeros(nstim,3);

S_tr(:,1)=DATA.kx;
S_tr(:,2)=DATA.ky;
S_tr(:,3)=DATA.sign;

Xmon=RF_size(1);
Ymon=RF_size(2);

LAMBDA=cat(3,DATA.signal{:});
LAMBDA=LAMBDA(:,select,:);

[T,N,~]=size(LAMBDA);


if flag==1
    RFST=zeros(T,N,Xmon,Ymon);

    for t=1:T

        Lambda_tr=squeeze(LAMBDA(t,:,:));
        RFS=RF_ESTS(Lambda_tr,S_tr,x_E,y_E,RF_size,A);

        RFST(t,:,:,:)=RFS;

    end

elseif(flag==0)
    Lambda_tr=squeeze(mean(LAMBDA(2:7,:,:),1));
    RFST=RF_ESTS(Lambda_tr,S_tr,x_E,y_E,RF_size,A);


end

end



%%

function E = rieszEnvelope(X)
% RIESZENVELOPE  Local envelope via the monogenic signal
%   E = rieszEnvelope(X) returns sqrt( X.^2 + R1.^2 + R2.^2 ),
%   where [R1,R2] are the Riesz transform of X.
% compute Riesz components
[M,N] = size(X);
F = fft2(X);
ux = ifftshift((-floor(N/2):ceil(N/2)-1)/N);
uy = ifftshift((-floor(M/2):ceil(M/2)-1)/M).';
[U,V] = meshgrid(ux,uy);
K = sqrt(U.^2+V.^2)+eps;
H1 = -1i*(U./K);
H2 = -1i*(V./K);
R1 = ifft2(F.*H1);
R2 = ifft2(F.*H2);
% envelope = magnitude of [X, R1, R2]
E = sqrt( abs(X).^2 + abs(R1).^2 + abs(R2).^2 );
end

%% 


function [Rt,norms]=krnls(log)
% This function gets the log data and returns the frequency domain RFs (fRFs) and their norms.

% "norms" is a two dimensional array. 
% First index -> number of cells
% Second index -> time


% Rt are kernels for all conditions 

maxk = max(log.kx);           % the maximum number of cycles per screen

sz = size(log.signal{1});     % size of the signal [time bins x number of cells] 

Rt = zeros([2*maxk+1 2*maxk+1 sz]);  % responses all 

for i = 1:size(log,1)                % lets accumulate the resposnes and compute the mean at the end 
    
    kx = log.kx(i)+maxk+1;           % pick (kx,ky) for this trial 
    ky = log.ky(i)+maxk+1;           % shift the index so it starts at 1 
    Rt(kx,ky,:,:) = squeeze(Rt(kx,ky,:,:)) + log.signal{i};
   
end




ncells = size(log.signal{1},2);

% filter 

sigma = 3;  % Little bit of smoothing the kernels 
norms=zeros(ncells,10);
for n=1:ncells
    for i=1:10

        Rt(:,:,i,n) = imgaussfilt(squeeze(Rt(:,:,i,n)),sigma);
        krnl=squeeze(Rt(:,:,i,n));
        norms(n,i)=norm(krnl(:));
       
        

    end
end


end


%% Extracting preferred frequency and orientation of the RFs


function [poris,freqs,frfm]=freq_extract(Rt)
% Tis function gets the fRFs as input and returns
% poris = preferred orientations
% freqs = preferred spatial frequencies
% frfm = mean of the fRFs

[~,~,N]=size(Rt);
poris=zeros(N,1);
freqs=zeros(N,1);
frfm=zeros(N,1);


t=1;
for i=1:N

    krnl=squeeze(Rt(:,:,i));
    krnl=0.5*(krnl+flip(flip(krnl,1),2));
    [rs, cs]=find(krnl==max(krnl(:)));


    if(length(cs)>1)
        thetaa=(180/pi)*atan( (cs(2)-cs(1))/(rs(2)-rs(1)));
        nuu=0.5*sqrt( (cs(2)-cs(1))^2+(rs(2)-rs(1))^2);
    else
        thetaa=2*pi;
        nuu=0;
    end


    
    poris(i,t)=thetaa;
    freqs(i,t)=nuu;
    frfm(i,t)=mean(krnl(:));


end



end

