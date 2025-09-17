
%% Loading the data
experiment='np01_003_000';
data=load(strcat('./data/',experiment,'.db'),'-mat'); % loading the log data

%% selecting top 40 based on the envlope sharpness
m=0.5;
[mask,index,env,Narea,Nmax]=getmasks(data.log,m); % envelopes are returned as env
select=index(1:40);
disp('selected=')
disp(select)

%% running the optimization algorithm
A0=[4.5,0 ; 0,4.5];
[A,fval,exitflag,output,grad]=EYE_CORRECT(data.log,select,mask,A0);

%% Evaluating RFs before and after correction
rfsb=RF(data.log,select,[],1); % before correction
rfsa=RF(data.log,select,A,1); % after correction

% rearranging the data
rfsa= permute(rfsa,[1 2 4 3]);
rfsb= permute(rfsb,[1 2 4 3]);
env= permute(env,[1 3 2]); 


%% visualizing the RFs
ftsz=18;
N=length(select);
disp(select);

for i=1:N

    ii=select(i);
    rfb=squeeze(mean(squeeze(rfsb(:,i,:,:)),1));
    rfa=squeeze(mean(squeeze(rfsa(:,i,:,:)),1));
    bmax=max(rfb(:));
    amax=max(rfa(:));
    cmax=max(bmax,amax);

    bmin=min(rfb(:));
    amin=min(rfa(:));
    cmin=min(bmin,amin);
    figure(1)
    ax1=subplot(2,1,1);
    imagesc(rfb);
    axis equal
    axis tight
    % axis off
    set(gca,'xticklabel',[])
    set(gca,'yticklabel',[])
    set(gca,'XTick',[])
    set(gca,'YTick',[])
    colormap(ax1,parula)
    clim([cmin cmax])
    title('Naive','fontweight','bold','fontsize',ftsz*1.2)

    colorbar
    ax2=subplot(2,1,2);
    imagesc(rfa);
    axis equal
    axis tight
    % axis off
    set(gca,'xticklabel',[])
    set(gca,'yticklabel',[])
    set(gca,'XTick',[])
    set(gca,'YTick',[])
    colormap(ax1,parula)
    clim([cmin cmax])
    title('Corrected','fontweight','bold','fontsize',ftsz*1.2)
    colorbar
    pause;
end





%% Main Functions 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%
function [A,fval,exitflag,output,grad,mask]=EYE_CORRECT(DATA,select,mask,A0)

% The cost function is based on the norm^2 of the receptive
% fields. To evaluate the cost function the receptive field of each neuron i 
% is normalized by mean(r_{i,t})+0.001, then the norm is calculated and
% summed across neurons.


% Inputs are:

% DATA=log

% "select" is the neurons selected to do the optimization based on their
% receptive fields. If select = [] all the neurons will be selected 

% A0 is the initial guess for 2by2 matrix A. if A=[] then A0=[4,0 ; 0,4]
% will be used


% Outputs are:
% A, fval , exitflag , output , grad at the optimum point.

% Putting data to the format to pass to optimization algorithm
disp('Putting data to the format to pass to optimization algorithm');
ncells = size(DATA.signal{1},2);
nstim = size(DATA,1);
x_E=DATA.xpos;
y_E=DATA.ypos;

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
mask1=mask(select,:,:);
RF_size=[480,135];

% running the minimization algorithm
X0=[ A0(1,1) , A0(1,2) , A0(2,1) , A0(2,2) ];
disp('starting optimization')
tic
opt = optimoptions('fminunc','SpecifyObjectiveGradient',true,'Display','iter','PlotFcns',@optimplotfval,'MaxIter',150,'TolFun',1e-5,'TolX',1e-5,'OutputFcn', @outfun);
[X,fval,exitflag,output,grad] = fminunc( @(X)EYE_CORRECT_COSTSD(X,x_E,y_E,S_tr,r,RF_size,mask1) , X0 , opt );
toc

A=[X(1),X(2);X(3),X(4)];

end

%% Cost function and its gradient
function [ELL,DELL]=EYE_CORRECT_COSTSD(X,x_E,y_E,S_tr,Lambda_tr,RF_size,mask)

 
A(1,1)=X(1);
A(1,2)=X(2);
A(2,1)=X(3);
A(2,2)=X(4);

[RFS,DRFS11,DRFS12,DRFS21,DRFS22]=RF_ESTSD(Lambda_tr,S_tr,x_E,y_E,RF_size,A);

RFS=RFS.*mask;
DELL = -2 * [sum(RFS.*DRFS11,'all') , sum(RFS.*DRFS12,'all') , sum(RFS.*DRFS21,'all') , sum(RFS.*DRFS22,'all')]/(RF_size(1)*RF_size(2));
ELL = -norm(RFS(:))^2 /(RF_size(1)*RF_size(2));

end



%% Auxiliary Functions
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%

function U=KXLY(mon_size,k,l)
% this function evaluates \alpha_{x,y}
Xmon=mon_size(1);
Ymon=mon_size(2);
ix=(1:Xmon);
iy=(1:Ymon);
[X,Y]=meshgrid(ix,iy);
U=(2*pi/Xmon)*(k*X+l*Y);
U=U';
end



%%


function [Rcos,Rsin,DaRcos,DaRsin,DbRcos,DbRsin,DcRcos,DcRsin,DdRcos,DdRsin]=reorganize(Lambda_tr,S_tr,x_E,y_E,A,RF_size)

% This function reorganizes the cas and cas' to decouple the trial
% dependent part and evaluates different gradient parts as explained in the
% paper.
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

[Rcos,Rsin,DaRcos,DaRsin,DbRcos,DbRsin,DcRcos,DcRsin,DdRcos,DdRsin]=reorganize(Lambda_tr,S_tr,x_E,y_E,A,RF_size);
[N,~]=size(Lambda_tr);
Xmon=RF_size(1);
Ymon=RF_size(2);
k=unique(S_tr(:,1));
l=unique(S_tr(:,2));
p=unique(S_tr(:,3));


Nk=length(k);
Nl=length(l);
Np=length(p);


DRFS11=zeros(N,Xmon,Ymon);
DRFS12=zeros(N,Xmon,Ymon);
DRFS21=zeros(N,Xmon,Ymon);
DRFS22=zeros(N,Xmon,Ymon);
RFS=zeros(N,Xmon,Ymon);

for ik=1:Nk
    kk=k(ik);
    for il=1:Nl
        ll=l(il);
        U=KXLY(RF_size,kk,ll);
        for ip=1:Np
           

            Rcos1=squeeze(Rcos(ik,il,ip,:));
            Rsin1=squeeze(Rsin(ik,il,ip,:));
            DaRcos1=squeeze(DaRcos(ik,il,ip,:));
            DaRsin1=squeeze(DaRsin(ik,il,ip,:));
            DbRcos1=squeeze(DbRcos(ik,il,ip,:));
            DbRsin1=squeeze(DbRsin(ik,il,ip,:));
            DcRcos1=squeeze(DcRcos(ik,il,ip,:));
            DcRsin1=squeeze(DcRsin(ik,il,ip,:));
            DdRcos1=squeeze(DdRcos(ik,il,ip,:));
            DdRsin1=squeeze(DdRsin(ik,il,ip,:));
            

            RFS=RFS+ (squeeze(tensorprod(Rcos1,cas(U))) + squeeze(tensorprod(Rsin1,cas1(U))));

            DRFS11=DRFS11+ (squeeze(tensorprod(DaRcos1,cas(U))) + squeeze(tensorprod(DaRsin1,cas1(U))));
            DRFS12=DRFS12+ (squeeze(tensorprod(DbRcos1,cas(U))) + squeeze(tensorprod(DbRsin1,cas1(U))));
            DRFS21=DRFS21+ (squeeze(tensorprod(DcRcos1,cas(U))) + squeeze(tensorprod(DcRsin1,cas1(U))));
            DRFS22=DRFS22+ (squeeze(tensorprod(DdRcos1,cas(U))) + squeeze(tensorprod(DdRsin1,cas1(U))));



        end

    end

end




end


%%
function [mask,index,envelopes,Narea,N1]=getmasks(log,m)
rfs=RF(log,[],[],0);

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


function stop = outfun(x, optimValues, state)
stop = false;
A=[x(1), x(2) ; x(3) , x(4)];
disp('A=')
disp(A)
drawnow
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
