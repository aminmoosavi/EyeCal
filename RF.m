


function RFST=RF(DATA,select,A,flag)

% Inputs are:

% DATA=log

% "select" is the neurons selected to do the optimization based on their
% receptive fields. If select = [] all the neurons will be selected


% Outputs are:
% Receptive fields for select neurons at all times with corrections based on A

% flag=0 RFs are calculated for all times. falg=1 RFs are calculated for
% the average 2:7 times.

nstim = size(DATA,1);
ncells = size(DATA.signal{1},2);
if(isempty( select ))
    select=1:ncells;
end

if(isempty( A ))
    A=[0,0 ; 0,0];
end

x_E=DATA.xpos;
y_E=DATA.ypos;

S_tr=zeros(nstim,3);

S_tr(:,1)=DATA.kx;
S_tr(:,2)=DATA.ky;
S_tr(:,3)=DATA.sign;
RF_size=[480,135];
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
function RFS=RF_ESTS(Lambda_tr,S_tr,x_E,y_E,RF_size,A)

[Rcos,Rsin,~,~,~,~,~,~,~,~]=reorganize(Lambda_tr,S_tr,x_E,y_E,A,RF_size);
[N,~]=size(Lambda_tr);
Xmon=RF_size(1);
Ymon=RF_size(2);
k=unique(S_tr(:,1));
l=unique(S_tr(:,2));
p=unique(S_tr(:,3));


Nk=length(k);
Nl=length(l);
Np=length(p);


RFS=zeros(N,Xmon,Ymon);

for ik=1:Nk
    kk=k(ik);
    for il=1:Nl
        ll=l(il);
        U=KXLY(RF_size,kk,ll);
        for ip=1:Np
            Rcos1=squeeze(Rcos(ik,il,ip,:));
            Rsin1=squeeze(Rsin(ik,il,ip,:));            
            RFS=RFS+ (squeeze(tensorprod(Rcos1,cas(U))) + squeeze(tensorprod(Rsin1,cas1(U))));

        end

    end

end

end


%%

function [Rcos,Rsin,DaRcos,DaRsin,DbRcos,DbRsin,DcRcos,DcRsin,DdRcos,DdRsin]=reorganize(Lambda_tr,S_tr,x_E,y_E,A,RF_size)
a=A(1,1);
b=A(1,2);
c=A(2,1);
d=A(2,2);

idxnanx=isnan(x_E);
idxnany=isnan(y_E);
x_E(idxnanx)=0;
y_E(idxnany)=0;



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
function y=cas1(x)

y=sin(x)-cos(x);

end

%%

function y=cas(x)

y=sin(x)+cos(x);

end
%%
function U=KXLY(mon_size,k,l)

Xmon=mon_size(1);
Ymon=mon_size(2);
ix=(1:Xmon);
iy=(1:Ymon);
[X,Y]=meshgrid(ix,iy);
U=(2*pi/Xmon)*(k*X+l*Y);
U=U';
end
