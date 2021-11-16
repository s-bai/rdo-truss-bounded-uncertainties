% A non-probabilistic reliable based truss optimization(NRBTO) MATLAB code
% Written by Song BAI, 2010.08.28
%
% This program aims at finding the maximum reliable index 
% with uncertain parameters under the constraint of total weight 
% using aggregation function method  
%
% Modified by Song BAI, 2011.01.04
% Uncertain Loads Transformation is added
% Modified by Song BAI, 2011.03.28
% This version considers the uncertainty of Young's modulus
% 
% Name explanation:
% tenbar102: "tenbar" represents the structure to optimize
%            "102" represents the 1st project and the 2nd stage
% 
% Input: u_c--> Vertical displacement constraint on the right bottom node
%        v_c--> Total weight constraint
% Modified by Song BAI, 2011.04.02
% This version considers the uncertainty of design variables, 
% namely the uncertainties of bar sectional areas.
% The five-bar truss is considered here.
% Modified by Song Bai, 2011.4.22
% This modified version is for verification of analytical sensitivity


function fivebar105(agg_P, u1_c , u2_c, v_c)
global m nn i0 j0 l T...
       e F var MatTran UncertLNum...
       u1_constr u2_constr volume_constr...
       U_nom ; % These global variables are read-only
global a; % These global variables are REWRITABLE

tic % Timing start

u1_constr=u1_c; 
u2_constr=u2_c ; 
volume_constr=v_c; 
P=agg_P ; % Penalty factor for aggregate function

% READ DATA FROM INPUT FILE
[n,m,nc,nn,x,y,ihl,ihr,e,a,px,py,UncertLNum,W,var,rho]=readfile();

%%% Parameters for GCMMA %%%
constr_num=1; design_num=m;
epsimin=0.00001;
amin=0.01*ones(design_num,1);
amax=40*ones(design_num,1);
low=amin;
upp=amax;
raa0=0.01;
raa=0.01*ones(constr_num,1);
a0=1;
para_a=zeros(constr_num,1);
c=100*ones(constr_num,1);  
d=1.0*ones(constr_num,1);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Define loads
F=zeros(nn,1); 
F(1:2:2*(n-nc),1)=px(1:n-nc);
F(2:2:2*(n-nc),1)=py(1:n-nc);

% Preprocess for FEA
[i0,j0,l,T]=preprocess(m,x,y,ihl,ihr,nc);

% Transform Matrix
[Q,Lam]=eig(W);
LamSqtInv=diag(1./diag(sqrt(Lam)));
MatTran=Q*LamSqtInv;

change=1.0; etau1_old=1.0; etau2_old=1.0; 
iter=0; itermax=50; % Upper limit for iteration steps
etahistory=zeros(1,itermax);
etahistoryu1=zeros(1,itermax);
etahistoryu2=zeros(1,itermax);
A=a;

% Start iteration
while change>0.001 && iter< itermax
    iter=iter + 1;
    
    % Calculate disp. under nominal sectional areas
    [U_nom,K,dKda,dUda]=fea(a,F); 
    
    sgn1=sign( u1_constr-U_nom(1) );
    if sgn1==0 
        sgn1=1; end; % Signum function for U1
    sgn2=sign( u2_constr-U_nom(2) );
    if sgn2==0 
        sgn2=1; end; % Signum function for U2
    
    % Calculate eta of u1
    [etau1sqr, lagmult1]=getindexetau1;
    etau1= sgn1*etau1sqr;
    
    % Calculate eta of u2
    [etau2sqr, lagmult2]=getindexetau2;
    etau2= sgn2*etau2sqr;
    
    % Pick out smaller eta using aggregation function method
    eta=1/P * log( exp(-P*etau1) + exp(-P*etau2) ); % eta is negtive
    f0val= eta;
    
    etahistory(iter)=-eta;    % Aggregated value
    etahistoryu1(iter)=etau1; % Eta of U1
    etahistoryu2(iter)=etau2; % Eta of U2
    fprintf('change=%f   eta1=%.3f   eta2=%.3f\n',change,etau1,etau2);
   
    df0dx  =zeros(m,1); 
    df0dxu1=zeros(m,1);
    df0dxu2=zeros(m,1);
    
    % Sensitivity analysis by FDM
%     for ii=1:m   
%         a=A;
%         a_turb=0.001*a(ii,1); % Turbulence 
%         a(ii,1)=a(ii,1)+a_turb;
%         [detau1sqr]=getindexetau1;
%         detau1=sgn1*detau1sqr;
%         [detau2sqr]=getindexetau2;
%         detau2=sgn2*detau2sqr;
%         df0dxu1(ii,1)=( detau1-etau1 ) / a_turb;
%         df0dxu2(ii,1)=( detau2-etau2 ) / a_turb;
%         DetaDeta1= -( exp(-P*etau1) / ( exp(-P*etau1) + exp(-P*etau2) ) );
%         DetaDeta2= -1-DetaDeta1;
%         df0dx(ii,1)= DetaDeta1*df0dxu1(ii,1) + DetaDeta2*df0dxu2(ii,1);
%     end;
    
    % Analytical sensitivity analysis (Debug)
%     df0dxu1_a=-lagmult1.ineqnonlin*dUda(1,:)';
%     df0dxu2_a=-lagmult2.ineqnonlin*dUda(2,:)';
%     DetaDeta1_a= -( exp(-P*etau1) / ( exp(-P*etau1) + exp(-P*etau2) ) );
%     DetaDeta2_a= -1-DetaDeta1;
%     df0dx_a= DetaDeta1*df0dxu1 + DetaDeta2*df0dxu2;
    
    % Analytical sensitivity analysis
    df0dxu1=-lagmult1.ineqnonlin*dUda(1,:)';
    df0dxu2=-lagmult2.ineqnonlin*dUda(2,:)';
    DetaDeta1= -( exp(-P*etau1) / ( exp(-P*etau1) + exp(-P*etau2) ) );
    DetaDeta2= -1-DetaDeta1;
    df0dx= DetaDeta1*df0dxu1 + DetaDeta2*df0dxu2;
        

    % Constraint sensitivity analysis for upper problem 
    fval=rho*l*A-volume_constr;
    % Objective sensitivity analysis for upper problem
    dfdx(1,:)=rho*l;

    % Call GCMMA solver 
    [amma,ymma,zmma,lam,xsi,eta_gcmma,mu,zet,s,f0app,fapp] = ...
         gcmmasub(constr_num,design_num,iter,epsimin,a,amin,amax,low,upp, ...
                  raa0,raa,f0val,df0dx,fval,dfdx,a0,para_a,c,d);
    % Update design variables
    a=amma;
    A=a;
    % Calculate change from previous iteration step 
    change=max(abs( (etau1_old-etau1)/etau1_old ),...
               abs( (etau2_old-etau2)/etau2_old ));

    etau1_old=etau1;
    etau2_old=etau2;
  
end;

%%%% Print results %%%%
fprintf('\n\nSolution converged after %d iteration steps\n',iter);
fprintf('\nThe optimal eta equals %f\n',-eta);
fprintf('\nThe optimal sectional areas are\n');
disp(a); % Output optimal sectional areas

fval=rho*l*A-volume_constr;
fprintf('The value of constraint function equals \n');
disp(fval); % Output constraint value at optimum

toc % Timing end

% Plot iterarion history curves
xlab=1:1:iter;
plot(xlab,etahistory(1,1:iter),'ko-',...
     xlab,etahistoryu1(1,1:iter),'b*-',...
     xlab,etahistoryu2(1,1:iter),'rv-');
legend('Total Eta','Eta of U1','Eta of U2');
xlabel ('Number of iteration ');
ylabel ('Eta');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%% End of main function %%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%% READ DATA FROM FILE %%%
function [n,m,nc,nn,x,y,ihl,ihr,e,a,px,py,UncertLNum,W,var,rho]=...
    readfile()
fid=fopen('data_2d5_a.dat','r');
infor=fscanf(fid,'%d',3);
n=infor(1); % NUM OF NODES
m=infor(2); % NUM OF BARS, ALSO THE NUM OF DESIGN VARIABLES
nc=infor(3); % NUM OF FIXED NODES
nn=2*(n-nc); % NUM OF NON-FIXED NODES
x=fscanf(fid,'%f',n); % NODAL COORIDINATES IN X DIRECTION
y=fscanf(fid,'%f',n); % NODAL COORIDINATES IN Y DIRECTION
ihl=fscanf(fid,'%d',m); % NUM OF THE LEFT NODE OF EACH BAR
ihr=fscanf(fid,'%d',m); % NUM OF THE RIGHT NODE OF EACH BAR
e=fscanf(fid,'%f',m); % YOUNG'S MODULUS
a=fscanf(fid,'%f',m); % SECTIONAL AREA OF INITIAL DESIGN
px=fscanf(fid,'%f',n-nc); % Fx
py=fscanf(fid,'%f',n-nc); % Fy
UncertLNum=fscanf(fid,'%d',1);
W=zeros(UncertLNum);
for i=1:UncertLNum
    W(i,:)=fscanf(fid,'%f',UncertLNum);
end
var=fscanf(fid,'%f',1);
rho=fscanf(fid,'%f',1);
fclose(fid); % CLOSE FILE

%%%%%%%%%% Function q2p %%%%%%%%%%
function [p]=q2p(q,nomValue,var,MatTran)
% This function transforms normalized parameters back into true ones
p = nomValue+var * diag( MatTran * q ) * nomValue; 
%%%%%%%%% function end %%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [f,g]=objeta(q)
% Objective function of index eta for fmincon
f= sum( q.^2 ) ;
g= 2*q ;
%%%%%%%%% function end %%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [etau1sqr,lagmult]=getindexetau1
% Calculate reliability index 
global UncertLNum ;
q0 = zeros( UncertLNum,1 ); % Initial value
options=optimset('LargeScale','off','Algorithm','active-set',...
                 'Display','off','GradObj','on','GradConstr','on',...
                 'Derivative','off');
[q,fval,ExitFlag,Output,Lamda]=...
    fmincon(@objeta,q0,[],[],[],[],[],[],@conetau1,options);
etau1sqr=fval;
lagmult=Lamda;
%%%%%%%%% function end %%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [c,ceq,GC,GCeq]=conetau1(q)
% Constraint function of index eta for fmincon
global nn F var MatTran...
       u1_constr ...
       U_nom...
       UncertLNum; % These global variables are read-only
global a ; % These global variables are REWRITABLE ?
           % The sectional area at the current step serve as the 
           % nominal value
dUdaii=zeros(nn,UncertLNum);

sgn=sign( U_nom(1)-u1_constr );
if sgn==0 
    sgn=-1; 
end;

aPerturb=q2p(q, a, var, MatTran); % Transform normalized parameters back

[U, K]=fea(aPerturb,F);

% Sensitivity Analysis
for ii=1:UncertLNum
    [dKdaii]=GetdKda(ii);
    dUdaii(1:nn,ii)=-K\(dKdaii*U);
end

c= sgn*( (U(1)-u1_constr) );
GC=sgn*var*MatTran*diag(a)* dUdaii(1,:)';
ceq=[];
GCeq=[];
%%%%%%%%% function end %%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [etau2sqr,Lagmult]=getindexetau2
% Calculate reliability index
global UncertLNum ;
q0 = zeros( UncertLNum,1 ); % Initial value
options=optimset('LargeScale','off','Algorithm','active-set',...
                 'Display','off','GradObj','on','GradConstr','on',...
                 'Derivative','off');
[q,fval,ExitFlag,Output,Lamda]=...
    fmincon(@objeta,q0,[],[],[],[],[],[],@conetau2,options);
etau2sqr=fval;
Lagmult=Lamda;
%%%%%%%%% function end %%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [c,ceq,GC,GCeq]=conetau2(q)
% Constraint function of index eta for fmincon
global nn F var MatTran...
       u2_constr ...
       U_nom...
       UncertLNum; % These global variables are read-only
global a ; % These global variables are REWRITABLE ?
           % The sectional area at the current step serve as the 
           % nominal value
dUdaii=zeros(nn,UncertLNum);

sgn=sign( U_nom(2)-u2_constr );

if sgn==0 
    sgn=-1; 
end;

aPerturb=q2p(q, a, var, MatTran); % Transform normalized parameters back

[U, K]=fea(aPerturb,F);

% Sensitivity Analysis
for ii=1:UncertLNum
    [dKdaii]=GetdKda(ii);
    dUdaii(1:nn,ii)=-K\(dKdaii*U);
end

c= sgn*( (U(2)-u2_constr) );
GC=sgn*var*MatTran*diag(a)* dUdaii(2,:)';
ceq=[];
GCeq=[];
%%%%%%%%% function end %%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%     FEA SUBROUTINE    %%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [U,K,dKda,dUda,sigma,dg_idx_j]=fea(aPerturb,F)
global m nn i0 j0 l T...
       e ; % These global variables are read-only.

K=zeros(nn,nn);
dKda=zeros(nn,nn);
dKdxi=zeros(nn,nn);
dUda=zeros(nn,m);
dUdxi=zeros(nn,m);

% Solve FE problem for displacement
for loop=1:m
    [KE]=esm(e(loop),aPerturb(loop),l(loop),T(2*loop-1,1),T(2*loop,1));
    if i0(loop)<0
        K(j0(loop)+(1:2),j0(loop)+(1:2))=...
            K(j0(loop)+(1:2),j0(loop)+(1:2))+KE(1:2,1:2);
    else
        K(i0(loop)+(1:2),i0(loop)+(1:2))=...
            K(i0(loop)+(1:2),i0(loop)+(1:2))+KE(1:2,1:2);
        K(j0(loop)+(1:2),j0(loop)+(1:2))=...
            K(j0(loop)+(1:2),j0(loop)+(1:2))+KE(1:2,1:2);
        K(i0(loop)+(1:2),j0(loop)+(1:2))=...
            K(i0(loop)+(1:2),j0(loop)+(1:2))-KE(1:2,1:2);
        K(j0(loop)+(1:2),i0(loop)+(1:2))=...
            K(j0(loop)+(1:2),i0(loop)+(1:2))-KE(1:2,1:2);
    end
end 
U=K\F; 

% Output derivative of K with respect to a, namely dKda
if nargout>2
    for loop=1:m
        [dKE]=desm(e(loop),l(loop),T(2*loop-1,1),T(2*loop,1));
        
        if i0(loop)<0
            dKda(j0(loop)+(1:2),j0(loop)+(1:2))=...
                dKda(j0(loop)+(1:2),j0(loop)+(1:2))+dKE(1:2,1:2);
        else
            dKda(i0(loop)+(1:2),i0(loop)+(1:2))=...
                dKda(i0(loop)+(1:2),i0(loop)+(1:2))+dKE(1:2,1:2);
            dKda(j0(loop)+(1:2),j0(loop)+(1:2))=...
                dKda(j0(loop)+(1:2),j0(loop)+(1:2))+dKE(1:2,1:2);
            dKda(i0(loop)+(1:2),j0(loop)+(1:2))=...
                dKda(i0(loop)+(1:2),j0(loop)+(1:2))-dKE(1:2,1:2);
            dKda(j0(loop)+(1:2),i0(loop)+(1:2))=...
                dKda(j0(loop)+(1:2),i0(loop)+(1:2))-dKE(1:2,1:2);
        end
    end
end

if nargout>3
    for loop=1:m
        [dKE]=desm(e(loop),l(loop),T(2*loop-1,1),T(2*loop,1));
        
        if i0(loop)<0
            dKdxi(j0(loop)+(1:2),j0(loop)+(1:2))=...
                dKdxi(j0(loop)+(1:2),j0(loop)+(1:2))+dKE(1:2,1:2);
        else
            dKdxi(i0(loop)+(1:2),i0(loop)+(1:2))=...
                dKdxi(i0(loop)+(1:2),i0(loop)+(1:2))+dKE(1:2,1:2);
            dKdxi(j0(loop)+(1:2),j0(loop)+(1:2))=...
                dKdxi(j0(loop)+(1:2),j0(loop)+(1:2))+dKE(1:2,1:2);
            dKdxi(i0(loop)+(1:2),j0(loop)+(1:2))=...
                dKdxi(i0(loop)+(1:2),j0(loop)+(1:2))-dKE(1:2,1:2);
            dKdxi(j0(loop)+(1:2),i0(loop)+(1:2))=...
                dKdxi(j0(loop)+(1:2),i0(loop)+(1:2))-dKE(1:2,1:2);
        end
        dUda(:,loop)=-K\(dKdxi*U);
        dKdxi=zeros(nn,nn);
    end
end

% Output bar stress
if nargout>4 
    [sigma]=axf(m,e,l,T,U,i0,j0); % EVALUATE BAR STRESS
end

% Output 1st derivative of stress with respect to a
if nargout>5 
    dg_idx_j=zeros(m,m);
    for ii=1:m
        for jj=1:m
            dg_idx_j(ii,jj)=daxf(ii,e(ii),l(ii),T,dUdxi(:,jj),...
                i0(ii),j0(ii));
            dg_idx_j(m+ii,jj)=-dg_idx_j(ii,jj);
        end;
    end;
end

%%% Calculate derivative of K with respect to sectional area %%%
function [dKdaii]=GetdKda(ii)
global nn i0 j0 l T e ; % These global variables are read-only.
dKdaii=zeros(nn,nn);
[dKE]=desm(e(ii), l(ii), T(2*ii-1,1), T(2*ii,1));
if i0(ii)<0
    dKdaii(j0(ii)+(1:2),j0(ii)+(1:2))=...
        dKdaii(j0(ii)+(1:2),j0(ii)+(1:2))+dKE(1:2,1:2);
else
    dKdaii(i0(ii)+(1:2),i0(ii)+(1:2))=...
        dKdaii(i0(ii)+(1:2),i0(ii)+(1:2))+dKE(1:2,1:2);
    dKdaii(j0(ii)+(1:2),j0(ii)+(1:2))=...
        dKdaii(j0(ii)+(1:2),j0(ii)+(1:2))+dKE(1:2,1:2);
    dKdaii(i0(ii)+(1:2),j0(ii)+(1:2))=...
        dKdaii(i0(ii)+(1:2),j0(ii)+(1:2))-dKE(1:2,1:2);
    dKdaii(j0(ii)+(1:2),i0(ii)+(1:2))=...
        dKdaii(j0(ii)+(1:2),i0(ii)+(1:2))-dKE(1:2,1:2);
end
    
%%% DERIVATIVE OF ELEMENT STIFFNESS MATRIX WITH RESPECT TO Ei %%%
function [dKE]=desm(e,l,calpha,salpha)
TE=[calpha^2,calpha*salpha;calpha*salpha,salpha^2];
[dKE]=e/l*TE*TE';

%%% BAR STRESS %%%
function [sigma]=axf(m,eStress,l,T,U,i0,j0)
sigma=zeros(m,1);
for loop=1:m
    if i0(loop)<0 
        u1=[0 0]';
    else u1=[U(i0(loop)+1) U(i0(loop)+2)]';
    end;
    u2=[U(j0(loop)+1) U(j0(loop)+2)]';
    sigma(loop,1)=eStress(loop)/l(loop)*T((2*loop-1):(2*loop))'*(u2-u1);
end

%%% BAR STRESS SENSITIVITY %%%
function d_ij=daxf(i,eStressSen,l,T,dUdxi,i0,j0)
if i0<0 
    du1=[0 0]';
else du1=[dUdxi(i0+1) dUdxi(i0+2)]';
end;
du2=[dUdxi(j0+1) dUdxi(j0+2)]';
d_ij=eStressSen/l*T((2*i-1):(2*i))'*(du2-du1);

%%% Preprocess for FEA %%%
function [i0,j0,l,T]=preprocess(m,x,y,ihl,ihr,nc)
i0=2*(ihl-nc)-2;
j0=2*(ihr-nc)-2;
l=zeros(1,m);
for i=1:m
    lx=x(ihr(i))-x(ihl(i)); 
    ly=y(ihr(i))-y(ihl(i));  
    l(i)=sqrt(lx*lx+ly*ly);
    c_alpha=lx/l(i);s_alpha=ly/l(i);
    T((2*i-1):(2*i),1)=[c_alpha s_alpha];
end

%%% ELEMENT STIFFNESS MATRIX %%%
function [KE]=esm(e,a,l,calpha,salpha)
TE=[calpha^2,calpha*salpha;calpha*salpha,salpha^2];
[KE]=e*a/l*TE*TE';

