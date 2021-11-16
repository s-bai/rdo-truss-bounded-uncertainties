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
% Modified by Song, BAI, 2011.4.23
% This version is meant to verify the analytical sensitivity


function tenbar105(agg_P, u2_c , u8_c, v_c)
global m nn i0 j0 l T...
       F e_nom var MatTran...
       u2_constr u8_constr volume_constr...
       id...
       U_nom sgn2 sgn8; % Fixed global variables
global a ; % REWRITABLE global variables
tic

u2_constr=u2_c; 
u8_constr=u8_c ; 
volume_constr=v_c; 
P=agg_P ; % Penalty factor for aggregate function

% READ DATA FROM INPUT FILE
[n,m,nc,nn,x,y,ihl,ihr,e,a,px,py,e_nom,W,var,rho]=readfile();

% The matrix id is used for distinguishing E1 and E2 group
id=[1,1,1,1,1,1,0,0,0,0;...
    0,0,0,0,0,0,1,1,1,1]; 

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

% Applied loads
F=zeros(nn,1); 
F(1:2:2*(n-nc),1)=px(1:n-nc);
F(2:2:2*(n-nc),1)=py(1:n-nc);
% ALLOCATE GLOBAL COORDINATES FOR EACH BAR
[i0,j0]=I0J0(ihl,ihr,nc);
% EVALUATE BAR LENGTH AND TRANSFORMATION MATRIX T
[l,T]=bar_len(m,x,y,ihl,ihr);

% Transform Matrix
[Q,Lam]=eig(W);
LamSqtInv=diag(1./diag(sqrt(Lam)));
MatTran=Q*LamSqtInv;

change=1.0; etau2_old=1.0; etau8_old=1.0; 
iter=0; itermax=50;
etahistory=zeros(1,itermax);
etahistoryu2=zeros(1,itermax);
etahistoryu8=zeros(1,itermax);
A=a;

% Start iteration
while change>0.001 && iter< itermax
    iter=iter + 1;
    
    % Calculate disp. under nominal stiffness matrix
    [U_nom,K,dKdE1,dKdE2,dUdxi]=fea(e,a,F,id); 
    % "e" is the nominal value of Young's modulus
    
    sgn2=sign( u2_constr-U_nom(2) );
    if sgn2==0 
        sgn2=1; 
    end % Signum function for U2
    
    sgn8=sign( u8_constr-U_nom(8) );
    if sgn8==0 
        sgn8=1; 
    end % Signum function for U8
    
    % Calculate eta of u2
    [etau22,lagmult1,q_opt_u2]=getindexetau2;
    etau2= sgn2*etau22;
    
    % Calculate eta of u8
    [etau82,lagmult2,q_opt_u8]=getindexetau8;
    etau8= sgn8*etau82;
    
    % Pick out smaller eta using aggregation function method
    eta=1/P * log( exp(-P*etau2) + exp(-P*etau8) ); % eta is negtive
    f0val= 0.1*eta;
    
    etahistory(iter)  =sign(-eta)*sqrt(abs(-eta));    % Aggregated value
    etahistoryu2(iter)=sign(etau2)*sqrt(abs(etau2)); % Eta of U2
    etahistoryu8(iter)=sign(etau8)*sqrt(abs(etau8)); % Eta of U8
    fprintf('change=%f   et2=%.3f   et8=%.3f\n',change,etau2,etau8);
   
    df0dx=zeros(m,1); 
    df0dxu2=zeros(m,1);
    df0dxu8=zeros(m,1);
    
%     % Sensitivity analysis by FDM
%     for ii=1:m   
%         a=A;
%         a_turb=0.001*a(ii,1); % Turbulence 
%         a(ii,1)=a(ii,1)+a_turb;
%         [detau22]=getindexetau2;
%         detau2=sgn2*detau22;
%         [detau82]=getindexetau8;
%         detau8=sgn8*detau82;
%         df0dxu2_FDM(ii,1)=( detau2-etau2 ) / a_turb;
%         df0dxu8_FDM(ii,1)=( detau8-etau8 ) / a_turb;
% %         DetaDeta2= -( exp(-P*etau2) / ( exp(-P*etau2) + exp(-P*etau8) ) );
% %         DetaDeta8= -1-DetaDeta2;
% %         df0dx(ii,1)= DetaDeta2*df0dxu2(ii,1) + DetaDeta8*df0dxu8(ii,1);
%     end;
%     a=A;
    
    % Analytical sensitivity analysis
    eGroup=q2p(q_opt_u2,e_nom,var,MatTran);
    eTemp(1,1:6) =eGroup(1);
    eTemp(1,7:10)=eGroup(2);
    [U_temp,K_temp,dKdE1_temp,dKdE2_temp,dUdxi_q_opt]=fea(eTemp,a,F,id);
    df0dxu2=-lagmult1.ineqnonlin*dUdxi_q_opt(2,:)';
    
    eGroup=q2p(q_opt_u8,e_nom,var,MatTran);
    eTemp(1,1:6) =eGroup(1);
    eTemp(1,7:10)=eGroup(2);
    [U_temp,K_temp,dKdE1_temp,dKdE2_temp,dUdxi_q_opt]=fea(eTemp,a,F,id);
    df0dxu8=-lagmult2.ineqnonlin*dUdxi_q_opt(8,:)';
    
    DetaDeta2= -( exp(-P*etau2) / ( exp(-P*etau2) + exp(-P*etau8) ) );
    DetaDeta8= -1-DetaDeta2;
    df0dx= (DetaDeta2*df0dxu2 + DetaDeta8*df0dxu8)*0.1;

    % Constraint sensitivity analysis for upper problem 
    fval=rho*l*A-volume_constr;
    % Objective sensitivity analysis for upper problem
    dfdx(1,:)=rho*l;
    
%     disp([sign(df0dxu2(1,1)*dUdxi(1,1)),sign(df0dxu2(1,1)*dUdxi(1,1))]);

    aLastStep=a;
    % Call GCMMA solver 
    [amma,ymma,zmma,lam,xsi,eta_gcmma,mu,zet,s,f0app,fapp] = ...
         gcmmasub(constr_num,design_num,iter,epsimin,a,amin,amax,low,upp, ...
                  raa0,raa,f0val,df0dx,fval,dfdx,a0,para_a,c,d);
    % Update design variables
    a=amma;
    A=a;
    
    % Calculate change from previous iteration step 
    change=max(abs( (etau2_old-etau2)/etau2_old ),...
        abs( (etau8_old-etau8)/etau8_old ));

    etau2_old=etau2;
    etau8_old=etau8;
  
end;

%%%% Print results %%%%
fprintf('\n\nSolution converged after %d iteration steps\n',iter);
fprintf('\nThe optimal eta equals %f\n',-eta);
fprintf('\nThe optimal sectional areas are\n');
disp(aLastStep);
fval=rho*l*aLastStep-volume_constr;
fprintf('The value of constraint function equals \n');
disp(fval);
toc
xlab=1:1:iter;
plot(xlab,etahistory(1,1:iter),'ko-',...
     xlab,etahistoryu2(1,1:iter),'b*-',...
     xlab,etahistoryu8(1,1:iter),'rv-');
legend('Total Eta','Eta of U2','Eta of U8');
xlabel ('Number of iteration ');
ylabel ('Eta');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%% End of main function %%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%% READ DATA FROM FILE %%%
function [n,m,nc,nn,x,y,ihl,ihr,e,a,px,py,e_nom,W,var,rho]=readfile()
fid=fopen('data_2d10_e.dat','r');
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
e_nom=fscanf(fid,'%f',UncertLNum);
W=zeros(UncertLNum);
for i=1:UncertLNum
    W(i,:)=fscanf(fid,'%f',UncertLNum);
end
var=fscanf(fid,'%f',1);
rho=fscanf(fid,'%f',1);
fclose(fid); % CLOSE FILE

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%    Lower Problem Sub    %%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%% Function q2p %%%%%%%%%%%
function [p]=q2p(q,e_nom,var,MatTran)
% This function transforms normalized loads back into true loads
p = e_nom+var * diag( MatTran * q ) * e_nom; 
%%%%%%%%% function end %%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [f,g]=objeta(q)
% Objective function of index eta for fmincon
f= (q(1)^2+q(2)^2) ;
g= [ 2*q(1), 2*q(2)]' ;
%%%%%%%%% function end %%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [eta1,lagmult1,q]=getindexetau2
% Calculate 2nd index eta during iteration process
q0=[ 0., 0. ]'; % Initial value
options=optimset('LargeScale','off','Algorithm','active-set','Display','off',...
    'GradObj','on','GradConstr','on','Derivative','off');
[q,fval,ExitFlag,Output,Lamda]=...
    fmincon(@objeta,q0,[],[],[],[],[],[],@conetau2,options);
eta1=fval;
lagmult1=Lamda;
%%%%%%%%% function end %%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [c,ceq,GC,GCeq]=conetau2(q)
% Constraint function of index eta for fmincon
global e_nom var MatTran...
       sgn2 u2_constr...
       F id; % Fixed global variables
       
   
global a ; % REWRITABLE global variables

% sgn=sign( U_nom(2)-u2_constr );
% if sgn==0 
%     sgn=-1; 
% end;

eGroup=q2p(q, e_nom, var, MatTran); % Transform normalized E back

eTemp(1,1:6) =eGroup(1);
eTemp(1,7:10)=eGroup(2);

[U, K, dKdE1, dKdE2]=fea(eTemp,a,F,id);

% Sensitivity Analysis
dUdE1=-K\(dKdE1*U);
dUdE2=-K\(dKdE2*U);

c= -sgn2*( (U(2)-u2_constr) );
GC=-sgn2 * var * ([ dUdE1(2), dUdE2(2) ] * diag(e_nom) * MatTran)';
ceq=[];
GCeq=[];
%%%%%%%%% function end %%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [eta2,lagmult2,q]=getindexetau8
% Calculate index eta during iteration process
q0=[ 0., 0. ]'; % Initial value
options=optimset('LargeScale','off','Algorithm','active-set','Display','off',...
    'GradObj','on','GradConstr','on','Derivative','off');
[q,fval,ExitFlag,Output,Lamda]=...
    fmincon(@objeta,q0,[],[],[],[],[],[],@conetau8,options);
eta2=fval;
lagmult2=Lamda;
%%%%%%%%% function end %%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [c,ceq,GC,GCeq]=conetau8(q)
% Constraint function of index eta for fmincon
global e_nom var MatTran...
       sgn8 u8_constr...
       F id; % Fixed global variables
       
global a ; % REWRITABLE global variables

% sgn=sign( U_nom(8)-u8_constr );
% if sgn==0 
%     sgn=-1; 
% end;

eGroup=q2p(q, e_nom, var, MatTran); % Transform normalized E back

eTemp(1,1:6) =eGroup(1);
eTemp(1,7:10)=eGroup(2);

[U, K, dKdE1, dKdE2]=fea(eTemp,a,F,id);

% Sensitivity Analysis
dUdE1=-K\(dKdE1*U);
dUdE2=-K\(dKdE2*U);

c= -sgn8*( U(8)-u8_constr );
GC=-sgn8 * var * ([ dUdE1(8), dUdE2(8) ] * diag(e_nom) * MatTran)';
ceq=[];
GCeq=[];
%%%%%%%%% function end %%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%     FEA SUBROUTINE    %%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [U,K,dKdE1,dKdE2,dUdxi,sigma,dg_idx_j]=fea(eFea,a,F,id)
global m nn i0 j0 l T; % Fixed global variables

K=zeros(nn,nn);
dKdE1=zeros(nn,nn);
dKdE2=zeros(nn,nn);
dUdxi=zeros(nn,m);
dKdxi=zeros(nn,nn);

for loop=1:m
    [KE]=esm(eFea(loop),a(loop),l(loop),T(2*loop-1,1),T(2*loop,1));
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
U=K\F; % Solve

if nargout>2 % Output dKdE
    for loop=1:m
        [dKE]=desm(a(loop),l(loop),T(2*loop-1,1),T(2*loop,1));
        if id(1,loop)>0
            if i0(loop)<0
                dKdE1(j0(loop)+(1:2),j0(loop)+(1:2))=...
                    dKdE1(j0(loop)+(1:2),j0(loop)+(1:2))+dKE(1:2,1:2);
            else
                dKdE1(i0(loop)+(1:2),i0(loop)+(1:2))=...
                    dKdE1(i0(loop)+(1:2),i0(loop)+(1:2))+dKE(1:2,1:2);
                dKdE1(j0(loop)+(1:2),j0(loop)+(1:2))=...
                    dKdE1(j0(loop)+(1:2),j0(loop)+(1:2))+dKE(1:2,1:2);
                dKdE1(i0(loop)+(1:2),j0(loop)+(1:2))=...
                    dKdE1(i0(loop)+(1:2),j0(loop)+(1:2))-dKE(1:2,1:2);
                dKdE1(j0(loop)+(1:2),i0(loop)+(1:2))=...
                    dKdE1(j0(loop)+(1:2),i0(loop)+(1:2))-dKE(1:2,1:2);
            end
        end
    end
    
    for loop=1:m
        [dKE]=desm(a(loop),l(loop),T(2*loop-1,1),T(2*loop,1));
        if id(2,loop)>0
            if i0(loop)<0
                dKdE2(j0(loop)+(1:2),j0(loop)+(1:2))=...
                    dKdE2(j0(loop)+(1:2),j0(loop)+(1:2))+dKE(1:2,1:2);
            else
                dKdE2(i0(loop)+(1:2),i0(loop)+(1:2))=...
                    dKdE2(i0(loop)+(1:2),i0(loop)+(1:2))+dKE(1:2,1:2);
                dKdE2(j0(loop)+(1:2),j0(loop)+(1:2))=...
                    dKdE2(j0(loop)+(1:2),j0(loop)+(1:2))+dKE(1:2,1:2);
                dKdE2(i0(loop)+(1:2),j0(loop)+(1:2))=...
                    dKdE2(i0(loop)+(1:2),j0(loop)+(1:2))-dKE(1:2,1:2);
                dKdE2(j0(loop)+(1:2),i0(loop)+(1:2))=...
                    dKdE2(j0(loop)+(1:2),i0(loop)+(1:2))-dKE(1:2,1:2);
            end
        end
    end
end

if nargout>4 % Output dUdxi
    for loop=1:m
        [dKE]=desm_da(eFea(loop),l(loop),T(2*loop-1,1),T(2*loop,1));
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
        dUdxi(:,loop)=-K\(dKdxi*U);
        dKdxi=zeros(nn,nn);
    end
end


if nargout>5 % Output bar stress
    [sigma]=axf(m,eFea,l,T,U,i0,j0); % EVALUATE BAR STRESS
end

if nargout>6 % Output 1st derivative of stress with respect to a
    dg_idx_j=zeros(m,m);
    for ii=1:m
        for jj=1:m
            dg_idx_j(ii,jj)=daxf(ii,eFea(ii),l(ii),T,dUdxi(:,jj),...
                i0(ii),j0(ii));
            dg_idx_j(m+ii,jj)=-dg_idx_j(ii,jj);
        end
    end
end
    
%%% Sensitivity of K WITH RESPECT TO Ei %%%
function [dKE]=desm(a,l,calpha,salpha)
TE = [calpha^2, calpha*salpha; ...
      calpha*salpha, salpha^2];
[dKE] = a / l * TE * TE';
%%% End of subroutine %%%

function [dKE]=desm_da(e,l,calpha,salpha)
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

%%% GLOBAL COORDINATES FOR EACH BAR %%%
function [i0,j0]=I0J0(ihl,ihr,nc)
i0=2*(ihl-nc)-2;
j0=2*(ihr-nc)-2;

%%% BAR LENGTH AND MATRIX T %%%
function [l,T]=bar_len(m,x,y,ihl,ihr)
l=zeros(1,m);
for i=1:m
    lx=x(ihr(i))-x(ihl(i)); 
    ly=y(ihr(i))-y(ihl(i));  
    l(i)=sqrt(lx*lx+ly*ly);
    c_alpha=lx/l(i);s_alpha=ly/l(i);
    T((2*i-1):(2*i),1)=[c_alpha s_alpha];
end

%%% ELEMENT STIFFNESS MATRIX %%%
function [KE]=esm(eStifMat,a,l,calpha,salpha)
TE=[calpha^2,calpha*salpha;calpha*salpha,salpha^2];
[KE]=eStifMat*a/l*TE*TE';