% A non-probabilistic reliable based truss optimization(NRBTO) MATLAB code
% Written by Song BAI, 2010.09.16
% This program aims at finding the maximun reliable index under the
% constraint of total weight using coherency function method
% This code was designed for 25-bars truss
%
% Input: u_c--> Vertical displacement constraint on the right bottom node
%        v_c--> Total weight constraint

%%%%%%% Main Function %%%%%%%%%
function rdo25(agg_P, u1_c , u2_c, v_c)
% Test input: rdo(10,3.5,-0.28,18)
global m nn ;
global a i0 j0 T l id ;
global F ; % Only some components are uncertain, the deterministic ones remain unchanged
global Z var MatTran; % Z is nominal uncertain parameter vector, var is the variation range
global U_nom ;
global u14_constr u15_constr volume_constr;
global sgn1 sgn2;
tic % Timing start

P = agg_P ; % Penalty factor of aggregate function
u14_constr = u1_c;
u15_constr = u2_c ;
volume_constr = v_c;
 
% READ DATA FROM FILE
[n,m,nc,nn,x,y,z,ihl,ihr,...
    pxNom, pyNom, pzNom, Z, a, W, var] = readfile(); 

[Q,Lam] = eig(W);
LamSqtInv = diag(1./diag(sqrt(Lam)));
MatTran = Q*LamSqtInv; % Transform Matrix "T"

F = zeros(nn,1); % Only some components of F are uncertain parameters
F(1:3:3*(n-nc),1) = pxNom(1:n-nc);
F(2:3:3*(n-nc),1) = pyNom(1:n-nc);
F(3:3:3*(n-nc),1) = pzNom(1:n-nc);

% The matrix id is used for distinguishing E1 and E2 group
id = [1,1,1,1,1,1,1,1,1,1,1,1,1,0,0,0,0,0,0,0,0,0,0,0,0;...
      0,0,0,0,0,0,0,0,0,0,0,0,0,1,1,1,1,1,1,1,1,1,1,1,1]; 

% ALLOCATE GLOBAL LOCATION FOR EACH BAR
[i0,j0] = I0J0(ihl,ihr,nc);

% EVALUATE BAR LENGTH AND TRANSFORMATION MATRIX T
[l,T] = bar_len(m,x,y,z,ihl,ihr);

%%% Parameters for GCMMA %%%
rho = 0.001; % Scaling factor for objective function
constr_num = 1;  % Constraint num = 1 (Volume constraint)
design_num = m;  % Design variable num
epsimin = 0.00001;
amin = 0.01*ones(design_num,1);
amax = 40*ones(design_num,1);
low = amin;
upp = amax;
raa0 = 0.01;
raa = 0.01*ones(constr_num,1);
aGCMMA = zeros(constr_num,1);
a0GCMMA = 1;
cGCMMA = 100*ones(constr_num,1); 
dGCMMA=ones(constr_num,1);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%

change = 1.0;
etau14_old = 1.0; etau15_old = 1.0;
iter = 0; itermax = 100;
etahistory    = zeros(1,itermax);
etahistoryu14 = zeros(1,itermax);
etahistoryu15 = zeros(1,itermax);
LowerIter1 = 0;
LowerIter2 = 0;

while change > 0.001 && iter < itermax
    iter = iter + 1;
    
    % Solution for displacement under nominal-valued uncertain parameters
    [U_nom] = fea3D(a, Z);
    
    sgn1 = sign( u14_constr - U_nom(14) );
    if sgn1 == 0
        sgn1 = 1;
    end
    
    sgn2 = sign( - (u15_constr - U_nom(15)) );
    if sgn2 == 0
        sgn2 = 1;
    end
    
    % Calculate eta of u14
    % The subroutine has the same structure as that considers load uncertainty only
    [etau14sqr, q_opt_u14, Lambda_u14, iterStep14] = getindexeta1(); % "q_opt_u14" is the MPF point
    etau14 = sgn1 * etau14sqr;
      
    % Calculate eta of u17
    % The subroutine has the same structure as that considers load uncertainty only
    [etau15sqr, q_opt_u15, Lambda_u15, iterStep15] = getindexeta2(); % "q_opt_u15" is the MPF point
    etau15 = sgn2 * etau15sqr;
    
    % Calculate objective function (eta with minimum value) using coherency function
    eta = 1 / P * log( exp(-P * etau14) + exp(-P * etau15) ); % Notice that eta is negtive
    f0val = eta; % Objective function variable for GCMMA
    
    etahistory(iter) = -eta;
    etahistoryu14(iter) = etau14;
    etahistoryu15(iter) = etau15; % Record objective function value at each iteration step
    
    fprintf('Iter. %d   change=%f   eta14=%.3f   eta15=%.3f\n',...
             iter, change, etau14, etau15); % Print resluts to screen
    fprintf('Lower iter. steps:  etau14=%d    etau15=%d\n', iterStep14,iterStep15);
    LowerIter1 = LowerIter1 + iterStep14;
    LowerIter2 = LowerIter2 + iterStep15;
    
    Popt14 = q2p(q_opt_u14, Z, var, MatTran);
    [Utemp, dUdE1temp, dUdE2temp, dUdFtemp, dUdx] = fea3D(a, Popt14);
    df0dxu14 =  - Lambda_u14.ineqnonlin * dUdx(14,:)';
    % Sensitivity of etaU14 with respect to desvars
    
    Popt15 = q2p(q_opt_u15, Z, var, MatTran);
    [Utemp, dUdE1temp, dUdE2temp, dUdFtemp, dUdx] = fea3D(a, Popt15);
    df0dxu15 =  Lambda_u15.ineqnonlin * dUdx(15,:)';
    % Sensitivity of etaU15 with respect to desvars
    
    DetaDeta14 = -( exp(-P*etau14) / ( exp(-P*etau14) + exp(-P*etau15) ) );
    DetaDeta15 = -1 - DetaDeta14;
    df0dx = DetaDeta14 * df0dxu14 + DetaDeta15 * df0dxu15;
    
%     % Sensitivity of eta with respect to desvars obtained by aggregate function
%         A = a;
%         for ii=1:m   % Sensitivity analysis by FDM
%             a = A;
%             aTurb = 0.01 * a(ii,1); % A turbulence
%             a(ii,1) = a(ii,1) + aTurb;
%             [detau14sqrFDM] = getindexeta1();
%             detau14FDM = sgn1 * detau14sqrFDM;
%             [detau15sqrFDM] = getindexeta2();
%             detau15FDM = sgn2 * detau15sqrFDM;
%             df0dxu14FDM(ii,1) = ( detau14FDM - etau14 ) / aTurb;
%             df0dxu15FDM(ii,1) = ( detau15FDM - etau15 ) / aTurb;
%             DetaDeta14 = -( exp(-P*etau14) / ( exp(-P*etau14) + exp(-P*etau15) ) );
%             DetaDeta15 = -1 - DetaDeta14;
%             df0dx(ii,1) = DetaDeta14 * df0dxu14FDM(ii,1) + DetaDeta15 * df0dxu15FDM(ii,1);
%         end
%         a = A;
    
    % Calculate constraint function (volume) value
    fval = rho * l * a - volume_constr;
    
    % Calculate sensitivities of constraint function 
    dfdx(1,:) = rho * l;
    
    % Call GCMMA solver
    [amma,ymma,zmma,lam,xsi,eta_gcmma,mu,zet,s,f0app,fapp] = ...
        gcmmasub(constr_num, design_num, iter, ...
                 epsimin, a, amin, amax, low, upp, raa0, raa, ...
                 f0val, df0dx, fval, dfdx, ...
                 a0GCMMA, aGCMMA, cGCMMA, dGCMMA); %#ok<NASGU>
   
    a = amma; % Update design variables
 
    % Calculate objective increasement 
    change = max(abs( (etau14_old - etau14) / etau14_old ),...
                 abs( (etau15_old - etau15) / etau15_old ));
    
    etau14_old = etau14;
    etau15_old = etau15;
    
%     disp([df0dxu14,df0dxu14FDM,df0dxu15,df0dxu15FDM]);
    
end

%%%% Print results %%%%
fprintf('\n\nSolution converged after %d iteration steps\n',iter);
fprintf('\nThe optimal eta equals %f\n',-eta);
fprintf('\nThe optimal sectional areas are\n');
for ii=1:m
    fprintf('Bar # %d --> %.4f\n',ii,a(ii));
end

fval = rho * l * a - volume_constr;
fprintf('The value of constraint function equals \n');
disp(fval);
toc % Timing ends
fprintf('Average lower iteration steps number: eta14=%f    eta15=%f\n',LowerIter1/iter,LowerIter2/iter);

% Plot iteration history curve
xlab=1:1:iter;
plot(xlab,etahistory(1,1:iter),'ko-',...
    xlab,etahistoryu14(1,1:iter),'b*-',...
    xlab,etahistoryu15(1,1:iter),'rv-');
legend('Total Eta','Eta of U14','Eta of U15');
xlabel ('Iteration step');
ylabel ('Eta');
%%%%%%%%%%%%%% End of main function %%%%%%%%%%%%%%%

%%% READ DATA FROM FILE %%%
function [n,m,nc,nn,x,y,z,ihl,ihr, pxNom, pyNom, pzNom, Z, a, W, var]=readfile()
fid=fopen('data_3d_ea.dat','r');

infor=fscanf(fid,'%d',3);
n=infor(1); % NUM OF NODES
m=infor(2); % NUM OF BARS
nc=infor(3); % NUM OF FIXED NODES
nn=3*(n-nc); % NUM OF NON-FIXED NODES

UncertainPara = fscanf(fid, '%d', 1);
UncertainE = fscanf(fid, '%d', 1);

x = fscanf(fid,'%f',n); % X NODAL COORIDINATES 
y = fscanf(fid,'%f',n); % Y NODAL COORIDINATES
z = fscanf(fid,'%f',n); % Z NODAL COORIDINATES

ihl = fscanf(fid,'%d',m); % NUM OF THE LEFT NODE OF EACH BAR
ihr = fscanf(fid,'%d',m); % NUM OF THE RIGHT NODE OF EACH BAR

e = fscanf(fid,'%f', UncertainE); % YOUNG'S MODULUS
a = fscanf(fid,'%f', m); % SECTIONAL AREA OF INITIAL DESIGN

pxNom = fscanf(fid,'%f',n-nc); % Fx
pyNom = fscanf(fid,'%f',n-nc); % Fy
pzNom = fscanf(fid,'%f',n-nc); % Fz

Z = [e(1); e(2); pyNom(5); pzNom(5)];

W = zeros(UncertainPara);
for i = 1 : UncertainPara
    W(i,:) = fscanf(fid, '%f', UncertainPara);
end

var = fscanf(fid,'%f',1);

fclose(fid); % CLOSE FILE

%%% FEA SUBROUTINE FOR SPACIAl TRUSS %%%
function [U, dUdE1, dUdE2, dUdF, dUdx] = fea3D(a, Z)
% Input: 
% a: sectional areas
% Z: uncertain parameter vector
% Output:
% U: displacement
% dUdEi: sensitivity of displacement with respect to Young's modulus
% dUdF: sensitivity of displacement with respect to force
% sigma: bar stress
% dg_idx_j: sensitivity of bar stress with respect to sectional areas
global m nn i0 j0 l T id ;
global F ; 

K = zeros(nn, nn); 
dKdx = zeros(nn, nn);
dKdE1 = zeros(nn, nn);
dKdE2 = zeros(nn, nn);
dUdx = zeros(nn, m);

E = zeros(m,1);
E(1:13,1)  = Z(1);
E(14:25,1) = Z(2); % Uncertain Young's modulus

Ffea = F;
Ffea(14,1) = Z(3); 
Ffea(15,1) = Z(4); % Reassign value to uncertain force components

for loop = 1:m % Calculate K
    [KE] = Ke( E(loop), a(loop), l(loop), ...
               T(3*loop-2,1), T(3*loop-1,1), T(3*loop,1) );
    if i0(loop) < 0
        K(j0(loop)+(1:3),j0(loop)+(1:3))=...
            K(j0(loop)+(1:3),j0(loop)+(1:3))+KE(1:3,1:3);
    else
        K(i0(loop)+(1:3),i0(loop)+(1:3))=...
            K(i0(loop)+(1:3),i0(loop)+(1:3))+KE(1:3,1:3);
        K(j0(loop)+(1:3),j0(loop)+(1:3))=...
            K(j0(loop)+(1:3),j0(loop)+(1:3))+KE(1:3,1:3);
        K(i0(loop)+(1:3),j0(loop)+(1:3))=...
            K(i0(loop)+(1:3),j0(loop)+(1:3))-KE(1:3,1:3);
        K(j0(loop)+(1:3),i0(loop)+(1:3))=...
            K(j0(loop)+(1:3),i0(loop)+(1:3))-KE(1:3,1:3);
    end
end

U = K \ Ffea; % Calculate DISPLACEMENT

if nargout > 1 
    for loop = 1:m   
        [dKE] = dKedE(a(loop), l(loop), T(3*loop-2,1), T(3*loop-1,1), T(3*loop,1));
        
        % Calculate dKedE1
        if id(1,loop) > 0
            if i0(loop) < 0
                dKdE1(j0(loop)+(1:3),j0(loop)+(1:3)) = ...
                    dKdE1(j0(loop)+(1:3),j0(loop)+(1:3))+dKE(1:3,1:3);
            else
                dKdE1(i0(loop)+(1:3),i0(loop)+(1:3)) = ...
                    dKdE1(i0(loop)+(1:3),i0(loop)+(1:3))+dKE(1:3,1:3);
                dKdE1(j0(loop)+(1:3),j0(loop)+(1:3)) = ...
                    dKdE1(j0(loop)+(1:3),j0(loop)+(1:3))+dKE(1:3,1:3);
                dKdE1(i0(loop)+(1:3),j0(loop)+(1:3)) = ...
                    dKdE1(i0(loop)+(1:3),j0(loop)+(1:3))-dKE(1:3,1:3);
                dKdE1(j0(loop)+(1:3),i0(loop)+(1:3)) = ...
                    dKdE1(j0(loop)+(1:3),i0(loop)+(1:3))-dKE(1:3,1:3);
            end
        end
        
        % Calculate dKedE2
        if id(2,loop) > 0
            if i0(loop) < 0
                dKdE2(j0(loop)+(1:3),j0(loop)+(1:3)) = ...
                    dKdE2(j0(loop)+(1:3),j0(loop)+(1:3))+dKE(1:3,1:3);
            else
                dKdE2(i0(loop)+(1:3),i0(loop)+(1:3)) = ...
                    dKdE2(i0(loop)+(1:3),i0(loop)+(1:3))+dKE(1:3,1:3);
                dKdE2(j0(loop)+(1:3),j0(loop)+(1:3)) = ...
                    dKdE2(j0(loop)+(1:3),j0(loop)+(1:3))+dKE(1:3,1:3);
                dKdE2(i0(loop)+(1:3),j0(loop)+(1:3)) = ...
                    dKdE2(i0(loop)+(1:3),j0(loop)+(1:3))-dKE(1:3,1:3);
                dKdE2(j0(loop)+(1:3),i0(loop)+(1:3)) = ...
                    dKdE2(j0(loop)+(1:3),i0(loop)+(1:3))-dKE(1:3,1:3);
            end
        end
        
    end
end

dUdE1 = - K \ (dKdE1 * U);
dUdE2 = - K \ (dKdE2 * U);

I = zeros(nn, 1);
I(14,1) = 1.0;
invKC14 = K \ I;

I = zeros(nn, 1);
I(15,1) = 1.0;
invKC15 = K \ I;

dUdF = [invKC14(14,1), invKC15(14,1); ...
        invKC14(15,1), invKC15(15,1)];
    
if nargout>3 % dUdx
    for loop = 1:m
        [dKE] = dKedx(E(loop), l(loop), T(3*loop-2,1), T(3*loop-1,1), T(3*loop,1));
        if i0(loop)<0
            dKdx(j0(loop)+(1:3),j0(loop)+(1:3))=...
                dKdx(j0(loop)+(1:3),j0(loop)+(1:3))+dKE(1:3,1:3);
        else
            dKdx(i0(loop)+(1:3),i0(loop)+(1:3))=...
                dKdx(i0(loop)+(1:3),i0(loop)+(1:3))+dKE(1:3,1:3);
            dKdx(j0(loop)+(1:3),j0(loop)+(1:3))=...
                dKdx(j0(loop)+(1:3),j0(loop)+(1:3))+dKE(1:3,1:3);
            dKdx(i0(loop)+(1:3),j0(loop)+(1:3))=...
                dKdx(i0(loop)+(1:3),j0(loop)+(1:3))-dKE(1:3,1:3);
            dKdx(j0(loop)+(1:3),i0(loop)+(1:3))=...
                dKdx(j0(loop)+(1:3),i0(loop)+(1:3))-dKE(1:3,1:3);
        end
        
        dUdx(:,loop) = - K \ (dKdx * U);
        dKdx = zeros(nn, nn);
    end
end

%%% GLOBAL LOCATION FOR EACH BAR %%%
function [i0,j0] = I0J0(ihl,ihr,nc)
i0 = 3*(ihl-nc)-3;
j0 = 3*(ihr-nc)-3;

%%% BAR LENGTH AND MATRIX T %%%
function [l,T] = bar_len(m,x,y,z,ihl,ihr)
l=zeros(1,m);
for i=1:m
    lx=x(ihr(i))-x(ihl(i));
    ly=y(ihr(i))-y(ihl(i));
    lz=z(ihr(i))-z(ihl(i));
    l(i)=sqrt(lx*lx+ly*ly+lz*lz);
    c_alpha=lx/l(i); c_beta=ly/l(i); c_gamma=lz/l(i);
    T((3*i-2):(3*i),1)=[c_alpha c_beta c_gamma];
end

%%% ELEMENT STIFFNESS MATRIX %%%
function [KE] = Ke( e, a, l, c_alpha, c_beta, c_gamma )
TE = [c_alpha c_beta c_gamma]';
[KE] = e * a / l * TE * TE';

%%% Sensitivity of Ke with respect to E %%%
function [dKedE] = dKedE( a, l, c_alpha, c_beta, c_gamma)
TE = [c_alpha, c_beta, c_gamma]';
[dKedE] = a / l * TE * TE';

%%% Sensitivity of Ke with respect to E %%%
function [dKedx] = dKedx( e, l, c_alpha, c_beta, c_gamma)
TE = [c_alpha, c_beta, c_gamma]';
[dKedx] = e / l * TE * TE';

%%%%%%%%%% Function q2p %%%%%%%%%%%
function [p] = q2p(q, Z, var, MatTran)
% This function transforms normalized loads back into true loads
p = Z + var * diag( MatTran * q ) * Z; 
%%%%%%%%% function end %%%%%%%%%%%

%%%%%%% Subroutine getindexeta1 %%%%%%%%
% Calculate index eta #2 
function [eta1sqr, q, Lambda, iterStep] = getindexeta1()
q0 = [ 0, 0, 0, 0 ]'; % Initial value
options = optimset('LargeScale','off','Algorithm','active-set','Display','off',...
                   'GradObj','on','GradConstr','on','Derivative','off');
[q, fval, ExitFlag, Output, Lambda] = ...
    fmincon(@objeta, q0, [],[],[],[],[],[], @coneta1, options);
eta1sqr = fval;
iterStep = Output.iterations;
%%%%%%%%% function end %%%%%%%%%%%

%%%%%%%% Subroutine getindexeta2 %%%%%%%%
function [eta2sqr, q, Lambda, iterStep] = getindexeta2()
% Calculate index eta during iteration process
q0 = [ 0, 0, 0, 0 ]'; % Initial value
options = optimset('LargeScale','off','Algorithm','active-set','Display','off',...
                   'GradObj','on','GradConstr','on','Derivative','off');
[q, fval, ExitFlag, Output, Lambda]=...
    fmincon(@objeta, q0, [],[],[],[],[],[], @coneta2, options);
eta2sqr = fval;
iterStep = Output.iterations;
%%%%%%%%% function end %%%%%%%%%%%

%%%%%%% Objective subroutine for GCMMA %%%%%%%%%
% Objective function of index eta for fmincon
function [f,g]=objeta(q)
f = sum(q.^2);
g = 2*q';
%%%%%%%%% function end %%%%%%%%%%%

%%%% Lower-level Constraint subroutine #1 for fmincon %%%%
function [c, ceq, GC, GCeq] = coneta1(q)
global a;
global Z var MatTran;
global u14_constr ; % Concerned displacement component constraint #14
global sgn1;

p = q2p(q, Z, var, MatTran);

[U, dUdE1, dUdE2, dUdF] = fea3D(a, p);

c =  - sgn1 * (U(14) - u14_constr);
% GC = - sgn1 * var * MatTran * diag(Z) * ...
%      [ dUdE1(14); dUdE2(14); dUdF(1,1); dUdF(2,1)];
GC = - sgn1 * var * ([ dUdE1(14), dUdE2(14), dUdF(1,1), dUdF(2,1)] * diag(Z) * MatTran)';
 
ceq = [];
GCeq = [];
%%%%%%%%% function end %%%%%%%%%%%

%%%% Constraint subroutine #2 for fmincon %%%%
function [c, ceq, GC, GCeq] = coneta2(q)
global a ;
global Z var MatTran;
global u15_constr ; % Concerned displacement component constraint #15
global sgn2;

p = q2p(q, Z, var, MatTran);

[U, dUdE1, dUdE2, dUdF] = fea3D(a, p);

c = - sgn2 * ( - (U(15) - u15_constr));
% GC = - sgn2 * var * MatTran * diag(Z) * ...
%      [ dUdE1(15); dUdE2(15); dUdF(1,2); dUdF(2,2)];
GC = - sgn2 * var * (-[ dUdE1(15), dUdE2(15), dUdF(1,2), dUdF(2,2)] * diag(Z) * MatTran)';

ceq = [];
GCeq = [];
%%%%%%%%% function end %%%%%%%%%%%
