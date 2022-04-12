function Top_MFSE_EM(targetXY,dVElmIdx,dSElmIdx,nElX,nElY,refine,eps_r,lambda)
%% THE INPUT PARAMETER MEANS
% targetXY: the focus point location;
% dVElmIdx: the element indexes of the design domain;
% dSElmIdx: the element indexes of the substrate domain;
% nElX: the number of x-axis elements in the rectangular model domain;
% nElY: the number of y-axis elements in the rectangular model domain;
% refine : the ratio of adjacent observation point distance to the FE mesh size; 
% eps_r: the relative permittivity of the solid material;
% lambda: the incident wavelength.
%% SETUP PARAMETERS
[nely,nelx] = size(dVElmIdx); elelen = 1;
corlen_x = 12; corlen_y = 5;
ptdist = elelen*refine;
nptx = nelx/refine; npty = nely/refine;
tolne = nelx*nely; tolnpt = nptx*npty;
%% BUILD CORRELATION MATRIX
[Xpt, Ypt] = meshgrid((0.5:1:nptx)*ptdist, (npty-0.5:-1:0.5)*ptdist);
Xpt = Xpt(:); Ypt = Ypt(:);
corMat = zeros(tolnpt,tolnpt);
for i = 1:size(corMat,1)
    for j = i + 1:size(corMat,2)
        corMat(i,j) = exp(-(((Xpt(j)-Xpt(i))^2/corlen_x^2+(Ypt(j)-Ypt(i))^2/corlen_y^2)));
    end
end
corMat = corMat + corMat';
for i = 1:size(corMat, 1)
    corMat(i,i) = 1;
end
%% DO SERIES EXPANSION OF THE MATERIAL FIELD
if size(corMat,1)<1e5
    [eigfunMat, eigvalMat] = eig(corMat);
else
    [eigfunMat, eigvalMat] = eigs(corMat,1500);
end
eigvalVec = diag(eigvalMat);
[eigvalVec, eigsortind] = sort(eigvalVec, 'descend');
neig = 0; tmpsum = 0.;
while tmpsum < (1-1e-4)*sum(abs(eigvalVec))
    neig = neig + 1;
    tmpsum = tmpsum + eigvalVec(neig);
end
EXPANMat = sparse(1:neig, 1:neig, eigvalVec(1:neig).^(-1/2), neig, neig)...
    *eigfunMat(:,eigsortind(1:neig))';
clear eigfunMat;
%% COMPUTE PHI ON ELEMENTS AND MATERIAL-FIELD POINTS 
[Xe, Ye] = meshgrid((0.5:1:nelx)*elelen, (nely-0.5:-1:0.5)*elelen);
Xe = Xe(:); Ye = Ye(:);
eIntopMat = zeros(neig, tolne);
grsize = min(round(tolnpt/20), tolne);
ngr = ceil(tolne/grsize);
for igr = 1:ngr
    eind = (igr-1)*grsize+1:min(igr*grsize, tolne);
    Xe_sub = Xe(eind); Ye_sub = Ye(eind);
    eptvals = exp(-((((repmat(Xpt',length(eind),1)-repmat(Xe_sub, 1, tolnpt)).^2)/corlen_x^2 ...
        +((repmat(Ypt',length(eind),1)-repmat(Ye_sub, 1, tolnpt)).^2)/corlen_y^2)))';
    eptvals(abs(eptvals) < 1e-9) = 0;
    eIntopMat(:,eind) = EXPANMat*eptvals;
end
ptIntopMat = EXPANMat*corMat'; clear corMat;
%% SETUP OF PHYSICS PARAMETERS
phy.scale = 1e-9;
phy.eps_r = eps_r;
phy.k = 2*pi/(lambda*phy.scale);
%% SETUP OF ALL INDEX SETS, ELEMENT MATRICES AND RELATED QUANTITIES
dis.nElX = nElX;
dis.nElY = nElY;
dis.tElmIdx = targetXY(1)*(nElY+1) + (nElY + 1) - targetXY(2);
dis.dVElmIdx = dVElmIdx; dis.dSElmIdx = dSElmIdx;
[dis.LEM,dis.MEM] = ELEMENT_MATRICES(phy.scale);
[dis] = INDEX_SETS_SPARSE(dis);
%% INITIALIZE ITERATION
filThr.beta = 0.5; volfrac = 0.5; 
temp_y = -log(1/volfrac-1)/filThr.beta;
x = (temp_y)*ones(1,tolnpt)/ptIntopMat; x = x';
loop = 0; obj = 0;
change = 1.; ichange = 1; n = neig;
xmin = -1000*ones(n,1); xmax = 1000*ones(n,1);
low = xmin; upp = xmax;
xold1 = x;  xold2 = x; clf; history = [];
while (change>=0.005 || filThr.beta<20) && loop<150
    loop = loop + 1;
    objold = obj;
    [obj,dcdx] = FEM(x,dis,phy,filThr,eIntopMat);
    %% Update Design Variables
    m = 1; cc = 10000*ones(m,1);
    d = zeros(m,1); a0 = 1; a = zeros(m,1);
    fval = zeros(m, 1); fval(1) = 1; dfdx = zeros(m, n);
    [xmma,~,~,~,~,~,~,~,~,low,upp] = ...
        mmasub(m,n,loop,x,xmin,xmax,xold1,xold2, ...
        obj,dcdx,fval,dfdx,low,upp,a0,a,cc,d);
    xold2 = xold1; xold1 = x; x = xmma;
    %% TUNE PROJECTION PARAMETER 
    change = abs(obj-objold)/-obj;
    if change < 0.005 && loop > 30
        ichange = ichange+1;
    else
        ichange = 1;
    end
    if mod(ichange,3) == 0
        filThr.beta = min(filThr.beta * 1.5,20);
    end
    %% PRINT RESULTS
    fprintf([' It.:%5i Obj.:%9.4f numdesvars :%5i' ...
        ' beta:%5.1f ch.:%6.3f\n'],...
        loop,-obj,neig,filThr.beta,change);
    history = cat(2,history,-obj);
end
%% PRINT ITERATION HISTORY
figure(4)
plot(1:loop,history,'-k','LineWidth',2);xlabel('Iteration');ylabel('Objective');
set(gca, 'fontsize', 15, 'LineWidth', 1.5,'FontName','Times New Roman','FontWeight','bold');
end

function [LaplaceElementMatrix,MassElementMatrix] = ELEMENT_MATRICES(scaling) 
%% ELEMENT MATRICES
aa = scaling/2; bb = scaling/2;
k1 = (aa^2+bb^2)/(aa*bb); k2 = (aa^2-2*bb^2)/(aa*bb); k3 = (bb^2-2*aa^2)/(aa*bb);
LaplaceElementMatrix = [k1/3 k2/6 -k1/6 k3/6 ; k2/6 k1/3 k3/6 -k1/6; ...
    -k1/6  k3/6  k1/3  k2/6; k3/6 -k1/6  k2/6  k1/3];
MassElementMatrix = aa*bb*[4/9 2/9 1/9 2/9 ; 2/9 4/9 2/9 1/9 ; ...
    1/9 2/9 4/9 2/9; 2/9 1/9 2/9 4/9];
end

function dis = INDEX_SETS_SPARSE(dis)
%% CONNECTIVITY AND INDEX SETS
% INDEX SETS FOR SYSTEM MATRIX
nEX = dis.nElX; nEY = dis.nElY;
nodenrs = reshape(1:(1+nEX)*(1+nEY),1+nEY,1+nEX);
edofVec = reshape(nodenrs(1:end-1,1:end-1)+1,nEX*nEY,1);
dis.edofMat = repmat(edofVec,1,4)+repmat([0 nEY+[1 0] -1],nEX*nEY,1);
dis.iS = reshape(kron(dis.edofMat,ones(4,1))',16*nEX*nEY,1);
dis.jS = reshape(kron(dis.edofMat,ones(1,4))',16*nEX*nEY,1);
dis.idxDSdx = reshape(dis.edofMat',1,4*nEX*nEY);
% INDEX SETS FOR BOUNDARY CONDITIONS
TMP = repmat([1:nEY;2:nEY+1],2,1);
dis.iB1 = reshape(TMP,4*nEY,1);
dis.jB1 = reshape([TMP(2,:);TMP(1,:);TMP(3,:);TMP(4,:)],4*nEY,1);
TMP = repmat([1:(nEY+1):(nEY+1)*nEX;(nEY+1)+1:(nEY+1):(nEY+1)*nEX+1],2,1);
dis.iB2 = reshape(TMP,4*nEX,1);
dis.jB2 = reshape([TMP(2,:);TMP(1,:);TMP(3,:);TMP(4,:)],4*nEX,1);
TMP = repmat([(nEY+1)*(nEX)+1:(nEY+1)*(nEX+1)-1;(nEY+1)*(nEX)+2:(nEY+1)*(nEX+1)],2,1);
dis.iB3 = reshape(TMP,4*nEY,1);
dis.jB3 = reshape([TMP(2,:);TMP(1,:);TMP(3,:);TMP(4,:)],4*nEY,1);
TMP = repmat([2*(nEY+1):nEY+1:(nEY+1)*(nEX+1);(nEY+1):nEY+1:(nEY+1)*(nEX)],2,1);
dis.iB4 = reshape(TMP,4*nEX,1);
dis.jB4 = reshape([TMP(2,:);TMP(1,:);TMP(3,:);TMP(4,:)],4*nEX,1);
dis.iRHS = TMP;
% INDEX SETS FOR INTEGRATION OF ALL ELEMENTS
ima0 = repmat([1,2,3,4,1,2,3,4,1,2,3,4,1,2,3,4],1,nEX*nEY)';
jma0 = repmat([1,1,1,1,2,2,2,2,3,3,3,3,4,4,4,4],1,nEX*nEY)';
addTMP = repmat(4*(0:nEX*nEY-1),16,1);
addTMP = addTMP(:);
dis.iElFull = ima0 + addTMP;
dis.jElFull = jma0 + addTMP;
% INDEX SETS FOR SENSITIVITY COMPUTATIONS
dis.iElSens = (1:4*nEX*nEY)';
jElSens = repmat(1:nEX*nEY,4,1);
dis.jElSens = jElSens(:);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% The code is supplementary to the corresponding paper:                    %
% Anisotropic material-field series expansion for topological design       %
% of optical metalens,                                                     %
% Zhaoyou Sun, Pai Liu, Yangjun Luo                                        %
% Submitted to Optics express, 2022                                        %
%                                                                          %
% This code is based on                                                    %
% Compact 200 line MATLAB code for inverse design in photonics by topology %
% optimization: tutorial, by Christiansen et al., J OPT SOC AM B (2021)    %                                             %
%                                                                          %
% Disclaimer:                                                              %
% The authors reserve all rights but do not guaranty that the code is      %
% free from errors. Furthermore, we shall not be liable in any event       %
% caused by the use of the program.                                        %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
