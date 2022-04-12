function [FOM,sensFOM] = FEM(x,dis,phy,filThr,eIntopMat)
% This function is modified on the basis of the code presented in the following paper:
% Christiansen, R. E., & Sigmund, O. (2021). Compact 200 line MATLAB code for inverse design in photonics by topology optimization: tutorial. 
% Journal of the Optical Society of America - B, 38(2), 510-520. https://doi.org/10.1364/JOSAB.405955
% The original code in the above paper can be download from: http://www.topopt.dtu.dk.
%% OBJECTIVE FUNCTION AND GRADIENT EVALUATION
% DISTRIBUTE MATERIAL IN MODEL DOMAIN BASED ON DESIGN FIELD
dFP(1:dis.nElY,1:dis.nElX) = 0;
dFP(dis.dSElmIdx(:)) = 1; 
dFP(dis.dVElmIdx(:)) = x'*eIntopMat;

% THE MATERIAL FIELD PROJECTION AND MATERIAL INTERPOLATION
dFPST = dFP; dFPST(dis.dVElmIdx(:)) = THRESHOLD(dFP(dis.dVElmIdx(:)), filThr.beta);
DdFSTDFS = DERIVATIVE_OF_THRESHOLD(dFP, filThr.beta);
[A,dAdx] = MATERIAL_INTERPOLATION(phy.eps_r,dFPST,1.0);

% CONSTRUCT THE SYSTEM MATRIX
[dis,F] = BOUNDARY_CONDITIONS_RHS(phy.k,dis,phy.scale);
dis.vS = reshape(dis.LEM(:)-phy.k^2*dis.MEM(:)*(A(:).'),16*dis.nElX*dis.nElY,1);
S = sparse([dis.iS(:);dis.iBC(:)],[dis.jS(:);dis.jBC(:)],[dis.vS(:);dis.vBC(:)]);

% SOLVING THE STATE SYSTEM: S * Ez = F
[L,U,Q1,Q2] = lu(S);
Ez = Q2 * (U\(L\(Q1 * F)));  Ez = full(Ez);

% FIGURE OF MERIT
P = sparse(dis.tElmIdx,dis.tElmIdx,1,...
    (dis.nElX+1)*(dis.nElY+1),(dis.nElX+1)*(dis.nElY+1));
FOM = Ez' * P * Ez;

% ADJOINT RIGHT HAND SIDE
AdjRHS = P*(2*real(Ez) - 1i*2*imag(Ez));

% SOLVING THE ADJOING SYSTEM: S.' * AdjLambda = AdjRHS
AdjLambda = (Q1.') * ((L.')\((U.')\((Q2.') * (-1/2*AdjRHS)))); % Solving

% COMPUTING SENSITIVITIES
dis.vDS = reshape(-phy.k^2*dis.MEM(:)*(dAdx(:).'),16*dis.nElX*dis.nElY,1);
DSdx = sparse(dis.iElFull(:),dis.jElFull(:),dis.vDS(:));
DSdxMulV = DSdx * Ez(dis.idxDSdx); 
DsdxMulV = sparse(dis.iElSens,dis.jElSens,DSdxMulV);
sens = 2*real(AdjLambda(dis.idxDSdx).' * DsdxMulV);
sens = full(reshape(sens,dis.nElY,dis.nElX));
sensFOM = eIntopMat * (reshape(sens(dis.dVElmIdx),[],1).*reshape(DdFSTDFS(dis.dVElmIdx),[],1));
FOM = -FOM; sensFOM = -sensFOM(:);

% PLOTTING AND PRINTING
figure(1);
imagesc((reshape(Ez.*conj(Ez),dis.nElY+1,dis.nElX+1))); colormap('jet');colorbar; axis equal;axis off;axis tight;
figure(2);
imagesc([reshape(1-dFPST(dis.dVElmIdx(:)),[],dis.nElX);reshape(1-dFPST(dis.dSElmIdx(:)),[],dis.nElX)]); colormap(gray); axis equal; caxis([0,1]);axis equal;axis off;axis tight;
figure(3);
ePhi = reshape(x'*eIntopMat,[],dis.nElX);
contourf(flip(ePhi,1),[0,0]); colormap([0 0 0; 0 0 1;1 0 0; 0 1 0; 1 1 1]); axis equal;axis off;axis tight;drawnow;
end

function [A,dAdx] = MATERIAL_INTERPOLATION(eps_r,x,alpha_i)
%% MATERIAL PARAMETER INTERPOLATION
A = 1 + x*(eps_r-1) - 1i * alpha_i * x .* (1 - x); % Interpolation
dAdx = (eps_r-1)*(1+0*x) - 1i * alpha_i * (1 - 2*x); % Derivative of interpolation
end

function [dis,F] = BOUNDARY_CONDITIONS_RHS(waveVector,dis,scaling)
%% SER THE BOUNDARY CONDITIONS
AbsBCMatEdgeValues = 1i*waveVector*scaling*[1/6 ; 1/6 ; 1/3 ; 1/3];
% ALL BOUNDARIES HAVE ABSORBING BOUNDARY CONDITIONS
dis.iBC = [dis.iB1(:);dis.iB2(:);dis.iB3(:);dis.iB4(:)];
dis.jBC = [dis.jB1(:);dis.jB2(:);dis.jB3(:);dis.jB4(:)];
dis.vBC = repmat(AbsBCMatEdgeValues,2*(dis.nElX+dis.nElY),1);
% BOTTOM BOUNDARY HAS INCIDENT PLANE WAVE
F = zeros((dis.nElX+1)*(dis.nElY+1),1);
F(dis.iRHS(1,:)) = F(dis.iRHS(1,:))+1i*waveVector;
F(dis.iRHS(2,:)) = F(dis.iRHS(2,:))+1i*waveVector;
F = scaling*F;
end

function ePhiProj = THRESHOLD(ePhi,beta)
%% PROJECTION
ePhiProj = 1./(1+exp(-beta*ePhi));
end

function edproj = DERIVATIVE_OF_THRESHOLD(ePhi,beta)
%% DERIVATIVE OF PROJECTION
edproj = exp(-beta*ePhi)./((1+exp(-beta*ePhi)).^2);
end