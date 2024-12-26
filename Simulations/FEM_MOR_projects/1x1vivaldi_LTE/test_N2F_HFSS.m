% test near fields to far fields transformation of fields extracted from
% HFSS field calculator tool. HFSS project = "1x1vivaldi_LTE"
clc; close all;
addpath('..\_matlabSolvers');

[z0, k0, lambda0] = getFreeSpaceElectricalParams(7.5e9);
[boxPos, boxN, dS] = buildBox([1 1 1 1 1 1],...
  -.01, .015, -.01, .035, -.003, .087, ...
  floor(.025*20/lambda0), ... % xPts
  floor(.045*20/lambda0), ... % yPts
  floor(.090*20/lambda0), ... % zPts
  .999, 1, 0);

%%--- extracting points for HFSS exportation
dlmwrite('Points.pts', boxPos.', 'delimiter', ' ', 'precision', 17);
disp('Waiting for eField.reg and hField.reg...')
pause;

%%--- load extracted fields
A= dlmread('eField.reg',' ',1,0);
E(1,:) = A(:,5) + 1i* A(:,6);
E(2,:) = A(:,7) + 1i* A(:,8);
E(3,:) = A(:,9) + 1i* A(:,10);
clear A;
A= dlmread('hField.reg',' ',1,0);
H(1,:) = A(:,5) + 1i* A(:,6);
H(2,:) = A(:,7) + 1i* A(:,8);
H(3,:) = A(:,9) + 1i* A(:,10);
clear A;

%%--- near fields to far electric field transformation
sizetref = 300;
sizet = 900;
nbrCoeffs = 31;
[phiFF, thetaFF] = meshgrid(deg2rad([0 90]), ...
  linspace(0, 2*pi*(sizetref-1)/sizetref, sizetref));
[Opt,Opp] =  vf_n2fOpFieldsFFT(k0, z0, boxPos, boxN, dS,...
  thetaFF, phiFF, nbrCoeffs);
figure; surf(abs(Opt(:,:,1)),'EdgeColor','none'); axis('tight');view(90,0);
[ifftOp, phiFFT, thetaFFT] = getLookAnglePoly(phiFF(1,:), sizet, ...
  sizetref, nbrCoeffs);
nf = [H(:); E(:)];

EtFF = zeros(sizet, size(phiFF,2));
EpFF = zeros(sizet, size(phiFF,2));
for i=1:size(phiFF,2)
  EtFF(:,i) = ifftOp*(Opt(:,:,i)*nf);
  EpFF(:,i) = ifftOp*(Opp(:,:,i)*nf);
end

Pr = computeRadiatedPower(nf, boxN, dS);
[gaint, gainp] = vf_computeGain(z0, EtFF, EpFF, Pr); % directivity

vf_plotFFPolarCutPlanes(figure(), 25.8, 10, thetaFFT, gaint, gainp, 'N2F', ...
    '\bf Directivity pattern');