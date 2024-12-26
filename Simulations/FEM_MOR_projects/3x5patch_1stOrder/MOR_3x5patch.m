%% MOR on 15 patches for planar scanning
clc; close all; %clear all;
addpath(genpath('..\_matlabSolvers'));
addpath(genpath('..\_matlabSolvers\PardisoInterface'));

flag = 2; % 1: solve all 2: load and build output matrices
nbrPorts = 15;
nbrCoeffs = 21;
nbrSmpls = 200;
freq = 2.35e9;
phi = [0 pi/2];
%% get system matrices
if flag == 1 || flag == 2
  [z0, k0, lambda0] = getFreeSpaceElectricalParams(freq);
  [boxPos, boxN, boxdS] = ...
    buildBox([1 0 1 1 1 1], -.14, .14, -.2, .2, 0, .07, ...
    20, 30, 5, .999, 1, 0); % scale the box to ensure that the sampling 
                            % points are part of the FEM domain
  [fullA, fullB, fullX] = getInputSystemMatrices(flag, nbrPorts);
  [fullCt, fullCp, fullF, sizetref, nbrCoeffs] = ...
    getOutputSystemMatrices(flag, ...
    nbrCoeffs, nbrSmpls, k0, z0, boxPos, boxN, boxdS, phi);
end
%% Array topology
equiRowX = [0 1 2 0 1 2 0 1 2 0 1 2 0 1 2];
equiRowY = [0 0 0 1 1 1 2 2 2 3 3 3 4 4 4];
% the FEM eigenmode solver can arbitrarily switch the phase of the
% excitations of 180° -> use of +singleEnded to recover the phases
phase = [-1 -1 -1 1 -1 1 1 1 -1 1 -1 -1 -1 -1 1];
%% MOR space selection 
phiPlane = 90; % [°]
nbrSpanAng = 5; % number of spanning angles
expSelCoeff = 1.5;
spanT = 90/(expSelCoeff^(nbrSpanAng-1)-1) * ...
  (expSelCoeff.^((1:nbrSpanAng)-1)-1);
spanP = phiPlane*ones(size(spanT)); % phi=90° plane

spanningV = zeros(size(fullA,1),size(spanT,2));
for j = 1:length(spanT)
  for k=1:nbrPorts
    spanningV(:,j) = spanningV(:,j) + fullX(:,k) .* phase(k) * ...
      exp(-1i*pi* ( (equiRowX(k)-1).* sin(deg2rad(spanT(j))) .* ...
      cos(deg2rad(spanP(j))) +  (equiRowY(k)-1) .* ...
      sin(deg2rad(spanT(j))) .* sin(deg2rad(spanP(j))) ) );
  end
end
%% model order reduction
if false % check orthonormalization timings
  tic; [V, R] = qr(spanningV,0); a = toc;
  tic; [V, s, v] = svd(spanningV,0); b = toc;
  fprintf('QR time = %2.4g s, SVD time = %2.4g s\n', a, b);
end
[V, R] = qr(spanningV,0);
romA = V'*fullA*V;
romB = V'*fullB;
romCt = zeros(nbrCoeffs,size(V,2),size(phi,2));
romCp = romCt;
for i=1:size(phi,2);
  romCt(:,:,i) = fullCt{i}*V;
  romCp(:,:,i) = fullCp{i}*V; 
end
romF = fullF*V;
sizet = 5000;
[ifftOp, phiFF, thetaFF] = getLookAnglePoly(phi, sizet, sizetref, ...
  nbrCoeffs);
%% testing for chosen angle
testTime = tic;
handle = figure();
testT = linspace(0,45,10);
testP = phiPlane*ones(1,10);
% plotSelectedAngles(spanT, spanP, testT, testP);
romTime = zeros(size(testT));
fullTime = romTime;
error = romTime;

romX = romA\romB;
for j=1:size(testT,2)
  
  %-----
  tic;
  romXscan = zeros(size(V,2),1);
  for k=1:nbrPorts
    romXscan = romXscan + romX(:,k) .* phase(k) * ...
      exp(-1i*pi* ( (equiRowX(k)-1).* sin(deg2rad(testT(j))) .* ...
      cos(deg2rad(testP(j))) +  (equiRowY(k)-1) .* ...
      sin(deg2rad(testT(j))) .* sin(deg2rad(testP(j))) ) );
  end
  
  EtFF=zeros(size(thetaFF,1),size(phiFF,2));
  EpFF=zeros(size(thetaFF,1),size(phiFF,2));    
  for i=1:size(thetaFF,2)
    EtFF(:,i) = ifftOp*(romCt(:,:,i)*romXscan);
    EpFF(:,i) = ifftOp*(romCp(:,:,i)*romXscan);
  end

  romTime(j) = toc;
  %-----
  tic;

  fullXscan = zeros(size(fullA,1),1);
  for k=1:nbrPorts
    fullXscan = fullXscan + fullX(:,k) .* phase(k) * ...
      exp(-1i*pi* ( (equiRowX(k)-1).* sin(deg2rad(testT(j))) .* ...
      cos(deg2rad(testP(j))) +  (equiRowY(k)-1) .* ...
      sin(deg2rad(testT(j))) .* sin(deg2rad(testP(j))) ) );
  end

  refEtFF=zeros(size(thetaFF,1),size(phiFF,2));
  refEpFF=zeros(size(thetaFF,1),size(phiFF,2));

  for i=1:size(thetaFF,2)
    refEtFF(:,i) = ifftOp*(fullCt{i}*fullXscan);
    refEpFF(:,i) = ifftOp*(fullCp{i}*fullXscan);
  end

  fullTime(j) = toc;
  %-----
  
  Pr = computeRadiatedPower(romF * romXscan, boxN, boxdS);
  PrRef = computeRadiatedPower(fullF * fullXscan, boxN, boxdS);
  fprintf('##> Pr relative difference = %2.4g\n', (PrRef-Pr)/PrRef);
  
  [gaint, gainp] = vf_computeGain(z0, EtFF, EpFF, Pr);
  [refGaint, refGainp] = vf_computeGain(z0, refEtFF, refEpFF, PrRef);

  vf_plotFFPolarCutPlanes(handle, 56.5, 20, thetaFF, gaint, gainp, 'ROM', ...
    thetaFF(1:40:end,:), refGaint(1:40:end,:), refGainp(1:40:end,:), ...
    'Full model', ['\bf \theta_s = ', sprintf('%2.3g°', testT(j)), ...
    ', \phi_s = ', sprintf('%2.3g°', testP(j))]);

  error(j) = getL2error(EtFF+EpFF, refEtFF+refEpFF);
   
end

fprintf('#> Testing time : %2.4g s\n', toc(testTime));