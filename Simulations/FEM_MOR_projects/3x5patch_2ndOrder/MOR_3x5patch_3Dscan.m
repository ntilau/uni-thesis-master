%% MOR on 15 patches for 3D scanning
clc; close all; clear all;
tStart = tic;
addpath(genpath('..\_matlabSolvers'));
addpath(genpath('..\_matlabSolvers\PardisoInterface'));

flag = 2; % 1: solve all 2: load and build output matrices 3: load System
nbrPorts = 15;
nbrCoeffs = 21;
nbrSmpls = 200; % need to oversample the operator
freq = 2.35e9;
phi = linspace(0, pi, 61);
if flag == 4
  tic;
  fprintf('#> Loading saved System matrices ... ');
  load SystemMatrices;
  fprintf('%2.4g s\n', toc);
end
%% get system matrices
if flag == 1 || flag == 2
  solTime = tic();
  [z0, k0, lambda0] = getFreeSpaceElectricalParams(freq);
  [boxPos, boxN, boxdS] = ...
    buildBox([1 0 1 1 1 1], -.14, .14, -.2, .2, 0, .07, ...
    22, 32, 6, .999, 0, 0);
  [fullA, fullB, fullX] = getInputSystemMatrices(flag, nbrPorts, true);
  [fullCt, fullCp, fullF, sizetref, nbrCoeffs] = ...
    getOutputSystemMatrices(flag, ...
    nbrCoeffs, nbrSmpls, k0, z0, boxPos, boxN, boxdS, phi);
  solvingTime = toc(solTime);
  save SystemMatrices fullA fullB fullCt fullCp fullF fullX ...
    sizetref nbrCoeffs z0 k0 boxPos boxN boxdS solvingTime;
end
%% Array topology
equiRowX = [0 1 2 0 1 2 0 1 2 0 1 2 0 1 2];
equiRowY = [0 0 0 1 1 1 2 2 2 3 3 3 4 4 4];
% the FEM eigenmode solver can arbitrarily switch the phase of the
% excitations of 180° -> use of +singleEnded to recover the phases
phase = [-1 -1 1 -1 -1 -1 1 -1 1 1 -1 -1 -1 1 -1]; % 3x5patch_2ndOrder
%% MOR space selection
[spanT, spanP] = getSpiralingHelicoidalTrajectory(15, 0, false);
plotSelectedAngles(spanT, spanP);

spanningV = zeros(size(fullA,1),size(spanT,2));
for j = 1:length(spanT)
  for k=1:nbrPorts
    spanningV(:,j) = spanningV(:,j) + fullX(:,k) .* phase(k) * ...
      exp(-1i*pi* ( (equiRowX(k)-1).* sin(deg2rad(spanT(j))) .* ...
      cos(deg2rad(spanP(j))) +  (equiRowY(k)-1) .* ...
      sin(deg2rad(spanT(j))) .* sin(deg2rad(spanP(j))) ) );
  end
end
%% Model order reduction
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
%% Testing ROM
testTime = tic();
% pattern shaping test
weightingX = ones(1,3); %cos(2*pi*(-1:1)/6);
weightingY = ones(1,5); %cos(2*pi*(-2:2)/10);

clear frames;
handle = figure();

nbrScan = 100;
[testT, testP] = getSpiralingHelicoidalTrajectory(nbrScan, 1, false, 0.5);

romTime = zeros(size(testT));
fullTime = romTime;
error = romTime;

rect = [100 100 560 420];
set(handle,'Position',rect);
set(handle,'Color',[1 1 1]);

romX = romA\romB;

for j=1:size(testT,2)
  
  angleTime = tic;
    
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

  handle = vf_plotFF3d(handle, thetaFF, phiFF, 40, gaint, gainp, ...
    ['\theta_s = ', sprintf('%2.3g', testT(j)), ...
    '°, \phi_s = ', sprintf('%2.3g', testP(j)), '°']);
  frames(j) = getframe(handle);

  error(j) = getL2error(EtFF+EpFF, refEtFF+refEpFF);
  
  fprintf('#> Scan angle time = %2.4g\n', toc(angleTime));
end
close(handle);
fprintf('#> Testing time = %2.4g s.\n',toc(testTime));
%% build movie
aviobj = avifile('scan.avi','Compression','none');
aviobj = addframe(aviobj,frames);
for i=1:length(frames)
  aviobj = addframe(aviobj,frames(length(frames)-i+1));
end
aviobj = close(aviobj);

fprintf('#> Total computation time = %2.4g s.\n',toc(tStart));