%%% check the number of spanning angles required in MOR for several
%%% steering operations
clc; close all; %clear all;
addpath(genpath('..\_matlabSolvers'));
addpath(genpath('..\_matlabSolvers\PardisoInterface'));

flag = 3; % 1: solve all 2: load and build output matrices
nbrPorts = 15;
nbrCoeffs = 21;
nbrSmpls = 200;
freq = 2.35e9;
phi = [0 pi/2];
% options for testing
scan3d = true; % false for phi constant planes scanning (2D scan)
N=19; % max number of vectors
phiPlane = 90; % [°] 2D scan
%% get system matrices
if flag == 1 || flag == 2
  [z0, k0, lambda0] = getFreeSpaceElectricalParams(freq);
  [boxPos, boxN, boxdS] = ...
    buildBox([1 0 1 1 1 1], -.14, .14, -.2, .2, 0, .07, ...
    20, 30, 5, .999, 1, 0);
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
spanError = zeros(1,N);
for q=1:N % number of vectors
  fprintf('#> Number of spanning vectors = %d\n', q);
  if scan3d
    [spanT, spanP, testT, testP] = ...
      getSpiralingHelicoidalTrajectory(q, 2*N, 0);
  else
    spanT = linspace(0, 90, q);
    spanP = phiPlane*ones(1, q);
  end
  spanningV = zeros(size(fullA,1),size(spanT,2));
  for j = 1:length(spanT)
    for k=1:nbrPorts
      spanningV(:,j) = spanningV(:,j) + fullX(:,k) .* phase(k) * ...
        exp(-1i*pi* ( (equiRowX(k)-1).* sin(deg2rad(spanT(j))) .* ...
        cos(deg2rad(spanP(j))) +  (equiRowY(k)-1) .* ...
        sin(deg2rad(spanT(j))) .* sin(deg2rad(spanP(j))) ) );
    end
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
  if ~scan3d
    if length(spanT)==1
      testT = 0;
    else
      testT = (spanT(2)-spanT(1))/2;
    end
    testP = phiPlane;
  end
  

  error = zeros(size(testT));
  romX = romA\romB;
  for j=1:size(testT,2)

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

    error(j) = getL2error(EtFF+EpFF, refEtFF+refEpFF);

  end
  
  spanError(q) = mean(error);
  
end

figProp = getFigureProperties();
figure; semilogy(1:length(spanError), spanError, '-r*', ...
  'LineWidth', figProp.lw, ...
  'MarkerSize',figProp.ms);
xlabel('Scan angle space dimension', 'Fontsize',figProp.fs);
ylabel('Relative error', 'Fontsize',figProp.fs);
% printEPS('','errorYZ')