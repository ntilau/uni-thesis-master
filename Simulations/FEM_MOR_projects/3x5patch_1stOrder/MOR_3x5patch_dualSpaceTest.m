%% test of MOR on 15 patches
clc; close all;% clear all;
addpath(genpath('..\_matlabSolvers'));
addpath(genpath('..\_matlabSolvers\PardisoInterface'));

flag = 2; % 1: solve all 2: load and build output matrices 
          % flag <= 3: computes the dual space vectors
nbrPorts = 15;
nbrCoeffs = 21;
nbrSmpls = 200;
freq = 2.35e9;
phi = [0 pi/2]; % constant phi planes
N = 3; % number of spanning vectors for model order reduction
%% get system matrices
if flag == 1 || flag == 2
  [z0, k0, lambda0] = getFreeSpaceElectricalParams(freq);
  [boxPos, boxN, boxdS] = ...
    buildBox([1 0 1 1 1 1], -.14, .14, -.2, .2, 0, .07, ...
    20, 30, 5, .999, 0, 0);
  [fullA, fullB, fullX] = getInputSystemMatrices(flag, nbrPorts);
  [fullCt, fullCp, fullF, sizetref, nbrCoeffs] = ...
    getOutputSystemMatrices(flag, ...
    nbrCoeffs, nbrSmpls, k0, z0, boxPos, boxN, boxdS, phi);
end
%% Array topology
equiRowX = [0 1 2 0 1 2 0 1 2 0 1 2 0 1 2];
equiRowY = [0 0 0 1 1 1 2 2 2 3 3 3 4 4 4];
phase = [-1 -1 -1 1 -1 1 1 1 -1 1 -1 -1 -1 -1 1];
%% MOR space selection
spanT = linspace(0,90,N);
spanP = 90*ones(1,N);
% plotSelectedAngles(spanT, spanP)

spanningV = zeros(size(fullA,1),size(spanT,2));
for j = 1:length(spanT)
  for k=1:nbrPorts
    spanningV(:,j) = spanningV(:,j) + fullX(:,k) .* phase(k) * ...
      exp(-1i*pi* ( (equiRowX(k)-1).* sin(deg2rad(spanT(j))) .* ...
      cos(deg2rad(spanP(j))) +  (equiRowY(k)-1) .* ...
      sin(deg2rad(spanT(j))) .* sin(deg2rad(spanP(j))) ) );
  end
end
%%
sizet = 1000;
[ifftOp, phiFF, thetaFF] = getLookAnglePoly(phi, sizet, sizetref, ...
  nbrCoeffs);

%% compute dual space solutions
if flag <= 3
  % output coefficients for N look angles
  idx = floor(linspace(sizet/(N+1),sizet*(N)/(N+1),N)); % select equally 
                                                        % spaced angles in 
                                                        % the range of
                                                        % sizet to ensure
                                                        % matching in those
                                                        % look angles
  dualSpanT = (thetaFF(idx,1)).'; % [rad]
  outputCoeffs = ifftOp(idx,:);
  % compute solution vectors
  spanningW = callParDiSo(fullA', (outputCoeffs*fullCt{2})');
              % must choose a plane for the matching test 2: phi(2)=90°
  
end
%% model order reduction
[V, R] = qr(spanningV,0);
[W, R] = qr(spanningW,0);
romA = W'*fullA*V;
romB = W'*fullB;
romCt = zeros(nbrCoeffs,size(V,2),size(phi,2));
romCp = romCt;
for i=1:size(phi,2);
  romCt(:,:,i) = fullCt{i}*V;
  romCp(:,:,i) = fullCp{i}*V; 
end
romF = fullF*V;
%% Testing ROM
testTime = tic;
testT = linspace(0,90,(N-1)*4+1);
testP = 45*ones(1,size(testT,2)); % notice we have matching even if the
                                  % scanning plane is different from the
                                  % dual space's one (test phi plane: 
                                  % phi = 45°)

% testT = 10;
% testP = 90;

error = zeros(1,size(testT,2));
errT = error;
errP = error;
errTxz = error;
errPxz = error;
errTyz = error;
errPyz = error;

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

    
  Pr = computeRadiatedPower(romF * romXscan, boxN, boxdS);
  PrRef = computeRadiatedPower(fullF * fullXscan, boxN, boxdS);
  fprintf('##> Pr relative difference = %2.4g\n', (PrRef-Pr)/PrRef);
  
  [gaint, gainp] = vf_computeGain(z0, EtFF, EpFF, PrRef);
                  % set the same power to see matching (radiated power may 
                  % not be the same)
  [refGaint, refGainp] = vf_computeGain(z0, refEtFF, refEpFF, PrRef);

    
  errTxz(j) = getL2error(EtFF(idx,1), refEtFF(idx,1));
  errPxz(j) = getL2error(EpFF(idx,1), refEpFF(idx,1));
  errTyz(j) = getL2error(EtFF(idx,2), refEtFF(idx,2));
  errPyz(j) = getL2error(EpFF(idx,2), refEpFF(idx,2));
  
  error(j) = getL2error(EtFF+EpFF, refEtFF+refEpFF);
  errT(j) = getL2error(EtFF, refEtFF);
  errP(j) = getL2error(EpFF, refEpFF);
    
  fprintf('t = %2.4g° p = %2.4g° => Error = \n\t',testT(j), testP(j));
  fprintf('tXZ %2.4g, pXZ %2.4g, | tYZ %2.4g |, pYZ %2.4g\n', errTxz(j), ...
      errPxz(j), errTyz(j), errPyz(j) );    

  vf_plotFFCutPlanes(thetaFF, gaint, gainp, 60, 2, ...
    thetaFF(1:2:end,:), refGaint(1:2:end,:), refGainp(1:2:end,:), ...
    dualSpanT);
  
%   figure(gcf); printEPS('','dual_pPol');
%   figure(gcf-1); printEPS('','dual_tPol');
      
  pause(1);
  close all;

end

fprintf('#> Testing time : %2.4g s\n', toc(testTime));

%% error plots
figure; semilogy(testT, error, '.-',testT, errT, '.-',testT, errP, '.-');
legend('total', '\theta_{pol}', '\phi_{pol}');
printEPS('','errorDirect');

figure; semilogy(testT, errTxz, '.-',testT, errPxz, '.-');
title('\phi = 0° plane');
legend('\theta_{pol}', '\phi_{pol}');
printEPS('','errorInverse00');


figure; semilogy(testT, errTyz, '.-',testT, errPyz, '.-');
title('\phi = 90° plane');
legend('\theta_{pol}', '\phi_{pol}');
printEPS('','errorInverse90');
