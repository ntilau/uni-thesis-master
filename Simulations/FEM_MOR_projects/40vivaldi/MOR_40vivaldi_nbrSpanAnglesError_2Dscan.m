%% MOR on Vivaldi antennas array
clc; close all;% clear all;
addpath(genpath('..\_matlabSolvers'));
addpath(genpath('..\_matlabSolvers\PardisoInterface'));

flag = 2; % 1: solve all 2: load and build output matrices
nbrPorts = 40;
nbrCoeffs = 25;
nbrSmpls = 200;
freq = 3e9;
phi = [0 pi/2];
phiPlane = 0; % [°]
spanAnglesSpace = 1:20; % nbr of spanning vectors per iteration
%% get system matrices
if flag == 1 || flag == 2
  [z0, k0, lambda0] = getFreeSpaceElectricalParams(freq);
  [boxPos, boxN, dS] =...
    buildBox([1 0 1 1 1 1], -.09, .09, -.09, .09, 0, .135,...
    40, 40, 30, 1, 0, 0);
  [fullA, fullB, fullX] = getInputSystemMatrices(flag, nbrPorts);
  [fullCt, fullCp, fullF, sizetref, nbrCoeffs] = ...
    getOutputSystemMatrices(flag, ...
    nbrCoeffs, nbrSmpls, k0, z0, boxPos, boxN, dS, phi);
end
%% Lumped ports order
equiRowX = [9 9 9 9 6 8 2 4 3 3 3 3 5 5 5 5 7 7 7 7 1 1 1 1 6 8 2 4 6 8 ...
  2 4 6 8 2 4 6 8 2 4 ];
equiRowY = [6 8 2 4 9 9 9 9 6 8 2 4 6 8 2 4 6 8 2 4 6 8 2 4 1 1 1 1 3 3 ...
  3 3 5 5 5 5 7 7 7 7 ];
phase = ones(1, nbrPorts);
alpha = 2*pi*0.01/lambda0; % physical spacing of 1 cm
%% MOR space selection

meanError = zeros(size(spanAnglesSpace));
maxError = meanError;

handle = figure;

for N = spanAnglesSpace;
  
  spanT = linspace(0,90,N);%% model order reduction
  spanP = phiPlane*ones(size(spanT));

  spanningV = zeros(size(fullA,1),size(spanT,2));
  for j = 1:length(spanT)
    for k=1:nbrPorts
      spanningV(:,j) = spanningV(:,j) + fullX(:,k) .* phase(k) * ...
        exp(-1i*alpha* ( (equiRowX(k)-1).* sin(deg2rad(spanT(j))) .* ...
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
  
  %% testing for chosen angle
  testTime = tic;

  testT = linspace(0,90, max(spanAnglesSpace));
  testP = phiPlane*ones(size(testT));
  % plotSelectedAngles(spanT, spanP, testT, testP);

  error = zeros(size(testT));

  romX = romA\romB;
  for j=1:size(testT,2)

    romXscan = zeros(size(V,2),1);
    for k=1:nbrPorts
      romXscan = romXscan + romX(:,k) .* phase(k) * ...
        exp(-1i*alpha* ( (equiRowX(k)-1).* sin(deg2rad(testT(j))) .* ...
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
        exp(-1i*alpha* ( (equiRowX(k)-1).* sin(deg2rad(testT(j))) .* ...
        cos(deg2rad(testP(j))) +  (equiRowY(k)-1) .* ...
        sin(deg2rad(testT(j))) .* sin(deg2rad(testP(j))) ) );
    end

    refEtFF=zeros(size(thetaFF,1),size(phiFF,2));
    refEpFF=zeros(size(thetaFF,1),size(phiFF,2));

    for i=1:size(thetaFF,2)
      refEtFF(:,i) = ifftOp*(fullCt{i}*fullXscan);
      refEpFF(:,i) = ifftOp*(fullCp{i}*fullXscan);
    end

    Pr = computeRadiatedPower(romF * romXscan, boxN, dS);
    PrRef = computeRadiatedPower(fullF * fullXscan, boxN, dS);
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

%   figure; semilogy(testT, error, '.-', 'LineWidth', 1.5, 'MarkerSize',7);
%   xlabel('Tested angle', 'FontSize', 12);
%   ylabel('Relative error', 'FontSize', 12);
%   printEPS('',['VivaldiErrorOrder',num2str(N)]);

  meanError(N) = mean(error);
  maxError(N) = max(error);
  
end
%%
figProp = getFigureProperties();
figure; semilogy(spanAnglesSpace, meanError(spanAnglesSpace), '*-', ...
  spanAnglesSpace, maxError(spanAnglesSpace), 'k*-',...
  'LineWidth',figProp.lw, 'MarkerSize', figProp.ms);
axis('tight');
xlabel('ROM order', 'FontSize', figProp.fs);
ylabel('Relative error', 'FontSize', figProp.fs);
legend('Avg error', 'Max error');