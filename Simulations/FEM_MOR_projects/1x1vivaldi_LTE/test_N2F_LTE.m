%%% LTE code for N2F
clc; close all;
addpath(genpath('..\_matlabSolvers'));
addpath(genpath('..\_matlabSolvers\PardisoInterface'));

nbrPorts = 1;
flag = 2; % 1: solve input system, 2: load previously solved
nbrCoeffs = 31;
nbrSmpls = 200;
freq = 7.5e9;
phi = [0,pi/2];
% phi = linspace(0,pi,61);

%% get system matrices
if flag == 1 || flag == 2
  [z0, k0] = getFreeSpaceElectricalParams(freq);
  [boxPos, boxN, dS] = ...
    buildBox([1 1 1 1 1 1], -.01, .015, -.01, .035, -0.003, .087, ...
    12, 22, 45, .99, 1, 0);
  [A, b, x] = getInputSystemMatrices(flag, nbrPorts);
  [ct, cp, F, sizetref, nbrCoeffs] = getOutputSystemMatrices(flag, ...
    nbrCoeffs, nbrSmpls, k0, z0, boxPos, boxN, dS, phi);
end

sizet = 3000;
[ifftOp, phiFF, thetaFF] = getLookAnglePoly(phi, sizet, sizetref, ...
  nbrCoeffs);

%% computing radiation pattern
handle = figure();

Pr = computeRadiatedPower(F*x, boxN, dS);
EtFF=zeros(size(thetaFF,1),size(phiFF,2));
EpFF=zeros(size(thetaFF,1),size(phiFF,2));

for i=1:size(phi,2)
  EtFF(:,i) = ifftOp*(ct{i}*x);
  EpFF(:,i) = ifftOp*(cp{i}*x);
end

[gaint, gainp] = vf_computeGain(z0, EtFF, EpFF, Pr);
vf_plotFFPolarCutPlanes(handle, 25.8, 10, thetaFF, gaint, gainp, 'N2F', ...
  '\bf Vivaldi pattern \rm (directivity)');
% vf_plotFF3d(handle, thetaFF, phiFF, 25.8, gaint, gainp);