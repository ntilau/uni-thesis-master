% Returns the output system matrices and the solutions-to-near-fields
% functional to compute the radiated power (directivity)
%
% [ct, cp, F, nbrSmpls, nbrCoeffs] = getOutputSystemMatrices(flag, ...
%   nbrCoeffs, nbrSmpls, k0, z0, boxPos, boxN, dS, phi)
%
% IN: flag = 1: dumps sampling points coordinates for the
%               solutions-to-fields operator computation
%            2: only computes the N2F operators with DFT-truncation
%     nbrCoeffs = number of DFT coefficients to retain
%     nbrSmpls = number of theta angles to consider for the DFT computation
%     k0, z0 = free-space wavenumber and impedance
%     surfPos = coordinates of the sampling points
%     n = normal unit vectors outwardly directed from the region of the
%         antennas
%     dS = surface patches areas (midpoint integration)
%     phi = phi angles of the constant phi planes to consider
%
% OUT: ct = solutions-to-EtFF (far electric field in theta pol.) operator
%      cp = solutions-to-EpFF (far electric field in phi pol.) operator
%      F = solutions-to-near-fields operator for radiated power computation
%      nbrSmpls = initial number of theta angles for trigometric 
%                 polynomials values computation
%      nbrCoeffs = effective number of DFT coefficients retained (odd
%                  number)
%
% Laurent Ntibarikure
function [ct, cp, F, nbrSmpls, nbrCoeffs] = getOutputSystemMatrices(flag, ...
  nbrCoeffs, nbrSmpls, k0, z0, surfPos, n, dS, phi)
%% dump fieldEvaluationPoints
if flag == 1
  mP = surfPos.';
  writeMatFull(mP, 'lte_fileset\FieldEvaluationPoints.fmat');
  beep;
  disp('#> Sampling coordinates dumped...');
  pause;
end
%% Output matrices
if flag == 1 || flag == 2
  disp('#> Loading Output System matrices ...');
  comp2Ext = mmread('lte_fileset\comp2Ext1.mm');
  num2UnNum = mmread('lte_fileset\num2UnNum1.mm');
  mFE = mmread('lte_fileset\functionalE.mm');
  mFH = mmread('lte_fileset\functionalH.mm');
  fieldFunctional  = [mFH;mFE];
  [phiFF,thetaFF] = meshgrid(phi,...
    linspace(0,2*pi*(nbrSmpls-1)/nbrSmpls,nbrSmpls));
  [Opt, Opp, nbrCoeffs] = vf_n2fOpFieldsFFT(k0, z0, surfPos, n, dS, ...
    thetaFF, phiFF, nbrCoeffs);
  coeffLevel = max(max(abs(Opt(floor(nbrCoeffs/2+1),:,:)))) / ...
    (max(max(abs(Opt(1,:,:)))));
  fprintf('##> nbrCoeffs = %d, level = %2.4g\n', ...
    nbrCoeffs, coeffLevel);
  ct = cell(size(phi));
  cp = ct;
  for i=1:size(phi,2)
    ct{i} = sparse(Opt(:,:,i)*fieldFunctional*num2UnNum*comp2Ext);
    cp{i} = sparse(Opp(:,:,i)*fieldFunctional*num2UnNum*comp2Ext);
  end
end
F = sparse(fieldFunctional*num2UnNum*comp2Ext);