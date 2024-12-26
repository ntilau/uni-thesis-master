% Returns the effective power radiated by the structure computing half of
% the real part of the Poynting vector flowing through the surface with
% radiation boundary conditions
%
% Pr = computeRadiatedPower(Fields, boxN, dS)
%
% IN: Fields = near fields collected in a vector |Hx(1)| sampling point 1
%                                                |Hy(1)|
%                                                |Hz(1)|
%                                                |Hx(2)| sampling point 2
%                                                |  .  |
%                                                |  .  |
%                                                |  .  |
%                                                |Hz(N)| sampling point N
%                                                |Ex(1)| sampling point 1
%                                                |Ey(1)|
%                                                |Ez(1)|
%                                                |Ex(2)| sampling point 2
%                                                |  .  |
%                                                |  .  |
%                                                |  .  |
%                                                |Ez(N)| sampling point N
%     n = outwardly directed from the antenna region normal unit vectors to
%         the surface at the sampling points
%     dS = surface patch area
%
% OUT: Pr = power radiated in [W]
function Pr = computeRadiatedPower(Fields, n, dS)

mH = Fields(1:floor(size(Fields,1)/2),1);
mE = Fields(floor(size(Fields,1)/2)+1:end,1);
E(1,:) = mE(1:3:end);
E(2,:) = mE(2:3:end);
E(3,:) = mE(3:3:end);
H(1,:) = mH(1:3:end);
H(2,:) = mH(2:3:end);
H(3,:) = mH(3:3:end);
S = cross(E,conj(H));
Sr = dot(S,n);
Pr=1/2*real(sum(Sr.*dS));

fprintf('##> Radiated power : %2.4g W\n', Pr);

