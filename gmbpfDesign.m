% Design GMBPF
more off
% Read parameters
gmbpfPara;

xiSL  = pi * 0.6 * ((-2):2) / sigma;   % Gussian center angular frequency (0.24 is obtained hylistic.)
Ex         = 5;                        % Extension for error calculation
passV      = power(10, -3/20);         % -3 dB

% Set the scan interval of K and stop band attenuation
if P < 16
  KL    = floor(2.6 * sigma):floor(3.8 * sigma);
  stopV = 1.0e-1;                        % 10e-3 for 6dB  10e-5 for 100dB
elseif P == 16
  KL    = floor(3.6 * sigma):floor(3.8 * sigma);
  stopV = 1.0e-3;                        % 10e-3 for 6dB  10e-5 for 100dB
elseif P < 27
  KL    = floor(3.6 * sigma):floor(5.0 * sigma);
  stopV = 1.0e-3;                        % 10e-3 for 6dB  10e-5 for 100dB
elseif P == 27
  KL    = floor(4.8 * sigma):floor(5.0 * sigma);
  stopV = 1.0e-5;                   % 10e-3 for 6dB  10e-5 for 100dB
else
  print("Search by youself \n");
  exit
end

scan = 1; 
nsPL = 11; % # of sp scan
sPB  = 0;  % sP scan bias

% Scan mu
quadL  = [1.0 sqrt(sqrt(10.0)) sqrt(10.0) sqrt(10.0 * sqrt(10.0))];
muBase = 1.0;
muL    = [0.0 muBase * quadL];

fName = sprintf("Data/errorListSft%d_D_P%02d.txt", shiftPixel, P);
fpE = fopen(fName, "a");
fprintf(fpE, "# xiC   sigma        totErr(percent)   K    shift   sP   eP  mu  edgeErr \n");

%Org filter
xiL         = xiSL + xiC;
KTgt        = floor(6  * sigma); % Long length is used for original filter.
NTgt        = 2 * KTgt + 1;
filterTgtL  = zeros(1, NTgt);
beta        = 1.0 / (2.0 * sigma * sigma);
nL          = (-KTgt):KTgt;

% Make target GMBPF 
for xi = xiL
  filterTgtL = filterTgtL + exp(- beta * (nL .* nL)) .* (cos((xi) * nL));
endfor
hightFilter = filterTgtL * cos((xiC) * nL)';

while 1
  gmbpfApprox
  if size(wlst, 2) >= 1
    K   = minResult(1);
    sP  = minResult(3);
    muL = minResult(7);
    if K > KL(2) && K < KL(end - 1)
      sPS   = floor(xiC * K / pi - P / 2) + sPB;
      nsPLH = (nsPL - 1) / 2;
      if sP > sPS - nsPLH + 1 && sP < sPS + nsPLH - 1
	break;
      else
	sPB = sPB + nsPLH;
      end
    else
      KL = (K-5):(K+5);
    end
  else
    muBase = muBase * 10.0;
    muL    = [muBase * quadL];
  end
end
pResult = minResult;
fprintf(fpE, "%9.6f  %1d %5d  %7d %4d %6d %4d %6.2e  %6.2e %6.2e\n", pResult(5), 3, pResult(6), pResult(1), pResult(2), pResult(3), pResult(4) - pResult(3) + 1, pResult(7), pResult(10), pResult(11));
fprintf("Min: %9.6f  %1d %5d  %7d %4d %6d %4d %6.2e  %6.2e %6.2e\n", pResult(5), 3, pResult(6), pResult(1), pResult(2), pResult(3), pResult(4) - pResult(3) + 1, pResult(7), pResult(10), pResult(11));

% Reculculate minimum filter
scan = 0; % Not scan
KL   = [K];
sPL  = [sP];
muL  = [mu];
gmbpfApprox

% Coefficients of ASFT filter
coefFilterSg = zeros(1, 2 * P);
if rem(sP, 2) == 0
  sSg = 1;
else
  sSg = -1;
end

sg = sSg;
for posC = 1:P
  coefFilterSg(posC)     =   sg * coefFilter(posC);
  coefFilterSg(posC + P) = - sg * coefFilter(posC + P);
  sg = -sg;
end

				% Save data necessary for fiter
thetaPL = (sP:eP) * piByK;
cosPara = cos(thetaPL) * exp(- alphaPixel);
sinPara = sin(thetaPL) * exp(- alphaPixel);
subPara = exp(- 2 * K * alphaPixel);

fName = sprintf("Data/save_%5.3f_%.0f_Sft%d_P%02d.mat", xiC / pi, sigma, shiftPixel, P);
save(fName, 'K', 'subPara', 'cosPara', 'sinPara', 'coefFilterSg');

% Output feature of the obtained filter if flagPlot == 1
if flagPlot == 1
  gmbpfFeature
end

