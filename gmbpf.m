% Calculate GMBPF
more off
% Read parameters
gmbpfPara

% Input signal
% inSig        = zeros(1, 100000);
% inSig(50000) = 1;
% inSig        = cos((xiC + 0.001 * pi) * (1:100000));

% Load coefficients of filter
fName = sprintf("Data/save_%5.3f_%.0f_Sft%d_P%02d.mat", xiC / pi, sigma, shiftPixel, P);
load(fName)

% Allocate memory for output
K2     = K + K;
L      = size(inSig, 2);
outSig = zeros(1, L + K2);
if L < K2
  inSig = [inSig zeros(1, K2 - L)];
  L = K2;
end

% Apply Filter
inteCos = zeros(1, P);
inteSin = zeros(1, P);
for pos = 1:K2
  add  = inSig(pos);

  inteCosTmp = inteCos;
  inteCos = cosPara .* inteCosTmp - sinPara .* inteSin + add;
  inteSin = sinPara .* inteCosTmp + cosPara .* inteSin;
  outSig(pos) = sum(coefFilterSg .* [inteCos inteSin]);
end

for pos = (K2 + 1):L
  add = inSig(pos);
  sub = subPara * inSig(pos - K2);

  inteCosTmp = inteCos;
  inteCos = cosPara .* inteCosTmp - sinPara .* inteSin + add;
  inteSin = sinPara .* inteCosTmp + cosPara .* inteSin;
  outSig(pos) = sum(coefFilterSg .* [inteCos inteSin]);
  inteCos = inteCos - sub;
end  

for pos = (L + 1):(L + K2)
  sub = subPara * inSig(pos - K2);

  inteCosTmp = inteCos;
  inteCos = cosPara .* inteCosTmp - sinPara .* inteSin + add;
  inteSin = sinPara .* inteCosTmp + cosPara .* inteSin;
  outSig(pos) = sum(coefFilterSg .* [inteCos inteSin]);
  inteCos = inteCos - sub;
end  

