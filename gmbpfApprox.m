%Calc Fitters

% KL, sPL, muL are necessary as well as P, xiC, 

if scan == 1
  resultTbl  = zeros(16, size(KL, 2) * nsPL * size(muL, 2));
else
  resultTbl  = zeros(16, 1);
end

pos = 1;
for K = KL
  KEx     = Ex * K;
  piByK   = pi / K;
  w       = sigma * piByK;
  sPC     = xiC / piByK;
  N       = 2 * K + 1;
  NEx     = 2 * KEx + 1;
  posL    = (K * (Ex - 1) + 1):(K * (Ex + 1) + 1);
  posExtL = [1:(K * (Ex - 1))  (K * (Ex + 1) + 2):NEx];  % Position for extension

  thetaL   = piByK * [-K:K];
  thetaExL = piByK * [-KEx:KEx];

  beta        = 1.0 / (2.0 * w * w);
  betaw       = 2.0 * beta;
  shiftTheta  = piByK * shiftPixel;
  alphaPixel  = shiftPixel / (sigma * sigma);
  thetaExSftL = thetaExL + shiftTheta;

  xiL = xiSL + xiC;
  filterExL  = zeros(1, NEx);
  for xi = xiL
    xiP       = xi / piByK;
    filterExL = filterExL + exp(- beta * (thetaExSftL .* thetaExSftL)) .* (cos(xiP * thetaExSftL));
  end
  filterExL    = filterExL / heightFilter;
  sqFilter     = filterExL * filterExL';   % Norm of extended fiter
  filterL      = filterExL(posL);          % Inter [-K, K] data
  filterExtL   = filterExL(posExtL);       % Outer [-K, K] data
  trErrFilter  = filterExtL * filterExtL'; % Trancation Error

  if scan == 1
    sPL  = (floor(sPC - P / 2) - (nsPL - 1) / 2 ):(floor(sPC - P / 2) + (nsPL - 1) / 2 ) + sPB;
  end
  
  for PN = 1:size(sPL, 2)
    sP     = sPL(PN);
    eP     = sP + P - 1;
    cosLL  = zeros(eP - sP + 1, N);
    sinLL  = zeros(eP - sP + 1, N);
    for p = sP:eP
      cosLL(p - sP + 1, :) = cos(p * thetaL);
      sinLL(p - sP + 1, :) = sin(p * thetaL);
    end
% Calc. ASFT
%    coefExp  = exp(- alphaPixel * ((-N + 1):0));
    coefExp  = exp(- alphaPixel * (0:(N - 1)));
    cosExpLL = zeros(eP - sP + 1, N);
    sinExpLL = zeros(eP - sP + 1, N);
    for p = sP:eP
      cosExpLL(p - sP + 1, :) = cosLL(p - sP + 1, :) .* coefExp;
      sinExpLL(p - sP + 1, :) = sinLL(p - sP + 1, :) .* coefExp;
    end
    if shiftPixel ~= 0 % ASFT
      triExpLL = [cosExpLL ; sinExpLL];
    else % SFT
      triExpLL = cosExpLL;
    end
    
    for mu = muL
      edgeCond   = triExpLL(:, [1 N]);                         % Edge height at 1 and 2K + 1
      coefFilter = inv(triExpLL * triExpLL' + mu * (edgeCond * edgeCond')) * triExpLL * filterL'; % Calc. filter coefficent
      
      filterAL    = coefFilter' * triExpLL;                       % Approximated data
      apErrFilter = (filterAL - filterL) * (filterAL - filterL)'; % Order error
      numerS      = 0.5 * (abs(filterAL(1)) + abs(filterAL(N)));  % Leap at edge point

      resultTbl(1, pos)  = K;
      resultTbl(2, pos)  = shiftPixel;
      resultTbl(3, pos)  = sP;
      resultTbl(4, pos)  = eP;
      resultTbl(5, pos)  = xiC;
      resultTbl(6, pos)  = sigma;
      resultTbl(7, pos)  = mu;
      resultTbl(8, pos)  = 100 * sqrt(trErrFilter / sqFilter);
      resultTbl(9, pos)  = 100 * sqrt(apErrFilter  / sqFilter); 
      resultTbl(10, pos) = 100 * sqrt((trErrFilter + apErrFilter) / sqFilter);
      resultTbl(11, pos) = abs(filterAL(N) / (filterL(N) - filterL(N - 1)) );
      %	resultTbl(12, pos) = numer / abs(mtSftL(2) - mtSftL(1));
      %	resultTbl(13, pos) = numer / abs(mtSftA(2) - mtSftA(1));
      %	resultTbl(14, pos) = invType;
      
      pResult = resultTbl(:, pos);
      % printf("%9.6f  %1d %5d  %7d %4d %6d %4d %6.2e  %6.2e %6.2e\n", pResult(5), 3, pResult(6), pResult(1), pResult(2), pResult(3), pResult(4) - pResult(3) + 1, pResult(7), pResult(10), pResult(11));
      pos = pos + 1;
      clear mtCSftA mtSSftA mtCSftL mtSSftL coefMtC coefMtS 
    end
  end
end

wlst = find(abs(resultTbl(11, :)) <= 1.2);
if size(wlst, 2) >= 1
  wResultTbl        = resultTbl(:, wlst);
  sTarget           = 10;       
  [minval minIndex] = min(wResultTbl(sTarget, :));  % Min error
  minResult         = wResultTbl(:, minIndex);      % Min data
end

clear resultTbl triExpLL cosLL sinLL sinExpLL cosExpLL 
