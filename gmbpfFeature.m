% Output Feature GMBPF
more off

% File of data to plot functions (of the last w). (Approximated and True) 
figure(1)
plot(thetaL, filterL, thetaL, filterAL);

dtftT         = (-KEx):KEx;
dtftTA        = (-K):K;
dtftTADelay   = 0:(2 * K);

dtftOmegaL = 0:(0.04/sigma):pi;
nDtftOmega = size(dtftOmegaL, 2);
nDC        = floor((nDtftOmega - 1) * xiC / pi ) + 1;
nDL        = (-500:500) + nDC;
dtftDOmegaL = dtftOmegaL(nDL); % For extracted passBand part
		     %  dtftDOmegaL = 0:(1.0/sigma):pi; % For plot all

dtftC       = filterExL * (cos(dtftT'  * dtftDOmegaL));
dtftS       = filterExL * (sin(dtftT'  * dtftDOmegaL));
dtft        = sqrt(dtftC  .* dtftC  + dtftS  .* dtftS);

dtftCA      = filterAL  * (cos(dtftTA' * dtftDOmegaL));
dtftSA      = filterAL  * (sin(dtftTA' * dtftDOmegaL));
dtftA       = sqrt(dtftCA .* dtftCA + dtftSA .* dtftSA);

dtftCADelay = filterAL  * (cos(dtftTADelay' * dtftDOmegaL));
dtftSADelay = filterAL  * (sin(dtftTADelay' * dtftDOmegaL));

figure(2)
semilogy(dtftDOmegaL / pi, abs(dtft), dtftDOmegaL / pi, abs(dtftA));

angleLA = atan2(dtftSA, dtftCA);
figure(3)
plot(dtftDOmegaL / pi, angleLA);

figure(4)
angleLADelay = atan2(dtftSADelay, dtftCADelay);
gDelay       = angleLADelay(2:1001) - angleLADelay(1:1000);
gDelay       = (gDelay - 2 * pi * (gDelay > 4) + 2 * pi * (gDelay < - 4) ) / (0.04/sigma);
figure(4)
plot(dtftDOmegaL(1:1000) / pi, gDelay);


fName = sprintf("Data/dtft_%5.3f_%.0f_Sft%d_P%02d.txt", xiC / pi, sigma, shiftPixel, P);
fpDTFT = fopen(fName, "w");
for pos = 1:(size(dtftDOmegaL, 2) - 1)
  fprintf(fpDTFT, "%e %e %e  %f  %e %e \n", dtftDOmegaL(pos) / pi, 20 * log10(abs(dtftC(pos))), 20 * log10(abs(dtftCA(pos))), angleLA(pos), (dtftDOmegaL(pos + 1) + dtftDOmegaL(pos) )/ (2 * pi), gDelay(pos) );
end
				% Cal. bandwidth
dtftOmegaBWL  = 0:(0.04/sigma):pi;
nDtftOmegaBW  = size(dtftOmegaBWL, 2);
nCBW          = floor((nDtftOmegaBW - 1) * xiC / pi ) + 1;
nBWL          = (-500:500) + nCBW;
dtftOmegaBWSL =  dtftOmegaBWL(nBWL);
dtftBWSCA     = filterAL * (cos(dtftTA' * dtftOmegaBWSL));
dtftBWSSA     = filterAL * (sin(dtftTA' * dtftOmegaBWSL));
dtftBWSA      = sqrt(dtftBWSCA .* dtftBWSCA + dtftBWSSA  .* dtftBWSSA);
dtftBWSA      = dtftBWSA / max(dtftBWSA);
for pos = 2:size(nBWL, 2)
  if (dtftBWSA(pos - 1) < passV && dtftBWSA(pos) >= passV)
    pbwL = dtftOmegaBWSL(pos);
  end
  if (dtftBWSA(pos - 1) > passV && dtftBWSA(pos) <= passV)
    pbwH = dtftOmegaBWSL(pos);
  end
  if (dtftBWSA(pos - 1) < stopV && dtftBWSA(pos) >= stopV)
    sbwL = dtftOmegaBWSL(pos);
  end
  if (dtftBWSA(pos - 1) > stopV && dtftBWSA(pos) <= stopV)
    sbwH = dtftOmegaBWSL(pos);
  end
end
fprintf("sigma = %f xiC = %f K = %d sP = %d P = %d shift = %d toterror = %f\n", sigma, xiC, K, sP, P, shiftPixel, minResult(10));
fprintf("Pass band width (%3.0fdB) = %f pi, Stop band width (%3.0fdB) = %f pi  %15.10e pi  %15.10e pi  %15.10e pi  %15.10e pi\n", 20 * log10(passV), (pbwH - pbwL) / pi, 20 * log10(stopV), (sbwH - sbwL) / pi, pbwL/pi, pbwH/pi, sbwL/pi, sbwH/pi);
fprintf(fpDTFT, "# sigma = %f xiC = %f K = %d sP = %d P = %d shift = %d toterror = %f\n", sigma, xiC, K, sP, P, shiftPixel, minResult(10));
fprintf(fpDTFT, "# Pass band width (%3.0fdB) = %f pi, Stop band width (%3.0fdB) = %f pi  %15.10e pi  %15.10e pi  %15.10e pi  %15.10e pi\n", 20 * log10(passV), (pbwH - pbwL) / pi, 20 * log10(stopV), (sbwH - sbwL) / pi, pbwL/pi, pbwH/pi, sbwL/pi, sbwH/pi);
fclose(fpDTFT);

% Output Filter  
fName = sprintf("Data/filter_%5.3f_%.0f_Sft%d_P%02d.txt", xiC / pi, sigma, shiftPixel, P);
fp = fopen(fName, "w");
pos = -K + 1;
stepA = 0.0;
stepO = 0.0;
for posE = (K * (Ex - 2) + 1):(K * (Ex - 1));
  stepO = stepO + filterExL(posE);
  fprintf(fp, "%e %d %d %e %e %e %e \n", thetaExL(posE), posE - K * Ex - 1, pos - 1,  0.0, filterExL(posE), stepA, stepO);
  pos = pos + 1;
end
for posE = (K * (Ex - 1) + 1):(K * (Ex + 1) + 1);
  stepA = stepA + filterAL(pos);
  stepO = stepO + filterExL(posE);
  fprintf(fp, "%e %d %d %e %e %e %e \n", thetaExL(posE), posE - K * Ex - 1, pos - 1, filterAL(pos), filterExL(posE), stepA, stepO);
  pos = pos + 1;
end
for posE = (K * (Ex + 1) + 2):(K * (Ex + 2) + 1);
  stepO = stepO + filterExL(posE);
  fprintf(fp, "%e %d %d %e %e %e %e \n", thetaExL(posE), posE - K * Ex - 1, pos - 1, 0.0, filterExL(posE), stepA, stepO);
  pos = pos + 1;
end
fclose(fp);

% Output Coefficients of ASFT filter
fName = sprintf("Data/coef_%5.3f_%.0f_Sft%d_P%02d.txt", xiC / pi, sigma, shiftPixel, P);
fp    = fopen(fName, "w");

if shiftPixel ~= 0 % ASFT
for posC = 1:(2 * P)
  fprintf(fp, "%e\n", coefFilterSg(posC));
end
else % SFT
for posC = 1:P
  fprintf(fp, "%e\n", coefFilterSg(posC));
end
end

fclose(fp);

