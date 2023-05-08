function p = Phase_Unwrap(p,cutoff)
%LocalUnwrap   Unwraps column vector of phase values.

m = length(p);

% Unwrap phase angles.  Algorithm minimizes the incremental phase variation 
% by constraining it to the range [-pi,pi]
dp = diff(p,1,1);                % Incremental phase variations

% Compute an integer describing how many times 2*pi we are off:
% dp in [-pi, pi]: dp_corr = 0,
% elseif dp in [-3*pi, 3*pi]: dp_corr = 1,
% else if dp in [-5*pi, 5*pi]: dp_corr = 2, ...
dp_corr = dp./(2*pi);

% We want to do round(dp_corr), except that we want the tie-break at n+0.5
% to round towards zero instead of away from zero (that is, (2n+1)*pi will
% be shifted by 2n*pi, not by (2n+2)*pi):
roundDown = abs(rem(dp_corr, 1)) <= 0.5;
dp_corr(roundDown) = fix(dp_corr(roundDown));

dp_corr = round(dp_corr);

% Stop the jump from happening if dp < cutoff (no effect if cutoff <= pi)
dp_corr(abs(dp) < cutoff) = 0;

% Integrate corrections and add to P to produce smoothed phase values
p(2:m,:) = p(2:m,:) - (2*pi)*cumsum(dp_corr,1);
