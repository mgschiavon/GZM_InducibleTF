%%   [GZM] Inducible Transcription Factors   %%
% ------------------------------------------- %
% FUNCTION: Calculate fitting error           %

% Created by Mariana Gómez-Schiavon
% August 2019

% FN_FitError : Calculate the sum of square errors between model's 
%               steady state (given a set of biophysical paramers and 
%               inducer concentration) and observed data.
%
%   Ef = FN_FitError(X,H,p,D)
%   X : Vector (array) of TF concentration
%   H : Vector of inducer (hormone) concentration
%   p : Structure with the kinetic parameters
%    .mY : 
%    .kb : 
%    .nO : 
%    .KO : 
%    .KX : 
%    .aX : 
%    .Im : 
%    .gY :
%   D : Measured output (data) matrix [length(H) x length(X)]
%
%   OUTPUT Ef : Sum of square errors
%
%   See also FN_SS_SimpleHill.m
%   See also FN_SS_HillxBasal.m
%   See also FN_SS_Mechanistic.m

function Ef = FN_FitError(X,H,p,D)
    Ye = FN_SS_Mechanistic(X,H,p);
    Ef = sum(sum((log10(D)-log10(Ye)).^2)./(2*var(log10(D))));
end
