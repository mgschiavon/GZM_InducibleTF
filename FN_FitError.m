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
%   M : Transcriptional model to consider
%       ('Mechanistic','HillxBasal','SimpleHill')
%   D : Measured output (data) matrix [length(H) x length(X)]
%
%   OUTPUT Ef : Sum of square errors
%
%   See also FN_SS_SimpleHill.m
%   See also FN_SS_HillxBasal.m
%   See also FN_SS_Mechanistic.m

function Ef = FN_FitError(X,H,p,M,D)
    if(strcmp(M,'SimpleHill'))
        Ye = FN_SS_SimpleHill(X,H,p);
    elseif(strcmp(M,'HillxBasal'))
        Ye = FN_SS_HillxBasal(X,H,p);
    elseif(strcmp(M,'Mechanistic'))
        Ye = FN_SS_Mechanistic(X,H,p);
    else
        'ERROR: Transcriptional model not defined. Options: SimpleHill, HillxBasal, Mechanistic.'
    end
    Ef = sum(sum((log10(D)-log10(Ye)).^2)./(2*var(log10(D))));
end
