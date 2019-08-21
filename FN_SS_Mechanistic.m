%%   [GZM] Inducible Transcription Factors   %%
% ------------------------------------------- %
% FUNCTION: Finding steady state solution     %
%           for the mechanistic model         %

% Created by Mariana Gómez-Schiavon
% August 2019

% FN_SS_Mechanistic : Finding steady state solution for the mechanistic
%                     model given a set of biophysical paramers and
%                     inducer (hormone) concentration.
%
%   Xs = FN_SS_Mechanistic(X,H,p)
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
%
%   OUTPUT Ye : Matrix [length(H) x length(X)]
%
%   See also FN_SS_SimpleHill.m
%   See also FN_SS_HillxBasal.m
%   See also FN_FitError.m

function Ye = FN_SS_Mechanistic(X,H,p)
    OHFn = @(x,n,k,a) [a+((1-a)*(x.^n)./((x.^n)+(k^n)))];
    OImx = @(x,I,b,m) [I*(1-(m*exp(-b*x/I)))];
    Ye = zeros(length(H),length(X));
    for h = 1:length(H)
        for i = 1:length(X)
            XT = X(i);
            Xa = roots([1,-(H(h)+XT+p.KX),H(h)*XT]);
            Xa = Xa([Xa<XT]);
            if(length(Xa)~=1)
                'error -- multiple solutions'
            end
            Xo = OHFn(Xa+(p.aX*(XT-Xa)),p.nO,p.KO,0);
            Ye(h,i) = OImx(Xo,p.Im,p.kb,p.mY)/p.gY;
        end
    end
    clear h i
end