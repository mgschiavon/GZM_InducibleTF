%%   [GZM] Inducible Transcription Factors   %%
% ------------------------------------------- %
% FUNCTION: Perform MCMC parameter estimation %

% Created by Mariana Gómez-Schiavon
% August 2019

% FN_FitMRW : Find the set of parameters that best fit the data.
%
%   [] = FN_FitMRW(X,H,p,D,s,fp,ExID)
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
%   s : Random number generator seed
%   f : Structure (array) with information to fit parameters
%    .par : Parameter to fit (e.g. 'mY')
%    .cov : Covariance to calculate parameter perturbations (e.g. 1e-3)
%    .lim : Range of acceptable values (e.g. [0,1])
%   I : Number of iterations
%   ExID : Code for output file name
%
%   OUTPUT Ef : Sum of square errors
%
%   See also FN_SS_SimpleHill.m
%   See also FN_SS_HillxBasal.m
%   See also FN_SS_Mechanistic.m
%   See also FN_FitError.m

function [] = FN_FitMRW(X,H,p,D,s,f,I,ExID)
    mrw.s = s;
    mrw.f = f;
    mrw.P = zeros(I,length(f));     % OUTPUT: Parameters.
    mrw.e = zeros(I,1);             % OUTPUT: Error function values.

    % (1) Define random number generator:
    rng(s,'twister');
    r.tL = rand(I,1);               % To evaluate proposal acceptance.
    r.Pe = zeros(I,length(f));      % Parameter perturbations.
    for i = 1:length(f)
        r.Pe(:,i) = mvnrnd(zeros(I,1),f(i).cov);
        % (2) Initialize system:
        mrw.P(1,i) = (rand()*(f(i).lim(2) - f(i).lim(1))) + f(i).lim(1);
        p.(f(i).par) = mrw.P(1,i);
    end
    mrw.e(1)   = FN_FitError(X,H,p,D);

    % (3) Iterate:
    for j = 2:I
        mrw.P(j,:) = mrw.P(j-1,:);
        mrw.e(j)   = mrw.e(j-1,:);
        % Alternative parameter set:
        for i = 1:length(f)
            p.(f(i).par) = min(max(f(i).lim(1),mrw.P(j,i) + r.Pe(j,i)),f(i).lim(2));
        end
        % Generate proposal:
        myE = FN_FitError(X,H,p,D);
        % If proposal is accepted, update system:
        if(r.tL(j) < exp(mrw.e(j)-myE))
            for i = 1:length(f)
                mrw.P(j,i) = p.(f(i).par);
            end
            mrw.e(j) = myE;
        end
        % Save progress:
        if(mod(j,10000)==0)
            j0 = j + 1
            save('TEMP_MRW.mat','mrw','j0','r');
        end
    end
    clear j myP myE

    % (4) Save:
    save(cat(2,'MRW_',ExID,'_s',num2str(s),'.mat'),'mrw');
    delete('TEMP_MRW.mat');
end
