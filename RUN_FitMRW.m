%%   [GZM] Inducible Transcription Factors   %%
% ------------------------------------------- %
% RUN: Find parameters that best fit the data %

% Created by Mariana Gómez-Schiavon
% August 2019

%   See also FN_SS_SimpleHill.m
%   See also FN_SS_HillxBasal.m
%   See also FN_SS_Mechanistic.m
%   See also FN_FitError.m
%   See also FN_FitMRW.m

clear

%% Inputs & conditions
% Transcriptional model:
M = 'SimpleHill';      % Options: 'Mechanistic', 'HillxBasal', 'SimpleHill'
% Data
ExID = 'GEMc';	% Experiment (TF) to consider
    load('DATA_synTF.mat','xd');
    X = mean(xd.(ExID).X,1);
    H = xd.(ExID).H;
    D = xd.(ExID).Y;
    clear xd
S = [1:1000];     % Random number seed(s) (1 per run)
I = 20000;      % Iterations per fitting run
printAll = 0;   % Flag for printing full random walk
% Kinetic parameters:
    p.nM = 10;
    p.m  = 0.1;
    p.a  = 0.003;
    p.n  = 1.6;
    p.K  = 1;
    p.g  = 0.01;
% Parameters to fit:
    i = 0;
    i = i + 1;
    f(i).par = 'nM';
    f(i).cov = 0.1;
    f(i).lim = [1e-3,1000];
    i = i + 1;
    f(i).par = 'm';
    f(i).cov = 0.1;
    f(i).lim = [2e-6,2];
    i = i + 1;
    f(i).par = 'a';
    f(i).cov = 0.1;
    f(i).lim = [2e-7,0.2];
    i = i + 1;
    f(i).par = 'n';
    f(i).cov = 0.1;
    f(i).lim = [1e-5,10];
    i = i + 1;
    f(i).par = 'K';
    f(i).cov = 0.1;
    f(i).lim = [1e-4,100];
    clear i

%% Run fitting:
bestP = zeros(length(S),length(f));
minE  = zeros(length(S),1);
for s = S
    cat(2,'Running seed #',num2str(s))
    [bP,mE] = FN_FitMRW(X,H,p,M,D,s,f,I,ExID,printAll)
    bestP(s,:) = bP;
    minE(s) = mE;
    if(mod(s,10)==0)
        save(cat(2,'TEMP_MRW_',M,'_',ExID,'.mat'));
    end
end
clear s bP mE ans
save(cat(2,'MRW_',M,'_',ExID,'.mat'));
delete(cat(2,'TEMP_MRW_',M,'_',ExID,'.mat'));
% load(cat(2,'MRW_',M,'_',ExID,'.mat'));

%% Figures
if(printAll)
    fig = figure();
    fig.Units = 'inches';
    fig.PaperPosition = [1 0 18 10];
    fig.Position = fig.PaperPosition;
    C = colormap('parula');
    for s = S
        load(cat(2,'MRW_',M,'_',ExID,'_s',num2str(s),'.mat'),'mrw');
        for i = 1:size(mrw.P,2)
            subplot(2,4,i)
            hold on;
            plot(mrw.P(:,i),'LineWidth',2,'Color',C(s*6,:))
                xlabel('Iterations')
                ylabel(f(i).par)
                xlim([0,I])
                set(gca,'YScale','log')
                box on
        end
        subplot(2,4,8)
        hold on;
        plot(mrw.e,'LineWidth',2,'Color',C(s*6,:))
            xlabel('Iterations')
            ylabel('Error')
            xlim([0,I])
            set(gca,'YScale','log')
            box on
    end
    clear s i a b
    print(gcf,cat(2,'MRW_',ExID,'_Runs.png'),'-dpng','-r300')
else
    fig = figure();
    fig.Units = 'inches';
    fig.PaperPosition = [0 0 5 8];
    fig.Position = fig.PaperPosition;
    
    hist(log10(minE),20)
        xlabel('log_{10}(min(error))')
        ylabel('Count')
        title(cat(2,'Experiment: ',ExID))
        box on
        
        axes('Position',[0.25 0.5 0.3 0.3])
            hist(log10(minE([minE<10])),20)
            xlabel('log_{10}(min(error))')
            ylabel('Count')
            title(cat(2,'min(min(error)) = ',num2str(min(minE))))
            box on
        print(gcf,cat(2,'MRW_',ExID,'_minE.png'),'-dpng','-r300')
        
    fig = figure();
    fig.Units = 'inches';
    fig.PaperPosition = [0 0 18 8];
    fig.Position = fig.PaperPosition;
    [a b] = sort(minE);
    for i = 1:length(f)
        subplot(2,ceil(length(f)/2),i)
        hold on;
            scatter(minE,bestP(:,i))
            for ii = 1:5
                scatter(a(ii),bestP(b(ii),i),'filled')
            end
                ylabel(f(i).par)
                xlabel('min(error)')
                set(gca,'XScale','log')
                box on
                grid on
    end
        print(gcf,cat(2,'MRW_',ExID,'_minExPar.png'),'-dpng','-r300')
end
    fig = figure();
    fig.Units = 'inches';
    fig.PaperPosition = [0 0 18 8];
    fig.Position = fig.PaperPosition;
    C = colormap('jet');
    Hi = log2(H); Hi(length(H)) = Hi(length(Hi)-1) - 1;
    [a b] = sort(minE);
    for ii = 1:10
        for i = 1:length(f)
            p.(f(i).par) = bestP(b(ii),i);
        end
        if(strcmp(M,'SimpleHill'))
            Ye = FN_SS_SimpleHill(X*p.nM,H,p);
        elseif(strcmp(M,'HillxBasal'))
            Ye = FN_SS_HillxBasal(X*p.nM,H,p);
        elseif(strcmp(M,'Mechanistic'))
            Ye = FN_SS_Mechanistic(X*p.nM,H,p);
        else
            'ERROR: Transcriptional model not defined. Options: SimpleHill, HillxBasal, Mechanistic.'
        end
        subplot(2,5,ii)
        hold on
        for i = 1:length(X)
            plot(Hi,Ye(:,i),'Color',C(i*10,:),'LineWidth',2)
            plot(Hi,D(:,i)*p.nM,'Color',C(i*10,:),'LineStyle','none','Marker','o')
        end
                xlabel('Hormone')
                ylabel('Output')
                title(cat(2,'Error = ',num2str(FN_FitError(X*p.nM,H,p,M,D*p.nM))))
                xlim([min(Hi)-0.5 max(Hi)+0.5])
                set(gca,'YScale','log','YGrid','on',...
                    'XTick',Hi([length(H):-3:1]),'XTickLabel',H([length(H):-3:1]),...
                    'XTickLabelRotation',45)
                box on
    end
    clear Hi s i 
    print(gcf,cat(2,'MRW_',M,'_',ExID,'_BestFits.png'),'-dpng','-r300')

%% END