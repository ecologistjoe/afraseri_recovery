%     % Data from Observations
%     D.A2ARD = [  4     3     3
%                  7    48     3
%                 19     2     9
%                 10     7     5
%                  4     5     9
%                  8     3    32
%                  1     1    11
%                  1    16     4];
%     D.R2R = [3 2 0 8 3 0]';
%     D.sumR = [17 5 2 11 3 1]';
%     D.S0 = [350,400,75,550,325,900,125,100]';
%     D.S1 = [334,134,0,0,267,200,67,34]';
%     D.J0 = [13,213,38,550,13,175,88,150]';
%     D.J1 = [534,550,50,367,200,267,0,50]';
%     D.A0 = [30,22,5,57,9,22,23,47]';
%     D.R0 = [17,5,2,0,11,3,1,0]';

load fir_dbh
load count_data
m = {'CD', 'LC', 'MC', 'MG', 'MS'};
%s = {};
X(30,:) = X(29,:);


% Build a dataframe D for each mountain in m
for i = 1:5
    q = 0;
    clear D;
    for k = 1:8
        
        for t = 0:1
            j = t*8+k;
            
            L = strcmp(Mount, m{i}) & (PlotNum == k);
            if t == 0
                Y = DBH(L,1:2);
            else
                Y = DBH(L,2:3);
            end

            A0 = (Y(:,1)>0) & (Y(:,1)<15);
            A1 = (Y(:,2)>0) & (Y(:,2)<15);
            R0 = Y(:,1) >= 15;
            R1 = Y(:,2) >= 15;
            D1 = Y(:,2)==0;

            D.A2ARD(j,1) = sum(A0 & A1) ;
            D.A2ARD(j,2) = sum(A0 & R1) ;
            D.A2ARD(j,3) = sum(A0 & D1) ;

            D.sumR(j) = sum(R0);
            D.R2R(j)  = sum(R0 & R1);

            idx = strcmp(plots, [m{i} num2str(k)]);
            if all(X(idx, t+[1 2 4 5]) > -1) 
                q = q +1;
                
                D.S0(q) = X(idx,1+t);
                D.S1(q) = X(idx,2+t);
                D.J0(q) = X(idx,4+t);
                D.J1(q) = X(idx,5+t);
                D.A0(q) = sum(A0);
                D.A1(q) = sum(A1);
                D.R0(q) = sum(R0);
                D.R1(q) = sum(R1);
            end
            
        end
    end

    D.R2R(D.sumR == 0) = [];
    D.sumR(D.sumR == 0) = [];
    D.A2ARD(sum(D.A2ARD,2)==0,:) = [];
    
    % Data Constants
    D.alpha3 = [1 1 1];
    D.alpha4 = [1 1 1 1];
    D.NUM_SAMPLES = length(D.S0);
    D.NUM_ADULT_SAMP = size(D.A2ARD,1);
    D.NUM_REPROD_SAMP = length(D.R2R);

    I(1).Ps = [0.1 0.1 0.1 0.7];
    I(1).Pj = [0.1 0.1 0.1 0.7];

    I(1).Pa = [.6 .3 .1];
    I(1).Prr = .8;
    I(1).Bsr = 10;
    I(1).Bjr = 10;
    
    I(1).tau_S = .1;
    I(1).tau_J = .1;
    I(1).tau_A = .1;
    I(1).tau_R = .1;
    
    for c = 2:5
        r = rand(1,4);
        I(c).Ps = r/sum(r);
        
        r = rand(1,4);
        I(c).Pj = r/sum(r);
        
        r = rand(1,3);
        I(c).Pa = r/sum(r);
        
        I(c).Prr = rand(1);
        
        I(c).Bsr = rand(1)*20;
        I(c).Bjr = rand(1)*20;

        I(c).tau_S = rand(1)/100;
        I(c).tau_J = rand(1)/100;
        I(c).tau_A = rand(1)/100;
        I(c).tau_R = rand(1)/100;
    
    end
    
    
    workingDir = ['c:/projects/firmodel/' m{i}];
    [~,~] = mkdir(workingDir);
    
    % Run the BUGS model
    [samples, stats] = matbugs(D, 'c:/projects/firmodel/bugs_model.txt', ...
            'init', I, ...
            'view', 0, 'nburnin', 100000, 'nsamples', 20000, ...
            'thin', 1, 'nChains', 5, ...
            'monitorParams', {'z','G', 'tau_S', 'tau_J', 'tau_A', 'tau_R'}, ...
            'workingDir', workingDir, ...
            'openBugs', 1, ...
            'DICstatus',1, ...
            'Bugdir', 'c:/Program Files (x86)/OpenBUGS');
        
    s{i}.stats = stats;
    s{i}.samp = samples;
    


   % stats.mean.Q, [stats.mean.z, stats.mean.thin_s, stats.mean.thin_j]

   
    load abies_stage_by_plots
   
    % Generate summary statistics
    G = reshape(s{i}.samp.G,[], 4,4);
    for y = 1:3
       
        x = c{i}(:,:,y);
        x(:, any(x==-1)) = [];
        
        rr = zeros([length(G) 7]);
        tt = zeros([length(G) 7]);
        for k = 1:length(G)
            g = squeeze(G(k,:,:));
            for j = 0:7-y
               gjx = g^j*mean(x,2);
                              
               rr(k, j+y) = (sum(gjx(3:4,:)));
               tt(k, j+y) = (sum(gjx));
            end
        end

        
        r = rr(:,2:5) ./ rr(:,1:4);
        t = tt(:,2:5) ./ tt(:,1:4);
        
        s{i}.t(:,:,y) = t;
        s{i}.r(:,:,y) = r;
        s{i}.tt(:,:,y) = tt;
        s{i}.rr(:,:,y) = rr;
        
    end
        
        

    
    
end
%for i = 1:5; zs(i) = s{i}.stats.mean.z; end; zs

%% REFERENCE:
% MATBUGS a Matlab interface for WinBugs, similar to R2WinBUGS
% [samples, stats] = matbugs(dataStruct,  bugsModelFileName, ...)
%
% This generates a file called 'script.txt' in the directory where
% winbugs14.exe is kept (so you must have write access to this directory!),
% then calls winbugs with this script, then reads the resulting samples from
% files into matlab structs.
%
% INPUT:
% dataStruct contains values of observed variables.
% bugsModel is the name of the model file; MUST END IN .txt
% 
% Note: variables with names 'a.b' in the bugs model file
% should be called 'a_b' in the matlab data structure.
%
% Optional arguments passed as 'string', value pairs [default in brackets, case
% insensitive]:
%
% 'init' - init(i).v is a struct containing initial values for variable 'v'
%          for chain i.  Uninitialized variables are given random initial
%          values by WinBUGS.  It is highly recommended that you specify the
%          initial values of all root (founder) nodes.
% 'monitorParams' - cell array of field names (use 'a_b' instead of 'a.b')
%                   [defaults to *, which currently does nothing...]
% 'nChains'  - number of chains [3]
% 'nBurnin'  - num samples for burn-in per chain [1000]
% 'nSamples' - num samples to keep after burn-in [5000]
% 'thin'     - keep every n'th step [1]
% 'blocking' - do block updating for GLMs (using IRLS proposal)?
%              [1 - doesn't yet work]
% 'view'     - set to 1 if you want to view WinBUGS output (then close the WinBUGS
%                window to return control to matlab)
%              set to 0 to close WinBUGS automatically without pausing [default 0]
% 'openBugs' - set to 1 to use openBugs file format [0]
% 'Bugdir'   - location of winbugs executable
%               Default is 'C:/Program Files/WinBUGS14' if not openBugs
%               Default is 'C:/Program Files/OpenBUGS' if OpenBugs.
% 'workingDir' - directory to store temporary data/init/coda files [pwd/tmp]
%
% 'DICstatus' - takes value 1 to set the DIC tool and 0 otherwise
% 'refreshrate' - sets the refresh rate for the updater. Default is 100. Values
%               of 10 prevent the computer from hanging too much if the model





