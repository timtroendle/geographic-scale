function uq(IFName1, IFName2, IFName3, BFName, OutputFolder)
%% SCRIPT CONTROL
rng(100);

CALCULATE = true; % Repeat the entire calculation or not. Warning: will overwrite the existing session file
SFNAME = 'uq_GeoScale_Session_MF.mat'; % Session filename (~500MB)

% Choose which plots you want
REG_PLOTSOBOL = false; % Regional Sobol indices (new outputs)
PLOTSOBOL = false;     % Multi-fidelity indices
PLOTDIAG  = false;     % Several diagnostics plots
PLOTED    = false;    % ED of the LF model (not that useful)
PLOTHIST  = false;     % Histogram of some selected QoI


% Import the previous session if availble, recalculate everything otherwise
if ~CALCULATE
    uqlab(SFNAME);
else
    uqlab
end


%% Import the files
M = 12;
XYLF = dlmread(IFName1,',',1,1);
NEDLF = size(XYLF,1);
XYHF = dlmread(IFName2,',',1,1);
NEDHF = size(XYHF,1);

XYReg = dlmread(IFName3,',',1,1);
NEDReg = size(XYReg,1);

% Get the I/O variable names from the input files
fid = fopen(IFName1);
tt = textscan(fid,'%s',34,'delimiter',',');
fclose(fid);
InVarNames = tt{1}(2:13);
OutVarNames = tt{1}(14:end-1);

% Get them also for the regional outputs
fid = fopen(IFName3);
tt = textscan(fid,'%s',22,'delimiter',',');
fclose(fid);
RegOutVarNames = tt{1}(14:end-1);


% Extract the EDs
% LF
%llLF = load(EDFName1);
%XEDLF = llLF.XED;
XEDLF = XYLF(:,1:M);
YEDLF = XYLF(:,M+1:end-1);

% HF
%llHF = load(EDFName2);
%XEDHF = llHF.XED;
XEDHF = XYHF(:,1:M);
YEDHF = XYHF(:,M+1:end-1);

%% Now extract also the extra outputs for the second analysis
XEDReg = XYReg(:,1:M);
YEDReg = XYReg(:,M+1:end);



%% Get the parameter Bounds (for the inputs)
fid = fopen(BFName);
tt = textscan(fid,'%s%s%s%s%s%s','delimiter',',');
M = size(tt{3}(2:end),1);
Bounds = zeros(2,M);
for ii = 1:M
    Bounds(1,ii) = str2double(tt{3}{ii+1});
    Bounds(2,ii) = str2double(tt{4}{ii+1});
end


%% All inputs are going to be uniform within the provided bounds
for ii = 1:M
    IOpts.Marginals(ii).Name = InVarNames{ii};
    IOpts.Marginals(ii).Type = 'Uniform';
    IOpts.Marginals(ii).Parameters = Bounds(:,ii)';
end

if CALCULATE
    myInput = uq_createInput(IOpts);
else
    myInput = uq_getInput();
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%               REGIONAL DATA ANALYSIS
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Create the PCE for the regional data (simpler, just once)
PCEOpts.Type = 'Metamodel';
PCEOpts.Name = 'Reg PCE';
PCEOpts.MetaType = 'PCE';
PCEOpts.Degree = 1:20;
PCEOpts.TruncOptions.qNorm = 0.55:.05:1;
PCEOpts.TruncOptions.MaxInteraction = 3;
PCEOpts.Bootstrap.Replications = 50;
PCEOpts.Display = 3;

PCEOpts.ExpDesign.X = XEDReg;
PCEOpts.ExpDesign.Y = YEDReg;
if CALCULATE
    % force overwriting in case of repeat trials
    myPCE_Reg = uq_createModel(PCEOpts,'-overwrite');
else
    myPCE_Reg = uq_getModel(PCEOpts.Name);
end

%% Sensitivity analysis on the Reg model
SOpts.Type = 'Sensitivity';
SOpts.Method = 'Sobol';
SOpts.Sobol.Order = 2;
SOpts.Model = myPCE_Reg;
SOpts.Input = myInput;

mySens_Reg = uq_createAnalysis(SOpts);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%               MULTIFIDELITY DATA ANALYSIS
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Calculate the LF PCE

PCEOpts.ExpDesign.X = XEDLF;
PCEOpts.ExpDesign.Y = YEDLF;
if CALCULATE
    % force overwriting in case of repeat trials
    myPCE_LF = uq_createModel(PCEOpts,'-overwrite');
else
    myPCE_LF = uq_getModel(PCEOpts.Name);
end

%% And the PCE of the residual
% Low fidelity version of the HF simulations
YEDHF_LF = uq_evalModel(myPCE_LF,XEDHF);
YEDRes = YEDHF - YEDHF_LF;
PCEOpts.Name = 'Res PCE';
PCEOpts.Degree = 1:5;
PCEOpts.TruncOptions.MaxInteraction = 2;
PCEOpts.Bootstrap.Replications = 50;
% MF uses the X of the HF set, and the difference of their responses
PCEOpts.ExpDesign.X = XEDHF;
PCEOpts.ExpDesign.Y = YEDRes;

if CALCULATE
    myPCE_Res = uq_createModel(PCEOpts,'-overwrite');
else
    myPCE_Res = uq_getModel(PCEOpts.Name);
end


%% Assemble the coefficients of the sparse MF PCE
% Let's be cautious: we want to match the coefficients carefully for each
% of the outputs. Some basis elements may not be common between the two
for ii = 1:15
    % get coefficients, bases, and their indices
    CLF = myPCE_LF.PCE(ii).Coefficients;
    CRes = myPCE_Res.PCE(ii).Coefficients;

    % we only want the non-zero ones
    ALF = myPCE_LF.PCE(ii).Basis.Indices(CLF~=0,:);
    ARes = myPCE_Res.PCE(ii).Basis.Indices(CRes~=0,:);
    CLF = CLF(CLF~=0);
    CRes = CRes(CRes~=0);

    % take their union
    MF(ii).Indices = union(ALF,ARes,'rows');

    % Now loop over the bases
    P = size(MF(ii).Indices,1);

    % Sparse expansions, so there should be no problems with the size of MF
    MF(ii).Coefficients = zeros(P,1);
    for pp = 1:size(MF(ii).Indices,1)
        idx1 = find(ismember(ALF,MF(ii).Indices(pp,:),'rows'));
        idx2 = find(ismember(ARes,MF(ii).Indices(pp,:),'rows'));
        if ~isempty(idx1)
            MF(ii).Coefficients(pp) = MF(ii).Coefficients(pp) + CLF(idx1);
        end
        if ~isempty(idx2)
            MF(ii).Coefficients(pp) = MF(ii).Coefficients(pp) + CRes(idx2);
        end
    end
end

%% Now combine everything in a single MF PCE
MFOpts.Type = 'Metamodel';
MFOpts.MetaType = 'PCE';
MFOpts.Name = 'MF PCE';
MFOpts.Method = 'Custom';
% Copy over most of the data structure from the LF PCE (contains all
% basis definitions and so on, saves some lines of code)
MFOpts.PCE = myPCE_LF.PCE;
% Now overwrite the indices with the actual ones from the MF assembly
for oo = 1:15
    MFOpts.PCE(oo).Basis.Indices = MF(oo).Indices;
    MFOpts.PCE(oo).Coefficients = MF(oo).Coefficients;
end
MFOpts.Input = myInput;
%% And create the MF PCE
if CALCULATE
    MFPCE = uq_createModel(MFOpts,'-overwrite');
else
    MFPCE = uq_getModel(MFOpts.Name);
end


%% Sensitivity analysis on the MF model
SOpts.Type = 'Sensitivity';
SOpts.Method = 'Sobol';
SOpts.Sobol.Order = 2;
SOpts.Model = MFPCE;
SOpts.Input = myInput;

mySens = uq_createAnalysis(SOpts);



%% BONUS: CROSS-VALIDATION ANALYSIS
% idea: let's only take 9/10 HF values at a time, and predict the remaining
% one to get an idea about the actual generalization properties (and to see
% if the LOO makes sense indeed

YPCECV = zeros(size(YEDHF));
for cc = 1:NEDHF
    idx = true(NEDHF,1);
    idx(cc) = false; % exclude the ith data point
    PCEOpts.Name = sprintf('CV PCE %i',cc);
    PCEOpts.ExpDesign.X = XEDHF(idx,:);
    PCEOpts.ExpDesign.Y = YEDRes(idx,:);
    if CALCULATE
        CVPCE = uq_createModel(PCEOpts,'-overwrite');
    else
        CVPCE = uq_getModel(PCEOpts.Name);
    end
    YPCECV(~idx,:) = uq_evalModel(CVPCE, XEDHF(~idx,:))+uq_evalModel(myPCE_LF,XEDHF(~idx,:));
end

%% Write results

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%               WRITE RESULTS
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Sobol indices

dlmwrite(strcat(OutputFolder, 'total-sobol-sf.csv'),mySens_Reg.Results.Total,'delimiter',',','precision','%.16e');
dlmwrite(strcat(OutputFolder, 'first-sobol-sf.csv'),mySens_Reg.Results.FirstOrder,'delimiter',',','precision','%.16e');
dlmwrite(strcat(OutputFolder, 'total-minus-first-sobol-sf.csv'),mySens_Reg.Results.Total - mySens_Reg.Results.FirstOrder,'delimiter',',','precision','%.16e');

dlmwrite(strcat(OutputFolder, 'total-sobol-mf.csv'),mySens.Results.Total,'delimiter',',','precision','%.16e');
dlmwrite(strcat(OutputFolder, 'first-sobol-mf.csv'),mySens.Results.FirstOrder,'delimiter',',','precision','%.16e');
dlmwrite(strcat(OutputFolder, 'total-minus-first-sobol-mf.csv'),mySens.Results.Total - mySens.Results.FirstOrder,'delimiter',',','precision','%.16e');

% Samples

XX = uq_getSample(myInput, 1e5);
YY = uq_evalModel(myPCE_Reg, XX);
dlmwrite(strcat(OutputFolder, 'pce-samples-sf.csv'),YY,'delimiter',',','precision','%.16e');

YY = uq_evalModel(MFPCE, XX);
dlmwrite(strcat(OutputFolder, 'pce-samples-mf.csv'),YY,'delimiter',',','precision','%.16e');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%               PLOTTING
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%% CREATE SOME PLOTS FOR DIAGNOSTICS AND PRELIMINARY ANALYSIS
%  get the output names
NOut = length(OutVarNames);


%% PLOT THE EXPERIMENTAL DESIGN (RAW DATA)
if PLOTED
    uq_figure('FileName', 'Figures/20190719_ED.fig');
    for ii = 1:NOut Y_i{ii} = sprintf('Y_{%i}',ii); end
    gplotmatrix(XEDLF, YEDLF,[],'k','o',3,[],[],InVarNames,Y_i);
end

if PLOTDIAG
    %% PLOT THE EFFECT OF THE MF SURROGATE
    uq_figure; uq_plot(YEDHF(:),YEDHF_LF(:),'+r')
    xl = xlim;
    axis equal
    hold on; uq_plot(YEDHF(:),YPCECV(:),'og')
    xlim(xl);
    ylim(xl);
    hold on; uq_plot(xl,xl,'--k','linewidth',1);
    set(gca,'XScale', 'log','YScale', 'log');

    for ii = 1:NOut
        uq_figure;
        histogram((YEDHF_LF(:,ii)-YEDHF(:,ii))/std(YEDHF_LF(:,ii)),10);
        hold on
        histogram((YPCECV(:,ii) - YEDHF(:,ii))/std(YEDHF_LF(:,ii)),10);
        legend('LF-HF', 'MF - HF')
        xl = xlim;
        xlim([-max(abs(xl)) max(abs(xl))])
    end

    %% MF QUALITY ESTIMATION
    SelectedVar = 1;
    figure; plot(YEDHF(:,SelectedVar),YEDHF_LF(:,SelectedVar),'+r')
    hold on; plot(YEDHF(:,SelectedVar),YPCECV(:,SelectedVar),'og')
    axis tight equal
    xl = xlim;
    xlim(xl),ylim(xl)
    hold on; plot(xl ,xl,'--k')
end

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%               SOBOL INDICES REGIONAL
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if REG_PLOTSOBOL
    %% PLOT TOTAL INDICES (REGIONAL)
    uq_figure('filename', 'Figures/20190719_Sobol_Total.fig');
    imagesc(mySens_Reg.Results.Total');
    caxis([0 1]); colormap summer; colorbar;axis tight;
    set(gca,'ytick', 1:NOut, 'yticklabel',RegOutVarNames,'YTickLabelRotation',25)
    set(gca,'xtick', 1:M, 'xticklabel',mySens_Reg.Results.VariableNames,'XTickLabelRotation',65)
    set(gca,'fontsize',14)
    title('Total Sobol'' indices', 'FontWeight','normal')
    ss = bone(512);
    colormap(ss(end:-1:1,:))

    %% PLOT FIRST ORDER INDICES
    uq_figure('filename', 'Figures/20190719_Sobol_FirstOrder.fig');
    imagesc(mySens_Reg.Results.FirstOrder');
    caxis([0 1]); colormap summer; colorbar;axis tight;
    set(gca,'ytick', 1:NOut, 'yticklabel',RegOutVarNames,'YTickLabelRotation',25)
    set(gca,'xtick', 1:M, 'xticklabel',mySens_Reg.Results.VariableNames,'XTickLabelRotation',65)
    set(gca,'fontsize',14)
    title('First Order Sobol'' indices', 'FontWeight','normal')
    ss = bone(512);
    colormap(ss(end:-1:1,:))

    %% PLOT NON-LINEAR EFFECTS (Tot - first order)
    uq_figure('filename', 'Figures/20190719_Sobol_TotFirstDifference.fig');
    imagesc(mySens_Reg.Results.Total' - mySens_Reg.Results.FirstOrder');
    caxis([0 1]); colormap summer; colorbar; axis tight;
    set(gca,'ytick', 1:NOut, 'yticklabel',RegOutVarNames,'YTickLabelRotation',25)
    set(gca,'xtick', 1:M, 'xticklabel',mySens_Reg.Results.VariableNames,'XTickLabelRotation',65)
    set(gca,'fontsize',14)
    title('Total - First Order indices', 'FontWeight','normal')
    ss = bone(512);
    colormap(ss(end:-1:1,:))
end

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%               SOBOL INDICES MULTIFIDELITY
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if PLOTSOBOL
    %% PLOT TOTAL INDICES (MULTIFIDELITY)
    uq_figure('filename', 'Figures/20190719_Sobol_Total.fig');
    imagesc(mySens.Results.Total');
    caxis([0 1]); colormap summer; colorbar;axis tight;
    set(gca,'ytick', 1:NOut, 'yticklabel',OutVarNames,'YTickLabelRotation',25)
    set(gca,'xtick', 1:M, 'xticklabel',mySens.Results.VariableNames,'XTickLabelRotation',65)
    set(gca,'fontsize',14)
    title('Total Sobol'' indices', 'FontWeight','normal')
    ss = bone(512);
    colormap(ss(end:-1:1,:))

    %% PLOT FIRST ORDER INDICES
    uq_figure('filename', 'Figures/20190719_Sobol_FirstOrder.fig');
    imagesc(mySens.Results.FirstOrder');
    caxis([0 1]); colormap summer; colorbar;axis tight;
    set(gca,'ytick', 1:NOut, 'yticklabel',OutVarNames,'YTickLabelRotation',25)
    set(gca,'xtick', 1:M, 'xticklabel',mySens.Results.VariableNames,'XTickLabelRotation',65)
    set(gca,'fontsize',14)
    title('First Order Sobol'' indices', 'FontWeight','normal')
    ss = bone(512);
    colormap(ss(end:-1:1,:))

    %% PLOT NON-LINEAR EFFECTS (Tot - first order)
    uq_figure('filename', 'Figures/20190719_Sobol_TotFirstDifference.fig');
    imagesc(mySens.Results.Total' - mySens.Results.FirstOrder');
    caxis([0 1]); colormap summer; colorbar; axis tight;
    set(gca,'ytick', 1:NOut, 'yticklabel',OutVarNames,'YTickLabelRotation',25)
    set(gca,'xtick', 1:M, 'xticklabel',mySens.Results.VariableNames,'XTickLabelRotation',65)
    set(gca,'fontsize',14)
    title('Total - First Order indices', 'FontWeight','normal')
    ss = bone(512);
    colormap(ss(end:-1:1,:))

    %% PLOTTING THE SOBOL' INDICES AS BARS FOR DIFFERENT OUTPUTS
    % Let's say we want to see outputs 4 and 13, then:
    SelectedVars = [3 4];
    uq_figure('FileName','Figures/20190719_Sobol_SideBySide')
    uq_bar([mySens.Results.Total(:,SelectedVars)])
    legend(OutVarNames(SelectedVars))
    title('Total Comparison')

    uq_figure('FileName','Figures/20190719_Sobol_SideBySide')
    uq_bar([mySens.Results.FirstOrder(:,SelectedVars)])
    legend(OutVarNames(SelectedVars))
    title('First Order Comparison')
    legend(OutVarNames(SelectedVars))
end


%% PLOTTING A HISTOGRAM OF THE PCE MODEL
% Make a histogram out of a million resampled values
XX = uq_getSample(1e6);
YY = uq_evalModel(MFPCE,XX); % MF model
YYReg = uq_evalModel(myPCE_Reg,XX); % Reg model

SelectedVar = 4;
SelectedVarReg = 4;
if PLOTHIST
    % Now make an histogram, possibly normalized as a pdf
    uq_figure('FileName','Figures/20190719_Histo');
    histogram(YY(:,SelectedVar),'Normalization','pdf','EdgeColor','none','FaceColor',[1 0 1]) % because fuchsia
    xlabel(OutVarNames{SelectedVar})
    uq_figure('FileName','Figures/20190719_Histo');
    histogram(YYReg(:,SelectedVarReg),'Normalization','pdf','EdgeColor','none','FaceColor',[1 0 1]) % because fuchsia
    xlabel(RegOutVarNames{SelectedVarReg})
end
end
