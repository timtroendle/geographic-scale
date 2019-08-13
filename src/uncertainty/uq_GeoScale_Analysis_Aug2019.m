%% Initialize the input
clearvars
uqlab
NED = 150;
rng(100);
OFName = 'GeoScaleED_20190719';
IFName = 'xy_20190719.csv';
%% Import the file (automatic Matlab generated code)
opts = delimitedTextImportOptions("NumVariables", 6);

% Specify range and delimiter
opts.DataLines = [2, Inf];
opts.Delimiter = ",";

% Specify column names and types
opts.VariableNames = ["Var1", "parameternameshort", "min", "max", "Var5", "Var6"];
opts.SelectedVariableNames = ["parameternameshort", "min", "max"];
opts.VariableTypes = ["string", "string", "double", "double", "string", "string"];
opts = setvaropts(opts, [1, 2, 5, 6], "WhitespaceRule", "preserve");
opts = setvaropts(opts, [1, 2, 5, 6], "EmptyFieldRule", "auto");
opts.ExtraColumnsRule = "ignore";
opts.EmptyLineRule = "read";

% Import the data
GeoScale_Bounds = readtable("uncertain-parameters.csv", opts);

clear opts

%% All inputs are going to be uniform within the provided bounds
M = size(GeoScale_Bounds,1);
for ii = 1:M
    IOpts.Marginals(ii).Name = char(GeoScale_Bounds.parameternameshort(ii));
    IOpts.Marginals(ii).Type = 'Uniform';
    IOpts.Marginals(ii).Parameters = [GeoScale_Bounds.min(ii) GeoScale_Bounds.max(ii)];
end

myInput = uq_createInput(IOpts);

%% Read the ED
opts = delimitedTextImportOptions("NumVariables", 26);

opts.DataLines = [2, Inf];
opts.Delimiter = ",";

% Specify column names and types
opts.VariableNames = ["Var1", "Var2", "Var3", "Var4", "Var5", "Var6", "Var7", "Var8", "Var9", "Var10", "Var11", "Var12", "Var13", "ylargescalecosteur", "ysmallscalecosteur", "ycostdiffeur", "ycostdiffrelative", "ylargescalepvgw", "ysmallscalepvgw", "ylargescalewindgw", "ysmallscalewindgw", "ylargescalebiofuelgw", "ysmallscalebiofuelgw", "ylargescalestoragegw", "ysmallscalestoragegw", "ylargescalestoragegwh", "ysmallscalestoragegwh", "ylargescaletransmissiongwkm", "Var29"];
opts.SelectedVariableNames = ["ylargescalecosteur", "ysmallscalecosteur", "ycostdiffeur", "ycostdiffrelative", "ylargescalepvgw", "ysmallscalepvgw", "ylargescalewindgw", "ysmallscalewindgw", "ylargescalebiofuelgw", "ysmallscalebiofuelgw", "ylargescalestoragegw", "ysmallscalestoragegw", "ylargescalestoragegwh", "ysmallscalestoragegwh", "ylargescaletransmissiongwkm"];
opts.VariableTypes = ["string", "string", "string", "string", "string", "string", "string", "string", "string", "string", "string", "string", "string", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "string"];
opts = setvaropts(opts, [1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 29], "WhitespaceRule", "preserve");
opts = setvaropts(opts, [1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 29], "EmptyFieldRule", "auto");
opts.ExtraColumnsRule = "ignore";
opts.EmptyLineRule = "read";

% Import the data
xy = readtable(IFName, opts);

YED = table2array(xy);
% Clear temporary variables
clear opts

% Now let's get the experimental design x from the output that I sent them
% in the first place
load GeoScaleED_20190719;


%% Now create the PCE

PCEOpts.Type = 'Metamodel';
PCEOpts.MetaType = 'PCE';
PCEOpts.Degree = 1:20;
PCEOpts.TruncOptions.qNorm = 0.55:.05:1;
PCEOpts.TruncOptions.MaxInteraction = 3;
PCEOpts.Bootstrap.Replications = 50;
PCEOpts.ExpDesign.X = XED;
PCEOpts.ExpDesign.Y = YED;
PCEOpts.Display = 3;

myPCE = uq_createModel(PCEOpts);
uq_print(myPCE);

%% Sensitivity analysis
SOpts.Type = 'Sensitivity';
SOpts.Method = 'Sobol';
SOpts.Sobol.Order = 2;
SOpts.Model = myPCE;
SOpts.Input = myInput;
mySens = uq_createAnalysis(SOpts);

%% CREATE SOME PLOTS FOR DIAGNOSTICS AND PRELIMINARY ANALYSIS
%  get the output names
OutNames = xy.Properties.VariableNames;
NOut = length(OutNames);

%% PLOT THE EXPERIMENTAL DESIGN (RAW DATA)
uq_figure('FileName', 'Figures/20190719_ED.fig'); 
VarNames = {myInput.Marginals.Name};
for ii = 1:NOut Y_i{ii} = sprintf('Y_{%i}',ii); end
gplotmatrix(XED, YED,[],'k','o',3,[],[],VarNames,Y_i);

%% PLOT THE PCE CONVERGENCE (LEAVE-ONE-OUT ERROR)
LOOs = [myPCE.Error(:).LOO];
uq_figure('Filename','Figures/20190719_PCE_loo.fig','Name','20190719_PCE_loo'); 
uq_plot(LOOs,NOut:-1:1,'ok')
xlabel('$\epsilon_{LOO}$')
set(gca,'ytick',1:NOut,'yticklabel',OutNames(end:-1:1),'xscale','log','xtick',[1e-2 1e-1])
ylim([0,NOut+1])
xlim([1e-3 2e-1])
title('PCE accuracy')
hold on; plot([5e-2, 5e-2],[0 NOut+1],'--g',[1e-1, 1e-1],[0 NOut+1],'--r','linewidth',2)

%% PLOT TOTAL INDICES
uq_figure('filename', 'Figures/20190719_Sobol_Total.fig'); 
imagesc(mySens.Results.Total');
caxis([0 1]); colormap summer; colorbar;axis tight;
set(gca,'ytick', 1:NOut, 'yticklabel',OutNames,'YTickLabelRotation',25)
set(gca,'xtick', 1:M, 'xticklabel',mySens.Results.VariableNames,'XTickLabelRotation',65)
set(gca,'fontsize',14)
title('Total Sobol'' indices', 'FontWeight','normal')

%% PLOT FIRST ORDER INDICES
uq_figure('filename', 'Figures/20190719_Sobol_FirstOrder.fig'); 
imagesc(mySens.Results.FirstOrder');
caxis([0 1]); colormap summer; colorbar;axis tight;
set(gca,'ytick', 1:NOut, 'yticklabel',OutNames,'YTickLabelRotation',25)
set(gca,'xtick', 1:M, 'xticklabel',mySens.Results.VariableNames,'XTickLabelRotation',65)
set(gca,'fontsize',14)
title('First Order Sobol'' indices', 'FontWeight','normal')

%% PLOT NON-LINEAR EFFECTS (Tot - first order)
uq_figure('filename', 'Figures/20190719_Sobol_TotFirstDifference.fig'); 
imagesc(mySens.Results.Total' - mySens.Results.FirstOrder');
caxis([0 1]); colormap summer; colorbar; axis tight;
set(gca,'ytick', 1:NOut, 'yticklabel',OutNames,'YTickLabelRotation',25)
set(gca,'xtick', 1:M, 'xticklabel',mySens.Results.VariableNames,'XTickLabelRotation',65)
set(gca,'fontsize',14)
title('Total - First Order indices', 'FontWeight','normal')

%% PLOT UNIVARIATE EFFECTS
% syntax: uq_PCE_displayEE(myPCE,1,NOUT), where NOUT is the output number.
% Example: uq_PCE_displayEE(myPCE,1,4) will plot the univariate effects
% for the relative price difference.
uq_PCE_displayEE(myPCE,1,4)



%% PLOTTING THE SOBOL' INDICES AS BARS FOR DIFFERENT OUTPUTS
% Let's say we want to see outputs 4 and 13, then:
SelectedVars = [4 13];
uq_figure('FileName','Figures/20190719_Sobol_SideBySide')
uq_bar(mySens.Results.Total(:,SelectedVars))
legend(OutNames(SelectedVars))


%% PLOTTING A HISTOGRAM OF THE PCE MODEL
% Make a histogram out of a million resampled values
XX = uq_getSample(1e6);
YY = uq_evalModel(XX);

% Now make an histogram, possibly normalized as a pdf
uq_figure('FileName','Figures/20190719_Histo'); 
histogram(YY(:,4),'Normalization','pdf','EdgeColor','none','FaceColor',[1 0 1]) % because fuchsia
xlabel(OutNames{4})