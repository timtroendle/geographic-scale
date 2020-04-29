function experimental_design(IFName, OFName, NED)
%% Initialize the input
uqlab
rng(100);
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
GeoScale_Bounds = readtable(IFName, opts);

clear opts

%% All inputs are going to be uniform within the provided bounds
M = size(GeoScale_Bounds,1);
for ii = 1:M
    IOpts.Marginals(ii).Name = char(GeoScale_Bounds.parameternameshort(ii));
    IOpts.Marginals(ii).Type = 'Uniform';
    IOpts.Marginals(ii).Parameters = [GeoScale_Bounds.min(ii) GeoScale_Bounds.max(ii)];
end

myInput = uq_createInput(IOpts);

%% Generate an highly optimized LHS
tic;XED = uq_getSample(NED,'LHS','iterations',10000);toc

%% Write out the results
% CSV file
fout = fopen(OFName,'w+t');
% Write the header first
for mm = 1:M
    fprintf(fout,'"%s"',GeoScale_Bounds.parameternameshort(mm));
    if mm<M
        fprintf(fout,',');
    else
        fprintf(fout,'\n');
    end
end

% And the data
for ii = 1:NED
    for mm = 1:M
    fprintf(fout,'%.16e',XED(ii,mm));
    if mm<M
        fprintf(fout,',');
    else
        fprintf(fout,'\n');
    end
    end
end

fclose(fout);
end
