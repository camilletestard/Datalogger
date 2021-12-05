%  AUTHORSHIP
%   Primary Developer: Stephen Meehan <swmeehan@stanford.edu> 
%   Math Lead & Secondary Developer:  Connor Meehan <connor.gw.meehan@gmail.com>
%   Bioinformatics Lead:  Wayne Moore <wmoore@stanford.edu>
%   Provided by the Herzenberg Lab at Stanford University 
%   License: BSD 3 clause
%

function [ umapOrEppOrMatch, reduction, clusterIdentifiers, extras]=...
    suh_pipelines(csvFileOrData, varargin)
%Brokers ALL input and output argument of either run_umap or run_epp
%or SuhMatch.Run depending on the named argument 'pipeline' which 
%defaults to 'umap'.
%For usage of umap see umap/run_umap.m
%For usage of epp  see epp/run_epp.m

if nargin<1
    csvFileOrData=[];
end

if ~isdeployed
    try
        edu.stanford.facs.swing.CpuInfo.isMac;
    catch
        initPaths
        try
            isMac=edu.stanford.facs.swing.CpuInfo.isMac;
        catch ex
            msg('Problem loading SUH jars');
            return;
        end
    end
end
if isempty(csvFileOrData)
    curPath=fileparts(mfilename('fullpath'));
    if ~isdeployed
        mexNN=fullfile(curPath, 'umap', UmapUtil.LocateMex);
        mexEPP=fullfile(curPath, 'epp', SuhEpp.LocateMex);
        if ~exist(mexEPP, 'file')
            mexEPP=fullfile(curPath, 'util', SuhEpp.LocateMex);
        end
        if ~exist(mexEPP, 'file') || ~exist(mexNN, 'file')
            UmapUtil.OfferFullDistribution(true)
            BasicMap.Global.save;
            mexEPP=fullfile(curPath, 'epp', SuhEpp.LocateMex);
            if ~exist(mexNN, 'file')
                msg('Must have UMAP''s MEX files to continue');
                return;
            end
            if ~exist(mexEPP, 'file')  
                if ~askYesOrNo('Proceed without EPP pipeline?')
                    return;
                end
            end
        end
    end
    ArgumentClinic;
    return;
end

args = inputParser;
addParameter(args, 'pipeline', 'epp', @ischar);
args=Args.NewKeepUnmatched(args, varargin{:});
varArgIn=Args.RemoveArg(varargin, 'pipeline');
umapOrEppOrMatch=[];
reduction=[];
clusterIdentifiers=[];
extras=[];
MatBasics.WarningsOff;
if ~strcmpi(args.pipeline, 'epp') && ~strcmpi(args.pipeline, 'umap')
    msg('first arg must be ''epp'',  ''umap'' or ''match''');
    return;
end
if strcmpi(args.pipeline, 'umap')
    if isdeployed
        if length(varArgIn)>1
            umapArgs=Args(UmapUtil.DefineArgs);
            args=Args.Str2NumOrLogical(umapArgs.p.Results, varArgIn);
            [reduction, umapOrEppOrMatch, clusterIdentifiers, extras]=...
                run_umap(csvFileOrData, args{:});
        else
            [reduction, umapOrEppOrMatch, clusterIdentifiers, extras]=...
                run_umap(csvFileOrData);
        end
    else
        [reduction, umapOrEppOrMatch, clusterIdentifiers, extras]=...
            run_umap(csvFileOrData, varArgIn{:});
    end
elseif strcmpi(args.pipeline, 'epp')
    umapOrEppOrMatch=run_epp(csvFileOrData, varArgIn{:});
    if isempty(umapOrEppOrMatch)
        fprintf('\n\nThe EPP hierarchy was NOT built!!\n\n');
        return;
    end
    
else
    [match, matchTable, trainingQfTreeMatch,  trainingQfTree, ...
        testQfTreeMatch, testQfTree]=SuhMatch.Run('training_set', ...
        csvFileOrData, varargin{:});
    umapOrEppOrMatch={match, matchTable, trainingQfTreeMatch,   ...
        trainingQfTree, testQfTreeMatch, testQfTree};
end

    function initPaths
        pPth=fileparts(mfilename('fullpath'));
        utilPath=fullfile(pPth, 'util');
        addpath(utilPath);
        MatBasics.WarningsOff
        if ~initJava
            error('Can not find suh.jar');
        end
        eppPath=fullfile(pPth, 'epp');
        addpath(eppPath);
        addpath(utilPath);
        umapPath=fullfile(pPth, 'umap');
        FileBasics.AddNonConflictingPaths({eppPath, utilPath, umapPath});
    end
end
