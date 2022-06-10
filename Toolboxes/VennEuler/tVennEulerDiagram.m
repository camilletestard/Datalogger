classdef tVennEulerDiagram < matlab.unittest.TestCase
    % tVennEulerDiagram Tests for vennEulerDiagram.
    %   Copyright 2021 The MathWorks, Inc.

    properties
        Figure
        ChartConstructor = @hTestVennEulerDiagram
    end

    properties(TestParameter)
        % Data for testing if correct errors are produced
        testErrorData = createTestErrorData;

        % Data for testing SetMembershipData
        testSMDData = createTestSMDData;

        % Data for testing SetCounts
        testSetCountsData = createTestSetCountsData;

        % Data for testing vennEulerDiagram.computeCircleAreas
        testCircleAreasData = createTestCircleAreasData;

        % Data for testing vennEulerDiagram.computeInitialCenters
        testInitialCentersData = createTestInitialCentersData;

        % Data for testing vennEulerDiagram.findDisjointIntersectionAreas
        testFindDIAreasData = createTestFindDIAreasData;

        % Data for testing vennEulerDiagram.computeStress
        testStressData = createTestStressData;

        % Data for testing vennEulerDiagram.computeGradients
        testComputeGradientsData = createTestComputeGradientsData;

        % Data for testing valid syntaxes
        testSyntaxesData = createTestSyntaxesData
    end
    
    methods (TestClassSetup)
        function createFigure(testCase)
            testCase.Figure = figure;
            testCase.addTeardown(@() close(testCase.Figure))
        end
    end
    
    methods (TestMethodTeardown)
        function clearFigure(testCase)
            clf(testCase.Figure);
        end
    end
    
    methods (Test)

        function testSyntaxes(testCase, testSyntaxesData)
            % Test that valid syntaxes correctly set properties
            args = testSyntaxesData.Args;
            h = hTestVennEulerDiagram(args{:});
            h.callUpdate();

            setListData = testSyntaxesData.SetListData;
            setMembershipData = testSyntaxesData.SetMembershipData;
            setLabels = testSyntaxesData.SetLabels;

            testCase.verifyEqual(h.SetListData, setListData);
            testCase.verifyEqual(h.SetMembershipData, setMembershipData);

            testCase.verifyEqual(h.SetLabels, setLabels);
        end

        function testErrors(testCase, testErrorData)
            % Test that invalid input produces the correct errors
            
            setListData = testErrorData.SetListData;
            propertyName = testErrorData.PropertyName;
            propertyValue = testErrorData.PropertyValue;
            errorID = testErrorData.ErrorID;

            % Error if incorrect number of set labels
            h = hTestVennEulerDiagram(setListData, propertyName, propertyValue);
            testCase.verifyError(@()h.callUpdate(), errorID);
        end
      
        function testSetMembershipData(testCase, testSMDData)
            % Test whether SetMembershipData is correctly set when using 
            % SetListData input
            expSetMembershipData = testSMDData.SetMembershipData;
            actSetMembershipData = hTestVennEulerDiagram.createSetMembershipData(testSMDData.SetListData);
            
            testCase.verifyEqual(actSetMembershipData, expSetMembershipData, ...
                'setMembershipData');
        end

        function testSetCounts(testCase, testSetCountsData)
            % Test whether SetCounts is correctly set when using SetListData
            % input
            expSetCounts = testSetCountsData.SetCounts;
            actSetCounts = hTestVennEulerDiagram.createSetCounts(testSetCountsData.SetListData);
            
            testCase.verifyEqual(actSetCounts, expSetCounts, 'setCounts');
        end

        function testCircleAreas(testCase, testCircleAreasData)
            % Test vennEulerDiagram.computeCircleAreas
            expCircleAreas = testCircleAreasData.CircleAreas;
            actCircleAreas = hTestVennEulerDiagram.callComputeCircleAreas(...
                testCircleAreasData.NumSets, ...
                testCircleAreasData.SetCounts);

            testCase.verifyEqual(actCircleAreas, expCircleAreas, 'circleAreas');
        end

        function testComputeInitialCenters(testCase, testInitialCentersData)
            % Test vennEulerDiagram.computeInitialCenters
            actCircleCenters = hTestVennEulerDiagram.callComputeInitialCircleCenters(...
                testInitialCentersData.SetCounts, ...
                testInitialCentersData.NumSets, ...
                testInitialCentersData.CircleAreas);

            actPairwiseDistances = zeros(1, testInitialCentersData.NumSets);

            for i = 1:(testInitialCentersData.NumSets - 1)
                actPairwiseDistances(i) = norm(actCircleCenters(i, :) - actCircleCenters(i + 1, :));
            end

            actPairwiseDistances(end) = norm(actCircleCenters(end, :) - actCircleCenters(1, :));

            % Verify that the center coordinates are roughly equidistant
            minPairwiseDistances = min(actPairwiseDistances) * ones(size(actPairwiseDistances));
            testCase.verifyEqual(actPairwiseDistances, minPairwiseDistances, 'AbsTol', 0.05, 'initialCircleCenters');
        end

        function testComputeCircleDiameters(testCase)
            % Test vennEulerDiagram.computeCircleDiameters, i,e, the final
            % total area should be 1/n of the old total area
            testAreas = 1:10;
            expSumScaledAreas = sum(testAreas) / numel(testAreas);
            scaledDiameters = hTestVennEulerDiagram.callComputeCircleDiameters(testAreas);
            actSumScaledAreas = sum(pi * (scaledDiameters / 2).^2);

            testCase.verifyEqual(actSumScaledAreas, expSumScaledAreas, 'computeCircleDiameters');
        end

        function testFindDisjointIntersectionAreas(testCase, testFindDIAreasData)
            % Test vennEulerDiagram.findDisjointIntersectionAreas
            expIntersectionAreas = testFindDIAreasData.IntersectionAreas;
            expShownIntersectionIndices = testFindDIAreasData.ShownIntersectionIndices;

            [actIntersectionAreas, ~, actShownIntersectionIndices] = ...
                hTestVennEulerDiagram.callFindDisjointIntersectionAreas(testFindDIAreasData.NumSets, ...
                    testFindDIAreasData.CirclePolyshapes); 

            testCase.verifyEqual(actIntersectionAreas, expIntersectionAreas, 'AbsTol', 0.05, 'intersectionAreas');
            testCase.verifyEqual(actShownIntersectionIndices, expShownIntersectionIndices, 'shownIntersectionIndices');
        end

        function testComputeStress(testCase, testStressData)
            % Test vennEulerDiagram.computeStress. Test for cases of high
            % stress/low stress.
            expStressUpperBound = testStressData.StressUpperBound;
            expStressLowerBound = testStressData.StressLowerBound;
            actStress = hTestVennEulerDiagram.callComputeStress(testStressData.IntersectionAreas, ...
                testStressData.SetCounts);

            % Verify that the stress is between the expected bounds of
            % stress
            testCase.verifyGreaterThanOrEqual(actStress, expStressLowerBound, 'lowerBoundStress');
            testCase.verifyLessThanOrEqual(actStress, expStressUpperBound, 'upperBoundStress');
        end
 
        function testComputeGradients(testCase, testComputeGradientsData)
            % Test vennEulerDiagram.computeGradients
            expGradientDirection = testComputeGradientsData.GradientDirection;

            areas_hat = hTestVennEulerDiagram.computeAreasHat(testComputeGradientsData.IntersectionAreas, testComputeGradientsData.SetCounts);
            actGradient = hTestVennEulerDiagram.callComputeGradients(testComputeGradientsData.NumSets, ...
                testComputeGradientsData.CurrCircleCenters, testComputeGradientsData.IntersectionAreas, areas_hat);

            % Verify that the direction of the gradient is approximately
            % correct by checking that the gradient has the expected sign
            % in each component
            testCase.verifyGreaterThanOrEqual(actGradient .* expGradientDirection, 0, 'gradientDirection');
        end

    end
end

function testErrorData = createTestErrorData
        setListData = {[1:3] [4:6] [7]}; %#ok<NBRAK> 
        numSets = 3;

        % Error if incorrect number of set labels
        setLabels = ["A", "B"];

        % Error if incorrect number of intersection colors
        intersectionColors = 0.5 * ones(2^numSets - 2, 3);

        % Error if incorrect number of intersection transparencies
        intersectionTransparencies = 0.5 * ones(2^numSets - 2, 1);
        
        % Error if incorrect number of circle face colors
        circleFaceColors = 0.5 * ones(numSets - 1, 3);

        % Error if incorrect number of circle face transparencies
        circleFaceTransparencies = ones(numSets - 1, 1);

        % Error if incorrect number of circle edge colors
        circleEdgeColors = ones(numSets - 1, 3);

        % Error if incorrect number of circle edge widths
        circleEdgeWidths = ones(numSets - 1, 1);

        testErrorData = struct( ...
            'SetLabels', struct('SetListData', setListData, ...
                'PropertyName', 'SetLabels', 'PropertyValue', setLabels, 'ErrorID', 'incorrectSize:SetLabels'), ...
            'IntersectionColors', struct('SetListData', setListData, ...
                'PropertyName', 'IntersectionColors', 'PropertyValue', intersectionColors, 'ErrorID', 'incorrectSize:IntersectionColors'), ...
            'IntersectionTransparencies', struct('SetListData', setListData, ...
                'PropertyName', 'IntersectionTransparencies', 'PropertyValue', intersectionTransparencies, 'ErrorID', 'incorrectSize:IntersectionTransparencies'), ...
            'CircleFaceColors', struct('SetListData', setListData, ...
                'PropertyName', 'CircleFaceColors', 'PropertyValue', circleFaceColors, 'ErrorID', 'incorrectSize:CircleFaceColors'), ...
            'CircleFaceTransparencies', struct('SetListData', setListData, ...
                'PropertyName', 'CircleFaceTransparencies', 'PropertyValue', circleFaceTransparencies, 'ErrorID', 'incorrectSize:CircleFaceTransparencies'), ...
            'CircleEdgeColors', struct('SetListData', setListData, ...
                'PropertyName', 'CircleEdgeColors', 'PropertyValue', circleEdgeColors, 'ErrorID', 'incorrectSize:CircleEdgeColors'), ...
            'CircleEdgeWidths', struct('SetListData', setListData, ...
                'PropertyName', 'CircleEdgeWidths', 'PropertyValue', circleEdgeWidths, 'ErrorID', 'incorrectSize:CircleEdgeWidths'));
end

function testSMDData = createTestSMDData
    % Create data for testing SetMembershipData is correct

    fewSetsSetListData = {[1 2 3] [1 2] [3 4 5] [1 5]};
    manySetsSetListData = {[1:3] [4:6, 3] [1:5] [6:9, 3] [10:12, 3:4] [4]}; %#ok<NBRAK> 
    duplicateSetsSetListData = {[1:3] [1:3] [3:5] [3:5]}; %#ok<NBRAK> 

    testSMDData = struct( ...
        'FewSets', struct('SetListData', {fewSetsSetListData}, ...
            'SetMembershipData', [1 1 0 1; 
                                  1 1 0 0; 
                                  1 0 1 0; 
                                  0 0 1 0; 
                                  0 0 1 1]), ...
        'ManySets', struct('SetListData', {manySetsSetListData}, ...
            'SetMembershipData', [1 0 1 0 0 0;
                                  1 0 1 0 0 0;
                                  1 1 1 1 1 0;
                                  0 1 1 0 1 1;
                                  0 1 1 0 0 0;
                                  0 1 0 1 0 0;
                                  0 0 0 1 0 0;
                                  0 0 0 1 0 0;
                                  0 0 0 1 0 0;
                                  0 0 0 0 1 0;
                                  0 0 0 0 1 0;
                                  0 0 0 0 1 0]), ...
        'DuplicateSets', struct('SetListData', {duplicateSetsSetListData}, ...
            'SetMembershipData', [1 1 0 0;
                                  1 1 0 0;
                                  1 1 1 1;
                                  0 0 1 1;
                                  0 0 1 1]));
end

function testCountsData = createTestSetCountsData
    % Create data for testing SetCounts is correct

    fewSetsSetListData = {[1 2 3] [1 2] [3 4 5] [1 5]};
    manySetsSetListData = {[1:3] [4:6, 3] [1:5] [6:9, 3] [10:12, 3:4] [4]}; %#ok<NBRAK> 
    duplicateSetsSetListData = {[1:3] [1:3] [3:5] [3:5]}; %#ok<NBRAK> 

    % SetCounts for each of the cell arrays of sets
    fewSetsCounts = zeros(2^(numel(fewSetsSetListData)), 1);
    manySetsCounts = zeros(2^(numel(manySetsSetListData)), 1);
    duplicateSetsCounts = zeros(4, 1);

    fewSetsCounts([4 5 6 12 13]) = 1;
    duplicateSetsCounts(4) = 1;
    duplicateSetsCounts([2 3]) = 2;
    manySetsCounts([7 11 32 55]) = 1;
    manySetsCounts(6) = 2;
    manySetsCounts([9 17]) = 3;

    testCountsData = struct(...
        'FewSets', struct('SetListData', {fewSetsSetListData}, ...
            'SetCounts', fewSetsCounts), ...
        'ManySets', struct('SetListData', {manySetsSetListData}, ...
            'SetCounts', manySetsCounts), ...
        'DuplicateSets', struct('SetListData', {duplicateSetsSetListData}, ...
            'SetCounts', duplicateSetsCounts));
end

function testCircleAreasData = createTestCircleAreasData
    % Create data for testing vennEulerDiagram.computeCircleAreas
    
    % Circle areas and counts for equal area circles
    EA_NumSets = 5;
    EA_CircleAreas = EA_NumSets * ones(EA_NumSets, 1);
    EA_SetCounts = zeros(2^EA_NumSets, 1);
    EA_SetCounts([4 5 6 7 8 29]) = 1;
    EA_SetCounts([2 3 9 17 25]) = 2;

    % Circle areas and counts for circles which do not intersect
    NI_NumSets = 5;
    NI_CircleAreas = [8; 4; 2; 1; 1];
    NI_SetCounts = zeros(2^NI_NumSets, 1);
    NI_SetCounts(2) = 8;
    NI_SetCounts(3) = 4;
    NI_SetCounts(5) = 2;
    NI_SetCounts([9 17]) = 1;

    % Circle areas and counts for circles with multiple disjoint
    % intersections
    MI_NumSets = 4;
    MI_CircleAreas = 8 * ones(MI_NumSets, 1);
    MI_SetCounts = ones(2^MI_NumSets, 1);
    MI_SetCounts([1 6 11])= 0;
    MI_SetCounts([2 3 5 9])= 2;

    % Circle areas and counts for circles with extremely different sizes
    CSS_NumSets = 4;
    CSS_CircleAreas = [3; 1; 1; 1000];
    CSS_SetCounts = zeros(2^CSS_NumSets, 1);
    CSS_SetCounts([11 13])= 1;
    CSS_SetCounts(10) = 3;
    CSS_SetCounts(9) = 995;

    testCircleAreasData = struct(...
        'EqualAreas', struct('CircleAreas', EA_CircleAreas, ...
            'NumSets', EA_NumSets, 'SetCounts', EA_SetCounts), ...
        'NoIntersections', struct('CircleAreas', NI_CircleAreas, ...
            'NumSets', NI_NumSets,'SetCounts', NI_SetCounts), ...
        'MultipleIntersections', struct('CircleAreas', MI_CircleAreas, ...
            'NumSets', MI_NumSets, 'SetCounts', MI_SetCounts), ...
        'ContrastingSetSizes', struct('CircleAreas', CSS_CircleAreas, ...
            'NumSets', CSS_NumSets, 'SetCounts', CSS_SetCounts));
end

function testInitialCentersData = createTestInitialCentersData
    % Create data for testing vennEulerDiagram.computeInitialCircleCenters

    % Equally sized sets which don't intersect
    NI_NumSets = 4;
    NI_SetCounts = zeros(12^NI_NumSets, 1);
    NI_SetCounts([2 3 5 9]) = 4;
    NI_SetCounts = NI_SetCounts / sum(NI_SetCounts);
    NI_CircleAreas = 4 * ones(NI_NumSets, 1);

    % Sets which have an equal number of elements in each disjoint
    % intersection for the same number of sets
    SymmetricNumSets = 4;
    SymmetricSetCounts = ones(2^SymmetricNumSets, 1);
    SymmetricSetCounts([1 6 11])= 0;
    SymmetricSetCounts([2 3 5 9])= 2;
    SymmetricSetCounts = SymmetricSetCounts / sum(SymmetricSetCounts);
    SymmetricCircleAreas = 8 * ones(SymmetricNumSets, 1);

    testInitialCentersData = struct(...
        'NoIntersections', struct('SetCounts', NI_SetCounts, ...
            'NumSets', NI_NumSets, 'CircleAreas', NI_CircleAreas), ...
        'SymmetricSets', struct('SetCounts', SymmetricSetCounts, ...
            'NumSets', SymmetricNumSets, 'CircleAreas', SymmetricCircleAreas));
end

function testFindDIAreasData = createTestFindDIAreasData
    % Create data for testing vennEulerDiagram.findDisjointIntersectionAreas

    numSides = 200;

    % Circles which are non-intersecting
    c1 = nsidedpoly(numSides, 'Center', [1 1], 'Radius', 1);
    c2 = nsidedpoly(numSides, 'Center', [1 -1], 'Radius', 1);
    c3 = nsidedpoly(numSides, 'Center', [-1 1], 'Radius', 1);
    c4 = nsidedpoly(numSides, 'Center', [-1 -1], 'Radius', 1);
    NI_CirclePolyshapes = [c1 c2 c3 c4];

    NI_NumSets = 4;
    NI_IntersectionAreas = zeros(2^NI_NumSets, 1);
    NI_IntersectionAreas([2 3 5 9]) = 25;
    NI_ShownIntersectionIndices = [2 3 5 9];

    % Circles with maximum number of disjoint intersections
    c1 = nsidedpoly(numSides, 'Center', [1 1], 'Radius', 1.5);
    c2 = nsidedpoly(numSides, 'Center', [1 -1], 'Radius', 1.5);
    c3 = nsidedpoly(numSides, 'Center', [-1 1], 'Radius', 1.5);
    c4 = nsidedpoly(numSides, 'Center', [-1 -1], 'Radius', 1.5);
    maxDI_CirclePolyshapes = [c1 c2 c3 c4];

    maxDI_NumSets = 4;
    maxDI_IntersectionAreas = zeros(2^maxDI_NumSets, 1);
    maxDI_IntersectionAreas([2 3 5 9]) = 18.2936;
    maxDI_IntersectionAreas([4 6 11 13]) = 6.4806;
    maxDI_IntersectionAreas([8 12 14 15]) = 0.1935;
    maxDI_IntersectionAreas(16) = 0.1291;
    maxDI_ShownIntersectionIndices = [2:6, 8:9, 11:16];

    % Test case where some circles are intersecting and other circles are
    % non intersecting (not all disjoint intersections are shown, only some)
    c1 = nsidedpoly(numSides, 'Center', [1 1], 'Radius', 1);
    c2 = nsidedpoly(numSides, 'Center', [1 -1], 'Radius', 1);
    c3 = nsidedpoly(numSides, 'Center', [-2 1], 'Radius', 1.5);
    c4 = nsidedpoly(numSides, 'Center', [-2 -1], 'Radius', 1.5);
    SI_CirclePolyshapes = [c1 c2 c3 c4];
    
    SI_NumSets = 4;
    SI_IntersectionAreas = zeros(2^SI_NumSets, 1);
    SI_IntersectionAreas([2 3]) = 16.6459;
    SI_IntersectionAreas([5 9]) = 29.2550;
    SI_IntersectionAreas(13) = 8.1982;
    SI_ShownIntersectionIndices = [2 3 5 9 13];

    % Test case where circles are nested
    c1 = nsidedpoly(numSides, 'Center', [0 0], 'Radius', 4);
    c2 = nsidedpoly(numSides, 'Center', [0 0], 'Radius', 3);
    c3 = nsidedpoly(numSides, 'Center', [0 0], 'Radius', 2);
    c4 = nsidedpoly(numSides, 'Center', [0 0], 'Radius', 1);
    NestedCirclePolyshapes = [c1 c2 c3 c4];

    NestedNumSets = 4;
    NestedIntersectionAreas = zeros(2^NestedNumSets, 1);
    NestedIntersectionAreas(2) = 43.75;
    NestedIntersectionAreas(4) = 31.25;
    NestedIntersectionAreas(8) = 18.75;
    NestedIntersectionAreas(16) = 6.25;
    NestedShownIntersectionIndices = [2 4 8 16];

    testFindDIAreasData = struct(...
        'NoIntersections', struct('NumSets', NI_NumSets, ...
            'CirclePolyshapes', NI_CirclePolyshapes, ...
            'IntersectionAreas', NI_IntersectionAreas, ...
            'ShownIntersectionIndices', NI_ShownIntersectionIndices), ...
        'MaxDisjointIntersections', struct('NumSets', maxDI_NumSets, ...
            'CirclePolyshapes', maxDI_CirclePolyshapes, ...
            'IntersectionAreas', maxDI_IntersectionAreas, ...
            'ShownIntersectionIndices', maxDI_ShownIntersectionIndices), ...
        'SomeIntersections', struct('NumSets', SI_NumSets, ...
            'CirclePolyshapes', SI_CirclePolyshapes, ...
            'IntersectionAreas', SI_IntersectionAreas, ...
            'ShownIntersectionIndices', SI_ShownIntersectionIndices), ...
        'NestedCircles', struct('NumSets', NestedNumSets, ...
            'CirclePolyshapes', NestedCirclePolyshapes, ...
            'IntersectionAreas', NestedIntersectionAreas, ...
            'ShownIntersectionIndices', NestedShownIntersectionIndices));
end

function testStressData = createTestStressData
    % Create data for testing vennEulerDiagram.computeStress

    % No stress when the intersection areas and set counts are
    % scalar multiples of one another
    NS_IntersectionAreas = [20 0 10 5 0 32 0 18 15]';
    NS_SetCounts = 3 * NS_IntersectionAreas;
    NS_StressLowerBound = 0;
    NS_StressUpperBound = 0;

    % Low stress when the intersection areas and set counts are close to
    % being scalar multiples of one another
    LS_IntersectionAreas = [20 0 10 5 0 32 0 18 15]';
    LS_SetCounts = [21 0 15 0 0 31 0 20 13]';
    LS_StressLowerBound = 0;
    LS_StressUpperBound = 0.5;

    % High stress when the intersection areas and set counts are far from
    % being scalar multiples of one another
    HS_IntersectionAreas = [20 0 10 5 0 32 0 18 15]';
    HS_SetCounts = flip([20 0 10 5 0 32 0 18 15]');
    HS_StressLowerBound = 0.5;
    HS_StressUpperBound = 1;

    % Max stress when the intersection areas and set counts have no indices
    % where they share positive values
    MS_IntersectionAreas = [20 0 10 5 0 32 0 18 15]';
    MS_SetCounts = [0 25 0 0 25 0 50 0 0]';
    MS_StressLowerBound = 1;
    MS_StressUpperBound = 1;

    testStressData = struct(...
        'NoStress', struct('IntersectionAreas', NS_IntersectionAreas, 'SetCounts', NS_SetCounts, ...
            'StressLowerBound', NS_StressLowerBound, 'StressUpperBound', NS_StressUpperBound), ...
        'LowStress', struct('IntersectionAreas', LS_IntersectionAreas, 'SetCounts', LS_SetCounts, ...
            'StressLowerBound', LS_StressLowerBound, 'StressUpperBound', LS_StressUpperBound), ...
        'HighStress', struct('IntersectionAreas', HS_IntersectionAreas, 'SetCounts', HS_SetCounts, ...
            'StressLowerBound', HS_StressLowerBound, 'StressUpperBound', HS_StressUpperBound), ...
        'MaxStress', struct('IntersectionAreas', MS_IntersectionAreas, 'SetCounts', MS_SetCounts, ...
            'StressLowerBound', MS_StressLowerBound, 'StressUpperBound', MS_StressUpperBound));
end

function testComputeGradientsData = createTestComputeGradientsData
    % Create data to test vennEulerDiagram.computeGradients

    numSides = 200;

    % When the sets have non-empty intersections and the circles are far apart, 
    % the gradients should be such that the corresponding circles move closer together
    c1 = nsidedpoly(numSides, 'Center', [1 1], 'Radius', 0.2);
    c2 = nsidedpoly(numSides, 'Center', [1 -1], 'Radius', 0.2);
    c3 = nsidedpoly(numSides, 'Center', [-1 -1], 'Radius', 0.2);
    c4 = nsidedpoly(numSides, 'Center', [-1 1], 'Radius', 0.2);
    CT_CirclePolyshapes = [c1 c2 c3 c4];

    CT_NumSets = 4;
    [CT_IntersectionAreas, ~, ~] = hTestVennEulerDiagram.callFindDisjointIntersectionAreas(CT_NumSets, CT_CirclePolyshapes);

    CT_SetCounts = ones(2^CT_NumSets, 1);
    CT_SetCounts([1 6 11])= 0;
    CT_SetCounts([2 3 5 9])= 2;
    CT_SetCounts = CT_SetCounts / sum(CT_SetCounts);
    CT_SetCounts = CT_SetCounts / sum(CT_SetCounts);
    
    CT_CurrCircleCenters = [1 1; 1 -1; -1 -1; -1 1];
    CT_GradientDirection = -1 * [1 1; 1 -1; -1 -1; -1 1];

    % When the sets have empty intersections and the circles are close together, 
    % the gradients should be such that the corresponding circles move farther apart
    c1 = nsidedpoly(numSides, 'Center', [0.1 0.1], 'Radius', 1);
    c2 = nsidedpoly(numSides, 'Center', [0.1 -0.1], 'Radius', 1);
    c3 = nsidedpoly(numSides, 'Center', [-0.1 -0.1], 'Radius', 1);
    c4 = nsidedpoly(numSides, 'Center', [-0.1 0.1], 'Radius', 1);
    FA_CirclePolyshapes = [c1 c2 c3 c4];

    FA_NumSets = 4;
    [FA_IntersectionAreas, ~, ~] = hTestVennEulerDiagram.callFindDisjointIntersectionAreas(FA_NumSets, FA_CirclePolyshapes);
    
    FA_SetCounts = zeros(2^FA_NumSets, 1);
    FA_SetCounts([2 3 5 9]) = 1;
    FA_SetCounts = FA_SetCounts / sum(FA_SetCounts);

    FA_CurrCircleCenters = 0.1 * [1 1; 1 -1; -1 -1; -1 1];
    FA_GradientDirection = [1 1; 1 -1; -1 -1; -1 1];

    testComputeGradientsData = struct(...
        'CloserTogether', struct('GradientDirection', CT_GradientDirection, ...
            'IntersectionAreas', CT_IntersectionAreas, ...
            'SetCounts', CT_SetCounts, ...
            'NumSets', CT_NumSets, ...
            'CurrCircleCenters', CT_CurrCircleCenters), ...
        'FartherApart', struct('GradientDirection', FA_GradientDirection, ...
            'IntersectionAreas', FA_IntersectionAreas, ...
            'SetCounts', FA_SetCounts, ...
            'NumSets', FA_NumSets, ...
            'CurrCircleCenters', FA_CurrCircleCenters));
end

function testSyntaxesData = createTestSyntaxesData
    % Create data to test valid syntaxes and property setting for vennEulerDiagram 
    setListData = {[1:3]; [3:5]; [3:6]}; %#ok<NBRAK> 
    setMembershipData = [1 0 0;
                          1 0 0;
                          1 1 1;
                          0 1 1;
                          0 1 1;
                          0 0 1];
    setLabels = ["A"; "B"; "C"];

    args1 = {setListData};
    args2 = {setListData setLabels};
    args3 = {setMembershipData};
    args4 = {setMembershipData setLabels};

    testSyntaxesData = struct( ...
        'SetListDataSyntax', struct('Args', {args1}, ...
            'SetListData', {setListData}, 'SetMembershipData', setMembershipData, 'SetLabels', strings(0, 1)), ...
        'SetListDataAndSetLabelsSyntax', struct('Args', {args2}, ...
            'SetListData', {setListData}, 'SetMembershipData', setMembershipData, 'SetLabels', setLabels), ...
        'SetMembershipDataSyntax', struct('Args', {args3}, ...
            'SetListData', {setListData}, 'SetMembershipData', setMembershipData, 'SetLabels', strings(0, 1)), ...
        'SetMembershipDataAndSetLabelsSyntax', struct('Args', {args4}, ...
            'SetListData', {setListData}, 'SetMembershipData', setMembershipData, 'SetLabels', setLabels), ...
        'EmptySyntax', struct('Args', {{{}}}, 'SetListData', {cell(0, 1)}, 'SetMembershipData', [], 'SetLabels', strings(0, 1)));
end
