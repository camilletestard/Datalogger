classdef hTestVennEulerDiagram < vennEulerDiagram
    % hTestVennEulerDiagram Class for testing vennEulerDiagram properties
    % and methods. Used in tVennEulerDiagram.
    % 
    %   Copyright 2021 The MathWorks, Inc.

    methods
        function callSetup(obj)
            setup(obj)
        end
        
        function callUpdate(obj)
            update(obj)
        end
        
        function visible = callVerifyDataProperties(obj)
            visible = verifyDataProperties(obj);
        end
    end
    
    methods (Static)
        function setMembershipData = createSetMembershipData(setListData)
            h = hTestVennEulerDiagram(setListData);
            h.callUpdate();
            setMembershipData = h.SetMembershipData;
        end

        function setCounts = createSetCounts(setListData)
            h = hTestVennEulerDiagram(setListData);
            h.callUpdate();
            setCounts = h.SetCounts;
        end

        function scaledDiameters = callComputeCircleDiameters(circleAreas)
            scaledDiameters = hTestVennEulerDiagram.computeCircleDiameters(numel(circleAreas), circleAreas);
        end

        function circleAreas = callComputeCircleAreas(numSets, setCounts)
            circleAreas = hTestVennEulerDiagram.computeCircleAreas(numSets, setCounts);
        end

        function initialCenters = callComputeInitialCircleCenters(setCounts, numSets, circleAreas)
            initialCenters = hTestVennEulerDiagram.computeInitialCircleCenters(setCounts, numSets, circleAreas);
        end

        function [intersectionAreas, intersectionPolyshapes, shownIntersectionIndices] = ...
                callFindDisjointIntersectionAreas(numSets, circlePolyshapes) 
            [intersectionAreas, intersectionPolyshapes, shownIntersectionIndices] = ...
                hTestVennEulerDiagram.findDisjointIntersectionAreas(numSets, circlePolyshapes); 
        end

        function stress = callComputeStress(intersectionAreas, setCounts)
            [stress, ~] = hTestVennEulerDiagram.computeStress(intersectionAreas, setCounts);
        end

        function areas_hat = computeAreasHat(intersectionAreas, setCounts)
            [~, areas_hat] = hTestVennEulerDiagram.computeStress(intersectionAreas, setCounts);
        end

        function gradients = callComputeGradients(numSets, currCircleCenters, areas, areas_hat)
            gradients = hTestVennEulerDiagram.computeGradients(numSets, currCircleCenters, areas, areas_hat);
        end
    end
end
