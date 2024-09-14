classdef ConstraintInfo < handle
    %----------------------------------------------------------------------
    % Simple structure that allows to store information for the constraints
    % of the problem
    %----------------------------------------------------------------------
    
    properties
        tag     % tag of the constraint 
        dualVar     % dual variable index
        dualVarTag  % dual variable tag (eq or in)
        BindScen_dualVars   % indices of binding scenario dual variables
        BindScen_dualVarsTag    % tags (eq or in)
        BindScen_dualIdxs   % indices for position in the w vector
        dualVar_Val     % value of the dualVar
        BindScen_dualVars_Val   % values for the binding scenario dual variables
    end
    
     methods
        function obj = ConstraintInfo()
            obj.tag = 'None';
            obj.dualVar = [];
            obj.dualVarTag = [];
            obj.BindScen_dualVars = [];
            obj.BindScen_dualVarsTag = [];
            obj.BindScen_dualIdxs = [];
            
            obj.dualVar_Val = nan;
            obj.BindScen_dualVars_Val = nan;
        end
     end
    
end