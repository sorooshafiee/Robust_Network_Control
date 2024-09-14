classdef ScenarioNode < handle
    %----------------------------------------------------------------------
    % Simple structure that allows to store information for the constraints
    % of the problem
    %----------------------------------------------------------------------
    
    properties
        w
        u
        childVec
        parentNode
    end
    
     methods
        function obj = ScenarioNode(w,nu,pN)
            if nargin < 1
                obj.w = [];
                obj.u = [];
                obj.parentNode = [];
            else
                if size(w,1) == 1
                    % dimension check
                    w = w';
                end
                obj.w = w;
                obj.u = sdpvar(nu,1);
                obj.parentNode = pN;
            end
            
            obj.childVec = [];
        end
        
        function add_child(self,object)
            self.childVec = [self.childVec,object];
        end
     end
    
end