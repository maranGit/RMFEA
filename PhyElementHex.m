classdef PhyElementHex < PhyElement
    properties
        
    end
    methods
        function obj = PhyElementHex()
            
        end
        
        % detail of virtual function in PhyElement.m
        setGeometry(obj)
        
        setInternalMaterialProperties(obj)
        
        Calculate_ElementStiffness_Force(obj)
        
        SpecificOutput(obj)
        
    end
end