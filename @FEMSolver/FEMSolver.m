classdef FEMSolver < handle
    properties
        dim % int
        ndofpn % int
        nNodes % (int) number of nodes in the domain
        ne % number of elements ( int )
        nodes % vector<PhyNode>
        elements % vector<PhyElement*>
        matAll % collection of material models
        dofs % vector of all dofs (pass by reference)
        nf % number of free dofs
        np % number of prescribed dofs
        ndof % ndof = nf + np
        K % stiffness matrix
        F % Force vector (external load)
        Fp % vector of prescribed dof forces
        nmats % number of materials
        mults = 1;
        nstep = 1;
    end
    properties (Hidden)
        EBC % essential BC: [#dof, displacement]
        FBC % all natural BC: [#dof]
        NBC % natural BC with loading: {#dof, nodal loading}
        % edofs % position of prescriped dofs in dofs
        % ndofs % position of loaded dofs in dofs
    end
    methods
        function obj = FEMSolver(dim)		% works as function with one or zero argument
            if nargin == 0
                obj.dim = 2;
                obj.ndofpn = 2;
            else
                obj.dim = dim;
                obj.ndofpn = dim;
            end
        end
        
        % All steps together:
        function FEMSolve(obj, runName)
            %{
            % Linear elastic small strain solver
            obj.Input(runName);
            obj.Calculate_ElementStiffness_Force(); 
            obj.Assemble();
            d = obj.Solve_Dofs();
            obj.Assign_dof(d);
            obj.UpdateFpNodalPrescribedForces;
            %}
            %=========================
            % small strain nonlinear material solver
            obj.Input(runName);
            fid = fopen('Output','w');
            for iStep = 1:obj.nstep
                % update EBC and NBC of certain dof for current step
                % since they are stored in element and node as reference
                % update here is enough
                EBC2 = obj.EBC;
                incEBC = obj.mults(iStep) * EBC2(:,2);
                NBC2 = obj.NBC;
                incNBC = obj.mults(iStep) * NBC2(:,2);
                for temp = 1:size(EBC2, 1)
                    currdof = obj.dofs(EBC2(temp,1));
                    currdof.v = currdof.v+incEBC(temp); % increase EBC
                end
                for temp = 1:size(NBC2, 1)
                    currdof = obj.dofs(NBC2(temp,1));
                    currdof.f = currdof.f+incNBC(temp); % increase NBC
                end
                %
                % update residual and stiffness of each element
                for temp = 1:obj.ne
                    obj.elements(temp).Calculate_ElementStiffness_Force(obj.matAll);
                end
                % obj.elements.Calculate_ElementStiffness_Force(obj.matAll); % should include one material model as input @@@@@
                %
                % assemble Fint and Fext, then compute residual
                Fext = [obj.dofs(obj.FBC).f]';
                Fint = zeros(obj.nf, 1);
                for iElem = 1:obj.ne
                    currFint = obj.elements(iElem).Fint;
                    currMap = obj.elements(iElem).dofMap;
                    freeFint = currFint(currMap>0);
                    freeMap = currMap(currMap>0);
                    Fint(freeMap) = Fint(freeMap) + freeFint;
                end
                Res = Fext - Fint;
                Res_0 = norm(Res);
                iter = 0;
                % Newton-Raphson loop for displacement solution
                while (norm(Res) > (1e-12)*Res_0) % hard-coded tolerance for now
                    % Assemble stiffness matrix
                    obj.F = Res;
                    obj.Assemble();
                    % compute disp
                    deltaD = obj.K \ obj.F;
                    % update disp (consider replacing for loop in future version!)
                    for iFDof = 1:obj.nf % update displacement of ith free dof
                        upDOF = obj.FBC(iFDof);
                        obj.dofs(upDOF).vUpdate(deltaD(iFDof));
                    end
                    % compute new residual and stiffness
                    for temp = 1:obj.ne
                        obj.elements(temp).Calculate_ElementStiffness_Force(obj.matAll);
                    end
                    %obj.elements.Calculate_ElementStiffness_Force(obj.matAll); % should include one material model as input @@@@@
                    % assemble Fint and Fext, then compute residual
                    Fint = zeros(obj.nf, 1);
                    for iElem = 1:obj.ne
                        currFint = obj.elements(iElem).Fint;
                        currMap = obj.elements(iElem).dofMap;
                        freeFint = currFint(currMap>0);
                        freeMap = currMap(currMap>0);
                        Fint(freeMap) = Fint(freeMap) + freeFint;
                    end
                    Res = Fext - Fint;
                    iter = iter + 1;
                    norm(Res);
                end
                %
                % update hardening variables
                obj.elements.Update;
                for temp = 1:obj.ndof
                    obj.dofs(temp).Update;
                end
                %
                % output state variables of each element
                % format: [ele1_gp1; ele1_gp2; ... ; ele(i)_gp(j)]
                stress = [obj.elements.stress];
                strain = [obj.elements.strain];
                disp=[obj.dofs.v];
                if obj.dim == 2
                    fprintf(fid, 'step %d\r\n', iStep);
                    fprintf(fid, '  displacement\r\n');
                    format = strcat(repmat('%5.2e, ', 1, 7), '%5.2e\r\n'); % work now only for one element model with 4 nodes
                    fprintf(fid, format, disp);
                    fprintf(fid, '  stress\r\n');
                    fprintf(fid, '  %5.6e, %5.6e, %5.6e\r\n', stress);
                    fprintf(fid, '  strain\r\n');
                    fprintf(fid, '  %5.6e, %5.6e, %5.6e\r\n', strain);
                elseif obj.dim == 3
                    fprintf(fid, '%5.2e, %5.2e, %5.2e, %5.2e, %5.2e, %5.2e\n', stress);
                    fprintf(fid, '%5.2e, %5.2e, %5.2e, %5.2e, %5.2e, %5.2e\n', strain);
                end
            end
            fclose(fid);
        end
        
        plotMesh2(obj, MNE_ID)
    end
    
    methods (Access = private)
        Input(obj, runName)
        
        % Step 3: compute nf from ndof and np, initialize stiffness matrix and force vector
        setSizes(obj)
        
        % Step 4: set prescribed dofs: already done when reading the input file
        % Step 5: Set global free nodal dof: already done when reading the input file
        % Step 6 and Step 7: dof positions; Step 7: Set F
        setPositions_F(obj)
        
        % Step 8: Element dof maps Me
        % Step 9: Set element dofs ae
        setElementDofMap_ae(obj)
        
        % Step 10: Compute element stiffness
        Calculate_ElementStiffness_Force(obj)
        
        % Step 11: Assembly from local to global system
        Assemble(obj)
        
        % Step 12: Solve global (free) dof a from Ka = F
        % successful solution returns true
        d = Solve_Dofs(obj)
        
        % Step 13: Assign a to nodes and elements
        Assign_dof(obj, d)
        
        % Step 14: Compute prescribed dof forces
        UpdateFpNodalPrescribedForces(obj)
        
    end
end