classdef vehicle < dynamicprops
    %VEHICLE Summary of this class goes here
    %   Detailed explanation goes here
    properties (SetAccess = private)
        length
        boatTailLength
        diameter
        baseDiameter
        structMass
        inertia_CM
        
        rCMpad_B
        rCMair_B
        rCP_B
        
        numFins
        finLE
        finThickness
        finLength
        finRootChord
        finTipChord
        finSweep
        finIncidence

        motor
        
        initPosVecGnd
        initVelVecBdy
        initEulAng
        initAngVelVec
    end
    
    properties (Dependent)
        padMass
        airMass
        
        finPlanformArea
        M6x6_B
    end
    
    methods
        %% Constructor
        function obj = vehicle
            % mass, volume and inertia
            obj.length              = SIM.parameter('Unit','m','Description','Vehicle length');
            obj.boatTailLength      = SIM.parameter('Unit','m','Description','Boattail length');
            obj.diameter            = SIM.parameter('Unit','m','Description','Vehicle diameter');
            obj.baseDiameter        = SIM.parameter('Unit','m','Description','Diameter of the base of the rocket');
            obj.structMass          = SIM.parameter('Unit','kg','Description','Vehicle mass without motor');
            obj.inertia_CM          = SIM.parameter('Unit','kg*m^2','Description','Inertia Matrix');
                        
            %Important Point Locations
            obj.rCMpad_B            = SIM.parameter('Value',[0;0;0],'Unit','m','Description','Vector going from the nose to the CM before motor burn');
            obj.rCMair_B            = SIM.parameter('Value',[0;0;0],'Unit','m','Description','Vector going from the nose to the CM before motor burn');
            obj.rCP_B               = SIM.parameter('Value',[0;0;0],'Unit','m','Description','Vector going from the nose to the center of pressure');
            
            % Overall Wing Properties (Used to create portWing and stbdWing
            obj.numFins       = SIM.parameter('Description','Number of fins','NoScale',true);
            obj.finLE         = SIM.parameter('Unit','m','Description','Distance from nose to fin leading edge');
            obj.finThickness  = SIM.parameter('Unit','m','Description','Thickness of the fin');
            obj.finLength     = SIM.parameter('Unit','m','Description','Length of the fin');
            obj.finRootChord  = SIM.parameter('Unit','m','Description','Fin root chord');
            obj.finTipChord   = SIM.parameter('Unit','m','Description','Fin tip chord');
            obj.finSweep      = SIM.parameter('Unit','deg','Description','Fin sweep angle');
            obj.finIncidence  = SIM.parameter('Unit','deg','Description','Fin flow incidence angle');
                        
            obj.motor = VEH.motor;
                        
            % initial conditions 
            obj.initPosVecGnd           = SIM.parameter('Unit','m','Description','Initial CM position represented in the inertial frame');
            obj.initVelVecBdy           = SIM.parameter('Unit','m/s','Description','Initial CM velocity represented in the body frame ');
            obj.initEulAng              = SIM.parameter('Unit','rad','Description','Initial Euler angles');
            obj.initAngVelVec           = SIM.parameter('Unit','rad/s','Description','Initial angular velocity vector');
        end
        function obj = buildMotor(obj)
            str = input('What is the name of the motor?\n','s');
            obj.motor.setMotorName(str,'');
            fileLoc = which(obj.motor.motorName.Value);
            if ~isfile(fileLoc)
                fprintf('The file containing the motor data file does not exist.\n');
                mass = input('Total motor mass in kg?\n');
                fuel = input('Fuel mass in kg?\n');
                T = input('Average thrust in N?\n');
                t = input('Burn time?\n');
                obj.motor.motorMass = SIM.parameter('Value',mass,'Unit','kg','Description','Total motor mass');
                obj.motor.fuelMass = SIM.parameter('Value',fuel,'Unit','kg','Description','Fuel mass');
                obj.motor.thrustVec = SIM.parameter('Value',[0 linspace(T,T,20) 0 0]','Unit','N','Description','Motor thrust curve');
                obj.motor.burnTimeVec = SIM.parameter('Value',[linspace(0,t,21) t+0.1 t+0.2]','Unit','s','Description','Motor thrust curve time reference');
            else
                temp = load(fileLoc);
                obj.motor.motorMass = SIM.parameter('Value',temp.m_m,'Unit','kg','Description','Total motor mass');
                obj.motor.fuelMass = SIM.parameter('Value',temp.m_f,'Unit','kg','Description','Fuel mass');
                obj.motor.thrustVec = SIM.parameter('Value',temp.F_t,'Unit','N','Description','Motor thrust curve');
                obj.motor.burnTimeVec = SIM.parameter('Value',temp.T_b,'Unit','s','Description','Motor thrust curve time reference');
            end
        end
        %% setters
        function setFluidCoeffsFileName(obj,val,units)
            if ~endsWith(val,'.mat')
                val = [val '.mat'] ;
            end
            obj.fluidCoeffsFileName.setValue(val,units);
        end

        function setD6x6_LE(obj,val,units)
            obj.D6x6_LE.setValue(val,units);
        end

        function setRB_LE(obj,val,units)
            obj.rB_LE.setValue(val(:),units);
        end

        function setRCM_LE(obj,val,units)
            obj.rCM_LE.setValue(val(:),units);
        end

        function setInitPosVecGnd(obj,val,units)
            obj.initPosVecGnd.setValue(val(:),units);
        end

        function setInitVelVecBdy(obj,val,units)
            obj.initVelVecBdy.setValue(val(:),units);
        end

        function setInitEulAng(obj,val,units)
            obj.initEulAng.setValue(val(:),units);
        end

        function setInitAngVelVec(obj,val,units)
            obj.initAngVelVec.setValue(val(:),units);
        end
        
        %% getters
        % mass
        function val = get.padMass(obj)
            val = SIM.parameter('Value',obj.structMass.Value+obj.motor.motorMass.Value,...
                'Unit','kg','Description','Vehicle mass on the launch pad');
        end
        
        function val = get.airMass(obj)
            val = SIM.parameter('Value',obj.structMass.Value+obj.motor.motorMass.Value-obj.motor.fuelMass.Value,...
                'Unit','kg','Description','Vehicle mass after motor burn');
        end
                            
        % aerodynamic reference area
        function val = get.finPlanformArea(obj)
            Sref = 1/2*obj.finLength.Value*(obj.finRootChord.Value+obj.finTipChord.Value);
            val = SIM.parameter('Value',Sref,'Unit','m^2',...
                'Description','Fin planform area');
        end
        
        function val = get.M6x6_B(obj)
            S = @(v) [0 -v(3) v(2);v(3) 0 -v(1);-v(2) v(1) 0];
            M = zeros(6,6);
            M(1,1) = obj.airMass.Value;
            M(2,2) = obj.airMass.Value;
            M(3,3) = obj.airMass.Value;
            M(1:3,4:6) = -obj.airMass.Value*S(obj.rCMair_B.Value);
            M(4:6,1:3) = obj.airMass.Value*S(obj.rCMair_B.Value);
            x = obj.rCMair_B.Value(1);
            y = obj.rCMair_B.Value(2);
            z = obj.rCMair_B.Value(3);
            M(4:6,4:6) = obj.inertia_CM.Value + (obj.airMass.Value * ...
                        [y^2 + z^2, -x*y     , -x*z;...
                         -x*y     , x^2 + z^2, -y*z;...
                         -x*z     , -y*z     , x^2 + y^2]);
            val = SIM.parameter('Value',M,'Unit','','Description',...
                '6x6 Mass-Inertia Matrix with origin at Wing LE Mid-Span');
        end
        
        %% other methods
        % Function to scale the object
        function obj = scale(obj,lengthScaleFactor,densityScaleFactor)
            props = findAttrValue(obj,'SetAccess','private');
            for ii = 1:numel(props)
                obj.(props{ii}).scale(lengthScaleFactor,densityScaleFactor);
            end
        end
        
        %Sets initial conditions on the path at the specified pathVariable
        function setICsOnPath(obj,initPathVar,pathFunc,geomParams,pathCntrPt,speed) %#ok<INUSL>
            % Sets initial conditions of the vehicle to be on the path
            [initPos,initVel] = eval(sprintf('%s(initPathVar,geomParams,pathCntrPt)',pathFunc));
            obj.setInitPosVecGnd(initPos,'m');
            obj.setInitVelVecBdy([-speed 0 0],'m/s');
            % Initial body z points radially out
            bdyZ = (initPos(:)-pathCntrPt(:))./sqrt(sum((initPos(:)-pathCntrPt(:)).^2));
            % Initial body x points backwards (opposite velocity(
            bdyX = -initVel;
            % Initial body y is calculated from the cross product of z & x
            bdyY = cross(bdyZ,bdyX);
            % Calculate euler angles from the rotation matrix
            obj.setInitEulAng(flip(rotm2eul([bdyX(:)'; bdyY(:)'; bdyZ(:)']')),'rad')
            % Initial angular velocity is zero
            obj.setInitAngVelVec([0 0 0],'rad/s');
        end
                
        function h = plot(obj,varargin)
            p = inputParser;
            addParameter(p,'FigHandle',[],@(x) isa(x,'matlab.ui.Figure'));
            addParameter(p,'AxHandle',[],@(x) isa(x,'matlab.graphics.axis.Axes'));
            addParameter(p,'EulerAngles',[0 0 0],@isnumeric)
            addParameter(p,'Position',[0 0 0]',@isnumeric)
            addParameter(p,'fuseRings',8,@isnumeric);
            addParameter(p,'Basic',false,@islogical) % Only plots aero surfaces if true
            parse(p,varargin{:})
            
            R = rotation_sequence(p.Results.EulerAngles);
            
            if isempty(p.Results.FigHandle) && isempty(p.Results.AxHandle)
                h.fig = figure;
                h.fig.Name ='Design';
            else
                h.fig = p.Results.FigHandle;
            end
            
            if isempty(p.Results.AxHandle)
                h.ax = gca;
            else
                h.ax = p.Results.AxHandle;
            end
            
            fs = obj.getPropsByClass("VEH.aeroSurf");
            % Aero surfaces (and fuselage)
            for ii = 1:4
                pts = R*obj.(fs{ii}).outlinePtsBdy.Value;
                h.surf{ii} = plot3(h.ax,...
                    pts(1,:)+p.Results.Position(1),...
                    pts(2,:)+p.Results.Position(2),...
                    pts(3,:)+p.Results.Position(3),...
                    'LineWidth',1.2,'Color','k','LineStyle','-',...
                    'DisplayName','Fluid Dynamic Surfaces');
                hold on
            end
            if p.Results.fuseRings == 0
                fusepts = [obj.fuse.rNose_LE.Value obj.fuse.rEnd_LE.Value];
                h.surf{5} = plot3(h.ax,fusepts(1,:),fusepts(2,:),fusepts(3,:),...
                                  'LineWidth',1.2,'Color','k','LineStyle','-',...
                                  'DisplayName','Fluid Dynamic Surfaces');
            else
                x=linspace(obj.fuse.rNose_LE.Value(1)+obj.fuse.diameter.Value/2,obj.fuse.rEnd_LE.Value(1)-obj.fuse.diameter.Value/2,p.Results.fuseRings);
                perSlice = 10;
                x = reshape(repmat(x,perSlice,1),[1 numel(x)*perSlice]);
                th=linspace(0,2*pi,perSlice);
                d=obj.fuse.diameter.Value;
                y=repmat(d/2*cos(th),1,p.Results.fuseRings);
                z=repmat(d/2*sin(th),1,p.Results.fuseRings);
                numextra=(perSlice-1)*2;
                xend=x(end);
                for i = 0:numextra-1
                    if mod(i,4)==0 || mod(i,4)==3
                        x(end+1)=x(1);
                    else
                        x(end+1)=xend;
                    end
                end
                y(end+1:end+numextra) = reshape(repmat(y(2:perSlice),2,1),[1 numextra]);
                z(end+1:end+numextra) = reshape(repmat(z(2:perSlice),2,1),[1 numextra]);
                
                [sx, sy, sz]=sphere;
                sx = reshape(obj.fuse.diameter.Value/2*sx,[1 numel(sx)]);
                sy = reshape(obj.fuse.diameter.Value/2*sy,[1 numel(sy)]);
                sz = reshape(obj.fuse.diameter.Value/2*sz,[1 numel(sz)]);
                nosex = sy(1:ceil(numel(sx)/2))+obj.fuse.rNose_LE.Value(1)+obj.fuse.diameter.Value/2;
                nosey = sx(1:ceil(numel(sx)/2));
                nosez = sz(1:ceil(numel(sx)/2));
%                 x(end+1:end+numel(nosex))=nosex;
%                 y(end+1:end+numel(nosey))=nosey;
%                 z(end+1:end+numel(nosez))=nosez;
                
                endx = sy(ceil(numel(sx)/2):end)+obj.fuse.rEnd_LE.Value(1)-obj.fuse.diameter.Value/2;
                endy = sx(ceil(numel(sx)/2):end);
                endz = sz(ceil(numel(sx)/2):end);
%                 x(end+1:end+numel(endx))=endx;
%                 y(end+1:end+numel(endy))=endy;
%                 z(end+1:end+numel(endz))=endz;
                
                
                h.surf{5}=plot3(h.ax,x,y,z,'LineWidth',.2,'Color','k','LineStyle','-',...
                      'DisplayName','Fluid Dynamic Surfaces');
                h.surf{6}=plot3(h.ax,nosex,nosey,nosez,'LineWidth',.2,'Color','k','LineStyle','-',...
                      'DisplayName','Fluid Dynamic Surfaces');
                h.surf{7}=plot3(h.ax,endx,endy,endz,'LineWidth',.2,'Color','k','LineStyle','-',...
                      'DisplayName','Fluid Dynamic Surfaces');
            end
                         
            if ~p.Results.Basic
                % Tether attachment points
                for ii = 1:obj.numTethers.Value
                    pts = R*obj.thrAttchPts_B(ii).posVec.Value;
                    h.thrAttchPts{ii} = plot3(h.ax,...
                        pts(1)+p.Results.Position(1),...
                        pts(2)+p.Results.Position(2),...
                        pts(3)+p.Results.Position(3),...
                        'r+','DisplayName','Tether Attachment Point');
                end
                % Turbines
                for ii = 1:obj.numTurbines.Value
                    pts = eval(sprintf('R*obj.turb%d.attachPtVec.Value',ii));
                    h.turb{ii} = plot3(h.ax,...
                        pts(1)+p.Results.Position(1),...
                        pts(2)+p.Results.Position(2),...
                        pts(3)+p.Results.Position(3),...
                        'm+','DisplayName','Turbine Attachment Point');
                end
                
                for ii = 1:4
                    pts = R*obj.fluidMomentArms.Value(:,ii);
                    h.momArms{ii} = plot3(h.ax,...
                        pts(1)+p.Results.Position(1),...
                        pts(2)+p.Results.Position(2),...
                        pts(3)+p.Results.Position(3),...
                        'b+','DisplayName','Fluid Dynamic Center');
                    
                end
                % Center of mass
                h.centOfMass = plot3(h.ax,...
                                    obj.rCM_LE.Value(1)+p.Results.Position(1),...
                                    obj.rCM_LE.Value(2)+p.Results.Position(2),...
                                    obj.rCM_LE.Value(3)+p.Results.Position(3),...
                                    'r*','DisplayName','Center of Mass');
                % Coordinate origin
                h.origin = plot3(h.ax,p.Results.Position(1),p.Results.Position(2),p.Results.Position(3),'kx','DisplayName','Body Frame Origin/Leading Edge');
                legend(h.ax,[h.surf{1} h.thrAttchPts{1} h.turb{1} h.momArms{2} h.centOfMass h.origin],'Location','northeast')
            end
            grid on
            axis equal
            xlabel('X (m)')
            ylabel('Y (m)')
            zlabel('Z (m)')
            view(-45,30)
            
            set(gca,'DataAspectRatio',[1 1 1])
        end   
        
        function plotCoeffPolars(obj)
            fh = findobj( 'Type', 'Figure', 'Name', 'Partitioned Aero Coeffs');
            
            if isempty(fh)
                fh = figure;
                fh.Position =[102 92 3*560 2*420];
                fh.Name ='Partitioned Aero Coeffs';
            else
                figure(fh);
            end
            
            % left wing
            ax1 = subplot(4,4,1);
            plot(obj.portWing.alpha.Value,obj.portWing.CL.Value);
            hWingCL_ax = gca;
            
            xlabel('$\alpha$ [deg]')
            ylabel('$C_{L}$')
            title('Port Wing')
            grid on
            hold on
            
            ax5 = subplot(4,4,5);
            plot(obj.portWing.alpha.Value,obj.portWing.CD.Value);
            xlabel('$\alpha$ [deg]')
            ylabel('$C_{D}$')
            grid on
            hold on
            hWingCD_ax = gca;
            
            ax9 = subplot(4,4,9);
            plot(obj.portWing.alpha.Value,obj.portWing.CL.Value(:)./obj.portWing.CD.Value(:))
            xlabel('$\alpha$ [deg]')
            ylabel('$\frac{C_{L}}{C_D}$')
            grid on
            hold on
            
            ax13 = subplot(4,4,13);
            plot(obj.portWing.alpha.Value,(obj.portWing.CL.Value(:).^3)./(obj.portWing.CD.Value(:).^2))
            xlabel('$\alpha$ [deg]')
            ylabel('$\frac{C_{L}^3}{C_D^2}$')
            grid on
            hold on
                       
            linkaxes([ax1,ax5,ax9,ax13],'x');
            
            % right wing
            ax2 = subplot(4,4,2);
            plot(obj.stbdWing.alpha.Value,obj.stbdWing.CL.Value);
            
            xlabel('$\alpha$ [deg]')
            ylabel('$C_{L}$')
            title('Starboard Wing')
            grid on
            hold on
            
            ax6 = subplot(4,4,6);
            plot(obj.stbdWing.alpha.Value,obj.stbdWing.CD.Value);
            xlabel('$\alpha$ [deg]')
            ylabel('$C_{D}$')
            grid on
            hold on
            
            ax10 = subplot(4,4,10);
            plot(obj.stbdWing.alpha.Value,obj.stbdWing.CL.Value(:)./obj.stbdWing.CD.Value(:))
            xlabel('$\alpha$ [deg]')
            ylabel('$\frac{C_{L}}{C_D}$')
            grid on
            hold on
            
            ax14 = subplot(4,4,14);
            plot(obj.stbdWing.alpha.Value,(obj.stbdWing.CL.Value(:).^3)./(obj.stbdWing.CD.Value(:).^2))
            xlabel('$\alpha$ [deg]')
            ylabel('$\frac{C_{L}^3}{C_D^2}$')
            grid on
            hold on
            
            linkaxes([ax2,ax6,ax10,ax14],'x');
            
            % HS
            ax3 = subplot(4,4,3);
            plot(obj.hStab.alpha.Value,obj.hStab.CL.Value);
            hhStabCL_ax = gca;
            
            xlabel('$\alpha$ [deg]')
            ylabel('$C_{L}$')
            title('Horizontal stabilizer')
            grid on
            hold on
            
            ax7 = subplot(4,4,7);
            plot(obj.hStab.alpha.Value,obj.hStab.CD.Value);
            hhStabCD_ax = gca;
            xlabel('$\alpha$ [deg]')
            ylabel('$C_{D}$')
            grid on
            hold on
            
            ax11 = subplot(4,4,11);
            plot(obj.hStab.alpha.Value,obj.hStab.CL.Value(:)./obj.hStab.CD.Value(:))
            xlabel('$\alpha$ [deg]')
            ylabel('$\frac{C_{L}}{C_D}$')
            grid on
            hold on
            
            ax15 = subplot(4,4,15);
            plot(obj.hStab.alpha.Value,(obj.hStab.CL.Value(:).^3)./(obj.hStab.CD.Value(:).^3))
            xlabel('$\alpha$ [deg]')
            ylabel('$\frac{C_{L}^3}{C_D^2}$')
            grid on
            hold on
            
            linkaxes([ax3,ax7,ax11,ax15],'x');
            
            % VS
            ax4 = subplot(4,4,4);
            plot(obj.vStab.alpha.Value,obj.vStab.CL.Value);
            hvStabCL_ax = gca;
            xlabel('$\alpha$ [deg]')
            ylabel('$C_{L}$')
            title('Vertical stabilizer')
            grid on
            hold on
            
            ax8 = subplot(4,4,8);
            plot(obj.vStab.alpha.Value,obj.vStab.CD.Value);
            hvStabCD_ax = gca;
            xlabel('$\alpha$ [deg]')
            ylabel('$C_{D}$')
            grid on
            hold on
            
            ax12 = subplot(4,4,12);
            plot(obj.vStab.alpha.Value,obj.vStab.CL.Value(:)./obj.vStab.CD.Value(:))
            xlabel('$\alpha$ [deg]')
            ylabel('$\frac{C_{L}}{C_D}$')
            grid on
            hold on
            
            ax16 = subplot(4,4,16);
            plot(obj.vStab.alpha.Value,(obj.vStab.CL.Value(:).^3)./(obj.vStab.CD.Value(:).^3))
            xlabel('$\alpha$ [deg]')
            ylabel('$\frac{C_{L}^3}{C_D^2}$')
            grid on
            hold on
            
            linkaxes([ax4,ax8,ax12,ax16],'x');
            
            %             axis([ax1 ax2 ax3 ax4],[-20 20 ...
%                 min([hWingCL_ax.YLim(1),hhStabCL_ax.YLim(1),hvStabCL_ax.YLim(1)])...
%                 max([hWingCL_ax.YLim(2),hhStabCL_ax.YLim(2),hvStabCL_ax.YLim(2)])]);
%             axis([ax5 ax6 ax7 ax8],[-20 20 ...
%                 min([hWingCD_ax.YLim(1),hhStabCD_ax.YLim(1),hvStabCD_ax.YLim(1)])...
%                 max([hWingCD_ax.YLim(2),hhStabCD_ax.YLim(2),hvStabCD_ax.YLim(2)])]);
            
        end
        
        %Get a struct of parameters of the desired class
        [output,varargout] = struct(obj,className);
        
        %returns a cell array of properties of the desired class
        output = getPropsByClass(obj,className);
               
        % calculate max equilibrium speed at a given azimuth, elevation, and
        % flow velocity vector
        [sp,tanRoll,velAng] = eqSpeed(obj,vf,az,el)
        
        % Functions to animate the vehicle
        val = animateSim(obj,tsc,timeStep,varargin)
        val = animateBody(obj,tsc,timeStep,varargin)
    end % methods
end