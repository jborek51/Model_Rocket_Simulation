classdef vehicle < dynamicprops
    %VEHICLE Summary of this class goes here
    %   Detailed explanation goes here
    properties (SetAccess = private)
        length
        diameter
        structMass
        inertia_CM
        
        rCMpad_B
        rCMair_B
        rCP_B
        selectCP
        
        numFins
        finLE
        finThickness
        finLength
        finRootChord
        finTipChord
        finSweep
        finIncidence
        
        K
        noseShape
        noseLength
        CnNose
        bodyLength
        boatLE
        boatTailLength
        baseDiameter

        motor
        
        initPosVecGnd
        initVelVecBdy
        initEulAng
        initAngVelVec
    end
    
    properties (Dependent)
        padMass
        airMass
        
        CnFin
        CnBoat
        CnRocket
        XcpNose
        XcpFin
        XcpBoat
        XcpBody
        rCP_Bx
        
        Ar
        Ap
        
        finOutline
        finPlanformArea
        M6x6_B
        initPosVecNose
    end
    
    methods
        %% Constructor
        function obj = vehicle
            % size mass, volume and inertia
            obj.length              = SIM.parameter('Unit','m','Description','Vehicle length');
            obj.bodyLength          = SIM.parameter('Unit','m','Description','Body-tube length');
            obj.boatTailLength      = SIM.parameter('Unit','m','Description','Boattail length');
            obj.noseLength          = SIM.parameter('Unit','m','Description','Nosecone length');
            obj.diameter            = SIM.parameter('Unit','m','Description','Vehicle diameter');
            obj.baseDiameter        = SIM.parameter('Unit','m','Description','Diameter of the base of the rocket');
            obj.structMass          = SIM.parameter('Unit','kg','Description','Vehicle mass without motor');
            obj.inertia_CM          = SIM.parameter('Unit','kg*m^2','Description','Inertia Matrix');
                        
            %Important Point Locations
            obj.rCMpad_B            = SIM.parameter('Value',[0;0;0],'Unit','m','Description','Vector going from the nose to the CM before motor burn');
            obj.rCMair_B            = SIM.parameter('Value',[0;0;0],'Unit','m','Description','Vector going from the nose to the CM before motor burn');
            obj.rCP_B               = SIM.parameter('Value',[0;0;0],'Unit','m','Description','Vector going from the nose to the center of pressure');
            obj.selectCP            = SIM.parameter('Value',0,'Unit','','Description','CP method selection: 0 = set;  1 = dynamic');
            
            % Overall fin properties 
            obj.numFins         = SIM.parameter('Unit','','Description','Number of fins','NoScale',true);
            obj.finLE           = SIM.parameter('Unit','m','Description','Distance from nose to fin leading edge');
            obj.finThickness    = SIM.parameter('Unit','m','Description','Thickness of the fin');
            obj.finLength       = SIM.parameter('Unit','m','Description','Length of the fin');
            obj.finRootChord    = SIM.parameter('Unit','m','Description','Fin root chord');
            obj.finTipChord     = SIM.parameter('Unit','m','Description','Fin tip chord');
            obj.finSweep        = SIM.parameter('Unit','deg','Description','Fin sweep angle');
            obj.finIncidence    = SIM.parameter('Unit','deg','Description','Fin flow incidence angle');
            
            obj.K               = SIM.parameter('Value',1,'Unit','','Description','Rocket body lift correction constant');
            obj.noseShape       = SIM.parameter('Description','Nosecone shape (Conical, Ogive, Parabolic, Haack)','NoScale',true);
            obj.CnNose          = SIM.parameter('Value',2,'Unit','','Description','Nosecone stability derivative');
            obj.boatLE          = SIM.parameter('Unit','m','Description','Distance from nose to boattail leading edge');
            
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
                          
        function val = get.CnFin(obj)
            N = obj.numFins.Value;
            df = obj.diameter.Value;
            Ls = obj.finLength.Value;
            Lr = obj.finRootChord.Value;
            Lt = obj.finTipChord.Value;
            sw = obj.finSweep.Value;
            V1 = [Ls,Ls/cosd(sw)*sind(sw)+Lt/2];
            V2 = [0,Lr/2];
            Lm = norm(V1-V2);
            Kfb = 1+(df/2)/(Ls+df/2);
            Cn = Kfb*(4*N*(Ls/df)^2)/(1+sqrt(1+(2*Lm/(Lr+Lt))^2));
            val = SIM.parameter('Value',Cn,'Unit','','Description','Fin stability derivative');
        end
                                  
        function val = get.CnBoat(obj)
            df = obj.diameter.Value;
            dn = obj.baseDiameter.Value;
            Cn = 2*((df/dn)^2-1);
            val = SIM.parameter('Value',Cn,'Unit','','Description','Boattail stability derivative');
        end
                                  
        function val = get.CnRocket(obj)
            Cn = obj.CnBoat.Value+obj.CnFin.Value+obj.CnNose.Value;
            val = SIM.parameter('Value',Cn,'Unit','','Description','Rocket base stability derivative');
        end
                                  
        function val = get.Ar(obj)
            ar = pi/4*obj.diameter.Value^2;
            val = SIM.parameter('Value',ar,'Unit','m^2','Description','Rocket base stability derivative');
        end
                                  
        function val = get.Ap(obj)
            an = pi*obj.diameter.Value*(sqrt(obj.diameter.Value^2+obj.noseLength.Value^2));
            ab = pi*obj.diameter.Value*(obj.bodyLength.Value+obj.boatTailLength.Value);
            val = SIM.parameter('Value',an+ab,'Unit','m^2','Description','Rocket base stability derivative');
        end
        
        function val = get.XcpNose(obj)
            if strcmpi(obj.noseShape.Value,'conical')
                X = 2/3*obj.noseLength.Value;
            elseif strcmpi(obj.noseShape.Value,'ogive')
                X = 0.466*obj.noseLength.Value;
            elseif strcmpi(obj.noseShape.Value,'parabolic') || strcmpi(obj.noseShape.Value,'haack')
                X = 1/2*obj.noseLength.Value;
            else
                error('noseShape must be Conical, Ogive, Parabolic, or Haack');
            end
            val = SIM.parameter('Value',X,'Unit','m','Description','Nosecone center of pressure location');
        end
        
        function val = get.XcpFin(obj)
            XLE = obj.finLE.Value;
            Ls = obj.finLength.Value;
            Lr = obj.finRootChord.Value;
            Lt = obj.finTipChord.Value;
            sw = obj.finSweep.Value;
            V1 = [Ls,Ls/cosd(sw)*sind(sw)+Lt/2];
            V2 = [0,Lr/2];
            Lm = norm(V1-V2);
            X = XLE+(Lm*(Lr+2*Lt))/(3*(Lr+Lt))+1/6*(Lr+Lt-(Lr*Lt)/(Lr+Lt));
            val = SIM.parameter('Value',X,'Unit','m','Description','Nosecone center of pressure location');
        end
        
        function val = get.XcpBoat(obj)
            XLE = obj.boatLE.Value;
            Lb = obj.boatTailLength.Value;
            df = obj.diameter.Value;
            dn = obj.baseDiameter.Value;
            X = XLE+Lb/3*(1+(1-(df/dn))/(1-(df/dn)^2));
            val = SIM.parameter('Value',X,'Unit','m','Description','Boattail center of pressure location');
        end
        
        function val = get.XcpBody(obj)
            XLE = obj.noseLength.Value;
            Lb = obj.bodyLength.Value;
            X = XLE+1/2*Lb;
            val = SIM.parameter('Value',X,'Unit','m','Description','Body center of pressure location');
        end
        
        function val = get.initPosVecNose(obj)
            xNose = obj.initPosVecGnd.Value+[0,0,obj.rCMpad_B.Value(1)]';
            val = SIM.parameter('Value',xNose,'Unit','m','Description','Initial nose position represented in the inertial frame');
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
        
        function val = get.finOutline(obj)
            %Returns the points of points in the body frame, relative to
            %the Wing leading edge. The points go clockwise looking down
            %from the surface positive z axis.
            %For 2 symmetric trapazoids, it starts in the center LE
            pts = zeros(obj.numFins.Value*3,5); 
            R = obj.diameter.Value/2;
            LE = obj.finLE.Value;
            sw = obj.finSweep.Value;
            t = obj.finTipChord.Value;
            r = obj.finRootChord.Value;
            L = obj.finLength.Value;
            %  point 1
            pts(1,1) = LE;  pts(4,1) = LE;  pts(7,1) = LE;  pts(10,1) = LE;
            pts(2,1) = R;   pts(5,1) = -R;  pts(9,1) = R;   pts(12,1) = -R;
            %  point 2
            pts(1,2) = LE+L/cosd(sw)*sind(sw);  pts(4,2) = LE+L/cosd(sw)*sind(sw);  
            pts(7,2) = LE+L/cosd(sw)*sind(sw);  pts(10,2) = LE+L/cosd(sw)*sind(sw);  
            pts(2,2) = L+R;  pts(5,2) = -L-R;  pts(9,2) = L+R;  pts(12,2) = -L-R; 
            %  point 3
            pts(1,3) = LE+L/cosd(sw)*sind(sw)+t;  pts(4,3) = LE+L/cosd(sw)*sind(sw)+t;  
            pts(7,3) = LE+L/cosd(sw)*sind(sw)+t;  pts(10,3) = LE+L/cosd(sw)*sind(sw)+t;  
            pts(2,3) = L+R;  pts(5,3) = -L-R;  pts(9,3) = L+R;  pts(12,3) = -L-R; 
            %  point 4
            pts(1,4) = LE+r; pts(4,4) = LE+r; pts(7,4) = LE+r;  pts(10,4) = LE+r;  
            pts(2,4) = R;    pts(5,4) = -R;   pts(9,4) = R;     pts(12,4) = -R; 
            %  point 5
            pts(1,5) = LE;  pts(4,5) = LE;  pts(7,5) = LE;  pts(10,5) = LE;
            pts(2,5) = R;   pts(5,5) = -R;  pts(9,5) = R;   pts(12,5) = -R;

            val = SIM.parameter('Value',pts,'Unit','m','Description',...
                'Fin outer postion points.');
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
            addParameter(p,'noseRings',3,@isnumeric);
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
            orig = obj.rCMpad_B.Value(1);
            %  Plot fins
            for ii = 1:4
                i1 = 3*(ii-1)+1;
                i2 = 3*(ii-1)+3;
                pts = obj.finOutline.Value(i1:i2,:);
                h.surf{ii} = plot3(h.ax,...
                    pts(1,:)-orig,...
                    pts(2,:),...
                    pts(3,:),...
                    'LineWidth',1.2,'Color','k','LineStyle','-',...
                    'DisplayName','Fluid Dynamic Surfaces');
                hold on
            end
            %  Plot body
            x = linspace(obj.noseLength.Value,obj.noseLength.Value+obj.bodyLength.Value,p.Results.fuseRings);
            perSlice = 16;
            x = reshape(repmat(x,perSlice,1),[1 numel(x)*perSlice]);
            th = linspace(0,2*pi,perSlice);
            y = repmat(obj.diameter.Value/2*cos(th),1,p.Results.fuseRings);
            z = repmat(obj.diameter.Value/2*sin(th),1,p.Results.fuseRings);
            numextra = (perSlice-1)*2;
            xend = x(end);
            for i = 0:numextra-1
                if mod(i,4) == 0 || mod(i,4) == 3
                    x(end+1) = x(1);
                else
                    x(end+1) = xend;
                end
            end
            y(end+1:end+numextra) = reshape(repmat(y(2:perSlice),2,1),[1 numextra]);
            z(end+1:end+numextra) = reshape(repmat(z(2:perSlice),2,1),[1 numextra]);
                        
            h.surf{5} = plot3(h.ax,x-orig,y,z,'LineWidth',.2,'Color','k','LineStyle','-',...
                'DisplayName','Fluid Dynamic Surfaces');
            %  Plot nose
            L = obj.noseLength.Value;
            N = p.Results.noseRings;
            r = obj.diameter.Value/2;
            x = linspace(0,L,N);
            th = linspace(0,2*pi,16);
            for i = 1:p.Results.noseRings
                rVec(i) = x(i)/L*r;
                C(i,1) = x(i);  C(i,2) = 0;  C(i,3) = 0;
                X(i,:) = C(i,1)+zeros(size(th))-orig;
                Y(i,:) = C(i,2)+rVec(i)*cos(th);
                Z(i,:) = C(i,3)+rVec(i)*sin(th);
            end
            h.surf{6} = patch(X,Y,Z,'k');
            %  Plot boattail
            xb = [obj.boatLE.Value,obj.boatLE.Value+obj.boatTailLength.Value]-orig;
            rVecb = [r obj.baseDiameter.Value/2];
            for i = 1:2
                C(i,1) = xb(i);  C(i,2) = 0;  C(i,3) = 0;
                Xb(i,:) = C(i,1)+zeros(size(th));
                Yb(i,:) = C(i,2)+rVecb(i)*cos(th);
                Zb(i,:) = C(i,3)+rVecb(i)*sin(th);
            end
            h.surf{7} = patch(Xb,Yb,Zb,'k');
            xCenter = obj.length.Value;
            yCenter = 0;
            zCenter = 0;
            yr = obj.baseDiameter.Value/2 * cos(th) + yCenter;
            zr = obj.baseDiameter.Value/2 * sin(th) + zCenter;
            xr = zeros(1, numel(th)) + xCenter;
            h.surf{8} = plot3(xr-orig, yr, zr, 'k-', 'LineWidth', 1);

            if ~p.Results.Basic                
                % Center of mass
                h.centOfMass = plot3(h.ax,...
                                    obj.rCMpad_B.Value(1)+p.Results.Position(1)-orig,...
                                    obj.rCMpad_B.Value(2)+p.Results.Position(2),...
                                    obj.rCMpad_B.Value(3)+p.Results.Position(3),...
                                    'r*','DisplayName','Center of Mass');
                % Center of mass
                h.centOfPress = plot3(h.ax,...
                                    obj.rCP_B.Value(1)+p.Results.Position(1)-orig,...
                                    obj.rCP_B.Value(2)+p.Results.Position(2),...
                                    obj.rCP_B.Value(3)+p.Results.Position(3),...
                                    'b*','DisplayName','Center of Pressure');
                % Coordinate origin
                legend(h.ax,[h.surf{1} h.centOfMass h.centOfPress],'Location','northwest')
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