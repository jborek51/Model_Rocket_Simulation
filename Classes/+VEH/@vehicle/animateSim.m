function animateSim(obj,tsc,timeStep,varargin)
%ANIMATESIM Method to animate a simulation using the provided timeseries

%% Input parsing
p = inputParser;


% ---Fundamental Animation Requirements---
% Timeseries collection structure with results from the simulation
addRequired(p,'tsc',@(x) or(isa(x,'signalcontainer'),isa(x,'struct')));
% Time step used in plotting
addRequired(p,'timeStep',@isnumeric);
% Time to start viewing
addParameter(p,'startTime',0,@isnumeric);
% Time to start viewing
addParameter(p,'endTime',tsc.positionVec.Time(end),@isnumeric);

% Vector of time stamps used to crop data
addParameter(p,'CropTimes',[],@isnumeric)

% ---Parameters for saving a gif---
% Switch to enable saving 0 = don't save
addParameter(p,'SaveGif',false,@islogical)
% Path to saved file, default is ./output
addParameter(p,'GifPath',fullfile(fileparts(which('Model_Rocket_Simulation.prj')),'Results'));
% Name of saved file, default is animation.gif
addParameter(p,'GifFile','animation.gif');
% Time step between frames of gif, default is time step (real time plot)
addParameter(p,'GifTimeStep',timeStep,@isnumeric)

% ---Parameters to save MPEG---
addParameter(p,'SaveMPEG',false,@islogical) % Boolean switch to save a MPEG
addParameter(p,'MPEGPath',fullfile(fileparts(which('OCTProject.prj')),'output'));
addParameter(p,'MPEGFile','animation');

% ---Plot Features---
% X limits on the plot
addParameter(p,'XLim',[],@isnumeric)
% Y limits on the plot
addParameter(p,'YLim',[],@isnumeric)
% Z limits on the plot
addParameter(p,'ZLim',[],@isnumeric)
% Name of the path geometry used
addParameter(p,'PathFunc',[],@ischar);
% Plot ground coordinate system axes
addParameter(p,'PlotAxes',true,@islogical);
% Set camera view angle [azimuth, elevation]
addParameter(p,'View',[71,22],@isnumeric)
% Set camera view angle [azimuth, elevation]
addParameter(p,'FigPos',[488 342 560 420],@isnumeric)
% Set font size
addParameter(p,'FontSize',get(0,'defaultAxesFontSize'),@isnumeric)
% Tracer (streaming red line behind the model)
addParameter(p,'PlotTracer',true,@islogical)
% Plot the ground station
addParameter(p,'GroundStation',[],@(x) isa(x,'OCT.sixDoFStation'))
% Plot the glider
addParameter(p,'Glider',[],@(x) isa(x,'OCT.vehicle'))
% Color tracer according to power production/consumption
addParameter(p,'ColorTracer',false,@islogical)
% Change Color tracer variable structure
% Must have properties: timesignal, min, max, minColor, and maxColor
% the timesignal must be singular in the non time dimention.
% the colors should be 3 by 1 vectors with values from 0 to 1.
addParameter(p,'ColorTracerVariableStruct',false)
% How long (in seconds) to keep the tracer on for
addParameter(p,'TracerDuration',5,@isnumeric)
% Plot a red dot on the closest point on the path
addParameter(p,'PathPosition',false,@islogical)
% Plot normal, tangent and desired vectors
addParameter(p,'NavigationVecs',false,@islogical)
% Plot tether nodal net force vecs
addParameter(p,'TetherNodeForces',false,@islogical)
% Plot velocity vector
addParameter(p,'VelocityVec',false,@islogical)
% Plot local aerodynamic force vectors on surfaces
addParameter(p,'LocalAero',false,@islogical)
% Add resulting net moment in the body frame to the table readout
addParameter(p,'FluidMoments',false,@islogical)
% Pause after each plot update (to go frame by frame)
addParameter(p,'Pause',false,@islogical)
% Zoom in the plot axes to focus on the body
addParameter(p,'ZoomIn',false,@islogical)
% Colorbar on right showing iteration power
addParameter(p,'ZoomInMove',false,@islogical)
% Colorbar on right showing iteration power
addParameter(p,'PowerBar',false,@islogical)
% Plot the tangent coordinate system
addParameter(p,'TangentCoordSys',false,@islogical);
% Optional scrolling plots on the side
addParameter(p,'ScrollPlots',{}, @(x) isa(x,'cell') && all(isa([x{:}],'timeseries'))); % Must be a cell array of timeseries objects
% Plot bedrock or not
addParameter(p,'Bedrock',true,@islogical)
% Plot bedrock or not
addParameter(p,'LineAngleEst',false,@islogical)
% Plot Flow Velocity Vector
addParameter(p,'FlowVec',false,@islogical)

% ---Parse the output---
parse(p,tsc,timeStep,varargin{:})


%% Setup some infrastructure type things
% If the user wants to save something and the specified directory does not
% exist, create it
if p.Results.SaveGif && ~exist(p.Results.GifPath, 'dir')
    mkdir(p.Results.GifPath)
end
if p.Results.SaveMPEG && ~exist(p.Results.MPEGPath, 'dir')
    mkdir(p.Results.MPEGPath)
end

if p.Results.SaveMPEG
    vidWriter = VideoWriter(fullfile(p.Results.MPEGPath,p.Results.MPEGFile), 'MPEG-4');
    open(vidWriter);
end

% Crop to the specified times
tscTmp = tsc.crop(p.Results.startTime,p.Results.endTime);
% Resample the timeseries to the specified framerate
tscTmp = tscTmp.resample(p.Results.timeStep);

%% Plot things
% Plot the aerodynamic surfaces
h = obj.plot('Basic',true);
h.ax = gca;
hold on

axes(h.ax)
% Get the "nominal" positions of the aerodynamic surfaces from that plot
for ii = 1:length(h.surf)
    hStatic{ii}.x = h.surf{ii}.XData;
    hStatic{ii}.y = h.surf{ii}.YData;
    hStatic{ii}.z = h.surf{ii}.ZData;
end
hold on

% Plot x, y and z ground fixed axes
if p.Results.PlotAxes
    posData = squeeze(tscTmp.positionVec.Data)';
    r = sqrt(sum(posData.^2,2));
    len = 0.1*max(r);
    plot3([0 len],[0 0],[0 0],...
        'Color','r','LineStyle','-');
    plot3([0 0],[0 len],[0 0],...
        'Color','g','LineStyle','-');
    plot3([0 0],[0 0],[0 len],...
        'Color','b','LineStyle','-');
end

% Set the plot limits to zoom in on the body
if p.Results.ZoomIn
    xlim(tscTmp.positionVec.Data(1,:,1)+obj.length.Value*[-.75 .75])
    ylim(tscTmp.positionVec.Data(2,:,1)+obj.length.Value*[-.75 .75])
    zlim(tscTmp.positionVec.Data(3,:,1)+obj.length.Value*[-.75 .75])
end

% set(gcf, 'Units', 'Normalized', 'OuterPosition', [0 0 1 1]);
% Get rid of tool bar and pulldown menus that are along top of figure.
% set(gcf, 'Toolbar', 'none', 'Menu', 'none');

if p.Results.VelocityVec
    h.velVec = plot3(...
        [tscTmp.positionVec.Data(1,:,1) tscTmp.positionVec.Data(1,:,1)+tscTmp.velocityVec.Data(1,:,1)*obj.length.Value./norm(tscTmp.velocityVec.Data(:,:,1))],...
        [tscTmp.positionVec.Data(2,:,1) tscTmp.positionVec.Data(2,:,1)+tscTmp.velocityVec.Data(2,:,1)*obj.length.Value./norm(tscTmp.velocityVec.Data(:,:,1))],...
        [tscTmp.positionVec.Data(3,:,1) tscTmp.positionVec.Data(3,:,1)+tscTmp.velocityVec.Data(3,:,1)*obj.length.Value./norm(tscTmp.velocityVec.Data(:,:,1))],...
        'Color','k','LineWidth',1.5,'LineStyle','--');
end

% Set the font size
set(gca,'FontSize',p.Results.FontSize);

% Set the viewpoint
view(p.Results.View)

% Set plot limits
setLimsToQuartSphere(gca,squeeze(tscTmp.positionVec.Data)',...
    'PlotAxes',true);

% Attempt to set plot axes limits automatically
allPlots = allchild(gca);
% Find min and maxes over position
minX = min(tscTmp.positionVec.Data(1,:));
maxX = max(tscTmp.positionVec.Data(1,:));
minY = min(tscTmp.positionVec.Data(2,:));
maxY = max(tscTmp.positionVec.Data(2,:));
minZ = min(tscTmp.positionVec.Data(3,:));
maxZ = max(tscTmp.positionVec.Data(3,:));
% Find min and max over all plotted data
for ii = 1:numel(allPlots)
    minX = min([minX allPlots(ii).XData(:)']);
    maxX = max([maxX allPlots(ii).XData(:)']);
    minY = min([minY allPlots(ii).YData(:)']);
    maxY = max([maxY allPlots(ii).YData(:)']);
    minZ = min([minZ allPlots(ii).ZData(:)']);
    maxZ = max([maxZ allPlots(ii).ZData(:)']);
end
% If one is not zero, make X and Y symmetric
xlim([minX maxX+5])
YlimVal = max(abs(minY),abs(maxY));
ylim([-YlimVal-5 YlimVal+5])
if abs(minZ)>abs(maxZ)
    zlim([minZ-5 maxZ])
else
    zlim([minZ maxZ+5])
end


% Set the custom x, y and z limits
if ~isempty(p.Results.XLim)
    xlim(p.Results.XLim)
end
if ~isempty(p.Results.YLim)
    ylim(p.Results.YLim)
end
if ~isempty(p.Results.ZLim)
    zlim(p.Results.ZLim)
end

% Set data aspect ratio to realistic (not skewed)
daspect([1 1 1])
% Set figure position 
h.fig.Position = p.Results.FigPos;
% Create a title
h.title = title({strcat(sprintf('Time = %.1f s',0),',',...
    sprintf(' Speed = %.1f m/s',norm(tscTmp.velocityVec.Data(:,:,1))))});

%% Update the plots
for ii = 1:numel(tscTmp.positionVec.Time)
%     timeStamp = tscTmp.positionVec.Time(ii);
    eulAngs   = tscTmp.eulerAngles.getsamples(ii).Data;
    posVec    = tscTmp.positionVec.getsamples(ii).Data;
    
    for jj = 1:numel(hStatic)
        % Rotate and translate all aero surfaces
        pts = rotation_sequence(eulAngs)...
            *[...
            hStatic{jj}.x(:)';...
            hStatic{jj}.y(:)';...
            hStatic{jj}.z(:)']+...
            tscTmp.positionVec.Data(:,:,ii);
        % Update the Vehicle outline
        h.surf{jj}.XData = pts(1,:);
        h.surf{jj}.YData = pts(2,:);
        h.surf{jj}.ZData = pts(3,:);
    end
    
    % Update the tracer
    if p.Results.PlotTracer
        delete(h.tracer(1));
        
        newLine = line(...
            [h.tracer(end).XData(end) posVec(1)],...
            [h.tracer(end).YData(end) posVec(2)],...
            [h.tracer(end).ZData(end) posVec(3)],...
            'Color',0.5*[1 1 1],'LineWidth',2);
        
        h.tracer(end).XData(end+1) = newLine.XData(1);
        h.tracer(end).YData(end+1) = newLine.YData(1);
        h.tracer(end).ZData(end+1) = newLine.ZData(1);

        h.tracer = [h.tracer(2:end) newLine];
        uistack(h.tracer(end),'top');
    end    
    
    % Update the title
    h.title.String = {strcat(...
        sprintf('Time = %.1f s',tscTmp.velocityVec.Time(ii)),',',...
        sprintf(' Speed = %.1f m/s',norm(tscTmp.velocityVec.getsamples(ii).Data)))};
            
    if p.Results.VelocityVec
        pt = tscTmp.positionVec.getsamples(ii).Data;
        velVec = tscTmp.velocityVec.getsamples(ii).Data;
        speed = sqrt(sum(velVec.^2));
        
        h.velVec.XData = [pt(1) pt(1)+velVec(1)*obj.length.Value./speed];
        h.velVec.YData = [pt(2) pt(2)+velVec(2)*obj.length.Value./speed];
        h.velVec.ZData = [pt(3) pt(3)+velVec(3)*obj.length.Value./speed];
        
    end
    
    % Set the plot limits to zoom in on the body
    if p.Results.ZoomIn
        pt = tscTmp.positionVec.getsamples(ii).Data;
        xlim(pt(1)+obj.length.Value*[-.75 .75])
        ylim(pt(2)+obj.length.Value*[-.75 .75])
        zlim(pt(3)+obj.length.Value*[-.75 .57])
    end
    
    if p.Results.ZoomInMove
        pt = tscTmp.positionVec.getsamples(ii).Data;
        [~,b] = size(tscTmp.gndStnPositionVec.Data);
        if b == 1
            pt2 = tscTmp.gndStnPositionVec.Data;
        else
            pt2 = tscTmp.gndStnPositionVec.getsamples(ii).Data;
        end
        r = pt-pt2;
        midpt =(pt+pt2)/2;
        limbound = norm(r)*[-1 1];
        xlim(midpt(1)+limbound)
        ylim(midpt(2)+limbound)
        zlim(midpt(3)+limbound)
    end
    
    % Update scrolling plots
    if ~isempty(p.Results.ScrollPlots)
        for jj = 1:numel(h.TimeLine)
            h.TimeLine(jj).XData = tscTmp.positionVec.Time(ii)*[1 1];
        end
    end
    
    drawnow
    
    %% Save gif of results
    if p.Results.SaveGif
        frame       = getframe(h.fig);
        im          = frame2im(frame);
        [imind,cm]  = rgb2ind(im,256);
        if ii == 1
            imwrite(imind,cm,fullfile(p.Results.GifPath,p.Results.GifFile),'gif', 'Loopcount',inf);
        else
            imwrite(imind,cm,fullfile(p.Results.GifPath,p.Results.GifFile),'gif','WriteMode','append','DelayTime',p.Results.GifTimeStep)
        end
    end
    
    % Save gif of results
    if p.Results.SaveMPEG
        frame = getframe(gcf);
        writeVideo(vidWriter,frame)
    end
    if p.Results.Pause
        pause
    end
end

if p.Results.SaveMPEG
    close(vidWriter);
end

end


