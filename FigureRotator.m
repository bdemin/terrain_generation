classdef FigureRotator < handle

% FigureRotator
%
% Apply this to a figure to use a mouse to rotate around a target position
% and zoom in and out. This is similar to the Rotate 3D capability found in
% the standard figure menu, but this allows greater flexibility and
% fluidity of movement. An example is provided in example_figure_rotator.m,
% and additional examples are provided below.
%
% Left-click and drag to rotate about the target.
% Scroll the mouse wheel to move towards/away from the target.
% Right-click and drag to zoom the camera in/out, changing the view angle.
%
% Example:
%
% figure();
% plot3(randn(1, 10), randn(1, 10), randn(1, 10));
% drawnow();
% f = FigureRotator(gca);
%
% The FigureRotator can later be stopped by calling the Stop() function.
%
% f.Stop();
%
% It's often helpful to specify the initial camera parameters, like
% position, target, up vector, and view angle. These can all be passed to
% the constructor.
%
% f = FigureRotator(gca, 'CameraTarget',    [0 0 0], ...
%                        'CameraPosition',  [15 0 0], ...
%                        'CameraUpVector',  [0 0 1], ...
%                        'CameraViewAngle', 60);
%
% The FigureRotator allows complete 3D rotation, so if you start losing
% track of "up", you can always re-align the camera's up vector with the 
% axes' [0 0 1] by calling RestoreUp().
%
% This object uses the figure's WindowButtonUpFcn, WindowButtonDownFcn,
% WindowButtonMotionFcn, and WindowScrollWheelFcn callbacks. If those are
% necessary for other tasks as well, callbacks can be attached to the
% FigureRotator, which will pass all arguments on to the provided callback
% function.
%
% Example:
%
% f = FigureRotator(gca);
% f.AttachCallback('WindowButtonDownFcn', 'disp(''clicked'');');
% 
% Or multiple callbacks can be set with a single call:
%
% f.AttachCallback('WindowButtonDownFcn',   'disp(''down'');', ...
%                  'WindowButtonUpFcn',     'disp(''up'');', ...
%                  'WindowButtonMotionFcn', 'disp(''moving'');', ...
%                  'WindowScrollWheelFcn',  @(~, ~) disp('scrolling'));
%
% A single FigureRotator can control multiple axes, even axes across
% multiple figures.
%
% Example:
% 
% figure(1);
% clf();
% ha1 = subplot(2, 1, 1);
% peaks;
% ha2 = subplot(2, 1, 2);
% peaks;
% 
% figure(2);
% clf();
% peaks;
% ha3 = gca();
% 
% f = FigureRotator([ha1 ha2 ha3]);
%
% Tucker McClure
% Copyright 2012, The MathWorks, Inc.
    
    properties
        
        % Figure and axes handles
        h_f;
        h_a;
        
        % Current states
        rotating = false;
        zooming  = false;

        rotate_start_point = [0 0]; % Mouse position on button down
        zoom_start_point   = [0 0]; % Mouse position on button down
        
        r_t0;       % Position of target wrt 0 on button down
        r_tc;       % Position of target wrt camera on button down
        up_hat;     % Up direction on button down
        right_hat;  % Right direction on button down
        view_angle; % View angle on button down
        
        wbdf; % Pass-through WindowButtonDownFcn
        wbuf; % Pass-through WindowButtonUpFcn
        wbmf; % Pass-through WindowButtonMotionFcn
        wswf; % Pass-through WindowScrollWheelFcn
        
    end
    
    methods
        
        % Construct a FigureRotator for the given axes.
        function o = FigureRotator(axes_handle, varargin)

            % Record the axes and figure.
            o.h_a = axes_handle;
            o.h_f = get(axes_handle, 'Parent');
            if iscell(o.h_f)
                o.h_f = cell2mat(o.h_f);
            end

            % Pass any arguments on to the axes object.
            if nargin >= 2
                set(o.h_a, varargin{:});
            end

            % Set the figure callbacks to register with this object.
            set(o.h_f, ...
                'WindowButtonDownFcn',   @o.ButtonDown, ...
                'WindowButtonUpFcn',     @o.ButtonUp, ...
                'WindowButtonMotionFcn', @o.Move, ...
                'WindowScrollWheelFcn',  @o.Wheel);
            
            % Set up the axes object for what we need. We get the last
            % word.
            set(o.h_a, ...
                'CameraPositionMode',  'manual', ...
                'CameraTargetMode',    'manual', ...
                'CameraUpVectorMode',  'manual', ...
                'CameraViewAngleMode', 'manual', ...
                'XLimMode',            'manual', ...
                'YLimMode',            'manual', ...
                'ZLimMode',            'manual', ...
                'DataAspectRatioMode', 'manual');
            
        end
        
        % Called when a user clicks
        function ButtonDown(o, h, event, varargin)

            % Get the button type.
            switch get(h, 'SelectionType')
                
                % Rotate around.
                case {'normal', 'extend'}

                    % Record the starting point and that we're rotating.
                    o.rotate_start_point = get(gcf(), 'CurrentPoint');
                    o.rotating           = true;

                    % Get positions (and account for non-equal aspect 
                    % ratios).
                    o.r_t0 = get(gca(), 'CameraTarget')';
                    r_c0   =    get(gca(), 'CameraPosition')' ...
                             ./ get(gca(), 'DataAspectRatio')';
                    o.r_tc = o.r_t0 - r_c0;

                    % Calculate "up" and "right" for camera.
                    o.up_hat    = get(gca(), 'CameraUpVector')';
                    o.right_hat = normalize(crs(o.r_tc)*o.up_hat);
                    o.up_hat    = normalize(crs(o.right_hat)*o.r_tc);
                
                % Zoom.
                case 'alt'
                
                    % Record the starting point and that we're zooming.
                    o.zoom_start_point = get(gcf(), 'CurrentPoint');
                    o.zooming          = true;
                    
                    % Record the current view angle.
                    o.view_angle = get(gca(), 'CameraViewAngle');
                    
            end

            % If there's a callback attachment, execute it.
            execute_callback(o.wbdf, h, event, varargin{:});
            
        end
        
        % Called when user releases a click
        function ButtonUp(o, h, event, varargin)
            
            % Get the button type.
            switch get(h, 'SelectionType')
                
                % Stop rotating.
                case {'normal', 'extend'}
                    o.rotating = false;
                    
                % Stop zooming.
                case 'alt'
                    o.zooming = false;
                    
            end
            
            % If there's a callback attachment, execute it.
            execute_callback(o.wbuf, h, event, varargin{:});
            
        end
        
        % Called when mouse moves in figure
        function Move(o, h, event, varargin)
            
            if o.rotating
                
                % Get the mouse position in the window.
                s = feval(@(x) x(3:4), get(h, 'Position'));
                p = get(gcf(), 'CurrentPoint');
                r = (p - o.rotate_start_point)./s;
                   
                % Calculate where the mouse is.
                r_mc = r(2)*o.up_hat + r(1)*o.right_hat;
                
                % Calculate the rotation axis.
                a_hat = normalize(cross(o.r_tc - r_mc, o.r_tc));
                
                % Calculate the rotation matrix.
                Q = aa2dcm(a_hat, min(norm(r_mc), 1)*pi);
                
                % Calculate the new camera position, accounting for
                % non-equal aspect ratios.
                r_n0 = -Q*o.r_tc + o.r_t0;
                r_n0 = r_n0 .* get(gca(), 'DataAspectRatio')';
                
                % Update the relevant quantities.
                set(gca(), 'CameraPosition', r_n0, ...
                           'CameraUpVector', Q*o.up_hat);
                          
            end
            
            if o.zooming
                
                % Get the mouse position in the window.
                s = feval(@(x) x(4), get(h, 'Position'));
                p = feval(@(x) x(2), get(gcf(), 'CurrentPoint'));
                r = -(p - o.zoom_start_point(2))/s;
                new_view_angle = min(2^r*o.view_angle, 180-eps);
                set(gca(), 'CameraViewAngle', new_view_angle);
                
            end
            
            % If there's a callback attachment, execute it.
            execute_callback(o.wbmf, h, event, varargin{:});
            
        end
        
        % Called for scroll wheel
        function Wheel(o, h, event, varargin)
            
            % Scalar to increase/decrease distance to target.
            s = 1.2^event.VerticalScrollCount;
            
            % Update the stored vector by increasing/decreasing its length
            % (in case we're moving).
            o.r_tc = s * o.r_tc;
            
            % Update what we're currently seeing by the appropriate amount.
            t0 = get(gca(), 'CameraTarget');
            c0 = get(gca(), 'CameraPosition'); 
            r_n0 = s * (c0 - t0) + t0;
            set(gca(), 'CameraPosition', r_n0);
            
            % If there's a callback attachment, execute it.
            execute_callback(o.wswf, h, event, varargin{:});
            
        end
        
        % Sometime users like to return "up" to [0 0 1], so we'll give them
        % a function to call.
        function RestoreUp(o)
            set(o.h_a, 'CameraUpVector', [0 0 1]);
        end

        % Add a pass-through callback for one of the callbacks
        % FigureRotator hogs to itself. This way, a user can still get all
        % the info he needs from a figure's callbacks *and* use the 
        % rotator.
        function AttachCallback(o, varargin)
            
            for k = 2:2:length(varargin)
                switch varargin{k-1}
                    case 'WindowButtonDownFcn'
                        o.wbdf = varargin{k};
                    case 'WindowButtonUpFcn'
                        o.wbuf = varargin{k};
                    case 'WindowButtonMotionFcn'
                        o.wbmf = varargin{k};
                    case 'WindowScrollWheelFcn'
                        o.wswf = varargin{k};
                    otherwise
                        warning('Invalid callback attachment.');
                end
            end
            
        end
        
        % We're done. Get rid of the callbacks. If there were pass-through
        % callbacks, replace our callbacks with those.
        function Stop(o)
            set(o.h_f, ...
                'WindowButtonDownFcn',   o.wbdf, ...
                'WindowButtonUpFcn',     o.wbuf, ...
                'WindowButtonMotionFcn', o.wbmf, ...
                'WindowScrollWheelFcn',  o.wswf);
        end
        
    end
    
end

% Safely normalize an input vector.
function x_hat = normalize(x)
    n = norm(x);
    if n > eps
        x_hat = x/n;
    else
        x_hat = x;
    end
end

% Convert the specified axis and angle of rotation to a direction cosine
% matrix.
function M = aa2dcm(ax, an)
    M = eye(3)*cos(an) + (1-cos(an))*(ax*ax') + crs(ax)*sin(an);
end

% Returns a skew-symmetric "cross product" matrix from 3-by-1 vector, v,
% such that cross(v, b) == crs(v)*b.
function M = crs(v)
    M = [ 0    -v(3)  v(2); ...
          v(3)  0    -v(1); ...
         -v(2)  v(1)  0];
end

% Execute whatever callback was requested.
function execute_callback(cb, h, event, varargin)
    if ~isempty(cb)
        if isa(cb, 'function_handle')
            cb(h, event)
        elseif iscell(cb)
            cb(h, event, varargin{:})
        else
            eval(cb);
        end
    end
end


