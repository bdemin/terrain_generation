%% Animate Scene
% 
% Caution! This script may not work consistently (or at all) on all
% machines and may take a very long time to run. Please see below for more.
%
% This script will cause the camera to "fly" around a generated scene as 
% time changes from early morning to late at night. It requires that a 
% scene is already generated with example_scene.m and and uses variables 
% that should be in the workspace.
%
% This script worked well on two of the author's machines, but is not
% likely to work on all machines, as rendering performance varies from
% machine to machine, depending on OpenGL drivers, MATLAB version, etc.
%
% If it *does* work, note it is likely to take longer than 10 minutes to
% complete and has been observed to run for more than an hour!
%
% Tucker McClure
% Copyright 2012, The MathWorks, Inc.

fprintf('Creating movie...\n');
tic();

% Select the figure so it's visible for the user in case it's hidden.
figure(1);
set(1, 'Position', [10 50 1920 1080], ...  % 1080p resolution
       'PaperPositionMode',     'auto', ...
       'WindowButtonDownFcn',   [], ...
       'WindowButtonUpFcn',     [], ...
       'WindowButtonMotionFcn', [], ...
       'WindowScrollWheelFcn',  []);

% Switch to a software OpenGL renderer.
% original_opengl = opengl();
% opengl('software');

% Specify how many frames we want.
n_frames = 30 * 30; % x seconds * 30 frames/second

% Add an Output directory if one doesn't exist.
if ~exist('Output', 'dir')
    mkdir('Output');
end

% Rotation functions
R_z = @(theta) [cos(theta) -sin(theta) 0; sin(theta) cos(theta) 0; 0 0 1];
R_y = @(theta) [cos(theta) 0 sin(theta); 0 1 0; -sin(theta) 0 cos(theta)];

% Create a VideoWriter names full_day_<datestamp>.mp4 and open it.
video_name = ['full_day_' datestr(now(), 'yyyy-mm-dd-HH-MM-SS') '.mp4'];
vw = VideoWriter(['Output' filesep video_name], 'MPEG-4');
vw.Quality = 90;
open(vw);

% Loop over time from morning to night.
for k = 1:n_frames
    
    % Sweep from t = 0.1 to 0.9.
    t = 0.8*k/n_frames + 0.1;
    
    % Put the sun in the right place.
    lightangle(h_light, 90, 360*(t - 0.25));
    
    % Set the sun color.
    sun_color = sun_tones(t);
    set(h_light, 'Color', sun_color);
    set(gca(), 'AmbientLightColor', sun_color);
    
    % Move the camera around.
    az = 2*pi/3 * (sin(  pi*k/n_frames - pi/2) + 1);
    el = pi/20  * (sin(2*pi*k/n_frames - pi/2) + 1) + 0.05;
    camera_position = camera_target' + R_z(-az) * R_y(el) * [-1 0 0]';
    set(gca, 'CameraPosition', camera_position', ...
             'CameraUpVector', R_z(-az) * R_y(el) * [0 0 1]');
    
    % Capture a frame and add it to the video.
    writeVideo(vw, getframe(1));
    
end

% Close the video file.
close(vw);

% Return the original renderer.
% opengl(original_opengl);

fprintf('Done.\n');
toc();
