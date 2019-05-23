% Example of Figure Rotator
%
% This scipt demonstrates how a FigureRotator can be applied to a figure,
% enabling a user to move around a target easily in 3D space. Only one line
% is necessary (37) to apply the FigureRotator, and it works despite that
% the figure is constantly being updated.
% 
% Additional examples can be found by typing 'help FigureRotator' at the
% command line.
%
% Tucker McClure
% Copyright 2012, The MathWorks, Inc.

% Create the figure.
figure(1)
clf();

% Draw some stuff.
n = 1000;
vertices = sign(randn(3*n, 3)).*randn(3*n, 3).^2;
faces    = reshape(1:3*n, n, 3);
colors   = (1-rand(3*n, 1).^5) * rand(1, 3);
h = patch('Vertices',        vertices, ...
          'Faces',           faces, ...
          'FaceVertexCData', colors);
      
% Set up shading and lighting.
shading interp;
set(1, 'Color', 'w');
set(gca, 'Position', [0 0 1 1], 'Projection', 'Perspective');
axis off;
camlight headlight;
lighting gouraud;

% Apply the FigureRotator. Any inputs beyond the current axes are optional
% and are passed directly to the axes object.
f = FigureRotator(gca, 'CameraTarget',    [0 0 0], ...
                       'CameraPosition',  [12 12 12], ...
                       'CameraViewAngle', 60);
                   

% Make it interesting to look out while we move around.
x1 = randn(3*n, 3);
x2 = randn(3*n, 3);
x3 = 1-rand(3*n, 3).^2;
f1 = rand(3*n, 3);

% Let the user stop the animation by pressing escape.
stop = false;
stop_fcn = @(~, k) assignin('base', 'stop', strcmp(k.Key, 'escape'));
set(1, 'KeyPressFcn', stop_fcn);

h_text = text(10, 10, 10, 'Press Escape to Stop Animation', ...
              'HorizontalAlignment', 'center', ...
              'Color', 'w');

% Loop forever.
tic();
while ~stop && ishandle(1)
    
    t = toc();
    if t > 3 && ishandle(h_text), delete(h_text); end % Delete the text.
    
    vertices = vertices + 0.003*(sin(0.5*t)*x1 + 0.6*sin(0.7*t)*x2);
    a1 = 0.75 + 0.25*sin(0.52*pi*t);
    a2 = 0.66 + 0.34*cos(0.56*t);
    x3 = max(min(x3 + 0.01*sin(f1*t), 1), 0);
    set(h, 'Vertices', vertices, ...
           'FaceVertexCData', a1*(a2.*colors + (1-a2)*x3));
    drawnow();
end
