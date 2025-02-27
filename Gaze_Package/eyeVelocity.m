function [speed, acceleration] = eyeVelocity(gazeError,framesPerSec)
% calculates frame-by-frame gaze velocity
% gazeError has to be an n * 2 matrix, with rows representing frames in
% time and columns representing x,y positions

    assert(size(gazeError,2) == 2);


    x = gazeError(:,1);
    y = gazeError(:,2);
    dt = 1/framesPerSec;%sec
    
    % position changes
    dx = diff(x);
    dy = diff(y);

    
    % velocity components
    vx = dx ./ dt;
    vy = dy ./ dt;

    % velocity magnitude
    speed = sqrt(vx.^2 + vy.^2);

    % Calculate the acceleration of eye movements
    ax = diff(vx) ./ dt;
    ay = diff(vy) ./ dt;
    acceleration = sqrt(ax.^2 + ay.^2);


end