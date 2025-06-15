function ewtheta = exponential(w,theta)

% w = theta * ||w||
% 
% omega_mag = norm(w);
% theta = omega_mag;
% 
% w = w / theta;


% theta is be defined using syms
% sometimes theta is variable and using w be unit vector


w1 = w(1,1);
w2 = w(2,1);
w3 = w(3,1);

w = [
    0, -w3, w2;
    w3, 0, -w1;
    -w2, w1, 0
    ];

I = [1,0,0;0,1,0;0,0,1];

ewtheta = I + w * sin(theta) + w^2 * (1 - cos(theta));
end