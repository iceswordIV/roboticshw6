clearvars, clc, close all

syms l0 l1 l2 L1 L2 L3 r1 r2 real
syms theta1 theta2 theta3 theta4 real
h1 = 0.3;
h2 = 0.4;
h3 = 0.5;
% theta1 = 0;
% theta2 = 0;
% theta3 = 0;

q1 = [0;0;h1];
q2 = [0;l1;h2];
q3 = [0;l1 + l2;h3];
q4 = [0;l1 + l2;l0];

w1 = [0;0;1];
w2 = [0;0;1];
w3 = [0;0;1];


%using unit vector for w
w1 = w1/norm(w1);
w2 = w2/norm(w2);
w3 = w3/norm(w3);
w4 = [0;0;0];

twist1 = revolute_twist(w1,q1);
twist2 = revolute_twist(w2,q2);
twist3 = revolute_twist(w3,q3);

v1 = [twist1(1,1);twist1(2,1);twist1(3,1)];
v2 = [twist2(1,1);twist2(2,1);twist2(3,1)];
v3 = [twist3(1,1);twist3(2,1);twist3(3,1)];
v4 = [0;0;1];

twist4 = [v4;w4];


center_point1 = [0;r1;h1];
center_point2 = [0;l1 + r2;h2];
center_point3 = [0;l1 + l2;h3];
center_point4 = [0;l1 + l2;l0];



I = [
    1,0,0;
    0,1,0;
    0,0,1
    ];

gst_1 = gst(q1);
gst_2 = gst(q2);
gst_3 = gst(q3);
gst_4 = gst(q4);


exponential_epison_hat_theta1 = exponential_epison_hat_theta(w1,v1,theta1);
exponential_epison_hat_theta2 = exponential_epison_hat_theta(w2,v2,theta2);
exponential_epison_hat_theta3 = exponential_epison_hat_theta(w3,v3,theta3);
exponential_epison_hat_theta4 = exponential_epison_hat_theta(w4,v4,theta4);





gsl1_0 = gsl_0(center_point1);
gsl2_0 = gsl_0(center_point2);
gsl3_0 = gsl_0(center_point3);
gsl4_0 = gsl_0(center_point4);


gst1 = exponential_epison_hat_theta1 * gsl1_0;
gst2 = exponential_epison_hat_theta1 * exponential_epison_hat_theta2 * gsl2_0;
gst3 = exponential_epison_hat_theta1 * exponential_epison_hat_theta2 * exponential_epison_hat_theta3 * gsl3_0;
gst4 = exponential_epison_hat_theta1 * exponential_epison_hat_theta2 * exponential_epison_hat_theta3 * exponential_epison_hat_theta4 * gsl4_0;

% J1 = [adjmginv(gst1) * twist1, zeros(6,1), zeros(6,1)]
% J2 = [adjmginv(gst2) * twist1, adjmginv(exponential_epison_hat_theta2 * gsl2_0) * twist2, zeros(6,1)]
% J3 = [adjmginv(gst3) * twist1, adjmginv(exponential_epison_hat_theta2 * exponential_epison_hat_theta3 * gsl3_0) * twist2, adjmginv(exponential_epison_hat_theta3 * gsl3_0) * twist3]

% suppose twist1, twist2, twist3 are 6×1 vectors in the spatial frame,
% and you have your exponential‐map function exponential_se3_hat(ω,v,θ)
%
%%

% --- after you compute twist1, twist2, twist3, and gsl1_0, gsl2_0, gsl3_0 ---

% 1) compute the *constant* adjoint-inverses of each COM frame:
Ad1 = adjmginv(exponential_epison_hat_theta1 * gsl1_0 );

E1 = exponential_epison_hat_theta(w1,v1,theta1);
E2 = exponential_epison_hat_theta(w2,v2,theta2);
E3 = exponential_epison_hat_theta(w3,v3,theta3);
E4 = exponential_epison_hat_theta(w4,v4,theta4);

% Link 2: column 1 sees E1·E2, column 2 sees E2 alone
Ad2_1 = adjmginv( E1*E2 * gsl2_0 );
Ad2_2 = adjmginv(      E2 * gsl2_0 );

% Link 3: col 1 sees E1·E2·E3, col 2 sees E2·E3, col 3 sees E3
Ad3_1 = adjmginv( E1*E2*E3 * gsl3_0 );
Ad3_2 = adjmginv(      E2*E3 * gsl3_0 );
Ad3_3 = adjmginv(           E3 * gsl3_0 );

% Link 4: col 1 sees E1·E2·E3·E4, col 2 sees E2·E3·E4,
% col 3 sees E3·E4, col 4 sees E4 alone
Ad4_1 = adjmginv( E1*E2*E3*E4 * gsl4_0 );
Ad4_2 = adjmginv(      E2*E3*E4 * gsl4_0 );
Ad4_3 = adjmginv(           E3*E4 * gsl4_0 );
Ad4_4 = adjmginv(                E4 * gsl4_0 );



% 2) assemble the three 6×3 Jacobians at zero angles:
J1 = [ Ad1 * twist1,  zeros(6,1), zeros(6,1), zeros(6,1)];
J2 = [ Ad2_1 * twist1, Ad2_2 * twist2, zeros(6,1), zeros(6,1)];
J3 = [ Ad3_1 * twist1, Ad3_2 * twist2,  Ad3_3 * twist3, zeros(6,1)];
J4 = [ Ad4_1 * twist1, Ad4_2 * twist2, Ad4_3 * twist3, Ad4_4 * twist4];


J1 = simplify(J1)
J2 = simplify(J2)
J3 = simplify(J3)
J4 = simplify(J4)

% inertia_matrices


%% --- inertial parameters ---
syms m1 m2 m3 m4               % masses
syms I1x I1y I1z            % link 1 principal inertias
syms I2x I2y I2z            % link 2 principal inertias
syms I3x I3y I3z            % link 3 principal inertias
syms I4x I4y I4z
syms dtheta1 dtheta2 dtheta3 dtheta4 % joint velocities

% 6×6 link inertia matrices at each COM
M1 = [ m1*eye(3),           zeros(3);
       zeros(3),      diag([I1x,I1y,I1z]) ];

M2 = [ m2*eye(3),           zeros(3);
       zeros(3),      diag([I2x,I2y,I2z]) ];

M3 = [ m3*eye(3),           zeros(3);
       zeros(3),      diag([I3x,I3y,I3z]) ];

M4 = [ m4*eye(3),           zeros(3);
       zeros(3),      diag([I4x,I4y,I4z]) ];
%% --- inertia matrix M(θ) ---
M = simplify( J1.'*M1*J1 ...
            + J2.'*M2*J2 ...
            + J3.'*M3*J3 ...
            + J4.'*M4*J4 ...
            );

%% --- Coriolis/centrifugal matrix C(θ,θ̇) ---
C = sym(zeros(4,4));
thetas  = [theta1;theta2;theta3;theta4];
dthetas = [dtheta1;dtheta2;dtheta3;dtheta4];


% Christoffel symbols Γ_{i,j,k}
for i = 1:3
  for j = 1:3
    for k = 1:3
      Gamma = 1/2*( ...
        diff(M(i,j), thetas(k)) + ...
        diff(M(i,k), thetas(j)) - ...
        diff(M(j,k), thetas(i)) );
      C(i,j) = C(i,j) + Gamma * dthetas(k);
    end
  end
end

C = simplify(C);

%% --- display results ---
M_str = char(M);
disp(M_str)

C_str = char(C);
disp(C_str)


%% --- numeric parameters ---
symVars = { ...
  l0,l1,l2, r1,r2, ...
  m1,m2,m3,m4, ...
  I1x,I1y,I1z, I2x,I2y,I2z, I3x,I3y,I3z, I4x,I4y,I4z, ...
  theta1,theta2,theta3,theta4, ...
  dtheta1,dtheta2,dtheta3,dtheta4 ...
};
M_fun = matlabFunction(M, 'Vars', symVars);
C_fun = matlabFunction(C, 'Vars', symVars);


%% --- NOW assign your numeric values ---
l0 = 0.2;    l1 = 0.8;    l2 = 0.6;
r1 = l1/2;   r2 = l2/2;

m1=50; m2=30; m3=20; m4=5;
I1x=m1*l1^2/12; I1z=I1x; I1y=0.05*I1x;
I2x=m2*l2^2/12; I2z=I2x; I2y=0.05*I2x;
I3z=0.10*I2z;   I3x=0.20*I3z; I3y=0.20*I3z;
I4z=0.60*I3z;   I4x=0.20*I4z; I4y=0.20*I4z;

theta1=pi/4; theta2=pi/6; theta3=pi/2.5; theta4=-0.4;
dtheta1=0.05; dtheta2=0.10; dtheta3=0.20; dtheta4=0.03;


%% --- finally call them to get numeric M and C ---
args = { ...
  l0,l1,l2, r1,r2, ...
  m1,m2,m3,m4, ...
  I1x,I1y,I1z, I2x,I2y,I2z, I3x,I3y,I3z, I4x,I4y,I4z, ...
  theta1,theta2,theta3,theta4, ...
  dtheta1,dtheta2,dtheta3,dtheta4 ...
};

M_val = M_fun(args{:});
C_val = C_fun(args{:});

disp('Numeric Mass Matrix M ='), disp(M_val)
disp('Numeric Coriolis   C ='), disp(C_val)




