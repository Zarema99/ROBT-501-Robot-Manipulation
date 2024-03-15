%%
%Exercise 1
disp('----------Exercise1----------')

syms k_x k_y k_z theta

assume(theta, 'real')
%Assume K is a versor
assume(k_x == 1/sqrt(3))
assume(k_y == 1/sqrt(3))
assume(k_z == 1/sqrt(3))

%Rotation about a generic axis K passing through the origin
R_K = [k_x*k_x*(1-cos(theta)) + cos(theta), k_x*k_y*(1-cos(theta)) - k_z*sin(theta), k_x*k_z*(1-cos(theta)) + k_y*sin(theta);
       k_x*k_y*(1-cos(theta)) + k_z*sin(theta), k_y*k_y*(1-cos(theta)) + cos(theta), k_y*k_z*(1-cos(theta)) - k_x*sin(theta);
       k_x*k_z*(1-cos(theta)) - k_y*sin(theta), k_y*k_z*(1-cos(theta)) + k_x*sin(theta), k_z*k_z*(1-cos(theta)) + cos(theta)];
   
disp("Dot products:");
X_dot_Y = simplify(dot(R_K(:, 1), R_K(:, 2)));
disp(X_dot_Y);
X_dot_Z = simplify(dot(R_K(:, 1), R_K(:, 3)));
disp(X_dot_Z);
Y_dot_Z = simplify(dot(R_K(:, 2), R_K(:, 3)));
disp(Y_dot_Z);

disp("Length of vectors:");
X_length = simplify(norm(R_K(:, 1)));
disp(X_length)
Y_length = simplify(norm(R_K(:, 2)));
disp(Y_length)
Z_length = simplify(norm(R_K(:, 3)));
disp(Z_length)

%Plot with theta = 30 deg
k_x = 1/sqrt(3);
k_y = 1/sqrt(3);
k_z = 1/sqrt(3);
theta = deg2rad(30);

%Rotation about a generic axis K passing through the origin
R_K_30 = [k_x*k_x*(1-cos(theta)) + cos(theta), k_x*k_y*(1-cos(theta)) - k_z*sin(theta), k_x*k_z*(1-cos(theta)) + k_y*sin(theta);
       k_x*k_y*(1-cos(theta)) + k_z*sin(theta), k_y*k_y*(1-cos(theta)) + cos(theta), k_y*k_z*(1-cos(theta)) - k_x*sin(theta);
       k_x*k_z*(1-cos(theta)) - k_y*sin(theta), k_y*k_z*(1-cos(theta)) + k_x*sin(theta), k_z*k_z*(1-cos(theta)) + cos(theta)];

figure(1)
%Fixed frame
quiver3(0,0,0,1,0,0, 'k')
text(1,0,0, 'X')
hold on
quiver3(0,0,0,0,1,0, 'k')
text(0,1,0, 'Y')
quiver3(0,0,0,0,0,1, 'k')
text(0,0,1, 'Z')

%Rotation axis
quiver3(0,0,0,1/sqrt(3),1/sqrt(3),1/sqrt(3), "r")
text(1/sqrt(3), 1/sqrt(3),1/sqrt(3), 'K')

%Moving frame (calculated with theta = 30 deg)
quiver3(0,0,0,0.9107,0.3333,-0.2440, 'b')
text(0.9107,0.3333,-0.2440, "X'")
quiver3(0,0,0,-0.2440,0.9107,0.3333, 'b')
text(-0.2440, 0.9107, 0.3333, "Y'")
quiver3(0,0,0,0.3333, -0.2440, 0.9107, 'b')
text(0.3333, -0.2440, 0.9107, "Z'")

title('Exercise 1')
%%
%Exercise 2
disp('----------Exercise2----------')

syms alpha beta gamma 
assume(alpha,'real')
assume(beta,'real')
assume(gamma,'real')

%Principal rotations
R_X = [1 0 0; 0 cos(alpha) -sin(alpha); 0 sin(alpha) cos(alpha)];
R_Y = [cos(beta) 0 sin(beta); 0 1 0; -sin(beta) 0 cos(beta)];
R_Z = [cos(gamma) -sin(gamma) 0; sin(gamma) cos(gamma) 0; 0 0 1];

%Post-multiplication rule
R_XYZ = R_X*R_Y*R_Z;

disp("Dot products:");
X_dot_Y = simplify(dot(R_XYZ(:, 1), R_XYZ(:, 2)));
disp(X_dot_Y);
X_dot_Z = simplify(dot(R_XYZ(:, 1), R_XYZ(:, 3)));
disp(X_dot_Z);
Y_dot_Z = simplify(dot(R_XYZ(:, 2), R_XYZ(:, 3)));
disp(Y_dot_Z);

disp("Length of vectors:");
X_length = simplify(norm(R_XYZ(:, 1)));
disp(X_length)
Y_length = simplify(norm(R_XYZ(:, 2)));
disp(Y_length)
Z_length = simplify(norm(R_XYZ(:, 3)));
disp(Z_length)

%Plot with alpha = 30 deg, beta = 45 deg, gamma = 60 deg
alpha = deg2rad(30);
beta = deg2rad(45);
gamma = deg2rad(60);

R_X = [1 0 0; 0 cos(alpha) -sin(alpha); 0 sin(alpha) cos(alpha)];
R_Y = [cos(beta) 0 sin(beta); 0 1 0; -sin(beta) 0 cos(beta)];
R_Z = [cos(gamma) -sin(gamma) 0; sin(gamma) cos(gamma) 0; 0 0 1];

%Post-multiplication rule
R_XY = R_X*R_Y;
R_xyz = R_X*R_Y*R_Z;

figure(2)
%Fixed frame
quiver3(0,0,0,1,0,0, 'k')
text(1,0,0, 'X')
hold on
quiver3(0,0,0,0,1,0, 'k')
text(0,1,0, 'Y')
quiver3(0,0,0,0,0,1, 'k')
text(0,0,1, 'Z')

%Rotation about X axis (X' coincides with X)
%calculated with alpha = 30 deg
quiver3(0,0,0,1,0,0, 'b')
text(1,0,0, "X, X'")
quiver3(0,0,0,0, cos(alpha), sin(alpha), 'b')
text(0, cos(alpha), sin(alpha), "Y'")
quiver3(0,0,0,0, -sin(alpha), cos(alpha), 'b')
text(0, -sin(alpha), cos(alpha), "Z'")

%Rotation about Y' axis (Y'' coincides with Y')
%calculated with alpha = 30 deg, beta = 45 deg
quiver3(0,0,0,0.7071, 0.3536, -0.6124, 'm')
text(0.7071, 0.3536, -0.6124, "X''")
quiver3(0,0,0,0, 0.8660, 0.5, 'm')
text(0, 0.8660, 0.5, "Y', Y''")
quiver3(0,0,0,0.7071, -0.3536, 0.6124, 'm')
text(0.7071, -0.3536, 0.6124, "Z''")

%Rotation about Z'' axis (Z''' coincides with Z'')
%calculated with alpha = 30 deg, beta = 45 deg, gamma = 60 deg
quiver3(0,0,0,0.3536, 0.9268, 0.1268, 'r')
text(0.3536, 0.9268, 0.1268, "X'''")
quiver3(0,0,0,-0.6124,0.1268,0.7803, 'r')
text(-0.6124,0.1268,0.7803, "Y'''")
quiver3(0,0,0,0.7071, -0.3536, 0.6124, 'r')
text(0.7071, -0.3536, 0.6124, "Z'', Z'''")


title('Exercise 2')
%%
%Exercise 3
disp('----------Exercise3----------')

alpha = deg2rad(60);
beta = deg2rad(45);

%Principal rotations
R_X = [1 0 0; 0 cos(alpha) -sin(alpha); 0 sin(alpha) cos(alpha)];
R_Y = [cos(beta) 0 sin(beta); 0 1 0; -sin(beta) 0 cos(beta)];

%Compound rotation (post-multiplication)
R_compound = R_X*R_Y;

disp("Theta:");
theta = rad2deg(acos((trace(R_compound) - 1)/2));
disp(theta)

disp("Screw axis:");
K = 1/(2*sin(deg2rad(theta)))*[R_compound(3,2) - R_compound(2,3); R_compound(1,3) - R_compound(3,1); R_compound(2,1) - R_compound(1,2)];
disp(K)

disp("Length of a screw axis:")
K_length = norm(K);
disp(K_length)

%Plot
figure(3)
%Fixed frame
quiver3(0,0,0,1,0,0, 'k')
text(1,0,0, 'X')
hold on
quiver3(0,0,0,0,1,0, 'k')
text(0,1,0, 'Y')
quiver3(0,0,0,0,0,1, 'k')
text(0,0,1, 'Z')

%Rotation about X axis (X' coincides with X)
%calculated with alpha = 60 deg
quiver3(0,0,0,1,0,0, 'b')
text(1,0,0, "X, X'")
quiver3(0,0,0,0, cos(alpha), sin(alpha), 'b')
text(0, cos(alpha), sin(alpha), "Y'")
quiver3(0,0,0,0, -sin(alpha), cos(alpha), 'b')
text(0, -sin(alpha), cos(alpha), "Z'")

%Rotation about Y' axis (Y'' coincides with Y')
%calculated with alpha = 60 deg, beta = 45 deg
quiver3(0,0,0,0.7071, 0.6124, -0.3536, 'm')
text(0.7071, 0.6124, -0.3536, "X''")
quiver3(0,0,0,0, 0.5, 0.8660, 'm')
text(0, 0.5, 0.8660, "Y', Y''")
quiver3(0,0,0,0.7071, -0.6124, 0.3536, 'm')
text(0.7071, -0.6124, 0.3536, "Z''")

%Screw axis
quiver3(0,0,0,0.7701,0.5525,0.3190, "r")
text(0.7096,0.5091,0.2939, 'K')

title('Exercise 3')