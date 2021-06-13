%%Project matlab file
%created by Kevin Pogrund

%% Create symbolic variables and set up constants
clearvars
clear ArmDynamics
clc

%Setting up variables
syms th1 th2 th3 'real'
syms dth1 dth2 dth3 'real'
syms ddth1 ddth2 ddth3 'real'
syms tau1 tau2 tau3 'real'
% syms m L I 'real'

m = 1;
L = 0.5;
I = 0.1;

q = [th1; th2; th3];
dq = [dth1; dth2; dth3];
ddq = [ddth1; ddth2; ddth3];

%% Positions, velocities and energies
%Positions
p1 = [L*cos(th1)/2; L*sin(th1)/2];
p2 = p1 + [L/2*cos(th1); L/2*sin(th1)] + [L/2*cos(th1 + th2); L/2*sin(th1 + th2)];
p3 = p2 + [L*cos(th1+th2)/2; L*cos(th1+th2)/2] + [L*cos(th1+th2+th3)/2; L*cos(th1+th2+th3)/2];

%Velocities
v1 = simplify(jacobian(p1,q)*dq);
v2 = simplify(jacobian(p2,q)*dq);
v3 = simplify(jacobian(p3,q)*dq);

%Energies
%Body
T1 = 0.5*m*transpose(v1)*v1 + 0.5*I*(dth1)^2;
T2 = 0.5*m*transpose(v2)*v2 + 0.5*I*(dth1 + dth2)^2;
T3 = 0.5*m*transpose(v3)*v3 + 0.5*I*(dth1 + dth2 + dth3)^2;

Ttot = simplify(T1+T2+T3);

%There is no potential energy seeing as we are in space, therefore Vtot=0

%% Variables in the manipulator equation
%Mass matrix
M = simplify(hessian(Ttot, dq));

% Derivative of Mass Matrix
dM = sym(zeros(length(M),length(M)));
for i=1:length(M)
    for j=1:length(M)
        dM(i,j) = jacobian(M(i,j),q)*dq;
    end
end
dM = simplify(dM);

% C Matrix
% Contains the centrifugal and coriolis accelerations
C = dM*dq - transpose(jacobian(Ttot,q));
C = simplify(C);

% G Matrix
% We are in space so G = 0;

%% Run simulink to get q, dq and ddq values
%{
%remove one of the '%' signs to comment this section
disp('simming')
sim('armDynamicsV2.slx');
out = ans; %#ok<NOANS>
clear ans
disp('done simming')

% %theta outputs
th1Out = out.th(:,1);
th2Out = out.th(:,2);
th3Out = out.th(:,3);

%dtheta outputs
dth1Out = out.dth(:,1);
dth2Out = out.dth(:,2);
dth3Out = out.dth(:,3);

%ddtheta outputs
ddth1Out = out.ddth(:,1);
ddth2Out = out.ddth(:,2);
ddth3Out = out.ddth(:,3);
tau1 = out.torque(1);
tau2 = out.torque(2);
tau3 = out.torque(3);
%}
%% Solve for coefficient of friction
% the section above needs to be uncommented
%{
syms b
B = eye(3);
fric = b*eye(3);
tau = [tau1; tau2; tau3];

manipulator = simplify(M*ddq + C + -B*tau + fric*dq);

d = [0 0 0];
for i = 2:size(th3Out, 1)
    EqMan=subs(manipulator(3,:),{th1, th2, th3, dth1, dth2, dth3, ddth1, ddth2,ddth3},{th1Out(i),th2Out(i),th3Out(i),dth1Out(i),dth2Out(i),dth3Out(i),ddth1Out(i),ddth2Out(i),ddth3Out(i)});
    d(3) = d(3)+solve(EqMan,b);
end
d(3)=d(3)/size(th1Out,1);
d
smallB = d(3)

%{
e = 0;
for i = 2:size(th1Out, 1)
    %code
    EqMan=subs(manipulator(3,:),{th1, th2, th3, dth1, dth2, dth3, ddth1, ddth2,ddth3},{th1Out(i),th2Out(i),th3Out(i),dth1Out(i),dth2Out(i),dth3Out(i),ddth1Out(i),ddth2Out(i),ddth3Out(i)});
    e = e+solve(EqMan,b);
end
e = e/size(th1Out,1);
%}
%}
%% stuff
b = 5.56;
B = eye(3);
tau = [tau1;tau2;tau3];
Q = b*dq;
disp('here');
acc_eqns = simplify(M \ (-C + B*tau - Q))
