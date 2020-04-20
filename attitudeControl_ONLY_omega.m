%%MATLAB script developed for attitude control for a dual-sat system in
%%order to create a virtual telescope to study the Alpha Centauri system
close all; clc; clear;
I = [0.09412757, 0.00703819, 0.01209629;
     0.00703819, 0.07495581, 0.03026143;
     0.01209629, 0.03026143, 0.02950416];%moment of inertia matrix
 
%% direction cosine matrix for desired attitude from inertial to body frames 
dcm = [0, -1, 0;
       1,  0, 0;
       0,  0, 1]; % Reference Attitude 
[thetaRx, thetaRy, thetaRz] = dcm2angle(dcm); %calculate needed rotation for desired attitude   
thetaRef = [thetaRx; thetaRy; thetaRz];
qr  = dcm2quat(dcm);

qr_inv = [qr(1) -qr(2:4)];

%gains for controler
Kp = 4; Kd = 1; Ki = 0;

%% Initial Conditions
phi = randn; theta = randn; psi = randn;  % initial conditions for the attitude
attitude = [phi; theta; psi]; 
omega = randn(3,1);
qi = randn(4,1); qi = qi/norm(qi);

%% integration (Dynamics)
options = [];
[t,x] = ode45(@dynamics, [0,6], [qi; omega], options, I,  Kp, Kd, Ki, qr);

%% Convert quaternions to Euler's angles
for i =1:length(t)
    [theta(i,1), theta(i,2), theta(i,3)] = quat2angle(x(i,1:4));
    u(i,:) = - Kp*(x(i,2:4)-qr(2:4))- Kd*x(i,5:end);   
    if (theta(i,1) <= (-1.5551)) && (theta(i,1) >= (-1.5865))
        disp(t(i))
    end
end
maxU = max(abs(u))
%% Calculate the attitude error
% dq =  q (x) qr^-1
for i=1:length(t)
    dq(i,:) = quatmult(x(i,1:4), qr_inv);
%     dq(i,:) = dq(i,:)/norm(dq(i,:));
end
%% Plotting 
figure;
umax = .635;

for i = 1:length(u(:,1))
    if u(i,1) >= umax
        u(i,1) = umax;
    end
    if u(i,1) <= -umax
        u(i,1) = -umax;
    end
    if u(i,2) >= umax
        u(i,2) = umax;
    end
    if u(i,2) <= -umax
        u(i,2) = -umax;
    end
    if u(i,3) >= umax
        u(i,3) = umax;
    end
    if u(i,3) <= -umax
        u(i,3) = -umax;
    end
end


% subplot 412; plot(t,theta(:,1),'+-r', t,theta(:,2),'>-b', t,theta(:,3),'o-m'); hold on;
% xlabel('Time (s)')
% ylabel('Euler Angles [rad]'); grid;
% legend('roll', 'pitch', 'yaw');

subplot 411; plot(t,x(:,5),'+-r', t,x(:,6),'>-b', t,x(:,7),'o-m'); hold on;
xlabel('Time (s)')
ylabel('Angular Velocity \omega [rad/s]'); grid;
legend('\omega_x', '\omega_y', '\omega_z');

subplot 412; plot(t,x(:,1),'+-r', t,x(:,2),'>-b', t,x(:,3),'o-m', t,x(:,4),'.-c'); hold on;
xlabel('Time (s)')
ylabel('Quaternions'); grid;
legend('q_0', 'q_1', 'q_2','q_3');

subplot 413; plot(t,dq(:,1),'+-r', t,dq(:,2),'>-b', t,dq(:,3),'o-m', t,dq(:,4),'.-c'); hold on;
xlabel('Time (s)')
ylabel('Quaternions Error'); grid;
legend('\delta q_0', '\delta q_1', '\delta q_2','\delta q_3');

subplot 414; plot(t,u(:,1),'+-r', t,u(:,2),'>-b', t,u(:,3),'o-m'); hold on;
xlabel('Time (s)')
ylabel('Controll Input [Nm]'); grid;
legend('Roll Control', 'Pitch Control', 'Yaw Control');
ylim();

%% Functions used for controller
function xdot = dynamics(~, x, I, Kp, Kd, Ki, qr)   
    q     = x(1:4); % quaternions 
    omega = x(5:end);   
    B  = [0 -omega';
          omega -tilde(omega)];
    qdot = 1/2*B*q;
    u = - Kd*(omega) - Kp*(q(2:4)-qr(2:4)');
    OmegaDot = I^-1*(cross(-omega,I*omega) + u);
    xdot = [qdot;OmegaDot];
end

function T= tilde(a)
 T = [0 -a(3) a(2);
      a(3) 0 -a(1);
      -a(2) a(1) 0];
end

function dq = quatmult(a,b)
dq = [a(1)*b(1)-a(2:4)*b(2:4)' a(1)*b(2:4)+b(1)*a(2:4)+cross(a(2:4),b(2:4))];
end

