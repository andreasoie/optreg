DATA_DIR = "data/";
PLOT_DIR = "plots/";

% X_FILENAME = "x16.mat";
% U_FILENAME = "u16.mat";
X_FILE = XF + ".mat";
U_FILE = UF + ".mat";

x_mat = load(DATA_DIR + X_FILE);
u_mat = load(DATA_DIR + U_FILE);

% OPTIMAL
x_data_opt = X_opt.data;

% Extract from Simulink
x_data = x_mat.data;
u_data = u_mat.data;
u_time = u_data(1, :);
u_pitch = u_data(2, :);
u_elevation = u_data(3, :);
u_pitch_opt = Up_opt.data;
u_elevation_opt = Ue_opt.data;
x_time = x_data(1, :);
x_travel = x_data(2, :);
x_travel_rate = x_data(3, :);
x_pitch = x_data(4, :);
x_pitch_rate = x_data(5, :);
x_elevation = x_data(6, :);
x_elevation_rate = x_data(7, :);
x_travel_opt = x_data_opt(:, 1);
x_travel_rate_opt = x_data_opt(:, 2);
x_pitch_opt = x_data_opt(:, 3);
x_pitch_rate_opt = x_data_opt(:, 4);
x_elevation_opt = x_data_opt(:, 5);
x_elevation_rate_opt = x_data_opt(:, 6);

figure('visible','off');
% set(gcf,'position', [0,0, 400, 800]);

subplot(5, 1, 1);
hold on
grid on
plot(x_time, u_pitch_opt, "linewidth", 2, "linestyle","--");
plot(u_time, u_pitch, "linewidth", 2);
xlabel("Time $t$ (s)","Interpreter", "latex")
ylabel("$rad$" , "Interpreter", "latex")
legend('u*')
hold off
title("Input: Pitch")

subplot(5, 1, 2);
hold on
grid on
plot(x_time, u_elevation_opt, "linewidth", 2,"linestyle","--");
plot(x_time, u_elevation, "linewidth", 2);
xlabel("Time $t$ (s)","Interpreter", "latex")
ylabel("$rad$" , "Interpreter", "latex")
legend('$p^{*}$', 'p', "Interpreter", "latex")
hold off
title("Input: Elevation")

subplot(5, 1, 3);
hold on
grid on
plot(x_time, x_travel_opt, "linewidth", 2,"linestyle","--");
plot(x_time, x_travel, "linewidth", 2);
xlabel("Time $t$ (s)","Interpreter", "latex")
ylabel("$rad$" , "Interpreter", "latex")
legend('\lambda*', '\lambda') 
hold off
title("Travel")
% 
% subplot(8, 1, 4);
% hold on
% grid on
% plot(x_time, x_travel_rate_opt, "linewidth", 2,"linestyle","--");
% plot(x_time, x_travel_rate, "linewidth", 2);
% xlabel("Time $t$ (s)","Interpreter", "latex")
% ylabel("$rad$" , "Interpreter", "latex")
% legend('$\dot{\lambda^{*}}$', '$\dot{\lambda}$',"Interpreter", "latex")
% hold off
% title("Travel rate")

subplot(5, 1, 4);
hold on
grid on
plot(x_time, x_pitch_opt, "linewidth", 2,"linestyle","--");
plot(x_time, x_pitch, "linewidth", 2);
xlabel("Time $t$ (s)","Interpreter", "latex")
ylabel("$rad$" , "Interpreter", "latex")
legend('$p*$', '$p$', "Interpreter", "latex") 
hold off
title("Pitch")

% subplot(8, 1, 6);
% hold on
% grid on
% plot(x_time, x_pitch_rate_opt, "linewidth", 2,"linestyle","--");
% plot(x_time, x_pitch_rate, "linewidth", 2);
% xlabel("Time $t$ (s)","Interpreter", "latex")
% ylabel("$rad$" , "Interpreter", "latex")
% legend('$\dot{p*}$', '$\dot{p}$', "Interpreter", "latex") 
% hold off
% title("Pitch rate")

subplot(5, 1, 5);
hold on
grid on
plot(x_time, x_elevation_opt, "linewidth", 2,"linestyle","--");
plot(x_time, x_elevation, "linewidth", 2);
xlabel("Time $t$ (s)","Interpreter", "latex")
ylabel("$rad$" , "Interpreter", "latex")
legend('$e*$', '$e$',"Interpreter", "latex")
hold off
title("Elevation")

% subplot(8, 1, 8);
% hold on
% grid on
% plot(x_time, x_elevation_rate_opt, "linewidth", 2,"linestyle","--");
% plot(x_time, x_elevation_rate, "linewidth", 2);
% xlabel("Time $t$ (s)","Interpreter", "latex")
% ylabel("$rad$" , "Interpreter", "latex")
% legend('$\dot{e^{*}}$', '$\dot{e}$',"Interpreter", "latex")
% hold off
% title("Elevation rate")

% PLOT_TITLE = deblank(XF) + "  " + deblank(UF);
% sgtitle(PLOT_TITLE)

SAVE_NAME = PLOT_DIR + XF + "_" + UF + "_custom.eps";
saveas(gcf, SAVE_NAME)  
