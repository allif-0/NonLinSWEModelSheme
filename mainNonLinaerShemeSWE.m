function mainNonLinaerShemeSWE
% =========================================================================
% Програмний модуль інформаційно-вимірювальної технології визначення джерела 
% забруднюючої речовини з використанням неявної схеми 
% у системі рівнянь мілкої води
% =========================================================================
clear; clc; close all;

%% Параметри системи
Nx = 100; Ny = 100;  % Зменшено для швидкості
Lx = 100; Ly = 100;
dx = Lx/Nx; dy = Ly/Ny;  dt = 0.2;
g = 9.81;  T_end = 25.0;  

max_newt_iters = 10;  tolerance = 1e-5;

%% Вектор невідомих: 3 змінні (H, u, v)
total_dof = 3 * Nx * Ny;
x_old = zeros(total_dof, 1);
x0_point = 27.0; y0_point = 23.0;

% Початкові умови
for i = 1:Nx
    for j = 1:Ny
        I = (i-1)*Ny + j;  % глобальний індекс вузла
        x_coord = (i-1) * dx;
        y_coord = (j-1) * dy;
        
        drop = ((x_coord - x0_point)^2 + ...
                (y_coord - y0_point)^2) / 100.0;
        H_init = 2.0 - min(drop, 1.0);
        
        % Вектор невідомих: [H, u, v] для кожного вузла
        x_old(3*(I-1) + 1) = H_init;  % H
        x_old(3*(I-1) + 2) = 0.0;     % u
        x_old(3*(I-1) + 3) = 0.0;     % v
    end
end

x_current = x_old;

%% Допоміжна функція для отримання індексів вузлів
% Повертає глобальний індекс для вузла (i, j)
get_I = @(i, j) (i-1)*Ny + j;

%% Основний часовий цикл
t = 0.0;
frame_idx = 1;
frames = {};
frames_time = [];

% Налаштування для візуалізації
figure('Position', [100, 100, 900, 700]);

while t < T_end
    fprintf('Розрахунок часу t = %.2f\n', t);
    % Ітераційний цикл Ньютона
    for n = 1:max_newt_iters
        % Ініціалізація розрідженої матриці та вектора правої частини
        A = sparse(total_dof, total_dof);
        b = zeros(total_dof, 1);
        
        for i = 1:Nx
            for j = 1:Ny
                I = get_I(i, j);
                
                row_H = 3*(I-1) + 1;  % індекс для H
                row_u = 3*(I-1) + 2;  % індекс для u
                row_v = 3*(I-1) + 3;  % індекс для v
                
                % ---------------------------------------------------------
                % Граничні умови (непротікання)
                % ---------------------------------------------------------
                if i == 1  % Ліва грань
                    A(row_u, row_u) = 1.0;
                    b(row_u) = 0.0;
                    
                    A(row_H, row_H) = 1.0;
                    I_right = get_I(2, j);
                    A(row_H, 3*(I_right-1) + 1) = -1.0;
                    b(row_H) = 0.0;
                    
                    A(row_v, row_v) = 1.0;
                    A(row_v, 3*(I_right-1) + 3) = -1.0;
                    b(row_v) = 0.0;
                    continue;
                    
                elseif i == Nx  % Права грань
                    A(row_u, row_u) = 1.0; b(row_u) = 0.0;
                    
                    A(row_H, row_H) = 1.0;
                    I_left = get_I(Nx-1, j);
                    A(row_H, 3*(I_left-1) + 1) = -1.0;
                    b(row_H) = 0.0;
                    
                    A(row_v, row_v) = 1.0;
                    A(row_v, 3*(I_left-1) + 3) = -1.0;
                    b(row_v) = 0.0;
                    continue;
                    
                elseif j == 1  % Нижня грань
                    A(row_v, row_v) = 1.0; b(row_v) = 0.0;
                    
                    A(row_H, row_H) = 1.0;
                    I_top = get_I(i, 2);
                    A(row_H, 3*(I_top-1) + 1) = -1.0;
                    b(row_H) = 0.0;
                    
                    A(row_u, row_u) = 1.0;
                    A(row_u, 3*(I_top-1) + 2) = -1.0;
                    b(row_u) = 0.0;
                    continue;
                    
                elseif j == Ny  % Верхня грань
                    A(row_v, row_v) = 1.0; b(row_v) = 0.0;
                    
                    A(row_H, row_H) = 1.0;
                    I_bottom = get_I(i, Ny-1);
                    A(row_H, 3*(I_bottom-1) + 1) = -1.0;
                    b(row_H) = 0.0;
                    
                    A(row_u, row_u) = 1.0;
                    A(row_u, 3*(I_bottom-1) + 2) = -1.0;
                    b(row_u) = 0.0;
                    continue;
                end
                
                % ---------------------------------------------------------
                % Визначення глобальних індексів для сусідніх вузлів
                % ---------------------------------------------------------
                % Індекс початку блоку для кожного вузла (перша змінна - H)
                col_C = 3*(I-1) + 1;                          % центр
                col_L = 3*(get_I(i-1, j)-1) + 1;             % лівий
                col_R = 3*(get_I(i+1, j)-1) + 1;             % правий
                col_B = 3*(get_I(i, j-1)-1) + 1;             % нижній
                col_T = 3*(get_I(i, j+1)-1) + 1;             % верхній
                
                % ---------------------------------------------------------
                % Рівняння 1: нерозривності
                % dH/dt + d(Hu)/dx + d(Hv)/dy = 0
                % ---------------------------------------------------------
                % Коефіцієнт при H_{i,j}^{n+1}
                A(row_H, col_C) = 1.0 / dt;
                
                coef_x = 1.0 / (4.0 * dx);
                
                % Потоки по X
                A(row_H, col_R) = A(row_H, col_R) + coef_x * ...
                                  x_current(col_R + 1 - 1 + 1);
                A(row_H, col_R + 1) = A(row_H, col_R + 1) + coef_x * ...
                                      x_current(col_R);
                
                A(row_H, col_L) = A(row_H, col_L) - coef_x * ...
                                  x_current(col_L + 1 - 1 + 1);
                A(row_H, col_L + 1) = A(row_H, col_L + 1) - coef_x * ...
                                      x_current(col_L);
                
                coef_y = 1.0 / (4.0 * dy);
                
                % Потоки по Y
                A(row_H, col_T) = A(row_H, col_T) + coef_y * ...
                                  x_current(col_T + 2 - 1 + 1);
                A(row_H, col_T + 2) = A(row_H, col_T + 2) + ...
                                      coef_y * x_current(col_T);
                
                A(row_H, col_B) = A(row_H, col_B) - coef_y * ...
                                  x_current(col_B + 2 - 1 + 1);
                A(row_H, col_B + 2) = A(row_H, col_B + 2) - ...
                                      coef_y * x_current(col_B);
                
                % Права частина
                term_dt = x_old(col_C) / dt;
                
                % Потоки по X (старий та поточний шари)
                flux_x_old_R = x_old(col_R) * x_old(col_R + 1);
                flux_x_old_L = x_old(col_L) * x_old(col_L + 1);
                flux_x_cur_R = x_current(col_R) * x_current(col_R + 1);
                flux_x_cur_L = x_current(col_L) * x_current(col_L + 1);
                
                part_x = (flux_x_old_R - flux_x_old_L - ...
                          flux_x_cur_R + flux_x_cur_L) / dx;
                
                % Потоки по Y
                flux_y_old_T = x_old(col_T) * x_old(col_T + 2);
                flux_y_old_B = x_old(col_B) * x_old(col_B + 2);
                flux_y_cur_T = x_current(col_T) * x_current(col_T + 2);
                flux_y_cur_B = x_current(col_B) * x_current(col_B + 2);
                
                part_y = (flux_y_old_T - flux_y_old_B - ...
                          flux_y_cur_T + flux_y_cur_B) / dy;
                
                b(row_H) = term_dt - 0.25 * (part_x + part_y);
                
                % ---------------------------------------------------------
                % Рівняння 2: імпульсу по X
                % d(Hu)/dt + d(Hu^2 + gH^2/2)/dx + d(Huv)/dy = 0
                % ---------------------------------------------------------
                A(row_u, col_C) = x_current(col_C + 1) / dt;
                A(row_u, col_C + 1) = x_current(col_C) / dt;
                
                coef_x = 1.0 / (4.0 * dx);
                
                A(row_u, col_R) = A(row_u, col_R) + coef_x * ...
                                   (x_current(col_R + 1)^2 + ...
                                    g * x_current(col_R));
                A(row_u, col_R + 1) = A(row_u, col_R + 1) + coef_x * ...
                                      (2.0 * x_current(col_R) * x_current(col_R + 1));
                
                A(row_u, col_L) = A(row_u, col_L) - coef_x * ...
                                   (x_current(col_L + 1)^2 + g * x_current(col_L));
                A(row_u, col_L + 1) = A(row_u, col_L + 1) - coef_x * ...
                                       (2.0 * x_current(col_L) * x_current(col_L + 1));
                
                coef_y = 1.0 / (4.0 * dy);
                
                A(row_u, col_T) = A(row_u, col_T) + coef_y * ...
                                   (x_current(col_T + 1) * x_current(col_T + 2));
                A(row_u, col_T + 1) = A(row_u, col_T + 1) + coef_y * ...
                                   (x_current(col_T) * x_current(col_T + 2));
                A(row_u, col_T + 2) = A(row_u, col_T + 2) + coef_y * ...
                                   (x_current(col_T) * x_current(col_T + 1));
                
                A(row_u, col_B) = A(row_u, col_B) - coef_y * ...
                                   (x_current(col_B + 1) * x_current(col_B + 2));
                A(row_u, col_B + 1) = A(row_u, col_B + 1) - coef_y * ...
                                   (x_current(col_B) * x_current(col_B + 2));
                A(row_u, col_B + 2) = A(row_u, col_B + 2) - coef_y * ...
                                   (x_current(col_B) * x_current(col_B + 1));
                
                % Права частина для рівняння імпульсу по X
                term_dt_u = (x_old(col_C) * x_old(col_C + 1)) / dt + ...
                            (x_current(col_C) * x_current(col_C + 1)) / dt;
                
                flux_x_old_R = x_old(col_R) * x_old(col_R + 1)^2 + ...
                               0.5 * g * x_old(col_R)^2;
                flux_x_old_L = x_old(col_L) * x_old(col_L + 1)^2 + ...
                               0.5 * g * x_old(col_L)^2;
                flux_x_cur_R = 2.0 * x_current(col_R) * x_current(col_R + 1)^2 + ...
                               0.5 * g * x_current(col_R)^2;
                flux_x_cur_L = 2.0 * x_current(col_L) * x_current(col_L + 1)^2 + ...
                               0.5 * g * x_current(col_L)^2;
                
                part_x_u = (flux_x_old_R - flux_x_old_L - ...
                            flux_x_cur_R + flux_x_cur_L) / dx;
                
                flux_y_old_T = x_old(col_T) * x_old(col_T + 1) * x_old(col_T + 2);
                flux_y_old_B = x_old(col_B) * x_old(col_B + 1) * x_old(col_B + 2);
                flux_y_cur_T = 2.0 * x_current(col_T) * ...
                               x_current(col_T + 1) * x_current(col_T + 2);
                flux_y_cur_B = 2.0 * x_current(col_B) * ...
                               x_current(col_B + 1) * x_current(col_B + 2);
                
                part_y_u = (flux_y_old_T - flux_y_old_B - ...
                            flux_y_cur_T + flux_y_cur_B) / dy;
                
                b(row_u) = term_dt_u - 0.25 * (part_x_u + part_y_u);
                
                % ---------------------------------------------------------
                % Рівняння 3: імпульсу по Y
                % d(Hv)/dt + d(Huv)/dx + d(Hv^2 + gH^2/2)/dy = 0
                % ---------------------------------------------------------
                A(row_v, col_C) = x_current(col_C + 2) / dt;
                A(row_v, col_C + 2) = x_current(col_C) / dt;
                
                coef_x = 1.0 / (4.0 * dx);
                
                A(row_v, col_R) = A(row_v, col_R) + coef_x * ...
                                   (x_current(col_R + 1) * x_current(col_R + 2));
                A(row_v, col_R + 1) = A(row_v, col_R + 1) + coef_x *...
                                       (x_current(col_R) * x_current(col_R + 2));
                A(row_v, col_R + 2) = A(row_v, col_R + 2) + coef_x * ...
                                       (x_current(col_R) * x_current(col_R + 1));
                
                A(row_v, col_L) = A(row_v, col_L) - coef_x * ...
                                   (x_current(col_L + 1) * x_current(col_L + 2));
                A(row_v, col_L + 1) = A(row_v, col_L + 1) - coef_x * ...
                                       (x_current(col_L) * x_current(col_L + 2));
                A(row_v, col_L + 2) = A(row_v, col_L + 2) - coef_x * ...
                                       (x_current(col_L) * x_current(col_L + 1));
                
                coef_y = 1.0 / (4.0 * dy);
                
                A(row_v, col_T) = A(row_v, col_T) + coef_y * ...
                                   (x_current(col_T + 2)^2 + g * x_current(col_T));
                A(row_v, col_T + 2) = A(row_v, col_T + 2) + coef_y * ...
                                       (2.0 * x_current(col_T) * x_current(col_T + 2));
                
                A(row_v, col_B) = A(row_v, col_B) - coef_y * ...
                                   (x_current(col_B + 2)^2 + g * x_current(col_B));
                A(row_v, col_B + 2) = A(row_v, col_B + 2) - coef_y * ...
                                   (2.0 * x_current(col_B) * x_current(col_B + 2));
                
                % Права частина для рівняння імпульсу по Y
                term_dt_v = (x_old(col_C) * x_old(col_C + 2)) / dt + ...
                            (x_current(col_C) * x_current(col_C + 2)) / dt;
                
                flux_x_old_R = x_old(col_R) * x_old(col_R + 1) * x_old(col_R + 2);
                flux_x_old_L = x_old(col_L) * x_old(col_L + 1) * x_old(col_L + 2);
                flux_x_cur_R = 2.0 * x_current(col_R) * ...
                               x_current(col_R + 1) * x_current(col_R + 2);
                flux_x_cur_L = 2.0 * x_current(col_L) * ...
                               x_current(col_L + 1) * x_current(col_L + 2);
                
                part_x_v = (flux_x_old_R - flux_x_old_L - ...
                            flux_x_cur_R + flux_x_cur_L) / dx;
                
                flux_y_old_T = x_old(col_T) * x_old(col_T + 2)^2 + ...
                               0.5 * g * x_old(col_T)^2;
                flux_y_old_B = x_old(col_B) * x_old(col_B + 2)^2 + ...
                               0.5 * g * x_old(col_B)^2;
                flux_y_cur_T = 2.0 * x_current(col_T) * x_current(col_T + 2)^2 + ...
                               0.5 * g * x_current(col_T)^2;
                flux_y_cur_B = 2.0 * x_current(col_B) * x_current(col_B + 2)^2 + ...
                               0.5 * g * x_current(col_B)^2;
                
                part_y_v = (flux_y_old_T - flux_y_old_B - ...
                            flux_y_cur_T + flux_y_cur_B) / dy;
                
                b(row_v) = term_dt_v - 0.25 * (part_x_v + part_y_v);
            end
        end
        
        % Розв'язання СЛАР
        x_next = A \ b;
        
        % Перевірка збіжності Ньютона
        error = norm(x_next - x_current, inf);
        x_current = x_next;
        
        if error < tolerance
            fprintf('Кількість ітерації для збіжності метода Ньютона: %d. Похибка: %.2e\n', n, error);
            break;
        end
        
        if n == max_newt_iters
            fprintf('  Ньютон не зійшовся. Похибка: %.2e\n', error);
        end
    end
    
    % Оновлення часового шару
    x_old = x_current;
    t = t + dt;
    
    % Збереження кадру для візуалізації
    H_frame = zeros(Nx, Ny);
    for i = 1:Nx
        for j = 1:Ny
            I = get_I(i, j);
            H_frame(i, j) = x_old(3*(I-1) + 1);
        end
    end
    frames{frame_idx} = H_frame;
    frames_time(frame_idx) = t;
    frame_idx = frame_idx + 1;
    
    % Візуалізація поточного стану
    if mod(round(t/dt), 5) == 0 || t <= dt || t >= T_end - dt/2
        subplot(1, 2, 1);
        imagesc(H_frame');
        axis equal tight;
        colorbar;
        title(sprintf('H (t = %.2f с)', t));
        xlabel('X'); ylabel('Y');
        set(gca, 'YDir', 'normal');
        caxis([min(H_frame(:)), max(H_frame(:))]);
        
        subplot(1, 2, 2);
        surf(H_frame');
        axis tight;
        colorbar;
        title(sprintf('3D поверхня (t = %.2f с)', t));
        xlabel('X'); ylabel('Y'); zlabel('H');
        view(45, 30);
        zlim([min(H_frame(:))-0.1, max(H_frame(:))+0.1]);
        
        drawnow;
    end
end

%% Збереження даних у MAT-файл
save('shallow_water_results.mat', 'frames', 'frames_time', 'dx', 'dy', 'dt', 'Nx', 'Ny', 'Lx', 'Ly');
fprintf('Дані збережено у файл: shallow_water_results.mat\n');