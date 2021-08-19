clc; clear; close all;
warning('off', 'all');

%-------------------------------------------------------------------------%
% 1) Initialize Data
%-------------------------------------------------------------------------%

load Measuerment_LOS_NLOS

color = lines(20);

data = [];
data_dim = 3;

cluster = cell(1, 10);
sigma_init = diag([1.^2 1.^2 1.^2]);

sigma_y = cell(1, 10);
tau_y = cell(1, 10);


% Hyperparameter
alpha = 0.01;
mc = 200; % 200

%-------------------------------------------------------------------------%
% 2) DP clustering & Plot
%-------------------------------------------------------------------------%

for ti = 1 : mc
    fprintf('********************> DP clustering... (%d / %d)\n', ti, mc)
    
    % - NLOS -> LOS
    %     for wow = 1 : 2
    %         if wow == 1
    %             m = 2;
    %         else
    %             m = 1;
    %         end
    
    % - LOS -> NLOS
    for m = 1 : 2
        
        % ti번째에 들어온 m번째 데이터 (1 : LOS, 2 : NLOS)
        rv = v.Time(ti).measurement(m, [1 4 5]);
    
        % 지금까지 들어온 누적 데이터
        data = [data; rv];
        data_N = size(data,1);
        if (ti == 1)
            % ti = 1일 때, LOS와 NLOS 신호를 각각 cluster 1, 2라고 가정
            z = [1; 2];
            Nclust = max(z);
            n_k = countNum(z, Nclust);
            str_c = ["Cluster1", "Cluster2"];
            new_z = m;
            
            cluster{m} = v.Time(ti).measurement(m, [1 4 5]);
            sigma_y{m} = sigma_init;
            tau_y{m} = inv(sigma_y{m});
        else
            % 1. POSTERIOR 계산
            p = NaN(Nclust+1, 1);
            prior_new = alpha / (data_N - 1 + alpha);

            for k = 1 : Nclust
                prior_exist = n_k(k) / (data_N - 1 + alpha);
                p(k) = pp_exist(prior_exist, n_k, k, data, data_N, tau_0, tau_y, mu_0, sigma_y, cluster);
            end
            p(Nclust+1) = pp_new(prior_new, data, data_N, mu_0, sigma_0, sigma_init);
            % 2. Posterior에 따라 새로운 클러스터에 할당!
            new_z = randsample(Nclust+1, 1, true, p);

            % - 이때 새 클러스터가 추가될 경우
            if new_z == Nclust+1
                n_k = ([n_k 1]);
                Nclust = Nclust + 1;
                str_c = [str_c, "Cluster" + new_z];
                cluster{new_z} = rv;
                sigma_y{new_z} = sigma_init;
                
            % - 기존의 클러스터에 추가될 경우
            else
                n_k(new_z) = n_k(new_z) + 1;
                cluster{new_z} = [cluster{new_z}; rv];
                sigma_y{new_z} = diag(diag([cov(cluster{new_z})]));
            end     
            z = [z; new_z];
        end
        % 3. mu_0, sigma_0, sigma_y 업데이트
        mu_0 = mean(data);
        sigma_0 = diag(std(data).^2);
        tau_0 = inv(sigma_0);
        
        tau_y{new_z} = inv(sigma_y{new_z}); 
        
        disp(n_k)
        
        % Plot
        set(gcf,'units','normalized','outerposition',[0 1/4.5 0.6 1/1.5]);
        subplot(1, 2, 1)
%         scatter(v.Time(ti).measurement(m, 1), v.Time(ti).measurement(m, 4), 15, color(m,:), 'filled')
        scatter3(v.Time(ti).measurement(m, 1), v.Time(ti).measurement(m, 4), v.Time(ti).measurement(m, 5), 15, color(m,:), 'filled')
        hold on
        grid on
        axis([350 600 0 7 -0.03 0.04])
        xlabel('X-axis')
        ylabel('Y-axis')
        zlabel('Z-axis')
        title('Correct Answer')
        legend('LOS', 'NLOS', 'Location', 'northeast')
        subplot(1, 2, 2)
%         scatter(v.Time(ti).measurement(m, 1), v.Time(ti).measurement(m, 4), 15, color(new_z,:), 'filled')
        scatter3(v.Time(ti).measurement(m, 1), v.Time(ti).measurement(m, 4), v.Time(ti).measurement(m, 5), 15, color(new_z,:), 'filled')
        hold on
        grid on
        axis([350 600 0 7 -0.03 0.04])
        xlabel('X-axis')
        ylabel('Y-axis')
        zlabel('Z-axis')
        title('DP result')
        legend(str_c, 'Location', 'northeast')
        drawnow
    end
end

%-------------------------------------------------------------------------%
% 3) Function
%-------------------------------------------------------------------------%

function n_k = countNum(z, Nclust)
n_k = [];
for k = 1 : Nclust
    n_k = [n_k numel(z(z == k))];
end
end

function p_e = pp_exist(prior_exist, n_k, k, data, data_N, tau_0, tau_y, mu_0, sigma_y, cluster)
tau_p = tau_0 + n_k(k) * tau_y{k};
sigma_p = inv(tau_p);

if size(cluster{k}, 1) == 1;
    sum_data = cluster{k};
else
    sum_data = sum(cluster{k});
end

% sum_data = 0;
% for n = 1 : data_N-1
%     if z(n) == k
%         sum_data = sum_data + data(n, :)
%     end
% end


mu_p = transpose(sigma_p * (tau_y{k} * transpose(sum_data) + tau_0 * transpose(mu_0)));
p_e = (prior_exist * mvnpdf(data(data_N, :), mu_p, sigma_p+sigma_y{k}));
end

function p_n = pp_new(prior_new, data, data_N, mu_0, sigma_0, sigma_init)
p_n = (prior_new * mvnpdf(data(data_N, :), mu_0, sigma_0+sigma_init)); % sigma_y{k} 대신 sigma_init
end
