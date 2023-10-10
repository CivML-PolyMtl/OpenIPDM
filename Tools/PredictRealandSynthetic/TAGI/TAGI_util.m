classdef TAGI_util
    methods (Static)
%         function true_init_speed = get_true_speed(file_path)
%             load(file_path);
%             true_state = GeneratedData.TrueState;
% 
%             start_year = 1964; % equal to index 1
%             all_years = cellfun(@(x) x(1,3),GeneratedData.ShortObserved,'UniformOutput',false);
%             years_index = cell2mat(all_years) - start_year;
% 
%             all_true_speeds = cellfun(@(x) x(2,:),true_state,'UniformOutput',false);
%             all_true_speeds = cell2mat(all_true_speeds);
% 
%             true_init_speed = zeros(length(years_index),1);
%             for i = 1:length(years_index)
%                 yr_idx = years_index(i);
%                 true_init_speed(i) = all_true_speeds(i,yr_idx);
%             end
% 
%             train_indx = gather(MdataEngy.StructureIndS);
%             true_init_speed = true_init_speed(train_indx);
%         end

        function save_dir = create_dir(custom_folder_name)
            d = datetime('now', 'Format','yyyy_MM_dd');
            dt = datetime('now', 'Format','yyyy_MM_dd-HH_mm');
            save_dir = sprintf('%s/Tools/PredictRealandSynthetic/Graphs/%s', pwd, custom_folder_name);
            if ~exist(save_dir, 'dir')
                mkdir(save_dir);
            end
        end
        % Saving functions
        function save_model(netInfo, save_dir)
            if ~exist(save_dir, 'dir')
                mkdir(save_dir);
            end
            save(sprintf('%s/net_info.mat', save_dir), 'netInfo');
            fid = fopen(sprintf('%s/net_info.json', save_dir),'w');
            encoded = jsonencode(netInfo);
            fprintf(fid,'%s',encoded);
            fclose(fid);
        end
        function save_plot(fig, save_dir, name)
%             save_dir = sprintf('%s/Figures', save_dir);
            if ~exist(save_dir, 'dir')
                mkdir(save_dir);
            end
            if isnumeric(name)
                name = string(name);
            end
            saveas(fig, sprintf('%s/fig_%s.png', save_dir, name));
            saveas(fig, sprintf('%s/fig_%s.fig', save_dir, name));
        end
        function save_train_data(x, save_dir, train)
            figure('visible','off');
            fig = histogram(x);
            xlabel('Structural attributes');
            ylabel('Frequency');
            save_dir = sprintf('%s/Training_data', save_dir);
            if ~exist(save_dir, 'dir')
                mkdir(save_dir);
            end
            if train
                title('Training');
                saveas(fig, sprintf('%s/x_train.png', save_dir));
            else
                title('Validation');
                saveas(fig, sprintf('%s/x_valid.png', save_dir));
            end
        end
        % Ploting functions
        function plot_pred(x_pred, y_pred, y_pred_var, x_train, y_train, x_val, y_val, RegressLoop, display)
            if display
                figure;
            else
                figure('visible','off');
            end
            scatter(x_train, y_train, 10, 'magenta', 'd', 'DisplayName','training data')
            if ~isempty(x_val) && ~isempty(y_val)
                hold on
                scatter(x_val, y_val, 10, 'blue', 'x', 'DisplayName','validation data')
            end
            hold on
            pl.regression(x_pred, y_pred, y_pred_var, 'black', 'green', 1);
            xlabel('$z$ (structural attribute)', 'Interpreter','latex')
            ylabel('$\dot{\mu}_{z}$', 'Interpreter','latex')  
            if ~isempty(RegressLoop)
                title(strcat('Regress loop: ', string(RegressLoop)))
            else
                title('Final results')
            end
        end
        function plot_pred2(x_pred, y_pred, x_pred_var, x_train, y_train, x_val, y_val, RegressLoop, display)
            % note x will be speed and y will be struc att
            f=1;
            if display
                figure;
            else
                figure('visible','off');
            end
            scatter(x_train, y_train, 10, 'magenta', 'd', 'DisplayName','training data')
            if ~isempty(x_val) && ~isempty(y_val)
                hold on
                scatter(x_val, y_val, 10, 'blue', 'x', 'DisplayName','validation data')
            end
            hold on
            [y_pred, idx] = sort(y_pred);
            xpl = x_pred(idx);
            x_pl_std = sqrt(x_pred_var(idx));
            plot(xpl, y_pred, 'color', 'black', 'LineWidth', 1, 'DisplayName', 'TAGI pred.')
            hold on
            patch([(xpl-f*x_pl_std)' fliplr((xpl+f*x_pl_std)')],[y_pred',  fliplr(y_pred')],...
                'green', 'EdgeColor','none','FaceColor', 'green','FaceAlpha', 0.2, 'DisplayName', strcat('$\dot{\mu}_{x}\pm$', num2str(f), '$\sigma$'))
            xlabel('$\dot{x_0}$', 'Interpreter','latex')
            ylabel('$z$ (structural attribute)', 'Interpreter','latex')
            h = legend('Location','southwest');
            set(h,'Interpreter','latex');
            if ~isempty(RegressLoop)
                title(strcat('Regress loop: ', string(RegressLoop)))
            else
                title('Final results')
            end
        end
        function plot_true(x_true, y_true, hold_on, deterministic)
            if hold_on
                hold on
                plot(x_true, y_true, 'red', 'DisplayName','true', 'LineWidth', 1)
                if ~deterministic
                    hold on
                    var = (exp(0.1^2)-1)*exp(0.1^2);
                    if any((y_true' + sqrt(var)) > 0)
                        up_bound = y_true' + sqrt(var);
                        up_bound(up_bound>0) = 0;
                    end
                    patch([x_true' fliplr(x_true')],[up_bound,  fliplr(y_true' - sqrt(var))], 'red', 'EdgeColor','none','FaceColor', 'red','FaceAlpha', 0.2, 'DisplayName','true $\pm \sigma_{w_0}$')
                end
            else
                figure
                plot(x_true, y_true, 'red', 'DisplayName','true', 'LineWidth', 1)
                if ~deterministic
                    var = (exp(0.1^2)-1)*exp(0.1^2);
                    if any((y_true' + sqrt(var)) > 0)
                        up_bound = y_true' + sqrt(var);
                        up_bound(up_bound>0) = 0;
                    end
                    patch([x_true' fliplr(x_true')],[up_bound,  fliplr(y_true' - sqrt(var))], 'red', 'EdgeColor','none','FaceColor', 'red','FaceAlpha', 0.2, 'DisplayName','true $\pm \sigma_{w_0}$')
                end
                xlabel('$z$ (structural attribute)', 'Interpreter','latex')
                ylabel('$\dot{x}_{0}$', 'Interpreter','latex')
                legend
            end
        end
        function plot_true2(x_true, y_true, std, f)
            hold on
%             [x_true, idx] = sort(x_true);
%             y_true = y_true(idx);
            plot(x_true, y_true, 'red', 'DisplayName','$\mu_{z}$', 'LineWidth', 1)
            hold on
            patch([x_true fliplr(x_true)],[y_true + f*std,  fliplr(y_true - f*std)], 'red', 'EdgeColor','none','FaceColor', 'red','FaceAlpha', 0.2, 'DisplayName',strcat('$\mu_{z}\pm$', num2str(f), '$\sigma_{w_0}$'))
        end
        function scatter_true_speed(pred_init_speed, true_init_speed)
            figure('visible','off');
            scatter(pred_init_speed, true_init_speed);
            xlabel('predicted speed');
            ylabel('true speed');
            axis equal;
            hold on
            plot(linspace(-0.7,0), linspace(-0.7,0), 'LineWidth', 1);
            xlim([-0.7 0]);
            ylim([-0.7 0]);
        end
        function scatter_true(pred_init, true_init)
            figure('visible','off');
            scatter(pred_init, true_init);
            xlabel('predicted');
            ylabel('true');
            min_val = min([pred_init; true_init]);
            max_val = max([pred_init; true_init]);
            axis equal;
            hold on
            plot(linspace(min_val,max_val), linspace(min_val,max_val), 'LineWidth', 1);
            xlim([min_val max_val]);
            ylim([min_val max_val]);
        end
        function hist_true_speed(pred_init_speed, true_init_speed)
            figure('visible','off');
            histogram(true_init_speed - pred_init_speed);
            xlabel('true - predicted');
            ylabel('frequency');
        end
        function x_one_hot = one_hot_encode_inp(x_cat, categories)
            x_one_hot = onehotencode(x_cat,2,'ClassNames',categories);
        end
        function [x_cat, x_tr] = split_cat_and_cont_inp(x_tr, cat_inp_idx)
            x_cat = x_tr(:,1:cat_inp_idx);
            x_tr = x_tr(:,cat_inp_idx+1:end);
        end
    end
end