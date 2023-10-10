%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% File:         pl
% Description:  plot figurue for TAGI
% Authors:      Luong-Ha Nguyen & James-A. Goulet 
% Created:      November 03, 2019
% Updated:      July 29, 2021
% Contact:      luongha.nguyen@gmail.com & james.goulet@polymtl.ca 
% Copyright (c) 2021 Luong-Ha Nguyen & James-A. Goulet. All rights reserved. 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
classdef pl
    methods (Static)
        function diagnosticCARL(myhist, mQhist, baseline, stat, loop, nsteps, epsCount)
            lw       = 2;
            numSteps = cumsum(stat(1:epsCount, 2));
            avgRew   = stat(1:epsCount, 4);
            my       = reshape(myhist(1:loop, :)', [loop*nsteps, 1]);
            mQ       = reshape(mQhist(1:loop, :)', [loop*nsteps, 1]);
            xbaseline= [0:loop]'.*nsteps;
            
            figure(1)
            plot(numSteps/1e6, avgRew/100, 'r', 'linewidth', lw);
            hold on
            plot([1:loop*nsteps]'/(1e6), my, '.k');
            plot([1:loop*nsteps]'/(1e6), mQ, '.g');
            plot(xbaseline/1e6, [baseline(1);baseline(1:loop)], 'm', 'linewidth', lw);
            xlabel('Num. of steps (M)')
            ylabel('rewards') 
            legend('avg. rew/100', 'returns', 'deltaQ-value', 'baseline')
            xlim([0, loop*nsteps/1E6])
        end
        function processedImgRL(img, rows, columns, name, export)
            %Mnist latent
            set(gcf,'Units','centimeters');
            set(gcf, 'Position', [1, 1, 12, 12])
            tiledlayout(rows, columns,'TileSpacing','none','Padding','none')
            numPlots = size(img, 4);
            for i = 1:numPlots
                nexttile
                imshow(img);
            end
            if export == 1
                fig = gcf;
                pp  = fig.PaperPosition;
                pp(1:2) = -0.15;
                ps  = pp(3:4);
                set(gcf,'Units','centimeters');
                set( gcf,'PaperSize', ps, 'PaperPosition', pp)
                print(name, '-dpdf', '-r0')
            end
        end
        function imgRL(img, rows, columns, name, export)
            %Mnist latent
            set(gcf,'Units','centimeters');
            set(gcf, 'Position', [1, 1, 12, 15.4])
            tiledlayout(rows, columns,'TileSpacing','none','Padding','none')
            numPlots = size(img, 1);
            for i = 1:numPlots
                nexttile
                imshow(img{i}./255);
            end
            if export == 1
                fig = gcf;
                pp  = fig.PaperPosition;
                pp(1:2) = -0.15;
                ps  = pp(3:4);
                set(gcf,'Units','centimeters');
                set( gcf,'PaperSize', ps, 'PaperPosition', pp)
                print(name, '-dpdf', '-r0')
            end
        end
        function regression(x, y, Sy, c1, c2, f)
            [xpl, idx] = sort(x);
            plot(xpl, y(idx), 'color', c1, 'LineWidth', 1, 'DisplayName', '$\dot{\mu}_{z}$')
            hold on
            patch([xpl' fliplr(xpl')],[y(idx)' + f*sqrt(Sy(idx)'),  fliplr(y(idx)' - f*sqrt(Sy(idx)'))],...
                c2, 'EdgeColor','none','FaceColor', c2,'FaceAlpha', 0.2, 'DisplayName', strcat('$\dot{\mu}_{z}\pm$', num2str(f), '$\sigma$'))
            % xlabel('x', 'Interpreter','latex')
            % ylabel('y', 'Interpreter','latex')            
%             h = legend('$\dot{\mu}_{z}$', ['$\dot{\mu}_{z}\pm$', num2str(f), '$\sigma$']);
            h = legend('Location','southwest');
            set(h,'Interpreter','latex');
        end
        function inversedNet(y, x, Sx, c1, c2, f)
            [xpl, idx] = sort(x);
            plot(xpl, y(idx), '-+', 'color', c1)
            hold on
            patch([x(idx)' + f*sqrt(Sx(idx)') fliplr(x(idx)' - f*sqrt(Sx(idx)'))], [y(idx)' fliplr(y(idx)')], c2, 'EdgeColor','none','FaceColor', c2,'FaceAlpha', 0.2)
            xlabel('$x$', 'Interpreter','latex')
            ylabel('$y$', 'Interpreter','latex')
%             h = legend('$E[y]$', ['$E[y]\pm3$', num2str(f), '$\sigma$']);
%             set(h,'Interpreter','latex');
        end
        function multipleImges(img)
            numPlots = size(img, 4);
            ndiv     = sqrt(numPlots);
            idx      = 1:numPlots;%randperm(size(img, 4), numPlots);
            img      = img(:, :, :, idx);
            for i = 1:numPlots
                ax = subplot(ndiv,ndiv,i);
                imshow(img(:,:,:,i), 'Parent', ax);
            end
        end
        function class(img, rows, columns, name, export)  
%             %Mnist latent
%             set(gcf,'Units','centimeters');
%             set(gcf, 'Position', [1, 1, 2.3, 2.3])
            % CelebA plot
            set(gcf,'Units','centimeters');
            set(gcf, 'Position', [1, 1, 7.92, 7.92])
            %CelebA latent vatiables
%             set(gcf,'Units','centimeters');
%             set(gcf, 'Position', [1, 1, 7.9, 7.9/10*5])
            tiledlayout(rows, columns,'TileSpacing','none','Padding','none')
            numPlots = size(img, 4);
            idx      = 1:numPlots;
            img      = img(:, :, :, idx);
            for i = 1:numPlots
                nexttile
                imshow(img(:,:,:,i));
            end
            if export == 1
%                 % Mnist export
%                 set(gcf,'Units','centimeters');
%                 set( gcf,'PaperSize',[2.5 2.5], 'PaperPosition',[0 0 2.5 2.5])
%                 print('test', '-dpdf')
                % celebA export
                fig = gcf;
                pp  = fig.PaperPosition;
                pp(1:2) = 0;
                ps  = pp(3:4);
                set(gcf,'Units','centimeters');
                set( gcf,'PaperSize', ps, 'PaperPosition', pp)
                print(name, '-dpdf')
            end
        end
        function latenVariables(img)
            numPlots = size(img, 4);
            idx      = 1:numPlots;%randperm(size(img, 4), numPlots);
            img      = img(:, :, :, idx);
            for i = 1:numPlots
                ax = subplot(5,10,i);
                imshow(img(:,:,:,i), 'Parent', ax);
            end
        end
        function vflatenVariables(img, imgSize, nlv, ns)
            img = dp.cvr2img(img, imgSize);
            numPlots = size(img, 4);
            idx      = 1 : numPlots;
            img      = img(:, :, :, idx);
            tiledlayout(nlv, ns,'TileSpacing','none','Padding','none')
            for i = 1 : numPlots
                nexttile
                imshow(img(:,:,:,i));
            end
        end
        function imgComparison(imgP, imgF, imgR)
            S        = size(imgP);
            img      = zeros(S(1), S(2), S(3), S(4)*3, 'like', imgP);
            img(:,:,:,1:S(4)) = repmat(imgR, [1,1,1,S(4)]);
            img(:,:,:,S(4)+1:2*S(4)) = imgP;
            img(:,:,:,2*S(4)+1:end) = imgF;
            numPlots = size(img, 4);
            idx      = 1:numPlots;%randperm(size(img, 4), numPlots);
            img      = img(:, :, :, idx);
            for i = 1:numPlots
                ax = subplot(3,numPlots/3,i);
                imshow(img(:,:,:,i), 'Parent', ax);
            end
        end
        function imgComparison_V2(imgP, imgR)
            S        = size(imgP);
            img      = zeros(S(1), S(2), S(3), S(4)*2, 'like', imgP);
            img(:,:,:,1:S(4)) = imgR;
            img(:,:,:,S(4)+1:end) = imgP;
            numPlots = size(img, 4);
            idx      = 1:numPlots;%randperm(size(img, 4), numPlots);
            img      = img(:, :, :, idx);
            for i = 1:numPlots
                ax = subplot(2,numPlots/2,i);
                imshow(img(:,:,:,i), 'Parent', ax);
            end
        end
        function trainImage(img, imgSize)
            img = dp.cvr2img(img, imgSize);
            for i = 1:size(img, 4)
                ax = subplot(1,size(img, 4),i);
                imshow(img(:,:,:, i), 'Parent', ax)
            end
        end
        function plotClassProb_V2(prob_class, y_obs)
            export_plot=1;
            y_test=y_obs;      
            pr=sortrows([prob_class(:,1:10),y_test],11);
            idx_y=0;
            for i=0:9
                idx=find(pr(:,11)==i);
                pr(idx,:)=sortrows(pr(idx,:),-(i+1));
                idx_y=[idx_y,idx_y(end)+numel(idx)];
            end
            
            p_class=[0.1:0.01:0.99 0.995 0.997 0.998 0.999];
            class_summary=zeros(10000,length(p_class));
            for i=1:10000
                for j=1:length(p_class)
                    if and(pr(i,pr(i,end)+1)==max(pr(i,1:10)),pr(i,pr(i,end)+1)>=p_class(j))
                        class_summary(i,j)=1;
                    elseif all(pr(i,1:10)<p_class(j))
                        class_summary(i,j)=2;
                    elseif and(pr(i,pr(i,end)+1)~=max(pr(i,1:10)),max(pr(i,1:10))>=p_class(j))
                        class_summary(i,j)=3;
                    end
                end
            end
            stat=[mean(class_summary==1);mean(class_summary==2);mean(class_summary==3)];
            
            
            figure('Position', [100 100 200 100]);
            h=area(p_class, stat');
            h(1).FaceColor = 'green';
            h(2).FaceColor = [1,1,0.1];
            h(3).FaceColor = 'red';
            xlim([0.1,0.999]);
            ylim([0,0.999]);
            ax = gca;
            ax.XTick = [0.1 0.5 0.9];
            %ax.YTick = [0 1];
            
            %xticklabels({})
            %xlabel('Threshold Pr')
            ylabel('$\Pr(~)$','Interpreter','Latex')
            set(gcf,'Color',[1 1 1])
            opts=['scaled y ticks = false,',...
                'scaled x ticks = false,',...
                'x label style={font=\large},',...
                'y label style={font=\large},',...
                'z label style={font=\large},',...
                'mark size=5',...
                ]; %'legend style={font=\large},',...
            %'title style={font=\large},',...
            if export_plot==1
                matlab2tikz('figurehandle',gcf,'filename',[ 'saved_figures/' 'mnist_' '.tex'] ,'standalone', true,'showInfo', false,'floatFormat','%.5g','extraTikzpictureOptions','font=\large','extraaxisoptions',opts);
            end
        end
        function plotClassProb(prob_class, y_obs)
            export_plot=0;
            y_test=y_obs;
            %correct_class=prob_class(metric.errorRate==0,1:10);
            %wrong_class=prob_class(metric.errorRate==1,1:10);
            
            pr=sortrows([prob_class(:,1:10),y_test],11);
            idx_y=0;
            for i=0:9
                idx=find(pr(:,11)==i);
                pr(idx,:)=sortrows(pr(idx,:),-(i+1));
                idx_y=[idx_y,idx_y(end)+numel(idx)];
            end
            
            p_class=[0.1:0.01:0.99 0.995 0.997 0.998 0.999];
            class_summary=zeros(10000,length(p_class));
            for i=1:10000
                for j=1:length(p_class)
                    if and(pr(i,pr(i,end)+1)==max(pr(i,1:10)),pr(i,pr(i,end)+1)>=p_class(j))
                        class_summary(i,j)=1;
                    elseif all(pr(i,1:10)<p_class(j))
                        class_summary(i,j)=2;
                    elseif and(pr(i,pr(i,end)+1)~=max(pr(i,1:10)),max(pr(i,1:10))>=p_class(j))
                        class_summary(i,j)=3;
                    end
                end
            end
            stat=[mean(class_summary==1);mean(class_summary==2);mean(class_summary==3)];
            
            
            figure('Position', [0 0 450 100]);
            subplot(1,3,1:2)
            h=imagesc(pr(:,1:10)',[0,1])
            
            colormap(jet)
            %colorbar
            ax = gca;
            ax.XTick = idx_y;
            ax.YTick = [-0.5:1:10.5];
            
%             ax.XColor = [1 1 1];
%             ax.YColor = [1 1 1];
            
            
            ax.GridColor = [1 1 1];
            ax.GridAlpha = 1;
            xticklabels({0 1 2 3 4 5 6 7 8 9})
            yticklabels({0 1 2 3 4 5 6 7 8 9})
            ylabel('True labels (0-9)')
            xlabel('Test set labels (0-9)')
            grid on
            
            
            
            subplot(1,3,3)
            h=area(p_class,[stat]')
            h(1).FaceColor = 'green';
            h(2).FaceColor = [1,1,0.1];
            h(3).FaceColor = 'red';
            xlim([0.1,0.999]);
            ylim([0,0.999]);
            ax = gca;
            ax.XTick = [0.1 0.5 0.9];
            %ax.YTick = [0 1];
            
            %xticklabels({})
            %xlabel('Threshold Pr')
            ylabel('$\Pr(~)$','Interpreter','Latex')
            set(gcf,'Color',[1 1 1])
            opts=['scaled y ticks = false,',...
                'scaled x ticks = false,',...
                'x label style={font=\large},',...
                'y label style={font=\large},',...
                'z label style={font=\large},',...
                'mark size=5',...
                ]; %'legend style={font=\large},',...
            %'title style={font=\large},',...
            if export_plot==1
                matlab2tikz('figurehandle',gcf,'filename',[ 'mnist_' num2str(nb_train) '.tex'] ,'standalone', true,'showInfo', false,'floatFormat','%.5g','extraTikzpictureOptions','font=\large','extraaxisoptions',opts);
            end
        end
    end
end