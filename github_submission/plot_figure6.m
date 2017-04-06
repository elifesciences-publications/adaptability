% curr_dir = '/Users/skryazhi/Dropbox WORK/Dropbox/kruglyak_assays_and_qtl_detection/_SKANALYSIS/data';
curr_dir = '/Users/skryazhi/Dropbox WORK/Dropbox/kruglyak_adaptability_paper/data 2017-01-10';
% filename = [curr_dir '/krug_data_01.mat'];
filename = [curr_dir '/krug_data_02.mat'];
load(filename, 'krug_data');

IFLOADSTAT = true;
% stat_filename = [curr_dir 'krug_data_01_stat.mat'];
stat_filename = [curr_dir 'krug_data_02_stat.mat'];

kre_ix = 9596;
q = [0.025 0.975];


%% Calculating the slopes of KRE33-RM and KRE33-BY strains (PCA and Linear regression)

% if IFLOADSTAT & exist(stat_filename, 'file') > 0
%     load(stat_filename, 's_obs', 's_distrib', 'Nperm', 'TF');
% else
%     Nperm = 10000; % Number of permutations to estimate the confidence intervals on excess fitness effect of the KRE33-BY allele
%     TF = krug_data.fnd_gt( : , kre_ix);
%     [s_obs, s_distr] = get_stat_distrib( krug_data, TF, Nperm);
%     save(stat_filename, 's_obs', 's_distr', 'Nperm', 'TF');
% end
% 
% env_fields = {'evo1', 'evo2'};
% env = {'30°C', '37°C'};
% bg = {'KRE33-RM', 'KRE33-BY'};
% 
% fprintf('Pleiotropic slope and angle analysis\n');
% fprintf('Nperm = %d\n', Nperm);
% 
% for ienv = 0:2
%     
%     if ienv == 0
%         so = s_obs.fnd;
%         sd = s_distr.fnd;
%         fprintf('----\nFounders:\n----\n');
%     else
%         so = s_obs.(env_fields{ienv});
%         sd = s_distr.(env_fields{ienv});
%         fprintf('----\nEvolved at %s:\n----\n', env{ienv});
%     end        
%     
%     for ibg = 1:2
%         fprintf('%s:\n', bg{ibg});
%         
%         for imenv = 1:2        
%             s1 = so.m(ibg,imenv);
%             s1_distr = quantile(sd.m(ibg,imenv,:),q);
%             fprintf('\tMean fitness at %s: %.2f (%.2f, %.2f) %%\n',...
%                 env{imenv}, s1, s1_distr);
%             
%             if ienv == 0
%                 continue;
%             end
%             s1 = so.m(ibg,imenv) - s_obs.fnd.m(ibg,imenv);
%             s1_distr = quantile(sd.m(ibg,imenv,:)-s_distr.fnd.m(ibg,imenv,:), q);
%             fprintf('\tMean fitness gain at %s: %.2f (%.2f, %.2f) %%\n',...
%                 env{imenv}, s1, s1_distr);            
%         end
%         
%         s1 = so.borth(ibg,2);
%         s1_distr = quantile(sd.borth(ibg,2,:),q);
%         fprintf('\tPCA pleiotropy slope: %.2f (%.2f, %.2f)\n',...
%             s1, s1_distr);
%         
%         s1 = so.blin(ibg,2);
%         s1_distr = quantile(sd.blin(ibg,2,:),q);
%         fprintf('\tLR pleiotropy slope: %.2f (%.2f, %.2f)\n',...
%             s1, s1_distr);
%         
%         if ienv == 0
%             continue;
%         end
%         s1 = so.mpa(ibg,1);
%         s1_distr = quantile(sd.mpa(ibg,1,:),q);
%         fprintf('\tMean pleiotropic angle: %.3f (%.3f, %.3f)\n',...
%             s1, s1_distr);
%     end
% end
% clear ienv ibg imenv s1 s1_distr so sd;
% 
% for ibg = 1:2
%     s1 = s_obs.pam_diff_env(ibg,1);
%     s1_distr = quantile(s_distr.pam_diff_env(ibg,1,:),q);
%     fprintf('Mean pleiotropic angle diff btw 37°C and 30°C for %s: %.3f (%.3f, %.3f)\n',...
%         bg{ibg}, s1, s1_distr);
% end
% clear ienv ibg imenv s1 s1_distr;
% 
% for ienv = 1:2
%     s1 = s_obs.pam_diff_kre33(ienv,1);
%     s1_distr = quantile(s_distr.pam_diff_kre33(ienv,1,:),q);
%     fprintf('Mean pleiotropic angle diff btw KRE33-BY and KRE33-RM at %s: %.3f (%.3f, %.3f)\n',...
%         env{ienv}, s1, s1_distr);
% end
% clear ienv ibg imenv s1 s1_distr;











%% plot founder and population fitness YPD30 vs SC37 (4 panels); Figures 1B and 6

close all;

xmin = -17;
xd = 5;
xmax = 17;
ymin = -32;
yd = 10;
ymax = 25;

tickfontsize = 8;
labelfontsize = 10;

cc = [0,0,0;        % black
    88, 81, 147;    % Dark Blue (BY, Low T)
    150, 63, 61;    % Dark Red (BY, High T);
    130, 118, 235;  % Light Blue (RM, Low T)
    226, 115, 90;   % Light Red (RM, High T)    
    0, 114, 178;    % blue
    213, 94, 0;     % vermillion
    86, 180, 233;   % sky blue
    230 159, 0;     % orange
    204, 121, 167;   % raddish purple
    0, 158, 115;    % bluish green
    240, 228, 66   % yellow
    ]./256;

mycc = {cc(4,:), cc(5,:); cc(2,:), cc(3,:)};

KRE33TF{1} = krug_data.fnd_gt(:,kre_ix);
KRE33TF{2} = ~krug_data.fnd_gt(:,kre_ix);


%%% Figure 1B
figure;
figurewidth = 8;
figureheight = 9;

x0 = 0.14; % left margin
y0 = 0.12;

subplotwidth = 0.85;
subplotheight = 1*figurewidth * subplotwidth / figureheight;

set(gcf, 'Units', 'centimeters');
set(gcf, 'Position', [0 0 figurewidth figureheight]);

subplotposx = x0;% + subplotwidth * [0 1] + dx * [0 1]; 
subplotposy = y0;% + subplotheight * [1 0] + dy * [1 0]; 


subplot('Position', [subplotposx(1) subplotposy(1) subplotwidth subplotheight]),
hold on, box on;
set(gca, 'FontName', 'Helvetica', 'FontSize', tickfontsize, 'Layer', 'top');

X = krug_data.fnd_fit(:,1:2);
Y = krug_data.fnd_fit(:,3:4);

X = X*100;
Y = Y*100;

coeffs = pca( [X(:,1), Y(:,1)]);
b(1,2) = coeffs(2,1)/coeffs(1,1);
b(1,1) = nanmean( Y(:,1) ) - nanmean( X(:,1) ) * b(1,2);

attribs.s = 'o';
attribs.mfc = 0*[1 1 1];
attribs.mec = 'none';

h = ploterrxy(X(:,1), Y(:,1), X(:,2), Y(:,2), attribs.s, 'hhy', 0, 'hhx', 0);
set(h, 'Color', attribs.mfc, 'MarkerSize', 5, 'MarkerFaceColor', attribs.mfc, 'MarkerEdgeColor', attribs.mec);

% xrange = [min(X(:,1)), max(X(:,1))];
% plot(xrange, b(1) + b(2)*xrange, 'k', 'LineWidth', 2);


plot([-20 20], [-20 20], '-', 'Color', 0.6*[1 1 1]);

set(gca, 'XLim', [xmin, xmax], 'XTick', -15:xd:xmax, 'XGrid', 'on');
set(gca, 'YLim', [ymin, ymax], 'YTick', -30:yd:ymax, 'YGrid', 'on');

title('Founders', 'FontName', 'Arial', 'FontSize', 14);
xlabel('X, %', 'FontName', 'Arial', 'FontSize', labelfontsize);
ylabel('Y, %', 'FontName', 'Arial', 'FontSize', labelfontsize);







%%% Figure 6
figure;

figurewidth = 16;
figureheight = 15;

x0 = 0.075; % left margin
y0 = 0.065;

dx = 0.0975; % between-axes interval
dy = 0.075;

subplotwidth = 0.40;
subplotheight = 1*figurewidth * subplotwidth / figureheight;

set(gcf, 'Units', 'centimeters');
set(gcf, 'Position', [0 0 figurewidth figureheight]);

subplotposx = x0 + subplotwidth * [0 1] + dx * [0 1]; 
subplotposy = y0 + subplotheight * [1 0] + dy * [1 0]; 

% Founders:
subplot('Position', [subplotposx(1) subplotposy(1) subplotwidth subplotheight]),
hold on, box on;
set(gca, 'FontName', 'Helvetica', 'FontSize', tickfontsize, 'Layer', 'top');

for k = 1:2
    TF = KRE33TF{k};
    
    X = krug_data.fnd_fit(TF,1:2);
    Y = krug_data.fnd_fit(TF,3:4);
    
    X = X*100;
    Y = Y*100;
    
    coeffs = pca( [X(:,1), Y(:,1)]);
    b(1,2) = coeffs(2,1)/coeffs(1,1);
    b(1,1) = nanmean( Y(:,1) ) - nanmean( X(:,1) ) * b(1,2);
    
    if k == 1
        attribs.s = 'o';
        attribs.mfc = 0.7*[1 1 1];
        attribs.mec = 'none';
        xRM = nanmean( X(:,1) );
        yRM = nanmean( Y(:,1) );
    else
        attribs.s = 'o';
        attribs.mfc = 0.4*[1 1 1];
        attribs.mec = 'none';
        xBY = nanmean( X(:,1) );
        yBY = nanmean( Y(:,1) );
    end
           
    h = ploterrxy(X(:,1), Y(:,1), X(:,2), Y(:,2), attribs.s, 'hhy', 0, 'hhx', 0);
    set(h, 'Color', attribs.mfc, 'MarkerSize', 5, 'MarkerFaceColor', attribs.mfc, 'MarkerEdgeColor', attribs.mec);
    
    xrange = [min(X(:,1)), max(X(:,1))];
    % plot(xrange, b(1) + b(2)*xrange, 'k', 'LineWidth', 2);
end

plot([-20 20], [-20 20], '-', 'Color', 0.6*[1 1 1]);

set(gca, 'XLim', [xmin, xmax], 'XTick', -15:xd:xmax, 'XGrid', 'on');
set(gca, 'YLim', [ymin, ymax], 'YTick', -30:yd:ymax, 'YGrid', 'on');

text(-16, 22, 'Founders', 'FontName', 'Arial', 'FontSize', 12);
xlabel('X, %', 'FontName', 'Arial', 'FontSize', labelfontsize);
ylabel('Y, %', 'FontName', 'Arial', 'FontSize', labelfontsize);



%%% Evolved pops at 30°C and 37°C
for ienv = 1:2
    if ienv == 1
        subplot('Position', [subplotposx(1) subplotposy(2) subplotwidth subplotheight]);
    elseif ienv == 2
        subplot('Position', [subplotposx(2) subplotposy(1) subplotwidth subplotheight]);
    end
    hold on, box on;
    set(gca, 'FontName', 'Helvetica', 'FontSize', tickfontsize, 'Layer', 'top');
    
    clear xRM yRM xBY yBY bRM bBY;
    for k = 1:2
        TF = KRE33TF{k};
                
        if ienv == 1
            ix0 = 1;
        elseif ienv == 2
            ix0 = 9;
        end
        
        X = krug_data.pop_fit(TF,ix0:2:ix0+6);
        Y = krug_data.pop_fit(TF,(ix0+1):2:(ix0+7));
        X(:,1) = nanmean(X,2);
        X(:,2) = nanstd(X,0,2)./sqrt( sum(~isnan(X),2)  );
        Y(:,1) = nanmean(Y,2);
        Y(:,2) = nanstd(Y,0,2)./sqrt( sum(~isnan(Y),2)  );
        
        X = X*100;
        Y = Y*100;

        coeffs = pca( [X(:,1), Y(:,1)]);
        b(1,2) = coeffs(2,1)/coeffs(1,1);
        b(1,1) = nanmean( Y(:,1) ) - nanmean( X(:,1) ) * b(1,2);

        attribs.s = 'o';
        attribs.mfc = mycc{k,ienv};
        attribs.mec = 'none';

%         if k == 1
%             
%             xRM = nanmean( X(:,1) );
%             yRM = nanmean( Y(:,1) );
%         else
%             attribs.s = 's';
%             attribs.mfc = 0.7*[1 1 1];
%             attribs.mec = 'none';
%             xBY = nanmean( X(:,1) );
%             yBY = nanmean( Y(:,1) );
%         end

                
        h = ploterrxy(X(:,1), Y(:,1), X(:,2), Y(:,2), attribs.s, 'hhy', 0, 'hhx', 0);
        set(h, 'Color', attribs.mfc, 'MarkerSize', 5, 'MarkerFaceColor', attribs.mfc, 'MarkerEdgeColor', attribs.mec);
        
        xrange = [min(X(:,1)), max(X(:,1))];
        % plot(xrange, b(1) + b(2)*xrange, 'k', 'LineWidth', 2);
    end
    plot([-20 20], [-20 20], '-', 'Color', 0.6*[1 1 1]);
    
    set(gca, 'XLim', [xmin, xmax], 'XTick', -15:xd:xmax, 'XGrid', 'on');
    set(gca, 'YLim', [ymin, ymax], 'YTick', -30:yd:ymax, 'YGrid', 'on');
    
    if ienv == 1
        %set(gca, 'YTickLabel', {});
        text(-16, 22, 'Evolved at OT', 'FontName', 'Arial', 'FontSize', 12);
    elseif ienv == 2
        %set(gca, 'YTickLabel', {});
        text(-16, 22, 'Evolved at HT', 'FontName', 'Arial', 'FontSize', 12);
    end

    xlabel('X, %', 'FontName', 'Arial', 'FontSize', labelfontsize);
    ylabel('Y, %', 'FontName', 'Arial', 'FontSize', labelfontsize);
end


% DY vs DX (RM)

X0 = krug_data.fnd_fit(:,1:2) * 100;
Y0 = krug_data.fnd_fit(:,3:4) * 100 ;

X1 = krug_data.pop_fit(:,1:2:7) * 100;
D1X(:,1) = nanmean(X1,2) - X0(:,1);
D1X(:,2) = sqrt( nanstd(X1,0,2).^2./sum(~isnan(X1),2)  + X0(:,2).^2 );

Y1 = krug_data.pop_fit(:, 2:2:8) * 100;
D1Y(:,1) = nanmean(Y1,2) - Y0(:,1);
D1Y(:,2) = sqrt( nanstd(Y1,0,2).^2./sum(~isnan(Y1),2)  + Y0(:,2).^2 );

X2 = krug_data.pop_fit(:,9:2:15) * 100;
D2X(:,1) = nanmean(X2,2) - X0(:,1);
D2X(:,2) = sqrt( nanstd(X2,0,2).^2./sum(~isnan(X2),2)  + X0(:,2).^2 );

Y2 = krug_data.pop_fit(:, 10:2:16) * 100;
D2Y(:,1) = nanmean(Y2,2) - Y0(:,1);
D2Y(:,2) = sqrt( nanstd(Y2,0,2).^2./sum(~isnan(Y2),2)  + Y0(:,2).^2 );


subplot('Position', [subplotposx(2) subplotposy(2) subplotwidth subplotheight]);
hold on, box on;
set(gca, 'FontName', 'Helvetica', 'FontSize', tickfontsize, 'Layer', 'top');

for k = 1:2
    TF = KRE33TF{k};
    
    attribs.s = 'o';
    attribs.mec = 'none';
    attribs.ms = 5;
    
    %%%% Evolved at 30°C:
%    plot( [zeros(nnz(TF),1) D1X(TF,1)]', [zeros(nnz(TF),1) D1Y(TF,1)]', '-', 'Color', mycc{k,1}, 'LineWidth', 1);
    
   attribs.mfc = mycc{k,1};
    plot(D1X(TF,1), D1Y(TF,1), attribs.s, 'Color', attribs.mfc, 'MarkerSize', attribs.ms, 'MarkerFaceColor', attribs.mfc, 'MarkerEdgeColor', attribs.mec);
%     %         h = ploterrxy(D1X(TF,1), D1Y(TF,1), D1X(TF,2), D1Y(TF,2), attribs.s, 'hhy', 0, 'hhx', 0);
%     %         set(h, 'Color', attribs.mfc, 'MarkerSize', 4, 'MarkerFaceColor', attribs.mfc, 'MarkerEdgeColor', attribs.mec);
    
    
    %%%% Evolved at 37°C:
%    plot( [zeros(nnz(TF),1) D2X(TF,1)]', [zeros(nnz(TF),1) D2Y(TF,1)]', '-', 'Color', mycc{k,2}, 'LineWidth', 1);

   attribs.mfc = mycc{k,2};
    plot(D2X(TF,1), D2Y(TF,1), attribs.s, 'Color', attribs.mfc, 'MarkerSize', attribs.ms, 'MarkerFaceColor', attribs.mfc, 'MarkerEdgeColor', attribs.mec);
%     %         h = ploterrxy(D2X(TF,1), D2Y(TF,1), D2X(TF,2), D2Y(TF,2), attribs.s, 'hhy', 0, 'hhx', 0);
%     %         set(h, 'Color', attribs.mfc, 'MarkerSize', 4, 'MarkerFaceColor', attribs.mfc, 'MarkerEdgeColor', attribs.mec);
end

for k = 1:2
    TF = KRE33TF{k};
    plot( [0 mean(D1X(TF,1))], [0 mean(D1Y(TF,1))], '-', 'Color', mycc{k,1}, 'LineWidth', 4);
    plot( [0 mean(D2X(TF,1))], [0 mean(D2Y(TF,1))], '-', 'Color', mycc{k,2}, 'LineWidth', 4);
end
clear k TF attribs h mycc;

xrange = -2:1:30;
plot(xrange, xrange, '-', 'Color', 0.6*[1 1 1], 'LineWidth', 1);

set(gca, 'XLim', [-1, 21], 'XTick', 0:5:20, 'XGrid', 'on');
set(gca, 'YLim', [-12, 42], 'YTick', -10:10:40, 'YGrid', 'on');

xlabel('\DeltaX, %', 'FontName', 'Arial', 'FontSize', labelfontsize);
ylabel('\DeltaY, %', 'FontName', 'Arial', 'FontSize', labelfontsize);


annotation('textbox', [0   0.90 0.1 0.1], 'String', 'A', 'FontName', 'Arial', 'FontSize', 12, 'FontWeight', 'bold', 'LineStyle', 'none');
annotation('textbox', [0   0.40 0.1 0.1], 'String', 'B', 'FontName', 'Arial', 'FontSize', 12, 'FontWeight', 'bold', 'LineStyle', 'none');
annotation('textbox', [0.5 0.90 0.1 0.1], 'String', 'C', 'FontName', 'Arial', 'FontSize', 12, 'FontWeight', 'bold', 'LineStyle', 'none');
annotation('textbox', [0.5 0.40 0.1 0.1], 'String', 'D', 'FontName', 'Arial', 'FontSize', 12, 'FontWeight', 'bold', 'LineStyle', 'none');



clear xRM yRM xBY yBY d bBY bRM KRE33 xrange;
clear figurewidth figureheight tickfontsize labelfontsize cc;
clear x0 y0 dx dy subplotwidth subplotheight subplotposx subplotposy;
clear xmin xmax ymin ymax xd yd attribs X Y b bint r rint stats ienv k ix0 TF h;
clear X0 Y0 X1 X2 Y1 Y2 D1X D1Y D2X D2Y KRE33TF;






%% plot founder and population fitness YPD30 vs SC37 (6 panels)

% close all;
% 
% figure;
% 
% xmin = -17;
% xd = 5;
% xmax = 17;
% ymin = -32;
% yd = 10;
% ymax = 25;
% 
% 
% figurewidth = 18;
% figureheight = 12;
% 
% tickfontsize = 8;
% labelfontsize = 10;
% 
% x0 = 0.07; % left margin
% y0 = 0.10;
% 
% dx = 0.08; % between-axes interval
% dy = 0.12;
% 
% subplotwidth = 0.255;
% subplotheight = 1*figurewidth * subplotwidth / figureheight;
% 
% set(gcf, 'Units', 'centimeters');
% set(gcf, 'Position', [0 0 figurewidth figureheight]);
% 
% cc = [0,0,0;        % black
%     0, 114, 178;    % blue
%     213, 94, 0;     % vermillion
%     86, 180, 233;   % sky blue
%     230 159, 0;     % orange
%     204, 121, 167;   % raddish purple
%     0, 158, 115;    % bluish green
%     240, 228, 66   % yellow
%     ]./256;
% 
% mycc = {cc(4,:), cc(5,:); cc(2,:), cc(3,:)};
% 
% subplotposx = x0 + subplotwidth * [0 1 2] + dx * [0 1 2]; 
% subplotposy = y0 + subplotheight * [1 0] + dy * [1 0]; 
% 
% 
% KRE33TF{1} = krug_data.fnd_gt(:,kre_ix);
% KRE33TF{2} = ~krug_data.fnd_gt(:,kre_ix);
% 
% 
% % Founders:
% subplot('Position', [subplotposx(1) subplotposy(1) subplotwidth subplotheight]),
% hold on, box on;
% set(gca, 'FontName', 'Helvetica', 'FontSize', tickfontsize, 'Layer', 'top');
% 
% for k = 1:2
%     TF = KRE33TF{k};
%     
%     X = krug_data.fnd_fit(TF,1:2);
%     Y = krug_data.fnd_fit(TF,3:4);
%     
%     X = X*100;
%     Y = Y*100;
%     
%     coeffs = pca( [X(:,1), Y(:,1)]);
%     b(1,2) = coeffs(2,1)/coeffs(1,1);
%     b(1,1) = nanmean( Y(:,1) ) - nanmean( X(:,1) ) * b(1,2);
%     
%     if k == 1
%         attribs.s = 'o';
%         attribs.mfc = 0.7*[1 1 1];
%         attribs.mec = 'none';
%         xRM = nanmean( X(:,1) );
%         yRM = nanmean( Y(:,1) );
%     else
%         attribs.s = 'o';
%         attribs.mfc = 0.4*[1 1 1];
%         attribs.mec = 'none';
%         xBY = nanmean( X(:,1) );
%         yBY = nanmean( Y(:,1) );
%     end
%            
%     h = ploterrxy(X(:,1), Y(:,1), X(:,2), Y(:,2), attribs.s, 'hhy', 0, 'hhx', 0);
%     set(h, 'Color', attribs.mfc, 'MarkerSize', 4, 'MarkerFaceColor', attribs.mfc, 'MarkerEdgeColor', attribs.mec);
%     
%     xrange = [min(X(:,1)), max(X(:,1))];
%     plot(xrange, b(1) + b(2)*xrange, 'k', 'LineWidth', 2);
% end
% 
% plot([-20 20], [-20 20], '-', 'Color', 0.6*[1 1 1]);
% 
% set(gca, 'XLim', [xmin, xmax], 'XTick', -15:xd:xmax, 'XGrid', 'on');
% set(gca, 'YLim', [ymin, ymax], 'YTick', -30:yd:ymax, 'YGrid', 'on');
% 
% text(-16, 21, 'Founders', 'FontName', 'Arial', 'FontSize', 12);
% xlabel('Fitness at 30°C (X), %', 'FontName', 'Arial', 'FontSize', labelfontsize);
% ylabel('Fitness at 37°C (Y), %', 'FontName', 'Arial', 'FontSize', labelfontsize);
% 
% 
% 
% %%% Evolved pops at 30°C and 37°C
% for ienv = 1:2
%     subplot('Position', [subplotposx(1+ienv) subplotposy(1) subplotwidth subplotheight]);
%     hold on, box on;
%     set(gca, 'FontName', 'Helvetica', 'FontSize', tickfontsize, 'Layer', 'top');
%     
%     clear xRM yRM xBY yBY bRM bBY;
%     for k = 1:2
%         TF = KRE33TF{k};
%                 
%         if ienv == 1
%             ix0 = 1;
%         elseif ienv == 2
%             ix0 = 9;
%         end
%         
%         X = krug_data.pop_fit(TF,ix0:2:ix0+6);
%         Y = krug_data.pop_fit(TF,(ix0+1):2:(ix0+7));
%         X(:,1) = nanmean(X,2);
%         X(:,2) = nanstd(X,0,2)./sqrt( sum(~isnan(X),2)  );
%         Y(:,1) = nanmean(Y,2);
%         Y(:,2) = nanstd(Y,0,2)./sqrt( sum(~isnan(Y),2)  );
%         
%         X = X*100;
%         Y = Y*100;
% 
%         coeffs = pca( [X(:,1), Y(:,1)]);
%         b(1,2) = coeffs(2,1)/coeffs(1,1);
%         b(1,1) = nanmean( Y(:,1) ) - nanmean( X(:,1) ) * b(1,2);
% 
%         attribs.s = 'o';
%         attribs.mfc = mycc{k,ienv};
%         attribs.mec = 'none';
%                 
%         h = ploterrxy(X(:,1), Y(:,1), X(:,2), Y(:,2), attribs.s, 'hhy', 0, 'hhx', 0);
%         set(h, 'Color', attribs.mfc, 'MarkerSize', 4, 'MarkerFaceColor', attribs.mfc, 'MarkerEdgeColor', attribs.mec);
%         
%         xrange = [min(X(:,1)), max(X(:,1))];
%         plot(xrange, b(1) + b(2)*xrange, 'k', 'LineWidth', 2);
%     end
%     plot([-20 20], [-20 20], '-', 'Color', 0.6*[1 1 1]);
%     
%     set(gca, 'XLim', [xmin, xmax], 'XTick', -15:xd:xmax, 'XGrid', 'on');
%     set(gca, 'YLim', [ymin, ymax], 'YTick', -30:yd:ymax, 'YGrid', 'on');
%     
%     if ienv == 1
%         %set(gca, 'YTickLabel', {});
%         text(-16, 21, 'Evolved at 30°C', 'FontName', 'Arial', 'FontSize', 12);
%     elseif ienv == 2
%         %set(gca, 'YTickLabel', {});
%         text(-16, 21, 'Evolved at 37°C', 'FontName', 'Arial', 'FontSize', 12);
%     end
% 
%     xlabel('Fitness at 30°C (X), %', 'FontName', 'Arial', 'FontSize', labelfontsize);
%     ylabel('Fitness at 37°C (Y), %', 'FontName', 'Arial', 'FontSize', labelfontsize);
% end
% 
% 
% % DY vs DX (RM)
% 
% X0 = krug_data.fnd_fit(:,1:2) * 100;
% Y0 = krug_data.fnd_fit(:,3:4) * 100 ;
% 
% X1 = krug_data.pop_fit(:,1:2:7) * 100;
% D1X(:,1) = nanmean(X1,2) - X0(:,1);
% D1X(:,2) = sqrt( nanstd(X1,0,2).^2./sum(~isnan(X1),2)  + X0(:,2).^2 );
% 
% Y1 = krug_data.pop_fit(:, 2:2:8) * 100;
% D1Y(:,1) = nanmean(Y1,2) - Y0(:,1);
% D1Y(:,2) = sqrt( nanstd(Y1,0,2).^2./sum(~isnan(Y1),2)  + Y0(:,2).^2 );
% 
% X2 = krug_data.pop_fit(:,9:2:15) * 100;
% D2X(:,1) = nanmean(X2,2) - X0(:,1);
% D2X(:,2) = sqrt( nanstd(X2,0,2).^2./sum(~isnan(X2),2)  + X0(:,2).^2 );
% 
% Y2 = krug_data.pop_fit(:, 10:2:16) * 100;
% D2Y(:,1) = nanmean(Y2,2) - Y0(:,1);
% D2Y(:,2) = sqrt( nanstd(Y2,0,2).^2./sum(~isnan(Y2),2)  + Y0(:,2).^2 );
% 
% attribs.s = 'o';
% attribs.mec = 'none';
% attribs.ms = 4;
% 
% for k = 1:2
%     subplot('Position', [subplotposx(k) subplotposy(2) subplotwidth subplotheight]);
%     hold on, box on;
%     set(gca, 'FontName', 'Helvetica', 'FontSize', tickfontsize, 'Layer', 'top');
%     
%     TF = KRE33TF{k};
%         
%     %%%% Evolved at 30°C:
%     attribs.mfc = mycc{k,1};
%     plot(D1X(TF,1), D1Y(TF,1), attribs.s, 'Color', attribs.mfc, 'MarkerSize', attribs.ms, 'MarkerFaceColor', attribs.mfc, 'MarkerEdgeColor', attribs.mec);    
%     
%     %%%% Evolved at 37°C:    
%     attribs.mfc = mycc{k,2};
%     plot(D2X(TF,1), D2Y(TF,1), attribs.s, 'Color', attribs.mfc, 'MarkerSize', attribs.ms, 'MarkerFaceColor', attribs.mfc, 'MarkerEdgeColor', attribs.mec);
%     
%     plot( [0 mean(D1X(TF,1))], [0 mean(D1Y(TF,1))], '-', 'Color', mycc{k,1}, 'LineWidth', 4);
%     plot( [0 mean(D2X(TF,1))], [0 mean(D2Y(TF,1))], '-', 'Color', mycc{k,2}, 'LineWidth', 4);
%     
%     xrange = -2:1:30;
%     plot(xrange, xrange, '-', 'Color', 0.6*[1 1 1], 'LineWidth', 1);
%     
%     set(gca, 'XLim', [-1, 21], 'XTick', 0:5:20, 'XGrid', 'on');
%     set(gca, 'YLim', [-12, 42], 'YTick', -10:10:40, 'YGrid', 'on');
%     
%     xlabel('\DeltaX, %', 'FontName', 'Arial', 'FontSize', labelfontsize);
%     ylabel('\DeltaY, %', 'FontName', 'Arial', 'FontSize', labelfontsize);
% end
% 
% clear k TF attribs h mycc;
% 
% 
% % Distribution of pleiotropic angles
% angle_mesh = -3*pi/16:pi/8:pi;
% counts = nan( 2, length(angle_mesh)-1 );
% for k = 1:2
%     
%     TF = KRE33TF{k};
%     
%     alpha_vec = get_pleiotropic_angle( [D1X(TF,1), D1Y(TF,1)] , [D2X(TF,1), D2Y(TF,1)] );
%     
%     counts(k,:) = histcounts( alpha_vec, angle_mesh);
% end
% 
% subplot('Position', [subplotposx(3) subplotposy(2) subplotwidth subplotheight]);
% hold on, box on;
% set(gca, 'FontName', 'Helvetica', 'FontSize', tickfontsize, 'Layer', 'top');
% 
% x = (angle_mesh(1:end-1) + angle_mesh(2:end) )/2;
% 
% h = bar( x, counts', 0.8, 'EdgeColor', 'none' );
% set(h(1), 'FaceColor', 0.7*[1 1 1]);
% set(h(2), 'FaceColor', 0.4*[1 1 1]);
% set(gca, 'XLim', [-0.2 2.6], 'XTick', 0:pi/4:3/4*pi, 'TickLabelInterpreter', 'latex',...
%     'XTickLabel', {'$0$', '$\pi/4$', '$\pi/2$', '$3\pi/4$'} );
% 
% xlabel('Pleiotropic angle, rad', 'FontName', 'Arial', 'FontSize', labelfontsize);
% ylabel('Number of founders', 'FontName', 'Arial', 'FontSize', labelfontsize);
% legend({'RM', 'BY'});
% 
% clear angle_mesh counts k TF alpha_vec x h;
% 
% annotation('textbox', [0   0.9 0.1 0.1], 'String', 'A', 'FontName', 'Arial', 'FontSize', 12, 'FontWeight', 'bold', 'LineStyle', 'none');
% annotation('textbox', [0.335 0.9 0.1 0.1], 'String', 'B', 'FontName', 'Arial', 'FontSize', 12, 'FontWeight', 'bold', 'LineStyle', 'none');
% annotation('textbox', [0.67 0.9 0.1 0.1], 'String', 'C', 'FontName', 'Arial', 'FontSize', 12, 'FontWeight', 'bold', 'LineStyle', 'none');
% annotation('textbox', [0   0.4 0.1 0.1], 'String', 'D', 'FontName', 'Arial', 'FontSize', 12, 'FontWeight', 'bold', 'LineStyle', 'none');
% annotation('textbox', [0.335 0.4 0.1 0.1], 'String', 'E', 'FontName', 'Arial', 'FontSize', 12, 'FontWeight', 'bold', 'LineStyle', 'none');
% annotation('textbox', [0.67 0.4 0.1 0.1], 'String', 'F', 'FontName', 'Arial', 'FontSize', 12, 'FontWeight', 'bold', 'LineStyle', 'none');
% 
% 
% 
% clear xRM yRM xBY yBY d bBY bRM KRE33 xrange;
% clear figurewidth figureheight tickfontsize labelfontsize cc;
% clear x0 y0 dx dy subplotwidth subplotheight subplotposx subplotposy;
% clear xmin xmax ymin ymax xd yd attribs X Y b bint r rint stats ienv k ix0 TF h;
% clear X0 Y0 X1 X2 Y1 Y2 D1X D1Y D2X D2Y KRE33TF;






%% plot Delta 1 Y vs Y and Delta 2 X vs X

% alpha = 0.01;
%     
% clf;
% 
% figurewidth = 18;
% figureheight = 11;
% 
% tickfontsize = 8;
% labelfontsize = 10;
% 
% x0 = 0.07; % left margin
% y0 = 0.1;
% 
% dx = 0.095; % between-axes interval
% dy = 0.03;
% 
% subplotwidth = 0.24;
% subplotheight = 1*figurewidth * subplotwidth / figureheight;
% 
% set(gcf, 'Units', 'centimeters');
% set(gcf, 'Position', [0 0 figurewidth figureheight]);
% 
% cc = [0,0,0;        % black
%     0, 114, 178;    % blue
%     213, 94, 0;     % vermillion
%     86, 180, 233;   % sky blue
%     230 159, 0;     % orange
%     204, 121, 167;   % raddish purple
%     0, 158, 115;    % bluish green
%     240, 228, 66   % yellow
%     ]./256;
% 
% 
% subplotposx = x0 + subplotwidth * [0 1 2] + dx * [0 1 2]; 
% subplotposy = y0 + subplotheight * [1 0] + dy * [1 0]; 
% 
% 
% %%%% Doing the regrassion analyses separately for RM and BY KRE33:
% clear X1 X2 Y1 Y2 X0 Y0 D1X D2X D1Y D2Y;
% for k = 1:2
%     if k == 1
%         TF = krug_data.fnd_gt(:,kre_ix);
% %         attribs.s = 'o';
% %         attribs.mfc = 'k';
% %         attribs.mec = 'none';
%     elseif k == 2
%         TF = ~krug_data.fnd_gt(:,kre_ix);
% %         attribs.s = 's';
% %         attribs.mfc = 'w';
% %         attribs.mec = 'k';
%     end
%     
%     attribs.s = 'o';
%     attribs.mfc = 0.6*[1 1 1];
%     attribs.mec = 0.6*[1 1 1];
% 
%     
%     %%% Founders
%     subplot('Position', [subplotposx(1) subplotposy(k) subplotwidth subplotheight]);
%     hold on, box on;
%     set(gca, 'FontName', 'Helvetica', 'FontSize', tickfontsize, 'Layer', 'top');
% 
%     X0 = krug_data.fnd_fit(TF,1:2);
%     Y0 = krug_data.fnd_fit(TF,3:4);
%     [bfnd,bint,r,rint,stats] = regress(Y0(:,1), [ones(size(X0,1),1), X0(:,1)]);
% %     r2 = stats(1);
% %     [b2,bint,r,rint,stats] = regress(X0(:,1), [ones(size(Y0,1),1), Y0(:,1)]);
%     
%     h = ploterrxy(X0(:,1), Y0(:,1), X0(:,2), Y0(:,2), attribs.s, 'hhy', 0, 'hhx', 0);
%     set(h, 'Color', attribs.mec, 'MarkerSize', 4, 'MarkerFaceColor', attribs.mfc, 'MarkerEdgeColor', attribs.mec);
% 
%     xrange = [min(X0(:,1)), max(X0(:,1))];
%     plot(xrange , bfnd(1) + bfnd(2)*xrange, 'k', 'LineWidth', 2);
%         
%     plot([-0.2 0.2], [-0.2 0.2], '-', 'Color', 0.3*[1 1 1]);
% 
%     set(gca, 'XLim', [-0.17, 0.17], 'XTick', -0.15:0.05:0.15, 'XGrid', 'on');
%     set(gca, 'YLim', [-0.32, 0.25], 'YTick', -0.3:0.1:0.25, 'YGrid', 'on');
%     
%     if k == 1
%         set(gca, 'XTickLabels', {});
%         title('Founders', 'FontName', 'Arial', 'FontSize', 14);
%         text(-0.16, 0.25, 'RM', 'FontName', 'Arial', 'FontSize', labelfontsize, 'VerticalAlignment', 'top');
%     elseif k == 2
%         xlabel('Fitness at 30°C, X', 'FontName', 'Arial', 'FontSize', labelfontsize);
%         text(-0.16, 0.25, 'BY', 'FontName', 'Arial', 'FontSize', labelfontsize, 'VerticalAlignment', 'top');
%     end
%     % text(-0.16, 0.20-(k-1)*0.1, sprintf('R^2 = %.2f', r2), 'FontName', 'Arial', 'FontSize', labelfontsize, 'VerticalAlignment', 'top');
% 
%     ylabel('Fitness at 37°C, Y', 'FontName', 'Arial', 'FontSize', labelfontsize);
%     clear bfnd bint r rint stats;
%     
%     
%     %%% Evolved at 30°C
%     subplot('Position', [subplotposx(2) subplotposy(k) subplotwidth subplotheight]);
%     hold on, box on;
%     set(gca, 'FontName', 'Helvetica', 'FontSize', tickfontsize, 'Layer', 'top');
% 
%     X1 = krug_data.pop_fit(TF,1:2:7);
%     D1X(:,1) = nanmean(X1,2) - X0(:,1);
%     D1X(:,2) = sqrt( nanstd(X1,0,2).^2./sum(~isnan(X1),2)  + X0(:,2).^2 );
%     
%     Y1 = krug_data.pop_fit(TF, 2:2:8);
%     D1Y(:,1) = nanmean(Y1,2) - Y0(:,1);
%     D1Y(:,2) = sqrt( nanstd(Y1,0,2).^2./sum(~isnan(Y1),2)  + Y0(:,2).^2 );
%             
%     h = ploterrxy(Y0(:,1), D1Y(:,1), Y0(:,2), D1Y(:,2), attribs.s, 'hhy', 0, 'hhx', 0);
%     set(h, 'Color', attribs.mec, 'MarkerSize', 4, 'MarkerFaceColor', attribs.mfc, 'MarkerEdgeColor', attribs.mec);
% 
%     xrange = min(Y0(:,1)):0.01:max(Y0(:,1))+0.01;    
%     
%     [bfndXDX,bint,r,rint,stats] = regress(D1X(:,1), [ones(size(X0,1),1), X0(:,1)]);    
%     bfndXDX(1) = nanmean(D1Y(:,1) - bfndXDX(2)*Y0(:,1));
%     
%     % plot([-0.32 0.25], [-0.32 0.25], '-', 'Color', 0.3*[1 1 1]);        
%     plot(xrange, bfndXDX(1) + bfndXDX(2) * xrange, '--', 'Color', 'k', 'LineWidth', 2);
%     
% %     % Linear regreassion and confidence envelope (from Stuart, Ord, Arnold p. 740, expression (32.28) )
%     [b,bint,r,rint,stats] = regress(D1Y(:,1), [ones(size(Y0,1),1), Y0(:,1)]);
%     plot(xrange, b(1) + b(2) * xrange, '-', 'Color', 'k', 'LineWidth', 2);
% 
%     s2 = stats(4);
%     n = size(Y0,1);
%     Y0m = nanmean( Y0(:,1) );
%     z = s2 * sqrt(1/n + (xrange - Y0m).^2/nansum( (Y0(:,1) - Y0m).^2 ));
%     confenv = z * tinv( 1 - alpha/2 , 1); 
% 
%     plot(xrange, b(1) + b(2) * xrange + confenv, '-', 'Color', 'k', 'LineWidth', 1);
%     plot(xrange, b(1) + b(2) * xrange - confenv, '-', 'Color', 'k', 'LineWidth', 1);
% 
%     
%     set(gca, 'XLim', [-0.32, 0.25], 'XTick', -0.3:0.1:0.25, 'XGrid', 'on');
%     set(gca, 'YLim', [-0.12, 0.32], 'YTick', -0.10:0.1:0.3, 'YGrid', 'on');
% 
%     if k == 1
%         set(gca, 'XTickLabels', {});
%         title('Evolved at 30°C', 'FontName', 'Arial', 'FontSize', 14);
%         text(0.24, 0.32, 'RM', 'FontName', 'Arial', 'FontSize', labelfontsize, 'VerticalAlignment', 'top', 'HorizontalAlignment', 'right');
%     elseif k == 2
%         xlabel('Fitness at 37°C, Y', 'FontName', 'Arial', 'FontSize', labelfontsize);
%         text(0.24, 0.32, 'BY', 'FontName', 'Arial', 'FontSize', labelfontsize, 'VerticalAlignment', 'top', 'HorizontalAlignment', 'right');
%     end
%     ylabel('\Delta_1Y', 'FontName', 'Arial', 'FontSize', labelfontsize);
% 
%     clear s2 n Y0m z confenv bfndXDX bint r rint stats;
% 
% 
%     
%     
%     %%% Evolved at 37°C
%     subplot('Position', [subplotposx(3) subplotposy(k) subplotwidth subplotheight]);
%     hold on, box on;
%     set(gca, 'FontName', 'Helvetica', 'FontSize', tickfontsize, 'Layer', 'top');
%     
%     X2 = krug_data.pop_fit(TF,9:2:15);
%     D2X(:,1) = nanmean(X2,2) - X0(:,1);
%     D2X(:,2) = sqrt( nanstd(X2,0,2).^2./sum(~isnan(X2),2)  + X0(:,2).^2 );
% 
%     Y2 = krug_data.pop_fit(TF, 10:2:16);
%     D2Y(:,1) = nanmean(Y2,2) - Y0(:,1);
%     D2Y(:,2) = sqrt( nanstd(Y2,0,2).^2./sum(~isnan(Y2),2)  + Y0(:,2).^2 );
% 
%     h = ploterrxy(X0(:,1), D2X(:,1), X0(:,2), D2X(:,2), attribs.s, 'hhy', 0, 'hhx', 0);
%     set(h, 'Color', attribs.mec, 'MarkerSize', 4, 'MarkerFaceColor', attribs.mfc, 'MarkerEdgeColor', attribs.mec);
%         
%     xrange = min(X0(:,1)):0.01:max(X0(:,1))+0.01;    
%     
%     [bfndYDY,bint,r,rint,stats] = regress(D2Y(:,1), [ones(size(Y0,1),1), Y0(:,1)]);    
%     bfndYDY(1) = nanmean(D2X(:,1) - bfndYDY(2)*X0(:,1));
%     
%     % plot([-0.32 0.25], [-0.32 0.25], '-', 'Color', 0.3*[1 1 1]);        
%     plot(xrange, bfndYDY(1) + bfndYDY(2) * xrange, '--', 'Color', 'k', 'LineWidth', 2);
% 
%     % Linear regression and confidence envelope (from Stuart, Ord, Arnold p. 740, expression (32.28) )
%     [b,bint,r,rint,stats] = regress(D2X(:,1), [ones(size(X0,1),1), X0(:,1)]);
%     plot(xrange, b(1) + b(2) * xrange, '-', 'Color', 'k', 'LineWidth', 2);
% 
%     s2 = stats(4);
%     n = size(X0,1);
%     X0m = nanmean( X0(:,1) );
%     z = s2 * sqrt(1/n + (xrange - X0m).^2/nansum( (X0(:,1) - X0m).^2 ));
%     confenv = z * tinv( 1 - alpha/2 , 1);
%     
%     plot(xrange, b(1) + b(2) * xrange + confenv, '-', 'Color', 'k', 'LineWidth', 1);
%     plot(xrange, b(1) + b(2) * xrange - confenv, '-', 'Color', 'k', 'LineWidth', 1);
%      
%     set(gca, 'XLim', [-0.17, 0.17], 'XTick', -0.15:0.05:0.15, 'XGrid', 'on');
%     set(gca, 'YLim', [-0.02, 0.22], 'YTick', 0:0.05:0.2, 'YGrid', 'on');
% 
%     if k == 1
%         set(gca, 'XTickLabels', {});
%         title('Evolved at 37°C', 'FontName', 'Arial', 'FontSize', 14);
%         text(0.16, 0.22, 'RM', 'FontName', 'Arial', 'FontSize', labelfontsize, 'VerticalAlignment', 'top', 'HorizontalAlignment', 'right');
%     elseif k == 2
%         xlabel('Fitness at 30°C, X', 'FontName', 'Arial', 'FontSize', labelfontsize);
%         text(0.16, 0.22, 'BY', 'FontName', 'Arial', 'FontSize', labelfontsize, 'VerticalAlignment', 'top', 'HorizontalAlignment', 'right');
%     end
%     ylabel('\Delta_2X', 'FontName', 'Arial', 'FontSize', labelfontsize); 
%     
%     clear s2 n X0m z confenv bfndYDY bint r rint stats;
%     clear X1 X2 Y1 Y2 X0 Y0 D1X D2X D1Y D2Y;
% end
% %
% %     if ienv == 1
% %         xlabel('Fitness, 30°C', 'FontName', 'Arial', 'FontSize', labelfontsize);
% %         text(-0.16, 0.22, 'Evolved at 30°C', 'FontName', 'Arial', 'FontSize', 12);
% %     else
% %         text(-0.16, 0.22, 'Evolved at 37°C', 'FontName', 'Arial', 'FontSize', 12);
% %     end
% 
% 
% clear figurewidth figureheight tickfontsize labelfontsize cc;
% clear x0 y0 dx dy subplotwidth subplotheight subplotposx subplotposy;
% clear hplot k TF attribs X0 Y0 b1 bint r rint stats r2 b2 h Y1pred Y1obs X2pred X2obs;   











%% plot Delta X vs Delta Y

% alpha = 0.01;
%     
% clf;
% 
% figurewidth = 18;
% figureheight = 11;
% 
% tickfontsize = 8;
% labelfontsize = 10;
% 
% x0 = 0.07; % left margin
% y0 = 0.1;
% 
% dx = 0.095; % between-axes interval
% dy = 0.03;
% 
% subplotwidth = 0.24;
% subplotheight = 1*figurewidth * subplotwidth / figureheight;
% 
% set(gcf, 'Units', 'centimeters');
% set(gcf, 'Position', [0 0 figurewidth figureheight]);
% 
% cc = [0,0,0;        % black
%     0, 114, 178;    % blue
%     213, 94, 0;     % vermillion
%     86, 180, 233;   % sky blue
%     230 159, 0;     % orange
%     204, 121, 167;   % raddish purple
%     0, 158, 115;    % bluish green
%     240, 228, 66   % yellow
%     ]./256;
% 
% 
% subplotposx = x0 + subplotwidth * [0 1 2] + dx * [0 1 2]; 
% subplotposy = y0 + subplotheight * [1 0] + dy * [1 0]; 
% 
% 
% %%%% Doing the regrassion analyses separately for RM and BY KRE33:
% for k = 1:2
%     if k == 1
%         TF = krug_data.fnd_gt(:,kre_ix);
% %         attribs.s = 'o';
% %         attribs.mfc = 'k';
% %         attribs.mec = 'none';
%     elseif k == 2
%         TF = ~krug_data.fnd_gt(:,kre_ix);
% %         attribs.s = 's';
% %         attribs.mfc = 'w';
% %         attribs.mec = 'k';
%     end
%     
%     attribs.s = 'o';
%     attribs.mfc = 0.6*[1 1 1];
%     attribs.mec = 0.6*[1 1 1];
% 
%     
%     %%% Founders
%     subplot('Position', [subplotposx(1) subplotposy(k) subplotwidth subplotheight]);
%     hold on, box on;
%     set(gca, 'FontName', 'Helvetica', 'FontSize', tickfontsize, 'Layer', 'top');
% 
%     X0 = krug_data.fnd_fit(TF,1:2);
%     Y0 = krug_data.fnd_fit(TF,3:4);
%     [bfnd,bint,r,rint,stats] = regress(Y0(:,1), [ones(size(X0,1),1), X0(:,1)]);
% %     r2 = stats(1);
% %     [b2,bint,r,rint,stats] = regress(X0(:,1), [ones(size(Y0,1),1), Y0(:,1)]);
%     
%     h = ploterrxy(X0(:,1), Y0(:,1), X0(:,2), Y0(:,2), attribs.s, 'hhy', 0, 'hhx', 0);
%     set(h, 'Color', attribs.mec, 'MarkerSize', 4, 'MarkerFaceColor', attribs.mfc, 'MarkerEdgeColor', attribs.mec);
%     xrange = [min(X0(:,1)), max(X0(:,1))];
%     plot(xrange , bfnd(1) + bfnd(2)*xrange, 'k', 'LineWidth', 2);
%         
%     plot([-0.2 0.2], [-0.2 0.2], '-', 'Color', 0.3*[1 1 1]);
% 
%     set(gca, 'XLim', [-0.17, 0.17], 'XTick', -0.15:0.05:0.15, 'XGrid', 'on');
%     set(gca, 'YLim', [-0.32, 0.25], 'YTick', -0.3:0.1:0.25, 'YGrid', 'on');
%     
%     if k == 1
%         set(gca, 'XTickLabels', {});
%         title('Founders', 'FontName', 'Arial', 'FontSize', 14);
%         text(-0.16, 0.25, 'RM', 'FontName', 'Arial', 'FontSize', labelfontsize, 'VerticalAlignment', 'top');
%     elseif k == 2
%         xlabel('Fitness at 30°C, X', 'FontName', 'Arial', 'FontSize', labelfontsize);
%         text(-0.16, 0.25, 'BY', 'FontName', 'Arial', 'FontSize', labelfontsize, 'VerticalAlignment', 'top');
%     end
%     % text(-0.16, 0.20-(k-1)*0.1, sprintf('R^2 = %.2f', r2), 'FontName', 'Arial', 'FontSize', labelfontsize, 'VerticalAlignment', 'top');
% 
%     ylabel('Fitness at 37°C, Y', 'FontName', 'Arial', 'FontSize', labelfontsize);
%     
%     
%     
%     %%% Evolved at 30°C
%     subplot('Position', [subplotposx(2) subplotposy(k) subplotwidth subplotheight]);
%     hold on, box on;
%     set(gca, 'FontName', 'Helvetica', 'FontSize', tickfontsize, 'Layer', 'top');
%     
%     X1 = krug_data.pop_fit(TF,1:2:7);
%     D1X(:,1) = nanmean(X1,2) - X0(:,1);
%     D1X(:,2) = sqrt( nanstd(X1,0,2).^2./sum(~isnan(X1),2)  + X0(:,2).^2 );
%     
%     Y1 = krug_data.pop_fit(TF, 2:2:8);
%     D1Y(:,1) = nanmean(Y1,2) - Y0(:,1);
%     D1Y(:,2) = sqrt( nanstd(Y1,0,2).^2./sum(~isnan(Y1),2)  + Y0(:,2).^2 );
%             
%     h = ploterrxy(D1X(:,1), D1Y(:,1), D1X(:,2), D1Y(:,2), attribs.s, 'hhy', 0, 'hhx', 0);
%     set(h, 'Color', attribs.mec, 'MarkerSize', 4, 'MarkerFaceColor', attribs.mfc, 'MarkerEdgeColor', attribs.mec);
% 
%     xrange = min(D1X(:,1)):0.01:max(D1X(:,1))+0.01;    
%     bfnd(1) = nanmean(D1Y(:,1) - bfnd(2)*D1X(:,1));
%     
%     plot([-0.01 0.21], [-0.01 0.21], '-', 'Color', 0.3*[1 1 1]);        
%     plot(xrange, bfnd(1) + bfnd(2) * xrange, '--', 'Color', 'k', 'LineWidth', 2);
%     
%     % Linear regreassion and confidence envelope (from Stuart, Ord, Arnold p. 740, expression (32.28) )
%     [b,bint,r,rint,stats] = regress(D1Y(:,1), [ones(size(D1X,1),1), D1X(:,1)]);
%     plot(xrange, b(1) + b(2) * xrange, '-', 'Color', 'k', 'LineWidth', 2);
% 
%     s2 = stats(4);
%     n = size(D1X,1);
%     D1Xm = nanmean( D1X(:,1) );
%     z = s2 * sqrt(1/n + (xrange - D1Xm).^2/nansum( (D1X(:,1) - D1Xm).^2 ));
%     confenv = z * tinv( 1 - alpha/2 , 1); 
% 
%     plot(xrange, b(1) + b(2) * xrange + confenv, '-', 'Color', 'k', 'LineWidth', 1);
%     plot(xrange, b(1) + b(2) * xrange - confenv, '-', 'Color', 'k', 'LineWidth', 1);
% 
%     
%     set(gca, 'XLim', [-0.01, 0.21], 'XTick',  0:0.05:0.2, 'XGrid', 'on');
%     set(gca, 'YLim', [-0.12, 0.32], 'YTick', -0.10:0.1:0.3, 'YGrid', 'on');
% 
% %     text(-0.24, 0.30, 'Evolved at 30°C', 'FontName', 'Arial', 'FontSize', 12, 'VerticalAlignment', 'top');
% %     text(-0.24, 0.20-(k-1)*0.1, sprintf('R^2 = %.2f', r2), 'FontName', 'Arial', 'FontSize', labelfontsize, 'VerticalAlignment', 'top');
% 
%     if k == 1
%         set(gca, 'XTickLabels', {});
%         title('Evolved at 30°C', 'FontName', 'Arial', 'FontSize', 14);
%         text(0, 0.32, 'RM', 'FontName', 'Arial', 'FontSize', labelfontsize, 'VerticalAlignment', 'top');
%     elseif k == 2
%         xlabel('\Delta_1X', 'FontName', 'Arial', 'FontSize', labelfontsize);
%         text(0, 0.32, 'BY', 'FontName', 'Arial', 'FontSize', labelfontsize, 'VerticalAlignment', 'top');
%     end
%     ylabel('\Delta_1Y', 'FontName', 'Arial', 'FontSize', labelfontsize);
% 
%     clear X1 D1X D1Y s2 n D1Xm z confenv;
% 
%     
%     
%     %%% Evolved at 37°C
%     subplot('Position', [subplotposx(3) subplotposy(k) subplotwidth subplotheight]);
%     hold on, box on;
%     set(gca, 'FontName', 'Helvetica', 'FontSize', tickfontsize, 'Layer', 'top');
%     
%     X2 = krug_data.pop_fit(TF,9:2:15);
%     D2X(:,1) = nanmean(X2,2) - X0(:,1);
%     D2X(:,2) = sqrt( nanstd(X2,0,2).^2./sum(~isnan(X2),2)  + X0(:,2).^2 );
%     
%     Y2 = krug_data.pop_fit(TF, 10:2:16);
%     D2Y(:,1) = nanmean(Y2,2) - Y0(:,1);
%     D2Y(:,2) = sqrt( nanstd(Y2,0,2).^2./sum(~isnan(Y2),2)  + Y0(:,2).^2 );
%             
%     h = ploterrxy(D2X(:,1), D2Y(:,1), D2X(:,2), D2Y(:,2), attribs.s, 'hhy', 0, 'hhx', 0);
%     set(h, 'Color', attribs.mec, 'MarkerSize', 4, 'MarkerFaceColor', attribs.mfc, 'MarkerEdgeColor', attribs.mec);
%         
%     xrange = min(D2X(:,1)):0.01:max(D2X(:,1))+0.01;
%     bfnd(1) = nanmean(D2Y(:,1) - bfnd(2)*D2X(:,1));
% 
%     plot([-0.01 0.21], [-0.01 0.21], '-', 'Color', 0.3*[1 1 1]);
%     plot(xrange, bfnd(1) + bfnd(2) * xrange, '--', 'Color', 'k', 'LineWidth', 2);
% 
%     % Linear regression and confidence envelope (from Stuart, Ord, Arnold p. 740, expression (32.28) )
%     [b,bint,r,rint,stats] = regress(D2Y(:,1), [ones(size(D2X,1),1), D2X(:,1)]);
%     plot(xrange, b(1) + b(2) * xrange, '-', 'Color', 'k', 'LineWidth', 2);
% 
%     s2 = stats(4);
%     n = size(D2X,1);
%     D2Xm = nanmean( D2X(:,1) );
%     z = s2 * sqrt(1/n + (xrange - D2Xm).^2/nansum( (D2X(:,1) - D2Xm).^2 ));
%     confenv = z * tinv( 1 - alpha/2 , 1);
%     
%     plot(xrange, b(1) + b(2) * xrange + confenv, '-', 'Color', 'k', 'LineWidth', 1);
%     plot(xrange, b(1) + b(2) * xrange - confenv, '-', 'Color', 'k', 'LineWidth', 1);
%     
%     
%     set(gca, 'XLim', [-0.01, 0.21], 'XTick',  0:0.05:0.2, 'XGrid', 'on');
%     set(gca, 'YLim', [-0.02, 0.42], 'YTick', -0:0.1:0.4, 'YGrid', 'on');
% 
% %     text(-0.24, 0.30, 'Evolved at 30°C', 'FontName', 'Arial', 'FontSize', 12, 'VerticalAlignment', 'top');
% %     text(-0.24, 0.20-(k-1)*0.1, sprintf('R^2 = %.2f', r2), 'FontName', 'Arial', 'FontSize', labelfontsize, 'VerticalAlignment', 'top');
% 
%     if k == 1
%         set(gca, 'XTickLabels', {});
%         title('Evolved at 37°C', 'FontName', 'Arial', 'FontSize', 14);
%         text(0, 0.42, 'RM', 'FontName', 'Arial', 'FontSize', labelfontsize, 'VerticalAlignment', 'top');
%     elseif k == 2
%         xlabel('\Delta_2X', 'FontName', 'Arial', 'FontSize', labelfontsize);
%         text(0, 0.42, 'BY', 'FontName', 'Arial', 'FontSize', labelfontsize, 'VerticalAlignment', 'top');
%     end
%     ylabel('\Delta_2Y', 'FontName', 'Arial', 'FontSize', labelfontsize); 
%     
%     clear X2 D2X D2Y s2 n D2Xm z confenv;
% end
% %
% %     if ienv == 1
% %         xlabel('Fitness, 30°C', 'FontName', 'Arial', 'FontSize', labelfontsize);
% %         text(-0.16, 0.22, 'Evolved at 30°C', 'FontName', 'Arial', 'FontSize', 12);
% %     else
% %         text(-0.16, 0.22, 'Evolved at 37°C', 'FontName', 'Arial', 'FontSize', 12);
% %     end
% 
% 
% clear figurewidth figureheight tickfontsize labelfontsize cc;
% clear x0 y0 dx dy subplotwidth subplotheight subplotposx subplotposy;
% clear hplot k TF attribs X0 Y0 b1 bint r rint stats r2 b2 h Y1pred Y1obs X2pred X2obs;   






%% plot phase space
    
% clf;
% 
% figurewidth = 16;
% figureheight = 8;
% 
% tickfontsize = 8;
% labelfontsize = 10;
% 
% x0 = 0.09; % left margin
% y0 = 0.13;
% 
% dx = 0.025; % between-axes interval
% dy = 0.03;
% 
% subplotwidth = 0.43;
% subplotheight = 1*figurewidth * subplotwidth / figureheight;
% 
% set(gcf, 'Units', 'centimeters');
% set(gcf, 'Position', [0 0 figurewidth figureheight]);
% 
% cc = [0,0,0;        % black
%     0, 114, 178;    % blue
%     213, 94, 0;     % vermillion
%     86, 180, 233;   % sky blue
%     230 159, 0;     % orange
%     204, 121, 167;   % raddish purple
%     0, 158, 115;    % bluish green
%     240, 228, 66   % yellow
%     ]./256;
% 
% 
% subplotposx = x0 + subplotwidth * [0 1] + dx * [0 1]; 
% subplotposy = y0;% + subplotheight * [1 0] + dy * [1 0]; 
% 
% 
% for k = 1:2
%     clear X1 X2 Y1 Y2 X0 Y0 D1X D2X D1Y D2Y;
%     if k == 1
%         TF = krug_data.fnd_gt(:,kre_ix);
%     elseif k == 2
%         TF = ~krug_data.fnd_gt(:,kre_ix);
%     end
%     
%     subplot('Position', [subplotposx(k) subplotposy(1) subplotwidth subplotheight]);
%     hold on, box on;
%     set(gca, 'FontName', 'Helvetica', 'FontSize', tickfontsize, 'Layer', 'top');
% 
%     X0 = krug_data.fnd_fit(TF,1:2);
%     Y0 = krug_data.fnd_fit(TF,3:4);
%     
%     X1 = krug_data.pop_fit(TF,1:2:7);
%     D1X(:,1) = nanmean(X1,2) - X0(:,1);
%     D1X(:,2) = sqrt( nanstd(X1,0,2).^2./sum(~isnan(X1),2)  + X0(:,2).^2 );
%     
%     Y1 = krug_data.pop_fit(TF, 2:2:8);
%     D1Y(:,1) = nanmean(Y1,2) - Y0(:,1);
%     D1Y(:,2) = sqrt( nanstd(Y1,0,2).^2./sum(~isnan(Y1),2)  + Y0(:,2).^2 );
%     
%     X2 = krug_data.pop_fit(TF,9:2:15);
%     D2X(:,1) = nanmean(X2,2) - X0(:,1);
%     D2X(:,2) = sqrt( nanstd(X2,0,2).^2./sum(~isnan(X2),2)  + X0(:,2).^2 );
%     
%     Y2 = krug_data.pop_fit(TF, 10:2:16);
%     D2Y(:,1) = nanmean(Y2,2) - Y0(:,1);
%     D2Y(:,2) = sqrt( nanstd(Y2,0,2).^2./sum(~isnan(Y2),2)  + Y0(:,2).^2 );
%     
%     %     plot( [X0(:,1) X1(:,1)]', [Y0(:,1) Y1(:,2)]', '-', 'Color', cc(4,:));
%     %     plot( [X0(:,1) X2(:,1)]', [Y0(:,1) Y2(:,2)]', '-', 'Color', cc(5,:));
%     
%     %     plot( [zeros(nnz(TF),1) D1X(:,1)]', [zeros(nnz(TF),1) D1Y(:,1)]', '-', 'Color', cc(4,:));
%     %     plot( [zeros(nnz(TF),1) D2X(:,1)]', [zeros(nnz(TF),1) D2Y(:,1)]', '-', 'Color', cc(5,:));
%     
%     for ix = randperm(nnz(TF))
%         plot( [0 D1X(ix,1)], [0 D1Y(ix,1)], '-', 'Color', cc(4,:));
%         plot( [0 D2X(ix,1)], [0 D2Y(ix,1)], '-', 'Color', cc(5,:));
%     end
%     
%     plot( [0 nanmean(D1X(:,1))], [0 nanmean(D1Y(:,1))], '-', 'Color', cc(2,:), 'LineWidth', 4);
%     plot( [0 nanmean(D2X(:,1))], [0 nanmean(D2Y(:,1))], '-', 'Color', cc(3,:), 'LineWidth', 4);
%     
%     set(gca, 'XLim', [-0.01, 0.21], 'XTick', -0.05:0.05:0.20, 'XGrid', 'on');
%     set(gca, 'YLim', [-0.11, 0.36], 'YTick', -0.10:0.05:0.35, 'YGrid', 'on');
%     
%     if k == 1
% %        title('RM', 'FontName', 'Arial', 'FontSize', 14);
%         text(-0.005, 0.355, 'RM', 'FontName', 'Arial', 'FontSize', 14, 'VerticalAlignment', 'top');
%         ylabel('Fitness increment at 37°C, \DeltaY', 'FontName', 'Arial', 'FontSize', labelfontsize);
%     elseif k == 2
%         set(gca, 'YTickLabels', {});
%         title('BY', 'FontName', 'Arial', 'FontSize', 14);
%         xlabel('Fitness increment at 30°C, \DeltaX', 'FontName', 'Arial', 'FontSize', labelfontsize,...
%             'Position', [-0.02, -0.15]);
%         text(-0.005, 0.355, 'BY', 'FontName', 'Arial', 'FontSize', 14, 'VerticalAlignment', 'top');
%     end    
% end
% clear X1 X2 Y1 Y2 X0 Y0 D1X D2X D1Y D2Y;
% clear figurewidth figureheight tickfontsize labelfontsize cc;
% clear x0 y0 dx dy subplotwidth subplotheight subplotposx subplotposy;
% clear hplot k TF attribs X0 Y0 b1 bint r rint stats r2 b2 h Y1pred Y1obs X2pred X2obs;   





%% OBSOLETE

%% Are pleiotropic slopes different if we stratify by other loci?
%% NOT FINISHED AND NOT USING

% N = 10000;
% slopes_perm.dpca = nan( N , 3 );
% slopes_perm.dlin = nan( N , 3 );
% 
% fprintf('Generating the null distribution of slope differences\n');
% 
% for lix = 1:N
%     
%     TF = (rand(krug_data.n_pop,1) > 0.5 );
%     if nnz(TF) == 0 || nnz(~TF) == 0
%         continue;
%     end
%     
%     RMn = nnz( TF );
%     BYn = nnz( ~TF );
%         
%     RM.X = krug_data.fnd_fit(TF,[1 2]) * 100;
%     RM.Y = krug_data.fnd_fit(TF,[3 4]) * 100;
%     
%     BY.X = krug_data.fnd_fit(~TF,[1 2]) * 100;
%     BY.Y = krug_data.fnd_fit(~TF,[3 4]) * 100;
%     
%     RM.m(1,1) = nanmean( RM.X(:,1) );
%     RM.m(1,2) = nanstd( RM.X(:,1) )/sqrt( sum(~isnan( RM.X(:,1) )) );
%     RM.m(2,1) = nanmean( RM.Y(:,1) );
%     RM.m(2,2) = nanstd( RM.Y(:,1) )/sqrt( sum(~isnan( RM.Y(:,1) )) );
%     
%     BY.m(1,1) = nanmean( BY.X(:,1) );
%     BY.m(1,2) = nanstd( BY.X(:,1) )/sqrt( sum(~isnan(BY.X(:,1))) );
%     BY.m(2,1) = nanmean( BY.Y(:,1) );
%     BY.m(2,2) = nanstd( BY.Y(:,1) )/sqrt( sum(~isnan(BY.Y(:,1))) );
%     
%     
%     % PCA
%     coeffs = pca( [RM.X(:,1), RM.Y(:,1)]);
%     a = coeffs(2,1)/coeffs(1,1);
%     
%     coeffs = pca( [BY.X(:,1), BY.Y(:,1)]);
%     b = coeffs(2,1)/coeffs(1,1);
%     slopes_perm.dpca(lix,1) = abs( b - a );
% 
%     
%     
%     % Linear regression
%     blin = regress( RM.Y(:,1) , [ones(RMn,1), RM.X(:,1)]);
%     a = blin(2);
%     blin = regress( BY.Y(:,1) , [ones(BYn,1), BY.X(:,1)]);
%     b = blin(2);
%     slopes_perm.dlin(lix,1) = abs( b - a );
%         
%         
%     RM.X1 = krug_data.pop_fit(TF,1:2:7) * 100;
%     RM.Y1 = krug_data.pop_fit(TF,2:2:8) * 100;
%     RM.X2 = krug_data.pop_fit(TF,9:2:15) * 100;
%     RM.Y2 = krug_data.pop_fit(TF,10:2:16) * 100;
%     
%     BY.X1 = krug_data.pop_fit(~TF,1:2:7) * 100;
%     BY.Y1 = krug_data.pop_fit(~TF,2:2:8) * 100;
%     BY.X2 = krug_data.pop_fit(~TF,9:2:15) * 100;
%     BY.Y2 = krug_data.pop_fit(~TF,10:2:16) * 100;
%     
%     RM.m1(1,1) = nanmean( RM.X1(:) );
%     RM.m1(1,2) = nanstd( RM.X1(:) )/sqrt( sum(~isnan( RM.X1(:) )) );
%     RM.m1(2,1) = nanmean( RM.Y1(:) );
%     RM.m1(2,2) = nanstd( RM.Y1(:) )/sqrt( sum(~isnan( RM.Y1(:) )) );
%     
%     RM.m2(1,1) = nanmean( RM.X2(:) );
%     RM.m2(1,2) = nanstd( RM.X2(:) )/sqrt( sum(~isnan( RM.X2(:) )) );
%     RM.m2(2,1) = nanmean( RM.Y2(:) );
%     RM.m2(2,2) = nanstd( RM.Y2(:) )/sqrt( sum(~isnan( RM.Y2(:) )) );
%     
%     BY.m1(1,1) = nanmean( BY.X1(:) );
%     BY.m1(1,2) = nanstd( BY.X1(:) )/sqrt( sum(~isnan( BY.X1(:) )) );
%     BY.m1(2,1) = nanmean( BY.Y1(:) );
%     BY.m1(2,2) = nanstd( BY.Y1(:) )/sqrt( sum(~isnan( BY.Y1(:) )) );
%     
%     BY.m2(1,1) = nanmean( BY.X2(:) );
%     BY.m2(1,2) = nanstd( BY.X2(:) )/sqrt( sum(~isnan( BY.X2(:) )) );
%     BY.m2(2,1) = nanmean( BY.Y2(:) );
%     BY.m2(2,2) = nanstd( BY.Y2(:) )/sqrt( sum(~isnan( BY.Y2(:) )) );
%     
%     % PCA:
%     coeffs = pca( [nanmean(RM.X1,2), nanmean(RM.Y1,2)]);
%     a = coeffs(2,1)/coeffs(1,1);
%     
%     coeffs = pca( [nanmean(BY.X1,2), nanmean(BY.Y1,2)]);
%     b = coeffs(2,1)/coeffs(1,1);
%     slopes_perm.dpca(lix,2) = abs( b - a );
%     
%     coeffs = pca( [nanmean(RM.X2,2), nanmean(RM.Y2,2)]);
%     a = coeffs(2,1)/coeffs(1,1);
%     
%     coeffs = pca( [nanmean(BY.X2,2), nanmean(BY.Y2,2)]);
%     b = coeffs(2,1)/coeffs(1,1);
%     slopes_perm.dpca(lix,3) = abs( b - a );
%     
%     % Linear regression
%     blin = regress( nanmean(RM.Y1,2) , [ones(RMn,1), nanmean(RM.X1,2)]);
%     a = blin(2);
%     blin = regress( nanmean(BY.Y1,2) , [ones(BYn,1), nanmean(BY.X1,2)]);
%     b = blin(2);
%     slopes_perm.dlin(lix,2) = abs( b - a );
% 
%     
%     blin = regress( nanmean(RM.Y2,2) , [ones(RMn,1), nanmean(RM.X2,2)]);
%     a = blin(2);
%     blin = regress( nanmean(BY.Y2,2) , [ones(BYn,1), nanmean(BY.X2,2)]);
%     b = blin(2);
%     slopes_perm.dlin(lix,3) = abs( b - a );
% end
% 
% clear TF RMn BYn RM BY coeffs blin lix;




% fprintf('Pleiotropic slope analysis for all loci\n');
% 
% slopes.bpca = nan( krug_data.n_snp , 6 );
% slopes.blin = nan( krug_data.n_snp , 6 );
% 
% slopes.dpca = nan( krug_data.n_snp , 3 );
% slopes.dlin = nan( krug_data.n_snp , 3 );
% 
% slopes.ppca = nan( krug_data.n_snp , 3 );
% slopes.plin = nan( krug_data.n_snp , 3 );
% 
% for lix = 1:krug_data.n_snp
%     
%     TF = krug_data.fnd_gt( : , lix);
%     if nnz(TF) == 0 || nnz(~TF) == 0
%         continue;
%     end
%     RMn = nnz( TF );
%     BYn = nnz( ~TF );
%         
%     RM.X = krug_data.fnd_fit(TF,[1 2]) * 100;
%     RM.Y = krug_data.fnd_fit(TF,[3 4]) * 100;
%     
%     BY.X = krug_data.fnd_fit(~TF,[1 2]) * 100;
%     BY.Y = krug_data.fnd_fit(~TF,[3 4]) * 100;
%     
%     RM.m(1,1) = nanmean( RM.X(:,1) );
%     RM.m(1,2) = nanstd( RM.X(:,1) )/sqrt( sum(~isnan( RM.X(:,1) )) );
%     RM.m(2,1) = nanmean( RM.Y(:,1) );
%     RM.m(2,2) = nanstd( RM.Y(:,1) )/sqrt( sum(~isnan( RM.Y(:,1) )) );
%     
%     BY.m(1,1) = nanmean( BY.X(:,1) );
%     BY.m(1,2) = nanstd( BY.X(:,1) )/sqrt( sum(~isnan(BY.X(:,1))) );
%     BY.m(2,1) = nanmean( BY.Y(:,1) );
%     BY.m(2,2) = nanstd( BY.Y(:,1) )/sqrt( sum(~isnan(BY.Y(:,1))) );
%     
%     
%     % PCA
%     coeffs = pca( [RM.X(:,1), RM.Y(:,1)]);
%     slopes.bpca(lix,1) = coeffs(2,1)/coeffs(1,1);
%     
%     coeffs = pca( [BY.X(:,1), BY.Y(:,1)]);
%     slopes.bpca(lix,2) = coeffs(2,1)/coeffs(1,1);
%     slopes.dpca(lix,1) = abs(slopes.bpca(lix,2) - slopes.bpca(lix,1));
%     slopes.ppca(lix,1) = nnz( slopes_perm.dpca(:,1) >= slopes.dpca(lix,1)) / N;
%     
%     % Linear regression
%     blin = regress( RM.Y(:,1) , [ones(RMn,1), RM.X(:,1)]);
%     slopes.blin(lix,1) = blin(2);
%     blin = regress( BY.Y(:,1) , [ones(BYn,1), BY.X(:,1)]);
%     slopes.blin(lix,2) = blin(2);
%     slopes.dlin(lix,1) = abs(slopes.blin(lix,2) - slopes.blin(lix,1));
%     slopes.plin(lix,1) = nnz( slopes_perm.dlin(:,1) >= slopes.dlin(lix,1)) / N;
%         
%         
%     RM.X1 = krug_data.pop_fit(TF,1:2:7) * 100;
%     RM.Y1 = krug_data.pop_fit(TF,2:2:8) * 100;
%     RM.X2 = krug_data.pop_fit(TF,9:2:15) * 100;
%     RM.Y2 = krug_data.pop_fit(TF,10:2:16) * 100;
%     
%     BY.X1 = krug_data.pop_fit(~TF,1:2:7) * 100;
%     BY.Y1 = krug_data.pop_fit(~TF,2:2:8) * 100;
%     BY.X2 = krug_data.pop_fit(~TF,9:2:15) * 100;
%     BY.Y2 = krug_data.pop_fit(~TF,10:2:16) * 100;
%     
%     RM.m1(1,1) = nanmean( RM.X1(:) );
%     RM.m1(1,2) = nanstd( RM.X1(:) )/sqrt( sum(~isnan( RM.X1(:) )) );
%     RM.m1(2,1) = nanmean( RM.Y1(:) );
%     RM.m1(2,2) = nanstd( RM.Y1(:) )/sqrt( sum(~isnan( RM.Y1(:) )) );
%     
%     RM.m2(1,1) = nanmean( RM.X2(:) );
%     RM.m2(1,2) = nanstd( RM.X2(:) )/sqrt( sum(~isnan( RM.X2(:) )) );
%     RM.m2(2,1) = nanmean( RM.Y2(:) );
%     RM.m2(2,2) = nanstd( RM.Y2(:) )/sqrt( sum(~isnan( RM.Y2(:) )) );
%     
%     BY.m1(1,1) = nanmean( BY.X1(:) );
%     BY.m1(1,2) = nanstd( BY.X1(:) )/sqrt( sum(~isnan( BY.X1(:) )) );
%     BY.m1(2,1) = nanmean( BY.Y1(:) );
%     BY.m1(2,2) = nanstd( BY.Y1(:) )/sqrt( sum(~isnan( BY.Y1(:) )) );
%     
%     BY.m2(1,1) = nanmean( BY.X2(:) );
%     BY.m2(1,2) = nanstd( BY.X2(:) )/sqrt( sum(~isnan( BY.X2(:) )) );
%     BY.m2(2,1) = nanmean( BY.Y2(:) );
%     BY.m2(2,2) = nanstd( BY.Y2(:) )/sqrt( sum(~isnan( BY.Y2(:) )) );
%     
%     % PCA:
%     coeffs = pca( [nanmean(RM.X1,2), nanmean(RM.Y1,2)]);
%     slopes.bpca(lix,3) = coeffs(2,1)/coeffs(1,1);   
%     
%     coeffs = pca( [nanmean(BY.X1,2), nanmean(BY.Y1,2)]);
%     slopes.bpca(lix,4) = coeffs(2,1)/coeffs(1,1);    
%     
%     coeffs = pca( [nanmean(RM.X2,2), nanmean(RM.Y2,2)]);
%     slopes.bpca(lix,5) = coeffs(2,1)/coeffs(1,1);
%     
%     coeffs = pca( [nanmean(BY.X2,2), nanmean(BY.Y2,2)]);
%     slopes.bpca(lix,6) = coeffs(2,1)/coeffs(1,1);
% 
%     slopes.dpca(lix,2) = abs(slopes.bpca(lix,4) - slopes.bpca(lix,3));
%     slopes.dpca(lix,3) = abs(slopes.bpca(lix,6) - slopes.bpca(lix,5));
%     
%     slopes.ppca(lix,2) = nnz( slopes_perm.dpca(:,2) >= slopes.dpca(lix,2)) / N;
%     slopes.ppca(lix,3) = nnz( slopes_perm.dpca(:,3) >= slopes.dpca(lix,3)) / N;
%     
%     % Linear regression
%     blin = regress( nanmean(RM.Y1,2) , [ones(RMn,1), nanmean(RM.X1,2)]);
%     slopes.blin(lix,3) = blin(2);
%     blin = regress( nanmean(BY.Y1,2) , [ones(BYn,1), nanmean(BY.X1,2)]);
%     slopes.blin(lix,4) = blin(2);
%     
%     blin = regress( nanmean(RM.Y2,2) , [ones(RMn,1), nanmean(RM.X2,2)]);
%     slopes.blin(lix,5) = blin(2);
%     blin = regress( nanmean(BY.Y2,2) , [ones(BYn,1), nanmean(BY.X2,2)]);
%     slopes.blin(lix,6) = blin(2);
%     
%     slopes.dlin(lix,2) = abs(slopes.blin(lix,4) - slopes.blin(lix,3));
%     slopes.dlin(lix,3) = abs(slopes.blin(lix,6) - slopes.blin(lix,5));
%     
%     slopes.plin(lix,2) = nnz( slopes_perm.dlin(:,2) >= slopes.dlin(lix,2)) / N;
%     slopes.plin(lix,3) = nnz( slopes_perm.dlin(:,3) >= slopes.dlin(lix,3)) / N;
% end
% 
% clear a b d TF RMn BYn RM BY coeffs blin lix;





%% Calculating the pleiotropic effect of BY allele (PCA)

% Nperm = 1000; % Number of permutations to estimate the confidence intervals on excess fitness effect of the KRE33-BY allele
% 
% TF = krug_data.fnd_gt( : , kre_ix);
% RMn = nnz( TF );
% BYn = nnz( ~TF );
% 
% fprintf('Pleiotropic slope analysis\n');
% fprintf('Nperm = %d\n', Nperm);
% fprintf('-----\nFounders\n-----\n');
% 
% RM.X = krug_data.fnd_fit(TF,[1 2]) * 100;
% RM.Y = krug_data.fnd_fit(TF,[3 4]) * 100;
% 
% BY.X = krug_data.fnd_fit(~TF,[1 2]) * 100;
% BY.Y = krug_data.fnd_fit(~TF,[3 4]) * 100;
% 
% RM.m(1,1) = nanmean( RM.X(:,1) );
% RM.m(1,2) = nanstd( RM.X(:,1) )/sqrt( sum(~isnan( RM.X(:,1) )) );
% RM.m(2,1) = nanmean( RM.Y(:,1) );
% RM.m(2,2) = nanstd( RM.Y(:,1) )/sqrt( sum(~isnan( RM.Y(:,1) )) );
% 
% BY.m(1,1) = nanmean( BY.X(:,1) );
% BY.m(1,2) = nanstd( BY.X(:,1) )/sqrt( sum(~isnan(BY.X(:,1))) );
% BY.m(2,1) = nanmean( BY.Y(:,1) );
% BY.m(2,2) = nanstd( BY.Y(:,1) )/sqrt( sum(~isnan(BY.Y(:,1))) );
% 
% 
% % PCA
% 
% X = [RM.X(:,1) - RM.m(1,1); BY.X(:,1) - BY.m(1,1)];
% Y = [RM.Y(:,1) - RM.m(2,1); BY.Y(:,1) - BY.m(2,1)];
% 
% coeffs = pca( [X, Y]);
% b_fnd(1,2) = coeffs(2,1)/coeffs(1,1);
% b_fnd(1,1) = RM.m(2,1) - RM.m(1,1) * b_fnd(1,2);
% b_fnd(1,3) = BY.m(2,1) - BY.m(1,1) * b_fnd(1,2);
% 
% fnd_df = BY.m(2,1) - ( b_fnd(1)+b_fnd(2)*BY.m(1,1) );
% 
% perms.fnd_b = nan(Nperm, 4);
% perms.fnd_df = nan(Nperm, 1);
% 
% for iperm = 1:Nperm
%     RMXperm = RM.X(:,1) + RM.X(:,2) .* randn( RMn , 1 );
%     RMYperm = RM.Y(:,1) + RM.Y(:,2) .* randn( RMn , 1 );
%     
%     BYXperm = BY.X(:,1) + BY.X(:,2) .* randn( BYn , 1 );
%     BYYperm = BY.Y(:,1) + BY.Y(:,2) .* randn( BYn , 1 );
%         
%     coeffs = pca( [RMXperm, RMYperm]);
%     perms.fnd_b(iperm,2) = coeffs(2,1)/coeffs(1,1);
%     perms.fnd_b(iperm,1) = nanmean( RMYperm ) - nanmean( RMXperm ) * perms.fnd_b(iperm,2);
% 
%     coeffs = pca( [BYXperm, BYYperm]);
%     perms.fnd_b(iperm,4) = coeffs(2,1)/coeffs(1,1);
%     perms.fnd_b(iperm,3) = nanmean( BYYperm ) - nanmean( BYXperm ) * perms.fnd_b(iperm,4);
% 
%     perms.fnd_df(iperm) = nanmean( BYYperm ) - ( perms.fnd_b(iperm,1) + perms.fnd_b(iperm,2) * nanmean( BYXperm ) );
%     
%     % Linear regression:
%     perms.fnd_blin(iperm,[1 2]) = regress( RMYperm, [ones(RMn,1), RMXperm]);
%     perms.fnd_blin(iperm,[3 4]) = regress( BYYperm, [ones(BYn,1), BYXperm]);
% end
% clear coeffs iperm RMXperm RMYperm BYXperm BYYperm ;
% 
% 
% fprintf('Mean fitness of KRE33-RM founders\n\tat 30°C: %.2f +- %.2f %%\n\tat 37°C: %.2f +- %.2f %%\n',...
%     RM.m(1,:), RM.m(2,:));
% 
% fprintf('Mean fitness of KRE33-BY founders\n\tat 30°C: %.2f +- %.2f %%\n\tat 37°C: %.2f +- %.2f %%\n',...
%     BY.m(1,:), BY.m(2,:));
% 
% fprintf('Mean fitness advantage of KRE33-BY founders\n\tat 30°C: %.2f +- %.2f %%\n\tat 37°C: %.2f +- %.2f %%\n',...
%     BY.m(1,1) - RM.m(1,1), sqrt( BY.m(1,2)^2 + RM.m(1,2)^2 ),...
%     BY.m(2,1) - RM.m(2,1), sqrt( BY.m(2,2)^2 + RM.m(2,2)^2 ) );
% 
% fprintf('PCA Pleiotropy slope\n\tfor KRE33-RM founders: %.2f (CI_{95%%} = [%.2f, %.2f]\n',...
%     b_fnd(2), quantile(perms.fnd_b(:,2), [0.025,0.975]));
% 
% % fprintf('Expected fitness advantage of BY founders\n\tat 37°C: %.2f %%\n',...
% %     ( b_fnd(1)+b_fnd(2)*BY.m(1,1)) - RM.m(2,1) );
% % fprintf('Excess fitness advantage of BY founders\n\tat 37°C: %.2f %% (CI_{95%%} = [%.2f, %.2f]\n',...
% %     fnd_df, quantile(perms.fnd_df, [0.025,0.975]));
% 
% 
% fprintf('-----\nEvolved\n-----\n');
% 
% RM.X1 = krug_data.pop_fit(TF,1:2:7) * 100;
% RM.Y1 = krug_data.pop_fit(TF,2:2:8) * 100;
% RM.X2 = krug_data.pop_fit(TF,9:2:15) * 100;
% RM.Y2 = krug_data.pop_fit(TF,10:2:16) * 100;
% 
% BY.X1 = krug_data.pop_fit(~TF,1:2:7) * 100;
% BY.Y1 = krug_data.pop_fit(~TF,2:2:8) * 100;
% BY.X2 = krug_data.pop_fit(~TF,9:2:15) * 100;
% BY.Y2 = krug_data.pop_fit(~TF,10:2:16) * 100;
% 
% RM.m1(1,1) = nanmean( RM.X1(:) );
% RM.m1(1,2) = nanstd( RM.X1(:) )/sqrt( sum(~isnan( RM.X1(:) )) );
% RM.m1(2,1) = nanmean( RM.Y1(:) );
% RM.m1(2,2) = nanstd( RM.Y1(:) )/sqrt( sum(~isnan( RM.Y1(:) )) );
% 
% RM.m2(1,1) = nanmean( RM.X2(:) );
% RM.m2(1,2) = nanstd( RM.X2(:) )/sqrt( sum(~isnan( RM.X2(:) )) );
% RM.m2(2,1) = nanmean( RM.Y2(:) );
% RM.m2(2,2) = nanstd( RM.Y2(:) )/sqrt( sum(~isnan( RM.Y2(:) )) );
% 
% BY.m1(1,1) = nanmean( BY.X1(:) );
% BY.m1(1,2) = nanstd( BY.X1(:) )/sqrt( sum(~isnan( BY.X1(:) )) );
% BY.m1(2,1) = nanmean( BY.Y1(:) );
% BY.m1(2,2) = nanstd( BY.Y1(:) )/sqrt( sum(~isnan( BY.Y1(:) )) );
% 
% BY.m2(1,1) = nanmean( BY.X2(:) );
% BY.m2(1,2) = nanstd( BY.X2(:) )/sqrt( sum(~isnan( BY.X2(:) )) );
% BY.m2(2,1) = nanmean( BY.Y2(:) );
% BY.m2(2,2) = nanstd( BY.Y2(:) )/sqrt( sum(~isnan( BY.Y2(:) )) );
% 
% % PCA:
% % [b1,bint,r,rint,stats] = regress( nanmean( RM.Y1,2), [ones(size(RM.X1,1),1), nanmean(RM.X1,2)]);
% coeffs = pca( [nanmean(RM.X1,2), nanmean(RM.Y1,2)]);
% b1(1,2) = coeffs(2,1)/coeffs(1,1);
% b1(1,1) = RM.m1(2,1) - RM.m1(1,1) * b1(1,2);
% df1 = BY.m1(2,1) - ( b1(1,1) + b1(1,2) * BY.m1(1,1) );
% 
% coeffs = pca( [nanmean(BY.X1,2), nanmean(BY.Y1,2)]);
% b1(1,4) = coeffs(2,1)/coeffs(1,1);
% b1(1,3) = BY.m1(2,1) - BY.m1(1,1) * b1(1,4);
% 
% 
% % [b2,bint,r,rint,stats] = regress( nanmean( RM.Y2,2), [ones(size(RM.X2,1),1), nanmean(RM.X2,2)]);
% coeffs = pca( [nanmean(RM.X2,2), nanmean(RM.Y2,2)]);
% b2(1,2) = coeffs(2,1)/coeffs(1,1);
% b2(1,1) = RM.m2(2,1) - RM.m2(1,1) * b2(1,2);
% df2 = BY.m2(2,1) - ( b2(1,1) + b2(1,2) * BY.m2(1,1) );
% 
% coeffs = pca( [nanmean(BY.X2,2), nanmean(BY.Y2,2)]);
% b2(1,4) = coeffs(2,1)/coeffs(1,1);
% b2(1,3) = BY.m2(2,1) - BY.m2(1,1) * b2(1,4);
% 
% % Linear regression
% blin1(1,[1 2]) = regress( nanmean(RM.Y1,2) , [ones(RMn,1), nanmean(RM.X1,2)]);
% blin1(1,[3 4]) = regress( nanmean(BY.Y1,2) , [ones(BYn,1), nanmean(BY.X1,2)]);
% 
% blin2(1,[1 2]) = regress( nanmean(RM.Y2,2) , [ones(RMn,1), nanmean(RM.X2,2)]);
% blin2(1,[3 4]) = regress( nanmean(BY.Y2,2) , [ones(BYn,1), nanmean(BY.X2,2)]);
% 
% 
% perms.b1 = nan(Nperm, 4);
% perms.blin1 = nan(Nperm, 4);
% perms.df1 = nan(Nperm, 1);
% 
% perms.b2 = nan(Nperm, 4);
% perms.blin2 = nan(Nperm, 4);
% perms.df2 = nan(Nperm, 1);
% 
% sigma1 = nanmean( nanstd(krug_data.errest(:,:,1),1,2) )*100;
% sigma2 = nanmean( nanstd(krug_data.errest(:,:,2),1,2) )*100;
% 
% for iperm = 1:Nperm
%     data_perm = nan(krug_data.n_pop, 16);
%     data_perm(:,1:2:15) = krug_data.pop_fit(:,1:2:15) * 100 + sigma1 * randn( krug_data.n_pop, 8 );
%     data_perm(:,2:2:16) = krug_data.pop_fit(:,2:2:16) * 100 + sigma2 * randn( krug_data.n_pop, 8 );
%     
%     RMperm.X1 = data_perm(TF,1:2:7);
%     RMperm.Y1 = data_perm(TF,2:2:8);    
%     BYperm.X1 = data_perm(~TF,1:2:7);
%     BYperm.Y1 = data_perm(~TF,2:2:8);
% 
%     RMperm.X2 = data_perm(TF,9:2:15);
%     RMperm.Y2 = data_perm(TF,10:2:16);
%     BYperm.X2 = data_perm(~TF,9:2:15);
%     BYperm.Y2 = data_perm(~TF,10:2:16);
%     
%     % PCA:
%     coeffs = pca( [nanmean(RMperm.X1,2), nanmean(RMperm.Y1,2)]);
%     perms.b1(iperm,2) = coeffs(2,1)/coeffs(1,1);
%     perms.b1(iperm,1) = nanmean(RMperm.Y1(:)) - nanmean(RMperm.X1(:)) * perms.b1(iperm,2);
%     perms.df1(iperm) = nanmean( BYperm.Y1(:) ) - ( perms.b1(iperm,1) + perms.b1(iperm,2) * nanmean( BYperm.X1(:) ) );
% 
%     coeffs = pca( [nanmean(BYperm.X1,2), nanmean(BYperm.Y1,2)]);
%     perms.b1(iperm,4) = coeffs(2,1)/coeffs(1,1);
%     perms.b1(iperm,3) = nanmean(BYperm.Y1(:)) - nanmean(BYperm.X1(:)) * perms.b1(iperm,4);
%     
%     coeffs = pca( [nanmean(RMperm.X2,2), nanmean(RMperm.Y2,2)]);
%     perms.b2(iperm,2) = coeffs(2,1)/coeffs(1,1);
%     perms.b2(iperm,1) = nanmean(RMperm.Y2(:)) - nanmean(RMperm.X2(:)) * perms.b2(iperm,2);
%     perms.df2(iperm) = nanmean( BYperm.Y2(:) ) - ( perms.b2(iperm,1) + perms.b2(iperm,2) * nanmean( BYperm.X2(:) ) );
% 
%     coeffs = pca( [nanmean(BYperm.X2,2), nanmean(BYperm.Y2,2)]);
%     perms.b2(iperm,4) = coeffs(2,1)/coeffs(1,1);
%     perms.b2(iperm,3) = nanmean(BYperm.Y2(:)) - nanmean(BYperm.X2(:)) * perms.b2(iperm,4);
%     
%     % Linear regression
%     perms.blin1(iperm,[1 2]) = regress( nanmean(RMperm.Y1,2) , [ones(RMn,1), nanmean(RMperm.X1,2)]);
%     perms.blin1(iperm,[3 4]) = regress( nanmean(BYperm.Y1,2) , [ones(BYn,1), nanmean(BYperm.X1,2)]);
% 
%     perms.blin2(iperm,[1 2]) = regress( nanmean(RMperm.Y2,2) , [ones(RMn,1), nanmean(RMperm.X2,2)]);
%     perms.blin2(iperm,[3 4]) = regress( nanmean(BYperm.Y2,2) , [ones(BYn,1), nanmean(BYperm.X2,2)]);
% 
% end
% clear data_perm coeffs iperm RMperm.X1 RMperm.Y1 BYperm.X1 BYperm.Y1 RMperm.X2 RMperm.Y2 BYperm.X2 BYperm.Y2;
% 
% fprintf('Evolved at 30°C:\n');
% fprintf('Mean fitness of RM pops\n\tat 30°C: %.2f +- %.2f %%\n\tat 37°C: %.2f +- %.2f %%\n',...
%     RM.m1(1,:), RM.m1(2,:));
% 
% fprintf('Mean fitness gain achieved by RM pops\n\tat 30°C: %.2f +- %.2f %%\n',...
%     RM.m1(1,1)-RM.m(1,1), sqrt(RM.m1(1,2)^2 + RM.m(1,2)^2) );
% fprintf('\tat 37°C: %.2f +- %.2f %%\n',...
%     RM.m1(2,1)-RM.m(2,1), sqrt(RM.m1(2,2)^2 + RM.m(2,2)^2) );
% 
% fprintf('Mean fitness of BY pops\n\tat 30°C: %.2f +- %.2f %%\n\tat 37°C: %.2f +- %.2f %%\n',...
%     BY.m1(1,:), BY.m1(2,:));
% 
% fprintf('Mean fitness gain achieved by BY pops\n\tat 30°C: %.2f +- %.2f %%\n',...
%     BY.m1(1,1)-BY.m(1,1), sqrt(BY.m1(1,2)^2 + BY.m(1,2)^2) );
% fprintf('\tat 37°C: %.2f +- %.2f %%\n',...
%     BY.m1(2,1)-BY.m(2,1), sqrt(BY.m1(2,2)^2 + BY.m(2,2)^2) );
% 
% fprintf('PCA Pleiotropy slope\n\tfor KRE33-RM founders: %.2f (CI_{95%%} = [%.2f, %.2f]\n',...
%     b1(2), quantile(perms.b1(:,2), [0.025,0.975]));
% fprintf('\tfor KRE33-BY founders: %.2f (CI_{95%%} = [%.2f, %.2f]\n',...
%     b1(4), quantile(perms.b1(:,4), [0.025,0.975]));
% 
% fprintf('LR Pleiotropy slope\n\tfor KRE33-RM founders: %.2f (CI_{95%%} = [%.2f, %.2f]\n',...
%     blin1(2), quantile(perms.blin1(:,2), [0.025,0.975]));
% fprintf('\tfor KRE33-BY founders: %.2f (CI_{95%%} = [%.2f, %.2f]\n',...
%     blin1(4), quantile(perms.blin1(:,4), [0.025,0.975]));
% 
% 
% % fprintf('Expected fitness advantage of BY pops\n\tat 37°C: %.2f %%\n',...
% %     ( b1(1)+b1(2)*BY.m1(1,1)) - RM.m1(2,1) );
% % 
% % fprintf('Excess fitness advantage of BY founders\n\tat 37°C: %.2f %% (CI_{95%%} = [%.2f, %.2f]\n',...
% %     df1, quantile(perms.df1, [0.025,0.975]));
% 
% 
% fprintf('\nEvolved at 37°C:\n');
% fprintf('Mean fitness of RM pops\n\tat 30°C: %.2f +- %.2f %%\n\tat 37°C: %.2f +- %.2f %%\n',...
%     RM.m2(1,:), RM.m2(2,:));
% fprintf('Mean fitness gain achieved by RM pops\n\tat 30°C: %.2f +- %.2f %%\n',...
%     RM.m2(1,1)-RM.m(1,1), sqrt(RM.m2(1,2)^2 + RM.m(1,2)^2) );
% fprintf('\tat 37°C: %.2f +- %.2f %%\n',...
%     RM.m2(2,1)-RM.m(2,1), sqrt(RM.m2(2,2)^2 + RM.m(2,2)^2) );
% 
% fprintf('Mean fitness of BY pops\n\tat 30°C: %.2f +- %.2f %%\n\tat 37°C: %.2f +- %.2f %%\n',...
%     BY.m2(1,:), BY.m2(2,:));
% 
% fprintf('Mean fitness gain achieved by BY pops\n\tat 30°C: %.2f +- %.2f %%\n',...
%     BY.m2(1,1)-BY.m(1,1), sqrt(BY.m2(1,2)^2 + BY.m(1,2)^2) );
% fprintf('\tat 37°C: %.2f +- %.2f %%\n',...
%     BY.m2(2,1)-BY.m(2,1), sqrt(BY.m2(2,2)^2 + BY.m(2,2)^2) );
% 
% fprintf('Mean fitness advantage of BY pops\n\tat 30°C: %.2f +- %.2f %%\n\tat 37°C: %.2f +- %.2f %%\n',...
%     BY.m2(1,1) - RM.m2(1,1), sqrt( BY.m2(1,2)^2 + RM.m2(1,2)^2 ),...
%     BY.m2(2,1) - RM.m2(2,1), sqrt( BY.m2(2,2)^2 + RM.m2(2,2)^2 ) );
% 
% fprintf('PCA Pleiotropy slope\n\tfor KRE33-RM founders: %.2f (CI_{95%%} = [%.2f, %.2f]\n',...
%     b2(2), quantile(perms.b2(:,2), [0.025,0.975]));
% fprintf('\tfor KRE33-BY founders: %.2f (CI_{95%%} = [%.2f, %.2f]\n',...
%     b2(4), quantile(perms.b2(:,4), [0.025,0.975]));
% 
% fprintf('LR Pleiotropy slope\n\tfor KRE33-RM founders: %.2f (CI_{95%%} = [%.2f, %.2f]\n',...
%     blin2(2), quantile(perms.blin2(:,2), [0.025,0.975]));
% fprintf('\tfor KRE33-BY founders: %.2f (CI_{95%%} = [%.2f, %.2f]\n',...
%     blin2(4), quantile(perms.blin2(:,4), [0.025,0.975]));
% 
% 
% % fprintf('Expected fitness advantage of BY pops\n\tat 37°C: %.2f %%\n',...
% %     ( b2(1)+b2(2)*BY.m2(1,1)) - RM.m2(2,1) );
% % 
% % fprintf('Excess fitness advantage of BY founders\n\tat 37°C: %.2f %% (CI_{95%%} = [%.2f, %.2f]\n',...
% %     df2, quantile(perms.df2, [0.025,0.975]));
% 
% RM.dX1 = RM.X1 - repmat(RM.X(:,1), 1, 4);
% RM.dX2 = RM.X2 - repmat(RM.X(:,1), 1, 4);
% RM.dY1 = RM.Y1 - repmat(RM.Y(:,1), 1, 4);
% RM.dY2 = RM.Y2 - repmat(RM.Y(:,1), 1, 4);
% 
% RM.d1 = [nanmean( RM.dX1,2 ), nanmean( RM.dY1,2 )];
% RM.d2 = [nanmean( RM.dX2,2 ), nanmean( RM.dY2,2 )];
% 
% RM.cos_alpha = sum(RM.d1 .* RM.d2, 2) ./ sqrt(sum(RM.d1.^2,2))./ sqrt(sum(RM.d2.^2,2));
% RM.alpha = acos( RM.cos_alpha )/pi*180;
% 
% BY.dX1 = BY.X1 - repmat(BY.X(:,1), 1, 4);
% BY.dX2 = BY.X2 - repmat(BY.X(:,1), 1, 4);
% BY.dY1 = BY.Y1 - repmat(BY.Y(:,1), 1, 4);
% BY.dY2 = BY.Y2 - repmat(BY.Y(:,1), 1, 4);
% 
% BY.d1 = [nanmean( BY.dX1,2 ), nanmean( BY.dY1,2 )];
% BY.d2 = [nanmean( BY.dX2,2 ), nanmean( BY.dY2,2 )];
% 
% BY.cos_alpha = sum(BY.d1 .* BY.d2, 2) ./ sqrt(sum(BY.d1.^2,2))./ sqrt(sum(BY.d2.^2,2));
% BY.alpha = acos( BY.cos_alpha )/pi*180;
% 
% % RM.alpha1 = (RM.Y1 - repmat(RM.Y(:,1), 1, 4))./(RM.X1 - repmat(RM.X(:,1), 1, 4));
% % BY.alpha1 = (BY.Y1 - repmat(BY.Y(:,1), 1, 4))./(BY.X1 - repmat(BY.X(:,1), 1, 4));
% % 
% % RM.alpha2 = (RM.Y2 - repmat(RM.Y(:,1), 1, 4))./(RM.X2 - repmat(RM.X(:,1), 1, 4));
% % BY.alpha2 = (BY.Y2 - repmat(BY.Y(:,1), 1, 4))./(BY.X2 - repmat(BY.X(:,1), 1, 4));
% 
% clear b_fnd b1 b2 blin blin1 blin2 BY BYn df1 df2 fnd_df Nperm perms RM RMn sigma1 sigma2 TF;





