[Xchn,Zmpn,Xch,Zmp,y,ydot,xcdotdes, t] = getXcZmpHaradaStepFeedBackForward( 0.1, 0.2, 0, [-5, -10], [0, 0]);

% figure(1)
% plot(t, Xchn, t, Zmpn, t, Xch, t, Zmp, t, y, t, ydot, t, xcdotdes)
% legend('Xchn', 'Zmpn', 'Xch', 'Zmp', 'y', 'ydot', 'xcdotdes')

% figure(2)

% Xchn(1)
% error = calculate_error(y, Xch, ydot, xcdotdes);
% plot(t, Xch, t, xcdotdes, t, y, t, ydot, t, error);
% legend('Xch', 'xcdotdes', 'y', 'xcdotdes', 'quadratic error')
% 
% figure(3)
is = 10;
js = 10;

cmap = distinguishable_colors(100);
graph_num = 0;

% Eigenvalues

for i = -10:2:-1
    for j = -10:2:-1
        if i ~= j
            graph_num = graph_num + 1;
            DesEig = [i, j];
            [Xchn,Zmpn,Xch,Zmp,y,ydot,xcdotdes, t] = getXcZmpHaradaStepFeedBackForward( 0.1, 0.2, 0, DesEig, [0, 0]);
            error = calculate_error(y, Xch, ydot, xcdotdes);
            p = plot(t, error, 'Color', cmap(i*j,:));
           
            hold on
            legendInfo{graph_num} = [num2str([i, j])];
        end
    end
end
legend(legendInfo);
hold off

% Initial values

figure(4)
plot_num = 0;
for i = 0.2:0.2:0.8
    for j = 0.2:0.2:0.8
        plot_num = plot_num + 1;
        x0 = [i, j];
        [Xchn,Zmpn,Xch,Zmp,y,ydot,xcdotdes, t] = getXcZmpHaradaStepFeedBackForward( 0.1, 0.2, 0, [-5, -10], x0);
        error = calculate_error(y, Xch, ydot, xcdotdes);
        p = plot(t, error, 'Color', cmap(plot_num,:));
        legendInfo{plot_num} = [num2str([i, j])];
        hold on
    end
end
legend(legendInfo);
xlim([0,1]);
ylim([0,70]);
hold off