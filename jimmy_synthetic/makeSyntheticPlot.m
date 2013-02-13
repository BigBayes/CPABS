f = figure()
subplot(2,2,1);
imagesc(Znew);
ylabel('Actor');
title('Features Z');
set(gca, 'YTick', [1 20:20:100]);
set(gca, 'XTick', []);
subplot(2,2,2);
imagesc(Xnew);
set(gca, 'XTick', []);
set(gca, 'YTick', []);
title('Network Y');
subplot(2,2,3);
imagesc(W);
xlabel('Feature');
ylabel('Feature');
title('Feature Weights W');
subplot(2,2,4);
imagesc(Znew');
xlabel('Actor');
title('Features Z^T')
set(gca, 'XTick', [1 20:20:100]);
set(gca, 'YTick', []);
set(gca, 'FontSize', 10);
colormap(1 - gray)
set(f, 'Position', [0 0 350 350]);