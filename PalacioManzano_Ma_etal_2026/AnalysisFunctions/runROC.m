
function AUC = runROC(pTime, data1, data2, sigLvl, Nreps, hAx)

offset = 0;

AUC = zeros(3,numel(pTime));    % AUC and its 95% bootstrap percentiles
for samp = 1:numel(pTime)
    Progress=100*samp/numel(pTime);
    if samp==1
        fprintf(['\t'  'Computing: ' num2str(Progress, '%12.0f') ' %%'])
        TextLength = length(num2str(Progress, '%12.0f'));
    else
    fprintf(repmat('\b',1,TextLength+2));
    fprintf([num2str(Progress, '%12.0f') ' %%']);
    end
    TextLength = length(num2str(Progress, '%12.0f'));

    N1 = size(data1,1);
    N2 = size(data2,1);
    [AUCsamp, Px, Py] = getAUC(data1(:,samp)+offset, data2(:,samp)+offset, N1, N2);
    AUC(1,samp) = AUCsamp;

    AUCreps = zeros(Nreps,1);
    for rep=1:Nreps
        indRand1 = randi(N1,N1,1);
        indRand2 = randi(N2,N2,1);

        [AUCsamp, Px, Py] = getAUC(data1(indRand1,samp)+offset, data2(indRand2,samp)+offset, N1, N2);
        AUCreps(rep,1) = AUCsamp;
    end
    AUCreps = sort(AUCreps);
    AUC(2,samp) = AUCreps(round(Nreps*sigLvl),1);
    AUC(3,samp) = AUCreps(round(Nreps*(1-sigLvl)),1);
end

fprintf('\n')

% Create closed polygon coordinates
x_shade = [pTime, fliplr(pTime)];
y_shade = [AUC(2,:), fliplr(AUC(3,:))];

hold(hAx, 'on')
plot(hAx, pTime, AUC(1,:), 'k', 'LineWidth', 2);
fill(hAx, x_shade, y_shade, [0 0 0], 'EdgeColor', 'none', 'FaceAlpha', 0.2);
line(hAx, [pTime(1) pTime(end)], [0.5 0.5], 'Color', [1 0 0], 'LineStyle', '--', 'LineWidth', 1)
hold(hAx, 'off')

