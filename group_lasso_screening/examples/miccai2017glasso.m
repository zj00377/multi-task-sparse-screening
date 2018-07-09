ADCUg = thickness(1:160, :);
ADCUg(161:330, :) = thickness(836-170+1:836, :);

ADCUglabel(1:160, 1) = 1;
ADCUglabel(161:330, 1) = -1;

X = zscore(ADCUg);
y = zscore(ADCUglabel);

example_gLasso

ADCUori = thickness(1:199, :);
ADCUori(200:429, :) = thickness(607:836, :);

for i=1:5000
    ADCUdata(:, i) = ADCUori(:, I(i));
end

ADCUdata = ADCUdata';
save('thickBase.txt', 'ADCUdata', '-ascii');

