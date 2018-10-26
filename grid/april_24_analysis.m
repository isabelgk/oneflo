
figure
for i = 1:3
    D = results.Data{i}(:,2);
    plot(1:197,D)
    hold on
    L{i} = results.Images{i};
end
legend(L)

standard = results2.Data{1}(:,2);
plot(1:197,standard)

averageSaturation = [mean(standard(10:180)),mean(D(10:180,1)),...
    mean(D(10:180,2)), mean(D(10:180,3))]
