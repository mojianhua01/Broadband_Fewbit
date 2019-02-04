clc;

Np = 512 * [2:1:10];

runtime = [26.39	30.28	34.32	38.47	42.51	46.57	51.16	55.16	59.14
2.47	3.8	5.22	6.71	8.3	9.26	10.27	12.87	13.61
26.8	31.47	36.78	41.21	45.97	50.55	56.46	60.71	65.23
3.28	4.36	5.4	6.8	7.9	8.96	9.89	12.09	12.69
];

figure,
plot1 = plot(Np, runtime, 'linewidth', 1);
set(plot1(1),'Marker','o','LineStyle','--');
set(plot1(2),'Marker','square','LineStyle','-');
set(plot1(3),'Marker','diamond','LineStyle','--');
set(plot1(4),'Marker','v','LineStyle','-');
grid on;
xlabel('Training length', 'fontsize',14)
ylabel('Runtime (ms)', 'fontsize',14)
handle  = legend('$$\mathbf{A} \widehat{\mathbf{x}}^{(n)}$$', ...
    'Fast implementation of $$\mathbf{A} \widehat{\mathbf{x}}^{(n)}$$',...
    '$$\mathbf{A}^{*} \widehat{\mathbf{s}}^{(n)}$$',...
    'Fast implementation of $$\mathbf{A}^{*} \widehat{\mathbf{s}}^{(n)}$$')
set(handle, 'interpreter','latex','fontsize',10)
set(gca, 'XTick', Np)