clear


figure

subplot(3,1,1)
mu=load('Hough_MultiIFO_UL.50_1000.txt');
plot(mu(:,1), mu(:,2), 'b-')
hold on
mu2=load('Multi_W50_1000predictedUL.11.0');
plot(mu2(:,1), mu2(:,4), 'r-')
xlabel('Frequency (Hz)')
ylabel(' Hough Multi-IFO 95% UL')
legend('measured','predicted')

subplot(3,1,2)
h1=load('Hough_H1_UL.50_1000.txt');
plot(h1(:,1), h1(:,2), 'b-')
hold on
h12=load('H1_W50_1000predictedUL.11.5');
plot(h12(:,1), h12(:,4), 'r-')
xlabel('Frequency (Hz)')
ylabel(' Hough H1 95% UL')
legend('measured','predicted')

subplot(3,1,3)
l1=load('Hough_L1_UL.50_1000.txt');
plot(l1(:,1), l1(:,2), 'b-')
hold on
l12=load('L1_W50_1000predictedUL.11.1');
plot(l12(:,1), l12(:,4), 'r-')
xlabel('Frequency (Hz)')
ylabel(' Hough L1 95% UL')
legend('measured','predicted')
