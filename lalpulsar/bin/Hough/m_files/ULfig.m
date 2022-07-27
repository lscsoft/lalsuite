clear

figure

mu=load('Hough_MultiIFO_UL.50_1000.txt');
subplot(3,1,1)
plot(mu(:,1), mu(:,2), 'r.')
xlabel('Frequency (Hz)')
ylabel(' Hough Multi-IFO 95% UL')


h1=load('Hough_H1_UL.50_1000.txt');
subplot(3,1,2)
plot(h1(:,1), h1(:,2), 'b.')
xlabel('Frequency (Hz)')
ylabel(' Hough H1 95% UL')


l1=load('Hough_L1_UL.50_1000.txt');
subplot(3,1,3)
plot(l1(:,1), l1(:,2), 'g.')
xlabel('Frequency (Hz)')
ylabel(' Hough L1 95% UL')
