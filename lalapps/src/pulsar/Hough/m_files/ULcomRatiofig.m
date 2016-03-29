clear


figure

subplot(3,1,1)
mu=load('Hough_MultiIFO_UL.50_1000.txt');
mu2=load('Multi_W50_1000predictedUL.11.0');
semilogy(mu(:,1), mu(:,2)./mu2(:,4), 'b.')
axis([50 1000 0.1 100])
xlabel('Frequency (Hz)')
ylabel('Multi-IFO ratio')
%legend('measured','predicted')

subplot(3,1,2)
h1=load('Hough_H1_UL.50_1000.txt');
h12=load('H1_W50_1000predictedUL.11.5');
semilogy(h1(:,1), h1(:,2)./h12(:,4), 'b.')
axis([50 1000 0.1 100])
xlabel('Frequency (Hz)')
ylabel('H1 UL ratio')

subplot(3,1,3)
l1=load('Hough_L1_UL.50_1000.txt');
l12=load('L1_W50_1000predictedUL.11.1');
%note L1 incompleted

for k=1:length(l1)
  fr=l1(k,1);
  index= find (l12(:,1)==fr);
  ratio(k)=l1(k,2)/l12(index,4);
end
 
semilogy(l1(:,1), ratio(:), 'b.')
axis([50 1000 0.1 100])
xlabel('Frequency (Hz)')
ylabel('L1 UL ratio')
