P1= load('skypatchfile.P1');
P2= load('skypatchfile.P2');

addpath('P1L1/skypatch_1/');
addpath('P1H1/skypatch_1/');
addpath('P1H2/skypatch_1/');
addpath('P2L1/skypatch_1/');
addpath('P2H1/skypatch_1/');
addpath('P2H2/skypatch_1/');

P1L1M000000
kk=size(map);
pixels =kk(1);

start_x= P1(1)-P1(3)/2;
end_x= P1(1)+P1(3)/2;
step_x= (end_x-start_x)/(pixels-1);
x=start_x:step_x:end_x;

start_y= P1(2)-P1(4)/2;
end_y= P1(2)+P1(4)/2;
step_y= (end_y-start_y)/(pixels-1);
y=start_y:step_y:end_y;
yy=end_y:-step_y:start_y;

subplot(3,2,1)

imagesc(x,yy,map),colorbar
colormap(cool)
xlabel('Right ascension (radians)')
ylabel('Declination (radians)')
title('P1')

P1H1M000000
subplot(3,2,3)
imagesc(x,yy,map),colorbar
xlabel('Right ascension (radians)')
ylabel('Declination (radians)')

P1H2M000000
subplot(3,2,5)
imagesc(x,yy,map),colorbar
xlabel('Right ascension (radians)')
ylabel('Declination (radians)')

%%%%%%%%%%%%%%%%%%%%%%%%%%%
P2L1M000000
kk=size(map);
pixels =kk(1);

start_x= P2(1)-P2(3)/2;
end_x= P2(1)+P2(3)/2;
step_x= (end_x-start_x)/(pixels-1);
x=start_x:step_x:end_x;

start_y= P2(2)-P2(4)/2;
end_y= P2(2)+P2(4)/2;
step_y= (end_y-start_y)/(pixels-1);
y=start_y:step_y:end_y;
yy=end_y:-step_y:start_y;

subplot(3,2,2)

imagesc(x,yy,map),colorbar
colormap(cool)
xlabel('Right ascension (radians)')
ylabel('Declination (radians)')
title('P2')

P2H1M000001
subplot(3,2,4)
imagesc(x,yy,map),colorbar
xlabel('Right ascension (radians)')
ylabel('Declination (radians)')

P2H2M000001
subplot(3,2,6)
imagesc(x,yy,map),colorbar
xlabel('Right ascension (radians)')
ylabel('Declination (radians)')
