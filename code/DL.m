clc;clear;

OffsetPos = 200;
OffsetVal = 200;

FileName = 'DL20201103_142426_100HZ������';
fid=fopen([FileName,'.txt'],'r');
FormatString=repmat('%s ',1,4);
dat=textscan(fid,FormatString); 
len = length(dat{1});

dat1(len) = 0;
dat2(len) = 0;  

for i = 1:len
    d1 = hex2dec(dat{1}{i});
    d2 = hex2dec(dat{2}{i});
    d3 = hex2dec(dat{3}{i});
    d4 = hex2dec(dat{4}{i});
    
    dat1(i) = (d1 -160) * 256 + d2;
    dat2(i) = (d3 -176) * 256 + d4;
end

left = dat1';
right = dat2';
xlsTable = table(left, right);
writetable(xlsTable, [FileName,'.csv']);

startPos = 1;
endPos = len;

for i = 1:1:len
    if (dat1(i) > OffsetVal) || (dat2(i) > OffsetVal)
        if i > OffsetPos
            startPos = i - OffsetPos;
        end
        break;
    end
end

for i = len:-1:1
    if (dat1(i) > OffsetVal) || (dat2(i) > OffsetVal)
        if i < len - OffsetPos
            endPos = i + OffsetPos;
        end
        break;
    end
end

figure;
plot(dat2(startPos:endPos),'b');
hold on;
plot(dat1(startPos:endPos),'r');
hold on;


return;



angle1 = dat(:,1);
angle2 = dat(:,2);




x = angle1 - angle2;
for i=1:length(x)
    if x(i) > 270
        x(i) = x(i) - 360;
    end
    if x(i) < -270
        x(i) = x(i) + 360;
    end
    
    x(i) = x(i) * 30 + 180;
end

plot(angle1,'r');
hold on;
plot(angle2,'g');
hold on;
plot(x,'b');
title('תȦ���Ƕ�ƫ��Ա�');
legend('���Ƶ�оƬ����','������ģ�鷽��δ�����','�ǶȲ�*30+180');


return;

j = 1;
cha(1) = 0;
for i=2:length(angle)
    x = angle(i) - angle(i-1);
    
    if x >= 300
        x = x - 360;
    end
    if x <= -5
        cha(j) = x;
        j = j + 1;
    end
end

plot(cha);
title('תȦ�����У�200ms�ɼ�һ�Σ�����ÿ��200ms,ʵ��ת���ĽǶ�');
mean(cha)
sum(cha)/5

return;

t=[0:pi/180:2*pi];

%FileName = '��ѱ��ϲ���3';
%FileName = '�����ϲ���3';
%FileName = '���°�ѱ��ϲ���3';
FileName = '�ֶ�1';


% ��ȡ����
dat=importdata([FileName,'.txt']); 
% GYRO_x = dat(:,1);
% GYRO_y = dat(:,2);
% GYRO_z = dat(:,3);

MAG_x  = dat(:,1);
MAG_y  = dat(:,2);
MAG_z  = dat(:,3);

k = (max(MAG_y)-min(MAG_y)) / (max(MAG_x)-min(MAG_x))
%MAG_x = k * MAG_x;

m = [MAG_x MAG_y ones(size(MAG_x))]\[-(MAG_x.^2+MAG_y.^2)];
xc = -.5*m(1)
yc = -.5*m(2)
R  = sqrt((m(1)^2+m(2)^2)/4-m(3))

figure;
plot(MAG_x, MAG_y, '*');
hold on;
plot(xc,yc,'r-x',(xc+R*cos(t)),(yc+R*sin(t)),'r-');
axis equal;






return;

for i = 1 : length(Cha)
    if Cha(i) > 300
        Cha(i) = Cha(i) - 360;
    end
	if Cha(i) < -300
        Cha(i) = Cha(i) + 360;
    end
end

for i = 1 : length(Filter_Cha)
    if Filter_Cha(i) > 300
        Filter_Cha(i) = Filter_Cha(i) - 360;
    end
	if Filter_Cha(i) < -300
        Filter_Cha(i) = Filter_Cha(i) + 360;
    end
end

figure;
plot(Angle);
hold on;
plot(Filter_Angle);
title([FileName, '   δ�˲����˲��ǶȶԱ�']);

set(gcf,'outerposition',get(0,'screensize'));
saveas(gcf, [FileName, '   δ�˲����˲��ǶȶԱ�'], 'png')
%close;

figure;
plot(Cha);
hold on;
plot(Filter_Cha);
title([FileName, '   δ�˲����˲��ǶȲ�Ա�']);

set(gcf,'outerposition',get(0,'screensize'));
saveas(gcf, [FileName, '   δ�˲����˲��ǶȲ�Ա�'], 'png')
%close;

return;

m = [x y ones(size(x))]\[-(x.^2+y.^2)];
xc = -.5*m(1)
yc = -.5*m(2)
R  = sqrt((m(1)^2+m(2)^2)/4-m(3))

figure;
plot(x,y);
hold on;
plot(xc,yc,'r-x',(xc+R*cos(t)),(yc+R*sin(t)),'r-');
axis equal;

R2=sqrt(x.*x+y.*y)-R;
figure;
plot(R2)

