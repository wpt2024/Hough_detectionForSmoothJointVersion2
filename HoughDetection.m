% ���Ի���任
clc
clear
close all 
% =====��ȡͼ��======
%I  = imread('rdfn06.bmp');
%I= imread('PaperResult.bmp');
I= imread('test.bmp');%����ͼƬ
%I= imread('sample03.bmp');%����ͼƬ
%I= imread('model3.bmp');%����ͼƬ
%I= imread('model6.jpg');%����ͼƬ
%I= imread('sample01.bmp');%����ͼƬ
%I=imread('joints.jpg');%import the model:model1.jpg
%I=imread('model2.jpg');%import the model:model1.jpg
%I=imread('model4.jpg');%import the model:model1.jpg
%I=imread('joints.jpg');%import the model:model1.jpg
%I=imread('model5.jpg');%import the model:model1.jpg
%I=imread('model6.jpg');%import the model:model1.jpg
%I= imread('1.bmp');%����ͼƬ


%=========================
lengthofImage=length(I)
rotI1 = imrotate(I,0,'crop'); % ��ת����Ƕ� ��45�ȣ�����ԭͼƬ��С
fig1 = imshow(rotI1);

rotI=rgb2gray(rotI1)%�Ҷ�ͼ
rotI=255-rotI%��ɫת��
%===
newI=imbinarize(rotI);%ͼ���ֵ��
%===ͼ��Ǽܻ�=
%newThinI=bwmorph(newI, 'thin', Inf);%ͼ��Ǽܻ�
newThinI=bwmorph(newI, 'skel', Inf);%ͼ��Ǽܻ�
%newThinI=bwmorph(newI);%ͼ��Ǽܻ�
%newThinI=bwmorph(newI, 'MinBracnchLength', 1);%ͼ��Ǽܻ� ȥ���϶̵ķֲ�

%===
%rotI=im2bw(rotI1)%��ֵ��ͼ��
%rotI=1-rotI%��ֵ��ͼ��

%imput the parameters
number1=input('Please input the number of right-incline joints:  (��СΪ1)')
number2=input('Please input the number of left-incline joints:  (��СΪ1)')
[row, colume]=size(newThinI)%���I�Ĵ�С
rockmass2Dsize=row*colume%�������-��������
gapvalue=input('����Ͼ�(size=size(I)):  (����5)')
minilengthvalue=input('������С����:  (����20)')
thelengthofimage=input('����������ͼƬ�ĳ��ȣ� ��ͼƬʵ�ʳߴ� /m��')%
thewidthofimage=input('����������ͼƬ�Ŀ�ȣ� ��ͼƬʵ�ʳߴ�/m��')%
realarea=thelengthofimage*thewidthofimage;%ʵ�ʵ����

%tracelength1=input('Please input the test-length of right-incline joints:  (һ��Ϊ1)')
%tracelength2=input('Please input the test-length of left-incline joints:  (һ��Ϊ1)')
% rotI=im2bw(rotI1)%ͼƬ��ֵ��
% ��ȡ��
BW = edge(newThinI,'canny'); %��ȡ��Ե
%BW = edge(rotI,'prewitt');
figure
hold on;
subplot(1,2,1)
imshow(BW); 
subplot(1,2,2);
imshow(newThinI);
% ����任

%���������ʶ��
%[H1,theta1,rho1] = hough(BW,'RhoResolution',1,'Theta',0:0.01:89.999); % �����ֵͼ��ı�׼����任��HΪ����任����theta,rhoΪ�������任�ĽǶȺͰ뾶ֵ
[H1,theta1,rho1] = hough(newThinI,'RhoResolution',1,'Theta',0:0.01:89.999); % �����ֵͼ��ı�׼����任��HΪ����任����theta,rhoΪ�������任�ĽǶȺͰ뾶ֵ

%H1Ϊhough transform����ÿ��Ԫ�طֱ��ŵ���Ӧ������Ӧ���ۼ�
%Theta1Ϊ���theta�Ĳ���ֵ������
%Rho1Ϊ���rho����ֵ������
k_r=size(H1,1)%����H1������ ������������

figure, 
%subplot(3,1,2); 
imshow(imadjust(mat2gray(H1)),[],'XData',theta1,'YData',rho1,...
    'InitialMagnification','fit');
xlabel('\theta1 (degrees)'), ylabel('\rho1');
axis on, axis normal, hold on;
colormap(hot) 
% ��ʾ����任�����еļ�ֵ��
% peaks1= houghpeaks(H1, number1) %�Զ��弣������
peaks1 = houghpeaks(H1, k_r) %���ݼ�ֵ������������ֵ��һ������������йأ����������������Ŀ�����������

%P1 = houghpeaks(H1,100,'threshold',ceil(0.3*max(H1(:)))); % �ӻ���任����H����ȡ100����ֵ��
%P1 = houghpeaks(H1,number1,'threshold',ceil(0.3*max(H1(:)))); %
%�ӻ���任����H����ȡ100����ֵ��%�Զ����������ȡ
P1 = houghpeaks(H1,k_r,'threshold',ceil(0.3*max(H1(:)))); % �ӻ���任����H����ȡ100����ֵ��%��ȡ�����ļ�ֵ��

x = theta1(P1(:,2));%��ֵ���thetaֵ����P�ĵڶ��д�ŵ��Ǽ�ֵ���thetaֵ
y = rho1(P1(:,1));%��ֵ���rhoֵ����P�ĵڶ��д�ŵ��Ǽ�ֵ���rhoֵ
plot(x,y,'s','color','black'); 
% ��ԭͼ�е���
%lines1 = houghlines(BW,theta1,rho1,P1,'FillGap',gapvalue,'MinLength',minilengthvalue);%'FillGap'��������ֱ��֮�����С�ڸ���ֵʱ������ֱ�߱��ϲ�Ϊһ��ֱ�ߣ�'MinLength'��������ֱ�ߵ���̳��ȣ���С�ڸ���ֵ��ֱ�߽���ɾ����
lines1 = houghlines(newThinI,theta1,rho1,P1,'FillGap',gapvalue,'MinLength',minilengthvalue);%'FillGap'��������ֱ��֮�����С�ڸ���ֵʱ������ֱ�߱��ϲ�Ϊһ��ֱ�ߣ�'MinLength'��������ֱ�ߵ���̳��ȣ���С�ڸ���ֵ��ֱ�߽���ɾ����

figure 
imshow(newThinI), hold on
max_len = 0;
for k = 1:length(lines1)
    % ���Ƹ�����
    xy = [lines1(k).point1; lines1(k).point2];
    plot(xy(:,1),xy(:,2),'LineWidth',2,'Color','blue');   
    % �����ߵ���㣨��ɫ�����յ㣨��ɫ��
    plot(xy(1,1),xy(1,2),'x','LineWidth',1,'Color','green');
    plot(xy(2,1),xy(2,2),'x','LineWidth',1,'Color','red');   
    % �����ߵĳ��ȣ�����߶�
    len = norm(lines1(k).point1 - lines1(k).point2);
    if ( len > max_len)
        max_len = len;
        xy_long = xy;
    end
end
% �Ժ�ɫ�߸�����ʾ�����
if ~isempty(xy_long)
  plot(xy_long(:,1),xy_long(:,2),'LineWidth',2,'Color','red');
end


% %���������ʶ��
[H2,theta2,rho2] = hough(newThinI,'RhoResolution',1,'Theta',-0.001:-0.01:-89.99); % �����ֵͼ��ı�׼����任��HΪ����任����theta,rhoΪ�������任�ĽǶȺͰ뾶ֵ
figure, 
imshow(imadjust(mat2gray(H2)),[],'XData',theta2,'YData',rho2,...
    'InitialMagnification','fit');
xlabel('\theta (degrees)'), ylabel('\rho');
axis on, axis normal, hold on;
colormap(hot) 
% ��ʾ����任�����еļ�ֵ��
k_l=size(H2,1)%����H2������ ������������
peaks2 = houghpeaks(H2, k_l) 
%peaks2 = houghpeaks(H2, number2) 

%P2 = houghpeaks(H2,5,'threshold',ceil(0.3*max(H2(:)))); % �ӻ���任����H����ȡ5����ֵ��
%P2 = houghpeaks(H2,number2,'threshold',ceil(0.3*max(H2(:)))); % �ӻ���任����H����ȡ5����ֵ��
P2 = houghpeaks(H2,k_l,'threshold',ceil(0.3*max(H2(:)))); % �ӻ���任����H����ȡ5����ֵ��

x = theta2(P2(:,2));
y = rho2(P2(:,1));
plot(x,y,'s','color','black'); 
% ��ԭͼ�е���
lines2 = houghlines(newThinI,theta2,rho2,P2,'FillGap',gapvalue,'MinLength',minilengthvalue);
figure, 
imshow(newThinI), hold on
max_len = 0;
for k = 1:length(lines2)
    % ���Ƹ�����
    xy = [lines2(k).point1; lines2(k).point2];
    plot(xy(:,1),xy(:,2),'LineWidth',2,'Color','green');   
    % �����ߵ���㣨��ɫ�����յ㣨��ɫ��
    plot(xy(1,1),xy(1,2),'x','LineWidth',1,'Color','green');
    plot(xy(2,1),xy(2,2),'x','LineWidth',1,'Color','red');   
    % �����ߵĳ��ȣ�����߶�
    len = norm(lines2(k).point1 - lines2(k).point2);
    if ( len > max_len)
        max_len = len;
        xy_long = xy;
    end
end
% �Ժ�ɫ�߸�����ʾ�����
plot(xy_long(:,1),xy_long(:,2),'LineWidth',2,'Color','red');

%���ܼ���ʶ��ͼ
figure, imshow(rotI), hold on
max_len = 0;
linestotal=[lines1,lines2]
for k = 1:length(linestotal)
    % ���Ƹ�����
    xy = [linestotal(k).point1; linestotal(k).point2];
    plot(xy(:,1),xy(:,2),'LineWidth',2,'Color','blue');   
    % �����ߵ���㣨��ɫ�����յ㣨��ɫ��
    plot(xy(1,1),xy(1,2),'x','LineWidth',2,'Color','green');
    plot(xy(2,1),xy(2,2),'x','LineWidth',2,'Color','red');   
    % �����ߵĳ��ȣ�����߶�
    len = norm(linestotal(k).point1 - linestotal(k).point2);
    if ( len > max_len)
        max_len = len;%������Ľ����߳���
        xy_long = xy;
    end
	k
end
% �Ժ�ɫ�߸�����ʾ�����
if ~isempty(xy_long)
  plot(xy_long(:,1),xy_long(:,2),'LineWidth',3,'Color','blue');
end
%��������Ǽ���ͼ
figure
dipangle1=90-[lines1.theta];
dipdirec1=[lines1.theta]+90;
dipangle2=[lines2.theta]+90;
dipdirec2=270+[lines2.theta];
dipa=[dipangle1,dipangle2];
dipd=[dipdirec1,dipdirec2];
dipdrad=deg2rad(dipd);

%hhh = my_polar([0 2*pi],[0 90]);hold on% 
%polarplot(dipdrad,dipa,'r *')
my_polar(dipdrad,dipa,'b o')
%rose(dipdrad,dipa,'b *');hold on
%Data output to files������ǵ���
filename = 'dipsdata.xlsx';
datadips(:,1)=dipd(1,:);
datadips(:,2)=dipa(1,:);
%datadips=unique(datadips)
%xlswrite(filename,datadips)

%======��ʾ��Ҫ�Ľ��======
max_len%��������
k%���ʶ����Ľ�������Ŀ
rockmass2Dsize%2D area
realarea%real 2D area 
density2D=k/rockmass2Dsize%the 2D density of joint traces
realdensity2D=k/realarea%

%end of file



