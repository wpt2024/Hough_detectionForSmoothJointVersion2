% 测试霍夫变换
clc
clear
close all 
% =====读取图像======
%I  = imread('rdfn06.bmp');
%I= imread('PaperResult.bmp');
I= imread('test.bmp');%导入图片
%I= imread('sample03.bmp');%导入图片
%I= imread('model3.bmp');%导入图片
%I= imread('model6.jpg');%导入图片
%I= imread('sample01.bmp');%导入图片
%I=imread('joints.jpg');%import the model:model1.jpg
%I=imread('model2.jpg');%import the model:model1.jpg
%I=imread('model4.jpg');%import the model:model1.jpg
%I=imread('joints.jpg');%import the model:model1.jpg
%I=imread('model5.jpg');%import the model:model1.jpg
%I=imread('model6.jpg');%import the model:model1.jpg
%I= imread('1.bmp');%导入图片


%=========================
lengthofImage=length(I)
rotI1 = imrotate(I,0,'crop'); % 旋转任意角度 如45度，保持原图片大小
fig1 = imshow(rotI1);

rotI=rgb2gray(rotI1)%灰度图
rotI=255-rotI%颜色转置
%===
newI=imbinarize(rotI);%图像二值化
%===图像骨架化=
%newThinI=bwmorph(newI, 'thin', Inf);%图像骨架化
newThinI=bwmorph(newI, 'skel', Inf);%图像骨架化
%newThinI=bwmorph(newI);%图像骨架化
%newThinI=bwmorph(newI, 'MinBracnchLength', 1);%图像骨架化 去掉较短的分叉

%===
%rotI=im2bw(rotI1)%二值化图像
%rotI=1-rotI%二值化图像

%imput the parameters
number1=input('Please input the number of right-incline joints:  (最小为1)')
number2=input('Please input the number of left-incline joints:  (最小为1)')
[row, colume]=size(newThinI)%输出I的大小
rockmass2Dsize=row*colume%计算面积-根据像素
gapvalue=input('输入断距(size=size(I)):  (比如5)')
minilengthvalue=input('输入最小迹长:  (比如20)')
thelengthofimage=input('输入所导入图片的长度： （图片实际尺寸 /m）')%
thewidthofimage=input('输入所导入图片的宽度： （图片实际尺寸/m）')%
realarea=thelengthofimage*thewidthofimage;%实际的面积

%tracelength1=input('Please input the test-length of right-incline joints:  (一般为1)')
%tracelength2=input('Please input the test-length of left-incline joints:  (一般为1)')
% rotI=im2bw(rotI1)%图片二值化
% 提取边
BW = edge(newThinI,'canny'); %提取边缘
%BW = edge(rotI,'prewitt');
figure
hold on;
subplot(1,2,1)
imshow(BW); 
subplot(1,2,2);
imshow(newThinI);
% 霍夫变换

%右倾节理迹线识别
%[H1,theta1,rho1] = hough(BW,'RhoResolution',1,'Theta',0:0.01:89.999); % 计算二值图像的标准霍夫变换，H为霍夫变换矩阵，theta,rho为计算霍夫变换的角度和半径值
[H1,theta1,rho1] = hough(newThinI,'RhoResolution',1,'Theta',0:0.01:89.999); % 计算二值图像的标准霍夫变换，H为霍夫变换矩阵，theta,rho为计算霍夫变换的角度和半径值

%H1为hough transform矩阵，每个元素分别存放的相应参数对应的累加
%Theta1为存放theta的采样值的向量
%Rho1为存放rho采样值的向量
k_r=size(H1,1)%返回H1的行数 即代表迹线条数

figure, 
%subplot(3,1,2); 
imshow(imadjust(mat2gray(H1)),[],'XData',theta1,'YData',rho1,...
    'InitialMagnification','fit');
xlabel('\theta1 (degrees)'), ylabel('\rho1');
axis on, axis normal, hold on;
colormap(hot) 
% 显示霍夫变换矩阵中的极值点
% peaks1= houghpeaks(H1, number1) %自定义迹线条数
peaks1 = houghpeaks(H1, k_r) %根据极值点来搜索，极值点一般和线条像素有关，但是线条的最大数目不超过这个数

%P1 = houghpeaks(H1,100,'threshold',ceil(0.3*max(H1(:)))); % 从霍夫变换矩阵H中提取100个极值点
%P1 = houghpeaks(H1,number1,'threshold',ceil(0.3*max(H1(:)))); %
%从霍夫变换矩阵H中提取100个极值点%自定义的条数提取
P1 = houghpeaks(H1,k_r,'threshold',ceil(0.3*max(H1(:)))); % 从霍夫变换矩阵H中提取100个极值点%提取出最大的极值点

x = theta1(P1(:,2));%极值点的theta值，即P的第二列存放的是极值点的theta值
y = rho1(P1(:,1));%极值点的rho值，即P的第二列存放的是极值点的rho值
plot(x,y,'s','color','black'); 
% 找原图中的线
%lines1 = houghlines(BW,theta1,rho1,P1,'FillGap',gapvalue,'MinLength',minilengthvalue);%'FillGap'：当两条直线之间距离小于该阈值时，两条直线被合并为一条直线；'MinLength'：保留的直线的最短长度（即小于该阈值的直线将被删除）
lines1 = houghlines(newThinI,theta1,rho1,P1,'FillGap',gapvalue,'MinLength',minilengthvalue);%'FillGap'：当两条直线之间距离小于该阈值时，两条直线被合并为一条直线；'MinLength'：保留的直线的最短长度（即小于该阈值的直线将被删除）

figure 
imshow(newThinI), hold on
max_len = 0;
for k = 1:length(lines1)
    % 绘制各条线
    xy = [lines1(k).point1; lines1(k).point2];
    plot(xy(:,1),xy(:,2),'LineWidth',2,'Color','blue');   
    % 绘制线的起点（黄色）、终点（红色）
    plot(xy(1,1),xy(1,2),'x','LineWidth',1,'Color','green');
    plot(xy(2,1),xy(2,2),'x','LineWidth',1,'Color','red');   
    % 计算线的长度，找最长线段
    len = norm(lines1(k).point1 - lines1(k).point2);
    if ( len > max_len)
        max_len = len;
        xy_long = xy;
    end
end
% 以红色线高亮显示最长的线
if ~isempty(xy_long)
  plot(xy_long(:,1),xy_long(:,2),'LineWidth',2,'Color','red');
end


% %左倾节理迹线识别
[H2,theta2,rho2] = hough(newThinI,'RhoResolution',1,'Theta',-0.001:-0.01:-89.99); % 计算二值图像的标准霍夫变换，H为霍夫变换矩阵，theta,rho为计算霍夫变换的角度和半径值
figure, 
imshow(imadjust(mat2gray(H2)),[],'XData',theta2,'YData',rho2,...
    'InitialMagnification','fit');
xlabel('\theta (degrees)'), ylabel('\rho');
axis on, axis normal, hold on;
colormap(hot) 
% 显示霍夫变换矩阵中的极值点
k_l=size(H2,1)%返回H2的行数 即代表迹线条数
peaks2 = houghpeaks(H2, k_l) 
%peaks2 = houghpeaks(H2, number2) 

%P2 = houghpeaks(H2,5,'threshold',ceil(0.3*max(H2(:)))); % 从霍夫变换矩阵H中提取5个极值点
%P2 = houghpeaks(H2,number2,'threshold',ceil(0.3*max(H2(:)))); % 从霍夫变换矩阵H中提取5个极值点
P2 = houghpeaks(H2,k_l,'threshold',ceil(0.3*max(H2(:)))); % 从霍夫变换矩阵H中提取5个极值点

x = theta2(P2(:,2));
y = rho2(P2(:,1));
plot(x,y,'s','color','black'); 
% 找原图中的线
lines2 = houghlines(newThinI,theta2,rho2,P2,'FillGap',gapvalue,'MinLength',minilengthvalue);
figure, 
imshow(newThinI), hold on
max_len = 0;
for k = 1:length(lines2)
    % 绘制各条线
    xy = [lines2(k).point1; lines2(k).point2];
    plot(xy(:,1),xy(:,2),'LineWidth',2,'Color','green');   
    % 绘制线的起点（黄色）、终点（红色）
    plot(xy(1,1),xy(1,2),'x','LineWidth',1,'Color','green');
    plot(xy(2,1),xy(2,2),'x','LineWidth',1,'Color','red');   
    % 计算线的长度，找最长线段
    len = norm(lines2(k).point1 - lines2(k).point2);
    if ( len > max_len)
        max_len = len;
        xy_long = xy;
    end
end
% 以红色线高亮显示最长的线
plot(xy_long(:,1),xy_long(:,2),'LineWidth',2,'Color','red');

%画总迹线识别图
figure, imshow(rotI), hold on
max_len = 0;
linestotal=[lines1,lines2]
for k = 1:length(linestotal)
    % 绘制各条线
    xy = [linestotal(k).point1; linestotal(k).point2];
    plot(xy(:,1),xy(:,2),'LineWidth',2,'Color','blue');   
    % 绘制线的起点（黄色）、终点（红色）
    plot(xy(1,1),xy(1,2),'x','LineWidth',2,'Color','green');
    plot(xy(2,1),xy(2,2),'x','LineWidth',2,'Color','red');   
    % 计算线的长度，找最长线段
    len = norm(linestotal(k).point1 - linestotal(k).point2);
    if ( len > max_len)
        max_len = len;%输出最大的节理迹线长度
        xy_long = xy;
    end
	k
end
% 以红色线高亮显示最长的线
if ~isempty(xy_long)
  plot(xy_long(:,1),xy_long(:,2),'LineWidth',3,'Color','blue');
end
%画倾向倾角极点图
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
%Data output to files倾向倾角导出
filename = 'dipsdata.xlsx';
datadips(:,1)=dipd(1,:);
datadips(:,2)=dipa(1,:);
%datadips=unique(datadips)
%xlswrite(filename,datadips)

%======显示必要的结果======
max_len%最大节理长度
k%输出识别出的节理迹线数目
rockmass2Dsize%2D area
realarea%real 2D area 
density2D=k/rockmass2Dsize%the 2D density of joint traces
realdensity2D=k/realarea%

%end of file



