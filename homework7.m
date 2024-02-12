%常微分方程数值求解
%先给出解析解的曲线
h=10^-3;
t=0:h:20;
ut=-1.499875*exp(-0.5*t)+0.499875*exp(-2000.5*t)+1;%两个方程
vt=-2.99975*exp(-0.5*t)-0.00025*exp(-2000.5*t)+1;
f=@(y) [-2000*y(1)+999.75*y(2)+1000.25;y(1)-y(2)];
figure(1)
plot(t,ut,"r",t,vt,"b")
xlabel("t")
legend("ut","vt","Location","best")
title("解析解的曲线")
%在数值解法结近似解
%显式方法
%f=[-2000u+999.75v+1000.25;u-v]
y1=zeros(2,20/h+1);
y1(:,1)=[0;-2];
for n=1:20/h
    k11=-2000*y1(1,n)+999.75*y1(2,n)+1000.25;
    k12=y1(1,n)-y1(2,n);
    k21=-2000*(y1(1,n)+0.5*h*k11)+999.75*(y1(2,n)+0.5*h*k12)+1000.25;
    k22=y1(1,n)+0.5*h*k11-y1(2,n)-0.5*h*k12;
    k31=-2000*(y1(1,n)+0.5*h*k21)+999.75*(y1(2,n)+0.5*h*k22)+1000.25;
    k32=y1(1,n)+0.5*h*k21-y1(2,n)-0.5*h*k22;
    k41=-2000*(y1(1,n)+h*k31)+999.75*(y1(2,n)+h*k32)+1000.25;
    k42=y1(1,n)+h*k31-y1(2,n)-h*k32;
    k1=[k11;k12];
    k2=[k21;k22];
    k3=[k31;k32];
    k4=[k41;k42];
    y1(:,n+1)=y1(:,n)+h/6*(k1+2*k2+2*k3+k4);
end
figure(2)
plot(0:h:20,y1(1,:),"r",0:h:20,y1(2,:),"b")
xlabel("t")
legend("u","v","Location","best")
erroru1=abs(y1(1,:)-ut);
errorv1=abs(y1(2,:)-vt);
title("四级四阶显式Runge-Kutta")
figure(3)
semilogy(0:h:20,erroru1,"r",0:h:20,errorv1,"b")
legend("error u","error v","Location","best")
title("显式方法误差绝对值")
ylabel("绝对误差")

%隐式方法
h=0.001;
m=2;
ksi=(sort(roots(polyLegendre(m))))';
a=(1+ksi)/2;
L=((1:m).^(-1))';
A=[];
for p=1:m
    A=[A;a.^(p-1)];
end
c=A\L;
AL=[];
for p=1:m
    AL=[AL,((a)'.^p)/p];
end
b=AL/(A');
y2=zeros(2,20/h+1);
y2(:,1)=[0;-2];
D=[-2000,999.75;1,-1];
for p=1:20/h
    k=(eye(2*m)-h*kron(D,b))\(kron(D,[1;1])*y2(:,p)+[1000.25;1000.25;0;0]);%解线性方程组
    y2(1,p+1)=y2(1,p)+h*(k(1:2)')*c;
    y2(2,p+1)=y2(2,p)+h*(k(3:4)')*c;
end
figure(4)
plot(0:h:20,y2(1,:),"r",0:h:20,y2(2,:),"b")
xlabel("t")
legend("u","v","Location","best")
erroru2=abs(y2(1,:)-ut(1:h/0.001:end));
errorv2=abs(y2(2,:)-vt(1:h/0.001:end));
title("二级四阶隐式Runge-Kutta")
figure(5)
semilogy(0:h:20,erroru2,"r",0:h:20,errorv2,"b")
legend("error u","error v","Location","best")
title("隐式方法误差绝对值")
ylabel("绝对误差")

function p=polyLegendre(n)%生成勒让德多项式
p=1;
for k=1:n
    p=conv(p,[1,0,-1]);
end
for k=1:n
    p=p(1:end-1).*(length(p)-1:-1:1);
end
p=p/2^n/factorial(n);
end