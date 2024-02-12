%数值积分
a=1;
b=3;
x=a:0.001:b;
fx=1./(x.^2).*sin(2*pi./x);
figure(1)
plot(x,fx)
grid on
title("f(x)=1/x^2*sin(2\pi/x)")
xlabel("x")
ylabel("f(x)")

If=-3/4/pi;%该积分的精确值
f=@(x)(1./x.^2.*sin(2*pi./x));
%Gauss-Legendre求积公式
N=4;%分成4个区间
h=(b-a)/N;
I0=0;
for p=1:N
   I0=I0+GaussLegendre(a+(p-1)*h,a+p*h,f); 
end
error1=I0-If;
%Romberg求积算法
epsl=1;
n=1;
q=1/2;
count=1;
I=zeros(10);%用于数据存储
while (epsl>10^-7)
    h=(b-a)/n;
    I(count,1)=h/2*(f(a)+sum(2*f(a+h:h:b-h))+f(b));
    for p=2:count
        I(count,p)=(4^(p-1)*I(count,p-1)-I(count-1,p-1))/(4^(p-1)-1);
    end
    n=n/q;
    if count>=2
        epsl=abs(I(count,count)-I(count-1,count-1));
    end
    count=count+1;
end
error2=I-(I~=0)*If;

function I=GaussLegendre(a,b,f)
P5=polyLegendre(5);
x0=sort(roots(P5));
x1=(a+b)/2+(b-a)/2*x0;
P4=polyLegendre(4);
P51=P5(1:end-1).*(length(P5)-1:-1:1);%5阶勒让德多项式求导
Ak=2/5./(polyval(P4,x0).*polyval(P51,x0));
I=sum(Ak.*f(x1))*(b-a)/2;
end

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