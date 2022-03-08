format long
N=32;
M=32;
e=2^-6;
d=1;
T=2;
s_0=2.2;
s_1=min(1/2,(e*log(N)*s_0));
%s_2=min((1-d)/2,e*log(N)*s_0);
dt=T/M;
t=zeros(1,N+1);
for n=1:M+1
    t(n)=(n-1)*dt;
end
x=zeros(N+1,1);
for i=1:(N/4)+1
    h(i)=((1-s_1)*4)/N;
   x(i)=(i-1)*h(i);
end
for i=(N/4)+2:(N/2)+1
    h(i)=(s_1*4)/N;
    x(i)=(1-s_1)+(i-((N/4)+1))*h(i);
end
for i=(N/2)+2:((3*N)/4)+1
    h(i)=(s_1*4)/N;
    x(i)=1+(i-((N/2)+1))*h(i);
end
for i=((3*N)/4)+2:N+1
    h(i)=((1-s_1)*4)/N;
    x(i)=(1+s_1)+(i-(((3*N)/4)+1))*h(i);
end
for i=2:N
    H(i)=h(i)+h(i+1);
end
u=zeros(N+1,M+1);
for n=1:M
    u(1,n+1)=0;
    u(N+1,n+1)=0;
end
for i=1:N+1
    
    u(i,1)=0;
end
a=zeros(1,N+1);
for i=1:(N/2)+1
    a(i)=-(4+x(i)^2);
end
for i=(N/2)+2:N+1
    a(i)=(6-x(i)^2);
end
f=zeros(N+1,M+1);
for n=1:M
for i=1:(N/2)+1
    f(i,n+1)=4*x(i)*exp(-t(n+1))*t(n+1)^2;
end
for i=(N/2)+2:N+1
     f(i,n+1)=4*(2-x(i))*exp(-t(n+1))*t(n+1)^2;
end
end
b=zeros(1,N+1);
for i=1:N+1
    b(i)=5;
end
c=zeros(1,N+1);
for i=1:N+1
    c(i)=2;
end
r_1=zeros(1,N+1);
r_2=zeros(1,N+1);
r_3=zeros(1,N+1);    
p_1=zeros(1,N+1);
p_2=zeros(1,N+1);   
p_3=zeros(1,N+1);
m_1=zeros(1,N+1);
m_2=zeros(1,N+1);
m_3=zeros(1,N+1);
l_1=zeros(1,N+1);
l_2=zeros(1,N+1); 
for i=2:(N/4)+1
    r_1(i)=((2*e)/(h(i)*H(i)))-(a(i-1)+a(i))/(2*h(i))-((b(i-1)+b(i))/4)-1/(2*dt);
    r_2(i)=-(2*e)/(h(i)*h(i+1))+(a(i-1)+a(i))/(2*h(i))-((b(i-1)+b(i))/4)-1/(2*dt);
    r_3(i)=(2*e)/(h(i+1)*H(i));
    p_1(i)=1/(2*dt);
    p_2(i)=1/(2*dt);
    p_3(i)=0;
    m_1(i)=1/2;
    m_2(i)=1/2;
    m_3(i)=0;
    l_1(i)=0;
    l_2(i)=(c(i-1)+c(i))/2;
end
for i=(N/4)+2:N/2
    r_1(i)=((2*e)/(h(i)*H(i)))-(a(i)/H(i));
    r_2(i)=(-(2*e)/(h(i)*h(i+1)))-b(i)-1/dt;
    r_3(i)=((2*e)/(h(i+1)*H(i)))+(a(i)/H(i));
    p_1(i)=0;
    p_2(i)=1/(dt);
    p_3(i)=0;
    m_1(i)=0;
    m_2(i)=1;
    m_3(i)=0;
    l_1(i)=0;
    l_2(i)=c(i);
end
for i=(N/2)+2:(3*N)/4
     r_1(i)=((2*e)/(h(i)*H(i)))-(a(i)/H(i));
    r_2(i)=(-(2*e)/(h(i)*h(i+1)))-b(i)-1/dt;
    r_3(i)=((2*e)/(h(i+1)*H(i)))+(a(i)/H(i));
    p_1(i)=0;
    p_2(i)=1/(dt);
    p_3(i)=0;
    m_1(i)=0;
    m_2(i)=1;
    m_3(i)=0;
    l_1(i)=c(i);
    l_2(i)=0;
end
for i=((3*N)/4)+1:N
     r_1(i)=(2*e)/(h(i)*H(i));
    r_2(i)=(-(2*e)/(h(i)*h(i+1)))-(a(i+1)+a(i))/(2*h(i+1))-((b(i+1)+b(i))/4)-1/(2*dt);
    r_3(i)=((2*e)/(h(i+1)*H(i)))+(a(i+1)+a(i))/(2*h(i+1))-((b(i+1)+b(i))/4)-1/(2*dt);
    p_1(i)=0;
    p_2(i)=1/(2*dt);
    p_3(i)=1/(2*dt);
    m_1(i)=0;
    m_2(i)=1/2;
    m_3(i)=1/2;
    l_1(i)=(c(i)+c(i+1))/2;
    l_2(i)=0;
end
q_1=-N/(8*s_1);
q_2=N/(2*s_1);
q_3=-3/2*(N/(4*s_1)+N/(4*s_1));
q_4=N/(2*s_1);
q_5=-N/(8*s_1);
A=zeros(N+1,N+1);
A(1,1)=1;
A(N+1,N+1)=1;
for n=1:M
    
for i=2:N/2
    A(i,i-1)=r_1(i);
    A(i,i)=r_2(i);
    A(i,i+1)=r_3(i);
end
for i=(N/2)+2:N
    A(i,i-1)=r_1(i);
    A(i,i)=r_2(i);
    A(i,i+1)=r_3(i);
end
A((N/2)+1,(N/2)-1)=q_1;
A((N/2)+1,N/2)=q_2;
A((N/2)+1,(N/2)+1)=q_3;
A((N/2)+1,(N/2)+2)=q_4;
A((N/2)+1,(N/2)+3)=q_5;
B=zeros(N+1,N+1);
for i=2:N/2
    B(i,i-1)=-p_1(i);
    B(i,i)=-p_2(i);
    B(i,i+1)=-p_3(i);
end
for i=(N/2)+2:N
      B(i,i-1)=-p_1(i);
    B(i,i)=-p_2(i);
    B(i,i+1)=-p_3(i);
end
C=zeros(N+1,1);
for i=1:N+1
    C(i)=u(i,n);
end
D=B*C;
E=zeros(N+1,N+1);
for i=2:N/2
    E(i,i-1)=m_1(i);
    E(i,i)=m_2(i);
    E(i,i+1)=m_3(i);
end
for i=(N/2)+2:N
    E(i,i-1)=m_1(i);
    E(i,i)=m_2(i);
    E(i,i+1)=m_3(i);
end
F=zeros(N+1,1);
for i=1:N+1
    F(i)=f(i,n+1);
end
J=zeros(N+1,1);
J(1)=0;
for i=2:N/2
    J(i)=(l_2(i))*0;
end
for i=(N/2)+2:N
    J(i)=(l_1(i))*(u(i-((N/2)+1),n+1));
end
G=E*F;
K=D+G+J;
 [L,R]=lu(A);
   u(:,n+1)=R\(L\K);
end

   surfc(t,x,u);
ylabel(' x');
xlabel(' t');
zlabel(' solution (u)');
