
function [S NN MM N M XX YY]=Surface_3(X_i_3,y_path_i_3,n_cand_3,m_cand_3,SIGMA_i_3)

syms x y
BETA_i_1=atan((n_cand_3(1)-X_i_3)/(m_cand_3(1)-y_path_i_3))+SIGMA_i_3(1);
BETA_i_2=atan((n_cand_3(1)-X_i_3)/(m_cand_3(1)-y_path_i_3))-SIGMA_i_3(1);
BETA_i_3=atan((n_cand_3(2)-X_i_3)/(m_cand_3(2)-y_path_i_3))+SIGMA_i_3(2);
BETA_i_4=atan((n_cand_3(2)-X_i_3)/(m_cand_3(2)-y_path_i_3))-SIGMA_i_3(2);
BETA_i_5=atan((n_cand_3(3)-X_i_3)/(m_cand_3(3)-y_path_i_3))+SIGMA_i_3(3);
BETA_i_6=atan((n_cand_3(3)-X_i_3)/(m_cand_3(3)-y_path_i_3))-SIGMA_i_3(3);

BETA_i_1=tan(BETA_i_1);
BETA_i_2=tan(BETA_i_2);
BETA_i_3=tan(BETA_i_3);
BETA_i_4=tan(BETA_i_4);
BETA_i_5=tan(BETA_i_5);
BETA_i_6=tan(BETA_i_6);

L_1=x-(BETA_i_1*(y-m_cand_3(1))+n_cand_3(1));
if BETA_i_1==inf
L_1=y-m_cand_3(1);
end

L_2=x-(BETA_i_2*(y-m_cand_3(1))+n_cand_3(1));
if BETA_i_2==inf
    L_2=y-m_cand_3(1);
end

L_3=x-(BETA_i_3*(y-m_cand_3(2))+n_cand_3(2));
if BETA_i_3==inf
    L_3=y-m_cand_3(2);
end

L_4=x-(BETA_i_4*(y-m_cand_3(2))+n_cand_3(2));
if BETA_i_4==inf
    L_4=y-m_cand_3(2);
end

L_5=x-(BETA_i_5*(y-m_cand_3(3))+n_cand_3(3));
if BETA_i_5==inf
    L_5=y-m_cand_3(3);
end

L_6=x-(BETA_i_6*(y-m_cand_3(3))+n_cand_3(3));
if BETA_i_6==inf
    L_6=y-m_cand_3(3);
end

%% first sensor,first line connected to other lines to other sensor lines
A1=solve(L_1,L_3);
A2=solve(L_1,L_4);
A3=solve(L_1,L_5);
A4=solve(L_1,L_6);
% x of the first line of first sensor
X1=sym2poly(A1.x);
X2=sym2poly(A2.x);
X3=sym2poly(A3.x);
X4=sym2poly(A4.x);
% y of the first line of first sensor
y1=sym2poly(A1.y);
y2=sym2poly(A2.y);
y3=sym2poly(A3.y);
y4=sym2poly(A4.y);
%% first sensor, second line connected to other lines of other sensor lines
A5=solve(L_2,L_3);
A6=solve(L_2,L_4);
A7=solve(L_2,L_5);
A8=solve(L_2,L_6);
% x of second line of first sensor
X5=sym2poly(A5.x);
X6=sym2poly(A6.x);
X7=sym2poly(A7.x);
X8=sym2poly(A8.x);
% x of second line of first sensor
y5=sym2poly(A5.y);
y6=sym2poly(A6.y);
y7=sym2poly(A7.y);
y8=sym2poly(A8.y);

%% second sensor, first line connected to other lines of other sensors
A9=solve(L_3,L_1);
A10=solve(L_3,L_2);
A11=solve(L_3,L_5);
A12=solve(L_3,L_6);
% x of second sensor, first line
X9=sym2poly(A9.x);
X10=sym2poly(A10.x);
X11=sym2poly(A11.x);
X12=sym2poly(A12.x);
% y of second sensor, first line
y9=sym2poly(A9.y);
y10=sym2poly(A10.y);
y11=sym2poly(A11.y);
y12=sym2poly(A12.y);

%% second sensor, second line connected to other lines
A13=solve(L_4,L_1);
A14=solve(L_4,L_2);
A15=solve(L_4,L_5);
A16=solve(L_4,L_6);

% x of second sensor, second line
X13=sym2poly(A13.x);
X14=sym2poly(A14.x);
X15=sym2poly(A15.x);
X16=sym2poly(A16.x);

% y of second sensor,second sline
y13=sym2poly(A13.y);
y14=sym2poly(A14.y);
y15=sym2poly(A15.y);
y16=sym2poly(A16.y);

%% Third sensor, first line connected to other lines
A17=solve(L_5,L_1);
A18=solve(L_5,L_2);
A19=solve(L_5,L_3);
A20=solve(L_5,L_4);

% x of third sensor, first line
X17=sym2poly(A17.x);
X18=sym2poly(A18.x);
X19=sym2poly(A19.x);
X20=sym2poly(A20.x);

% y of third sensor, second line
y17=sym2poly(A17.y);
y18=sym2poly(A18.y);
y19=sym2poly(A19.y);
y20=sym2poly(A20.y);

%% third sensor, second line connected to other lines
A21=solve(L_6,L_1);
A22=solve(L_6,L_2);
A23=solve(L_6,L_3);
A24=solve(L_6,L_4);

% x of third sensor,second line
X21=sym2poly(A21.x);
X22=sym2poly(A22.x);
X23=sym2poly(A23.x);
X24=sym2poly(A24.x);

% y of third sensor, second line
y21=sym2poly(A21.y);
y22=sym2poly(A22.y);
y23=sym2poly(A23.y);
y24=sym2poly(A24.y);


%% seperating the X of each line of each sensor connected points

%% first sensor,first line
XX1=[X1 X2 X3 X4];
YY1=[y1 y2 y3 y4];

%% first sensor, second line
XX2=[X5 X6 X7 X8];
YY2=[y5 y6 y7 y8];

%% second sensor,first line
XX3=[X9 X10 X11 X12];
YY3=[y9 y10 y11 y12];

%% Second sensor,second line
XX4=[X13 X14 X15 X16];
YY4=[y13 y14 y15 y16];


%% Third sensor, first line
XX5=[X17 X18 X19 X20];
YY5=[y17 y18 y19 y20];

%% Third sensor, second line
XX6=[X21 X22 X23 X24];
YY6=[y21 y22 y23 y24];

%% the whole connected points 
XX=[X1 X2 X3 X4 X5 X6 X7 X8 X9 X10 X11 X12 X13 X14 X15 X16 X17 X18 X19 X20 X21 X22 X23 X24];
YY=[y1 y2 y3 y4 y5 y6 y7 y8 y9 y10 y11 y12 y13 y14 y15 y16 y17 y18 y19 y20 y21 y22 y23 y24 ];

%% computing distance for middle point selecting

for i=1:length(XX1)
 a1(i)=abs(XX1(i)-n_cand_3(1));
a_signx_1(i)=(XX1(i)-n_cand_3(1));
    
 a2(i)=abs(XX2(i)-n_cand_3(1));
 a_signx_2(i)=(XX2(i)-n_cand_3(1));
    
a3(i)=abs(XX3(i)-n_cand_3(2));
a_signx_3(i)=(XX3(i)-n_cand_3(2));

a4(i)=abs(XX4(i)-n_cand_3(2));
a_signx_4(i)=(XX4(i)-n_cand_3(2));

a5(i)=abs(XX5(i)-n_cand_3(3));
a_signx_5(i)=(XX5(i)-n_cand_3(3));

a6(i)=abs(XX6(i)-n_cand_3(3));
a_signx_6(i)=(XX6(i)-n_cand_3(3));

b1(i)=abs(YY1(i)-m_cand_3(1));
b_signy_1(i)=(YY1(i)-m_cand_3(1));

b2(i)=abs(YY2(i)-m_cand_3(1));
b_signy_2(i)=(YY2(i)-m_cand_3(1));

b3(i)=abs(YY3(i)-m_cand_3(2));
b_signy_3(i)=(YY3(i)-m_cand_3(2));

b4(i)=abs(YY4(i)-m_cand_3(2));
b_signy_4(i)=(YY4(i)-m_cand_3(2));

b5(i)=abs(YY5(i)-m_cand_3(3));
b_signy_5(i)=(YY5(i)-m_cand_3(3));

b6(i)=abs(YY6(i)-m_cand_3(3));
b_signy_6(i)=(YY6(i)-m_cand_3(3));

ab1(i)=(a1(i)^2+b1(i)^2)^(1/2);
ab2(i)=(a2(i)^2+b2(i)^2)^(1/2);
ab3(i)=(a3(i)^2+b3(i)^2)^(1/2);
ab4(i)=(a4(i)^2+b4(i)^2)^(1/2);
ab5(i)=(a5(i)^2+b5(i)^2)^(1/2);
ab6(i)=(a6(i)^2+b6(i)^2)^(1/2);

end

%% first sens, first line, noghte poshti 

p=0;
k=0;
    if (prod(a_signx_1)>0)
        XX1_c=XX1;
        YY1_c=YY1;
    end
    if (prod(a_signx_1)<0)
        for j=1:length(XX1)
            
if a_signx_1(j)>0
    p=p+1;
end
if a_signx_1(j)<0
    k=k+1;
end
        end
    end
    
    j=0;
    for i=1:length(XX1) 
  if p==3 && a_signx_1(i)>0
      j=j+1;
    XX1_c(j)=XX1(i);
    YY1_c(j)=YY1(i);
  end
    end
    
  j=0;
  for i=1:length(XX1)
      if k==3 && a_signx_1(i)<0
            j=j+1;
    XX1_c(j)=XX1(i);
    YY1_c(j)=YY1(i);
      end
  end
    


%%first sens,first line
for i=1:length(XX1_c)
    c1(i)=abs(XX1_c(i)-n_cand_3(1));
    d1(i)=abs(YY1_c(i)-m_cand_3(1));
    ab1(i)=(c1(i)^2+d1(i)^2)^(1/2);
end
k=0;
for i=1:length(XX1_c)
    c1(i)=abs(XX1_c(i)-n_cand_3(1));
    d1(i)=abs(YY1_c(i)-m_cand_3(1));
    cd1(i)=(c1(i)^2+d1(i)^2)^(1/2);
    if length(XX1_c)==4
    if cd1(i)<max(ab1) && cd1(i)>min(ab1)
        k=k+1;
    n1(k)=XX1_c(i);
    m1(k)=YY1_c(i);
    end
    end
    if length(XX1_c)==3
    if cd1(i)<max(ab1) 
        k=k+1;
    n1(k)=XX1_c(i);
    m1(k)=YY1_c(i);    
    end    
    end
    end

        
        
  
%% first sensor,second line
p=0;
k=0;
    if (prod(a_signx_2)>0)
        XX2_c=XX2;
        YY2_c=YY2;
    end
    if (prod(a_signx_2)<0)
     for i=1:length(XX2)       
if a_signx_2(i)>0
    p=p+1;
end
if a_signx_2(i)<0
    k=k+1;
end
end
        end
        

    j=0;
    for i=1:length(XX2) 
  if p==3 && a_signx_2(i)>0
      j=j+1;
    XX2_c(j)=XX2(i);
    YY2_c(j)=YY2(i);
  end
  
  j=0;
  for i=1:length(XX2)
      if k==3 && a_signx_2(i)<0
            j=j+1;
    XX2_c(j)=XX2(i);
    YY2_c(j)=YY2(i);
      end
  end

end
k=0;    
for i=1:length(XX2_c)
c2(i)=abs(XX2_c(i)-n_cand_3(1));
d2(i)=abs(YY2_c(i)-m_cand_3(1));
ab2(i)=(c2(i)^2+d2(i)^2)^(1/2);
end
for i=1:length(XX2_c)
    c2(i)=abs(XX2_c(i)-n_cand_3(1));
d2(i)=abs(YY2_c(i)-m_cand_3(1));
cd2(i)=(c2(i)^2+d2(i)^2)^(1/2);
    
if length(XX2_c)==4
if cd2(i)<max(ab2) && cd2(i)>min(ab2)
    k=k+1;
    n2(k)=XX2_c(i);
    m2(k)=YY2_c(i);
end
end
if length(XX2_c)==3
if cd2(i)<max(ab2) 
    k=k+1;
    n2(k)=XX2_c(i);
    m2(k)=YY2_c(i);
end    
end
end



p=0;
k=0;
    if (prod(a_signx_3)>0)
        XX3_c=XX3;
        YY3_c=YY3;
    end
    if (prod(a_signx_3)<0)
     for i=1:length(XX3)       
if a_signx_3(i)>0
    p=p+1;
end
if a_signx_3(i)<0
    k=k+1;
end
        end
        end

    j=0;
    for i=1:length(XX3) 
  if p==3 && a_signx_3(i)>0
      j=j+1;
    XX3_c(j)=XX3(i);
    YY3_c(j)=YY3(i);
  end
    end
  j=0;
  for i=1:length(XX3)
      if k==3 && a_signx_3(i)<0
            j=j+1;
    XX3_c(j)=XX3(i);
    YY3_c(j)=YY3(i);
      end
  end




k=0;
for i=1:length(XX3_c)
c3(i)=abs(XX3_c(i)-n_cand_3(2));
d3(i)=abs(YY3_c(i)-m_cand_3(2));
ab3(i)=(c3(i)^2+d3(i)^2)^(1/2);
end
for i=1:length(XX3_c)
c3(i)=abs(XX3_c(i)-n_cand_3(2));
d3(i)=abs(YY3_c(i)-m_cand_3(2));
cd3(i)=(c3(i)^2+d3(i)^2)^(1/2);
if length(XX3_c)==4
if cd3(i)<max(ab3) && cd3(i)>min(ab3)
    k=k+1;
n3(k)=XX3_c(i);
m3(k)=YY3_c(i);
end
end
if length(XX3_c)==3
if cd3(i)<max(ab3) 
    k=k+1;
n3(k)=XX3_c(i);
m3(k)=YY3_c(i) ;
end
end
end
    


p=0;
k=0;
    if (prod(a_signx_4)>0)
        XX4_c=XX4;
        YY4_c=YY4;
    end
    if (prod(a_signx_4)<0)
     for i=1:length(XX4)       
if a_signx_4(i)>0
    p=p+1;
end
if a_signx_4(i)<0
    k=k+1;
end
        end
        end

    j=0;
    for i=1:length(XX4) 
  if p==3 && a_signx_4(i)>0
      j=j+1;
    XX4_c(j)=XX4(i);
    YY4_c(j)=YY4(i);
  end
    end
  j=0;
  for i=1:length(XX4)
      if k==3 && a_signx_4(i)<0
            j=j+1;
    XX4_c(j)=XX4(i);
    YY4_c(j)=YY4(i);
      end
  end


k=0;

for i=1:length(XX4_c)
 c4(i)=abs(XX4_c(i)-n_cand_3(2));
 d4(i)=abs(YY4_c(i)-m_cand_3(2));
 ab4(i)=(c4(i)^2+d4(i)^2)^(1/2);
end
for i=1:length(XX4_c)
 c4(i)=abs(XX4_c(i)-n_cand_3(2));
 d4(i)=abs(YY4_c(i)-m_cand_3(2));
 cd4(i)=(c4(i)^2+d4(i)^2)^(1/2);
 if length(XX4_c)==4
   if cd4(i)<max(ab4) && cd4(i)>min(ab4)
       k=k+1;
       n4(k)=XX4_c(i);
       m4(k)=YY4_c(i);
   end
 end
 if length(XX4_c)==3
   if cd4(i)<max(ab4) 
   k=k+1;
       n4(k)=XX4_c(i);
       m4(k)=YY4_c(i);
   end
   end
 end  



 
p=0;
k=0;
    if (prod(a_signx_5)>0)
        XX5_c=XX5;
        YY5_c=YY5;
    end
    if (prod(a_signx_5)<0)
     for i=1:length(XX5)       
if a_signx_5(i)>0
    p=p+1;
end
if a_signx_5(i)<0
    k=k+1;
end
        end
        end

    j=0;
    for i=1:length(XX5) 
  if p==3 && a_signx_5(i)>0
      j=j+1;
    XX5_c(j)=XX5(i);
    YY5_c(j)=YY5(i);
  end
    end
  j=0;
  for i=1:length(XX5)
      if k==3 && a_signx_5(i)<0
            j=j+1;
    XX5_c(j)=XX5(i);
    YY5_c(j)=YY5(i);
      end
  end




for i=1:length(XX5_c)
    c5(i)=abs(XX5_c(i)-n_cand_3(3));
d5(i)=abs(YY5_c(i)-m_cand_3(3));
ab5(i)=(c5(i)^2+d5(i)^2)^(1/2);
end
k=0;
for i=1:length(XX5_c)
c5(i)=abs(XX5_c(i)-n_cand_3(3));
d5(i)=abs(YY5_c(i)-m_cand_3(3));
cd5(i)=(c5(i)^2+d5(i)^2)^(1/2);
if length(XX5_c)==4
   if cd5(i)<max(ab5) && cd5(i)>min(ab5)
       k=k+1;
       n5(k)=XX5(i);
       m5(k)=YY5(i);
   end
   end
   if length(XX5_c)==3
       if cd5(i)<max(ab5) 
       k=k+1;
       n5(k)=XX5(i);
       m5(k)=YY5(i);
       
   end
end
end
    


p=0;
k=0;
    if (prod(a_signx_6)>0)
        XX6_c=XX6;
        YY6_c=YY6;
    end
    if (prod(a_signx_6)<0)
     for i=1:length(XX6)       
if a_signx_6(i)>0
    p=p+1;
end
if a_signx_6(i)<0
    k=k+1;
end
        end
        end

    j=0;
    for i=1:length(XX6) 
  if p==3 && a_signx_6(i)>0
      j=j+1;
    XX6_c(j)=XX6(i);
    YY6_c(j)=YY6(i);
  end
    end
  j=0;
  for i=1:length(XX6)
      if k==3 && a_signx_6(i)<0
            j=j+1;
    XX6_c(j)=XX6(i);
    YY6_c(j)=YY6(i);
      end
  end





k=0;
for i=1:length(XX6_c)
c6(i)=abs(XX6_c(i)-n_cand_3(3));
d6(i)=abs(YY6_c(i)-m_cand_3(3));
ab6(i)=(c6(i)^2+d6(i)^2)^(1/2);
end
for i=1:length(XX6_c)
c6(i)=abs(XX6_c(i)-n_cand_3(3));
d6(i)=abs(YY6_c(i)-m_cand_3(3));
cd6(i)=(c6(i)^2+d6(i)^2)^(1/2);
if length(XX6_c)==4
       if cd6(i)<max(ab6) && cd6(i)>min(ab6)
           k=k+1;
           n6(k)=XX6(i);
           m6(k)=YY6(i);
       end
end
if length(XX6_c)==3       
       if cd6(i)<max(ab6) 
           k=k+1;
           n6(k)=XX6(i);
           m6(k)=YY6(i);       
     end
end

end




%% selected nodes after calculating distance 
N=[n1 n2 n3 n4 n5 n6];
M=[m1 m2 m3 m4 m5 m6];
% NN=N(1:6)
% MM=M(1:6);

% % seprating nodes to avoid choosing same nodes
 k=0;
 p=0;
for i=1:length(N)
    for j=i+1:length(N)
    if (N(i)==N(j))&&(M(i)==M(j)) 
        k=k+1;
NN(k)=N(i);
MM(k)=M(i);
    end
    end
    end    
    
  NN;
  MM;
% n_v = length(N);
% i=2;
% j=1;
% while j<n_v
%     while i<=n_v
% 
%         if N(j)==N(i)
%             NN = [N(1:i-1) N(i+1:n_v)];
%             n_v = n_v-1;
%             i = i-1;
%         end
%         i=i+1;
%     end
%     j=j+1;
%     i = j+1;
% end
% 
% i=2;
% j=1;
% n_v=length(M);
% while j<n_v
%     while i<=n_v
% 
%         if N(j)==N(i)
%             MM = [M(1:i-1) M(i+1:n_v)];
%             n_v = n_v-1;
%             i = i-1;
%         end
%         i=i+1;
%     end
%     j=j+1;
%     i = j+1;
% end
% NN



if length(NN)==6;
    
for i=1:length(NN)
    for j=1:length(NN)
        R(i,j)=(NN(i)-NN(j))/(MM(i)-MM(j));
    L(i,j)=x-(R(i,j)*(y-MM(j))+NN(j));
    if L(i,j)==inf
        L(i,j)=y-NN(j);
    end
    end
end
%% first line and calculating surface
%substituting nodes in equations
F1=subs(L(1,2),{x,y},{NN(3),MM(3)});
F2=subs(L(1,2),{x,y},{NN(4),MM(4)});
F3=subs(L(1,2),{x,y},{NN(5),MM(5)});
F4=subs(L(1,2),{x,y},{NN(6),MM(6)});
nn=[NN(3) NN(4) NN(5) NN(6)];
mm=[MM(3) MM(4) MM(5) MM(6)];
f1=[F1 F2 F3 F4 ];

k=0;
j=0;
for i=1:length(f1)
    if f1(i)>0
        k=k+1;
        NN_U(k)=nn(i);
        MM_U(k)=mm(i);
    end
        if f1(i)<0
j=j+1;
NN_D(j)=nn(i);
MM_D(j)=mm(i);
    end
end
if k==2 && j==2
    N_U=[NN_U NN(1) NN(2)];
    M_U=[MM_U MM(1) MM5(2)];
    N_D=[NN_D NN(1) NN(2)];
    M_D=[MM_D MM(1) MM(2)];
    
    % up surface
for i=1:length(N_U)
    for j=1:length(N_U)
        RR(i,j)=(N_U(i)-N_U(j))/(M_U(i)-M_U(j));
    LL(i,j)=x-(RR(i,j)*(y-M_U(j))+N_U(j));
    if LL(i,j)==inf
        LL(i,j)=y-N_U(j);
    end
    end
end
    FF1=subs(LL(1,2),{x,y},{N_U(3),M_U(3)});
    FF2=subs(LL(1,2),{x,y},{N_U(4),M_U(4)});
    
    FF3=subs(LL(1,3),{x,y},{N_U(2),M_U(2)});
    FF4=subs(LL(1,3),{x,y},{N_U(4),M_U(4)});
    
    FF5=subs(LL(1,4),{x,y},{N_U(2),M_U(2)});
    FF6=subs(LL(1,4),{x,y},{N_U(3),M_U(3)});
    
if FF1*FF2<0
    X_side=[N_U(3) N_U(4)];
    Y_side=[M_U(3) M_U(4)];
    X_front=[N_U(1) N_U(2)];
    Y_front=[M_U(1) M_U(2)];
    R_M=RR(1,2);
    [SS]=CS(R_M,X_side,Y_side,X_front,Y_front);
    SS_1=SS;
end    

if FF3*FF4<0
    X_side=[N_U(2) N_U(4)];
    Y_side=[M_U(2) M_U(4)];
    
    X_front=[N_U(1) N_U(3)];
    Y_front=[M_U(1) M_U(3)];
    R_M=RR(1,3);
    [SS]=CS(R_M,X_side,Y_side,X_front,Y_front);
     SS_1=SS;
end
    

if FF5*FF6<0
    X_side=[N_U(2) N_U(3)];
    Y_side=[M_U(2) M_U(3)];
    
    X_front=[N_U(1) N_U(4)];
    Y_front=[M_U(1) M_U(4)];
    R_M=RR(1,4);
    [SS]=CS(R_M,X_side,Y_side,X_front,Y_front);
     SS_1=SS;
end









% down surface

for i=1:length(N_D)
    for j=1:length(N_D)
        RR(i,j)=(N_D(i)-N_D(j))/(M_D(i)-M_D(j));
    LL(i,j)=x-(RR(i,j)*(y-M_D(j))+N_D(j));
    if LL(i,j)==inf
        LL(i,j)=y-N_D(j);
    end
    end
end
    FF1=subs(LL(1,2),{x,y},{N_D(3),M_D(3)});
    FF2=subs(LL(1,2),{x,y},{N_D(4),M_D(4)});
    
    FF3=subs(LL(1,3),{x,y},{N_D(2),M_D(2)});
    FF4=subs(LL(1,3),{x,y},{N_D(4),M_D(4)});
    
    FF5=subs(LL(1,4),{x,y},{N_D(2),M_D(2)});
    FF6=subs(LL(1,4),{x,y},{N_D(3),M_D(3)});
    
if FF1*FF2<0
    X_side=[N_D(3) N_D(4)];
    Y_side=[M_D(3) M_D(4)];
    X_front=[N_D(1) N_D(2)];
    Y_front=[M_D(1) M_D(2)];
    R_M=RR(1,2);
    [SS]=CS(R_M,X_side,Y_side,X_front,Y_front);
    SS_2=SS;
end    

if FF3*FF4<0
    X_side=[N_D(2) N_D(4)];
    Y_side=[M_D(2) M_D(4)];
    
    X_front=[N_D(1) N_D(3)];
    Y_front=[M_D(1) M_D(3)];
    R_M=RR(1,3);
    [SS]=CS(R_M,X_side,Y_side,X_front,Y_front);
     SS_2=SS;
end
    

if FF5*FF6<0
    X_side=[N_D(2) N_D(3)];
    Y_side=[M_D(2) M_D(3)];
    
    X_front=[N_D(1) N_D(4)];
    Y_front=[M_D(1) M_D(4)];
    R_M=RR(1,4);
    [SS]=CS(R_M,X_side,Y_side,X_front,Y_front);
     SS_2=SS;
end
%% final surface calculating
S=SS_1+SS_2;
end














%% second line... calculating surface 


F5=subs(L(1,3),{x,y},{NN(2),MM(2)});
F6=subs(L(1,3),{x,y},{NN(4),MM(4)});
F7=subs(L(1,3),{x,y},{NN(5),MM(5)});
F8=subs(L(1,3),{x,y},{NN(6),MM(6)});
nn=[NN(2) NN(4) NN(5) NN(6)];
mm=[MM(2) MM(4) MM(5) MM(6)];
f2=[F5 F6 F7 F8];




k=0;
j=0;
for i=1:length(f1)
    if f2(i)>0
        k=k+1;
        NN_U(k)=nn(i);
        MM_U(k)=mm(i);
    end
        if f2(i)<0
j=j+1;
NN_D(j)=nn(i);
        MM_D(j)=mm(i);
    end
end
if k==2 && j==2
    N_U=[NN_U NN(1) NN(3)];
    M_U=[MM_U MM(1) MM(3)];
    N_D=[NN_D NN(1) NN(3)];
    M_D=[MM_D MM(1) MM(3)];
    
    %% upper surface
for i=1:length(N_U)
    for j=1:length(N_U)
        RR(i,j)=(N_U(i)-N_U(j))/(M_U(i)-M_U(j));
    LL(i,j)=x-(RR(i,j)*(y-M_U(j))+N_U(j));
    if LL(i,j)==inf
        LL(i,j)=y-N_U(j);
    end
    end
end
    FF1=subs(LL(1,2),{x,y},{N_U(3),M_U(3)});
    FF2=subs(LL(1,2),{x,y},{N_U(4),M_U(4)});
    
    FF3=subs(LL(1,3),{x,y},{N_U(2),M_U(2)});
    FF4=subs(LL(1,3),{x,y},{N_U(4),M_U(4)});
    
    FF5=subs(LL(1,4),{x,y},{N_U(2),M_U(2)});
    FF6=subs(LL(1,4),{x,y},{N_U(3),M_U(3)});
    
if FF1*FF2<0
    X_side=[N_U(3) N_U(4)];
    Y_side=[M_U(3) M_U(4)];
    X_front=[N_U(1) N_U(2)];
    Y_front=[M_U(1) M_U(2)];
    R_M=RR(1,2);
    [SS]=CS(R_M,X_side,Y_side,X_front,Y_front);
    SS_1=SS;
end    

if FF3*FF4<0
    X_side=[N_U(2) N_U(4)];
    Y_side=[M_U(2) M_U(4)];
    
    X_front=[N_U(1) N_U(3)];
    Y_front=[M_U(1) M_U(3)];
    R_M=RR(1,3);
    [SS]=CS(R_M,X_side,Y_side,X_front,Y_front);
     SS_1=SS;
end
    

if FF5*FF6<0
    X_side=[N_U(2) N_U(3)];
    Y_side=[M_U(2) M_U(3)];
    
    X_front=[N_U(1) N_U(4)];
    Y_front=[M_U(1) M_U(4)];
    R_M=RR(1,4);
    [SS]=CS(R_M,X_side,Y_side,X_front,Y_front);
     SS_1=SS;
end










%% downer surface
for i=1:length(N_D)
    for j=1:length(N_D)
        RR(i,j)=(N_D(i)-N_D(j))/(M_D(i)-M_D(j));
    LL(i,j)=x-(RR(i,j)*(y-M_D(j))+N_D(j));
    if LL(i,j)==inf
        LL(i,j)=y-N_D(j);
    end
    end
end
    FF1=subs(LL(1,2),{x,y},{N_D(3),M_D(3)});
    FF2=subs(LL(1,2),{x,y},{N_D(4),M_D(4)});
    
    FF3=subs(LL(1,3),{x,y},{N_D(2),M_D(2)});
    FF4=subs(LL(1,3),{x,y},{N_D(4),M_D(4)});
    
    FF5=subs(LL(1,4),{x,y},{N_D(2),M_D(2)});
    FF6=subs(LL(1,4),{x,y},{N_D(3),M_D(3)});
    
if FF1*FF2<0
    X_side=[N_D(3) N_D(4)];
    Y_side=[M_D(3) M_D(4)];
    X_front=[N_D(1) N_D(2)];
    Y_front=[M_D(1) M_D(2)];
    R_M=RR(1,2);
    [SS]=CS(R_M,X_side,Y_side,X_front,Y_front);
    SS_2=SS;
end    

if FF3*FF4<0
    X_side=[N_D(2) N_D(4)];
    Y_side=[M_D(2) M_D(4)];
    
    X_front=[N_D(1) N_D(3)];
    Y_front=[M_D(1) M_D(3)];
    R_M=RR(1,3);
    [SS]=CS(R_M,X_side,Y_side,X_front,Y_front);
     SS_2=SS;
end
    

if FF5*FF6<0
    X_side=[N_D(2) N_D(3)];
    Y_side=[M_D(2) M_D(3)];
    
    X_front=[N_D(1) N_D(4)];
    Y_front=[M_D(1) M_D(4)];
    R_M=RR(1,4);
    [SS]=CS(R_M,X_side,Y_side,X_front,Y_front);
     SS_2=SS;
end
%% final surface calculating
S=SS_1+SS_2;
end



%% third line... calculating surface 

F9=subs(L(1,4),{x,y},{n2,m2});
F10=subs(L(1,4),{x,y},{n3,m3});
F11=subs(L(1,4),{x,y},{n5,m5});
F12=subs(L(1,4),{x,y},{n6,m6});
f3=[F9 F10 F11 F12];

nn=[NN(2) NN(3) NN(5) NN(6)];
mm=[MM(2) MM(3) MM(5) MM(6)];

k=0;
j=0;
for i=1:length(f1)
    if f3(i)>0
        k=k+1;
        NN_U(k)=nn(i);
        MM_U(k)=mm(i);
    end
        if f3(i)<0
j=j+1;
NN_D(j)=nn(i);
        MM_D(j)=mm(i);
    end
end
if k==2 && j==2
    N_U=[NN_U NN(1) NN(4)];
    M_U=[MM_U MM(1) MM(4)];
    N_D=[NN_D NN(1) NN(4)];
    M_D=[MM_D MM(1) MM(4)];
    
    %% upper surface
for i=1:length(N_U)
    for j=1:length(N_U)
        RR(i,j)=(N_U(i)-N_U(j))/(M_U(i)-M_U(j));
    LL(i,j)=x-(RR(i,j)*(y-M_U(j))+N_U(j));
    if LL(i,j)==inf
        LL(i,j)=y-N_U(j);
    end
    end
end
    FF1=subs(LL(1,2),{x,y},{N_U(3),M_U(3)});
    FF2=subs(LL(1,2),{x,y},{N_U(4),M_U(4)});
    
    FF3=subs(LL(1,3),{x,y},{N_U(2),M_U(2)});
    FF4=subs(LL(1,3),{x,y},{N_U(4),M_U(4)});
    
    FF5=subs(LL(1,4),{x,y},{N_U(2),M_U(2)});
    FF6=subs(LL(1,4),{x,y},{N_U(3),M_U(3)});
    
if FF1*FF2<0
    X_side=[N_U(3) N_U(4)];
    Y_side=[M_U(3) M_U(4)];
    X_front=[N_U(1) N_U(2)];
    Y_front=[M_U(1) M_U(2)];
    R_M=RR(1,2);
    [SS]=CS(R_M,X_side,Y_side,X_front,Y_front);
    SS_1=SS;
end    

if FF3*FF4<0
    X_side=[N_U(2) N_U(4)];
    Y_side=[M_U(2) M_U(4)];
    
    X_front=[N_U(1) N_U(3)];
    Y_front=[M_U(1) M_U(3)];
    R_M=RR(1,3);
    [SS]=CS(R_M,X_side,Y_side,X_front,Y_front);
     SS_1=SS;
end
    

if FF5*FF6<0
    X_side=[N_U(2) N_U(3)];
    Y_side=[M_U(2) M_U(3)];
    
    X_front=[N_U(1) N_U(4)];
    Y_front=[M_U(1) M_U(4)];
    R_M=RR(1,4);
    [SS]=CS(R_M,X_side,Y_side,X_front,Y_front);
     SS_1=SS;
end




%% downer surface
for i=1:length(N_D)
    for j=1:length(N_D)
        RR(i,j)=(N_D(i)-N_D(j))/(M_D(i)-M_D(j));
    LL(i,j)=x-(RR(i,j)*(y-M_D(j))+N_D(j));
    if LL(i,j)==inf
        LL(i,j)=y-N_D(j);
    end
    end
end
    FF1=subs(LL(1,2),{x,y},{N_D(3),M_D(3)});
    FF2=subs(LL(1,2),{x,y},{N_D(4),M_D(4)});
    
    FF3=subs(LL(1,3),{x,y},{N_D(2),M_D(2)});
    FF4=subs(LL(1,3),{x,y},{N_D(4),M_D(4)});
    
    FF5=subs(LL(1,4),{x,y},{N_D(2),M_D(2)});
    FF6=subs(LL(1,4),{x,y},{N_D(3),M_D(3)});
    
if FF1*FF2<0
    X_side=[N_D(3) N_D(4)];
    Y_side=[M_D(3) M_D(4)];
    X_front=[N_D(1) N_D(2)];
    Y_front=[M_D(1) M_D(2)];
    R_M=RR(1,2);
    [SS]=CS(R_M,X_side,Y_side,X_front,Y_front);
    SS_2=SS;
end    

if FF3*FF4<0
    X_side=[N_D(2) N_D(4)];
    Y_side=[M_D(2) M_D(4)];
    
    X_front=[N_D(1) N_D(3)];
    Y_front=[M_D(1) M_D(3)];
    R_M=RR(1,3);
    [SS]=CS(R_M,X_side,Y_side,X_front,Y_front);
     SS_2=SS;
end
    
if FF5*FF6<0
    X_side=[N_D(2) N_D(3)];
    Y_side=[M_D(2) M_D(3)];
    
    X_front=[N_D(1) N_D(4)];
    Y_front=[M_D(1) M_D(4)];
    R_M=RR(1,4);
    [SS]=CS(R_M,X_side,Y_side,X_front,Y_front);
     SS_2=SS;
end
%% final surface calculating
S=SS_1+SS_2;
end




%% forth line... calculating surface 

F13=subs(L(1,5),{x,y},{n2,m2});
F14=subs(L(1,5),{x,y},{n3,m3});
F15=subs(L(1,5),{x,y},{n4,m4});
F16=subs(L(1,5),{x,y},{n6,m6});
f4=[F13 F14 F15 F16];
nn=[NN(2) NN(3) NN(4) NN(6)];
mm=[MM(2) MM(3) MM(4) MM(6)];



k=0;
j=0;
for i=1:length(f1)
    if f4(i)>0
        k=k+1;
        NN_U(k)=nn(i);
        MM_U(k)=mm(i);
    end
        if f4(i)<0
j=j+1;
NN_D(j)=nn(i);
        MM_D(j)=mm(i);
    end
end
if k==2 && j==2
    N_U=[NN_U NN(1) NN(5)];
    M_U=[MM_U MM(1) MM(5)];
    N_D=[NN_D NN(1) NN(5)];
    M_D=[MM_D MM(1) MM(5)];
    
    %% upper surface
for i=1:length(N_U)
    for j=1:length(N_U)
        RR(i,j)=(N_U(i)-N_U(j))/(M_U(i)-M_U(j));
    LL(i,j)=x-(RR(i,j)*(y-M_U(j))+N_U(j));
    if LL(i,j)==inf
        LL(i,j)=y-N_U(j);
    end
    end
end
    FF1=subs(LL(1,2),{x,y},{N_U(3),M_U(3)});
    FF2=subs(LL(1,2),{x,y},{N_U(4),M_U(4)});
    
    FF3=subs(LL(1,3),{x,y},{N_U(2),M_U(2)});
    FF4=subs(LL(1,3),{x,y},{N_U(4),M_U(4)});
    
    FF5=subs(LL(1,4),{x,y},{N_U(2),M_U(2)});
    FF6=subs(LL(1,4),{x,y},{N_U(3),M_U(3)});
    
if FF1*FF2<0
    X_side=[N_U(3) N_U(4)];
    Y_side=[M_U(3) M_U(4)];
    X_front=[N_U(1) N_U(2)];
    Y_front=[M_U(1) M_U(2)];
    R_M=RR(1,2);
    [SS]=CS(R_M,X_side,Y_side,X_front,Y_front);
    SS_1=SS;
end    

if FF3*FF4<0
    X_side=[N_U(2) N_U(4)];
    Y_side=[M_U(2) M_U(4)];
    
    X_front=[N_U(1) N_U(3)];
    Y_front=[M_U(1) M_U(3)];
    R_M=RR(1,3);
    [SS]=CS(R_M,X_side,Y_side,X_front,Y_front);
     SS_1=SS;
end
    

if FF5*FF6<0
    X_side=[N_U(2) N_U(3)];
    Y_side=[M_U(2) M_U(3)];
    
    X_front=[N_U(1) N_U(4)];
    Y_front=[M_U(1) M_U(4)];
    R_M=RR(1,4);
    [SS]=CS(R_M,X_side,Y_side,X_front,Y_front);
     SS_1=SS;
end


%% downer surface
for i=1:length(N_D)
    for j=1:length(N_D)
        RR(i,j)=(N_D(i)-N_D(j))/(M_D(i)-M_D(j));
    LL(i,j)=x-(RR(i,j)*(y-M_D(j))+N_D(j));
    if LL(i,j)==inf
        LL(i,j)=y-N_D(j);
    end
    end
end
    FF1=subs(LL(1,2),{x,y},{N_D(3),M_D(3)});
    FF2=subs(LL(1,2),{x,y},{N_D(4),M_D(4)});
    
    FF3=subs(LL(1,3),{x,y},{N_D(2),M_D(2)});
    FF4=subs(LL(1,3),{x,y},{N_D(4),M_D(4)});
    
    FF5=subs(LL(1,4),{x,y},{N_D(2),M_D(2)});
    FF6=subs(LL(1,4),{x,y},{N_D(3),M_D(3)});
    
if FF1*FF2<0
    X_side=[N_D(3) N_D(4)];
    Y_side=[M_D(3) M_D(4)];
    X_front=[N_D(1) N_D(2)];
    Y_front=[M_D(1) M_D(2)];
    R_M=RR(1,2);
    [SS]=CS(R_M,X_side,Y_side,X_front,Y_front);
    SS_2=SS;
end    

if FF3*FF4<0
    X_side=[N_D(2) N_D(4)];
    Y_side=[M_D(2) M_D(4)];
    X_front=[N_D(1) N_D(3)];
    Y_front=[M_D(1) M_D(3)];
    R_M=RR(1,3);
    [SS]=CS(R_M,X_side,Y_side,X_front,Y_front);
     SS_2=SS;
end
    

if FF5*FF6<0
    X_side=[N_D(2) N_D(3)];
    Y_side=[M_D(2) M_D(3)];
    X_front=[N_D(1) N_D(4)];
    Y_front=[M_D(1) M_D(4)];
    R_M=RR(1,4);
    [SS]=CS(R_M,X_side,Y_side,X_front,Y_front);
     SS_2=SS;
end
%% final surface calculating
S=SS_1+SS_2;
end




%% fifth line... calculating surface 

F17=subs(L(1,6),{x,y},{n2,m2});
F18=subs(L(1,6),{x,y},{n3,m3});
F19=subs(L(1,6),{x,y},{n4,m4});
F20=subs(L(1,6),{x,y},{n5,m5});
nn=[NN(2) NN(3) NN(4) NN(5)];
mm=[MM(2) MM(3) MM(4) MM(5)];
f5=[F17 F18 F19 F20];




k=0;
j=0;
for i=1:length(f1)
    if f4(i)>0
        k=k+1;
        NN_U(k)=nn(i);
        MM_U(k)=mm(i);
    end
        if f4(i)<0
j=j+1;
        NN_D(j)=nn(i);
        MM_D(j)=mm(i);
    end
end
if k==2 && j==2
    N_U=[NN_U NN(1) NN(6)];
    M_U=[MM_U MM(1) MM(6)];
    N_D=[NN_D NN(1) NN(6)];
    M_D=[MM_D MM(1) MM(6)];
    
    %% upper surface
for i=1:length(N_U)
    for j=1:length(N_U)
        RR(i,j)=(N_U(i)-N_U(j))/(M_U(i)-M_U(j));
    LL(i,j)=x-(RR(i,j)*(y-M_U(j))+N_U(j));
    if LL(i,j)==inf
        LL(i,j)=y-N_U(j);
    end
    end
end
    FF1=subs(LL(1,2),{x,y},{N_U(3),M_U(3)});
    FF2=subs(LL(1,2),{x,y},{N_U(4),M_U(4)});
    
    FF3=subs(LL(1,3),{x,y},{N_U(2),M_U(2)});
    FF4=subs(LL(1,3),{x,y},{N_U(4),M_U(4)});
    
    FF5=subs(LL(1,4),{x,y},{N_U(2),M_U(2)});
    FF6=subs(LL(1,4),{x,y},{N_U(3),M_U(3)});
    
if FF1*FF2<0
    X_side=[N_U(3) N_U(4)];
    Y_side=[M_U(3) M_U(4)];
    X_front=[N_U(1) N_U(2)];
    Y_front=[M_U(1) M_U(2)];
    R_M=RR(1,2);
    [SS]=CS(R_M,X_side,Y_side,X_front,Y_front);
    SS_1=SS;
end    

if FF3*FF4<0
    X_side=[N_U(2) N_U(4)];
    Y_side=[M_U(2) M_U(4)];
    
    X_front=[N_U(1) N_U(3)];
    Y_front=[M_U(1) M_U(3)];
    R_M=RR(1,3);
    [SS]=CS(R_M,X_side,Y_side,X_front,Y_front);
     SS_1=SS;
end
    

if FF5*FF6<0
    X_side=[N_U(2) N_U(3)];
    Y_side=[M_U(2) M_U(3)];
    
    X_front=[N_U(1) N_U(4)];
    Y_front=[M_U(1) M_U(4)];
    R_M=RR(1,4);
    [SS]=CS(R_M,X_side,Y_side,X_front,Y_front);
     SS_1=SS;
end




%% downer surface
for i=1:length(N_D)
    for j=1:length(N_D)
        RR(i,j)=(N_D(i)-N_D(j))/(M_D(i)-M_D(j));
    LL(i,j)=x-(RR(i,j)*(y-M_D(j))+N_D(j));
    if LL(i,j)==inf
        LL(i,j)=y-N_D(j);
    end
    end
end
    FF1=subs(LL(1,2),{x,y},{N_D(3),M_D(3)});
    FF2=subs(LL(1,2),{x,y},{N_D(4),M_D(4)});
    
    FF3=subs(LL(1,3),{x,y},{N_D(2),M_D(2)});
    FF4=subs(LL(1,3),{x,y},{N_D(4),M_D(4)});
    
    FF5=subs(LL(1,4),{x,y},{N_D(2),M_D(2)});
    FF6=subs(LL(1,4),{x,y},{N_D(3),M_D(3)});
    
if FF1*FF2<0
    X_side=[N_D(3) N_D(4)];
    Y_side=[M_D(3) M_D(4)];
    X_front=[N_D(1) N_D(2)];
    Y_front=[M_D(1) M_D(2)];
    R_M=RR(1,2);
    [SS]=CS(R_M,X_side,Y_side,X_front,Y_front);
    SS_2=SS;
end    

if FF3*FF4<0
    X_side=[N_D(2) N_D(4)];
    Y_side=[M_D(2) M_D(4)];
    X_front=[N_D(1) N_D(3)];
    Y_front=[M_D(1) M_D(3)];
    R_M=RR(1,3);
    [SS]=CS(R_M,X_side,Y_side,X_front,Y_front);
     SS_2=SS;
end
    

if FF5*FF6<0
    X_side=[N_D(2) N_D(3)];
    Y_side=[M_D(2) M_D(3)];
    X_front=[N_D(1) N_D(4)];
    Y_front=[M_D(1) M_D(4)];
    R_M=RR(1,4);
    [SS]=CS(R_M,X_side,Y_side,X_front,Y_front);
     SS_2=SS;
end
%% final surface calculating
S=SS_1+SS_2;

end
elseif length(NN)<6
    S=0;
end
end



% 
% for i=1:length(NN)
% k=k+1;
% if k<length(NN) || k==length(NN) 
%     R(i)=(NN(i)-NN(k))/(MM(i)-MM(k));
% end
% end


        
        
        

















































