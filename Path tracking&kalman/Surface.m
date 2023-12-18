%% 2 0ut of n
% x1 is x
% x2 is y

function S_t=Surface(X_i,y_path_i,n_cand,m_cand,SIGMA_i)
syms x1 x2
BETA_i_1=atan((n_cand(1)-X_i)/(m_cand(1)-y_path_i))+SIGMA_i(1);
BETA_i_2=atan((n_cand(1)-X_i)/(m_cand(1)-y_path_i))-SIGMA_i(1);
BETA_i_3=atan((n_cand(2)-X_i)/(m_cand(2)-y_path_i))+SIGMA_i(2);
BETA_i_4=atan((n_cand(2)-X_i)/(m_cand(2)-y_path_i))-SIGMA_i(2);

BETA_i_1=tan(BETA_i_1);
BETA_i_2=tan(BETA_i_2);
BETA_i_3=tan(BETA_i_3);
BETA_i_4=tan(BETA_i_4);


L_1=x1-(BETA_i_1*(x2-m_cand(1))+n_cand(1));
L_2=x1-(BETA_i_2*(x2-m_cand(1))+n_cand(1));
L_3=x1-(BETA_i_3*(x2-m_cand(2))+n_cand(2));
L_4=x1-(BETA_i_4*(x2-m_cand(2))+n_cand(2));

A1=solve(L_1,L_3);
A2=solve(L_1,L_4);
A3=solve(L_2,L_3);
A4=solve(L_2,L_4);
X1=sym2poly(A1.x1);
X2=sym2poly(A2.x1);
X3=sym2poly(A3.x1);
X4=sym2poly(A4.x1);

y1=sym2poly(A1.x2);
y2=sym2poly(A2.x2);
y3=sym2poly(A3.x2);
y4=sym2poly(A4.x2);
X=[X1 X2 X3 X4];
Y=[y1 y2 y3 y4];
for i=1:length(X)
    for j=1:length(X) 
    R(i,j)=(X(i)-X(j))/(Y(i)-Y(j));
    L(i,j)=x1-(R(i,j)*(x2-Y(j))+X(j));
    if L(i,j)==inf
        L(i,j)=x2-X(j);
    end
    end
end
%substitiung nodes in line equation
F1=subs(L(1,2),{x1,x2},{X3,y3});
F2=subs(L(1,2),{x1,x2},{X4,y4});

F3=subs(L(1,3),{x1,x2},{X2,y2});
F4=subs(L(1,3),{x1,x2},{X4,y4});

F5=subs(L(1,4),{x1,x2},{X2,y2});
F6=subs(L(1,4),{x1,x2},{X3,y3});

if F1*F2<0
    X_side=[X3 X4];
    Y_side=[y3 y4];
    X_front=[X1 X2];
    Y_front=[y1 y2];
    R_M=R(1,2);
    [SS]=CS(R_M,X_side,Y_side,X_front,Y_front);
elseif F3*F4<0
    X_side=[X2 X4];
    Y_side=[y2 y4];
    X_front=[X1 X3];
    Y_front=[y1 y3];
    R_M=R(1,3);
    [SS]=CS(R_M,X_side,Y_side,X_front,Y_front);

elseif F5*F6<0
    X_side=[X2 X3];
    Y_side=[y2 y3];
    X_front=[X1 X4];
    Y_front=[y1 y4];
    R_M=R(1,4);
    [SS]=CS(R_M,X_side,Y_side,X_front,Y_front);
    
end
%% total surface
S_t=SS;
    
    









% end
% k=1;
% j=0;
% for i=1:length(X)
%     k=k+1;
%         if k<length(X) || k==length(X)
%         r(i)=(X(i)-X(k))/(Y(i)-Y(k));
%         if r(i)==max(R) || r(i)==min(R)
%             j=j+1;
%             XX(j)=X(i);
%             YY(j)=Y(i);
%         end
%         end
% end
% 
% 
% 
% 
% k=1;
% for i=1:length(X)
%         j=0;
%         k=k+1;
%         if k<length(X) || k==length(X)
%         r(i)=(X(i)-X(k))/(Y(i)-Y(k));
%         
%     if r(i)<max(R) && r(i)>min(R)
%         j=j+1;
%         l_1=x2-(r(i)*(x1-X(i))+Y(i));
%         l_2=x2-((-1/r(i))*(x1-XX(j))+YY(j));
%         l_3=x2-((-1/r(i))*(x1-XX(j+1))+YY(j+1));
%         AA1=solve(l_1,l_2);
%         AA2=solve(l_1,l_3);
%         xx2=sym2poly(AA1.x1);
%         xx3=sym2poly(AA2.x1);
%         yy2=sym2poly(AA1.x2);
%         yy3=sym2poly(AA2.x2);
%         dis_X_ghaede=abs(X(i)-X(k));
%         dis_Y_ghaede=abs(Y(i)-Y(k));
%         dis_X_Y_ghaede=((dis_X_ghaede)^2+(dis_Y_ghaede^2))^(1/2);
%         dis_x_1_ertefa=abs(XX(j)-xx2);
%         dis_x_2_ertefa=abs(XX(j+1)-xx3);
%         dis_y_1_ertefa=abs(YY(j)-y2);
%         dis_y_2_ertefa=abs(YY(j+1)-y3);
%         dis_xy_1_ertefa=((dis_x_1_ertefa)^2+(dis_y_1_ertefa)^2)^(1/2);
%         dis_xy_2_ertefa=((dis_x_2_ertefa)^2+(dis_y_2_ertefa)^2)^(1/2);
%         S=0.5*(dis_X_Y_ghaede*dis_xy_1_ertefa)+0.5*(dis_X_Y_ghaede*dis_xy_2_ertefa);
%     
%     end
%         end
% end

    
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        

















