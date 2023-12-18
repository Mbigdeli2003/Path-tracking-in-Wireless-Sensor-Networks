function [SS]=CS(R_M,X_side,Y_side,X_front,Y_front)
    syms x y
    %% lines equations
    L_M=x-(R_M*(y-Y_front(1))+X_front(1));
    L_S_1=x-(-1/R_M*(y-Y_side(1))+X_side(1));
    L_S_2=x-(-1/R_M*(y-Y_side(2))+X_side(2));
    %% connected points
    A1=solve(L_M,L_S_1);
    A2=solve(L_M,L_S_2);
    
    X1=sym2poly(A1.x);
    X2=sym2poly(A2.x);
    
    Y1=sym2poly(A1.y);
    Y2=sym2poly(A2.y);
    %%ghaede
    g_x=abs(X_front(1)-X_front(2));
    g_y=abs(Y_front(1)-Y_front(2));
    G=(g_x^2+g_y^2)^(1/2);
    
    %%ertefa1
    e_x_1=abs(X_side(1)-X1);
    e_y_1=abs(Y_side(1)-Y1);
    E_1=(e_x_1^2+e_y_1^2)^(1/2);
    %first triangular
    S1=0.5*(G*E_1);
    
    %%ertefa2
    e_x_2=abs(X_side(2)-X2);
    e_y_2=abs(Y_side(2)-Y2);
    E_2=(e_x_2^2+e_y_2^2)^(1/2);
    %second triand
    S2=0.5*(G*E_2);
    
    %% total surf
    SS=S1+S2;
end
    
    
    
    