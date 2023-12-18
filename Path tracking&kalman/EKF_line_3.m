clc
clear all
close all
load line_3
         y_path=y_path';
T=1;
A=[1 0 T 0;0 1 0 T;0 0 1 0;0 0 0 1];
B=[0.5*T^2 0;0 0.5*T^2;T 0;0 T];
p_h=H*H';
P(:,:,1)=[p_h(1,1) 0 0 0;0 p_h(2,2) 0 0;0 0 1e6 0;0 0 0 1e6 ];

inv_P(:,:,1)=1e-2*eye(4);
 X_HAT(:,1)=[X_hat(:,1)' 0 0]';

         
for i=2:length(y_path)
     clear r_1 r_2 H_1 H_2 H_k Z_k delta_1 delta_2 delta K 
 varN=0.01^(1/2);
 nn=M_fix(i)+10;
   r=(random('normal',0,varN,nn,2));
   r=r';
        H_k=zeros(2,1);
        Z_k=zeros(2,1);
        delta=[zeros(4,1) zeros(4,1)];
         Q=1e-4*eye(2);
       X_HAT(:,i)=A*X_HAT(:,i-1)+B*Q*[r(1,1) r(2,2)]';

      for j=1:M_fix(i)

    
             N_Best=N_CAND{i};
             M_Best=M_CAND{i};
             x_abs_1(j)=abs(X_HAT(1,i)-(N_Best(1)));
             x_abs_2(j)=abs(X_HAT(1,i)-(N_Best(2)));
             y_abs_1(j)=abs(X_HAT(2,i)-M_Best(1));
             y_abs_2(j)=abs(X_HAT(2,i)-M_Best(2));
             r_1(j)=(x_abs_1(j)^2+y_abs_1(j)^2)^(1/2);
             r_2(j)=((x_abs_2(j))^2+(y_abs_2(j)^2))^(1/2);
             
            H_1(j)=atan2((X_HAT(1,i)-(N_Best(1))),(X_HAT(2,i)-(M_Best(1))));
            H_2(j)=atan2((X_HAT(1,i)-N_Best(2)),(X_HAT(2,i)-M_Best(2)));

            h_1(j)=atan2((X(i)-(N_Best(1))),(y_path(i)-(M_Best(1))));
            h_2(j)=atan2((X(i)-(N_Best(2))),(y_path(i)-M_Best(2)));

             
            Z_1(j)=h_1(j)+r(1,j);
             Z_2(j)=h_2(j)+r(2,j);
         
         H_k=[H_k;H_1(j);H_2(j)];
          Z_k=[Z_k;Z_1(:,j);Z_2(:,j)];
         
      
         delta_1(:,j)=(1/r_1(j))*[cos(H_1(j)) -sin(H_1(j))  0 0 ]';
         delta_2(:,j)=(1/r_2(j))*[cos(H_2(j)) -sin(H_2(j)) 0 0 ]';
         delta=[delta delta_1(:,j) delta_2(:,j)];
      end
         SIZE_DELTA=size(delta);
         
         R=1e-3*eye(SIZE_DELTA(2));
         
         P(:,:,i)=(A*P(:,:,i-1)*A')+B*Q*B';
         inv_P(:,:,i+1)=inv(P(:,:,i))+delta*inv(R)*delta';
         P(:,:,i+1)=inv(inv_P(:,:,i+1));
         K=(P(:,:,i+1))*(delta)*inv(R);
          P(:,:,i)=P(:,:,i+1);  
%           X_HAT(:,i)=A*X_HAT(:,i-1);
           X_HAT(:,i+1)=X_HAT(:,i)+(K*(Z_k-H_k));
              X_HAT(:,i)=X_HAT(:,i+1);
              X_Hat(:,i)=X_HAT(:,i+1);

      end

%  end
 
 
 
% end

figure
plot(X,y_path,'b','linewidth',4)
hold on
plot(X_Hat(1,:),X_Hat(2,:),'y','linewidth',4);
hold on
plot(X_hat(1,:),X_hat(2,:));legend('Real Tatrget','EKF Estimation','first estimation')
plot(n,m,'r*','linewidth',2.5)
hold on
plot(NN(:,:),MM(:,:),'g*');grid

ERROR_R_X=abs(X_Hat(1,:)-X);
ERROR_R_Y=abs(X_Hat(2,:)-y_path');
ERROR_R_XY=(((ERROR_R_X.^2)+(ERROR_R_Y.^2)).^(1/2));
disp('error mean EKF')
mean(ERROR_R_XY)



ERROR_X=abs(X-X_hat(1,:));
ERROR_Y=abs(y_path'-X_hat(2,:));
ERROR_xy=((ERROR_X.^2+ERROR_Y.^2).^(1/2));
disp('error mean first estimation')
mean(ERROR_xy)
figure;
plot(ERROR_xy,'r--','linewidth',3);
hold on
plot(ERROR_R_XY,'b','linewidth',3);legend('first estimation and real','ERROR EKF and real')
% y_path=y_path';
% X_Y_real = {'X', 'y_path'; X;y_path}
% xlswrite('line.xls', X,X_hat(1,:),X_hat_2(1,:), 'X_estimation', 'E1');


% 
% 
%     for i=2:length(y_path)
%         if length(N_CAND{i})==2
%        clear P X_HAT inv_P  Q x_abs_1 x_abs_2 y_abs_1 y_abs_1
%        clear r_1 r_2 H_1 H_2 Z_1 Z_2 H_k Z_k delta _1 delta_2 delta
%         P(:,:,1)=1e6*eye(4);
%         X_HAT(:,1)=[X_hat(:,i)' 0 0]';
% %         X_HAT_d=6*ones(4,M_fix(i));
% %  X_HAT(:,1)=[0.1 0.1 0 0]';   
%          R=0.5;
%          Q=0.01;
%             
%             
%             
%         for j=2:M_fix(i)
%             nn=1000;
% varN=0.01^(1/2);
% r=random('normal',0,varN,nn,3);
%             
%             X_HAT(:,j)=A*X_HAT(:,j-1);
%               N_Best=N_CAND{i};
%               M_Best=M_CAND{i};
%               x_abs_1(j)=abs(X(i)-N_Best(1));
%               x_abs_2(j)=abs(X(i)-N_Best(2));
%               y_abs_1(j)=abs(y_path(i))-M_Best(1);
%               y_abs_2(j)=abs(y_path(i)-M_Best(2));
%               r_1(j)=(x_abs_1(j)^2+y_abs_1(j)^2)^(1/2);
%               r_2(j)=(x_abs_2(j)^2+y_abs_2(j)^2)^(1/2);
%               
%              H_1(j)=atan((y_path(i)-M_Best(1))/(X(i)-N_Best(1)));
%              Z_1(j)=H_1(j)+r(j,1);
%              H_2(j)=atan((y_path(i)-M_Best(2))/(X(i)-N_Best(2)));
%             Z_2(j)=H_2(j)+r(j,2);
%        %first node
% %        [phi_1,R_1] = cart2pol(N_Best(1),M_Best(1));
%         grad_1(j,:)=(1/(r_1(j)))*[-sin(H_1(j)) cos(H_1(j)) 0 0];
%         y_1(j)=((1/(r_1(j)))*grad_1(j,:))*((Z_1(j)-H_1(j))+(grad_1(j,:)'));
%         Y_1(j)=(1/(R))*grad_1(j,:)*grad_1(j,:)';
%         
%         %second node
%         
%         grad_2(j,:)=(1/(r_2(j)))*[-sin(H_2(j)) cos(H_2(j)) 0 0];
%         y_2(j)=((1/(r_2(j)))*grad_2(j,:))*((Z_2(j)-H_2(j))+(grad_2(j,:)'));
%         Y_2(j)=(1/(R))*grad_2(j,:)*grad_2(j,:)';
%         
%         
%         
%         % final estimation and updating cov matrix
% %         P(:,:,i-1)=cov((X_HAT(:,i-1)-X_HAT(:,i-1))*(X_HAT(:,i-1)-X_HAT(:,i-1))');
%         P(:,:,j)=A*P(:,:,j-1)*A';
%          inv_P(:,:,j)=(inv(P(:,:,j-1)))+(Y_1(j)+Y_2(j));
%         X_HAT(:,j+1)=inv(inv_P(:,:,j))*((inv(P(:,:,j-1))*X_HAT(:,j))+y_1(j)+y_2(j));
%         end
%         end
%     X_Hat(:,i)=X_HAT(:,j+1);
%     end
% 
%     figure
% plot(X,y_path,'b','linewidth',4)
% hold on
% plot(X_Hat(1,:),X_Hat(2,:),'y','linewidth',4);legend('Real Tatrget','EKF Estimation')
% plot(n,m,'r*','linewidth',2.5)
% hold on
% plot(NN(:,:),MM(:,:),'g*');grid
% 
% ERROR_R_X=abs(X_Hat(1,:)-X);
% ERROR_R_Y=abs(X_Hat(2,:)-y_path');
% ERROR_R_XY=(ERROR_R_X.^2+ERROR_R_Y.^2).^(1/2);
% 
% 
% ERROR_x=abs(X_Hat(1,:)-X_hat(1,:));
% ERROR_y=abs(X_Hat(2,:)-X_hat(2,:));
% ERROR_XY=(ERROR_x.^2+ERROR_y.^2).^(1/2);
% 
% 
% 
% ERROR_X=abs(X-X_hat(1,:));
% ERROR_Y=abs(y_path'-X_hat(2,:));
% ERROR_xy=(ERROR_X.^2+ERROR_Y.^2).^(1/2);
% figure;
% plot(ERROR_xy,'r--','linewidth',3);
% hold on
% plot(ERROR_R_XY,'b','linewidth',3);legend('first estimation and real','ERROR EKF and real')
% 
