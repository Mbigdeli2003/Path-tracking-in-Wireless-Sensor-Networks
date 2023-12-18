clc
clear all
close all
tic
%% defing surface
x=300;%length
y=100; %width
% defining reseloution
res=100;
X=linspace(0,x,res);
Y=linspace(0,y,res);
load X_daraje_3
load Y_daraje_3
    
%% fitting and evaluating
%% data selection randomly for path
% polynomial estimation and calculating the path 
%polynomial
% number of point selection
% N=20;
% randomly selection point for polynomial function estimation
% for i=1:N
% xdata(i)=round(rand*300);
% ydata(i)=round(rand*100);
% end
% % curve fiiting and setting the path 
% curve = fit( xdata', ydata', 'poly3' );
% y_path=abs(curve(X));
% X=linspace(0,300,100);
% y_path=linspace(0,100,100);

% for i=1:length(y_path)
%     if y_path(i)>y
%         y_path(i)=y;
%     end
% end
    
    
%% slecet random nodes
% node number

% % random node selection
% % n is x and m is y of observer nodes
n=[zeros(1,11) 50*ones(1,11) 100*ones(1,11) 150*ones(1,11) 200*ones(1,11) 250*ones(1,11) 300*ones(1,11)];
m=[0:10:100 0:10:100 0:10:100 0:10:100 0:10:100 0:10:100 0:10:100];
N=max(size(n));

%% calculating distance from observor nodes to path points
% SRE defining and selecting sensors
% NN is x of selected sensors and MM is y of selected sens
 SR=0.1;
% M_SR=zeros(length(y_path),N);
% N_SR=zeros(length(y_path),N);
for i=1:length(y_path)
for j=1:N
   dis_X(i,j)=abs(n(j)-X(i));
   dis_Y(i,j)=abs(m(j)-y_path(i));
  dis_X_Y(i,j)=(dis_X(i,j)^2+dis_Y(i,j)^2)^(1/2);
  
end
end


for i=1:length(y_path)
    k=0;p=0;
    for j=1:N
        if dis_X_Y(i,j)<SR*max(max(dis_X_Y))   
 N_SR(i,j)=n(j);
 M_SR(i,j)=m(j);
        else
            N_SR(i,j)=nan;
            M_SR(i,j)=nan;
        end
        
         if (N_SR(i,j)>0 || N_SR(i,j)==0)
%             && (N_SR(i,j)~=N_SR(i-1,j-1))
            k=k+1;
            if k<=3
         NN(i,k)=n(j);
            end
          end

              if (M_SR(i,j)>0 || M_SR(i,j)==0) 
%                  && (M_SR(i,j)~=M_SR(i-1,j-1))
                 p=p+1;
                 if p<=3
              MM(i,p)=m(j);
                 end
              end 
        
    end
end
% end
%% noise geeneration mean=0;variance=0.001;
nn=1000;
varN=0.01^(1/2);
r=random('normal',0,varN,nn,3);

%% calculating angles with 3 selected sensors for first estimation with
%% noise
for i=1:length(y_path)
for j=1:min(size(MM))  
    a(i,j)=NN(i,j)-X(i);
    b(i,j)=(MM(i,j)-y_path(i));
    if (NN(i,j)-X(i))==0
         a(i,j)=(NN(i,j)-(X(i)))+1;
%          
    end
    if (MM(i,j)-y_path(i))==0
       b(i,j)=(MM(i,j)-(y_path(i)))+1;
    end
BETA(i,j)=a(i,j)/b(i,j);
Beta(i,j)=atan2(a(i,j),b(i,j));
z(i,j)=Beta(i,j)+r(i,j);
end
end

%% first estimate with three sensors
for i=1:length(y_path) 
j=1;  
 Z=[(MM(i,j)*sin(z(i,j)))-(NN(i,j)*cos(z(i,j)));
     MM(i,j+1)*sin(z(i,j+1))-NN(i,j+1)*cos(z(i,j+1));
     MM(i,j+2)*sin(z(i,j+2))-NN(i,j+2)*cos(z(i,j+2))];
 H=[-cos(z(i,j)) sin(z(i,j));
    -cos(z(i,j+1)) sin(z(i,j+1));
    -cos(z(i,j+2)) sin(z(i,j+2))];
X_hat(:,i)=pinv(H)*Z;
H_imp{i}=H
end
   for j=1:length(y_path) 
if X_hat(1,j)>300
    X_hat(1,j)=300;
end
if X_hat(1,j)<0
    X_hat(1,j)=0;
end
if X_hat(2,j)>100
    X_hat(2,j)=100;
end
if X_hat(2,j)<0
    x_hat(2,j)=0;
   end
   end
    
    %first estimate with 2 sensors
for i=1:length(y_path) 
j=1;  
 Z=[(MM(i,j)*sin(z(i,j)))-(NN(i,j)*cos(z(i,j)));
     MM(i,j+1)*sin(z(i,j+1))-NN(i,j+1)*cos(z(i,j+1))];
 H=[-cos(z(i,j)) sin(z(i,j));
    -cos(z(i,j+1)) sin(z(i,j+1))];
X_hat_2(:,i)=pinv(H)*Z;
end

   for j=1:length(y_path) 
if X_hat(1,j)>300
    X_hat(1,j)=300;
end
if X_hat(1,j)<0
    X_hat(1,j)=0;
end
if X_hat(2,j)>100
    X_hat(2,j)=100;
end
if X_hat(2,j)<0
    X_hat(2,j)=0;
end
   end





   for j=1:length(y_path) 
if X_hat_2(1,j)>300
    X_hat_2(1,j)=300;
end
if X_hat_2(1,j)<0
    X_hat_2(1,j)=0;
end
if X_hat_2(2,j)>100
    X_hat_2(2,j)=100;
end
if X_hat_2(2,j)<0
    X_hat_2(2,j)=0;
end
   end

%% calculating sampling number of nodes
R=1;
for i=1:length(y_path)
for j=1:min(size(N_SR))
   Dis_X(i,j)=abs(N_SR(i,j)-X(i));
   Dis_Y(i,j)=abs(M_SR(i,j)-y_path(i));
  Dis_X_Y(i,j)=((Dis_X(i,j)^2)+(Dis_Y(i,j)^2))^(1/2);
  
end
end

for i=1:length(y_path)
for j=1:min(size(N_SR))
    if Dis_X_Y(i,j)==max(Dis_X_Y(i,:))
        M(i)=0.01/(((1/3)*(asin(R/Dis_X_Y(i,j))))^2);
        M_fix(i)=ceil(M(i));
    
end
end
end

for i=1:length(y_path)
for j=1:min(size(N_SR))
    if Dis_X_Y(i,j)==max(Dis_X_Y(i,:))
        Sigma(i,j)=0.01;
    end
    if Dis_X_Y(i,j)>0 || Dis_X_Y(i,j)==0
    Sigma(i,j)=((1/3)*(asin(R/Dis_X_Y(i,j))))^2*M(i);
    else
      Sigma(i,j)=NaN;
   
end
end
end
for i=1:length(y_path)
    k=0;
    for j=1:min(size(N_SR))
%      X_i=X(i);
%      y_path_i=y_path(i);
          
     if Sigma(i,j)>0 || Sigma(i,j)==0
         k=k+1;
         SIG(i,k)=Sigma(i,j);
         ncand(i,k)=N_SR(i,j);
         mcand(i,k)=M_SR(i,j);
     end
    end
end








% 
% % 2 out of n
for i=1:length(y_path)
    for j=1:min(size(SIG))
        k=j+1;
        for p=k:min(size(SIG))
            
       if SIG(i,j)>0 && SIG(i,p)>0
           SIGMA_i=[SIG(i,j) SIG(i,p)];
           n_cand=[ncand(i,j) ncand(i,p)];
           m_cand=[mcand(i,j) mcand(i,p)];
           X_i=X(i);
           y_path_i=y_path(i);
           try
           S_t=Surface(X_i,y_path_i,n_cand,m_cand,SIGMA_i);
           catch exception
            S_t=0;
           end
           s_2(j,p,i)=S_t;
           
       end
    end
end
end



for i=1:length(y_path)
    for j=1:min(size(SIG))
        k=j+1;
        for p=k:min(size(SIG))
            
       if SIG(i,j)>0 && SIG(i,p)>0
           SIGMA_i=[SIG(i,j) SIG(i,p)];
           n_cand=[ncand(i,j) ncand(i,p)];
           m_cand=[mcand(i,j) mcand(i,p)];
           X_i=X(i);
           y_path_i=y_path(i);
           try
           S_t=Surface(X_i,y_path_i,n_cand,m_cand,SIGMA_i);
           catch exception
            S_t=0;
           end
           S_2(j,p,i)=S_t;
           if S_2(j,p,i)==min(min(nonzeros(s_2(:,:,i)))) 
               N_CAND{i}=[ncand(i,j) ncand(i,p)];
                  M_CAND{i}=[mcand(i,j) mcand(i,p)];
                  SIGMA{i}=[SIG(i,j) SIG(i,p)];
                  end
               
           
       end
    end
end
end 



% 3 out of n

% for i=1:length(y_path)
%     for j=1:min(size(SIG))
%         k=j+1;
%         for p=k:min(size(SIG))
%             l=p+1;
%             for o=l:min(size(SIG))
%               if SIG(i,j)>0 && SIG(i,p)>0 && SIG(i,o)>0
%                   SIGMA_i_3=[SIG(i,j) SIG(i,p) SIG(i,o)];
%                   n_cand_3=[ncand(i,j) ncand(i,p) ncand(i,o)];
%                   m_cand_3=[mcand(i,j) mcand(i,p) mcand(i,o)];
%                   X_i_3=X(i);
%                   y_path_i_3=y_path(i);
%                   try
%                   [S NN MM N M]=Surface_3(X_i_3,y_path_i_3,n_cand_3,m_cand_3,SIGMA_i_3);
%                   catch exception
%                       S=0;
%                   end
%                   S_3(j,p,o,i)=S;
%                   
% %                    close all
% %                   plot(XX,YY,'g.','linewidth',4);
% %                   hold on
% %                    plot(N,M,'r*','linewidth',5);
% %                    hold on
% %                    plot(NN,MM,'b*','linewidth',6);
% %                    pause
%               end
%             end
%         end
%     end
% end



%% finding minimum surface and the best sensors to estimate location

% for i=1:length(y_path)
%     for j=1:min(size(SIG))
%         k=j+1;
%         for p=k:min(size(SIG))
%             
%               if SIG(i,j)>0 && SIG(i,p)>0
%            SIGMA_i=[SIG(i,j) SIG(i,p)];
%            n_cand=[ncand(i,j) mcand(i,p)];
%            m_cand=[mcand(i,j) mcand(i,p)];
%            X_i=X(i);
%            y_path_i=y_path(i);
%            try
%            S_t=Surface(X_i,y_path_i,n_cand,m_cand,SIGMA_i);
%            catch exception
%                S_t=0;
%            end
%            s_2(j,p,i)=S_t;
%            
%        end
%             
%             
%             l=k+1;
%             for o=l:min(size(SIG))
%               if SIG(i,j)>0 && SIG(i,p)>0 && SIG(i,o)>0
%                   SIGMA_i_3=[SIG(i,j) SIG(i,p) SIG(i,o)];
%                   n_cand_3=[ncand(i,j) ncand(i,p) ncand(i,o)];
%                   m_cand_3=[mcand(i,j) mcand(i,p) mcand(i,o)];
%                   X_i_3=X(i);
%                   y_path_i_3=y_path(i);
%                   try
%                   S=Surface_3(X_i_3,y_path_i_3,n_cand_3,m_cand_3,SIGMA_i_3);
%                   catch exception
%                       S=0;
%                   S_3(j,p,o,i)=S;
%                   end
%                   if S_3(j,p,o,i)==min(min(min(S_3(:,:,:,i)))) && S_3(j,p,o,i)==min(min(s_2(:,:,i))) 
%                       N_CAND{i}=[ncand(i,j) ncand(i,p) ncand(i,o)];
%                      M_CAND{i}=[mcand(i,j) mcand(i,p) mcand(i,o)];
%                      SIGMA{i}=[SIG(i,j) SIG(i,p) SIG(i,o)];
%                   elseif S_2(j,p,i)==min(min(min(S_3(:,:,:,i)))) && S_2(j,p,i)==min(min(s_2(:,:,i))) 
%                   N_CAND{i}=[ncand(i,j) mcand(i,p)];
%                   M_CAND{i}=[mcand(i,j) mcand(i,p)];
%                   SIGMA{i}=[SIG(i,j) SIG(i,p)];
%                   end
%                       
%               end
%             end
%         end
%     end
% end
           
           
%% extended kalman filter_final estimation
% defining matrix and cosntants
%% with constant spped movement 
% T is time between snapshots
% v is speed of target
%c is sound speed
% T=1;
% v=10;
% c=347;
% F=[1 0 T 0;0 1 0 T;0 0 1 0;0 0 0 1];
% nn=1000;
% varN=0.01^(1/2);
% r=random('normal',0,varN,nn,3);
% V_X=zeros(1,length(X));
% V_Y=zeros(1,length(X));
% 
% if length(N_CAND{i})==2
%     for i=2:length(y_path_i)
%         saw(i)=atan((y_path(i)-y_path(i-1))/(X(i)-X(i-1)));
%         for j=1:M_fix(i)
%             nn=1000;
% varN=0.01^(1/2);
% r=random('normal',0,varN,nn,3);
%             
%             N_Best=N_CAND{i};
%             M_Best=M_CAND{i};
%            H_1(i,j)=atan((y_path(i)-N_Best(1))/(X(i)-M_Best(1)));
%            Z_1(i,j)=atan((y_path(i)-N_Best(1))/(X(i)-M_Best(1)))+r(j);
%            H_2(i,j)=atan((y_path(i)-N_Best(2))/(X(i)-M_Best(2)));
%            Z_2(i,j)=atan((y_path(i)-N_Best(2))/(X(i)-M_Best(2)))+r(j);
%         H=[H_1(i,j);H_2(i,j)];
%         Z=[Z_1(i,j);Z_2(i,j)];
%         V_X(i)=v*saw(i);
%         V_Y(i)=v*saw(i);
%         X_HAT(:,i-1)=[X(i-1);y_path(i-1);V_X(i-1);V_Y(i-1)];
%        X_HAT(:,i)=F*X_hat(:,i-1);
%        
%        %first node
%        [phi_1,R_1] = cart2pol(N_Best(1),M_Best(1));
%         grad_1(:,i)=(1/(R_1))*[-sin(phi_1);cos(phi_1);0;0];
%         y_1(i)=(1/(r(j)))*grad_1(:,i)'*(Z_1(i,j)-H_1(i,j))+(grad_1(:,i)'*X_HAT(:,i));
%         Y_1(i)=1/(r(j))*grad_1(:,i)*grad_1(:,i)';
%         
%         %second node
%         [phi_2,R_2] = cart2pol(N_Best(2),M_Best(2));
%         grad_2(:,i)=(1/(R_2))*[-sin(phi_2);cos(phi_2);0;0];
%         y_2(i)=(1/(r(j)))*grad_2(:,i)'*(Z_2(i,j)-H_2(i,j))+(grad_2(:,i)'*X_HAT(:,i));
%         Y_2(i)=1/(r(j))*grad_2(:,i)*grad_2(:,i)';
%         
%         
%         
%         % final estimation and updating cov matrix
%         P(:,:,i-1)=cov((X_HAT(:,i-1)-X_HAT(:,i-1))*(X_HAT(:,i-1)-X_HAT(:,i-1))');
%         P(:,:,i)=inv(P(:,:,i-1))+(Y_1(i)+Y_2(i));
%         X_HAT(:,i)=P(:,:,i-1)*((inv(P(:,:,i-1))*X_HAT(:,i))+y_1(i)+y_2(i));
%         end
%     end
% end
%         
%         
% 
% 
% 
% 
% 
% if length(N_CAND{i})==3
%     for i=2:length(y_path_i)
%         saw(i)=atan((y_path(i)-y_path(i-1))/(X(i)-X(i-1)));
%         for j=1:M_fix(i)
%             nn=1000;
% varN=0.01^(1/2);
% r=random('normal',0,varN,nn,3);
%             
%             N_Best=N_CAND{i};
%             M_Best=M_CAND{i};
%            H_1(i,j)=atan((y_path(i)-N_Best(1))/(X(i)-M_Best(1)));
%            Z_1(i,j)=atan((y_path(i)-N_Best(1))/(X(i)-M_Best(1)))+r(j);
%            H_2(i,j)=atan((y_path(i)-N_Best(2))/(X(i)-M_Best(2)));
%            Z_2(i,j)=atan((y_path(i)-N_Best(2))/(X(i)-M_Best(2)))+r(j);
%            H_3(i,j)=atan((y_path(i)-N_Best(3))/(X(i)-M_Best(3)));
%            Z_3(i,j)=atan((y_path(i)-N_Best(2))/(X(i)-M_Best(2)))+r(j);
%            H=[H_1(i,j);H_2(i,j);H_3(i,j)];
%         Z=[Z_1(i,j);Z_2(i,j);Z_3(i,j)];
%         V_X(i)=v*saw(i);
%         V_Y(i)=v*saw(i);
%         X_HAT(:,i-1)=[X(i-1);y_path(i-1);V_X(i-1);V_Y(i-1)];
%         X_HAT(:,i)=F*X_hat(:,i-1);
%         
%         %first node
%         [phi_1,R_1] = cart2pol(N_Best(1),M_Best(1));
%         grad_1(:,i)=(1/(R_1))*[-sin(phi_1);cos(phi_1);0;0];
%         y_1(i)=(1/(r(j)))*grad_1(:,i)'*(Z_1(i,j)-H_1(i,j))+(grad_1(:,i)'*X_HAT(:,i));
%         Y_1(i)=1/(r(j))*grad_1(:,i)*grad_1(:,i)';
%         
%         %second node
%         [phi_2,R_2] = cart2pol(N_Best(2),M_Best(2));
%         grad_2(:,i)=(1/(R_2))*[-sin(phi_2);cos(phi_2);0;0];
%         y_2(i)=(1/(r(j)))*grad_2(:,i)'*(Z_2(i,j)-H_2(i,j))+(grad_2(:,i)'*X_HAT(:,i));
%         Y_2(i)=1/(r(j))*grad_2(:,i)*grad_2(:,i)';
%         
%        %third node
%         [phi_3,R_3] = cart2pol(N_Best(3),M_Best(3));
%         grad_3(:,i)=(1/(R_3))*[-sin(phi_3);cos(phi_3);0;0];
%         y_3(i)=(1/(r(j)))*grad_3(:,i)'*(Z_3(i,j)-H_3(i,j))+(grad_3(:,i)'*X_HAT(:,i));
%         Y_3(i)=1/(r(j))*grad_3(:,i)*grad_3(:,i)';
%         
%         
%         % final estimation and update cov matrix
%         P(:,:,i-1)=cov((X_HAT(:,i-1)-X_HAT(:,i-1))*(X_HAT(:,i-1)-X_HAT(:,i-1))');
%         P(:,:,i)=inv(P(:,:,i-1))+(Y_1(i)+Y_2(i)+Y_3(i));
%         X_HAT(:,i)=P(:,:,i-1)*((inv(P(:,:,i-1))*X_hat(:,i))+y_1(i)+y_2(i)+y_3(i));
%         end
%     end
% end




figure
plot(X,y_path,'b','linewidth',4);
hold on
plot(X_hat(1,:),X_hat(2,:),'y','linewidth',4);legend('Real Target','First Estimation')
plot(n,m,'r*','linewidth',2.5)
hold on
plot(NN(:,:),MM(:,:),'g*');grid




% figure
% plot(X,y_path,'b','linewidth',4)
% hold on
% plot(X_HAT(1,:),X_hat(2,:),'y','linewidth',4);legend('Real Tatrget','EKF Estimation')
% % plot(n,m,'r*','linewidth',2.5)
% % hold on
% % plot(NN(:,:),MM(:,:),'g*');grid
% 
% toc




% figure;
%  plot(M_fix,'r*','linewidth',3)