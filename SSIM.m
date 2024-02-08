function [ reconstructed_image_before_trans,iter_num  ] = FSAMP_Before_Trans( y_deresidual,theta,S, sub_pixels,sigma)
%CS_SAMP Summary of this function goes here
%Version: 1.0 written by jbb0523 @2015-05-08
%   Detailed explanation goes here
%   y_deresidual = Phi * x
%   x = Psi * theta
%	y_deresidual = Phi*Psi * theta
%   let     theta = Phi*Psi, then    y=theta*reconstructed_image
%   S is step
%   now we know     y_deresidual and theta?we are going to get     reconstructed_image
%   Reference:Thong T.Do?Lu Gan?Nam Nguyen?Trac D.Tran?Sparsity adaptive
%   matching pursuit algorithm for practical compressed sensing[C]?Asilomar
%   Conference on Signals?Systems?and Computers?Pacific Grove?California?
%   2008?10?581-587.
%   Available at:
%   http://dsp.rice.edu/sites/dsp.rice.edu/files/cs/asilomar08_final.pdf

enter = 1;
[y_rows,y_columns] = size(y_deresidual);
if y_rows<y_columns
    y_deresidual = y_deresidual';%y_deresidual should be a column vector
end
[M,N] = size(theta);%the measurement matrix theta is m*n
reconstructed_image_before_trans = zeros(N,1);%output
r_n = y_deresidual;%  r0=y
iter_num=0;
St = atan(S)*M*pi;%atan(S)*M*pi
Ds = M/4;
if(Ds < 1)
    Ds = 1;
end
if(S>St/4)
    size_of_step = linspace(S,St/4,Ds);
else
    size_of_step = linspace(St/4,S,Ds);%original
end
size_of_step1 =round(size_of_step);
L = size_of_step1(1,1);%(Size of the finalist in the first stage)
if norm(y_deresidual) == 0
    Pos_theta = floor(linspace(1,y_rows+1,y_rows));
    Pos_theta = Pos_theta';
    theta_ls = zeros(y_rows,1);
    enter = 0;
    iter_num=1;
else
    Pos_theta = [];% final list
end
r_lastIter = 0;

while(enter)
    % tic
    if  L>N %length(Ck)<=M (1,1)
        L = N;
    end
    %(1)Preliminary Test
    product = theta'*r_n;  %cannot understand
    [~,pos]=sort(abs(product),'descend');
    Sk = pos(1:L);% short list
    %(2)Make Candidate List
    Ck = union(Pos_theta,Sk); %candidate list  ,'stable'
    %(3)Final Test
    At = theta(:,Ck);
    [ck,~]=size(Ck);
    if(ck>M)
         theta_ls = At'*(At*At')^(-1)*y_deresidual;
    else
        theta_ls = (At'*At)^(-1)*At'*y_deresidual;
    end
   
    [~,pos] = sort(abs(theta_ls),'descend');
    F = Ck(pos(1:L));
    
    %(4)Compute Residue
    if(L>M||L==M)
         theta_ls =  theta(:,F)'*(theta(:,F)*theta(:,F)')^(-1)*y_deresidual;
    else
        theta_ls =  (theta(:,F)'*theta(:,F))^(-1)*theta(:,F)'*y_deresidual;
    end
     r_new = y_deresidual - theta(:,F)*theta_ls;
    iter_num = iter_num+  1;
    
    %(3)Final Test
    if(norm( norm(r_new) - norm(r_lastIter) )<sigma || norm(r_new) < sigma || L == N ) %norm( norm(r_new)-norm(r_n))<sigma
        Pos_theta = F;
        break;
    end
    
    r_lastIter = r_new;
    
    if norm(r_new)>=norm(r_n)
        times = norm(r_n)/norm(r_new);
        L =round( L- times*L);
    else
        Pos_theta = F;
        r_n = r_new;%where to put
        if iter_num>size(size_of_step)
            L = L+size_of_step1;%Update the size of finalist
        else
            L = L+size_of_step1(1,iter_num);%Update the size of finalist
        end 
    end
end
reconstructed_image_before_trans(Pos_theta)=theta_ls;

end
