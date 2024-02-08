function [ reconstructed_image,iter_num  ] = FSAMP( y_deresidual,theta,M,N,S, sub_pixels,sigma)
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

% A(:,pos)=zeros(M,1);                          %  ?????????????????????????

enter = 1;
[y_rows,y_columns] = size(y_deresidual);
if y_rows<y_columns
    y_deresidual = y_deresidual';%y_deresidual should be a column vector
end
% [M,N] = size(theta);%the measurement matrix theta is m*n
% S = M/2; %  attention---------------------
reconstructed_image = zeros(N,1);%output
r_n = y_deresidual;%  r0=y
iter_num=0;
St = atan(S)*M/4/pi;%atan(S)*M/4pi
Ds = M/4;
if(Ds < 1)
    Ds = 1;
end
% if(S>St/4)
%     size_of_step = linspace(S,St,Ds);
% else
    size_of_step = linspace(St,S,Ds);%original
% end
size_of_step1 =round(size_of_step);
L = size_of_step1(1,1);%(Size of the finalist in the first stage)
% step_canIndex = 2;
theta_ls_result = [];
if norm(y_deresidual) == 0
    Pos_theta = floor(linspace(1,y_rows+1,y_rows));
    Pos_theta = Pos_theta';
    theta_ls_result = zeros(y_rows,1);
    enter = 0;
    iter_num=1;
else
    Pos_theta = [];% final list
end



At = theta(:,1:M);
theta_ls = (At'*At)^(-1)*At'*y_deresidual; 
[~,pos] = sort(abs(theta_ls),'descend');
while(enter)
   if  L>M 
       L=M;
   end
        % L = M;
    F = pos(1:L);
    %(4)Compute Residue
    theta_ls =  (theta(:,F)'*theta(:,F))^(-1)*theta(:,F)'*y_deresidual;

    r_new = y_deresidual - theta(:,F)*theta_ls;
    iter_num = iter_num+  1;
    
    %(3)Final Test
    if( norm(r_new) < sigma || L == M ) 
        Pos_theta = F; 
        theta_ls_result = theta_ls;
%         iter_num
        break;
    end
    
    if norm(r_new)>norm(r_n)
        break; %%ifwht & hadamard not enter
    else
        Pos_theta = F;
        theta_ls_result = theta_ls;
        r_n = r_new;
        if iter_num>size(size_of_step)
            L = L+size_of_step1(1,1);%Update the size of finalist
        else
            L = L+size_of_step1(1,iter_num);%Update the size of finalist
        end 
    end 
end

% Pos_theta
% theta_ls_result
reconstructed_image(Pos_theta)=theta_ls_result;
reconstructed_image = roundn(reconstructed_image,-10);

%inverse transformation
final_image_transformation =ifwht(reconstructed_image);
reconstructed_image = reshape(final_image_transformation.',sub_pixels,[])';


% if floor(final_image_transformation1) ~= floor(reconstructed_image)
%     disp("not same");
% end


end
