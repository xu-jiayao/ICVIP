clear;
close all;
clc;

path(path, './images/1024&1024');
path(path, './images/result');

imageOriginalPath = './images/1024&1024';
imageFiles = [dir(fullfile(imageOriginalPath,'*png'));
              dir(fullfile(imageOriginalPath,'*tiff'));
              dir(fullfile(imageOriginalPath,'*jpg'));
              dir(fullfile(imageOriginalPath,'*bmp'));
              dir(fullfile(imageOriginalPath,'*mat'))];
numFiles = length(imageFiles);
measurement_matrix_construction = 'binary_hadamard';
image_reconstruction_algorithm  = 'FSAMP';
image_transformation_algorithm  = 'ifwht';

%initialization of parameterss
sub_pixels       = 8;%block size
n   = sub_pixels*sub_pixels;%n
sampling_rate = 0.25;
m  = round(n*sampling_rate);%m

%time record
calculation_time = 0;


for image_loop = 1:2 
    
    %load image
    load_image = imread(imageFiles(image_loop).name);
 
    %judge the image is color or gray
    if size(load_image,3)==3
        original_image = double(rgb2gray(load_image));
    else
        original_image = double(load_image);
    end
    
    
    imwrite(uint8(original_image),fullfile(strcat('./images/dot_product/', ...
        '_original_',  ...
        imageFiles(image_loop).name)));

    
%     get the measurement matrix
    [theta, phi,psi] =  GetMeasurementMatrix(measurement_matrix_construction,...
        image_transformation_algorithm,m,n);
    %seperate the image into cells
    N_1 = zeros(1, size(original_image,1)/sub_pixels) + sub_pixels;
    N_2 = zeros(1, size(original_image,2)/sub_pixels) + sub_pixels;
    C = mat2cell(double(original_image), N_1, N_2);
    
    %coding part
    %initial the numbers
    num_of_rows = size(original_image,1)/sub_pixels;
    num_of_columns = size(original_image,2)/sub_pixels;
    y_deresidual = cell(num_of_rows,num_of_columns);
    
    %coding
    for indexX = 1:num_of_rows
        for indexY = 1:num_of_columns
            one_block_image = reshape(C{indexX,indexY}.',1,[])';
            y_deresidual{indexX,indexY} = BCS_encoder(one_block_image, phi);
        end
    end
    
    
    
    %%%%%%%%%%%%%%%%%%decoding
    %original
    tic
    total_iter_num = 0;
    reconstructed_image= cell(num_of_rows,num_of_columns);


% %     open the mulithreading
    p = gcp('nocreate');
    if isempty(p)
        poolsize = 0;
        parpool('local',4);
    else
        poolsize = p.NumWorkers;
    end
    gd = gpuDevice();



%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%multithreading
%%%%%%%%%%%%%%
    tic
    %initialation of the parameters of the decoding
    if(num_of_rows<num_of_columns)
        Big_num = num_of_columns;
        Small_num = num_of_rows;
    else
        Big_num = num_of_rows;
        Small_num = num_of_columns;
    end
    
    Total_number = Big_num+Small_num+1;%test-->number
    total_iter_num4 = 0;
    total_iter_num3 = 0;
    total_iter_num2 = 0;
    total_iter_num1 = 0;
    total_iter_num = 0;
    
    sum_of_index1 = 2;
    sum_of_index2 = 3;
    sum_of_index3 = 4;
    sum_of_index4 = 5;
    
    reconstructed1 = cell(num_of_rows,num_of_columns);
    reconstructed2 = cell(num_of_rows,num_of_columns);
    reconstructed3 = cell(num_of_rows,num_of_columns);
    reconstructed4 = cell(num_of_rows,num_of_columns);
    
    columns1 = 0;rows1 =0;columns2=0;rows2 = 0;columns3=0;rows3 =0; columns4=0;rows4 = 0;
    

%     timeDifference = 0;
    flagAddThree = 0;
    spmd(4)
        %only for square images
        for blockY = 4:4:num_of_rows
            
            %calculate the index
            TotalNumber = num_of_columns+blockY;
            startRows = blockY-3;
            endRows = blockY;
            sum_of_index = 2+blockY-4;%test-->number
            if flagAddThree ~= 0
                sum_of_index = sum_of_index + 3;
                flagAddThree = 0;
            end
            while(sum_of_index <= TotalNumber)
                
                if(sum_of_index-startRows<num_of_columns+1)
                    rows = startRows;
                    columns = sum_of_index-startRows;
                else
                    columns = num_of_columns;
                    rows = sum_of_index-columns;
                end
                columns1 = columns;
                rows1 =rows;
                columns2=columns1-1;
                rows2 = rows1+1;
                columns3=columns2-1;
                rows3 = rows2+1;
                columns4=columns3-1;
                rows4 = rows3+1;
                
				%test->>why
                if TotalNumber - sum_of_index <3 
                    flagAddThree = 1;
                    if TotalNumber - sum_of_index == 2
                        rows4=blockY+1;
                        columns4 =2+blockY - rows4;
                    elseif TotalNumber - sum_of_index == 1
                        rows3=blockY+1;
                        columns3 =  2+blockY+1 - rows3;
                        columns4=columns3-1;
                        rows4 = rows3+1;
                    elseif TotalNumber - sum_of_index == 0
                        rows2=blockY+1;
                        columns2 =  2+blockY+2 - rows2;
                        columns3=columns2-1;
                        rows3 = rows2+1;
                        columns4=columns3-1;
                        rows4 = rows3+1;
                    end
                end
%                 
%                 biggestTime = 0;
%                 smallestTime = 10;
                
                %seperate into 4 threads
                if labindex == 1
                    
                        [reconstructed1{rows1, columns1},iter_num1] = FSAMP(y_deresidual{rows1,columns1},theta,m,n,n,sub_pixels,0.0001);
                        total_iter_num1 = total_iter_num1+ iter_num1;
                        
                        if columns1~=1
                            [reconstructed1{rows1, columns1},reconstructed1{rows1, columns1-1}] = PostProcess(reconstructed1{rows1, columns1},...
                            reconstructed1{rows1, columns1-1},1);
                        end
%                         time1 = toc;
%                     end
                    
                elseif labindex == 2
%                     tic
                    if((columns2>0)&&(rows2>0)&&(rows2<num_of_rows+1)&&(columns2<num_of_columns+1))
                        
                        [reconstructed2{rows2, columns2},iter_num2] = FSAMP(y_deresidual{rows2,columns2},theta,m,n,n,sub_pixels,0.0001);
                        total_iter_num2 = total_iter_num2+ iter_num2;
                        if columns2~=1
                            [reconstructed2{rows2, columns2},reconstructed2{rows2, columns2-1}] = PostProcess(reconstructed2{rows2, columns2},...
                            reconstructed2{rows2, columns2-1},1);
                        end
                    end
%                     time2 = toc;
                elseif labindex == 3
%                     tic
                    if((columns3>0)&&(rows3>0)&&(rows3<num_of_rows+1)&&(columns3<num_of_columns+1))
                        
                        [reconstructed3{rows3, columns3},iter_num3] = FSAMP(y_deresidual{rows3,columns3},theta,m,n,n,sub_pixels,0.0001);
                        total_iter_num3 = total_iter_num3+ iter_num3;
                        if columns3~=1
                            [reconstructed3{rows3, columns3},reconstructed3{rows3, columns3-1}] = PostProcess(reconstructed3{rows3, columns3},...
                            reconstructed3{rows3, columns3-1},1);
                        end
                    end
%                     time3 = toc;
                else
%                     tic
                    if((columns4>0)&&(rows4>0)&&(rows4<num_of_rows+1)&&(columns4<num_of_columns+1))
                        
                        [reconstructed4{rows4, columns4},iter_num4] = FSAMP(y_deresidual{rows4,columns4},theta,m,n,n,sub_pixels,0.001);
                        total_iter_num4 = total_iter_num4+ iter_num4;
                        if columns4~=1
                            [reconstructed4{rows4, columns4},reconstructed4{rows4, columns4-1}] = PostProcess(reconstructed4{rows4, columns4},...
                            reconstructed4{rows4, columns4-1},1);
                        end
                    end
%                     time4 = toc;
                    %multithreading
                end
                sum_of_index = sum_of_index+1;
            end
            %seperate the image into blocks
        end
        
        if blockY<num_of_rows
            blockY = blockY+1;
            %calculate the index
            TotalNumber = num_of_columns+num_of_rows;
            startRows = blockY;
            endRows = num_of_rows;
            sum_of_index = 1+blockY;
            if flagAddThree ~= 0
                sum_of_index = sum_of_index + 3;
                flagAddThree = 0;
            end
            while(sum_of_index <= TotalNumber)
                
                if(sum_of_index-startRows<num_of_columns+1)
                    rows = startRows;
                    columns = sum_of_index-startRows;
                else
                    columns = num_of_columns;
                    rows = sum_of_index-columns;
                end
                columns1 = columns;
                rows1 =rows;
                if num_of_rows-blockY>0
                    columns2=columns1-1;
                    rows2 = rows1+1;
                end
                if num_of_rows-blockY>1
                    columns3=columns2-1;
                    rows3 = rows2+1;
                end
      
                
                %seperate into 3 threads
                if labindex == 1
%                         tic
                        [reconstructed1{rows1, columns1},iter_num1] = FSAMP(y_deresidual{rows1,columns1},theta,m,n,n,sub_pixels,0.0001);
                        total_iter_num1 = total_iter_num1+ iter_num1;
                        if columns1~=1
                            [reconstructed1{rows1, columns1},reconstructed1{rows1, columns1-1}] = PostProcess(reconstructed1{rows1, columns1},...
                            reconstructed1{rows1, columns1-1},1);
                        end
%                         time1 = toc;
                    
                elseif labindex == 2
%                     tic
                    if((columns2>0)&&(rows2>0)&&(rows2<num_of_rows+1)&&(columns2<num_of_columns+1))
                        
                        [reconstructed2{rows2, columns2},iter_num2] = FSAMP(y_deresidual{rows2,columns2},theta,m,n,n,sub_pixels,0.0001);
                        total_iter_num2 = total_iter_num2+ iter_num2;
                        if columns2~=1
                            [reconstructed2{rows2, columns2},reconstructed2{rows2, columns2-1}] = PostProcess(reconstructed2{rows2, columns2},...
                            reconstructed2{rows2, columns2-1},1);
                        end
                    end
%                     time2 = toc;
                elseif labindex == 3
%                     tic
                    if((columns3>0)&&(rows3>0)&&(rows3<num_of_rows+1)&&(columns3<num_of_columns+1))
                        
                        [reconstructed3{rows3, columns3},iter_num3] = FSAMP(y_deresidual{rows3,columns3},theta,m,n,n,sub_pixels,0.0001);
                        total_iter_num3 = total_iter_num3+ iter_num3;
                        if columns3~=1
                            [reconstructed3{rows3, columns3},reconstructed3{rows3, columns3-1}] = PostProcess(reconstructed3{rows3, columns3},...
                            reconstructed3{rows3, columns3-1},1);
                        end
                    end

%                     
                    %multithreading
                end
                sum_of_index = sum_of_index+1;
                %while
            end
        end
        
        %spmd
    end
    
%     collect all of the data in different thread into one thread
    reconstructed_image = reconstructed1{1};
    reconstructed2 = reconstructed2{2};
    reconstructed3 = reconstructed3{3};
    reconstructed4 = reconstructed4{4};
    for i=1:num_of_rows
        for j=1:num_of_columns
            if  ~isempty(reconstructed2{i,j})
                reconstructed_image{i,j} = reconstructed2{i,j};
            elseif ~isempty(reconstructed3{i,j})
                reconstructed_image{i,j} = reconstructed3{i,j};
            elseif ~isempty(reconstructed4{i,j})
                reconstructed_image{i,j} = reconstructed4{i,j};
            end
        end
    end
    final_image_reconstruct = round(cell2mat(reconstructed_image)); 

    time_record = toc;
    
    

    % parpool close
     delete(gcp('nocreate'));
    
      image_psnr = PSNR(final_image_reconstruct, original_image)
      image_ssim = ssim(final_image_reconstruct, original_image)

   
    imwrite(uint8(final_image_reconstruct),fullfile(strcat('./images/result/', ...
        imageFiles(image_loop).name)));

    
end
