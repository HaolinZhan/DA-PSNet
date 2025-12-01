clc 
clear all
close all
%% 
sw=4400; % spectral width in the direct dimension spectral width in the direct dimension
sw1=50;  %spectral width in the indirect dimension

chunk_point=88; %chunk_point=SW/SW1, namely the number of sampling points in each chunk
ni=93;%The number of pure shift data chunks

point=8184;
fn=8184;
t21 = 1/sw * (0 : chunk_point-1);% 
t2=t21(1:chunk_point);
t_1 = 1/sw1* (0 :ni-1);% t_
t1=t_1';

peak=0;
num=1; %The number of training datasets

%%
ori_signal_2D=zeros(ni,chunk_point);
idl_signal_2D=zeros(ni,chunk_point);
% 
nus_chunk_blank1=zeros(point,1);
origins=zeros(fn,num);
ideals=zeros(fn,num);
noises=zeros(fn,num);
%%
for i=1:num    
    peak=0;
    while(peak<7)       
        for n=1:26
            f=1000+(3000-1000)*rand(1);
            A=0.05+(1-0.05)*rand(1);
            T2=0.08+(0.6-0.08)*rand(1);
            J=2+(12-2)*rand(1);
            J_num=round(rand(1)*9);
            p=(0+(1-0)*rand(1))<0.5;
                       
            %% 
            orgin_1D=(A*exp(1i*2*pi*f*(t2)).*exp(-t2/T2).*((cos(pi*J*(t2-1/(2*sw1))))).^J_num);
            ori_signal_2D(:,:,n)=(p*exp(1i * 2*pi*f* t1) .* exp(-t1/T2)) *orgin_1D;
            %% 
            ideal_1D=(A*exp(1i*2*pi*f*(t2)).*exp(-t2/T2));
            idl_signal_2D(:,:,n)=(p*exp(1i * 2*pi*f* t1) .* exp(-t1/T2) )*ideal_1D ;           
            if(p==1)
                peak=peak+1;
            end
        end       
    end     
    
    finalori_signal_2D =sum(ori_signal_2D,3);
    finalidl_signal_2D=sum(idl_signal_2D,3);
    
    ori_signal_part=finalori_signal_2D.';
    origin=ori_signal_part(:);
    
    
    idl_signal_part=finalidl_signal_2D.';
    ideal=idl_signal_part(:);
    

noise_sd=0.15;
% noise_sd=(2e-7)+((8e-7)-(2e-7))*rand(1); %2021-AC-DN-Unet
% noise_sd=(0)+((0.5/30)-(0))*rand(1);%2021-JPCL-PS-Resnet,A=1-30
% noise_sd=(0)+((8.5e-4/20)-(0))*rand(1);%2023-JPCL,A=1-20
% noise_sd=1e-4;%2019-Angewdte
% noise_sd=(0.001)+((0.15)-(0.001))*rand(1);%2021-JACS
% noise_sd=(0.001)+((0.03)-(0.001))*rand(1);%2021-JBNMR
noise=noise_sd*(randn(point,1)+1i*randn(point,1));%生成均值为0，方差为noise_sd平方，即标准差为noise_sd的高斯白噪声；

noise_ori=origin+noise;


 %%
    for c=1:ni  % fully sampled 
%     for c=[1,2,3,4,5,7,8,10,12,15,18,22,28,37,93] % sparsely sampled mask
        s=chunk_point*(c-1)+1;
        e=chunk_point*c;

        seg_orgin=noise_ori(s:e);
        nus_chunk_blank1(s:e)=seg_orgin;
        
        nus_origin=nus_chunk_blank1;
    end

%%
origin_sp=(fft(origin,fn));

origins(:,i)=origin_sp;
ideal_sp=(fft(ideal,fn));
ideals(:,i)=ideal_sp;


noise_origin_sp=(fft(noise_ori,fn));

noise_nus_ori_sp=(fft(nus_origin,fn));
noises(:,i)=noise_nus_ori_sp;


end    

  real_ideal_sp=real(ideals);
  max_idl_value=max(abs(real_ideal_sp));

  norm_ideal_sp_real=real_ideal_sp./max_idl_value;

  real_origin_sp=real(noises);
  max_ori_value=max(abs(real_origin_sp));
  norm_origin_sp_real=real_origin_sp./max_ori_value;

%%
figure(10);plot(norm_ideal_sp_real); 
title('Normalized ideal spectra');
% %% ****************************************      
figure(11);plot(((norm_origin_sp_real)));
title('Normalized original spectra');