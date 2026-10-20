T_fav_pred0=zeros(1096,6232);
z_values = linspace(0.00001, 0.99999, 10000);
Tv = spline(flamelet_Z,flamelet_T,z_values);


t0 = now;

for i=1:548
    disp(['row ',num2str(i),' of 548'])
    for k=1:6232
x_cor = i;
z_cor = k;
%try

if i==1
    tic
end

z_ex = Z_fav(x_cor,z_cor);
z_var_ex = Z_var(x_cor,z_cor);

if z_var_ex<10e-8  
    z_var_ex=0; 
    if (1-z_ex) > (z_ex)
        z_ex = 0;
    else
        z_ex = 1; 
    end
end

%disp(['Z = ', num2str(z_ex)]);
%disp(['Z_var = ', num2str(z_var_ex)]);
%disp(['X_coor = ', num2str(x_cor)]);
%disp(['Z_coor = ', num2str(z_cor)]);


a = (z_ex) *( ((z_ex) * (1-z_ex) / (z_var_ex)) -1);
b = (a / z_ex) - a;

if isnan(a) || isnan(b)
    % Degenerate case: Assign 1 to the PDF at z = z_ex
    % p_z = @(z) (z == z_ex); % This is exact if z_ex matches a z value perfectly
    % z_values = linspace(0, 1, 300);
    % % Calculate the PDF values for each z
    % pdf_values = p_z(z_values);
    % max_pdf = max(pdf_values); %To normalize probability
    % pdf_values = pdf_values/max_pdf; %Probability normalized btw 0 and 1

    %Tv = TofZ(z_values);
    % T_fav_pred0_ik = sum(Tv.*pdf_values)/sum(pdf_values); 
    %T_fav_pred0_ik = TofZ(z_ex);
    if z_ex==0
        T_fav_pred0_ik=flamelet_Z(1);
    end

    if z_ex==1
        T_fav_pred0_ik=flamelet_Z(end);
    end

    T_fav_pred0(i,k)=T_fav_pred0_ik;
else
    %B = beta(a,b);
    %p_z = @(z) (1 / B) * (z.^(a - 1)) .* ((1 - z).^(b - 1));
    
    % Calculate the PDF values for each z
    %pdf_values = p_z(z_values);
    pdf_values = betapdf(z_values,a,b);

    
    T_fav_pred0_ik = sum(Tv.*pdf_values)/sum(pdf_values); 

    T_fav_pred0(i,k)=T_fav_pred0_ik;
    T_fav_pred0(1096+1-i,k)=T_fav_pred0_ik;
end

% % Calculate the PDF values for each z
% pdf_values = p_z(z_values);
% max_pdf = max(pdf_values); %To normalize probability
% 
% Tv = TofZ(z_values);
% T_fav_pred0_ik = sum(Tv.*pdf_values)/sum(pdf_values); 
% 
% T_fav_pred0(i,k)=T_fav_pred0_ik;


    end
    if i==10
        toc
    end

    t_used = now-t0;
    vel = i/t_used;
    time_remaining = (548-i)/vel;
    disp(['Velocity: ',num2str(vel*24*3600),' seconds per row'])
    disp(['Expected end date ',datestr(now+time_remaining)])
    disp('-----------------------------------------------------------')
   
end

%% Visualization

% tic
% plot2DDataC(T_fav_pred0,1/728,'nws');
% toc

% T_fav_half = zeros(1096, 6232);  
% T_fav_half(1:548, :) = T_fav(1:548, :);  
% MASK=double([T_fav_half>0]);
% printResultPNG(T_fav_half,MASK,0,3500,'nws','T_fav.png');
tic
T_fav_half = zeros(1096, 6232);  
T_fav_half(1:548, :) = T_fav(1:548, :);  
MASK=double([T_fav_half>0]);
% printResultPNG(T_fav_pred0-T_fav_half,MASK,-3000,1500,'nws','T_fav_pred_minus_T_fav.png');

T_diff_pred_fav = T_fav_pred0-T_fav_half;
T_diff_pred_fav = T_diff_pred_fav(1:548,:);
plot2DDataC((T_fav_pred0-T_fav_half),1/728,'comp_50');
% caxis([-3000 1000]);  
toc

% T_diff_pred_fav = T_fav_pred0-T_fav_half;
% T_diff_pred_fav = T_diff_pred_fav(1:548,:);
% printResultPNG(T_diff_pred_fav,MASK,-3300,2000,'nws','T_fav_pred_minus_T_fav.png');