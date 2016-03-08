%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%	
%	This file is part of the Snoopy code.

%   Snoopy code is free software: you can redistribute it and/or modify
%    it under the terms of the GNU General Public License as published by
%    the Free Software Foundation, either version 3 of the License, or
%    (at your option) any later version.

%    Snoopy code is distributed in the hope that it will be useful,
%    but WITHOUT ANY WARRANTY; without even the implied warranty of
%    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%    GNU General Public License for more details.

%    You should have received a copy of the GNU General Public License
%    along with Snoopy code.  If not, see <http://www.gnu.org/licenses/>.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% Get the data 

rep='../';
full_file=importdata('../rates_r00.txt');

timevar=full_file.data;
nvar=size(timevar,2);
var_name=strread(full_file.textdata{1},'%s',nvar);

% Assign variables to columns of data 
for i=1:nvar
   assignin('base',var_name{i},timevar(:,i)); 
end

% Create a table with all variables
tblA = table(KX,KY,KZ,TIME,EM,GR);

% Sort the rows of the table
tblB = sortrows(tblA);

% Find how many different modes are present
% @mode_box - array that contains the location of last mode of the group 
% @n_modes  - number of various dynamo modes
i=0; 
k=0;
mode_box{1}=1;
for j=1:(size(KX,1)-1)
    if (tblB{j,1}==tblB{j+1,1} &tblB{j,2}==tblB{j+1,2} & tblB{j,3}==tblB{j+1,3})
        i=i+1;
    else
        k = k+1;
        mode_box{k+1}=j+1;
    end
end
mode_box{k+2}=size(KX,1);
n_modes = size(KX,1) - i; 

% Look at each mode and determine whether its growing exponentially 
% and if it's exserting an oscillating bahaviour
for i=1:n_modes
    % Check that a mode has data at more than one time 
    % for such modes, take the growth rate from Snoopy
   if (mode_box{i+1}-mode_box{i} == 1)
        gr_arr{i} = tblB{i,6};
        chi_arr{i} = blanks(1);
   else
       % Get the EM and time data to fit an exp
       for j=1:(mode_box{i+1}-mode_box{i})
           x{j} = tblB{j+mode_box{i}-1,4};
           y{j} = tblB{j+mode_box{i}-1,5};
       end
       x_arr = cell2mat(x);
       y_arr = cell2mat(y);
       
       if (i==2) 
            figure(3)
            semilogy(x_arr, y_arr);
            legend('I=2');  
       end
        
       % Get ln(EM) and convert cells into arrays
       for k=1:size(y,2)
            if (y{k} >= 0.0) 
                ln_y{k} = log(y{k}); 
            else
                'ERROR: negative energy'
            end   
       end
       x_arr = cell2mat(x);
       y_arr = cell2mat(y);
       ln_y_arr = cell2mat(ln_y);
       
       % Linear fit to TIME vs ln(EM), real growth rate is given 
       % by the slope of the linear curve
       f = polyfit(x_arr,ln_y_arr,1);
       gr_arr{i} = f(1);
       
       % Find how good of a fit the linear curve is 
       % by looking at R^2, which predicts % variance in EM 
       yfit = polyval(f,x_arr);
       yresid = y_arr - yfit;
       SSresid = sum(yresid.^2);
       SStotal = (length(y_arr)-1) * var(y_arr);
       RSq{i} = 1 - SSresid/SStotal;
       
       
       % Try linear fit using predicted function
       rng default                    % for reproducibility
       fun= @(x00)x00(1) + x00(2)*exp(x00(3)*x_arr)-y_arr;
       x0 = [0,0,0];
       options = optimoptions(@lsqnonlin,'Algorithm','levenberg-marquardt');
       x00 = lsqnonlin(fun,x0,[],[],options);
       %plot(x_arr,y_arr, x_arr,vestimated(1)+vestimated(2)*exp(cgr{i}*x_arr));
       lin_gr{i}=x00(3); 
       
       figure(2)
       if (i==1) 
       semilogy(x_arr, y_arr, x_arr, x00(1) + x00(2)*exp(x00(3)*x_arr));
       legend('simulation',num2str(i))  
       end
       hold on
       if (i==2 | i==3)
       semilogy(x_arr, y_arr, x_arr, x00(1) + x00(2)*exp(x00(3)*x_arr));
       legend('simulation',num2str(i))
       end
       
       % If R^2 > 1, the growth rate is most likely complex 
       % Trigonometric fit would provide us with the period 
       rng default                    % for reproducibility
       % Create a comlplex exponential model. The model computes 
       % a vector of differences between predicted values 
       % and observed values.
       objfcn = @(v)v(1)*exp((v(2)+sqrt(-1)*v(3))*x_arr)-y_arr;
       % Because the data is complex, 
       % set the Algorithm option to 'levenberg-marquardt'
       opts = optimoptions(@lsqnonlin,...
        'Algorithm','levenberg-marquardt','Display','off');
       x0 = (1+sqrt(-1))*[1;1;1]; % arbitrary initial guess
       [vestimated,resnorm,residuals,exitflag,output] = lsqnonlin(objfcn,x0,[],[],opts);
       vestimated,resnorm,exitflag,output.firstorderopt;
       cgr_arr{i}=vestimated(2) + sqrt(-1)*vestimated(3);
       %plot(x_arr,y_arr, x_arr,vestimated(1)+vestimated(2)*exp(cgr{i}*x_arr));
    end
end
% Save the results to the table 
for i=1:n_modes
   kx_arr{i} = tblB{mode_box{i},1}; 
   ky_arr{i} = tblB{mode_box{i},2}; 
   kz_arr{i} = tblB{mode_box{i},3}; 
end

tblC_title={'KX'; 'KY'; 'KZ';'ExpGR'; 'RSq';'LinGR';'CplxGR'};
tblC = table(transpose(kx_arr),transpose(ky_arr),transpose(kz_arr),transpose(gr_arr),transpose(RSq),transpose(lin_gr),transpose(cgr_arr),'VariableNames', transpose(tblC_title));








    