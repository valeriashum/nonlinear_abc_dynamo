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


% Example File: how to read the timevar file in Matlab
% This example reads all the variable in the timevar file
% and display the time history of several quantities...

rep='../';


mean_start=100;     % Where do we start our time-averages?

full_file=importdata(['../20x2x2_R_0plus/short_timevar_200_10']);
timevar=transpose(full_file.data);
nvar=size(timevar,1);
var_name=strread(full_file.textdata{1},'%s',nvar);

for i=1:nvar
   assignin('base',var_name{i},timevar(i,:)); 
end

%mag.helicity=timevar(22,:);

% Only when braginskii pressure is present
%brag.p=timevar(23,:);
%brag.stress=timevar(24,:);


%mode1.kx=timevar(23,:);
%mode1.wz=timevar(24,:);
%mode2.kx=timevar(25,:);
%mode2.wz=timevar(26,:);

%vel.shear=timevar(20,:);
%vel.vxmode=timevar(20,:);
%vel.vymode=timevar(21,:);


if exist('dvxvy') 
    vxvy=dvxvy;
end

if exist('vxvy')&&exist('bxby')

    alpha=vxvy-bxby;

    figure(8)
    plot(t,alpha,t,bxby,t,vxvy);
    legend('alpha coefficient','Maxwell Stress','Reynolds Stress');
    title('Turbulent Transport')
    
    clear mean_transport;
    
    if length(vxvy) > mean_start
        mean_vel_transport=cumsum(vxvy(mean_start:end))./(1:length(vxvy(mean_start:end)));
        mean_mag_transport=cumsum(bxby(mean_start:end))./(1:length(bxby(mean_start:end)));
        mean_tot_transport=cumsum(alpha(mean_start:end))./(1:length(alpha(mean_start:end)));

        e=mean(ev(mean_start:end))
        alpha_v=mean_vel_transport(end);
        alpha_b=mean_mag_transport(end);
        alpha=alpha_v-alpha_b
    else
        disp(['Caution, we are not removing initial conditions']);
        e=mean(ev)
        alpha_v=mean(vxvy);
        alpha_b=-mean(bxby);
        alpha=alpha_v-alpha_b
    end


elseif exist('vxvy')
    alpha=vxvy;

    figure(8)
    semilogy(t,alpha);
    title('Turbulent Transport')
    
    clear mean_transport;
    
    if length(vxvy) > mean_start
        mean_vel_transport=cumsum(vxvy(mean_start:end))./(1:length(vxvy(mean_start:end)));
        mean_tot_transport=cumsum(alpha(mean_start:end))./(1:length(alpha(mean_start:end)));

        e=mean(ev(mean_start:end))
        alpha_v=mean_vel_transport(end);
        alpha=alpha_v
    else
        disp(['Caution, we are not removing initial conditions']);
        e=mean(ev)
        alpha_v=mean(vxvy);
        alpha=alpha_v
    end
    
end


if exist('em')
    set(gca,'FontSize',30,'fontWeight','bold')
    set(findall(gcf,'type','text'),'FontSize',30,'fontWeight','bold')

    figure(2)
    semilogy(t,em)
    title('Magnetic energy')
    xlabel('Time')
    
    figure(20)
    plot(t,log(em))
    title('Magnetic energy')
end

if exist('em_gr')&&exist('Rm')
    figure(26)
    N=1/200.;
    %xx = 3.512407e-02:2.5e-04:2.883686e+01;                              %for N=1/40 
    %xx = 1.0006231e-01:0.000025:1.229837e-01;                            %for N=1/80 2
    %xx = 6.717514e-02:2.5e-04:1.414214e-01;                              %for N=1/200
    xx = 2.124265e-01:2.5e-04:2.906888e-01;                              %for N=1/200 10
    yy = spline(Rm,em_gr,xx);                                             
    %subplot(2,1,1), plot(Rm,em_gr,'o',xx,yy, xx, N.*(xx-N./xx));         %for N=1/40 
    plot(Rm,em_gr,'o',xx,yy, xx, N.*(xx-N./xx),xx,0.0);                  %for N=1/200
    xlabel('1/\eta','fontsize',16);
    ylabel('Growth Rate','fontsize',16);
    set(gca,'fontsize',16);
    title('Box: 400\pi x 10\pi x 10\pi','fontsize',16)
     
    
    %subplot(2,1,2), plot(Rm,em_gr,'o',xx,yy, xx,N.*(xx-N./xx),xx,0.0);   %for N=1/40 
    %xlim([0.0 2.0]);                                                      %for N=1/40 
    %ylim([-0.2235 0.05]);                                                 %for N=1/40 
    %subplot(2,1,2), plot(Rm,em_gr,'o',xx,yy, xx,N.*(xx-N./xx),xx,0.0);   %for N=1/40 
    %xlim([0.11 0.12]);                                                   %for N=1/80 
    
    legend('SNOOPY','SNOOPY interpolation','Asymptotic Approx','Orientation','horizontal' );
    legend boxoff;                                                
    %xlabel('1/\eta','fontsize',16);
    %ylabel('Growth Rate','fontsize',16);
    %set(gca,'fontsize',14);
end

if exist('k_par')&&exist('k_BxJ')&&exist('k_alfvenic')
    set(gca,'FontSize',15,'fontWeight','bold')
    set(findall(gcf,'type','text'),'FontSize',15,'fontWeight','bold')

    figure(25)
    subplot(2,1,1),  plot(t,k_par,t,k_BxJ ,t,k_alfvenic,'LineWidth',2 );legend('k_{||}/2\pi','k_{\perp}/2\pi','k_{alfvenic}/2\pi');title('Characteristic Wavenumbers');
    subplot(2,1,2),  plot(t,k_par./k_BxJ,t,1., 'LineWidth',2 );legend('k_{||}/k_{\perp}')
    xlabel('Time')
end
if exist('Bx_kxdom')&&exist('Bx_kydom')&&exist('Bx_kzdom')&&exist('By_kxdom')&&exist('By_kydom')&&exist('By_kzdom')&&exist('Bz_kxdom')&&exist('Bz_kydom')&&exist('Bz_kzdom')
    figure(21)
    subplot(3,1,1), plot(t,Bx_kxdom,t,Bx_kydom,t,Bx_kzdom, t, 2.);
    title('Magnetic field B_x: dominant wavenumber k_x,k_y,k_z ');
    legend('k_x','k_y','k_z','k_u','Orientation','horizontal');
    subplot(3,1,2), plot(t,By_kxdom,t,By_kydom,t,By_kzdom ,t, 2.) ;title('Magnetic field B_y: dominant wavenumber k_x,k_y,k_z '); 
    subplot(3,1,3), plot(t,Bz_kxdom,t,Bz_kydom,t,Bz_kzdom, t, 2.) ;title('Magnetic field B_z: dominant wavenumber k_x,k_y,k_z '); xlabel('time');
   
    %t vs k^2
    figure(22)
    plot(t,sqrt(Bx_kxdom.^2 + Bx_kydom.^2 + Bx_kzdom.^2),'d', t,sqrt(By_kxdom.^2 + By_kydom.^2 + By_kzdom.^2),'d',t ,sqrt(Bz_kxdom.^2 + Bz_kydom.^2 + Bz_kzdom.^2),'d', t, sqrt(2));
    legend('B_x','B_y','B_z', 'u_{ABC}');
    title('Time vs squareroot of |k^2|'); xlabel('time');ylabel('sqrt(k^2)');   
end

if exist('B_kxdom')&&exist('B_kydom')&&exist('B_kzdom')
    figure(21)
    subplot(2,1,1), plot(t,B_kxdom,'d',t,B_kydom,'d',t,B_kzdom,'d', t, 1.);
    title('Magnetic field <B>: dominant wavenumber k_x,k_y,k_z ');
    legend('k_x','k_y','k_z','k_u=1','Orientation','horizontal');
    subplot(2,1,2), plot(t,sqrt(B_kxdom.^2 +B_kydom.^2+B_kzdom.^2),'d', t, sqrt(2.)) ;title('Magnetic field <B>: dominant wavenumber k '); xlabel('time');
    legend('k','k_u','Orientation','horizontal');
end


if exist('vxmax')&&exist('vxmin')&&exist('vymax')&&exist('vymin')&&exist('vzmax')&&exist('vzmin')
    figure(3)
    subplot(3,1,1), plot(t,vxmax,t,vxmin);
    title('Velocity field; extrema vx,vy,vz');
    subplot(3,1,2), plot(t,vymax,t,vymin);
    subplot(3,1,3), plot(t,vzmax,t,vzmin);
end

if exist('bxmax')&&exist('bxmin')&&exist('bymax')&&exist('bymin')&&exist('bzmax')&&exist('bzmin')
    figure(4)
    subplot(3,1,1), plot(t,bxmax,t,bxmin);
    title('Magnetic field; extrema bx,by,bz');
    subplot(3,1,2), plot(t,bymax,t,bymin);
    subplot(3,1,3), plot(t,bzmax,t,bzmin);
end


if exist('j2')&&exist('w2')
    figure(11)
    semilogy(t,w2,t,j2,'linewidth',2);
    legend('Enstrophy','J^2')
    xlabel('t','fontsize',16)
    set(gca,'fontsize',16)
    
elseif exist('w2')
    figure(11)
    semilogy(t,w2,'linewidth',2);
    legend('Enstrophy')
    xlabel('t','fontsize',16)
    set(gca,'fontsize',16)
end

if exist('hm')
    figure(12)
    plot(t,hm,'linewidth',2);
    title('Magnetic helicity')
    xlabel('t','fontsize',16)
    set(gca,'fontsize',16)
end


if exist('dmax')&&exist('dmin')
    figure(13)
    plot(t,dmax,t,dmin);
    title('Density field extrema');
end

if exist('az2')
    figure(14)
    plot(t,az2);
    title('potential vector')
end

if exist('vxaz')
        figure(15)
        plot(t,vxaz);
        title('potential vector source')
end

if exist('hc')
    figure(16)
    plot(t,hc);
    title('Cross helicity')
end
