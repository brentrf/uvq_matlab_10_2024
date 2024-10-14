%Animation
function MAIN()
    clear;
    close all;
    thetas = [-90:90]';
 
    %INPUTS  ........
    Rmax            = 2.5;     	%[m] Irradiation pattern range 
    pattern.FWHM    = 100;      %[deg] FWHM of pattern
    head.radius     = 0.25;    	%[m] 
    Yclearance     = 1.5;      	%[m]  head_center to ceiling
    R_competitor    = Yclearance - head.radius;

    %Pattern Width LUT
    LUT.alphas = [0.5, 1,  2,  3,  4,  6, 10];
    LUT.HWHM = [75, 60, 45, 38, 33, 27, 21];
    figure; plot(LUT.alphas,LUT.HWHM);    xlabel('cospower'); ylabel('HWHM [deg]'); title('LUT')
    
    % LOOK UP TABLE FOR PATTERN
    pattern.alpha = interp1(LUT.HWHM,LUT.alphas,pattern.FWHM/2);
    y_temp = cosd(thetas);
    figure; polarplot(deg2rad(thetas),y_temp);
    hold on; polarplot(deg2rad(thetas),y_temp.^2);
    hold on; polarplot(deg2rad(thetas),y_temp.^3);
    hold on; polarplot(deg2rad(thetas),y_temp.^5);
    legend('60','45','38','33');
    thetalim([-90,90]);

    figure; plot(thetas,y_temp);
    hold on; plot(thetas,y_temp.^2);
    hold on; plot(thetas,y_temp.^3);
    hold on; plot(thetas,y_temp.^4);
    hold on; plot(thetas,y_temp.^6);
    hold on; plot(thetas,y_temp.^10);
    hold on; plot(thetas,y_temp.^0.5);
    legend('60','45','38','33','27','21','75');

    %====  PLOT PATTERNS =======
    hf1=figure(); polarplot(1,1); hf2=figure();polarplot(1,1);   hf3=figure();polarplot(1,1); 
    betas = [30,45,60,75]';
    R=0.95;
    for kb=1:length(betas)
        legnames{kb} = [num2str( betas(kb) ),'deg'];

        pattern.beta_deg = betas(kb);
        pattern.alpha    = interp1(LUT.HWHM,LUT.alphas,pattern.beta_deg);
        R = (1/(1-cosd(pattern.beta_deg)))^0.5;
        [thetas1, rho_out1] = generate_surface_Uniform(thetas,R,pattern);
        R = (1/(0.5*(sind(pattern.beta_deg)^2)))^0.5;
        [thetas2, rho_out2] = generate_surface_LambTrunc(thetas,R,pattern);    
        R = (1+pattern.alpha)^0.5;
        [thetas3, rho_out3] = generate_surface(thetas,R,pattern); 
        
        figure(hf1);  polarplot(deg2rad(thetas1),rho_out1); hold all; thetalim([-90,90]);
        figure(hf2);  polarplot(deg2rad(thetas2),rho_out2);  hold all;   thetalim([-90,90]);     
        figure(hf3);  polarplot(deg2rad(thetas2),rho_out3); hold all; thetalim([-90,90]);
    end 
    figure(hf1); legend(legnames);
    figure(hf2); legend(legnames);
    figure(hf3); legend(legnames);    
    
    %================ generate animation

    %define person/object
    persons_headctr_x = [-2:0.05:+2]';                  %[m]  motion vector (x positions)
    Nsteps = length(persons_headctr_x);
    persons_headctr_y = ones(Nsteps,1)*(-Yclearance);   %[m]  motion vector (y positions)

    %set up video
    v = VideoWriter('test.avi');
    open(v);

    %loop over frames
    hf = figure;
    for ks = 1:Nsteps
        
        %generate & place person
        head.x = persons_headctr_x(ks,1);
        head.y = persons_headctr_y(ks,1);
        head.theta = atan2d(head.y,head.x);
        poly_points = generate_polygon_person(head);    %generate points for "person"
        pgon = polyshape(poly_points.x,poly_points.y);  %polygon

        %Generate Surface of Max Allowed Intensity
        cosd(head.theta-90);
        R = -head.radius*0.2+min([Rmax, ...
                ((head.x^2 + head.y^2)^0.5 - head.radius ) / (cosd(head.theta+90)^pattern.alpha) ]);
        
        [thetas, rho_out] = generate_surface(thetas,R,pattern);
        [thetas, rho_out_competition] = generate_surface(thetas,R_competitor,pattern);
        
        figure(hf); lamp.x = rho_out.*cosd(thetas-90);  lamp.y = rho_out.*sind(thetas-90);
        plot(lamp.x,lamp.y); hold all;
        lamp_competition.x = rho_out_competition.*cosd(thetas-90); 
        lamp_competition.y = rho_out_competition.*sind(thetas-90); 
        plot(lamp_competition.x,lamp_competition.y,'--','Color','r'); 
        hold on; plot(2*persons_headctr_x,zeros(Nsteps,1),'Linewidth',3','color','k');
        hold on; plot(2*persons_headctr_x,ones(Nsteps,1)*(-Yclearance-5.5*head.radius),'Linewidth',2','color',[0.2,0.2,0.2]);
        xlim([ min(persons_headctr_x)-head.radius*2 , max(persons_headctr_x)+head.radius*2]);
        ylim([ -Yclearance-5*head.radius, 0]); axis equal;
        hold on; plot(pgon);  %add person
        hold off;
        title('Active Sensing allows *higher efficacy* while maintaining ACGIH TLVs');
        
        movie_F = getframe(hf);

        writeVideo(v,movie_F);

        pause(0.15);
        
        
        %polar plot
        %         figure(hf); polarplot(deg2rad(thetas),rho_out); rlim([0,10]);
        %         ax = gca; ax.ThetaZeroLocation = 'bottom'; thetalim([-90,90]);
    end
    close(v);

    
end


function poly_points = generate_polygon_person(head)

    rhead = head.radius;  %[m]    
    poly.x1 = [rhead*cosd([-85:5:265])]';
    poly.y1 = [rhead*sind([-85:5:265])]';
    poly.x2 = -[0.2*rhead,rhead*1.5,rhead*1.5,rhead*1.25,rhead*1.25,...
               -rhead*1.25,-rhead*1.25,-rhead*1.5,-rhead*1.5,-0.2*rhead]';
    poly.y2 = -1*(rhead+ [0.2*rhead,0.2*rhead, 2*rhead, 2*rhead,4.5*rhead,...
                   4.5*rhead,2*rhead,2*rhead,0.2*rhead,0.2*rhead]' );

    poly_points.x = head.x + [poly.x1;poly.x2];
    poly_points.y = head.y + [poly.y1;poly.y2];
end


function [thetas, rho_out] = generate_surface(thetas,R,pattern)
    
    rho_out = R * cosd(thetas).^pattern.alpha;
    
    
end


function [thetas, rho_out] = generate_surface_LambTrunc(thetas,R,pattern)
    
    rho_out = R * cosd(thetas).*(abs(thetas)<=pattern.beta_deg) ;
    
end

function [thetas, rho_out] = generate_surface_Uniform(thetas,R,pattern)
    
    rho_out = R .*(abs(thetas)<=pattern.beta_deg) ;
    
end


