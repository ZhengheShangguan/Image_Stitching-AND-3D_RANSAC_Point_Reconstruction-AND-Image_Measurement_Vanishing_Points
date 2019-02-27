function create_gif(X, center, avr_proj_Res)
% this function is for creating a GIF based on different views


% Display the camera center and reconstructed points in 3D, report the residuals
for i = 1:11
    
    view_vec = [15*(i-1), 2*(i-1)];
    h = figure('units','normalized','outerposition',[0 0 1 1]);
    hold on
    plot3(X(:,1),X(:,2),X(:,3),'b.');
    plot3(center(1),center(2),center(3),'r*');
    axis equal
    view(view_vec)
    legend(['Projected Residual\_image02 = ' num2str(avr_proj_Res)]);
    xlabel('X-axis')
    ylabel('Y-axis')
    zlabel('Z-axis')
    title('Reconstructed 3D points for image 02')

    drawnow 

    % Capture the plot as an image 
    frame = getframe(h); 
    im = frame2im(frame); 
    [imind,cm] = rgb2ind(im,256); 

    % Write to the GIF File 
    if i == 1 
        imwrite(imind,cm,'image.gif','gif', 'Loopcount',inf); 
    else 
        imwrite(imind,cm,'image.gif','gif','WriteMode','append'); 
    end

end


end