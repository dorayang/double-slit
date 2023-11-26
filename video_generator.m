% This 
function video_generator(file, x, y, t, psi, frame_rate, z_label, title_entry, vid_mode, v)
    figure(2)
    writerObj= VideoWriter(file);
    writerObj.FrameRate = frame_rate;
    
%     % For removing values before the slit
%     psi(:,:, 1:round((length(y) - 1)/4) + 1) = zeros(length(t), length(x), round((length(y) - 1)/4) + 1);
% 
%     % For removing values outside of potential well
%     psi(:,:, 1:round(0.55*length(y))) = zeros(length(t), length(x),round(0.55*length(y)));
%     psi(:,:, round(0.75*length(y)):end) = zeros(length(t), length(x),round(length(y) - 0.75*length(y))+1);
%     psi(:,1:round(0.4*length(x)), :) = zeros(length(t), round(0.4*length(x)), length(y));
%     psi(:,round(0.6*length(x)):end, :) = zeros(length(t), round(length(x)- 0.6*length(x))+1, length(y));

    
    % open the video writer
    open(writerObj);
    % write the frames to the video
    for i=1:2:length(t)

        S = strcat("frame ", string(i), " of ", string(length(t)));
        disp(S);

        if vid_mode == 1 % Surface plots
            
            surf(x, y, v + (reshape(psi(i,:,:), [length(x), length(y)])));
            colormap turbo
            shading interp
            zlim([min(psi,[],'all')-1,max(psi,[],'all') + 1])
            xlabel("y")
            ylabel("x")
            zlabel(z_label)
            title(title_entry)
            legend(strcat("Time: ", string(compose("%.5f", t(i)))), 'location', 'northoutside')
            writeVideo(writerObj, getframe(gcf));
        elseif vid_mode == 2 % Contour plots
             contourf(x, y, (reshape(psi(i,:,:), [length(x), length(y)])), 'edgecolor','none', 'LevelStep',(max(psi,[],'all') -  min(psi,[],'all'))/100);           
             colormap turbo
             c = colorbar;
             c.Label.String = "\surd{(\psi \psi^*)}";
             clim([min(psi,[],'all'), max(psi,[],'all')]) 
             xlabel("y")
             ylabel("x")
             title(title_entry)   
             annotation('textbox', [0.13, 0.07, 0.1, 0.1], 'String', strcat("Time: ", string(compose("%.5f", t(i)))), 'BackgroundColor','w')
             writeVideo(writerObj, getframe(gcf));
        end
    end
end