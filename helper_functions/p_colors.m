function phase_color = p_colors(palette)
% function phase_color = p_colors(palette)
switch palette
    case 'yellow_red'
        
        phase_color{1}=[255,201,70]./255;
        
        phase_color{2}=[253,141,33]./255;
        %phase_color = 'r';
        
        phase_color{3}=[227,26,28]./255;
        %phase_color = 'g';
        
        phase_color{4}=[142,23,15]./255;
        %phase_color = 'k';
        
    case 'green_blue'
        phase_color={[161 218 180]./255; [65 182 196]./255; [34 94 168]./255; [10 30 69]./255};
    case 'uni_blue'
        phase_color={[239,243,255]./255; [189,215,231]./255; [107,174,214]./255; [33,113,181]./255};
        % HSL 225, 100%, 97
    case 'myblue'
        phase_color=[107,174,214]./255;
     case 'mygray'
        phase_color=[0.3 0.3 0.3];
end
