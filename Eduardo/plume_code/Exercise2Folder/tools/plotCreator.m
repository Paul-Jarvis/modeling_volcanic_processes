function plotCreator(x, y, plot_xlabel, plot_ylabel, plot_title, color, style, hold_idx)
    
    global plot_counter
    
    plot_width = 3;
    plot_fontsize = 14;
    
    if (hold_idx==0)
        figure1 = figure(plot_counter);
        plot_counter = plot_counter+1;
        axes1 = axes('Parent',figure1,'FontSize',plot_fontsize);
        box(axes1,'on');
        hold(axes1,'all');
    end
    
    plot(x,y,'Color',color,'LineStyle',style,'LineWidth',plot_width);
    xlabel(plot_xlabel,'FontSize',plot_fontsize);
    ylabel(plot_ylabel,'FontSize',plot_fontsize);
    title(plot_title,'FontSize',plot_fontsize);
    
    if (hold_idx>0)
        hold off;
    end    

end