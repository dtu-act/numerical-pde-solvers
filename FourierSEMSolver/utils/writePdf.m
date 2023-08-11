function writePdf(fig, path)    
    set(fig,'Units','centimeters');
    screenposition = get(fig,'Position');
    set(fig,...
        'PaperPosition',[0 0 screenposition(3:4)],...
        'PaperSize', [screenposition(3:4)]);
    print(fig, '-dpdf', '-painters', path)
end