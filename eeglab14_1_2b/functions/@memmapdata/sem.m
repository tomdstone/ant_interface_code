function y = sem(x)
    if size(x, 2) > 1
        error('sem can only take row/column vectors')
    else
        y = std(x) / sqrt(length(x));
    end
end