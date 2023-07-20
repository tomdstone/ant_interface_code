function y=nansem(x)
    if size(x, 2) > 1
        error('sem can only take row/column vectors')
    else
        nonanx=x(~isnan(x));
        y=std(nonanx) / sqrt(length(nonanx));
    end
end