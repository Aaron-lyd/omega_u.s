function var_rms = root_mean_square(var)


var_rms = var - nanmean(var(:));
good = ~isnan(var_rms);
var_rms = rms(var_rms(good));

end


