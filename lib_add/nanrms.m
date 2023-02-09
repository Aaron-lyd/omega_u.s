function y = nanrms(x)

y = sqrt(mean(x .* conj(x), 'omitnan'));