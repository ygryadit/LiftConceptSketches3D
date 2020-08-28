function vals_out = combineLabels(vals1, vals2)

maks_average = find(~isnan(vals1) & ~isnan(vals2));
maks_dave = find(~isnan(vals1) & isnan(vals2));
maks_jerry = find(isnan(vals1) & ~isnan(vals2));

vals_out = NaN*ones(size(vals1));
vals_out(maks_average) = 0.5*(vals1(maks_average)+vals2(maks_average));
vals_out(maks_dave) = vals1(maks_dave);
vals_out(maks_jerry) = vals2(maks_jerry);


end