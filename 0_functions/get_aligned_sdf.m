function [sdf, raster] = get_aligned_sdf(spikes, alignTimes, ops)

spkTimes = []; spkTimes = round(spikes);
spkTimes(spkTimes == 0) = [];

sdf_session = [];
sdf_session = SpkConvolver (spkTimes, round(max(alignTimes)+5000), ops.sdf_filter);
sdf = nan(length(alignTimes),range(ops.timewin)+1);

for ii = 1:length(alignTimes)
    try
        sdf(ii,:) = sdf_session(alignTimes(ii)+ops.timewin(1):alignTimes(ii)+ops.timewin (end));
        raster{ii,1} = spkTimes(spkTimes>alignTimes(ii)+ops.timewin (1) & spkTimes<alignTimes(ii)+ops.timewin (end))-alignTimes(ii);
    catch
        sdf(ii,:) = nan(1,length(ops.timewin ));
        raster{ii,1} = [];
    end
end

sdf = single(sdf);

end

