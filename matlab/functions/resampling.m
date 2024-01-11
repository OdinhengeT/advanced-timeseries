function [particles, weights] = resampling(particles, weights, J)
    resampledIndices = randsample(1:length(weights), length(weights), true, weights);
    if length(resampledIndices) > J
        resampledIndices = resampledIndices(1:J);
    elseif length(resampledIndices) < J
        resampledIndices = randsample(resampledIndices, J, true);
    end
    index = mod(resampledIndices,length(particles));
    index(index==0) = length(particles);
    particles = particles(index,:);
    weights = weights(index);
end
