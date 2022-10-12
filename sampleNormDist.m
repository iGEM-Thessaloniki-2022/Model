function [sample] = sampleNormDist(Dist_mean, Dist_std, sample_size)

NormalDist = makedist('Normal','mu',Dist_mean,'sigma',Dist_std);

sample = random(NormalDist,sample_size,1);