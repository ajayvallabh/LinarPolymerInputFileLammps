function [randomv] = random(iseed)
%UNTITLED5 Summary of this function goes here
%   Detailed explanation goes here
aa = 16807;
mm = 2147483647;
sseed = iseed;
sseed = mod((aa*sseed),mm)
randomv = sseed/mm;
iseed=sseed;
end

