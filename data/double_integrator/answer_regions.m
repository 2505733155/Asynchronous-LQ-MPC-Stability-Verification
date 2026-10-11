function [reg, law] = answer_regions()
%ANSWER_REGIONS Input-constrained explicit MPC regions for this example.

persistent cached_reg cached_law
if isempty(cached_reg)
    cached_reg = answer_regions_input_only_legacy();
    cached_law = struct('mode', "input_only_legacy");
end
reg = cached_reg;
if nargout > 1
    law = cached_law;
end
end
