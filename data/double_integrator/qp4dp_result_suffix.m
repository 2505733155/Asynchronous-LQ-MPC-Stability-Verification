function suffix = qp4dp_result_suffix(requested)
%QP4DP_RESULT_SUFFIX Keep corrected-QP outputs separate from nominal records.
suffix = string(requested);
if strlength(suffix)==0
    suffix = "_qp4dp";
elseif ~contains(suffix,"qp4dp")
    suffix = suffix + "_qp4dp";
end
end
