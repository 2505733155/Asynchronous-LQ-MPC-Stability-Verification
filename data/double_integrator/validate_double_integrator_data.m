function summary = validate_double_integrator_data()
%VALIDATE_DOUBLE_INTEGRATOR_DATA Validate the public seven-region PWA data.
% The routine performs inexpensive structural checks before any solver call.

reg = answer_regions();
required = {'K', 'k', 'H', 'b'};
expectedRegionCount = 7;

if numel(reg) ~= expectedRegionCount
    error('Expected %d PWA regions, found %d.', expectedRegionCount, numel(reg));
end

originHits = 0;
for i = 1:numel(reg)
    for field = required
        if ~isfield(reg(i), field{1})
            error('Region %d is missing the required field "%s".', i, field{1});
        end
    end

    K = reg(i).K;
    k = reg(i).k;
    H = reg(i).H;
    b = reg(i).b;
    if ~isnumeric(K) || ~isequal(size(K), [1, 2]) || any(~isfinite(K), 'all')
        error('Region %d must have a finite 1-by-2 feedback gain K.', i);
    end
    if ~isnumeric(k) || ~isscalar(k) || ~isfinite(k)
        error('Region %d must have a finite scalar offset k.', i);
    end
    if ~isnumeric(H) || size(H, 2) ~= 2 || isempty(H) || any(~isfinite(H), 'all')
        error('Region %d must have a nonempty finite H matrix with two columns.', i);
    end
    if ~isnumeric(b) || ~isequal(numel(b), size(H, 1)) || any(~isfinite(b), 'all')
        error('Region %d must have a finite b vector with one entry per row of H.', i);
    end
    if all(H * [0; 0] <= b(:) + 1e-10)
        originHits = originHits + 1;
    end
end

if originHits == 0
    error('The origin is not covered by the public PWA partition.');
end

summary = struct('region_count', numel(reg), 'origin_region_count', originHits);
fprintf('Validated public double-integrator PWA data: %d regions; origin covered by %d region(s).\n', ...
    summary.region_count, summary.origin_region_count);
end
