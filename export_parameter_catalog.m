function T = export_parameter_catalog(filename)
% EXPORT_PARAMETER_CATALOG  Export parameter source/type table for the paper.
%
%   T = export_parameter_catalog()
%   T = export_parameter_catalog('export/phase_parameters_by_source.csv')

if nargin < 1 || isempty(filename)
    if ~exist('export','dir'), mkdir('export'); end
    filename = fullfile('export','phase_parameters_by_source.csv');
end

P1 = measles_params(1);
fields = fieldnames(P1);
rows = {};
for i = 1:numel(fields)
    nm = fields{i};
    if startsWith(nm,'param_') || any(strcmp(nm, {'name','profile_type','profile_interpretation','archetype_modifier_note'}))
        continue;
    end
    val = P1.(nm);
    if ~(isnumeric(val) || islogical(val) || ischar(val) || isstring(val))
        continue;
    end
    source_type = get_meta(P1.param_class, nm, 'not_classified');
    source_label = get_meta(P1.param_source_label, nm, 'not classified');
    value_reference = stringify_value(val);

    changed = false(1,6);
    values = strings(1,6);
    values(1) = string(value_reference);
    for w = 2:6
        Pw = measles_params(w);
        if isfield(Pw,nm)
            values(w) = string(stringify_value(Pw.(nm)));
            changed(w) = ~isequaln(Pw.(nm), val);
        else
            values(w) = "NA";
        end
    end
    if any(changed(2:6))
        paper_class = 'exploratory_archetype_modifier';
    else
        paper_class = source_type;
    end

    rows(end+1,:) = {nm, paper_class, source_type, source_label, value_reference, ...
        char(values(2)), char(values(3)), char(values(4)), char(values(5)), char(values(6))}; %#ok<AGROW>
end

T = cell2table(rows, 'VariableNames', {'Parameter','Paper_classification','Baseline_source_type', ...
    'Source_or_role','Healthy_adult_value','Child_value','Elderly_value', ...
    'Immunocompromised_value','Partial_vaccine_value','Fully_vaccinated_value'});
writetable(T, filename);
fprintf('  Saved parameter catalog: %s\n', filename);
end

function s = get_meta(meta, field, fallback)
if isfield(meta, field), s = meta.(field); else, s = fallback; end
end

function s = stringify_value(v)
if isnumeric(v) || islogical(v)
    if isscalar(v)
        s = sprintf('%.8g', double(v));
    else
        s = mat2str(v);
    end
elseif isstring(v) || ischar(v)
    s = char(v);
else
    s = '<non-scalar>';
end
end
