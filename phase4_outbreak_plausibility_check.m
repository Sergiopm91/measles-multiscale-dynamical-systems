function T = phase4_outbreak_plausibility_check()
% PHASE4_OUTBREAK_PLAUSIBILITY_CHECK
% Soft plausibility check against documented measles outbreak regimes.
%
% No clinical calibration is performed here. The table only places modeled
% high-susceptibility and high-coverage suppression regimes next to outbreak
% descriptors to support qualitative plausibility in the manuscript.

if ~exist('export','dir'), mkdir('export'); end
Outbreak = {'Samoa 2019'; 'Madagascar 2018-2019'; 'Disneyland 2014-2015'};
Use_in_paper = {'high-susceptibility severe outbreak plausibility'; ...
                'large outbreak in low-coverage context'; ...
                'localized outbreak despite high national coverage'};
Required_citation_placeholder = {'[Samoa 2019 official/WHO source]'; ...
                                 '[Madagascar 2018-2019 WHO/CDC source]'; ...
                                 '[Disneyland 2014-2015 CDC/MMWR source]'};
Model_comparison = {'Use 0-60% coverage regimes; compare large attack-rate potential only'; ...
                    'Use low/intermediate coverage regimes; qualitative size only'; ...
                    'Use clustered susceptibility/high average coverage; qualitative persistence only'};
Claim_allowed = {'plausibility check, not calibration'; 'plausibility check, not calibration'; 'plausibility check, not calibration'};
T = table(Outbreak, Use_in_paper, Required_citation_placeholder, Model_comparison, Claim_allowed);
writetable(T, fullfile('export','phase4_outbreak_plausibility_check.csv'));
fprintf('\nOutbreak plausibility check exported: export/phase4_outbreak_plausibility_check.csv\n');
end
