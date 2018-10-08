function phaseout = musdynContinous_FtildeState_GM_dampedHesWrap(input)

persistent splinestruct

if isempty(splinestruct) || size(splinestruct.MA,1) ~= length(input.phase.time.f)
    splinestruct = SplineInputData_GM(input.phase.time.f,input);
end

input.auxdata.splinestruct = splinestruct;

phaseout = musdynContinous_FtildeState_GM_dampedADiGatorHes(input);