function phaseout = musdynContinous_FtildeState_GM_dampedGrdWrap(input)

persistent splinestruct

if isempty(splinestruct) || size(splinestruct.MA,1) ~= length(input.phase.time.f) || size(splinestruct.LMT,2)~= input.auxdata.NMuscles
    splinestruct = SplineInputData_GM(input.phase.time.f,input);
end

input.auxdata.splinestruct = splinestruct;

phaseout = musdynContinous_FtildeState_GM_dampedADiGatorGrd(input);