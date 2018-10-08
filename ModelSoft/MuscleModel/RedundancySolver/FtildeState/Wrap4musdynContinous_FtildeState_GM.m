function phaseout = Wrap4musdynContinous_FtildeState_GM(input)

persistent splinestruct

if isempty(splinestruct) || size(splinestruct.MA,1) ~= length(input.phase.time) 
    splinestruct = SplineInputData_GM(input.phase.time,input);
end

input.auxdata.splinestruct = splinestruct;

phaseout = musdynContinous_FtildeState_GM_damped(input);