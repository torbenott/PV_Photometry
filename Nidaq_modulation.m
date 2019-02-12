function Modulated_LED=Nidaq_modulation(amp,freq,TaskParameters)
%Generates a sin wave for LED amplitude modulation.

duration = TaskParameters.GUI.NidaqDuration;
sample_rate = TaskParameters.GUI.NidaqSamplingRate;

if TaskParameters.GUI.Modulation==0 || freq==0 || amp==0
    Modulated_LED=(amp/2)*ones(duration*sample_rate,1);
else
DeltaT=1/sample_rate;
Time=0:DeltaT:(duration-DeltaT);
Modulated_LED=amp*(sin(2*pi*freq*Time)+1)/2;
Modulated_LED=Modulated_LED';
end
end