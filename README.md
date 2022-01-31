# NeuronalWaveformAnalysis

This is the basic code used for the analysis of neuronal electrophysiological recordings in the article in the journal "Brain", titled "The effect of striatal dopaminergic grafts on the neuronal activity in the substantia nigra pars reticulata and subthalamic nucleus in hemiparkinsonian rats".
https://doi.org/10.1093/brain/awr226

Start with "computeMeasures14v4.m", which computes all the measures for one neural recording file (in Spike2 format).  That code was then called repeatedly (using other code not currently included in this repository), to analyze each of the recording files, and the measurements were collected for graphing.  As shown in that paper, analyses were performed such as neuronal firing rate, burst index, coefficient of variation, proportion of spikes in burst, discharge density histogram range, and sample entropy. 

The code is unfortunately rather messy, and it will probably not run without errors at this time, since I haven't used it in quite a few years. Libraries such as the SON library and other libraries may need to be updated. However, this code did run successfully in its time, and provides the basic idea of the analyses used. I am making this available to help other younger researchers get started, as I have been helped by others before me.
