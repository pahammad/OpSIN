**********************************************
OpSIN: Opto-genetic spikesorting and identification of neurons
**********************************************

This package contains various Matlab routines for extracting spike events associated with Channelrhodopson (ChR2) targeted neurons of interest.

For bug reports/comments/feedback, e-mail: parvez@ieee.org

**********************************************
Main routines:
**********************************************

* OpSIN_SS_1Ch_batchWithGuisNew: semi-supervised spikesorting tool

This is the main entry point of the OpSIN algorithm. This tool is a combination of automated processing and manual annotation, to enable fast processing of single channel e-phys recordings. The automated processing allows the user to quickly divide the spike candidate pool into three categories: good/unsure/bad. Once this determination is made, the user is guided through the unsure set of spike candidates to finalize the spike annotation.

**********************************************
Useful sub-routines:
**********************************************

1. Benchmarker_MoMa_1Ch_batch: Mostly manual benchmarking tool

This tool is useful to generate expert-annotated ground-truth for single-channel e-phys recordings from a loose-patch or suction-electrode type preparation. The assumption is that one or a few neurons of interest are tagged with ChR2. It is the slowest of the available options in OpSIN package, but it is also the most objective approach to generate expert-guided ground-truth data.

Several other useful subroutines can be found in the Extra Goodies subfolder.

**********************************************
Example Walkthrough: “Using_OpSIN.m”
**********************************************

This is an example to walk through a case scenario for OpSIN based processing of optogenetic spike-sorting of ephys (electrophysiology) signals.

A) To be done before you run OpSIN

A1) The example datafile '100W_constant_light_stimulation_example_01' contains the recorded electrophysiology data (20kHz sample rate) in column 1, and the voltage trace used for the light stimulation in column 2.

A2) Copy the example datafile into a new folder and specify the directory by the variable name 'generalDataRootDirectory'

A3) You will have to specify the voltage value that was used during the pre-stimulation pulses to aid the spike selection process performed by the algorithm. The variable is called: 'lightChannelONlevel'

B) Run OpSIN

B1) Run'OpSIN_SS_1Ch_batchWithGuis.m' 
A window will open asking you to select the file you want to spike sort. Select the example file '100W_constant_light_stimulation_example_01'.

you will be asked to enter '1' to continue the analysis:
'Enter 1 to continue:' Enter '1' and press return.

B2)
Two graphs will popup as a result of the algorithm looking for possible spike candidates during pre-stimulation pulses.

One graph contains all candidate spike waveforms that were found by the algorithm during light stimulation periods.

The second graph shows candidate spike waveforms clustered together to allow the user to combine relevant clusters that contain the spike waveform of interest (representing ChR2 evoked activity). 

You will be asked to enter the number of CLEAN cluster template IDs to pool together"
'Enter the number of CLEAN cluster template IDs you want to pool together:'
Enter the number of clusters you would like to include and press return.

You will have to select clusters based on visual inspection of the spike shape. Your selection include the cluster containing the one containing the highest number of spike waveforms.

B3)
You will be asked to specify the cluster ID(s) on which you like to base the template matching:
'Enter the BEST clean template waveform ID, ideally one centered on first peak (for double-bump spikes):'

You will be asked whether the chosen cluster contains spikes which were centered on the thin peak or the broad peak:
 For double bumps: was your first template choice centered (~40-50) on the thin peak? (0 = No; 1 = Yes)
Enter '0' and press return.

B4) 

A GUI will pop up (SELECT SPIKE DETECTION THRESHOLD) asking you to specify the spike detection threshold specific to this recording.
Once you adjusted the 'SD section threshold' bar using the mouse, press the 'FINALIZE'
button.

B5)

You will be asked to specify whether the electrophysiology data contains double bump spikes, or not:
'Enter 1 if the data contains double bump spikes, 0 otherwise:'
Enter '0' and press return.

B6)

A GUI will pop up (SELECT THRESHOLDS FOR GOOD, UNSURE & BAD SPIKES), asking you to specify both the 'GOOD' and 'BAD' spike Threshold in the context of a 'Spike Probability Histogram' that indicates the probability distribution of all detected spike candidates (lower right graph). A high probability value reflects a good fit the   spike template.
On the upper right you can also see all 'Estimated spike probabilities' as they are distributed throughout the entire recording. The red trace shows the voltage signal of the light stimulation.

The upper left graph shows all the spikes that - with the currently selected probability threshold - are defined as 'Good Spikes' by the algorithm.

The middle left graph shows all spikes that - with the currently selected probability threshold - are considered as borderline (unsure) spikes. Once that the probability threshold is defined by the user, all spikes that are contained in here will each be manually checked by the user whether they fit the template of a ChR2 evoked spike, or not.

The lower left graph contains all spikes that - with the currently selected probability threshold - are considered 'Bad Spikes' by the algorithm and will accordingly not be included as ChR2 evoked spike candidates.

The probability threshold should be set such as to exclude all spike waveforms that are not in agreement with the template spike waveform. To this end, a thorough visual inspection of all three candidate spike pools is required with the goal to keep the 'Good Spikes' pool (upper left graph) as homogenous as possible. More spikes in the middle panel mean a long user assisted spike selection process.

Once the probability threshold bar by using the mouse press the 'finalize' button.

B7)

In case you kept some spike waveforms in the unsure spike waveform pool you will be asked to decide for every single uncertain spike candidate whether to include it or not.

'Enter your input (1 = highConf good spike, 11= lowConf good spike, 2= highConf bad spike, 22= lowConf bad spike (default), -1=undo):' 

enter '1' if you want to include a spike
enter '2' if you want to discard a spike

B8)

Once this is done for all uncertain spike candidates the algorithm is finished.

The final sorted dataset is saved under the original filename plus '_OpSIN' and the current date.

The final spike data (1 represeting a spike, 0 representing no spike) is saved in the fourth column of the final dataset.

