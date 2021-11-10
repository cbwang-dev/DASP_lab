# Answer for lab session 1

## Exercise 1.1

If you have done this correctly, a new mat-file `Computed RIRs.mat` should appear in the current matlab directory. Open this mat-file and check which variables it contains. Try to find out what each variable represents (based on the user-defined values in the simulation environment). The RIRs are stored in the variables `RIR_sources` and `RIR_noise`. Try to figure out what each dimension of these variables represents. How many RIRs have been generated, and why?

> it contains `fs_RIR`, `m_pos`, `rev_time`, `RIR_noise`, `RIR_sources`, `room_dim`, `s_pos`, and`v_pos`. The configuration is saved in `./GUI_setups/1_1.mat`.
>
> > `fs_RIR`: sample frequency, 44100 Hz. 
> > `m_pos`: an array whose size is `[3,2]`, denoting the position of mics. 
> > `rev_time`: T60 reverbration time, 1.5 seconds. 
> > `RIR_noise`: each mic's room impulse response of noise signal(1), size `[22050,3]`.
> > `RIR_sources`: each mic's room impulse response of audio signals(3), size `[22050,3,2]`.
> > `room_dim`: dimension of room, the value is `[10,10]`, meters respectively	.
> > `s_pos`: position of audio sources, size `[2,2]`. 
> > `v_pos`: position of noise source, size `[1,2]`. 
>
> There are overall 9 RIRs generated. `noise(1)->mics(3):3 RIRs` and `audio(2)->mics(3):6 RIRs`.

Create a new figure and plot the RIR from source 2 to microphone 3. Can you recognize the direct path component, the early reflections, and the reverberant tail?

> Use command `plot(RIR_sources(:,3,2))` to generate plot. This is saved in `./figures/1_1_11.bmp`. The direct path component, the early reflections, and the reverberant tail are noted in the file. Further details are provided in course slides: 2-4, of `H09P9A`: *Multimedia Technology and Coding*. 

## Exercise 1.2

What do you conclude in terms of the speech intelligibility in highly reverberant scenarios?

> In highly reverbrant scenarios, the speech tends to be not intelligible. 

Create a new scenario containing a 2-microphone array with d = 1 cm, and a single target speech source impinging on the microphone array at an angle of 45◦ (note that the microphone array has a vertical configuration). There is no noise source. Set the sampling frequency of the RIRs to 44.1 kHz, and set the reverberation time to 0 s (no reverberation). Try to predict what the RIRs for both microphones will look like, and in what sense the two microphone signals will differ from each other. Confirm by simulation.

