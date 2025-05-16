# ChangeDetectionPANN

This repository contains the offical code used in the Training-free AI for Earth Observation Change Detection using Physics Aware Neuromorphic Networks as publised in ... 

## Training-free AI for Earth Observation Change Detection using Physics Aware Neuromorphic Networks

> **Abstract:** Earth observations from low Earth orbit satellites provide vital information for decision makers to better manage time-sensitive events such as natural disasters. For the data to be most effective for first responders, low latency is required between data capture and its arrival to decision makers. A major bottleneck is in the bandwidth-limited downlinking of the data from satellites to ground stations. One approach to overcome this challenge is to process at least some of the data on-board and prioritise pertinent data to be downlinked. In this work we propose a Physics Aware Neuromorphic Network (PANN) to detect changes caused by natural disasters from a sequence of multi-spectral satellite images and produce a change map, enabling relevant data to be prioritised for downlinking.
>
> The PANN used in this study is motivated by physical neural networks comprised of nano-electronic circuit elements known as ``memristors'' (nonlinear resistors with memory).  The weights in the network are dynamic and update in response to varying input signals according to memristor equations of state and electrical circuit conservation laws. The PANN thus generates physics-constrained dynamical output features which are used to detect changes in a natural disaster detection task by applying a distance-based metric. Importantly, this makes the whole model training free, allowing it to be implemented with minimal computing resources. The PANN was benchmarked against a state-of-the-art AI model and achieved comparable to or better results in each natural disaster category. It thus presents a promising solution to the challenge of resource-constrained on-board processing.

### Example Usuage

Example scripts for running the model, evaluating the output and change maps and creating the feature space plots are provided in scripts. File paths will need to be changed in the example scripts. When evaluating the change maps or creating the feature space plots it is expected that the relevant files produced when running the model are in the same directory as the script. If you want to save the outputs from the PANN model to create the feature space visualisation plots, make sure `saveFiles = true` when running the model.

### Data

The evaluation data used in this project is from the [RaVAEn paper](https://github.com/spaceml-org/RaVAEn/tree/master). It is expected to be in the same format and directory structure as the original data.
