# EELS-simulator general overview
In this project, we model high-brightness beams within an electron microscope as 2D normal distributions within phase space and find the optimum configuration of optics for achieving maximum time and space resolution to be one that nearly maximizes the transverse slope of the distribution before entering the camera, and minimizes the longitudinal slope of the distribution before entering the specimen. 

## Introduction
Traditionally, electron microscopes have used long sequences of small packets of electrons to scan their specimen in order to minimize the expansion of the packets as the electrons repulse each other. This however has the drawback of having limited time-resolution, as in order to fully scan a material, the entire sequence of packets must be sent.

For observing systems that degrade under electron pulses (photoactive or nanoelectronic devices) or rapidly change (chemical reactions), all the information must be collected in as few pulses as possible in order to obtain accurate data.

The solution to this problem is high-brightness beams with adaptive optics. With more electrons, more data is collected, where the optics serve to either compress the pulse in time before hitting the specimen or to compress the pulse in space before hitting the camera in order to achieve greater time and space resolution respectively.

This simulator thus provides both a model of the concept of such a system as well as a way to determine the optimal configuration of the pre-camera phase space.

## Modeling the system

We model the electron pulse as a set of two 2D normal distributions in phase space where phase space is simply a coordinate system where the y axis represents velocity and the x axis represents position in both the transverse and longitudinal directions.

As the pulse travels through the microscope it undergoes several changes along the way. The first is free expansion- since the faster electrons will travel farther away from the slower electrons over time, the distribution either stretches or compresses over time depending on whether the faster electrons are in front of or behind the slower ones. We take advantage of this fact to orient the distribution to our liking.

This is shown below in fig a and fig c.

![alt text](https://github.com/ZovcIfzm/EELS-simulator/blob/master/images/manipulations123.png)

The second is lensing. As the distribution travels through the RF lenses, the field within the RF cavity is tuned so that electrons farther back receive a greater increase in energy than electrons further in front. The result of this is that the slope of the distribution is flipped, so that if previously the faster electrons were in the front of the distribution, they are now in the back, as shown in fig b above. As mentioned before, we can then take advantage of this to orient the distribution vertically if we want to as shown in fig c. The magnetic lenses work in the same way, just in the transverse direction rather than the RF's longitudinal manipulation.

The third is specimen energy loss. As the distribution travels through the specimen, different percentages of it lose different amounts of energy based on the composition of the material. This data is later recovered by the camera to allow us to determine exactly what that composition is. This shattering is modeled by creating hundreds of copies of the initial distribution, each with their own distinct velocities as a result of the energy loss.

This is shown in the figure below.

![alt text](https://github.com/ZovcIfzm/EELS-simulator/blob/master/images/shattering.png)

The last is an analyzer shift. As the distributions travel through the analyzer, which contains a magnetic field, they experience a Lorentz force in the transverse direction directly proportional to their velocity. This separates the distributions that experienced less energy loss from those that experienced more energy loss, allowing us to isolate the frequency of the different amounts of energy loss, allowing us to determine the composition of the material, as different elements of the composition cause different amounts of energy loss.


## Finding the optimum phase space configuration
While these are the basic concepts simulating the evolution of an electron pulse, now we must use them to solve a problem: What is the optimum configuration of the pre-specimen/pre-camera phase space for optimal time and space resolution? To do this we must first establish certain statistical metrics and algorithms.

## Deviation
The metric we will use to compare how closely our simulated spectrum is from the actual spectrum is standard deviation. While not a perfect statistical measurement for this situation, it will do as it increases in magnitude as the two spectrums deviate. This is calculated by comparing the populations at each corresponding point within each spectrum (the beginning lines, the second lines, the third lines, and so on). For example, the spectrum for the particular material we will be scanning, hexagonal Boron Nitride powder, is shown below and we would be comparing each of those distinct lines with that given by our simulation.

![alt text](https://github.com/ZovcIfzm/EELS-simulator/blob/master/images/spectrum.png)

## Pixel integration
In order to see whether or not our configuration yields greater resolution or not, we need to model the camera itself. The camera for our model is a set of pixels in the transverse directions of a certain width that register all the energy within that transverse width. For example if we had 10 pixels of equal length, we would split them among the length of all the distributions, and they would receive all the energy within that particular column of the microscope.

To calculate all the energy within that column, we must analytically integrate the portions of each distribution within that column. However since normal distributions extend infinitely, this would mean we'd have to integrate each distribution. Since we have limited computing power, we choose to approximate within 99.9999% accuracy by integrating only distributions whose center is within 5 standard deviations of either edge of a pixel.

This itself is slightly overkill as 3 standard deviations already holds 99.7% of all data, but since there is little difference in computational time for 3 deviations vs 5 as opposed to none at all, we will indulge ourselves.

## Optimization procedure
The optimization is split into two portions. Before hitting the specimen, we want to pulse to be as compact along the longitudinal lens as possible so that the time between when the pulse hits the specimen and when the pulse leaves the specimen is as small as possible, this is the time resolution. This is trivially verifiable, if the pulse is compact along the longitudinal lens we know for sure that the time resolution will be better as the distance between the front and end of the pulse is shorter.

However finding the optimum configuration of optics for space compression is not as clear cut. Does compacting the pulse along the transverse lens create more accurate space resolution as there is less overlap between distributions, allowing us to more clearly see the particular populations for each amount of energy loss? As opposed to time resolution, this is not trivially verifiable as instead of simply distance being the only dependent condition for better resolution, now we have to deal with how the distributions interact with each other, the specimen, the analyzer, and the camera, and then deal with any aberrations they may cause.

We begin with a simple 512 pixel setup, with specimen data provided for 1024 resultant mini-distributions as a result of shattering through the specimen. With this setup, we run simulations continuously for progressively smaller and smaller pulse widths as a result of free contraction, as well as eventually larger pulse widths as that free contraction becomes free expansion as the faster electrons cross over the slower electrons.

![alt text](https://github.com/ZovcIfzm/EELS-simulator/blob/master/images/SlopeVsDeviation50.png)
(fig 1)

![alt text](https://github.com/ZovcIfzm/EELS-simulator/blob/master/images/widthDeviation.png)
(fig 2)

As seen above in fig 2, as the width decreases in size where width is the transverse standard deviation of the distribution, the deviation between the simulated spectrum and the actual spectrum decreases. This is further and perhaps more clearly illustrated with fig 1, the slope graph, where slope is inversely proportional to width.

However, if the pixel count and shatter number are not cleanly related by an integer multiple, while deviation still tends to decrease as width decreases, after a certain point aberrations appear within the spectrum and cause increasing deviation, opposite to what we would expect.

An example is given below with a 500 pixel system. Where fig 4 is the graph of width vs. deviation, and fig 5 is the spectrum of a focused distribution (one with minimal width)

![alt text](https://github.com/ZovcIfzm/EELS-simulator/blob/master/images/abberatedDeviation.png)
(fig 4)
![alt text](https://github.com/ZovcIfzm/EELS-simulator/blob/master/images/abberatedSpecFocused.png)
(fig 5)

If we increase the pixel count (to say 5000) the overall graph has a reduced overall abberated shape, however the data between adjacent pixels is extremely variant, causing an overall unclean spectrum as seen below in fig 6. Similarly as the pixel count is now higher than the number of different frequencies given to us to compare to (1024), we can no longer measure deviation.

![alt text](https://github.com/ZovcIfzm/EELS-simulator/blob/master/images/5000pixAbb.png)
(fig 6)

## Conclusion
While increasing time resolution simply requires maximal time compression, maximizing space resolution is a little more tricky. As there are a discrete number of pixels on the camera taking information from a continous discribution, there is naturally an increasing amount of data loss with increasing space compression. As a result, it is safe to say that to maximize space resolution requires maximizing space compression up until data loss becomes detrimental.

## Possible sources of error
As a result of this being a simulation, there is an inherent source of error in not accounting for real-life effects that have been excluded. For example, uncertainty is not accounted for, nor are the effect of changes to the specimen from being hit by a high-intensity electron pulse that causes data collected from the front of the pulse to be different than that of the back.

Given the intense statistical nature of this simulator as well, errors may crop up in the implementation of statistical methods. For example, an edge case may not be accounted for in which some parts of a distribution may be double counted when summing, or perhaps an equation somewhere was implemented wrong.

Similarly given the sheer number of calculations and the precision required, the fact that almost all of the data is represented in a double data type format (a floating point that represents a value as a power of 2) resulting floating point error may be stacking up to create deviations within the data that may not be equal among all types of distributions.

Finally the fact that we have used some approximations must be accounted for. In the analytical integration, we have not calculated the intensity of every point within the distribution, rather we have split the distribution into a rectangular grid, calculating the value at the middle of each rectangle, having it represent the average value of that rectangle. Similarly we have not integrated every distribution within our pixel sum calculations- leaving out 0.0001% of data per distribution which could add up given that there are 1024 total shattered distributions in total.

## Future work
Future work for this project will consist of improving the accuracy of the simulator as currently it has just been freshly completed, and just with any freshly completed large project it probably contains some bugs that cause inaccuracies.

Similarly the code will be made more concise due to the sheer amount of repetition present due to the spaghetti nature of the code. Time and space (memory) optimizations will also be done to the code as well, in addition to the simulation.

Potential future work could be in creating an algorithm that calculates the exact configuration needed of optics needed for specific real-life electron microscope systems, as well as creating a GUI that allows easier usage of the program for repeated use.

## Acknowledgments and references
The author acknowledges C.-Y Ruan for providing helpful references as well as assistance with conceptual understanding, as well as Arham Jain for invaluable help in constructing and optimizing the early simulator.

1 J. Williams, F. Zhou, T. Sun, Z. Tao, K. Chang, K. Makino, M. Berz, P. M. Duxbury, and C.-Y. Ruan, Structural Dynamics 4, 044035 (2017)

2 F. Zhou, J. Williams, and C.-Y. Ruan, Elsevier 0009-2614 (2017)

3 C. -Y. Ruan, Science 354 (6310), 283-284


