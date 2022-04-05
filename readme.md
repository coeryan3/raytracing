# Assignment 1a: Getting Started with Ray Casting
#### Author: Ryan Coe (coe00003)

## Overall Description
###### *Similar to hw1a description. Added description of Blinn-Phong equation in shade-ray.*

My file, raycasting.ppm, takes in a file input from the user where
the file contains input information about the scene in which the program
is supposed to program. The program will search the file for specific
identifiers and assign the values associated with them to specified 
local variables. The program then re-loops through the data again for the 
data regarding spheres and material colors. 

Once this is complete, the program proceeds to check the data collected from
the input file and determines whether the input values satisfy the 
basic requirements.

Thus, moving forward, the program calculates the viewing window coordinates: 
u, v, w. Once calculated the four corners of the viewing window are then 
determined for further use within the program.

The values found from the four corners are then used to calculate step values
to map pixels to the viewing window. By mapping out the pixels to the viewing 
window, the program can then calculate the rays from points in the viewing 
window back to the viewing position (eye).

Then, the program calls the funciton ray_trace, which uses the parameters of 
a ray, the spheres within the scene, the viewing position, and the different 
colors. This then uses the ray formula and the sphere formula to determine
whether there is an intersection between the ray and a sphere. If there is 
an intersection the pixel which the ray is pointing to is then the function 
shade_ray is called. This function utilizes the Blinn-Phong Illumination
equation to calculate the color based on the input values provided for 
material colors (mtlcolors) and light sources. The function also calls a
helper function to help calculate the shadow ray if any existing objects 
are between the intersection point and the light source. After calculating
the values from the illumnation equation, a color is then assigned to a pixel
in the array img. 

To finish things off the program then takes the data of the pixels (img) and generates 
a ppm file. This file is then available for the user to open from the directory.


## How to Run
In order to run the program, type into command prompt:

make
./raycasting

Or type:

make test

which will then prompt the user to enter in a file name. 
For example, an input file my be ray_input1b.txt.


## Error with Code
Overall the program seems to function and return a ppm image. 
However, there was one downside to the current program which will need
to be fixed for further work. When the input file contains multiple 
spheres with different colors, it has been noticed that the output image 
displays each sphere as the same/similar base color. I might also 
hypothesize that some value/variable is not reseting when calling shade-ray.
This can be seen in the image provided multiple_spheres.ppm.

Due to timing reasons, the analysis was conducted only looking at one 
sphere rather than multiple spheres of different colors.
