# Generating-SEM-image-for-segmentation
A classical interactive method for homogeneous nano-paritcle cluster image generation.

Due to a terrific effort to obtain labels for some certain types of medical images, it is understandable to use some codes to reduce workload. All codes are based on Matlab R2018a and are still in progress. For a more detailed discussion, please feel free to open an issue.


## The whole flowchart

![Flowchart](https://github.com/AdamGreen95/Generating-SEM-image-for-segmentation/raw/master/20200904200408.png)

To summarise the whole project briefly, we drew a flowchart as an auxiliary demostration. 
To begin with, we get a image featuring with a lot of particle closely clustering. we adopt a build-in function: 
``` 
[impoly](https://de.mathworks.com/help/images/ref/impoly.html)
``` 
in matlab to segment a target particle as a seed which will be used for later iteration process. Then we set multiple sets of variables to determine the parameters which will effect how the particles are placed. We do this step because too many paremeters are inconvenient while not all parameters can significantly affect the characteristics of the image. The number of particles determines the number of iterations, for every iteration, we adopt random or fixed changes to current particle, such as size change, rotation, shading and so on. It must be mentioned that the position where we place the particle is followed 3 different strategies in our codes. If we treat the core part of flowchat as a module, the input of the module is a homogeneous cluster SEM image with manual interactions, a predifine the variables as well, the outputs are the synthesised iamge and its label. Obviously different sets of variables will gain different images. 
