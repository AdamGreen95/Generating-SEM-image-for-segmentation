# Generating-SEM-image-for-segmentation
A classical interactive method for homogeneous nano-paritcle cluster image generation.

Due to a terrific effort to obtain labels for some certain types of medical images, it is understandable to use some codes to reduce workload. All codes are based on Matlab R2018a and are still in progress. For a more detailed discussion, please feel free to open an issue.

## The whole flowchart

![Flowchart](https://github.com/AdamGreen95/Generating-SEM-image-for-segmentation/raw/master/20200904200408.png)

I'll summarise the whole progress briefly. To begin with, we get a image with a lot of particle clustering closely, we adopt build_in function: [impoly](https://de.mathworks.com/help/images/ref/impoly.html) in matlab to segment a target particle as an seed which will be used for later iteration process. 
If we treat the core part of flowchat as a module, the input of the module is a homogeneous cluster SEM image, while outputs are the synthesised iamge and its label. Obviously different sets of parameter will gain different images. 
