# Local-Feature-Matching

Based on the assignment developed by James Hays, a local feature matching algorithm is developed; starting with extracting interest points, to preparing local feature descriptors and finally matching the points between different images. For the step of extracting interest points, Harris corner detection algorithm was implemented. A SIFT-like local feature descriptor was chosen for the local feature description phase, and a nearest neighbor ratio test was executed in the matching step. The algorithm was tested on three pairs of images: Notre Dame de Paris, Mount Rushmore, and the Gaudi's Episcopal Palace, and the results were evaluated against ground truth.

# Running the Code

In order to test the best chosen combination of the different algorithms and parameters used, kindly call the function "proj2_averageAccuracy( )". It will display the precision and accuracy on the first 100 matches of each of the 3 image pairs. To display the precision and accuracy on the best 100 matches of each of the 3 image pairs when considering all the matches, change "maxPtsEval" variable on line 7 of the function to be 1000.

In order to display the matches as dots colored with green if they are good/correct and red if they are bad/wrong, call the function "proj2()" and:
 - set the "imagePair" variable on line 14 to be the number of the pair you want to test
 - set the "num_pts_to_evaluate" variable on line 68 to be as desired

For the detection step, the Harris Corner Detection algorithm is the one running by default. If you wish to change it to the SIFT detector, comment the lines titled "Harris Corner Detector" and uncomment those titled "SIFT detector" in the function "get_interest_points()".

For the description step, the simplified version of the SIFT descriptor is the one running by default. If you wish to change it to the relatively complicated version, comment the lines titled "Relatively Complicated Descriptor" and uncomment those titled "Simplified Descriptor" in the function "get_features()".

For the matching step, the Euclidean distance is manually calculated by default. If you wish to test the KNNSearch algorithm, comment the lines titled "Manually Computing Euclidean Distance" and uncomment those titled "KNNSearch Algorithm" in the function "match_features()".
