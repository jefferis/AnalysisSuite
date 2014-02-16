# runitAmiraFileFunctions.R

# Some test neurons to verify that all is well for the 
# AmiraFileFunctions rountines
# now going to use the RUnit package (see pdf in doc dir)

# runTestFile(file.path(TestDir,"IO","runitAmiraFileFunctions.R"))
require(RUnit)

AmiraTutorialDirectory=NULL
if(file.exists("/Applications")){
	# we're on a mac, find the latest version of amira
	AmiraRootDirectory=rev(list.files("/Applications",patt="Amira*",full=TRUE))[1]
	AmiraTutorialDirectory=file.path(AmiraRootDirectory,"/data/tutorials/")
	if(!file.exists(AmiraTutorialDirectory)) AmiraTutorialDirectory=NULL
}


test.ReadAmiramesh<-function(){
  neurites.data=ReadAmiramesh(file.path(TestDir,"Data","neurons","Neurites.am"))
  neurites.rds=readRDS(file.path(TestDir,"Data","neurons","Neurites.rds"))
  checkEquals(neurites.data,neurites.rds)
  
  neurites.header=ReadAmiramesh.Header(file.path(TestDir,"Data","neurons","Neurites.am"))
  neurites.hrds=readRDS(file.path(TestDir,"Data","neurons","Neurites.header.rds"))
  checkEquals(neurites.header,neurites.hrds)
}

test.ReadAM3D<-function(){
	# Test output of Felix Ever's skeletonize routines
	fieldsToCheck=c("NeuronName", "NumPoints", "StartPoint", "BranchPoints", "EndPoints", 
	"NumSegs", "SegList", "nTrees", "d")
	
	result=structure(list(NeuronName = "Neurites", NumPoints = 291L, StartPoint = 1L, 
	    BranchPoints = c(98L, 256L, 272L), EndPoints = c(1L, 54L, 
	    202L, 257L, 274L), NumSegs = 7L, SegList = seglist(c(1, 3, 4, 
	    5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 
	    21, 22, 2, 24, 25, 26, 27, 28, 29, 30, 31, 32, 33, 34, 35, 
	    36, 37, 38, 39, 40, 41, 42, 43, 44, 45, 46, 47, 48, 49, 50, 
	    51, 52, 53, 23, 55, 56, 57, 58, 59, 60, 61, 62, 63, 64, 65, 
	    66, 67, 68, 69, 70, 71, 72, 73, 74, 75, 76, 77, 78, 79, 80, 
	    81, 82, 83, 84, 85, 86, 87, 88, 89, 90, 91, 92, 93, 94, 95, 
	    96, 97, 98), c(98, 54), c(98, 100, 101, 102, 103, 104, 105, 
	    106, 107, 108, 109, 110, 111, 112, 113, 114, 115, 116, 117, 
	    118, 119, 120, 121, 122, 123, 124, 125, 126, 127, 128, 129, 
	    130, 131, 132, 133, 134, 135, 136, 137, 138, 139, 140, 141, 
	    142, 143, 144, 145, 146, 147, 99, 149, 150, 151, 152, 153, 
	    154, 155, 156, 157, 158, 159, 160, 161, 162, 163, 164, 165, 
	    166, 167, 168, 169, 170, 171, 172, 173, 174, 175, 176, 177, 
	    178, 179, 180, 181, 182, 183, 184, 185, 186, 187, 188, 189, 
	    190, 191, 192, 193, 194, 195, 196, 197, 198, 199, 200, 201, 
	    148, 203, 204, 205, 206, 207, 208, 209, 210, 211, 212, 213, 
	    214, 215, 216, 217, 218, 219, 220, 221, 222, 223, 224, 225, 
	    226, 227, 228, 229, 230, 231, 232, 233, 234, 235, 236, 237, 
	    238, 239, 240, 241, 242, 243, 244, 245, 246, 247, 248, 249, 
	    250, 251, 252, 253, 254, 255, 256), c(256, 202), c(256, 258, 
	    259, 260, 261, 262, 263, 264, 265, 266, 267, 268, 269, 270, 
	    271, 272), c(272, 273, 257), c(272, 275, 276, 277, 278, 279, 
	    280, 281, 282, 283, 284, 285, 286, 287, 288, 289, 290, 291, 
	    274)), nTrees = 1L, d = structure(list(PointNo = 1:291, Label = c(2, 
	    2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 
	    2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 
	    2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 
	    2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 
	    2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 
	    2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 
	    2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 
	    2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 
	    2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 
	    2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 
	    2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 
	    2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 
	    2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 
	    2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 
	    2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 
	    2, 2, 2, 2, 2), X = c(374.888, 367.069, 374.368, 373.722, 
	    373.54, 373.06, 372.566, 372.127, 371.78, 371.355, 371.035, 
	    370.838, 370.688, 370.279, 370.065, 369.657, 369.463, 369.148, 
	    368.773, 368.404, 367.915, 367.46, 353.928, 366.653, 366.212, 
	    365.852, 365.353, 364.904, 364.445, 363.922, 363.397, 362.923, 
	    362.431, 361.949, 361.538, 361.045, 360.716, 360.307, 359.857, 
	    359.463, 359.048, 358.636, 358.247, 357.909, 357.534, 357.09, 
	    356.569, 356.077, 355.622, 355.213, 354.876, 354.556, 354.252, 
	    333.369, 353.573, 353.202, 352.765, 352.335, 351.875, 351.421, 
	    350.954, 350.503, 350.076, 349.699, 349.23, 348.823, 348.355, 
	    347.894, 347.434, 346.964, 346.519, 346.088, 345.621, 345.191, 
	    344.733, 344.313, 343.89, 343.439, 342.975, 342.531, 342.051, 
	    341.585, 341.088, 340.637, 340.216, 339.763, 339.265, 338.77, 
	    338.271, 337.779, 337.29, 336.787, 336.276, 335.77, 335.291, 
	    334.826, 334.337, 333.855, 314.925, 333.431, 332.992, 332.523, 
	    332.064, 331.566, 331.112, 330.672, 330.289, 329.891, 329.427, 
	    328.943, 328.547, 328.139, 327.73, 327.263, 326.831, 326.391, 
	    325.979, 325.62, 325.355, 324.995, 324.667, 324.359, 324.053, 
	    323.735, 323.446, 323.087, 322.742, 322.442, 322.189, 321.906, 
	    321.57, 321.237, 320.855, 320.465, 320.063, 319.7, 319.292, 
	    318.879, 318.48, 318.082, 317.704, 317.319, 316.905, 316.475, 
	    316.045, 315.621, 315.262, 293.393, 314.521, 314.097, 313.667, 
	    313.274, 312.86, 312.426, 311.999, 311.585, 311.14, 310.716, 
	    310.291, 309.868, 309.473, 309.108, 308.732, 308.352, 307.986, 
	    307.601, 307.212, 306.839, 306.433, 306.074, 305.695, 305.344, 
	    305.029, 304.767, 304.485, 304.322, 303.999, 303.742, 303.387, 
	    303.024, 302.64, 302.239, 301.846, 301.347, 300.981, 300.541, 
	    300.122, 299.689, 299.236, 298.752, 298.364, 297.906, 297.534, 
	    297.064, 296.574, 296.167, 295.678, 295.294, 294.802, 294.327, 
	    293.862, 275.121, 292.979, 292.551, 292.082, 291.67, 291.281, 
	    290.922, 290.559, 290.221, 289.845, 289.478, 289.128, 288.804, 
	    288.504, 288.205, 287.898, 287.608, 287.347, 287.088, 286.884, 
	    286.679, 286.425, 286.247, 286.058, 285.855, 285.57, 285.272, 
	    284.972, 284.721, 284.576, 284.46, 284.221, 283.92, 283.517, 
	    283.187, 282.757, 282.293, 281.837, 281.371, 280.97, 280.587, 
	    280.223, 279.848, 279.444, 279.015, 278.557, 278.093, 277.641, 
	    277.164, 276.729, 276.318, 275.959, 275.645, 275.467, 275.273, 
	    268.316, 274.985, 274.684, 274.324, 273.95, 273.466, 273.011, 
	    272.56, 272.122, 271.703, 271.273, 270.948, 270.647, 270.339, 
	    269.921, 269.4, 268.847, 269.844, 269.25, 269.297, 269.457, 
	    269.693, 269.842, 269.914, 269.975, 269.998, 270.004, 270.004, 
	    270.016, 270.076, 270.148, 270.164, 270.134, 270.053, 269.952
	    ), Y = c(127.498, 134.529, 127.91, 128.142, 128.732, 129.072, 
	    129.291, 129.574, 129.973, 130.268, 130.671, 131.121, 131.576, 
	    131.884, 132.313, 132.596, 133.006, 133.344, 133.633, 133.918, 
	    134.038, 134.25, 141.004, 134.77, 134.832, 135.048, 135.129, 
	    135.284, 135.427, 135.327, 135.326, 135.565, 135.806, 136.055, 
	    136.402, 136.576, 136.97, 137.27, 137.475, 137.761, 137.995, 
	    138.223, 138.474, 138.76, 139.007, 139.161, 139.258, 139.454, 
	    139.667, 139.9, 140.188, 140.466, 140.748, 144.76, 141.07, 
	    141.079, 141.028, 141.074, 141.077, 141.207, 141.391, 141.627, 
	    141.889, 142.213, 142.355, 142.632, 142.791, 142.98, 143.181, 
	    143.362, 143.584, 143.825, 143.991, 144.222, 144.314, 144.336, 
	    144.445, 144.598, 144.775, 144.973, 145.076, 145.205, 145.151, 
	    145.128, 145.099, 145.166, 145.118, 145.114, 145.121, 145.163, 
	    145.167, 145.122, 145.106, 145.117, 144.984, 144.832, 144.796, 
	    144.78, 132.196, 144.586, 144.405, 144.287, 144.112, 144.067, 
	    143.929, 143.814, 143.593, 143.394, 143.345, 143.285, 143.041, 
	    142.826, 142.584, 142.413, 142.206, 142.019, 141.737, 141.398, 
	    140.983, 140.651, 140.291, 139.924, 139.559, 139.195, 138.793, 
	    138.459, 138.112, 137.718, 137.289, 136.889, 136.541, 136.192, 
	    135.892, 135.592, 135.299, 134.955, 134.669, 134.383, 134.109, 
	    133.858, 133.599, 133.35, 133.113, 132.926, 132.761, 132.623, 
	    132.417, 118.666, 132.05, 131.875, 131.696, 131.497, 131.345, 
	    131.181, 131.079, 130.955, 130.861, 130.727, 130.591, 130.372, 
	    130.103, 129.799, 129.519, 129.29, 129.051, 128.823, 128.543, 
	    128.244, 128.004, 127.744, 127.441, 127.107, 126.75, 126.371, 
	    126.005, 125.563, 125.197, 124.767, 124.393, 124.028, 123.677, 
	    123.327, 122.961, 122.78, 122.412, 122.141, 121.865, 121.648, 
	    121.473, 121.329, 121.03, 120.839, 120.505, 120.305, 120.158, 
	    119.851, 119.705, 119.382, 119.272, 119.094, 118.883, 104.161, 
	    118.403, 118.16, 118.009, 117.756, 117.471, 117.145, 116.817, 
	    116.469, 116.167, 115.855, 115.531, 115.195, 114.832, 114.473, 
	    114.196, 113.976, 113.711, 113.451, 113.188, 112.888, 112.588, 
	    112.217, 111.876, 111.55, 111.339, 111.122, 110.939, 110.842, 
	    110.675, 110.421, 110.144, 109.787, 109.507, 109.183, 109.025, 
	    108.881, 108.722, 108.559, 108.269, 107.953, 107.623, 107.325, 
	    107.117, 107.037, 106.931, 106.762, 106.549, 106.4, 106.154, 
	    105.875, 105.555, 105.238, 104.867, 104.524, 102.149, 104.142, 
	    103.75, 103.42, 103.083, 102.868, 102.637, 102.515, 102.415, 
	    102.263, 102.174, 102.22, 102.307, 102.279, 102.152, 102.093, 
	    102.156, 107.951, 102.476, 102.893, 103.305, 103.706, 104.127, 
	    104.53, 104.874, 105.156, 105.423, 105.705, 105.992, 106.278, 
	    106.555, 106.817, 107.076, 107.351, 107.644), Z = c(61.998, 
	    64.303, 62.069, 62.07, 62.284, 62.274, 62.147, 62.031, 62.065, 
	    62.097, 62.191, 62.366, 62.553, 62.613, 62.823, 62.974, 63.223, 
	    63.461, 63.668, 63.891, 64.009, 64.134, 66.716, 64.483, 64.731, 
	    65.029, 65.171, 65.387, 65.609, 65.631, 65.556, 65.504, 65.563, 
	    65.589, 65.638, 65.562, 65.587, 65.555, 65.493, 65.525, 65.57, 
	    65.667, 65.817, 66.041, 66.289, 66.543, 66.67, 66.729, 66.734, 
	    66.714, 66.715, 66.71, 66.72, 65.032, 67.055, 67.383, 67.611, 
	    67.856, 68.06, 68.244, 68.329, 68.308, 68.247, 68.2, 68.115, 
	    68.059, 68.01, 67.972, 68.015, 67.979, 67.895, 67.792, 67.724, 
	    67.613, 67.436, 67.166, 66.927, 66.787, 66.759, 66.641, 66.519, 
	    66.353, 66.251, 66.024, 65.73, 65.497, 65.389, 65.269, 65.182, 
	    65.067, 64.923, 64.843, 64.831, 64.829, 64.898, 64.969, 65.002, 
	    65.038, 58.312, 64.847, 64.68, 64.522, 64.393, 64.303, 64.122, 
	    63.901, 63.656, 63.416, 63.221, 63.095, 62.901, 62.704, 62.543, 
	    62.467, 62.313, 62.154, 62.073, 61.966, 61.845, 61.721, 61.589, 
	    61.437, 61.281, 61.153, 61.087, 60.977, 60.878, 60.814, 60.809, 
	    60.734, 60.631, 60.532, 60.451, 60.401, 60.371, 60.339, 60.281, 
	    60.232, 60.088, 59.909, 59.71, 59.524, 59.408, 59.272, 59.098, 
	    58.881, 58.604, 53.007, 58.066, 57.882, 57.721, 57.505, 57.292, 
	    57.14, 56.92, 56.681, 56.474, 56.235, 55.991, 55.803, 55.631, 
	    55.472, 55.322, 55.117, 54.89, 54.675, 54.543, 54.421, 54.275, 
	    54.058, 53.967, 53.864, 53.724, 53.527, 53.317, 53.11, 52.937, 
	    52.788, 52.691, 52.585, 52.491, 52.499, 52.536, 52.554, 52.562, 
	    52.604, 52.728, 52.881, 53.015, 53.033, 52.895, 52.773, 52.679, 
	    52.679, 52.644, 52.671, 52.655, 52.772, 52.883, 52.946, 52.985, 
	    44.643, 52.941, 52.88, 52.816, 52.706, 52.599, 52.505, 52.436, 
	    52.341, 52.238, 52.126, 51.995, 51.827, 51.662, 51.476, 51.189, 
	    50.844, 50.513, 50.18, 49.818, 49.492, 49.207, 48.95, 48.657, 
	    48.352, 48.009, 47.676, 47.317, 46.889, 46.43, 45.996, 45.624, 
	    45.383, 45.192, 44.947, 44.704, 44.549, 44.406, 44.321, 44.264, 
	    44.268, 44.218, 44.1, 43.907, 43.673, 43.509, 43.442, 43.418, 
	    43.389, 43.368, 43.415, 43.547, 43.766, 44.042, 44.342, 40.232, 
	    44.142, 43.977, 43.783, 43.618, 43.573, 43.43, 43.181, 42.899, 
	    42.608, 42.304, 41.875, 41.424, 40.954, 40.583, 40.343, 40.313, 
	    33.754, 40.05, 39.77, 39.528, 39.328, 39.081, 38.771, 38.395, 
	    37.969, 37.536, 37.113, 36.694, 36.279, 35.859, 35.424, 34.985, 
	    34.56, 34.152), W = c(1.388, 2.782, 1.562, 1.582, 1.804, 
	    1.994, 2.225, 2.517, 2.665, 2.612, 2.498, 2.337, 2.149, 2.115, 
	    2.096, 2.165, 2.186, 2.272, 2.405, 2.483, 2.564, 2.688, 4.078, 
	    2.772, 2.773, 2.756, 2.667, 2.443, 2.33, 2.182, 1.962, 1.948, 
	    1.756, 1.696, 1.756, 1.681, 1.718, 1.692, 1.626, 1.684, 1.798, 
	    1.992, 2.256, 2.664, 3.197, 3.621, 3.918, 4.21, 4.428, 4.558, 
	    4.546, 4.436, 4.239, 4.527, 3.84, 3.883, 4.072, 4.135, 4.137, 
	    4.001, 3.89, 3.81, 3.739, 3.701, 3.625, 3.569, 3.487, 3.468, 
	    3.505, 3.666, 3.919, 4.13, 4.228, 4.114, 4.072, 4.114, 4.201, 
	    4.28, 4.271, 4.184, 4.003, 3.854, 3.789, 3.754, 3.8, 3.896, 
	    4.303, 4.598, 4.716, 4.534, 4.41, 4.331, 4.336, 4.358, 4.395, 
	    4.455, 4.51, 4.52, 5.576, 4.37, 4.16, 4.027, 3.867, 3.817, 
	    3.759, 3.697, 3.683, 3.695, 3.66, 3.622, 3.634, 3.727, 3.868, 
	    4.121, 4.426, 4.606, 4.765, 4.664, 4.529, 4.289, 4.23, 4.249, 
	    4.33, 4.437, 4.81, 4.984, 5.127, 5.274, 5.5, 5.607, 5.615, 
	    5.627, 5.574, 5.557, 5.554, 5.608, 5.7, 5.884, 5.973, 5.966, 
	    5.97, 5.966, 5.918, 5.83, 5.763, 5.596, 5.571, 3.689, 5.512, 
	    5.423, 5.337, 5.5, 5.584, 5.72, 5.66, 5.818, 5.985, 6.223, 
	    6.43, 6.657, 6.647, 6.397, 6.018, 5.635, 5.317, 5.102, 5.039, 
	    4.92, 4.77, 4.518, 4.464, 4.485, 4.512, 4.53, 4.533, 4.452, 
	    4.375, 4.322, 4.169, 4.081, 3.965, 3.996, 4.086, 4.123, 4.113, 
	    4.046, 3.999, 3.893, 3.74, 3.619, 3.326, 3.179, 3.178, 3.233, 
	    3.142, 3.164, 3.14, 3.267, 3.401, 3.538, 3.613, 2.862, 3.715, 
	    3.743, 3.772, 3.803, 3.837, 3.874, 3.913, 3.956, 4, 4.047, 
	    4.098, 4.152, 4.208, 4.263, 4.313, 4.35, 4.369, 4.364, 4.33, 
	    4.267, 4.179, 4.069, 3.946, 3.818, 3.695, 3.587, 3.502, 3.444, 
	    3.414, 3.408, 3.416, 3.428, 3.431, 3.417, 3.383, 3.328, 3.256, 
	    3.176, 3.093, 3.015, 2.949, 2.897, 2.86, 2.839, 2.83, 2.829, 
	    2.833, 2.839, 2.844, 2.848, 2.851, 2.855, 2.858, 2.861, 4.71, 
	    3.006, 3.144, 3.273, 3.392, 3.505, 3.618, 3.736, 3.863, 3.999, 
	    4.139, 4.277, 4.406, 4.517, 4.604, 4.664, 4.699, 3.977, 4.668, 
	    4.662, 4.639, 4.594, 4.528, 4.442, 4.344, 4.242, 4.143, 4.057, 
	    3.99, 3.948, 3.929, 3.929, 3.942, 3.958, 3.972), Parent = c(-1, 
	    22, 1, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 
	    18, 19, 20, 21, 53, 2, 24, 25, 26, 27, 28, 29, 30, 31, 32, 
	    33, 34, 35, 36, 37, 38, 39, 40, 41, 42, 43, 44, 45, 46, 47, 
	    48, 49, 50, 51, 52, 98, 23, 55, 56, 57, 58, 59, 60, 61, 62, 
	    63, 64, 65, 66, 67, 68, 69, 70, 71, 72, 73, 74, 75, 76, 77, 
	    78, 79, 80, 81, 82, 83, 84, 85, 86, 87, 88, 89, 90, 91, 92, 
	    93, 94, 95, 96, 97, 147, 98, 100, 101, 102, 103, 104, 105, 
	    106, 107, 108, 109, 110, 111, 112, 113, 114, 115, 116, 117, 
	    118, 119, 120, 121, 122, 123, 124, 125, 126, 127, 128, 129, 
	    130, 131, 132, 133, 134, 135, 136, 137, 138, 139, 140, 141, 
	    142, 143, 144, 145, 146, 201, 99, 149, 150, 151, 152, 153, 
	    154, 155, 156, 157, 158, 159, 160, 161, 162, 163, 164, 165, 
	    166, 167, 168, 169, 170, 171, 172, 173, 174, 175, 176, 177, 
	    178, 179, 180, 181, 182, 183, 184, 185, 186, 187, 188, 189, 
	    190, 191, 192, 193, 194, 195, 196, 197, 198, 199, 200, 256, 
	    148, 203, 204, 205, 206, 207, 208, 209, 210, 211, 212, 213, 
	    214, 215, 216, 217, 218, 219, 220, 221, 222, 223, 224, 225, 
	    226, 227, 228, 229, 230, 231, 232, 233, 234, 235, 236, 237, 
	    238, 239, 240, 241, 242, 243, 244, 245, 246, 247, 248, 249, 
	    250, 251, 252, 253, 254, 255, 273, 256, 258, 259, 260, 261, 
	    262, 263, 264, 265, 266, 267, 268, 269, 270, 271, 272, 291, 
	    272, 275, 276, 277, 278, 279, 280, 281, 282, 283, 284, 285, 
	    286, 287, 288, 289, 290)), .Names = c("PointNo", "Label", "X", "Y", "Z", 
	    "W", "Parent"), class = "data.frame", row.names = c(NA, 
	    -291L))), .Names = c("NeuronName", "NumPoints", 
	"StartPoint", "BranchPoints", "EndPoints", "NumSegs", "SegList", 
	"nTrees", "d"))
	
	result.new=read.neuron(file.path(TestDir,"Data","neurons","Neurites.am"))
	checkEquals(result[fieldsToCheck],result.new[fieldsToCheck],tol=1e-4)
	
	nwm=read.neuron(file.path(TestDir,"Data","neurons","neuron_with_materials.am"))
	labels=c(1L, 7L, 7L, 7L, 7L, 7L, 7L, 7L, 7L, 7L, 7L, 7L, 7L, 7L, 7L, 
	7L, 7L, 7L, 7L, 7L, 7L, 7L, 7L, 7L, 7L, 7L, 7L, 7L, 7L, 7L, 7L, 
	7L, 7L, 7L, 7L, 7L, 7L, 7L, 7L, 7L, 7L, 7L, 7L, 7L, 7L, 3L, 4L, 
	4L, 4L, 4L, 4L, 4L, 4L, 4L, 4L, 4L, 4L, 4L, 4L, 4L, 4L, 4L, 4L, 
	4L, 4L, 4L, 4L, 4L, 4L, 4L, 4L, 0L, 7L, 7L, 7L, 7L, 7L, 7L, 4L, 
	4L, 4L, 4L, 4L, 4L, 4L, 4L, 4L, 4L, 4L, 4L, 4L, 4L, 4L, 4L, 4L, 
	3L, 3L, 3L, 3L, 3L, 3L, 3L, 3L, 3L, 3L, 3L, 3L, 3L, 3L, 3L, 3L, 
	3L, 3L, 3L, 3L, 2L, 2L, 2L, 2L, 2L, 2L, 2L, 2L, 2L, 2L, 2L, 2L, 
	2L, 2L, 2L, 2L, 2L, 2L, 2L, 2L, 2L, 2L, 2L, 2L, 2L, 2L, 2L, 2L, 
	2L, 2L, 2L, 2L, 2L, 2L, 2L, 2L, 2L, 2L, 2L, 2L, 2L, 2L, 2L, 2L, 
	2L, 2L, 2L, 2L, 2L, 2L, 2L, 2L, 2L, 2L, 2L, 2L, 2L, 2L, 2L, 2L, 
	2L, 2L, 2L, 2L, 2L, 2L, 2L, 2L, 2L, 2L, 2L, 2L, 2L, 2L, 2L, 2L, 
	2L, 2L, 2L, 2L, 2L, 2L, 2L, 2L, 2L, 2L, 2L, 2L, 2L, 2L, 2L, 2L, 
	2L, 2L, 2L, 2L, 2L, 2L, 2L, 2L, 2L, 2L, 2L, 2L, 2L, 2L, 2L, 2L, 
	2L, 2L, 2L, 2L, 2L, 2L, 2L, 2L, 2L, 2L, 2L, 2L, 2L, 2L, 2L, 2L, 
	2L, 2L, 2L, 2L, 2L, 2L, 2L, 2L, 2L, 2L, 2L, 2L, 2L, 2L, 2L, 2L, 
	2L, 2L, 2L, 2L, 2L, 2L, 2L, 2L, 2L, 2L, 2L, 2L, 2L, 2L, 2L, 2L, 
	2L, 2L, 2L, 2L, 2L, 2L, 2L, 2L, 2L, 2L, 2L, 2L, 2L, 2L, 2L, 2L, 
	2L, 2L, 2L, 2L, 2L, 2L, 2L, 2L, 2L, 2L, 2L, 2L, 2L, 2L, 2L, 2L, 
	2L, 2L, 2L, 2L, 2L, 2L, 2L, 2L, 2L, 2L, 2L, 2L, 2L, 2L, 2L, 3L, 
	3L, 3L, 3L, 3L, 3L, 3L, 3L, 3L, 3L, 3L, 3L, 3L, 3L, 3L, 3L, 3L, 
	3L, 3L, 3L, 3L, 3L, 3L, 3L, 3L, 3L, 3L, 3L, 3L, 3L, 3L, 3L, 3L, 
	3L, 3L, 3L, 3L, 3L, 3L, 3L, 3L, 3L, 3L, 3L, 3L, 3L, 3L, 3L, 3L, 
	3L, 3L, 3L, 3L, 3L, 3L, 3L, 3L, 3L, 3L, 3L, 3L, 3L, 3L, 3L, 3L, 
	3L, 3L, 3L, 3L, 3L, 3L, 3L, 3L, 3L, 3L, 3L, 3L, 3L, 3L, 3L, 3L, 
	3L, 3L, 3L, 3L)
	checkEquals(nwm$d$Label,labels)

	nwis = read.neuron(file.path(TestDir, "Data", "neurons", 
		"NeuritesWithIsolatedPoints_veryshort.am"))
		
	nwis_base=structure(list(NeuronName = "NeuritesWithIsolatedSegment", NumPoints = 14L, 
		    StartPoint = 1, BranchPoints = integer(0), EndPoints = c(1L, 
		    12L), NumSegs = 1L, SegList = list(c(1, 2, 3, 4, 5, 6, 7, 
		    8, 9, 10, 11, 12)), nTrees = 2L, d = structure(list(PointNo = 1:14, 
		        Label = c(2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2), 
		        X = c(374.926, 374.838, 374.779, 374.399, 373.922, 373.518, 
		        373.138, 372.689, 372.261, 371.859, 371.634, 371.257, 
		        373.549, 373.913), Y = c(127.815, 128.305, 128.772, 128.496, 
		        128.564, 128.755, 129.046, 129.264, 129.485, 129.77, 
		        130.204, 130.525, 128.38, 128.055), Z = c(62, 62.077, 
		        62.149, 62.068, 62.139, 62.218, 62.255, 62.174, 62.073, 
		        62.082, 62.137, 62.174, 63, 63), W = c(0.2, 0.411, 0.655, 
		        1.247, 1.564, 1.794, 1.988, 2.211, 2.461, 2.708, 2.717, 
		        2.678, 0.2, 0.2), Parent = c(-1, 1, 2, 3, 4, 5, 6, 7, 
		        8, 9, 10, 11, -1, 13), NeighbourCount = c(1L, 2L, 2L, 
		        2L, 2L, 2L, 2L, 2L, 2L, 2L, 2L, 1L, 1L, 1L), SubTree = c(1, 
		        1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 2, 2)), .Names = c("PointNo", 
		    "Label", "X", "Y", "Z", "W", "Parent", "NeighbourCount", 
		    "SubTree"), class = "data.frame", row.names = c(NA, -14L)), 
		    OrientInfo = structure(list(AxonOriented = TRUE, AVMPoint = NA, 
		        AVLPoint = NA, Scl = c(1, 1, 1), NewAxes = c(1, 2, 3)), .Names = c("AxonOriented", 
		    "AVMPoint", "AVLPoint", "Scl", "NewAxes"))), .Names = c("NeuronName", 
		"NumPoints", "StartPoint", "BranchPoints", "EndPoints", "NumSegs", 
		"SegList", "nTrees", "d", "OrientInfo"))
	
	checkEquals(nwis,nwis_base,"Compare parsed version of NeuritesWithIsolatedSegment")

	# FIXME - this still errors out because trees with only 1 point are
	# not read correctly 
	# nwip = ReadNeuronFromAM3D(file.path(TestDir, "Data", "neurons", 
	# 	"NeuritesWithIsolatedPoints_veryshort.am"))
}

test.ReadWriteNeuronFromAM<-function(){
  # Test reading of Amira Lineset format neurons using sample data from
  # Felix Evers' hxskeletonize plugin for Amira and
  # flycircuit.tw data
  fieldsToCheck=c("NumPoints", "StartPoint", "BranchPoints", "EndPoints", 
        "NumSegs", "SegList")
  fieldsToCheckLong=c("NeuronName", "NumPoints", "StartPoint", "BranchPoints", "EndPoints", 
    "NumSegs", "SegList", "d")
      
	tmpfile=tempfile()
	on.exit(unlink(tmpfile))
	am3d=read.neuron(file.path(TestDir,"Data","neurons","testneuron_am3d.am"))
	# converted to lineset in amira by hxskeletonize
	lineset=ReadNeuronFromAM(file.path(TestDir,"Data","neurons","testneuron_lineset.am"))
  # check seglist describes equivalent graph - they end up in different order
  g1=Neuron2Graph(lineset)
  g2=Neuron2Graph(am3d)
  checkEquals(g1,g2,msg='seglist graphs do not match')
#  checkEquals(lineset[fieldsToCheck],am3d[fieldsToCheck],tol=1e-6,
#    msg="Same tracing saved as AM3D and AM by Amira produces different results")
  
  WriteNeuronToAM(am3d,tmpfile)
  amfromam3d=ReadNeuronFromAM(tmpfile)
  checkEquals(amfromam3d[fieldsToCheck],am3d[fieldsToCheck],tol=1e-6,checkNames = FALSE,
    msg="AM3D and AM file saved by R produce different results")
  checkEquals(amfromam3d$d[,1:6],am3d$d[,1:6],tol=1e-3)
      
  fcam=ReadNeuronFromAM(file.path(TestDir,"Data","neurons","testneuron_fclineset.am"))
  WriteNeuronToAM(fcam,tmpfile,Force=TRUE)
  fcamresaved=ReadNeuronFromAM(tmpfile)
  checkEquals(fcamresaved[fieldsToCheck],fcam[fieldsToCheck],tol=1e-6,
    msg="AM file from flycircuit and AM file resaved by R produce different results")
  checkEquals(fcamresaved$d[,1:5],fcam$d[,1:5],tol=1e-6)
    
}


test.ReadWrite3DDensityAmiraBinary<-function(){
	testData=array(rnorm(10^3),dim=rep(10,3))
	tmpfile=tempfile()
	
	Write3DDensityToAmiraLattice(tmpfile,testData,ftype="binary",dtype='double')
	testData.new=Read3DDensityFromAmiraLattice(tmpfile)
	checkEquals(testData,testData.new,msg="Failed round trip read/write Amira binary data")
	unlink(tmpfile)

	Write3DDensityToAmiraLattice(tmpfile,testData,ftype="hxzip",dtype='double')
	testData.new=Read3DDensityFromAmiraLattice(tmpfile)
	checkEquals(testData,testData.new,msg="Failed round trip read/write Amira hxzip data")
	unlink(tmpfile)
	
	Write3DDensityToAmiraLattice(tmpfile,testData,ftype="binary",dtype='double',endian='big')
	testData.new=Read3DDensityFromAmiraLattice(tmpfile)
	checkEquals(testData,testData.new,msg="Failed round trip read/write Amira big endian data")
	unlink(tmpfile)

	Write3DDensityToAmiraLattice(tmpfile,testData,ftype="binary",dtype='double',endian='little')
	testData.new=Read3DDensityFromAmiraLattice(tmpfile)
	checkEquals(testData,testData.new,msg="Failed round trip read/write Amira big little data")
	unlink(tmpfile)
}

test.ReadHxZipData<-function(){
  amfile=file.path(TestDir,"IO","AL-a_M.am")
  x=Read3DDensityFromAmiraLattice(amfile)
  # count non-zero elements
  counts=c(6623L, 3304L, 2046L, 1529L, 1257L, 1054L, 907L, 706L, 657L, 
           557L, 458L, 444L, 399L, 325L, 307L, 247L, 269L, 224L, 185L, 196L, 
           186L, 147L, 150L, 146L, 120L, 138L, 122L, 105L, 94L, 95L, 90L, 
           86L, 88L, 81L, 88L, 66L, 63L, 78L, 53L, 76L, 53L, 43L, 47L, 49L, 
           54L, 46L, 42L, 43L, 33L, 38L, 37L, 28L, 26L, 19L, 29L, 33L, 21L, 
           30L, 21L, 25L, 12L, 19L, 20L, 19L, 8L, 14L, 14L, 15L, 10L, 16L, 
           14L, 15L, 9L, 9L, 13L, 10L, 8L, 9L, 12L, 6L, 6L, 8L, 7L, 9L, 
           9L, 7L, 7L, 3L, 6L, 9L, 9L, 5L, 7L, 6L, 6L, 3L, 6L, 5L, 7L, 10L, 
           8L, 4L, 5L, 2L, 3L, 4L, 4L, 6L, 1L, 6L, 4L, 1L, 6L, 0L, 6L, 10L, 
           0L, 6L, 2L, 4L, 3L, 4L, 4L, 4L, 3L, 6L, 1L, 3L, 2L, 3L, 9L, 2L, 
           1L, 3L, 3L, 1L, 2L, 4L, 2L, 3L, 6L, 2L, 4L, 2L, 2L, 4L, 0L, 1L, 
           4L, 3L, 1L, 3L, 2L, 1L, 5L, 2L, 3L, 1L, 4L, 3L, 1L, 3L, 0L, 0L, 
           1L, 0L, 2L, 1L, 2L, 3L, 2L, 3L, 2L, 1L, 2L, 4L, 1L, 0L, 0L, 2L, 
           1L, 5L, 2L, 7L, 1L, 2L, 5L, 3L, 2L, 0L, 0L, 1L, 1L, 1L, 1L, 3L, 
           0L, 1L, 0L, 1L, 1L, 0L, 1L, 2L, 0L, 3L, 2L, 2L, 0L, 3L, 1L, 1L, 
           2L, 1L, 2L, 0L, 2L, 2L, 1L, 0L, 1L, 0L, 1L, 0L, 0L, 2L, 0L, 0L, 
           2L, 3L, 0L, 0L, 0L, 0L, 2L, 0L, 3L, 0L, 1L, 0L, 2L, 0L, 0L, 2L, 
           2L, 2L, 1L, 2L, 0L, 0L, 0L, 0L, 0L, 0L, 1L)
  # nb 255 bins because 0 is ignored
  checkEquals(counts,tabulate(x,nbins=255))
}

test.ReadWrite3DDensityAmiraText<-function(){
	testData=array(rnorm(10^3),dim=rep(10,3))
	tmpfile=tempfile()
	
	Write3DDensityToAmiraLattice(tmpfile,testData,ftype="text")
	testData.new=Read3DDensityFromAmiraLattice(tmpfile)
	unlink(tmpfile)
	checkEquals(testData,testData.new,tol=1e-6)
}

test.ReadWriteAmiraLandmarksSingle<-function(){
	testData=matrix(rnorm(15),ncol=3)
	tmpfile=tempfile()
	WriteAmiraLandmarks(tmpfile,testData)
	testData.new=ReadAmiraLandmarks(tmpfile)
	unlink(tmpfile)
	names(testData.new)<-NULL
	checkEquals(testData,testData.new,tol=1e-6)	
}

test.ReadWriteAmiraLandmarksPaired<-function(){
	testData=replicate(2,matrix(rnorm(15),ncol=3),simplify=FALSE)
	tmpfile=tempfile()

	WriteAmiraLandmarks(tmpfile,testData)
	testData.new=ReadAmiraLandmarks(tmpfile)
	unlink(tmpfile)
	names(testData.new)<-NULL
	checkEquals(testData,testData.new,tol=1e-6)	
}
checkNonNullOutput<-function(expr,msg){
	t<-try(eval(expr))
	checkTrue(!inherits(t,'try-error') && !is.null(t),msg=msg)
	t
}

if(!is.null(AmiraTutorialDirectory)){
test.ReadAmiraTutorialLandmarks<-function(){
	landmarks<-checkNonNullOutput(
		ReadAmiraLandmarks(file.path(AmiraTutorialDirectory,"lobus.landmarks.am")),
		"Failed to read lobus.landmarks.am")
}
test.ReadAmiraTutorialVolumeData<-function(){
	lobus<-checkNonNullOutput(
		ReadAmiramesh(file.path(AmiraTutorialDirectory,"lobus.am")),
		"Failed to read lobus.landmarks.am")
}
test.ReadAmiraTutorialLabelsData<-function(){
	lobus.labels<-checkNonNullOutput(
		ReadAmiramesh(file.path(AmiraTutorialDirectory,"lobus.labels.am")),
		"Failed to read lobus.labels.am with ReadAmiramesh")
}
test.ReadAmiraTutorialLabelsDataWithRead3DDensityFromAmiraLattice<-function(){
	lobus.labels<-checkNonNullOutput(
		Read3DDensityFromAmiraLattice(file.path(AmiraTutorialDirectory,"lobus.labels.am")),
		"Failed to read lobus.labels.am with Read3DDensityFromAmiraLattice")
}
test.ReadAmiraTutorialSurfaceData<-function(){
	lobus.surf<-checkNonNullOutput(
		ParseAMSurfToContourList(file.path(AmiraTutorialDirectory,"lobus.surf")),
		"Failed to read lobus.surf")
	plot3dsurface(lobus.surf)
}
}

test.readGzipAmirameshNeuron<-function(){
  closeAllConnections()
  nopencons.start=length(showConnections())
  am3d=read.neuron(file.path(TestDir,"Data","neurons","testneuron_am3d.am"))
  checkException(read.neuron(file.path(TestDir,"IO","testneuron_am3d.am.gz")),
      silent=TRUE)
  am3da=read.neuron(file.path(TestDir,"Data","neurons","testneuron_am3d_ascii.am"))
  checkEquals(am3d,am3da,fieldsToExclude='NeuronName')
  am3daz=read.neuron(file.path(TestDir,"Data","neurons","testneuron_am3d_ascii.am.gz"))
  checkEquals(am3d,am3daz,fieldsToExclude='NeuronName')
  
  fcl=read.neuron(file.path(TestDir,"Data","neurons","testneuron_fclineset.am"))
  fclz=read.neuron(file.path(TestDir,"Data","neurons","testneuron_fclineset.am.gz"))
  checkEquals(fcl,fclz,fieldsToExclude='NeuronName')
  nopencons.end=length(showConnections())
  checkEquals(nopencons.end,nopencons.start,"Managed to leave a connection open")
}

test.ParseMaterials<-function(){
	amfile=file.path(TestDir,"Data","labels","ComplexLabelFieldHeader.am")
	h=ReadAmiramesh.Header(amfile,Verbose = FALSE,)
	materials_baseline=structure(list(id = c(1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 
							13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 
							29, 30, 31, 32, 33, 34, 35, 36, 37, 38, 39, 40, 41, 42, 43, 44, 
							45, 46, 47, 48, 49, 50, 51, 52, 53, 54, 55, 56, 57, 58, 59, 60, 
							61, 62, 63, 64, 65, 66, 67, 68, 69), col = structure(c("#B3CCCC", 
									"#1CFF0A", "#E10000", "#2929CC", "#13A8BF", "#87CC29", "#005F00", 
									"#4AC403", "#32EDFF", "#BD0000", "#A36500", "#48CC28", "#2B7BB0", 
									"#197E20", "#29CC67", "#87CC29", "#2974CC", "#CC4D29", "#FF0000", 
									"#FFCB00", "#FF8D00", "#661900", "#246614", "#3227CC", "#CC276B", 
									"#FFCB30", "#4300B2", "#47B279", "#0066FF", "#CD2895", "#E57B1A", 
									"#A425C3", "#CC8616", "#FFE102", "#5E28CC", "#28CC84", "#3B28CC", 
									"#CB1264", "#4AC403", "#31ECFF", "#BD0000", "#A36400", "#47CC27", 
									"#2A83B0", "#187D20", "#28CC66", "#87CC28", "#CC4C28", "#FF0000", 
									"#FFCB00", "#FF8D00", "#661900", "#236613", "#3226CC", "#CC266B", 
									"#FFCB2E", "#4200B1", "#47B178", "#0066FF", "#CD2793", "#E47B1A", 
									"#A425C3", "#CC8616", "#FFE102", "#5E27CC", "#3B27CC", "#CB1263", 
									"#28CC66", "#28CC66"), .Names = c("Exterior", "FB", "EB", "SAD", 
									"NO", "SOG", "PB", "CRE_R", "EPA_R", "VES_R", "ATL_R", "PLP_R", 
									"AVLP_R", "AL_R", "GOR_R", "SCL_R", "FLA", "ICL_R", "ME_R", "LOP_R", 
									"LO_R", "MB_R", "PVLP_R", "OTU_R", "WED_R", "SMP_R", "LH_R", 
									"SLP_R", "LB_R", "SIP_R", "IB_R", "IVLP_R", "IPS_R", "SPS_R", 
									"LAL_R", "PRW", "AME_R", "GA_R", "CRE_L", "EPA_L", "VES_L", "ATL_L", 
									"PLP_L", "AVLP_L", "AL_L", "GOR_L", "SCL_L", "ICL_L", "ME_L", 
									"LOP_L", "LO_L", "MB_L", "PVLP_L", "OTU_L", "WED_L", "SMP_L", 
									"LH_L", "SLP_L", "LB_L", "SIP_L", "IB_L", "IVLP_L", "IPS_L", 
									"SPS_L", "LAL_L", "AME_L", "GA_L", "PAN_L", "PAN_R"), class = "AsIs"), 
					level = c(0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 
							15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 
							30, 31, 32, 33, 34, 35, 36, 37, 38, 39, 40, 41, 42, 43, 44, 
							45, 46, 47, 48, 49, 50, 51, 52, 53, 54, 55, 56, 57, 58, 59, 
							60, 61, 62, 63, 64, 65, 66, 67, 68)), .Names = c("id", "col", 
					"level"), row.names = c("Exterior", "FB", "EB", "SAD", "NO", 
					"SOG", "PB", "CRE_R", "EPA_R", "VES_R", "ATL_R", "PLP_R", "AVLP_R", 
					"AL_R", "GOR_R", "SCL_R", "FLA", "ICL_R", "ME_R", "LOP_R", "LO_R", 
					"MB_R", "PVLP_R", "OTU_R", "WED_R", "SMP_R", "LH_R", "SLP_R", 
					"LB_R", "SIP_R", "IB_R", "IVLP_R", "IPS_R", "SPS_R", "LAL_R", 
					"PRW", "AME_R", "GA_R", "CRE_L", "EPA_L", "VES_L", "ATL_L", "PLP_L", 
					"AVLP_L", "AL_L", "GOR_L", "SCL_L", "ICL_L", "ME_L", "LOP_L", 
					"LO_L", "MB_L", "PVLP_L", "OTU_L", "WED_L", "SMP_L", "LH_L", 
					"SLP_L", "LB_L", "SIP_L", "IB_L", "IVLP_L", "IPS_L", "SPS_L", 
					"LAL_L", "AME_L", "GA_L", "PAN_L", "PAN_R"), class = "data.frame")
	checkEquals(h$Materials,materials_baseline)
}

test.is.amiramesh<-function(){
  amfile=file.path(TestDir,"Data","labels","ComplexLabelFieldHeader.am")
  checkTrue(is.amiramesh(amfile))
  # not enough to have a file ending
  tf=tempfile(fileext='.am')
  writeLines("#somethingelse",tf)
  on.exit(unlink(tf))
  checkTrue(!is.amiramesh(tf))
}

test.amiratype<-function(){
  amfile=file.path(TestDir,"Data","labels","ComplexLabelFieldHeader.am")
  h=ReadAmiramesh.Header(amfile,Verbose = FALSE)
  checkEquals(amiratype(amfile),'uniform.field')
}
  
neuron.isomorphic <- function(a, b) {
  ga=try(as.igraph(a),silent=T)
  gb=try(as.igraph(b),silent=T)
  isTRUE(try(graph.isomorphic(ga,gb),silent=T))
}

igraphvsoriginal.amiraneuron<-function(f,checkFn=checkTrue){
  ig=try(read.neuron.amiraskel(f,method='igraph',Verbose=FALSE))
    checkTrue(inherits(ig,'neuron'),
              msg='Failed to read neuron ",f," using igraph method')
    orig=try(read.neuron.amiraskel(f,method='original',Verbose=FALSE), silent=TRUE)
    if(inherits(orig,'try-error')) {
      message("Failled to read neuron ",f," using original method")
    } else {
      checkFn(neuron.isomorphic(ig, orig),
                msg=paste("hxskel neurons ",f,"are not equivalent by igraph and old methods"))
    }
}

find.am3d<-function(){  
  if(grepl('apple',R.version$platform)){
    ff=system('mdfind "AmiraMesh 3D"|grep "\\.am$"',intern=TRUE)
    # now check first lines to be sure that they are actually amiramesh
    Filter(function(f) isTRUE(amiratype(f)=='SkeletonGraph'),ff)
  } else character(0)
}

test.igraphvsoriginal<-function(){
  ff=file.path(TestDir,'Data','neurons',
               c("Neurites.am", "NeuritesWithIsolatedPoints.am", 
                 "NeuritesWithIsolatedSegment_veryshort.am",
                 "NeuritesWithNA.am",
                 "neuron_with_materials.am",
                 "testneuron_am3d_ascii.am","testneuron_am3d.am",
                 "UnbranchedNeurite.am"))
  ff=find.am3d()
  exceptions=c("NeuritesWithNA.am",'120818-JK56fill_skel.am')
  
  require(plyr)
  l_ply(ff[!basename(ff)%in%exceptions], igraphvsoriginal.amiraneuron,
        .progress='text')
}
