# Differential Geometry in Edge Detection: accurate estimation of position, orientation and curvature

Author of this Release Package: 
	Yuliang Guo (yuliang_guo@brown.edu)
	Benjamin Kimia (benjamin_kimia@brown.edu)

Reference: 
	"Differential Geometry in Edge Detection: accurate estimation of position, orientation and curvature"
	 B B Kimia, X Li, Y Guo, A Tamrakar
	IEEE Transactions on Pattern Analysis and Machine Intelligence 2018

### 1. Setup
In matlab, run
cd curvature_estimation/curvelet
mex  form_curvelet_mex.cxx CC_curve_model_3d.cxx curvelet.cxx curveletmap.cxx form_curvelet_process.cxx

### 2. Run

see demo.m

### 3. Improve a third-party edge detector

the last section of code in demo.m give an example of using this software to perform accurate edge detection with supixel localization and highly accurate orientation based on a third-party edge detector

Reference:
	"Fast Edge Detection Using Structured Forests"
	Piotr Dollar and C. Lawrence Zitnick
	IEEE Transactions on Pattern Analysis and Machine Intelligence 2015
