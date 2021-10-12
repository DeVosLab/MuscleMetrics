/*	**********************************************

 	MuscleMetrics.ijm

	**********************************************
	
	Author: 			Winnok H. De Vos
	Date Created: 		November 25, 2016
	Date Last Modified:	December 20, 2016

 	Interactive tools to facilitate assessment of muscle organisation
 	Includes a segmentation of connected structures and measurement of size, shape and local thickness
 	Also includes an option to measure directionality 
 	Segmentation of muscle fibers relies on a directional derivative for which the plugin FeatureJ is used (written by E. Meijiering, https://imagescience.org/meijering/software/featurej/)
 	Elliptic fourier descriptors requires an additional plugin written by Thomas Boudier and Ben Tupper (http://imagejdocu.tudor.lu/doku.php?id=plugin:analysis:fourier_shape_analysis:start)
 	
 	Written for Ineke D'Hondt and Bart Braeckman, Ugent
 	
	**********************************************
*/

/*
 	***********************

	Variable initiation

	***********************
*/

//	Numbers
var line_thickness 			= 	200;										//	thickness of selection for straightening
var gauss_scale				= 	1;											//	gauss/laplace scale
var min_size				= 	50;											//	minimum structure size
var max_size				= 	1000000;									//	maximum structure size
var pixel_size 				= 	0.15;										//	pixel size in micron	
var enhance_radius			=	10;											//	radius for local contrast enhancement
var bg_radius				= 	0;											//	radius for background correction
		
// Strings
var threshold_method		= 	"Yen";										//	autothreshold method
var micron					= 	getInfo("micrometer.abbreviation");			//	micron symbol
var enhancer				= 	"laplace";									//	muscle fiber enhancement algorithm

// Arrays
var thresholds 				= 	getList("threshold.methods");				//	autothreshold algorithms 
var enhancers				= 	newArray("gauss","laplace");				//	muscle fiber enhancement algorithms

// Booleans
var thickness				=	true;										//	include local thickness measurement
var skeleton				=	true;										//	include skeleton analysis
var efd						=	false;										//	include elliptic fourier descriptor extraction
var corner					= 	true;										//	include corner detection

/*
 	***********************

	Macro Tools

	***********************
*/

macro "Autorun"
{
	erase();
	setOptions();
}

macro "Settings Action Tool - C438 T5f16S"
{
	setOptions()
	setup();
}

macro "Sarcomere Selection Action Tool - C333 R33dd"
{
	erase();
	run("Set Measurements...", " fit redirect=None decimal=4");
	id = getImageID;
	title = getTitle;
	if(selectionType()==-1)
	{
		setTool("polygon");
		waitForUser("First Draw a Polygon Around One Sarcomer and Press OK"); 
	}
	run("Duplicate...","title=["+title+"] duplicate channels=1");
	sarcomere_id = getImageID;
	run("Clear Outside", "slice");
	run("Measure");
	angle = getResult("Angle",0);
	run("Select None");
	run("Rotate... ", "angle="+angle+" grid=1 interpolation=Bilinear enlarge");
	selectImage(id); close;
	erase();
}

macro "Straighten Action Tool - C333 P0347b8fc0 "
{
	line_thickness = getNumber("Line Thickness",line_thickness);
	id = getImageID;
	if(selectionType()==-1)
	{
		setTool("polyline");
		waitForUser("First Draw a Segmented line and Press OK"); 
	}
	run("Properties... ", "  width="+line_thickness);
	run("Straighten...", "title=Straightened line="+line_thickness+" process");
	selectImage(id);
	close;
}

macro "Analyze Directionality Action Tool - C333 T3e16< T3e16>"
{
	id = getImageID;
	selectImage(id);
	run("Directionality", "method=[Fourier components] nbins=90 histogram=-90");
}

macro "Analyze Objects Action Tool - C333 T3f16A"
{
	erase();
	setBatchMode(true);
	id = getImageID;
	decalibrateImage(id);
	obj_nr = segmentObjects(id);
	if(obj_nr > 0) analyzeObjects(id,obj_nr);
	calibrateImage(id);
	setBatchMode("exit and display");
	toggleOverlay();
}

macro "Save And Close Action Tool -  C333 T3f16C"
{
	dir 	= getDirectory("");
	id 		= getImageID;
	title	= getTitle; 
	prefix 	= getString("Prefix", title);
	// save image
	selectImage(id);
	run("Remove Overlay");
	saveAs(".tif", dir+prefix+"_Cropped.tif");
	close;
	// save results
	selectWindow("Object_Results");
	saveAs("Results", dir+prefix+"_Object_Results.txt");
	run("Close");
	// save log
	selectWindow("Log");
	saveAs("Text", dir+prefix+"_Log.txt");
	run("Close");
	// save rois
	roiManager("deselect");
	roiManager("Save", dir+prefix+"_RoiSet.zip");
	erase();
}

macro "Toggle Overlay Action Tool - Caaa O11ee"
{
	toggleOverlay();
}

macro "[t] Toggle Overlay"
{
	toggleOverlay();
}

/*
 	***********************

	Functions

	***********************
*/

function erase()
{
	run("Clear Results");
	roiManager("reset");
	run("Collect Garbage");
}

function setOptions()
{
	run("Options...", "iterations = 1 count = 1");
	run("Colors...", "foreground=white background=black selection=yellow");
	run("Appearance...", "  antialiased menu = 0");
	run("Overlay Options...", "stroke = red width = 1 fill = none");
	setOption("BlackBackground", false);
	run("Set Measurements...", "area mean min fit perimeter shape feret's redirect=None decimal=4");
	run("Misc...", "divide=Infinity reverse");
	run("Conversions...", " ");
	setBackgroundColor(0, 0, 0);
	setForegroundColor(255,255,255);
}

function getMoment()
{
     MonthNames = newArray("Jan","Feb","Mar","Apr","May","Jun","Jul","Aug","Sep","Oct","Nov","Dec");
     DayNames = newArray("Sun", "Mon","Tue","Wed","Thu","Fri","Sat");
     getDateAndTime(year, month, dayOfWeek, dayOfMonth, hour, minute, second, msec);
     time_string ="Date: "+DayNames[dayOfWeek]+" ";
     if (dayOfMonth<10) {time_string = time_string+"0";}
     time_string = time_string+dayOfMonth+"-"+MonthNames[month]+"-"+year+"\nTime: ";
     if (hour<10) {time_string = time_string+"0";}
     time_string = time_string+hour+":";
     if (minute<10) {time_string = time_string+"0";}
     time_string = time_string+minute+":";
     if (second<10) {time_string = time_string+"0";}
     time_string = time_string+second;
     return time_string;
}

function setup()
{
	setOptions();
	Dialog.create("MuscleMetrics - Object Analysis");
	Dialog.addNumber("Pixel Size (ignored if calibrated)",pixel_size,3,5,micron);
	Dialog.addNumber("Background Correction Radius",bg_radius,0,5,"");
	Dialog.addNumber("Local Contrast Enhance Radius",enhance_radius,0,5,"");	
	Dialog.addChoice("Object Enhancement Algorithm",enhancers,enhancer);
	Dialog.addNumber("Gaussian Smoothing Radius",gauss_scale,0,5,"");
	Dialog.addChoice("Autothreshold Algorithm",thresholds,threshold_method);
	Dialog.addNumber("Minimum Object Area ",min_size,0,5,"pixels");
	Dialog.addNumber("Maximum Object Area ",max_size,0,5,"pixels");
	labels 			= newArray(4);		
	defaults 		= newArray(4);
	labels[0] 		= "Incl. Local Thickness";	
	defaults[0] 	= thickness;
	labels[1] 		= "Incl. Skeleton Analysis";			
	defaults[1] 	= skeleton;
	labels[2] 		= "Incl. Elliptic Descriptors";		
	defaults[2] 	= efd;
	labels[3] 		= "Incl. Corner Detection";		
	defaults[3] 	= corner;
	Dialog.addCheckboxGroup(2,2,labels,defaults);
	Dialog.show;
	print("\\Clear");
	moment = getMoment();
	print(moment+"\n");
	print("*******************************************");
	print("               Settings");
	print("*******************************************");
	pixel_size				= Dialog.getNumber();	print("image pixel size:",pixel_size);
	bg_radius				= Dialog.getNumber();	print("backkground correction radius:",bg_radius);
	enhance_radius			= Dialog.getNumber();	print("local contrast enhancement radius:",enhance_radius);
	enhancer				= Dialog.getChoice();	print("object enhancer:",enhancer);
	gauss_scale				= Dialog.getNumber();	print("gauss smoothing scale:",gauss_scale);
	threshold_method		= Dialog.getChoice();	print("autothreshold method:",threshold_method);
	min_size 				= Dialog.getNumber();	print("min object size:",min_size);
	max_size 				= Dialog.getNumber();	print("max object size:",max_size);
	thickness				= Dialog.getCheckbox();	print("include local thickness:",thickness);
	skeleton				= Dialog.getCheckbox();	print("include skeletonization:",skeleton);
	efd						= Dialog.getCheckbox();	print("include EFD:",efd);
	corner					= Dialog.getCheckbox();	print("include corner detection:",corner);
	print("*******************************************");	
}

function calibrateImage(id)
{
	selectImage(id);
	getPixelSize(unit, pixelwidth, pixelheight);
	if(indexOf(unit,"pixel")>=0)run("Properties...", " unit=µm pixel_width="+pixel_size+" pixel_height="+pixel_size);
	else pixel_size = pixelwidth;
	print("Image Calibrated, Calibration =",pixel_size,micron,"per pixel");
}

function decalibrateImage(id)
{
	selectImage(id);
	getPixelSize(unit, pixel_width, pixel_height);
	if(indexOf(unit,"pixel")==0)pixel_size = pixel_width;
	run("Properties...", " unit=pixel pixel_width=1 pixel_height=1");
	print("Image Decalibrated, Calibration = 1 pixel per pixel");
}

function toggleOverlay()
{
	roiManager("Show All without labels");
	roiManager("Show None");
	run("Select None"); roiManager("deselect");
	if(Overlay.size == 0)run("From ROI Manager");
	else run("Remove Overlay");
}

function segmentObjects(id)
{
	// segment individual muscle fibers and return ROIset
	selectImage(id); 
	//	pre-processing
	if(enhance_radius>0)run("Enhance Local Contrast (CLAHE)", "blocksize="+enhance_radius+" histogram=256 maximum=2 mask=*None* fast_(less_accurate)");
	if(bg_radius>0)run("Subtract Background...", "rolling="+bg_radius);
	//	object enhancement
	if(enhancer == "gauss")
	{
		run("Duplicate...","title=Enhanced");
		enh_id = getImageID;
		selectImage(enh_id);
		run("Gaussian Blur...", "sigma="+gauss_scale);
		run("Invert");
	}
	else if(enhancer == "laplace")
	{
		run("FeatureJ Derivatives", "x-order=0 y-order=2 z-order=0 smoothing="+gauss_scale);
		enh_id = getImageID;
		selectImage(enh_id); 
		rename("Enhanced");
	}
	//	binarization and object detection
	selectImage(enh_id); 
	setAutoThreshold(""+threshold_method);
	setOption("BlackBackground", false);
	run("Convert to Mask");
	run("Analyze Particles...", "size="+min_size+"-"+max_size+"pixel show=Nothing clear include add");
	selectImage(enh_id); close;
	obj_nr = roiManager("count");
	return(obj_nr);
}	

function analyzeObjects(id,obj_nr)
{
	//	create binary image of objects
	selectImage(id);
	getDimensions(im_width, im_height, channels, slices, frames);
	newImage("Objects", "8-bit black", im_width, im_height, 1);
	obj_id = getImageID;
	selectImage(obj_id);
	roiManager("Deselect");
	roiManager("Fill");
	setThreshold(1,255);
	setOption("BlackBackground", false);
	run("Convert to Mask");
	//	elliptic fourier descriptors (limited version to detect deviation from very modest approximation, i.e. using 5 descriptors)
	if(efd)
	{
		run("Set Measurements...", " perimeter redirect=None decimal=4");
		efd_sum	= newArray(obj_nr);
		efd_perim = newArray(obj_nr);
		run("Clear Results");
		selectImage(obj_id);
		for(i = 0; i < obj_nr; i++)
		{
			roiManager("select",i);
			//run("Elliptic Fourier D.", "number=10 results");
			run("Elliptic Fourier D.", "number=5 reconstruction results");
			IJ.renameResults("Results-EFD-Objects","Results");
			result_nr = nResults;
			for(v = 2; v < result_nr; v++) 	// first two coeffs can be ignored
			{	
				efd_sum[i] += getResult("efd",v);
			}
			run("Clear Results"); 
			selectWindow("Objects-EFD");
			run("Measure"); 
			close;
			efd_perim[i] = getResult("Perim.",0);
			selectWindow("Results"); 
			run("Close"); 
			selectImage(obj_id); 
			run("Select None");		
		}
	}
	//	measure local thickness
	if(thickness)
	{
		run("Set Measurements...", "area standard mean redirect=None decimal=4");
		av_thickness = newArray(obj_nr); 	// average thickness
		sd_thickness = newArray(obj_nr); 	// standard deviation
		cv_thickness = newArray(obj_nr); 	// coeff. of variation	
		run("Clear Results");
		selectImage(obj_id);
		run("Local Thickness (complete process)", "threshold=1");
		selectWindow("Objects_LocThk"); 
		thick_id = getImageID;
		selectImage(thick_id);
		roiManager("Deselect");
		roiManager("Measure");
		for(i = 0; i < obj_nr; i++)
		{
			av_thickness[i] = getResult("Mean",i);
			sd_thickness[i]	= getResult("StdDev",i);
			cv_thickness[i]	= getResult("StdDev",i)/getResult("Mean",i);
		}
		selectImage(thick_id); close;
		run("Clear Results"); run("Collect Garbage");
	}
	// analyze skeleton
	if(skeleton)
	{
		skeleton_area 		= newArray(obj_nr); 	// skeleton length
		skeleton_mean 		= newArray(obj_nr); 	// mean intensity 
		branch_point_nr 	= newArray(obj_nr); 	// nr of branch points
		skeleton_feret		= newArray(obj_nr); 	// skeleton feret diameter
		straightness		= newArray(obj_nr); 	// straightness of the skeleton (1 is a straight line)
		selectImage(obj_id);
		run("Select None");
		run("Skeletonize");
		run("Analyze Skeleton (2D/3D)", "prune=none");
		tag_id = getImageID;
		selectImage(tag_id);
		run("Macro...", "code=if(v!=70)v=0"); // 70 = branch points
		run("Clear Results");
		run("Set Measurements...", "area mean feret's redirect=[Tagged skeleton] decimal=4");
		selectImage(obj_id);
		run("Analyze Particles...", "size=3-Inf pixel show=Nothing display clear include");
		for(i = 0; i < obj_nr; i++)
		{
			skeleton_area[i] 	= getResult("Area",i);
			branch_point_nr[i]	= getResult("Mean",i)*getResult("Area",i)/70; 
			skeleton_feret[i]	= getResult("Feret",i);
			straightness[i]		= getResult("Feret",i)/getResult("Area",i); 
		}
		selectImage(tag_id); close;
		run("Clear Results"); run("Collect Garbage");
	}
	// detect bending points (corners)
	if(corner)
	{
		bending_energy = newArray(obj_nr);
		corner_nr = newArray(obj_nr);
		run("Set Measurements...", "mean integrated redirect=None decimal=4");
		selectImage(obj_id); // now the objects have been skeletonized
		run("FeatureJ Structure", "smallest smoothing=1.0 integration=1");
		struct_id = getImageID;
		selectImage(struct_id);
		roiManager("Deselect");
		roiManager("Measure");
		for(i = 0; i < obj_nr; i++)
		{
			bending_energy[i] 	= getResult("RawIntDen",i);
		}
		run("Clear Results"); 
		selectImage(struct_id);
		run("Find Maxima...", "noise=50 output=[Single Points]");
		selectWindow("Objects smallest structure eigenvalues Maxima");
		max_id = getImageID;
		selectImage(struct_id); close;
		selectImage(max_id);
		roiManager("Deselect");
		roiManager("Measure");
		for(i = 0; i < obj_nr; i++)
		{
			corner_nr[i] 	= getResult("RawIntDen",i)/255;
		}
		selectImage(max_id); close;
		run("Clear Results"); run("Collect Garbage");
	}
	//	measure general object parameters
	run("Set Measurements...", "area mean integrated standard fit perimeter shape feret's redirect=None decimal=4");
	roiManager("Deselect");
	selectImage(id);
	roiManager("Measure");
	// measure objects and complete the results table
	for(i = 0; i < obj_nr; i++)
	{
		if(efd)
		{	 
			setResult("EFD_Sum",i,efd_sum[i]); 
			setResult("Curvature",i,efd_perim[i]/getResult("Perim.",i));
		}
		if(thickness)
		{
			setResult("Mean_Thickness",i,av_thickness[i]);  
			setResult("Std_Thickness",i,sd_thickness[i]);  
			setResult("CoV_Thickness",i,cv_thickness[i]); 		 
		}
		if(skeleton)
		{
			setResult("Skeleton_Length",i,skeleton_area[i]); 
			setResult("Skeleton_Mean",i,skeleton_mean[i]);
			setResult("Branch_Point_Nr",i,branch_point_nr[i]);  
			setResult("Skeleton_Feret",i,skeleton_feret[i]); 
			setResult("Straightness",i,straightness[i]); 
		}
		if(corner)
		{
			setResult("Bending_Energy",i,bending_energy[i]);
			setResult("Bend_Point_Nr",i,corner_nr[i]); 
		}
	}
	updateResults;
	IJ.renameResults("Results","Object_Results");
	selectImage(obj_id); close;	
}

