
Toolkits:

(1) ITK v4.5.1 x64[CMAKE]
ON: Module_ITKReview; ITK_BUILD_DEFAULT_MODULES
ON: ITK_USE_SYSTEM_HDF5
OTHERS: DEFAULT

(2) VTK v5.10.1 x64[CMAKE]
ON: VTK_USE_SYSTEM_HDF5
OTHERS: DEFAULT

(3) QT 4.8.2 x64


===========================================================					   


Configuration of LiverSegASM:

[CMAKE] OFF: BUILD_CUDA_LIB; ITK_USE_GPU


===========================================================
Change log 2014-05-19
1. Change HDF5 libraries dependency. 
Use ITK'S HDF5 instead of standalone HD5 libraries. So when compile the LiverSegmentation.exe, standalone HDF5 will no longer be necessary, because we use the HDF5 of ITK.
2. Change boundary profile classification using Adaboost instead of FLANN.
In this way, the speed of whole process improved 2~3 times. So using this version, it will take 50~80 seconds more or less instead of up to 200 seconds.
3. Extract profile including:
	(1)contrast
	(2)intensity moment
	(3)gradient moment
	(4)gradient
	(5)offset from liver centroid
Using the composed features above, the impact of tumours can be reduced to the minimum.
4. Using deformable mesh instead of boundary search strategy, and then using this deformed mesh to guide the fitting of statistical shape model.

===========================================================
Run command:
LiverSegmentation.exe [inputImage] [SSMFile] [OutputDir] [AdaboostClassifierFile] [AdaboostSegmentFile] [GeometryFile] [AtlasImage]

For example:
LiverSegmentation.exe liver-orig1.nii.gz Data\StatisticalShapeModel_20131212.h5 D:\Output Data\AdaboostClassifier-default.h5 Data\learn3-AdaBoostResults Data\geoImage.mha Data\atlasImage.nii.gz