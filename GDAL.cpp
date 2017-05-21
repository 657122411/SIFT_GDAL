// GDAL.cpp : 定义控制台应用程序的入口点。
//

#include "stdafx.h"
#include "gdal_priv.h"
#pragma comment( lib, "opencv_highgui249d.lib")
#pragma comment( lib, "opencv_core249d.lib")
 
#include "ASiftDetector.h"
#include "utils.h"

// 提取sift特征
Mat ExtractSIFTFeature(const string &imgfn, vector<KeyPoint> &points)  {
    auto img = imread(imgfn, true);    //imgfn: image file name
    SiftFeatureDetector detector;
    vector<KeyPoint> keypoints;
    detector.detect(img, keypoints);
    SiftDescriptorExtractor extractor;
    Mat descriptors;
    if (!keypoints.size()) {
        return Mat();
    }
    // extract sift feature
    extractor.compute(img, keypoints, descriptors);
    //points.resize(keypoints.size());
    for (int i = 0; i < keypoints.size(); i++){
        points.push_back(keypoints[i]);
    }
    return descriptors;
}

/*
    //注册文件格式
    GDALAllRegister();
   
    const char* pszFile = "C:\\Users\\Xyy\\Desktop\\SIFT\\GF_L2A.tif";
    GDALDataset *poDataset;
    //使用只读方式打开图像
    poDataset = (GDALDataset*) GDALOpen( pszFile,GA_ReadOnly );
    if( poDataset == NULL )
    {
        printf( "File: %s不能打开！\n",pszFile);
        return 0;
    }
 
    //输出图像的格式信息
    printf( "Driver:%s/%s\n",
        poDataset->GetDriver()->GetDescription(),
        poDataset->GetDriver()->GetMetadataItem( GDAL_DMD_LONGNAME) );
 
    //输出图像的大小和波段个数
    printf( "Size is%dx%dx%d\n",
        poDataset->GetRasterXSize(),poDataset->GetRasterYSize(),
        poDataset->GetRasterCount());
 
    //输出图像的投影信息
    if( poDataset->GetProjectionRef() != NULL )
        printf( "Projectionis `%s'\n", poDataset->GetProjectionRef() );
 
    //输出图像的坐标和分辨率信息
    double adfGeoTransform[6];
    if( poDataset->GetGeoTransform( adfGeoTransform) == CE_None )
    {
        printf( "Origin =(%.6f,%.6f)\n",
            adfGeoTransform[0], adfGeoTransform[3]);
 
        printf( "PixelSize = (%.6f,%.6f)\n",
            adfGeoTransform[1], adfGeoTransform[5]);
    }
 
    GDALRasterBand *poBand;
    int            nBlockXSize, nBlockYSize;
    int            bGotMin, bGotMax;
    double         adfMinMax[2];
 
    //读取第一个波段
    poBand = poDataset->GetRasterBand( 1 );
 
    //获取图像的块大小并输出
    poBand->GetBlockSize(&nBlockXSize, &nBlockYSize );
    printf( "Block=%dx%dType=%s, ColorInterp=%s\n",
        nBlockXSize, nBlockYSize,
        GDALGetDataTypeName(poBand->GetRasterDataType()),
        GDALGetColorInterpretationName(
        poBand->GetColorInterpretation()));
 
    //获取该波段的最大值最小值，如果获取失败，则进行统计
    adfMinMax[0] = poBand->GetMinimum( &bGotMin);
    adfMinMax[1] = poBand->GetMaximum( &bGotMax);
 
    if( ! (bGotMin&& bGotMax) )
        GDALComputeRasterMinMax((GDALRasterBandH)poBand, TRUE, adfMinMax);
 
    printf( "Min=%.3fd,Max=%.3f\n", adfMinMax[0], adfMinMax[1] );
 
    //输出图像的金字塔信息
    if( poBand->GetOverviewCount() > 0 )
        printf( "Band has%d overviews.\n", poBand->GetOverviewCount() );
 
    //输出图像的颜色表信息
    if( poBand->GetColorTable() != NULL)
        printf( "Band hasa color table with %d entries.\n",
        poBand->GetColorTable()->GetColorEntryCount() );
 
    float *pafScanline;
    int   nXSize = poBand->GetXSize();
   
    //读取图像的第一行数据
    pafScanline = (float*) CPLMalloc(sizeof(float)*nXSize);
    poBand->RasterIO(GF_Read, 0, 0, nXSize,1, 
        pafScanline, nXSize,1, GDT_Float32, 0, 0 );
 
    CPLFree(pafScanline);
 
    //关闭文件
    GDALClose((GDALDatasetH)poDataset);
	
	*/



int main(int argc, const char * argv[]) {
    // insert code here...
    string imgfn = "C:\\Users\\Xyy\\Desktop\\SIFT\\GF01_MS2_cut.tif";
    Mat queryImage, queryBlackImg, qDescriptor;
    queryImage = imread(imgfn);
    vector<KeyPoint> qKeypoints;
    qDescriptor = ExtractSIFTFeature(imgfn, qKeypoints);
    /*
	drawKeypoints(queryImage, qKeypoints, queryBlackImg);
    imshow("Sift", queryBlackImg);
    cvWaitKey(0);
	*/
    
    string objFileName = "C:\\Users\\Xyy\\Desktop\\SIFT\\15.tif";
    Mat objectImage, objBlackImg, objDesriptor;
    objectImage = imread(objFileName);
    vector<KeyPoint> objKeypoints;
    objDesriptor = ExtractSIFTFeature(objFileName, objKeypoints);
    /*
	drawKeypoints(objectImage, objKeypoints, objBlackImg);
    imshow("Sift", objBlackImg);
    cvWaitKey(0);
	*/
    
    //直接找1近邻，SIFT找匹配点
    FlannBasedMatcher matcher;
	
    vector< DMatch > matches1;
    matcher.match(qDescriptor, objDesriptor, matches1);
    findInliers(qKeypoints, objKeypoints, matches1, imgfn, objFileName); // 通过1近邻匹配的点对
    
    // 使用Lowe的raw feature方法匹配
    vector<Point2f> qeK, obK;
    vector<vector<DMatch>> matches;
    vector<DMatch> good_matches2;
    matcher.knnMatch(qDescriptor, objDesriptor, matches, 2);
    for (size_t i = 0; i < matches.size(); i++){
        if (matches[i][0].distance < 0.8*matches[i][1].distance){
            good_matches2.push_back(matches[i][0]);
            qeK.push_back(qKeypoints[matches[i][0].queryIdx].pt);
            obK.push_back(objKeypoints[matches[i][1].trainIdx].pt);
        }
    }
    findInliers(qKeypoints, objKeypoints, good_matches2, imgfn, objFileName); // 通过Lowe匹配的点对

    
    /*
	ASiftDetector asiftDetector;
    vector<KeyPoint> asiftKeypoints_query;
    Mat asiftDescriptors_query;
    asiftDetector.detectAndCompute(queryImage, asiftKeypoints_query, asiftDescriptors_query);
    vector<KeyPoint> asiftKeypoints_object;
    Mat asiftDescriptors_object;
    asiftDetector.detectAndCompute(objectImage, asiftKeypoints_object, asiftDescriptors_object);*/
    
    /*Matching descriptor vectors using FLANN matcher, ASIFT找匹配点
    std::vector< DMatch > asiftMatches;
    matcher.match(asiftDescriptors_query, asiftDescriptors_object, asiftMatches);
    findInliers(asiftKeypoints_query, asiftKeypoints_object, asiftMatches, imgfn, objFileName);*/
    
    // 使用内置函数画匹配点对  Mat img_matches;
     /*
	 drawMatches( queryImage, qKeypoints, objectImage, objKeypoints,  matches, img_matches, Scalar::all(-1), Scalar::all(-1),
     vector<char>(), DrawMatchesFlags::NOT_DRAW_SINGLE_POINTS );
     imshow("Good Matches & Object detection", img_matches);
     waitKey(0);
	 */
    
    return 0;
}

	
	
		
		   
	/*vector<Point2f> qeK, obK;
    vector<vector<DMatch>> matches;
    vector<DMatch> good_matches2;
    matcher.knnMatch(qDescriptor, objDesriptor, matches, 2);
    for (size_t i = 0; i < matches.size(); i++){
        if (matches[i][0].distance < 0.8*matches[i][1].distance){
            good_matches2.push_back(matches[i][0]);
            qeK.push_back(result[matches[i][0].queryIdx].pt);

			//qeK.push_back(qKeypoints[matches[i][0].queryIdx].pt);

            obK.push_back(objKeypoints[matches[i][1].trainIdx].pt);
        }
    }
    findInliers(result, objKeypoints, good_matches2, imgfn, objFileName); // 通过Lowe匹配的点对
	*/
	/*findInliers(qKeypoints, objKeypoints, good_matches2, imgfn, objFileName); // 通过Lowe匹配的点对
		vector< DMatch > matches1;
		Mat testdescriptors;
		SiftFeatureDetector detector;
		SiftDescriptorExtractor extractor;
		extractor.compute(queryImage, result, testdescriptors);
		matcher.match(testdescriptors, objDesriptor, matches1);
		*/
    	/*
    vector< DMatch > matches1;
    matcher.match(qDescriptor, objDesriptor, matches1);
    findInliers(qKeypoints, objKeypoints, matches1, imgfn, objFileName); // 通过1近邻匹配的点对
    */
    /*ASiftDetector asiftDetector;
    vector<KeyPoint> asiftKeypoints_query;
    Mat asiftDescriptors_query;
    asiftDetector.detectAndCompute(queryImage, asiftKeypoints_query, asiftDescriptors_query);
    vector<KeyPoint> asiftKeypoints_object;
    Mat asiftDescriptors_object;
    asiftDetector.detectAndCompute(objectImage, asiftKeypoints_object, asiftDescriptors_object);*/
    
    /*Matching descriptor vectors using FLANN matcher, ASIFT找匹配点
    std::vector< DMatch > asiftMatches;
    matcher.match(asiftDescriptors_query, asiftDescriptors_object, asiftMatches);
    findInliers(asiftKeypoints_query, asiftKeypoints_object, asiftMatches, imgfn, objFileName);*/
    
    // 使用内置函数画匹配点对  Mat img_matches;
     /*
	 drawMatches( queryImage, qKeypoints, objectImage, objKeypoints,  matches, img_matches, Scalar::all(-1), Scalar::all(-1),
     vector<char>(), DrawMatchesFlags::NOT_DRAW_SINGLE_POINTS );
     imshow("Good Matches & Object detection", img_matches);
     waitKey(0);
	 */
	 /*drawKeypoints(queryImage, qKeypoints, queryBlackImg);
    imshow("Sift", queryBlackImg);
    cvWaitKey(0);*/
		